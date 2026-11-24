#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions Configure);
Configure(qw(no_ignore_case no_auto_abbrev));
use File::Basename qw(dirname);
use File::Path qw(make_path);
use Cwd qw(abs_path);

# Cluster valid deletion-supporting rows directly from the annotated
# split_reads.tsv. No per-read assignment file is generated, and the input
# split_reads.tsv is never modified.

my ($config, $input, $clusters_out, $help);
GetOptions(
    'Config=s'      => \$config,
    'Input=s'       => \$input,
    'ClustersOut=s' => \$clusters_out,
    'Help'          => \$help,
) or die usage();

if ($help) {
    print usage();
    exit 0;
}

die usage() unless defined $config && defined $input && defined $clusters_out;
$config = abs_path($config) if -e $config;
$input  = abs_path($input)  if -e $input;

die "[ERROR] Config file not found: $config\n" unless defined $config && -s $config;
die "[ERROR] Input file not found: $input\n" unless defined $input && -s $input;
die "[ERROR] Config and input paths must be absolute\n"
    unless $config =~ m{^/} && $input =~ m{^/};
die "[ERROR] Output path must be absolute\n" unless $clusters_out =~ m{^/};
make_path(dirname($clusters_out)) unless -d dirname($clusters_out);

my %conf = read_config($config);
my $cluster_window = required_config(\%conf, 'SA_SPLIT_CLUSTER_WINDOW');
my $min_support_reads = required_config(\%conf, 'SA_SPLIT_MIN_SUPPORT_READS');
validate_nonnegative_integer('SA_SPLIT_CLUSTER_WINDOW', $cluster_window);
validate_positive_integer('SA_SPLIT_MIN_SUPPORT_READS', $min_support_reads);

my ($records_ref, $input_count, $filtered_count) = read_split_reads($input);
my @clusters = cluster_records($records_ref, $cluster_window);

for my $cluster (@clusters) {
    $cluster->{Cluster_Status} =
        $cluster->{Support_Reads} >= $min_support_reads ? 'Pass' : 'LowSupport';
}

write_clusters($clusters_out, \@clusters);

print STDERR "[INFO] Input split-read records : $input_count\n";
print STDERR "[INFO] Filtered non-DEL records: $filtered_count\n";
print STDERR "[INFO] Valid DEL records       : ", scalar(@$records_ref), "\n";
print STDERR "[INFO] Clusters                : ", scalar(@clusters), "\n";
print STDERR "[INFO] PASS clusters           : ",
    scalar(grep { $_->{Cluster_Status} eq 'Pass' } @clusters), "\n";
print STDERR "[INFO] LowSupport clusters     : ",
    scalar(grep { $_->{Cluster_Status} eq 'LowSupport' } @clusters), "\n";
print STDERR "[INFO] Cluster output          : $clusters_out\n";
exit 0;

# ------------------------------------------------------------------
# Input reader
# ------------------------------------------------------------------

sub read_split_reads {
    my ($file) = @_;
    open my $fh, '<', $file or die "[ERROR] Cannot read $file: $!\n";
    my $header_line = <$fh>;
    die "[ERROR] Empty split-read file: $file\n" unless defined $header_line;
    chomp $header_line;
    $header_line =~ s/\r$//;

    my @header = split /\t/, $header_line, -1;
    my %index;
    for my $i (0 .. $#header) {
        die "[ERROR] Duplicate input column '$header[$i]' in $file\n"
            if exists $index{$header[$i]};
        $index{$header[$i]} = $i;
    }

    for my $column (qw/
        SampleID Read_Name Chrom Start End Size_bp
        Junction_Type Breakpoint_Status Breakpoint_Reason
    /) {
        die "[ERROR] Missing column '$column' in $file\n" unless exists $index{$column};
    }

    my @records;
    my $input_count = 0;
    my $filtered_count = 0;
    my $line_no = 1;
    while (my $line = <$fh>) {
        $line_no++;
        chomp $line;
        $line =~ s/\r$//;
        next if $line eq '';
        $input_count++;

        my @values = split /\t/, $line, -1;
        die "[ERROR] Column count mismatch at $file:$line_no\n"
            unless @values == @header;
        my %record;
        @record{@header} = @values;

        # Only PASS deletion gaps enter clustering. All other rows remain in
        # split_reads.tsv as read-level audit evidence.
        unless ($record{Breakpoint_Status} eq 'PASS'
            && $record{Junction_Type} eq 'Deletion_Gap'
            && $record{Breakpoint_Reason} eq 'Valid_Deletion_Gap') {
            $filtered_count++;
            next;
        }

        unless (is_positive_integer($record{Start})
            && is_positive_integer($record{End})
            && $record{End} >= $record{Start}
            && is_positive_integer($record{Size_bp})
            && $record{Size_bp} == $record{End} - $record{Start} + 1
            && defined $record{Chrom} && $record{Chrom} ne '' && $record{Chrom} ne 'NA'
            && defined $record{SampleID} && $record{SampleID} ne '' && $record{SampleID} ne 'NA'
            && defined $record{Read_Name} && $record{Read_Name} ne '' && $record{Read_Name} ne 'NA') {
            die "[ERROR] Invalid PASS deletion record at $file:$line_no\n";
        }

        push @records, \%record;
    }
    close $fh or die "[ERROR] Failed while reading $file\n";

    @records = sort {
        $a->{SampleID} cmp $b->{SampleID}
            || chrom_key($a->{Chrom}) cmp chrom_key($b->{Chrom})
            || $a->{Start} <=> $b->{Start}
            || $a->{End} <=> $b->{End}
            || $a->{Read_Name} cmp $b->{Read_Name}
    } @records;

    return (\@records, $input_count, $filtered_count);
}

# ------------------------------------------------------------------
# Clustering
# ------------------------------------------------------------------

sub cluster_records {
    my ($records_ref, $window) = @_;
    return () unless @$records_ref;

    my (@clusters, @current);
    for my $record (@$records_ref) {
        if (!@current) {
            @current = ($record);
            next;
        }

        my $current_sample = $current[0]{SampleID};
        my $current_chrom  = $current[0]{Chrom};
        if ($record->{SampleID} ne $current_sample || $record->{Chrom} ne $current_chrom) {
            push @clusters, build_cluster(\@current);
            @current = ($record);
            next;
        }

        my $median_start = median_int(map { $_->{Start} } @current);
        my $median_end   = median_int(map { $_->{End} } @current);

        if (abs($record->{Start} - $median_start) <= $window
            && abs($record->{End} - $median_end) <= $window) {
            push @current, $record;
        }
        else {
            push @clusters, build_cluster(\@current);
            @current = ($record);
        }
    }

    push @clusters, build_cluster(\@current) if @current;

    for my $i (0 .. $#clusters) {
        $clusters[$i]{Cluster_ID} = sprintf('SR_CLUSTER_%06d', $i + 1);
    }
    return @clusters;
}

sub build_cluster {
    my ($records_ref) = @_;
    my @starts = map { $_->{Start} } @$records_ref;
    my @ends   = map { $_->{End} } @$records_ref;
    my @sizes  = map { $_->{Size_bp} } @$records_ref;
    my %read_names = map { $_->{Read_Name} => 1 } @$records_ref;

    return {
        SampleID       => $records_ref->[0]{SampleID},
        Chrom          => $records_ref->[0]{Chrom},
        Start          => median_int(@starts),
        End            => median_int(@ends),
        Support_Reads  => scalar(keys %read_names),
        Support_Records=> scalar(@$records_ref),
        Start_Min      => min_num(@starts),
        Start_Max      => max_num(@starts),
        End_Min        => min_num(@ends),
        End_Max        => max_num(@ends),
        Mean_Size_bp   => mean_num(@sizes),
        Size_Min       => min_num(@sizes),
        Size_Max       => max_num(@sizes),
    };
}

# ------------------------------------------------------------------
# Output
# ------------------------------------------------------------------

sub write_clusters {
    my ($file, $clusters_ref) = @_;
    my $tmp_file = "$file.tmp.$$";
    unlink $tmp_file if -e $tmp_file;
    open my $fh, '>', $tmp_file or die "[ERROR] Cannot write $tmp_file: $!\n";

    my @header = qw(
        SampleID Cluster_ID Chrom Start End Size_bp Support_Reads Support_Records
        Start_Range End_Range Mean_Size_bp Size_Range_bp Cluster_Status
    );
    print {$fh} join("\t", @header), "\n";

    for my $cluster (@$clusters_ref) {
        my @values = (
            $cluster->{SampleID},
            $cluster->{Cluster_ID},
            $cluster->{Chrom},
            $cluster->{Start},
            $cluster->{End},
            $cluster->{End} - $cluster->{Start} + 1,
            $cluster->{Support_Reads},
            $cluster->{Support_Records},
            $cluster->{Start_Min} . '-' . $cluster->{Start_Max},
            $cluster->{End_Min} . '-' . $cluster->{End_Max},
            sprintf('%.2f', $cluster->{Mean_Size_bp}),
            $cluster->{Size_Min} . '-' . $cluster->{Size_Max},
            $cluster->{Cluster_Status},
        );
        print {$fh} join("\t", @values), "\n";
    }
    close $fh or die "[ERROR] Cannot close $tmp_file: $!\n";
    rename $tmp_file, $file
        or die "[ERROR] Cannot replace cluster output $file: $!\n";
}

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------

sub median_int {
    my @values = sort { $a <=> $b } @_;
    return 0 unless @values;
    return int($values[int(@values / 2)] + 0.5) if @values % 2;
    return int(($values[@values / 2 - 1] + $values[@values / 2]) / 2 + 0.5);
}

sub mean_num {
    my @values = @_;
    return 0 unless @values;
    my $sum = 0;
    $sum += $_ for @values;
    return $sum / @values;
}

sub min_num {
    my @values = @_;
    return 0 unless @values;
    my $minimum = $values[0];
    for my $value (@values) { $minimum = $value if $value < $minimum; }
    return $minimum;
}

sub max_num {
    my @values = @_;
    return 0 unless @values;
    my $maximum = $values[0];
    for my $value (@values) { $maximum = $value if $value > $maximum; }
    return $maximum;
}

sub chrom_key {
    my ($chrom) = @_;
    $chrom =~ s/^chr//i;
    return sprintf('%09d', $chrom) if $chrom =~ /^\d+$/;
    return 'X' if uc($chrom) eq 'X';
    return 'Y' if uc($chrom) eq 'Y';
    return 'Z' . $chrom;
}

sub read_config {
    my ($file) = @_;
    my %config;
    open my $fh, '<', $file or die "[ERROR] Cannot read $file: $!\n";
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next if $line =~ /^\s*(?:#|$)/;
        $line =~ s/\s+#.*$//;
        next unless $line =~ /^\s*([^=\s]+)\s*=\s*(.*?)\s*$/;
        die "[ERROR] Duplicate config key $1\n" if exists $config{$1};
        $config{$1} = $2;
    }
    close $fh;
    return %config;
}

sub required_config {
    my ($config_ref, $key) = @_;
    die "[ERROR] Missing config $key\n"
        unless defined $config_ref->{$key} && $config_ref->{$key} ne '';
    return $config_ref->{$key};
}

sub validate_positive_integer {
    my ($name, $value) = @_;
    die "[ERROR] $name must be a positive integer\n"
        unless $value =~ /^\d+$/ && $value > 0;
}

sub validate_nonnegative_integer {
    my ($name, $value) = @_;
    die "[ERROR] $name must be a non-negative integer\n"
        unless $value =~ /^\d+$/;
}

sub is_positive_integer {
    return defined $_[0] && $_[0] =~ /^\d+$/ && $_[0] > 0;
}

sub usage {
    return "Usage: perl ClusterSASplitReads.pl --Config CONF --Input split_reads.tsv --ClustersOut split_clusters.tsv\n";
}
