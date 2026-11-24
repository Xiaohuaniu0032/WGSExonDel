#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions Configure);
Configure(qw(no_ignore_case no_auto_abbrev));
use File::Basename qw(dirname);
use File::Path qw(make_path);
use Cwd qw(abs_path);

my ($config, $input, $clusters_out, $help);
GetOptions(
    'Config=s'      => \$config,
    'Input=s'       => \$input,
    'ClustersOut=s' => \$clusters_out,
    'Help'             => \$help,
) or die usage();

if ($help) {
    print usage();
    exit 0;
}

die usage() unless defined $config && defined $input && defined $clusters_out;

$config = abs_path($config) if -e $config;
$input  = abs_path($input)  if -e $input;

die "[ERROR] Config not found: $config\n" unless defined $config && -s $config;
die "[ERROR] Input not found: $input\n" unless defined $input && -s $input;
die "[ERROR] Config and input paths must be absolute\n"
    unless $config =~ m{^/} && $input =~ m{^/};

die "[ERROR] Output path must be absolute\n"
    unless $clusters_out =~ m{^/};
make_path(dirname($clusters_out)) unless -d dirname($clusters_out);

my %conf = read_config($config);
my $cluster_distance = required_config(\%conf, 'DISCORDANT_CLUSTER_DISTANCE');
my $min_support_pairs = required_config(\%conf, 'MIN_DISCORDANT_PAIRS');
validate_nonnegative_integer('DISCORDANT_CLUSTER_DISTANCE', $cluster_distance);
validate_positive_integer('MIN_DISCORDANT_PAIRS', $min_support_pairs);

my ($pairs_ref, $input_count) = read_pairs($input);
my @clusters = cluster_pairs($pairs_ref, $cluster_distance);

for my $cluster (@clusters) {
    $cluster->{Cluster_Status} =
        $cluster->{Support_Pairs} >= $min_support_pairs ? 'Pass' : 'LowSupport';
}

write_clusters($clusters_out, \@clusters);

print STDERR "[INFO] Input filtered pairs: $input_count\n";
print STDERR "[INFO] Unique pairs: ", scalar(@$pairs_ref), "\n";
print STDERR "[INFO] Clusters: ", scalar(@clusters), "\n";
print STDERR "[INFO] PASS clusters: ",
    scalar(grep { $_->{Cluster_Status} eq 'Pass' } @clusters), "\n";
print STDERR "[INFO] LowSupport clusters: ",
    scalar(grep { $_->{Cluster_Status} eq 'LowSupport' } @clusters), "\n";
print STDERR "[INFO] Clusters: $clusters_out\n";
exit 0;

# ------------------------------------------------------------------
# Input reader
# ------------------------------------------------------------------

sub read_pairs {
    my ($file) = @_;

    open my $fh, '<', $file or die "[ERROR] Cannot read $file: $!\n";
    my $header_line = <$fh>;
    die "[ERROR] Empty input file: $file\n" unless defined $header_line;
    chomp $header_line;
    $header_line =~ s/\r$//;

    my @header = split /\t/, $header_line, -1;
    my %index = map { $header[$_] => $_ } 0 .. $#header;
    for my $column (qw/Pair_ID SampleID Read_Name Chrom Start End Size_bp Insert_Size Pair_Orientation/) {
        die "[ERROR] Missing column '$column' in $file\n"
            unless exists $index{$column};
    }

    my (@pairs, %seen_pair);
    my $input_count = 0;

    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next if $line eq '';
        $input_count++;

        my @values = split /\t/, $line, -1;
        my %record = map {
            $_ => defined $values[$index{$_}] ? $values[$index{$_}] : 'NA'
        } @header;

        die "[ERROR] Invalid interval for $record{Pair_ID}\n"
            unless $record{Start} =~ /^\d+$/
                && $record{End} =~ /^\d+$/
                && $record{End} >= $record{Start};

        next if $seen_pair{$record{Pair_ID}}++;
        push @pairs, \%record;
    }
    close $fh;

    @pairs = sort {
        chrom_key($a->{Chrom}) cmp chrom_key($b->{Chrom})
            || $a->{Start} <=> $b->{Start}
            || $a->{End} <=> $b->{End}
            || $a->{Read_Name} cmp $b->{Read_Name}
    } @pairs;

    return (\@pairs, $input_count);
}

# ------------------------------------------------------------------
# Clustering
# ------------------------------------------------------------------

sub cluster_pairs {
    my ($pairs_ref, $distance) = @_;
    return () unless @$pairs_ref;

    my (@clusters, @current);
    for my $pair (@$pairs_ref) {
        if (!@current) {
            @current = ($pair);
            next;
        }

        my $current_chrom = $current[0]{Chrom};
        if ($pair->{Chrom} ne $current_chrom) {
            push @clusters, build_cluster(\@current);
            @current = ($pair);
            next;
        }

        my $median_start = median_int(map { $_->{Start} } @current);
        my $median_end   = median_int(map { $_->{End} } @current);

        if (abs($pair->{Start} - $median_start) <= $distance
            && abs($pair->{End} - $median_end) <= $distance) {
            push @current, $pair;
        }
        else {
            push @clusters, build_cluster(\@current);
            @current = ($pair);
        }
    }

    push @clusters, build_cluster(\@current) if @current;

    for my $i (0 .. $#clusters) {
        $clusters[$i]{Cluster_ID} = sprintf('DR_CLUSTER_%06d', $i + 1);
    }

    return @clusters;
}

sub build_cluster {
    my ($records_ref) = @_;
    my @records = @$records_ref;

    my @starts = map { $_->{Start} } @records;
    my @ends = map { $_->{End} } @records;
    my @insert_sizes = map { $_->{Insert_Size} } @records;

    my %read_names = map { $_->{Read_Name} => 1 } @records;
    my %orientation_count;
    $orientation_count{$_->{Pair_Orientation}}++ for @records;

    return {
        SampleID => $records[0]{SampleID},
        Chrom => $records[0]{Chrom},
        Start => median_int(@starts),
        End => median_int(@ends),
        Support_Pairs => scalar(keys %read_names),
        Start_Min => min_num(@starts),
        Start_Max => max_num(@starts),
        End_Min => min_num(@ends),
        End_Max => max_num(@ends),
        Median_Insert_Size => median_num(@insert_sizes),
        Mean_Insert_Size => mean_num(@insert_sizes),
        Min_Insert_Size => min_num(@insert_sizes),
        Max_Insert_Size => max_num(@insert_sizes),
        Pair_Orientation_Count => format_count_hash(\%orientation_count),
    };
}

# ------------------------------------------------------------------
# Output writers
# ------------------------------------------------------------------

sub write_clusters {
    my ($file, $clusters_ref) = @_;
    open my $fh, '>', $file or die "[ERROR] Cannot write $file: $!\n";

    my @header = qw(
        SampleID Cluster_ID Chrom Start End Size_bp Support_Pairs
        Start_Range End_Range Median_Insert_Size Mean_Insert_Size
        Insert_Size_Range Pair_Orientation_Count Cluster_Status
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
            $cluster->{Support_Pairs},
            $cluster->{Start_Min} . '-' . $cluster->{Start_Max},
            $cluster->{End_Min} . '-' . $cluster->{End_Max},
            sprintf('%.2f', $cluster->{Median_Insert_Size}),
            sprintf('%.2f', $cluster->{Mean_Insert_Size}),
            $cluster->{Min_Insert_Size} . '-' . $cluster->{Max_Insert_Size},
            $cluster->{Pair_Orientation_Count},
            $cluster->{Cluster_Status},
        );
        print {$fh} join("\t", @values), "\n";
    }
    close $fh;
}

# ------------------------------------------------------------------
# Numeric and configuration helpers
# ------------------------------------------------------------------

sub median_int {
    my @values = sort { $a <=> $b } @_;
    return 0 unless @values;
    return int($values[int(@values / 2)] + 0.5) if @values % 2;
    return int(($values[@values / 2 - 1] + $values[@values / 2]) / 2 + 0.5);
}

sub median_num {
    my @values = sort { $a <=> $b } @_;
    return 0 unless @values;
    return $values[int(@values / 2)] if @values % 2;
    return ($values[@values / 2 - 1] + $values[@values / 2]) / 2;
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
    for my $value (@values) {
        $minimum = $value if $value < $minimum;
    }
    return $minimum;
}

sub max_num {
    my @values = @_;
    return 0 unless @values;
    my $maximum = $values[0];
    for my $value (@values) {
        $maximum = $value if $value > $maximum;
    }
    return $maximum;
}

sub format_count_hash {
    my ($counts_ref) = @_;
    return join(',', map { $_ . ':' . $counts_ref->{$_} } sort keys %$counts_ref) || 'NA';
}

sub chrom_key {
    my ($chrom) = @_;
    (my $value = $chrom) =~ s/^chr//i;
    return sprintf('%09d', $value) if $value =~ /^\d+$/;
    return 'X' if uc($value) eq 'X';
    return 'Y' if uc($value) eq 'Y';
    return 'Z' . $value;
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

sub usage {
    return "Usage: perl ClusterDiscordantReads.pl --Config CONF --Input discordant_pairs.tsv --ClustersOut discordant_clusters.tsv\n";
}
