#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions Configure);
Configure(qw(no_ignore_case no_auto_abbrev));
use File::Basename qw(dirname);
use File::Path qw(make_path);
use Cwd qw(abs_path);

my ($config, $input, $assignments_out, $clusters_out, $help);
GetOptions(
    'Config=s'         => \$config,
    'Input=s'          => \$input,
    'AssignmentsOut=s' => \$assignments_out,
    'ClustersOut=s'    => \$clusters_out,
    'Help'             => \$help,
) or die usage();

if ($help) {
    print usage();
    exit 0;
}

die usage() unless defined $config && defined $input
    && defined $assignments_out && defined $clusters_out;

$config = abs_path($config) if -e $config;
$input  = abs_path($input)  if -e $input;

die "[ERROR] Config file not found: $config\n" unless defined $config && -s $config;
die "[ERROR] Input file not found: $input\n" unless defined $input && -s $input;
die "[ERROR] Config and input paths must be absolute\n"
    unless $config =~ m{^/} && $input =~ m{^/};
die "[ERROR] Output paths must be absolute\n"
    unless $assignments_out =~ m{^/} && $clusters_out =~ m{^/};
die "[ERROR] Output files must be different\n"
    if $assignments_out eq $clusters_out;

for my $path ($assignments_out, $clusters_out) {
    make_path(dirname($path)) unless -d dirname($path);
}

my %conf = read_config($config);
my $cluster_window = required_config(\%conf, 'SA_SPLIT_CLUSTER_WINDOW');
my $min_support_reads = required_config(\%conf, 'SA_SPLIT_MIN_SUPPORT_READS');
validate_nonnegative_integer('SA_SPLIT_CLUSTER_WINDOW', $cluster_window);
validate_positive_integer('SA_SPLIT_MIN_SUPPORT_READS', $min_support_reads);

my ($events_ref, $input_count) = read_events($input);
my @clusters = cluster_events($events_ref, $cluster_window);

for my $cluster (@clusters) {
    $cluster->{Cluster_Status} =
        $cluster->{Support_Reads} >= $min_support_reads ? 'Pass' : 'LowSupport';
}

write_assignments($assignments_out, \@clusters);
write_clusters($clusters_out, \@clusters);

print STDERR "[INFO] Input deletion events: $input_count\n";
print STDERR "[INFO] Unique deletion events: ", scalar(@$events_ref), "\n";
print STDERR "[INFO] Clusters: ", scalar(@clusters), "\n";
print STDERR "[INFO] PASS clusters: ",
    scalar(grep { $_->{Cluster_Status} eq 'Pass' } @clusters), "\n";
print STDERR "[INFO] LowSupport clusters: ",
    scalar(grep { $_->{Cluster_Status} eq 'LowSupport' } @clusters), "\n";
print STDERR "[INFO] Assignments: $assignments_out\n";
print STDERR "[INFO] Clusters: $clusters_out\n";
exit 0;

# ------------------------------------------------------------------
# Input reader
# ------------------------------------------------------------------

sub read_events {
    my ($file) = @_;
    open my $fh, '<', $file or die "[ERROR] Cannot read $file: $!\n";
    my $header_line = <$fh>;
    die "[ERROR] Empty event file: $file\n" unless defined $header_line;
    chomp $header_line;
    $header_line =~ s/\r$//;

    my @header = split /\t/, $header_line, -1;
    my %index = map { $header[$_] => $_ } 0 .. $#header;
    for my $column (qw/Event_ID SampleID Read_Name Chrom Start End Size_bp Strand/) {
        die "[ERROR] Missing column '$column' in $file\n"
            unless exists $index{$column};
    }

    my (@events, %seen_event_id);
    my $input_count = 0;

    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        next if $line eq '';
        $input_count++;

        my @values = split /\t/, $line, -1;
        my %event = map {
            $_ => defined $values[$index{$_}] ? $values[$index{$_}] : 'NA'
        } @header;

        unless ($event{Start} =~ /^\d+$/ && $event{End} =~ /^\d+$/
            && $event{Start} >= 1 && $event{End} >= $event{Start}) {
            die "[ERROR] Invalid event coordinates for $event{Event_ID}\n";
        }

        die "[ERROR] Duplicate Event_ID '$event{Event_ID}' in $file\n"
            if $seen_event_id{$event{Event_ID}}++;
        push @events, \%event;
    }
    close $fh;

    @events = sort {
        chrom_key($a->{Chrom}) cmp chrom_key($b->{Chrom})
            || $a->{Start} <=> $b->{Start}
            || $a->{End} <=> $b->{End}
            || $a->{SampleID} cmp $b->{SampleID}
            || $a->{Read_Name} cmp $b->{Read_Name}
    } @events;

    return (\@events, $input_count);
}

# ------------------------------------------------------------------
# Clustering
# ------------------------------------------------------------------

sub cluster_events {
    my ($events_ref, $window) = @_;
    return () unless @$events_ref;

    my (@clusters, @current);
    for my $event (@$events_ref) {
        if (!@current) {
            @current = ($event);
            next;
        }

        my $current_chrom = $current[0]{Chrom};
        if ($event->{Chrom} ne $current_chrom) {
            push @clusters, build_cluster(\@current);
            @current = ($event);
            next;
        }

        my $median_start = median_int(map { $_->{Start} } @current);
        my $median_end   = median_int(map { $_->{End} } @current);

        if (abs($event->{Start} - $median_start) <= $window
            && abs($event->{End} - $median_end) <= $window) {
            push @current, $event;
        }
        else {
            push @clusters, build_cluster(\@current);
            @current = ($event);
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
        SampleID => $records_ref->[0]{SampleID},
        Chrom => $records_ref->[0]{Chrom},
        Start => median_int(@starts),
        End => median_int(@ends),
        Support_Reads => scalar(keys %read_names),
        Support_Records => scalar(@$records_ref),
        Start_Min => min_num(@starts),
        Start_Max => max_num(@starts),
        End_Min => min_num(@ends),
        End_Max => max_num(@ends),
        Mean_Size_bp => mean_num(@sizes),
        Size_Min => min_num(@sizes),
        Size_Max => max_num(@sizes),
        Primary_Strand_Count => count_strands($records_ref, 'Primary_Strand'),
        SA_Strand_Count => count_strands($records_ref, 'SA_Strand'),
        records => $records_ref,
    };
}

sub count_strands {
    my ($records_ref, $field) = @_;
    my %counts;
    for my $record (@$records_ref) {
        next unless defined $record->{$field};
        $counts{ $record->{$field} }++;
    }
    return join(',', map { $_ . ':' . $counts{$_} } sort keys %counts) || 'NA';
}

# ------------------------------------------------------------------
# Output writers
# ------------------------------------------------------------------

sub write_assignments {
    my ($file, $clusters_ref) = @_;
    open my $fh, '>', $file or die "[ERROR] Cannot write $file: $!\n";
    print {$fh} "Event_ID\tCluster_ID\tCluster_Status\n";

    for my $cluster (@$clusters_ref) {
        for my $event (@{ $cluster->{records} }) {
            print {$fh} join("\t",
                $event->{Event_ID},
                $cluster->{Cluster_ID},
                $cluster->{Cluster_Status},
            ), "\n";
        }
    }
    close $fh;
}

sub write_clusters {
    my ($file, $clusters_ref) = @_;
    open my $fh, '>', $file or die "[ERROR] Cannot write $file: $!\n";

    my @header = qw(
        SampleID Cluster_ID Chrom Start End Size_bp Support_Reads Support_Records
        Start_Range End_Range Mean_Size_bp Size_Range_bp
        Primary_Strand_Count SA_Strand_Count Cluster_Status
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
            $cluster->{Primary_Strand_Count},
            $cluster->{SA_Strand_Count},
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

sub usage {
    return "Usage: perl ClusterSASplitReads.pl --Config CONF --Input split_deletion_events.tsv --AssignmentsOut split_cluster_assignments.tsv --ClustersOut split_clusters.tsv\n";
}
