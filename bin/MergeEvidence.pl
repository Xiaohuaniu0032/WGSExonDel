#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions Configure);
Configure(qw(no_ignore_case no_auto_abbrev));
use File::Basename qw(dirname basename);
use File::Path qw(make_path);

my %opt = (min_overlap => 0.90);
GetOptions(
    'Sample=s'                  => \$opt{sample},
    'Depth=s'                   => \$opt{depth},
    'Split=s'                   => \$opt{split},
    'Discordant=s'              => \$opt{discordant},
    'DepthCandidatesOut=s'      => \$opt{depth_out},
    'SplitCandidatesOut=s'      => \$opt{split_out},
    'DiscordantCandidatesOut=s' => \$opt{discordant_out},
    'Out=s'                     => \$opt{out},
    'MinOverlap=f'              => \$opt{min_overlap},
    'Help'                      => \$opt{help},
) or die usage();

if ($opt{help}) {
    print usage();
    exit 0;
}

for my $key (qw/sample depth split discordant depth_out split_out discordant_out out/) {
    die "[ERROR] Missing required option: $key\n" . usage()
        unless defined $opt{$key} && $opt{$key} ne '';
}

die "[ERROR] MinOverlap must be greater than 0 and no greater than 1\n"
    unless $opt{min_overlap} > 0 && $opt{min_overlap} <= 1;
die "[ERROR] Invalid Sample ID\n" unless $opt{sample} =~ /^[A-Za-z0-9_.-]+$/;

for my $key (qw/depth split discordant/) {
    die "[ERROR] Input file not found: $opt{$key}\n" unless -s $opt{$key};
}

my %output_paths;
for my $key (qw/depth_out split_out discordant_out out/) {
    my $path = $opt{$key};
    die "[ERROR] Output path must be absolute: $path\n" unless $path =~ m{^/};
    die "[ERROR] Duplicate output path: $path\n" if $output_paths{$path}++;
    make_path(dirname($path)) unless -d dirname($path);
}

log_msg("Preparing evidence candidates for $opt{sample}");

my @depth_records = read_depth_records($opt{depth}, $opt{sample});
my @split_records = read_split_records($opt{split}, $opt{sample});
my @discordant_records = read_discordant_records($opt{discordant}, $opt{sample});

@depth_records = sort { interval_sort($a, $b) } @depth_records;
@split_records = sort { interval_sort($a, $b) } @split_records;
@discordant_records = sort { interval_sort($a, $b) } @discordant_records;

assign_candidate_ids(\@depth_records, 'DEPTH');
assign_candidate_ids(\@split_records, 'SPLIT');
assign_candidate_ids(\@discordant_records, 'DISCORDANT');

write_candidate_table($opt{depth_out}, 'depth', \@depth_records);
write_candidate_table($opt{split_out}, 'split', \@split_records);
write_candidate_table($opt{discordant_out}, 'discordant', \@discordant_records);

my %depth_by_chrom = index_by_chromosome(@depth_records);
my %discordant_by_chrom = index_by_chromosome(@discordant_records);
my @merged_candidates;

for my $split (@split_records) {
    my ($chrom, $start, $end) = @{$split}{qw/Chrom Start End/};
    my $candidate_length = interval_length($start, $end);
    die "[ERROR] Invalid split candidate interval: $chrom:$start-$end\n"
        unless $candidate_length > 0;

    # Depth support is based on the union of overlapping CUSUM DEL intervals.
    my @depth_hits = grep {
        overlap_length($start, $end, $_->{Start}, $_->{End}) > 0
    } records_for_chromosome(\%depth_by_chrom, $chrom);

    my ($depth_covered_bases, $depth_coverage_fraction) =
        covered_bases_by_union(\@depth_hits, $start, $end);
    my $depth_support = $depth_coverage_fraction >= $opt{min_overlap} ? 1 : 0;

    # Discordant support retains the existing interval-overlap rule.
    my @discordant_overlaps = grep {
        overlap_length($start, $end, $_->{Start}, $_->{End}) > 0
    } records_for_chromosome(\%discordant_by_chrom, $chrom);

    my (@discordant_hits, $best_discordant_fraction);
    $best_discordant_fraction = 0;
    for my $discordant (@discordant_overlaps) {
        my $overlap_bases = overlap_length(
            $start, $end, $discordant->{Start}, $discordant->{End}
        );
        my $overlap_fraction = $overlap_bases / $candidate_length;
        $best_discordant_fraction = $overlap_fraction
            if $overlap_fraction > $best_discordant_fraction;
        push @discordant_hits, $discordant
            if $overlap_fraction >= $opt{min_overlap};
    }

    my $discordant_support = @discordant_hits ? 1 : 0;

    # Keep the established Evidence_Level rule exactly as before.
    my $evidence_level = $depth_support
        ? ($discordant_support ? 'High' : 'Moderate')
        : 'Low';

    my @evidence_types = ('Split');
    push @evidence_types, 'Depth' if $depth_support;
    push @evidence_types, 'Discordant' if $discordant_support;

    my $discordant_pairs = sum_unique_source_values(
        \@discordant_hits,
        'Support_Pairs'
    );
    my $median_insert_size = join_unique_values(
        map { metadata_value($_, 'Median_Insert_Size') } @discordant_hits
    );

    push @merged_candidates, {
        SampleID                  => $opt{sample},
        Candidate_ID              => sprintf('CANDIDATE_%06d', scalar(@merged_candidates) + 1),
        Split_Cluster_ID          => $split->{Source_ID},
        Chrom                     => $chrom,
        Start                     => $start,
        End                       => $end,
        Size_bp                   => $candidate_length,
        Split_Reads               => metadata_value($split, 'Support_Reads'),
        Depth_Support             => $depth_support ? 'Yes' : 'No',
        Depth_Covered_Bases       => $depth_covered_bases,
        Depth_Coverage_Fraction   => sprintf('%.4f', $depth_coverage_fraction),
        Depth_Candidate_Count     => scalar(@depth_hits),
        Discordant_Support        => $discordant_support ? 'Yes' : 'No',
        Discordant_Pairs          => $discordant_pairs,
        Discordant_Overlap_Fraction => sprintf('%.4f', $best_discordant_fraction),
        Discordant_Cluster_Count  => scalar(@discordant_hits),
        Median_Insert_Size        => $median_insert_size,
        Evidence_Types            => join(',', @evidence_types),
        Evidence_Level            => $evidence_level,
    };
}

write_merged_table($opt{out}, \@merged_candidates);

my %level_counts;
$level_counts{$_->{Evidence_Level}}++ for @merged_candidates;

log_msg('Depth candidates: ' . scalar(@depth_records));
log_msg('PASS Split candidates: ' . scalar(@split_records));
log_msg('PASS Discordant candidates: ' . scalar(@discordant_records));
log_msg('Merged Split-centered candidates: ' . scalar(@merged_candidates));
log_msg('Evidence_Level: High=' . ($level_counts{High} || 0)
    . ', Moderate=' . ($level_counts{Moderate} || 0)
    . ', Low=' . ($level_counts{Low} || 0));
exit 0;

# ------------------------------------------------------------------
# Input readers
# ------------------------------------------------------------------

sub read_depth_records {
    my ($file, $sample) = @_;
    my ($header, $rows) = read_tsv($file);
    my %index = header_index($header);

    my %column = (
        sample => required_column(\%index, 'SampleID', 'Sample'),
        chrom  => required_column(\%index, 'Chrom', 'CHROM'),
        start  => required_column(\%index, 'Start', 'START'),
        end    => required_column(\%index, 'End', 'END'),
        type   => required_column(\%index, 'Type', 'TYPE'),
    );

    my @records;
    for my $row (@$rows) {
        next unless field_value($row, $column{sample}) eq $sample;
        next unless field_value($row, $column{type}) =~ /DEL/i;

        my ($start, $end) = valid_interval(
            field_value($row, $column{start}),
            field_value($row, $column{end}),
            $file,
        );

        push @records, {
            SampleID => $sample,
            Chrom    => normalize_chromosome(field_value($row, $column{chrom})),
            Start    => $start,
            End      => $end,
            Source_ID => 'DEPTH_EVENT_' . (scalar(@records) + 1),
            Meta     => $row,
            Header   => $header,
        };
    }
    return @records;
}

sub read_split_records {
    my ($file, $sample) = @_;
    my ($header, $rows) = read_tsv($file);
    my %index = header_index($header);

    my %column = (
        sample  => required_column(\%index, 'SampleID'),
        id      => required_column(\%index, 'Cluster_ID'),
        chrom   => required_column(\%index, 'Chrom'),
        start   => required_column(\%index, 'Start'),
        end     => required_column(\%index, 'End'),
        status  => required_column(\%index, 'Cluster_Status'),
    );

    my @records;
    for my $row (@$rows) {
        next unless field_value($row, $column{sample}) eq $sample;
        next unless field_value($row, $column{status}) =~ /^Pass$/i;

        my ($start, $end) = valid_interval(
            field_value($row, $column{start}),
            field_value($row, $column{end}),
            $file,
        );

        push @records, {
            SampleID => $sample,
            Chrom    => normalize_chromosome(field_value($row, $column{chrom})),
            Start    => $start,
            End      => $end,
            Source_ID => field_value($row, $column{id}),
            Meta     => $row,
            Header   => $header,
        };
    }
    return @records;
}

sub read_discordant_records {
    my ($file, $sample) = @_;
    my ($header, $rows) = read_tsv($file);
    my %index = header_index($header);

    my %column = (
        sample  => required_column(\%index, 'SampleID'),
        id      => required_column(\%index, 'Cluster_ID'),
        chrom   => required_column(\%index, 'Chrom'),
        start   => required_column(\%index, 'Start'),
        end     => required_column(\%index, 'End'),
        status  => required_column(\%index, 'Cluster_Status'),
    );

    my @records;
    for my $row (@$rows) {
        next unless field_value($row, $column{sample}) eq $sample;
        next unless field_value($row, $column{status}) =~ /^Pass$/i;

        my ($start, $end) = valid_interval(
            field_value($row, $column{start}),
            field_value($row, $column{end}),
            $file,
        );

        push @records, {
            SampleID => $sample,
            Chrom    => normalize_chromosome(field_value($row, $column{chrom})),
            Start    => $start,
            End      => $end,
            Source_ID => field_value($row, $column{id}),
            Meta     => $row,
            Header   => $header,
        };
    }
    return @records;
}

# ------------------------------------------------------------------
# Candidate tables
# ------------------------------------------------------------------

sub assign_candidate_ids {
    my ($records_ref, $prefix) = @_;
    for my $i (0 .. $#$records_ref) {
        $records_ref->[$i]{Candidate_ID} = sprintf('%s_%06d', $prefix, $i + 1);
    }
}

sub write_candidate_table {
    my ($file, $type, $records_ref) = @_;
    my @header;

    if ($type eq 'depth') {
        @header = qw(
            SampleID Candidate_ID Chrom Start End Size_bp N_Bins
            Baseline_Depth Mean_Depth Mean_Depth_Ratio Min_Depth_Ratio
            Max_CUSUM Depth_File
        );
    }
    elsif ($type eq 'split') {
        @header = qw(
            SampleID Candidate_ID Source_Cluster_ID Chrom Start End
            Size_bp Support_Reads
        );
    }
    else {
        @header = qw(
            SampleID Candidate_ID Source_Cluster_ID Chrom Start End Size_bp
            Support_Pairs Median_Insert_Size Mean_Insert_Size Insert_Size_Range
            Pair_Orientation_Count
        );
    }

    open my $fh, '>', $file or die "[ERROR] Cannot write $file: $!\n";
    print {$fh} join("\t", @header), "\n";

    for my $record (@$records_ref) {
        my %value = (
            SampleID              => $record->{SampleID},
            Candidate_ID          => $record->{Candidate_ID},
            Source_Cluster_ID     => $record->{Source_ID},
            Chrom                 => $record->{Chrom},
            Start                 => $record->{Start},
            End                   => $record->{End},
            Size_bp               => interval_length($record->{Start}, $record->{End}),
            Support_Reads         => metadata_value($record, 'Support_Reads'),
            N_Bins                => metadata_value($record, 'N_Bins'),
            Baseline_Depth        => metadata_value($record, 'Baseline'),
            Mean_Depth            => metadata_value($record, 'Mean_Depth'),
            Mean_Depth_Ratio      => metadata_value($record, 'Mean_Ratio'),
            Min_Depth_Ratio       => metadata_value($record, 'Min_Ratio'),
            Max_CUSUM             => metadata_value($record, 'Max_CUSUM'),
            Depth_File            => basename(metadata_value($record, 'Depth_File')),
            Support_Pairs         => metadata_value($record, 'Support_Pairs'),
            Median_Insert_Size    => metadata_value($record, 'Median_Insert_Size'),
            Mean_Insert_Size      => metadata_value($record, 'Mean_Insert_Size'),
            Insert_Size_Range     => metadata_value($record, 'Insert_Size_Range'),
            Pair_Orientation_Count => metadata_value($record, 'Pair_Orientation_Count'),
        );

        print {$fh} join("\t", map { clean_value($value{$_}) } @header), "\n";
    }
    close $fh;
}

sub write_merged_table {
    my ($file, $records_ref) = @_;
    my @header = qw(
        SampleID Candidate_ID Split_Cluster_ID Chrom Start End Size_bp Split_Reads
        Depth_Support Depth_Covered_Bases Depth_Coverage_Fraction Depth_Candidate_Count
        Discordant_Support Discordant_Pairs Discordant_Overlap_Fraction
        Discordant_Cluster_Count Median_Insert_Size Evidence_Types Evidence_Level
    );

    open my $fh, '>', $file or die "[ERROR] Cannot write $file: $!\n";
    print {$fh} join("\t", @header), "\n";
    for my $record (@$records_ref) {
        print {$fh} join("\t", map { clean_value($record->{$_}) } @header), "\n";
    }
    close $fh;
}

# ------------------------------------------------------------------
# TSV and interval helpers
# ------------------------------------------------------------------

sub read_tsv {
    my ($file) = @_;
    open my $fh, '<', $file or die "[ERROR] Cannot read $file: $!\n";
    my $line = <$fh>;
    die "[ERROR] Empty TSV: $file\n" unless defined $line;
    chomp $line;
    $line =~ s/\r$//;

    my @header = split /\t/, $line, -1;
    my @rows;
    while (my $row_line = <$fh>) {
        chomp $row_line;
        $row_line =~ s/\r$//;
        next if $row_line eq '';
        my @values = split /\t/, $row_line, -1;
        die "[ERROR] Column count mismatch in $file\n" unless @values == @header;
        push @rows, \@values;
    }
    close $fh;
    return (\@header, \@rows);
}

sub header_index {
    my ($header_ref) = @_;
    my %index;
    for my $i (0 .. $#$header_ref) {
        my $name = $header_ref->[$i];
        die "[ERROR] Duplicate header '$name'\n" if exists $index{$name};
        $index{$name} = $i;
    }
    return %index;
}

sub required_column {
    my ($index_ref, @names) = @_;
    for my $name (@names) {
        return $index_ref->{$name} if exists $index_ref->{$name};
    }
    die "[ERROR] Missing required column: " . join('/', @names) . "\n";
}

sub field_value {
    my ($row_ref, $index) = @_;
    return defined $row_ref->[$index] ? $row_ref->[$index] : '';
}

sub metadata_value {
    my ($record, $name) = @_;
    return 'NA' unless $record->{Meta} && $record->{Header};

    my %index = map { $record->{Header}[$_] => $_ } 0 .. $#{ $record->{Header} };
    return 'NA' unless exists $index{$name};

    my $value = $record->{Meta}[ $index{$name} ];
    return defined $value && $value ne '' ? $value : 'NA';
}

sub valid_interval {
    my ($start, $end, $file) = @_;
    die "[ERROR] Invalid interval in $file: $start-$end\n"
        unless $start =~ /^\d+$/ && $end =~ /^\d+$/ && $start > 0 && $end >= $start;
    return ($start, $end);
}

sub normalize_chromosome {
    my ($chrom) = @_;
    $chrom =~ s/^chr//i;
    return 'chr' . $chrom;
}

sub index_by_chromosome {
    my @records = @_;
    my %index;
    push @{ $index{$_->{Chrom}} }, $_ for @records;
    return %index;
}

sub records_for_chromosome {
    my ($index_ref, $chrom) = @_;
    return @{ $index_ref->{$chrom} || [] };
}

sub covered_bases_by_union {
    my ($hits_ref, $candidate_start, $candidate_end) = @_;
    return (0, 0) unless @$hits_ref;

    my @intervals = sort { $a->[0] <=> $b->[0] }
        map {
            [
                $_->{Start} > $candidate_start ? $_->{Start} : $candidate_start,
                $_->{End} < $candidate_end ? $_->{End} : $candidate_end,
            ]
        }
        grep {
            $_->{Start} <= $candidate_end && $_->{End} >= $candidate_start
        } @$hits_ref;

    return (0, 0) unless @intervals;

    my ($current_start, $current_end) = @{ $intervals[0] };
    my $covered_bases = 0;

    for my $i (1 .. $#intervals) {
        if ($intervals[$i][0] <= $current_end + 1) {
            $current_end = $intervals[$i][1]
                if $intervals[$i][1] > $current_end;
        }
        else {
            $covered_bases += $current_end - $current_start + 1;
            ($current_start, $current_end) = @{ $intervals[$i] };
        }
    }

    $covered_bases += $current_end - $current_start + 1;
    my $candidate_length = interval_length($candidate_start, $candidate_end);
    my $fraction = $candidate_length ? $covered_bases / $candidate_length : 0;
    return ($covered_bases, $fraction);
}

sub sum_unique_source_values {
    my ($records_ref, $field) = @_;
    my (%seen, $sum);
    $sum = 0;

    for my $record (@$records_ref) {
        next if $seen{ $record->{Source_ID} }++;
        my $value = metadata_value($record, $field);
        $sum += $value if $value =~ /^\d+(?:\.\d+)?$/;
    }
    return $sum;
}

sub join_unique_values {
    my @values = @_;
    my %seen;
    my @unique = grep {
        defined $_ && $_ ne 'NA' && $_ ne '' && !$seen{$_}++
    } @values;
    return @unique ? join(',', @unique) : 'NA';
}

sub interval_length {
    my ($start, $end) = @_;
    return $end >= $start ? $end - $start + 1 : 0;
}

sub overlap_length {
    my ($start1, $end1, $start2, $end2) = @_;
    my $start = $start1 > $start2 ? $start1 : $start2;
    my $end = $end1 < $end2 ? $end1 : $end2;
    return $end >= $start ? $end - $start + 1 : 0;
}

sub interval_sort {
    return chrom_rank($a->{Chrom}) <=> chrom_rank($b->{Chrom})
        || $a->{Start} <=> $b->{Start}
        || $a->{End} <=> $b->{End};
}

sub chrom_rank {
    my ($chrom) = @_;
    $chrom =~ s/^chr//i;
    return $1 if $chrom =~ /^(\d+)$/;
    return 23 if uc($chrom) eq 'X';
    return 24 if uc($chrom) eq 'Y';
    return 25 if uc($chrom) eq 'M' || uc($chrom) eq 'MT';
    return 26;
}

sub clean_value {
    my ($value) = @_;
    return 'NA' unless defined $value && $value ne '';
    $value =~ s/[\t\r\n]/ /g;
    return $value;
}

sub log_msg {
    print STDERR "[INFO] $_[0]\n";
}

sub usage {
    return "Usage: perl MergeEvidence.pl --Sample ID --Depth cusum_events.tsv --Split split_clusters.tsv --Discordant discordant_clusters.tsv --DepthCandidatesOut depth_candidates.tsv --SplitCandidatesOut split_candidates.tsv --DiscordantCandidatesOut discordant_candidates.tsv --Out merged_candidates.tsv [--MinOverlap 0.90]\n";
}
