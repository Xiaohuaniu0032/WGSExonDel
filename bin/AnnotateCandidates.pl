#!/usr/bin/env perl
# -*- coding: utf-8 -*-

use strict;
use warnings;
use Getopt::Long qw(GetOptions Configure);
Configure(qw(no_ignore_case no_auto_abbrev));
use File::Basename qw(dirname);
use File::Path qw(make_path);
use Cwd qw(abs_path);

my %opt;
GetOptions(
    'Config=s'         => \$opt{config},
    'GeneAnnotation=s' => \$opt{gene_annotation},
    'ExonAnnotation=s' => \$opt{exon_annotation},
    'Input=s'          => \$opt{input},
    'Out=s'            => \$opt{out},
    'Help'             => \$opt{help},
) or die usage();

if ($opt{help}) {
    print usage();
    exit 0;
}

for my $key (qw/config gene_annotation exon_annotation input out/) {
    die usage() unless defined $opt{$key} && $opt{$key} ne '';
}

$opt{config} = resolve_existing_file($opt{config}, 'Config');
$opt{gene_annotation} = resolve_existing_file($opt{gene_annotation}, 'GeneAnnotation');
$opt{exon_annotation} = resolve_existing_file($opt{exon_annotation}, 'ExonAnnotation');
$opt{input} = resolve_existing_file($opt{input}, 'Input');

die "[ERROR] Output path must be absolute: $opt{out}\n"
    unless $opt{out} =~ m{^/};

my $out_dir = dirname($opt{out});
make_path($out_dir) unless -d $out_dir;
die "[ERROR] Cannot create output directory: $out_dir\n" unless -d $out_dir;

my %conf = read_config($opt{config});
my $min_exon_fraction = required_number(\%conf, 'MIN_EXON_OVERLAP_FRACTION', 0, 1);
my $flank_distance = required_integer(\%conf, 'TARGET_REGION_FLANK', 0);
my $gene_list_file = $conf{HCM_CORE_GENE_LIST};
$gene_list_file = resolve_existing_file($gene_list_file, 'HCM_CORE_GENE_LIST');

my $gene_classification = read_gene_classification($gene_list_file);
my $genes_by_chrom = read_interval_annotation(
    $opt{gene_annotation},
    [qw/Gene Transcript Chrom Start End Strand/],
    'gene'
);
my $exons_by_chrom = read_interval_annotation(
    $opt{exon_annotation},
    [qw/Gene Transcript Exon Chrom Start End Strand/],
    'exon'
);

my ($candidate_header, $candidate_rows) = read_tsv($opt{input});
my %candidate_idx = header_index($candidate_header, $opt{input});
my @required_columns = qw/
    SampleID Candidate_ID Chrom Start End Split_Reads
    Depth_Support Depth_Coverage_Fraction
    Discordant_Support Discordant_Pairs Discordant_Overlap_Fraction
    Median_Insert_Size Evidence_Level
/;
for my $column (@required_columns) {
    die "[ERROR] Required column '$column' is missing from $opt{input}\n"
        unless exists $candidate_idx{$column};
}

my @output_header = qw/
    SampleID Candidate_ID Gene Gene_Classification Transcript Region_Type
    Chrom Start End Size_bp Affected_Exons Fully_Covered_Exons
    Split_Reads Depth_Support Depth_Coverage_Fraction
    Discordant_Support Discordant_Pairs Discordant_Overlap_Fraction Median_Insert_Size Evidence_Level
/;

print "[INFO] Candidate annotation started\n";
print "[INFO] Config: $opt{config}\n";
print "[INFO] Input: $opt{input}\n";
print "[INFO] Gene annotation: $opt{gene_annotation}\n";
print "[INFO] Exon annotation: $opt{exon_annotation}\n";
print "[INFO] HCM gene list: $gene_list_file\n";
print "[INFO] Minimum exon overlap fraction: $min_exon_fraction\n";
print "[INFO] Flanking distance: $flank_distance bp\n";
print "[INFO] Candidate rows: ", scalar(@$candidate_rows), "\n";
print "[INFO] Output: $opt{out}\n";

open my $out_fh, '>', $opt{out}
    or die "[ERROR] Cannot write output file $opt{out}: $!\n";
print {$out_fh} join("\t", @output_header), "\n";

my %region_count;
my %seen_candidate;
my $annotated_count = 0;

for my $row_number (0 .. $#$candidate_rows) {
    my $row = $candidate_rows->[$row_number];
    my %candidate;
    @candidate{@$candidate_header} = @$row;

    my $candidate_id = field_value(\%candidate, 'Candidate_ID');
    die "[ERROR] Empty Candidate_ID at data row " . ($row_number + 1) . "\n"
        if is_missing($candidate_id);
    die "[ERROR] Duplicate Candidate_ID '$candidate_id' in $opt{input}\n"
        if $seen_candidate{$candidate_id}++;

    my $chrom = field_value(\%candidate, 'Chrom');
    my $start = field_value(\%candidate, 'Start');
    my $end = field_value(\%candidate, 'End');
    die "[ERROR] Invalid candidate coordinates at data row " . ($row_number + 1) . ": $chrom:$start-$end\n"
        unless !is_missing($chrom) && is_positive_integer($start) && is_positive_integer($end);
    ($start, $end) = ($end, $start) if $start > $end;

    my $annotation = annotate_candidate(
        chrom => $chrom,
        start => $start,
        end => $end,
        genes_by_chrom => $genes_by_chrom,
        exons_by_chrom => $exons_by_chrom,
        gene_classification => $gene_classification,
        min_exon_fraction => $min_exon_fraction,
        flank_distance => $flank_distance,
    );

    $region_count{$annotation->{region_type}}++;
    $annotated_count++ if $annotation->{region_type} ne 'Intergenic';

    my $depth_support = normalize_support(field_value(\%candidate, 'Depth_Support'));
    my $discordant_support = normalize_support(field_value(\%candidate, 'Discordant_Support'));
    my $depth_fraction = format_fraction(field_value(\%candidate, 'Depth_Coverage_Fraction'));
    my $discordant_pairs = normalize_pair_count(
        field_value(\%candidate, 'Discordant_Pairs'),
        $discordant_support
    );
    my $discordant_fraction = format_fraction(field_value(\%candidate, 'Discordant_Overlap_Fraction'));
    my $median_insert_size = $discordant_support eq 'No'
        ? 'NA'
        : normalize_numeric_or_na(field_value(\%candidate, 'Median_Insert_Size'));

    my @output = (
        field_value(\%candidate, 'SampleID'),
        $candidate_id,
        $annotation->{gene},
        $annotation->{gene_classification},
        $annotation->{transcript},
        $annotation->{region_type},
        normalize_chrom($chrom),
        $start,
        $end,
        $end - $start + 1,
        $annotation->{affected_exons},
        $annotation->{fully_covered_exons},
        normalize_count(field_value(\%candidate, 'Split_Reads')),
        $depth_support,
        $depth_fraction,
        $discordant_support,
        $discordant_pairs,
        $discordant_fraction,
        $median_insert_size,
        field_value(\%candidate, 'Evidence_Level'),
    );

    @output = map { clean_tsv_value($_) } @output;
    print {$out_fh} join("\t", @output), "\n";
}

close $out_fh or die "[ERROR] Cannot close output file $opt{out}: $!\n";

print "[INFO] Annotated within/near gene regions: $annotated_count\n";
for my $region (qw/Exonic Intronic Flanking Intergenic/) {
    print "[INFO] $region candidates: ", ($region_count{$region} || 0), "\n";
}
print "[INFO] Candidate annotation finished\n";

exit 0;

# ------------------------------------------------------------------
# Candidate annotation
# ------------------------------------------------------------------

sub annotate_candidate {
    my (%args) = @_;
    my ($chrom, $start, $end) = @args{qw/chrom start end/};
    my $gene_records = find_overlapping_records($args{genes_by_chrom}, $chrom, $start, $end);
    my $exon_records = find_overlapping_records($args{exons_by_chrom}, $chrom, $start, $end);

    my $region_type;
    my @selected_genes;

    if (@$gene_records) {
        @selected_genes = @$gene_records;
        $region_type = @$exon_records ? 'Exonic' : 'Intronic';
    }
    else {
        my ($nearest_genes, $distance) = find_nearest_genes(
            $args{genes_by_chrom}, $chrom, $start, $end
        );
        if (@$nearest_genes && $distance <= $args{flank_distance}) {
            @selected_genes = @$nearest_genes;
            $region_type = 'Flanking';
        }
        else {
            @selected_genes = ();
            $region_type = 'Intergenic';
        }
    }

    my @gene_names = sort_unique(map { $_->{Gene} } @selected_genes);
    my (@transcripts, @classifications);
    for my $gene (@gene_names) {
        my @gene_transcripts = sort_unique(
            map { $_->{Transcript} }
            grep { $_->{Gene} eq $gene } @selected_genes
        );
        push @transcripts, @gene_transcripts ? join(';', @gene_transcripts) : 'NA';
        push @classifications, $args{gene_classification}{$gene} // 'NA';
    }

    my (@affected, @fully_covered);
    my %affected_keys;
    for my $exon (@$exon_records) {
        my $overlap = overlap_length($start, $end, $exon->{Start}, $exon->{End});
        next if $overlap <= 0;
        my $exon_length = $exon->{End} - $exon->{Start} + 1;
        my $fraction = $overlap / $exon_length;
        next if $fraction < $args{min_exon_fraction};

        my $key = join('|', @$exon{qw/Gene Transcript Exon/});
        next if $affected_keys{$key}++;
        push @affected, $exon;

        if ($start <= $exon->{Start} && $end >= $exon->{End}) {
            push @fully_covered, $exon;
        }
    }

    my $exon_label = sub {
        my ($records) = @_;
        my @labels = sort_exon_labels(map { 'EX' . $_->{Exon} } @$records);
        return @labels ? join(',', @labels) : 'NA';
    };

    return {
        gene => @gene_names ? join(',', @gene_names) : 'NA',
        gene_classification => @classifications ? join(',', @classifications) : 'NA',
        transcript => @transcripts ? join(',', @transcripts) : 'NA',
        region_type => $region_type,
        affected_exons => $exon_label->(\@affected),
        fully_covered_exons => $exon_label->(\@fully_covered),
    };
}

sub find_overlapping_records {
    my ($index, $chrom, $start, $end) = @_;
    my @records;
    my %seen;
    for my $key (chrom_keys($chrom)) {
        for my $record (@{ $index->{$key} || [] }) {
            last if $record->{Start} > $end;
            next if $record->{End} < $start;
            my $id = join('|', @$record{qw/Gene Transcript Chrom Start End/});
            push @records, $record unless $seen{$id}++;
        }
    }
    return \@records;
}

sub find_nearest_genes {
    my ($index, $chrom, $start, $end) = @_;
    my $best_distance;
    my @nearest;
    my %seen;
    for my $key (chrom_keys($chrom)) {
        for my $record (@{ $index->{$key} || [] }) {
            my $distance = $record->{Start} > $end
                ? $record->{Start} - $end
                : $start > $record->{End}
                    ? $start - $record->{End}
                    : 0;
            my $id = join('|', @$record{qw/Gene Transcript Chrom Start End/});
            next if $seen{$id}++;
            if (!defined $best_distance || $distance < $best_distance) {
                $best_distance = $distance;
                @nearest = ($record);
            }
            elsif ($distance == $best_distance) {
                push @nearest, $record;
            }
        }
    }
    return (\@nearest, defined $best_distance ? $best_distance : 1_000_000_000);
}

# ------------------------------------------------------------------
# Annotation table readers
# ------------------------------------------------------------------

sub read_interval_annotation {
    my ($file, $required, $type) = @_;
    my ($header, $rows) = read_tsv($file);
    my %idx = header_index($header, $file);
    for my $column (@$required) {
        die "[ERROR] Required column '$column' is missing from $file\n"
            unless exists $idx{$column};
    }

    my %by_chrom;
    my $line_number = 1;
    for my $row (@$rows) {
        $line_number++;
        my %record;
        @record{@$header} = @$row;
        my ($gene, $transcript, $chrom, $start, $end) = @record{qw/Gene Transcript Chrom Start End/};
        die "[ERROR] Invalid $type annotation at $file line $line_number\n"
            if is_missing($gene) || is_missing($transcript) || is_missing($chrom)
                || !is_positive_integer($start) || !is_positive_integer($end) || $start > $end;

        my $normalized_chrom = normalize_chrom($chrom);
        my $entry = {
            Gene => $gene,
            Transcript => $transcript,
            Chrom => $normalized_chrom,
            Start => 0 + $start,
            End => 0 + $end,
            Strand => $record{Strand} // 'NA',
            Exon => $record{Exon} // 'NA',
        };
        if ($type eq 'exon') {
            die "[ERROR] Missing exon number at $file line $line_number\n"
                if is_missing($entry->{Exon});
            $entry->{Exon} =~ s/^EX//i;
        }
        for my $key (chrom_keys($normalized_chrom)) {
            push @{ $by_chrom{$key} }, $entry;
        }
    }

    for my $chrom (keys %by_chrom) {
        @{ $by_chrom{$chrom} } = sort {
            $a->{Start} <=> $b->{Start} || $a->{End} <=> $b->{End}
                || $a->{Gene} cmp $b->{Gene}
        } @{ $by_chrom{$chrom} };
    }

    die "[ERROR] No valid $type annotation records found in $file\n"
        unless keys %by_chrom;
    return \%by_chrom;
}

sub read_gene_classification {
    my ($file) = @_;
    my ($header, $rows) = read_tsv($file);
    my %idx = header_index($header, $file);
    my $gene_index = exists $idx{Gene} ? $idx{Gene} : 0;
    my $class_index = exists $idx{Evidence} ? $idx{Evidence}
        : exists $idx{Classification} ? $idx{Classification} : undef;
    my %classification;
    for my $row (@$rows) {
        my $gene = $row->[$gene_index];
        next if is_missing($gene);
        my $class = defined $class_index ? $row->[$class_index] : 'NA';
        $classification{$gene} = $class unless is_missing($class);
    }
    return \%classification;
}

# ------------------------------------------------------------------
# TSV, config and value helpers
# ------------------------------------------------------------------

sub read_tsv {
    my ($file) = @_;
    open my $fh, '<', $file or die "[ERROR] Cannot open $file: $!\n";
    my $line = <$fh>;
    die "[ERROR] Empty TSV file: $file\n" unless defined $line;
    chomp $line;
    $line =~ s/\r$//;
    my @header = split /\t/, $line, -1;
    my %seen;
    for my $name (@header) {
        die "[ERROR] Empty header field in $file\n" if !defined $name || $name eq '';
        die "[ERROR] Duplicate header field '$name' in $file\n" if $seen{$name}++;
    }
    my @rows;
    my $line_number = 1;
    while (my $row_line = <$fh>) {
        $line_number++;
        chomp $row_line;
        $row_line =~ s/\r$//;
        next if $row_line =~ /^\s*$/;
        my @values = split /\t/, $row_line, -1;
        die "[ERROR] Column count mismatch in $file line $line_number: expected " . scalar(@header)
            . ', observed ' . scalar(@values) . "\n"
            unless @values == @header;
        push @rows, \@values;
    }
    close $fh;
    return (\@header, \@rows);
}

sub header_index {
    my ($header, $file) = @_;
    my %index;
    for my $i (0 .. $#$header) {
        $index{$header->[$i]} = $i;
    }
    return %index;
}

sub read_config {
    my ($file) = @_;
    open my $fh, '<', $file or die "[ERROR] Cannot open config $file: $!\n";
    my %config;
    my $line_number = 0;
    while (my $line = <$fh>) {
        $line_number++;
        chomp $line;
        $line =~ s/\r$//;
        next if $line =~ /^\s*$/ || $line =~ /^\s*#/;
        $line =~ s/\s+#.*$//;
        die "[ERROR] Invalid config line $line_number in $file: $line\n"
            unless $line =~ /^\s*([A-Z][A-Z0-9_]*)\s*=\s*(.*?)\s*$/;
        my ($key, $value) = ($1, $2);
        die "[ERROR] Duplicate config key '$key' in $file\n" if exists $config{$key};
        $value =~ s/^['"]|['"]$//g;
        die "[ERROR] Empty config value for '$key' in $file\n" if $value eq '';
        $config{$key} = $value;
    }
    close $fh;
    return %config;
}

sub resolve_existing_file {
    my ($path, $label) = @_;
    die "[ERROR] $label path must be absolute: " . (defined $path ? $path : 'NA') . "\n"
        unless defined $path && $path =~ m{^/};
    my $resolved = abs_path($path);
    die "[ERROR] $label file not found or empty: $path\n"
        unless defined $resolved && -s $resolved;
    return $resolved;
}

sub required_number {
    my ($config, $key, $min, $max) = @_;
    die "[ERROR] Missing config parameter: $key\n"
        unless exists $config->{$key};
    my $value = $config->{$key};
    die "[ERROR] $key must be numeric: $value\n"
        unless defined $value && $value =~ /^\d+(?:\.\d+)?$/;
    die "[ERROR] $key must be > $min and <= $max\n"
        unless $value > $min && $value <= $max;
    return 0 + $value;
}

sub required_integer {
    my ($config, $key, $min) = @_;
    die "[ERROR] Missing config parameter: $key\n" unless exists $config->{$key};
    my $value = $config->{$key};
    die "[ERROR] $key must be an integer >= $min: $value\n"
        unless defined $value && $value =~ /^\d+$/ && $value >= $min;
    return 0 + $value;
}

sub field_value {
    my ($record, $key) = @_;
    return 'NA' unless exists $record->{$key} && defined $record->{$key} && $record->{$key} ne '';
    return $record->{$key};
}

sub normalize_support {
    my ($value) = @_;
    return 'NA' if is_missing($value);
    return 'Yes' if $value =~ /^(?:1|yes|true|pass|supported)$/i;
    return 'No' if $value =~ /^(?:0|no|false|fail|none|unsupported)$/i;
    return $value;
}

sub normalize_pair_count {
    my ($value, $support) = @_;
    return 0 if $support eq 'No';
    return 0 + $value if defined $value && $value =~ /^\d+$/;
    return 'NA';
}

sub normalize_count {
    my ($value) = @_;
    return 0 + $value if defined $value && $value =~ /^\d+$/;
    return 'NA';
}

sub numeric_or_zero {
    my ($value) = @_;
    return 0 if is_missing($value);
    return 0 + $value if $value =~ /^\d+(?:\.\d+)?$/;
    die "[ERROR] Expected a non-negative numeric value, observed: $value\n";
}

sub format_fraction {
    my ($value) = @_;
    return '0.0000' if is_missing($value);
    die "[ERROR] Invalid fraction value: $value\n"
        unless $value =~ /^\d+(?:\.\d+)?$/ && $value >= 0 && $value <= 1;
    return sprintf('%.4f', $value);
}

sub normalize_numeric_or_na {
    my ($value) = @_;
    return 'NA' if is_missing($value);
    my @values = split /,/, $value, -1;
    die "[ERROR] Invalid numeric value: $value\n"
        unless @values && !grep { $_ !~ /^-?\d+(?:\.\d+)?$/ } @values;
    return join(',', @values);
}

sub normalize_chrom {
    my ($chrom) = @_;
    $chrom =~ s/^\s+|\s+$//g;
    return $chrom;
}

sub chrom_keys {
    my ($chrom) = @_;
    my @keys = ($chrom);
    if ($chrom =~ /^chr(.+)$/i) {
        push @keys, $1;
    }
    else {
        push @keys, 'chr' . $chrom;
    }
    my %seen;
    return grep { !$seen{$_}++ } @keys;
}

sub overlap_length {
    my ($start1, $end1, $start2, $end2) = @_;
    my $start = $start1 > $start2 ? $start1 : $start2;
    my $end = $end1 < $end2 ? $end1 : $end2;
    return 0 if $end < $start;
    return $end - $start + 1;
}

sub sort_unique {
    my @values = grep { defined $_ && $_ ne '' && $_ ne 'NA' } @_;
    my %seen;
    return sort grep { !$seen{$_}++ } @values;
}

sub sort_exon_labels {
    my @labels = @_;
    my %seen;
    @labels = grep { !$seen{$_}++ } @labels;
    return sort {
        my ($a_num) = $a =~ /^EX(\d+)$/;
        my ($b_num) = $b =~ /^EX(\d+)$/;
        defined $a_num && defined $b_num
            ? $a_num <=> $b_num
            : $a cmp $b;
    } @labels;
}

sub is_positive_integer {
    my ($value) = @_;
    return defined $value && $value =~ /^\d+$/ && $value > 0;
}

sub is_missing {
    my ($value) = @_;
    return 1 unless defined $value;
    return 1 if $value eq '' || $value eq 'NA' || $value eq '.';
    return 0;
}

sub clean_tsv_value {
    my ($value) = @_;
    $value = 'NA' unless defined $value && $value ne '';
    $value =~ s/[\r\n\t]/ /g;
    return $value;
}

sub usage {
    return <<'USAGE';
Usage:
  perl bin/AnnotateCandidates.pl \
    --Config /absolute/path/to/hcm_exondel.conf \
    --GeneAnnotation /absolute/path/to/MANERefSeqGenes.tsv \
    --ExonAnnotation /absolute/path/to/MANERefSeqExons.tsv \
    --Input /absolute/path/to/merged_candidates.tsv \
    --Out /absolute/path/to/annotated_candidates.tsv

Purpose:
  Annotate split-read-centered candidates using MANE gene and exon intervals.

Region_Type definitions:
  Exonic     Candidate overlaps at least one annotated exon.
  Intronic   Candidate overlaps a gene interval but no exon interval.
  Flanking   Candidate does not overlap a gene and is within TARGET_REGION_FLANK bp
             of the nearest MANE gene interval.
  Intergenic Candidate is outside all gene intervals and gene flanks.

Affected exon definition:
  An exon is affected when candidate overlap / exon length is at least
  MIN_EXON_OVERLAP_FRACTION. Exon labels are reported as EX1,EX2,...
USAGE
}
