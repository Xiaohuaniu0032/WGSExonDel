
#!/usr/bin/env perl

use strict;
use warnings;

@ARGV == 2 or die "Usage: perl $0 <VCF_List> <Output_TSV>\n";

my ($vcf_list, $output_file) = @ARGV;

open my $LIST, '<', $vcf_list
    or die "[ERROR] Cannot open '$vcf_list': $!\n";

open my $OUT, '>', $output_file
    or die "[ERROR] Cannot open '$output_file': $!\n";

# Output columns
my @columns = qw(
    SampleID
    Variant_Status
    Chrom
    Pos
    ID
    Ref
    Alt
    Qual
    Filter
    SVType
    End
    SVLen
    PE
    SR
    MAPQ
    SRMAPQ
    CIPOS
    CIEND
    CT
    GT
    GQ
    FT
    DR
    DV
    RR
    RV
);

print {$OUT} join("\t", @columns), "\n";

my ($file_count, $positive_count, $negative_count, $variant_count) = (0, 0, 0, 0);

while (my $vcf = <$LIST>) {

    chomp $vcf;
    $vcf =~ s/\r$//;
    $vcf =~ s/^\s+|\s+$//g;

    next if $vcf eq '';
    next if $vcf =~ /^#/;

    die "[ERROR] VCF file not found: $vcf\n"
        unless -f $vcf;

    open my $IN, '<', $vcf
        or die "[ERROR] Cannot open '$vcf': $!\n";

    my @samples;
    my @sample_indices;
    my $has_header = 0;

    my %has_del;

    while (my $line = <$IN>) {

        chomp $line;
        $line =~ s/\r$//;

        next if $line eq '';

        # Parse VCF header
        if ($line =~ /^#CHROM(?:\t|\s)/) {

            my @header = split_vcf_fields($line);

            die "[ERROR] No sample column in VCF: $vcf\n"
                unless @header >= 10;

            for my $i (9 .. $#header) {
                push @samples, $header[$i];
                push @sample_indices, $i;
                $has_del{$header[$i]} = 0;
            }

            $has_header = 1;
            next;
        }

        next if $line =~ /^#/;

        die "[ERROR] VCF header not found: $vcf\n"
            unless $has_header;

        my @fields = split_vcf_fields($line);

        die "[ERROR] Invalid VCF record: $vcf\n"
            unless @fields >= 8;

        # Parse INFO
        my %info;

        if ($fields[7] ne '.') {
            for my $item (split /;/, $fields[7]) {
                my ($key, $value) = split /=/, $item, 2;
                $info{$key} = defined $value ? $value : 1;
            }
        }

        # Keep DEL variants only
        my $svtype = $info{SVTYPE};

        $svtype = 'DEL'
            if !defined $svtype && $fields[4] =~ /<DEL>/;

        next unless defined $svtype && $svtype eq 'DEL';

        # Parse FORMAT
        my @format_keys = @fields >= 9
            ? split(/:/, $fields[8], -1)
            : ();

        my %format_index;

        for my $i (0 .. $#format_keys) {
            $format_index{$format_keys[$i]} = $i;
        }

        # Process each sample
        for my $i (0 .. $#samples) {

            my $sample_id = $samples[$i];
            my $sample_index = $sample_indices[$i];

            die "[ERROR] Missing sample data: $vcf\n"
                unless @fields > $sample_index;

            my @values = split /:/, $fields[$sample_index], -1;

            my %format;

            for my $key (keys %format_index) {
                $format{$key} = $values[$format_index{$key}]
                    if defined $values[$format_index{$key}];
            }

            # Keep samples with alternate genotype
            next unless has_alt_genotype($format{GT});

            my @row = (
                $sample_id,
                'DEL',
                value_or_na($fields[0]),
                value_or_na($fields[1]),
                value_or_na($fields[2]),
                value_or_na($fields[3]),
                value_or_na($fields[4]),
                value_or_na($fields[5]),
                value_or_na($fields[6]),
                'DEL',
                value_or_na($info{END}),
                value_or_na($info{SVLEN}),
                value_or_na($info{PE}),
                value_or_na($info{SR}),
                value_or_na($info{MAPQ}),
                value_or_na($info{SRMAPQ}),
                value_or_na($info{CIPOS}),
                value_or_na($info{CIEND}),
                value_or_na($info{CT}),
                value_or_na($format{GT}),
                value_or_na($format{GQ}),
                value_or_na($format{FT}),
                value_or_na($format{DR}),
                value_or_na($format{DV}),
                value_or_na($format{RR}),
                value_or_na($format{RV}),
            );

            print {$OUT} join("\t", @row), "\n";

            $has_del{$sample_id} = 1;
            $variant_count++;
        }
    }

    close $IN;

    die "[ERROR] VCF header not found: $vcf\n"
        unless $has_header;

    # Output negative samples
    for my $sample_id (@samples) {

        next if $has_del{$sample_id};

        my @row = ('NA') x scalar(@columns);

        $row[0] = $sample_id;
        $row[1] = 'Neg';

        print {$OUT} join("\t", @row), "\n";

        $negative_count++;
    }

    $positive_count += scalar grep { $has_del{$_} } @samples;

    $file_count++;
}

close $LIST;

die "[ERROR] No valid VCF files found.\n"
    unless $file_count > 0;

# Append column descriptions
my @descriptions = (
    'Sample identifier extracted from the VCF sample column.',
    'Variant status: DEL indicates deletion detected; Neg indicates no DEL detected.',
    'Chromosome.',
    'Variant start position (1-based).',
    'Variant identifier.',
    'Reference allele.',
    'Alternate allele.',
    'Phred-scaled variant quality score.',
    'Variant filtering status.',
    'Structural variant type.',
    'Variant end position.',
    'Structural variant length (bp). Deletions are usually negative.',
    'Number of paired-end read pairs supporting the deletion.',
    'Number of split reads supporting the deletion.',
    'Mapping quality reported by DELLY.',
    'Split-read mapping quality reported by DELLY.',
    'Confidence interval offsets around the start position.',
    'Confidence interval offsets around the end position.',
    'Breakpoint connection type.',
    'Sample genotype.',
    'Genotype quality.',
    'Genotype-level filtering status.',
    'Number of reference-supporting discordant read pairs.',
    'Number of variant-supporting discordant read pairs.',
    'Number of reference-supporting split reads.',
    'Number of variant-supporting split reads.'
);

print {$OUT} "\n# Column descriptions\n";

for my $i (0 .. $#columns) {
    print {$OUT} "# $columns[$i]: $descriptions[$i]\n";
}

close $OUT
    or die "[ERROR] Cannot close output file: $!\n";

print "Done!\n";
print "Input VCF files:  $file_count\n";
print "Positive samples: $positive_count\n";
print "Negative samples: $negative_count\n";
print "DEL records:      $variant_count\n";
print "Output file:      $output_file\n";


# Split VCF fields
sub split_vcf_fields {
    my ($line) = @_;

    return split /\t/, $line, -1 if index($line, "\t") >= 0;

    return split /\s+/, $line;
}


# Check whether genotype contains an alternate allele
sub has_alt_genotype {
    my ($gt) = @_;

    return 0 unless defined $gt && $gt ne '';

    my @alleles = split m{[/|]}, $gt;

    for my $allele (@alleles) {
        return 1 if $allele =~ /^\d+$/ && $allele > 0;
    }

    return 0;
}


# Convert missing values to NA
sub value_or_na {
    my ($value) = @_;

    return 'NA' unless defined $value;
    return 'NA' if $value eq '' || $value eq '.';

    return $value;
}


