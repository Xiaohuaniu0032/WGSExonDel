
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
    STRANDS
    SU
    PE
    SR
    BD
    CIPOS
    CIEND
    CIPOS95
    CIEND95
    GT
);

print {$OUT} join("\t", @columns), "\n";

my ($file_count, $positive_count, $negative_count, $variant_count) =
    (0, 0, 0, 0);

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

    my (@samples, @sample_indices);
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

                die "[ERROR] Empty sample ID in VCF: $vcf\n"
                    unless defined $header[$i] && $header[$i] ne '';

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

        my $single_sample = @samples == 1;

        # Process each sample
        for my $i (0 .. $#samples) {

            my $sample_id = $samples[$i];
            my $sample_index = $sample_indices[$i];

            die "[ERROR] Missing sample data: $vcf\n"
                unless @fields > $sample_index;

            my @values = split /:/, $fields[$sample_index], -1;

            my %format;

            for my $key (keys %format_index) {

                my $index = $format_index{$key};

                $format{$key} = $values[$index]
                    if $index <= $#values;
            }

            # Check sample-level support
            my $has_support_fields =
                grep { exists $format_index{$_} } qw(SU PE SR);

            my $has_support =
                positive_number($format{SU}) ||
                positive_number($format{PE}) ||
                positive_number($format{SR});

            # Check genotype
            my $has_alt = has_alt_genotype($format{GT});

            # For single-sample VCF, use INFO support
            # when sample-level support fields are unavailable
            if (!$has_support && !$has_alt &&
                !$has_support_fields && $single_sample) {

                $has_support =
                    positive_number($info{SU}) ||
                    positive_number($info{PE}) ||
                    positive_number($info{SR});
            }

            # Skip samples without DEL support
            next unless $has_support || $has_alt;

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
                value_or_na($info{STRANDS}),
                sample_or_info($format{SU}, $info{SU}, $single_sample),
                sample_or_info($format{PE}, $info{PE}, $single_sample),
                sample_or_info($format{SR}, $info{SR}, $single_sample),
                value_or_na($format{BD}),
                value_or_na($info{CIPOS}),
                value_or_na($info{CIEND}),
                value_or_na($info{CIPOS95}),
                value_or_na($info{CIEND95}),
                value_or_na($format{GT}),
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

# Column descriptions
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
    'Structural variant length (bp).',
    'Breakpoint orientation.',
    'Total number of supporting reads.',
    'Number of paired-end reads supporting the variant.',
    'Number of split reads supporting the variant.',
    'Amount of BED evidence supporting the variant.',
    'Confidence interval offsets around the start position.',
    'Confidence interval offsets around the end position.',
    '95% confidence interval offsets around the start position.',
    '95% confidence interval offsets around the end position.',
    'Sample genotype.'
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

    if (index($line, "\t") >= 0) {
        return split /\t/, $line, -1;
    }

    $line =~ s/^\s+|\s+$//g;

    return split /\s+/, $line;
}


# Check whether genotype contains an alternate allele
sub has_alt_genotype {

    my ($gt) = @_;

    return 0 unless defined $gt && $gt ne '' && $gt ne '.';

    my @alleles = split m{[/|]}, $gt;

    for my $allele (@alleles) {

        return 1 if $allele =~ /^\d+$/ && $allele > 0;
    }

    return 0;
}


# Check whether a support value is positive
sub positive_number {

    my ($value) = @_;

    return 0 unless defined $value;
    return 0 unless $value =~ /^\d+(?:\.\d+)?$/;

    return $value > 0;
}


# Get sample-level value, with INFO fallback for single-sample VCF
sub sample_or_info {

    my ($sample_value, $info_value, $single_sample) = @_;

    if (defined $sample_value &&
        $sample_value ne '' &&
        $sample_value ne '.') {

        return $sample_value;
    }

    return value_or_na($info_value) if $single_sample;

    return 'NA';
}


# Convert missing values to NA
sub value_or_na {

    my ($value) = @_;

    return 'NA' unless defined $value;
    return 'NA' if $value eq '' || $value eq '.';

    return $value;
}


