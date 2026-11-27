
#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename qw(basename dirname);

@ARGV == 2 or die "Usage: perl $0 <InputList> <OutputFile>\n";

my ($input_list, $output_file) = @ARGV;

open my $LIST, '<', $input_list
    or die "[ERROR] Cannot open '$input_list': $!\n";

my @files;

while (my $file = <$LIST>) {

    chomp $file;
    $file =~ s/\r$//;
    $file =~ s/^\s+|\s+$//g;

    next if $file eq '';
    next if $file =~ /^#/;

    die "[ERROR] Input file not found: $file\n"
        unless -f $file;

    push @files, $file;
}

close $LIST;

die "[ERROR] No input files found.\n" unless @files;

# Read header from the first non-empty input file
my $header;

for my $file (@files) {

    open my $IN, '<', $file
        or die "[ERROR] Cannot open '$file': $!\n";

    while (my $line = <$IN>) {
        chomp $line;
        $line =~ s/\r$//;

        next if $line eq '';

        $header = $line;
        last;
    }

    close $IN;

    last if defined $header;
}

die "[ERROR] Cannot find header in input files.\n"
    unless defined $header;

my @header_fields = split /\t/, $header, -1;

# Identify column indexes
my %column;

for my $i (0 .. $#header_fields) {
    $column{$header_fields[$i]} = $i;
}

for my $required (qw(SampleID Candidate_ID)) {
    die "[ERROR] Required column '$required' not found.\n"
        unless exists $column{$required};
}

open my $OUT, '>', $output_file
    or die "[ERROR] Cannot open '$output_file': $!\n";

# Write header once
print {$OUT} "$header\n";

my $file_count     = 0;
my $positive_count = 0;
my $negative_count = 0;
my $row_count      = 0;

for my $file (@files) {

    open my $IN, '<', $file
        or die "[ERROR] Cannot open '$file': $!\n";

    my $file_header;
    my $has_data = 0;

    while (my $line = <$IN>) {

        chomp $line;
        $line =~ s/\r$//;

        next if $line eq '';

        # Check header
        unless (defined $file_header) {

            $file_header = $line;

            die "[ERROR] Inconsistent header: $file\n"
                unless $file_header eq $header;

            next;
        }

        # Skip repeated header
        next if $line eq $header;

        # Validate column count
        my @fields = split /\t/, $line, -1;

        die "[ERROR] Column count mismatch: $file\n"
            unless @fields == @header_fields;

        print {$OUT} "$line\n";

        $has_data = 1;
        $row_count++;
    }

    close $IN;

    if ($has_data) {

        $positive_count++;

    } else {

        # Extract SampleID from:
        # .../SampleID/S6_Report/annotated_candidates.tsv

        my $sample_id = basename(dirname(dirname($file)));

        # Generate negative sample record
        my @fields = ('NA') x scalar(@header_fields);

        $fields[$column{SampleID}]     = $sample_id;
        $fields[$column{Candidate_ID}] = 'Neg';

        print {$OUT} join("\t", @fields), "\n";

        $negative_count++;
    }

    $file_count++;
}

close $OUT
    or die "[ERROR] Cannot close '$output_file': $!\n";

print "Done!\n";
print "Input files:      $file_count\n";
print "Positive samples: $positive_count\n";
print "Negative samples: $negative_count\n";
print "Candidate records: $row_count\n";
print "Total records:    " . ($row_count + $negative_count) . "\n";
print "Output file:      $output_file\n";


