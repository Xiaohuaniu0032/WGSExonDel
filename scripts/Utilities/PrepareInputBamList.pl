#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename qw(basename);

@ARGV == 2 or die "Usage: perl $0 <InputFile> <OutputFile>\n";

my ($input_file, $output_file) = @ARGV;

open my $IN, '<', $input_file
    or die "Error: Cannot open input file '$input_file': $!\n";

open my $OUT, '>', $output_file
    or die "Error: Cannot open output file '$output_file': $!\n";

while (my $bam_path = <$IN>) {

    chomp $bam_path;
    $bam_path =~ s/\r$//;
    $bam_path =~ s/^\s+|\s+$//;

    next if $bam_path eq '';

    my $filename = basename($bam_path);

    my ($sample_id) = $filename =~ /^(.+)\.final\.merge\.bam$/;

    die "Error: Invalid BAM filename at input line $.: $bam_path\n"
        unless defined $sample_id;

    print {$OUT} "$sample_id\t$bam_path\n"
        or die "Error: Failed to write output file: $!\n";
}

close $IN;
close $OUT or die "Error: Failed to close output file: $!\n";

print "Done! Output file: $output_file\n";

