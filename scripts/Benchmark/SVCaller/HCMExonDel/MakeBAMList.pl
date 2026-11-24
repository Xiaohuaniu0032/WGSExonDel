#!/usr/bin/env perl
use strict;
use warnings;

my $input = shift @ARGV or die "Usage: perl $0 <bam_path_list>\n";
open my $IN, "<", $input or die "Cannot open $input: $!\n";

while (my $bam = <$IN>) {
    chomp $bam;
    $bam =~ s/^\s+|\s+$//g;
    next if $bam eq "";
    next if $bam =~ /^#/;

    my ($filename) = $bam =~ m{([^/]+)$};
    my $sample = $filename // "";
    $sample =~ s/\.sorted\.markdup\.bam$//;
    $sample =~ s/\.bam$//;

    print join("\t", $sample, $bam), "\n";
}

close $IN;
