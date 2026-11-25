#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);

my $file = "$Bin/../../db/MANE/MANE.GRCh38.v1.4.Summary.txt";

print "Gene\tTranscript\tMANE_status\n";

open my $IN, "<", $file or die "[ERROR] Cannot open MANE summary: $file\n";
<$IN>;
while (my $line = <$IN>) {
    chomp $line;
    my @fields = split /\t/, $line, -1;
    next unless @fields >= 6;

    my $gene = $fields[3] // "";
    my $transcript = $fields[5] // "";
    my $status = $fields[-5] // "";

    $transcript =~ s/\.\d+$//;
    next if $gene eq "" || $transcript eq "";

    print join("\t", $gene, $transcript, $status), "\n";
}
close $IN;
