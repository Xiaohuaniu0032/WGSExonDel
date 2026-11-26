
#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename qw(basename);

@ARGV == 2 or die "Usage: perl $0 <InputFile> <OutputFile>\n";

my ($input_file, $output_file) = @ARGV;

open my $IN, '<', $input_file
    or die "Error: Cannot open '$input_file': $!\n";

open my $OUT, '>', $output_file
    or die "Error: Cannot open '$output_file': $!\n";

while (my $bam_path = <$IN>) {

    chomp $bam_path;
    $bam_path =~ s/\r$//;

    next if $bam_path =~ /^\s*$/;

    my ($sample_id) = split /\./, basename($bam_path);

    print {$OUT} "$sample_id\t$bam_path\n";
}

close $IN;
close $OUT;

print "Done! Output file: $output_file\n";

