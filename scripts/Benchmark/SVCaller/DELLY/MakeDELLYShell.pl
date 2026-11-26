
#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename qw(basename);
use File::Path qw(make_path);
use Cwd qw(abs_path);

my ($bam_list, $outdir, $ref, $delly);

GetOptions(
    "bam-list|b=s" => \$bam_list,
    "outdir|o=s"   => \$outdir,
    "ref|r=s"      => \$ref,
    "delly|d=s"    => \$delly,
) or die usage();

die usage() unless $bam_list && $outdir && $ref;

$delly ||= "delly";

die "[ERROR] BAM list not found: $bam_list\n" unless -f $bam_list;
die "[ERROR] Reference fasta not found: $ref\n" unless -f $ref;

# Resolve absolute paths
$bam_list = abs_path($bam_list);
$ref      = abs_path($ref);

make_path($outdir) unless -d $outdir;
$outdir = abs_path($outdir);

open my $IN, '<', $bam_list
    or die "[ERROR] Cannot open $bam_list: $!\n";

while (my $bam = <$IN>) {

    chomp $bam;
    $bam =~ s/^\s+|\s+$//g;

    next if $bam eq "";
    next if $bam =~ /^#/;

    die "[ERROR] BAM not found: $bam\n" unless -f $bam;

    my $bam_abs  = abs_path($bam);
    my $bam_base = basename($bam);

    # Extract sample ID from BAM filename
    # Example:
    # 25B09089386.final.merge.bam
    # => 25B09089386
    #
    # HCM_ACTN2_DEL_1000_2000_rep1.sorted.markdup.bam
    # => HCM_ACTN2_DEL_1000_2000_rep1

    my ($sample) = split /\./, $bam_base;

    die "[ERROR] Cannot extract sample ID from BAM: $bam\n"
        unless defined $sample && $sample ne "";

    # Extract gene name from simulated sample
    my $gene;

    if ($sample =~ /^HCM_([^_]+)_DEL_/) {
        $gene = $1;
    } elsif ($sample =~ /^HCM_([^_]+)_/) {
        $gene = $1;
    }

    # Set output directory
    # Simulated data: OUTDIR/GENE/SAMPLE
    # Real WGS data:  OUTDIR/SAMPLE

    my $sample_dir = defined $gene
        ? "$outdir/$gene/$sample"
        : "$outdir/$sample";

    make_path($sample_dir) unless -d $sample_dir;

    my $sh  = "$sample_dir/run_delly_$sample.sh";
    my $vcf = "$sample.delly.vcf";

    open my $SH, '>', $sh
        or die "[ERROR] Cannot write $sh: $!\n";

    print {$SH} <<"EOF";
#!/bin/bash
set -euo pipefail

cd "$sample_dir"

"$delly" call \\
    -t DEL \\
    -g "$ref" \\
    "$bam_abs" \\
    > "$vcf"

EOF

    close $SH;

    chmod 0755, $sh;

    print "[OK] $sh\n";
}

close $IN;

print "Done! DELLY shell scripts generated successfully.\n";

sub usage {
    return <<'USAGE';

Usage:
  perl MakeDELLYShell.pl \
    --bam-list all.bam \
    --outdir /path/to/delly_out \
    --ref ref.fa

Optional:
  --delly /path/to/delly

Example:
  perl MakeDELLYShell.pl \
    --bam-list all.bam \
    --outdir /path/to/SVCallerBenchmark/DELLY \
    --ref /path/to/ref.fa \
    --delly delly

USAGE
}


