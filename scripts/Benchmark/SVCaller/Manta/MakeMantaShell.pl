#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use File::Path qw(make_path);
use Cwd qw(abs_path);

if (@ARGV < 4) {
    die "Usage: perl $0 <all.bam> <OUTDIR> <REF.fa> <MANTA_DIR>\n";
}

my ($bam_list, $outdir, $ref, $manta) = @ARGV;

$outdir = abs_path($outdir) if -d $outdir;
$ref    = abs_path($ref)    if -e $ref;
$manta  = abs_path($manta)  if -d $manta;

open my $IN, "<", $bam_list or die "Cannot open $bam_list: $!\n";

while (my $bam = <$IN>) {
    chomp $bam;
    next if $bam =~ /^\s*$/;
    next if $bam =~ /^\s*#/;

    die "BAM file not found: $bam\n" unless -e $bam;

    my $bam_base = basename($bam);

    # Example input BAM filename:
    # Example: HCM_ACTN2_DEL_1000_2000_rep1.sorted.markdup.bam
    my $sample = $bam_base;
    $sample =~ s/\.sorted\.markdup\.bam$//;
    $sample =~ s/\.bam$//;

    # Extract the gene symbol from the sample name.
    # Example: HCM_ACTN2_DEL_1000_2000_rep1 -> ACTN2
    my $gene;
    if ($sample =~ /^HCM_([^_]+)_DEL_/) {
        $gene = $1;
    } else {
        die "Cannot parse gene name from sample: $sample\n";
    }

    my $run_dir = "$outdir/$gene/$sample";

    make_path($run_dir) unless -d $run_dir;

    my $sh_file = "$run_dir/run_manta_${sample}.sh";

    open my $SH, ">", $sh_file or die "Cannot write $sh_file: $!\n";

    print $SH <<"EOF";
#!/bin/bash
set -euo pipefail

MANTA=$manta
BAM=$bam
REF=$ref
SAMPLE=$sample
RUN_DIR=$run_dir

mkdir -p \${RUN_DIR}

if [ ! -f "\${BAM}.bai" ]; then
    samtools index \${BAM}
fi

if [ ! -f "\${REF}.fai" ]; then
    samtools faidx \${REF}
fi

\${MANTA}/bin/configManta.py \\
  --bam \${BAM} \\
  --referenceFasta \${REF} \\
  --runDir \${RUN_DIR}

\${RUN_DIR}/runWorkflow.py -j 8

bcftools view \\
  -i 'INFO/SVTYPE="DEL"' \\
  \${RUN_DIR}/results/variants/diploidSV.vcf.gz \\
  -Oz \\
  -o \${RUN_DIR}/\${SAMPLE}.manta.DEL.vcf.gz

tabix -p vcf \${RUN_DIR}/\${SAMPLE}.manta.DEL.vcf.gz

bcftools view \\
  \${RUN_DIR}/\${SAMPLE}.manta.DEL.vcf.gz \\
  > \${RUN_DIR}/\${SAMPLE}.manta.DEL.vcf

echo "Manta finished: \${RUN_DIR}/\${SAMPLE}.manta.DEL.vcf"
EOF

    close $SH;

    chmod 0755, $sh_file;

    print "Generated: $sh_file\n";
}

close $IN;

print "All Manta shell scripts generated successfully.\n";

