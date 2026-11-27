
#!/usr/bin/env perl

use strict;
use warnings;

use File::Basename qw(basename);
use File::Path qw(make_path);
use Cwd qw(abs_path);

@ARGV == 4 or die usage();

my ($bam_list, $outdir, $ref, $manta) = @ARGV;

# Validate input
die "[ERROR] BAM list not found: $bam_list\n"
    unless -f $bam_list;

die "[ERROR] Reference fasta not found: $ref\n"
    unless -f $ref;

die "[ERROR] Manta directory not found: $manta\n"
    unless -d $manta;

die "[ERROR] Manta configManta.py not found\n"
    unless -f "$manta/bin/configManta.py";

# Resolve absolute paths
$bam_list = abs_path($bam_list);
$ref      = abs_path($ref);
$manta    = abs_path($manta);

make_path($outdir) unless -d $outdir;
$outdir = abs_path($outdir);

open my $IN, '<', $bam_list
    or die "[ERROR] Cannot open $bam_list: $!\n";

while (my $bam = <$IN>) {

    chomp $bam;
    $bam =~ s/^\s+|\s+$//g;

    next if $bam eq "";
    next if $bam =~ /^#/;

    die "[ERROR] BAM file not found: $bam\n"
        unless -f $bam;

    # Resolve BAM path
    my $bam_abs  = abs_path($bam);
    my $bam_base = basename($bam);

    # Extract sample ID
    # Example:
    # 25B09089386.target.bam
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

    # Manta workflow directory
    my $run_dir = "$sample_dir/manta_run";

    # Output VCF
    my $vcf_gz = "$sample_dir/$sample.manta.DEL.vcf.gz";
    my $vcf    = "$sample_dir/$sample.manta.DEL.vcf";

    # Shell script
    my $sh_file = "$sample_dir/run_manta_$sample.sh";

    open my $SH, '>', $sh_file
        or die "[ERROR] Cannot write $sh_file: $!\n";

    # Quote shell arguments
    my $q_manta    = shell_quote($manta);
    my $q_bam      = shell_quote($bam_abs);
    my $q_ref      = shell_quote($ref);
    my $q_sample   = shell_quote($sample);
    my $q_run_dir  = shell_quote($run_dir);
    my $q_vcf_gz   = shell_quote($vcf_gz);
    my $q_vcf      = shell_quote($vcf);

    print {$SH} <<"EOF";
#!/usr/bin/env bash
set -euo pipefail

# Change to the directory containing this script
cd -- "\$(dirname -- "\${BASH_SOURCE[0]}")"

MANTA=$q_manta
BAM=$q_bam
REF=$q_ref
SAMPLE=$q_sample
RUN_DIR=$q_run_dir
VCF_GZ=$q_vcf_gz
VCF=$q_vcf

# Index BAM
if [[ ! -f "\${BAM}.bai" && ! -f "\${BAM%.bam}.bai" ]]; then
    samtools index "\${BAM}"
fi

# Index reference
if [[ ! -f "\${REF}.fai" ]]; then
    samtools faidx "\${REF}"
fi

# Configure Manta workflow
if [[ ! -f "\${RUN_DIR}/runWorkflow.py" ]]; then
    "\${MANTA}/bin/configManta.py" \\
        --bam "\${BAM}" \\
        --referenceFasta "\${REF}" \\
        --runDir "\${RUN_DIR}"
fi

# Run Manta
"\${RUN_DIR}/runWorkflow.py" -j 8

# Extract DEL variants
bcftools view \\
    -i 'INFO/SVTYPE="DEL"' \\
    "\${RUN_DIR}/results/variants/diploidSV.vcf.gz" \\
    -Oz \\
    -o "\${VCF_GZ}"

# Index VCF
tabix -p vcf "\${VCF_GZ}"

# Generate uncompressed VCF
bcftools view "\${VCF_GZ}" > "\${VCF}"

echo "Manta finished: \${VCF}"

EOF

    close $SH
        or die "[ERROR] Cannot close $sh_file: $!\n";

    chmod 0755, $sh_file
        or die "[ERROR] Cannot chmod $sh_file: $!\n";

    my $display_dir = defined $gene
        ? "$gene/$sample"
        : $sample;

    print "[OK] $display_dir/run_manta_$sample.sh\n";
}

close $IN;

print "Done! Manta shell scripts generated successfully.\n";

sub shell_quote {
    my ($s) = @_;

    $s =~ s/'/'"'"'/g;

    return "'$s'";
}

sub usage {
    return <<'USAGE';

Usage:
  perl MakeMantaShell.pl <BAM_LIST> <OUTDIR> <REF.fa> <MANTA_DIR>

Example:
  perl MakeMantaShell.pl \
    all.bam \
    Manta \
    hg19.fa \
    /path/to/manta-1.6.0.centos6_x86_64

USAGE
}


