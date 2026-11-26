
#!/usr/bin/env perl

use strict;
use warnings;

use File::Basename qw(basename);
use File::Path qw(make_path);
use Cwd qw(abs_path);
use Getopt::Long qw(GetOptions);

my ($outdir, $bam_list, $SAMTOOLS, $EXTRACT_SPLIT, $LUMPYEXPRESS, $help);

$SAMTOOLS      = "samtools";
$EXTRACT_SPLIT = "extractSplitReads_BwaMem";
$LUMPYEXPRESS  = "lumpyexpress";

GetOptions(
    "outdir=s"        => \$outdir,
    "bam-list=s"      => \$bam_list,
    "samtools=s"      => \$SAMTOOLS,
    "extract-split=s" => \$EXTRACT_SPLIT,
    "lumpyexpress=s"  => \$LUMPYEXPRESS,
    "help|h"          => \$help,
) or die usage();

if ($help) {
    print usage();
    exit 0;
}

die usage() unless defined $outdir && defined $bam_list;

die "[ERROR] BAM list not found: $bam_list\n"
    unless -s $bam_list;

# Resolve output directory
make_path($outdir) unless -d $outdir;
$outdir = abs_path($outdir);

# Resolve BAM list
$bam_list = abs_path($bam_list);

open my $IN, '<', $bam_list
    or die "[ERROR] Cannot open BAM list: $bam_list: $!\n";

while (my $bam = <$IN>) {

    chomp $bam;
    $bam =~ s/^\s+|\s+$//g;

    next if $bam eq "";
    next if $bam =~ /^#/;

    unless ($bam =~ /\.bam$/i) {
        warn "[WARN] Skip non-BAM line: $bam\n";
        next;
    }

    die "[ERROR] BAM not found: $bam\n"
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

    die "[ERROR] Cannot extract sample ID from: $bam\n"
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

    # Output files (relative to shell directory)
    my $discordants_bam = "$sample.discordants.bam";
    my $split_bam       = "$sample.split.bam";
    my $vcf             = "$sample.lumpy.vcf";

    my $shell = "$sample_dir/run_lumpyexpress_$sample.sh";

    open my $SH, '>', $shell
        or die "[ERROR] Cannot write shell script: $shell: $!\n";

    # Quote shell arguments
    my $q_bam         = shell_quote($bam_abs);
    my $q_discordants = shell_quote($discordants_bam);
    my $q_split       = shell_quote($split_bam);
    my $q_vcf         = shell_quote($vcf);

    my $q_samtools      = shell_quote($SAMTOOLS);
    my $q_extract_split = shell_quote($EXTRACT_SPLIT);
    my $q_lumpyexpress  = shell_quote($LUMPYEXPRESS);

    print {$SH} <<"EOF";
#!/usr/bin/env bash
set -euo pipefail

# Change to the directory containing this script
cd -- "\$(dirname -- "\${BASH_SOURCE[0]}")"

# Extract discordant read pairs
$q_samtools view -b -F 1294 $q_bam > $q_discordants

# Extract split reads
$q_samtools view -h $q_bam | \\
    $q_extract_split -i stdin | \\
    $q_samtools view -Sb - > $q_split

# Run LUMPY Express
$q_lumpyexpress \\
    -B $q_bam \\
    -S $q_split \\
    -D $q_discordants \\
    -o $q_vcf

EOF

    close $SH
        or die "[ERROR] Cannot close shell script: $shell: $!\n";

    chmod 0755, $shell
        or die "[ERROR] Cannot chmod shell script: $shell: $!\n";

    my $display_dir = defined $gene
        ? "$gene/$sample"
        : $sample;

    print "[OK] $display_dir/run_lumpyexpress_$sample.sh\n";
}

close $IN;

print "Done! LUMPY shell scripts generated successfully.\n";

sub shell_quote {
    my ($s) = @_;

    $s =~ s/'/'"'"'/g;

    return "'$s'";
}

sub usage {
    return <<'USAGE';

Usage:
  perl MakeLumpyExpressShellBatch.pl \
    --outdir output_directory \
    --bam-list bam_list.txt \
    [--samtools samtools] \
    [--extract-split extractSplitReads_BwaMem] \
    [--lumpyexpress lumpyexpress]

Example:
  perl MakeLumpyExpressShellBatch.pl \
    --outdir JX_2 \
    --bam-list JX2.target.bam.list \
    --samtools samtools \
    --extract-split extractSplitReads_BwaMem \
    --lumpyexpress lumpyexpress

USAGE
}

