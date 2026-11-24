#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename qw(basename);
use File::Path qw(make_path);
use Getopt::Long qw(GetOptions);

# Generate one LUMPY shell script per BAM listed in the input table.

my ($outdir, $bam_list, $SAMTOOLS, $EXTRACT_SPLIT, $LUMPYEXPRESS, $help);
$SAMTOOLS = "samtools";
$EXTRACT_SPLIT = "extractSplitReads_BwaMem";
$LUMPYEXPRESS = "lumpyexpress";

GetOptions(
    "outdir=s"       => \$outdir,
    "bam-list=s"     => \$bam_list,
    "samtools=s"     => \$SAMTOOLS,
    "extract-split=s"=> \$EXTRACT_SPLIT,
    "lumpyexpress=s" => \$LUMPYEXPRESS,
    "help|h"         => \$help,
) or die usage();

if ($help) {
    print usage();
    exit 0;
}

die usage() unless defined $outdir && defined $bam_list && -s $bam_list;

make_path($outdir) unless -d $outdir;

open my $IN, "<", $bam_list or die "Cannot open bam list: $bam_list\n";

while (my $bam = <$IN>) {
    chomp $bam;
    $bam =~ s/^\s+|\s+$//g;

    next if $bam eq "";
    next if $bam =~ /^#/;

    unless ($bam =~ /\.bam$/) {
        warn "[WARN] Skip non-BAM line: $bam\n";
        next;
    }

    my $bam_base = basename($bam);

    my $sample;
    if ($bam_base =~ /^(.+)\.sorted\.markdup\.bam$/) {
        $sample = $1;
    } elsif ($bam_base =~ /^(.+)\.bam$/) {
        $sample = $1;
    } else {
        warn "[WARN] Cannot parse sample name from: $bam\n";
        next;
    }

    my $gene;
    if ($sample =~ /^HCM_([^_]+)_DEL_/) {
        $gene = $1;
    } else {
        warn "[WARN] Cannot parse gene name from sample: $sample\n";
        next;
    }

    my $gene_dir   = "$outdir/$gene";
    my $sample_dir = "$gene_dir/$sample";

    make_path($sample_dir) unless -d $sample_dir;

    my $discordants_bam = "$sample_dir/$sample.discordants.bam";
    my $split_bam       = "$sample_dir/$sample.split.bam";
    my $vcf             = "$sample_dir/$sample.lumpy.vcf";

    my $shell = "$sample_dir/run_lumpyexpress_$sample.sh";

    open my $SH, ">", $shell or die "Cannot write shell script: $shell\n";

    my $q_bam         = shell_quote($bam);
    my $q_discordants = shell_quote($discordants_bam);
    my $q_split       = shell_quote($split_bam);
    my $q_vcf         = shell_quote($vcf);

    my $q_samtools      = shell_quote($SAMTOOLS);
    my $q_extract_split = shell_quote($EXTRACT_SPLIT);
    my $q_lumpyexpress  = shell_quote($LUMPYEXPRESS);

    print $SH <<"EOF";
#!/usr/bin/env bash
set -euo pipefail

$q_samtools view -b -F 1294 $q_bam > $q_discordants

$q_samtools view -h $q_bam | \\
$q_extract_split -i stdin | \\
$q_samtools view -Sb - > $q_split

$q_lumpyexpress \\
  -B $q_bam \\
  -S $q_split \\
  -D $q_discordants \\
  -o $q_vcf

EOF

    close $SH;

    chmod 0755, $shell or die "Cannot chmod shell script: $shell\n";

    print "[OK] $gene/$sample/run_lumpyexpress_$sample.sh\n";
}

close $IN;

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
    --bam-list bam_list.tsv \
    [--samtools samtools] \
    [--extract-split extractSplitReads_BwaMem] \
    [--lumpyexpress lumpyexpress]
USAGE
}

