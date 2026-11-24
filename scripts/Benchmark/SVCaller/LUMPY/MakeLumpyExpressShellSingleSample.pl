#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Path qw(make_path);

my ($bam, $sample, $outdir, $samtools, $extract_split, $lumpyexpress, $help);
$samtools = $ENV{SAMTOOLS} || "samtools";
$extract_split = $ENV{EXTRACT_SPLIT_READS} || "extractSplitReads_BwaMem";
$lumpyexpress = $ENV{LUMPYEXPRESS} || "lumpyexpress";

GetOptions(
    "bam=s"           => \$bam,
    "sample=s"        => \$sample,
    "outdir=s"        => \$outdir,
    "samtools=s"      => \$samtools,
    "extract-split=s" => \$extract_split,
    "lumpyexpress=s"  => \$lumpyexpress,
    "help|h"          => \$help,
) or die usage();

if ($help) {
    print usage();
    exit 0;
}

die usage() unless defined $bam && defined $sample && defined $outdir;
make_path($outdir) unless -d $outdir;

my $discordants_bam = "$outdir/$sample.discordants.bam";
my $split_bam = "$outdir/$sample.split.bam";
my $vcf = "$outdir/$sample.lumpy.vcf";
my $runsh = "$outdir/run_lumpyexpress_$sample.sh";

open my $SH, ">", $runsh or die "Cannot write $runsh: $!\n";
print $SH "#!/usr/bin/env bash\nset -euo pipefail\n\n";
print $SH "$samtools view -b -F 1294 '$bam' > '$discordants_bam'\n";
print $SH "$samtools view -h '$bam' | $extract_split -i stdin | $samtools view -Sb - > '$split_bam'\n";
print $SH "$lumpyexpress -B '$bam' -S '$split_bam' -D '$discordants_bam' -o '$vcf'\n";
close $SH;
chmod 0755, $runsh or die "Cannot chmod $runsh: $!\n";

sub usage {
    return <<'USAGE';
Usage:
  perl MakeLumpyExpressShellSingleSample.pl \
    --bam sample.bam --sample SAMPLE --outdir output_directory \
    [--samtools samtools] [--extract-split extractSplitReads_BwaMem] \
    [--lumpyexpress lumpyexpress]
USAGE
}
