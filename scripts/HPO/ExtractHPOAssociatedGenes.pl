
#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use FindBin qw($Bin);
use File::Path qw(make_path);
use Cwd qw(abs_path);

my ($gene_file, $hpo_file, $outdir, $help);
$outdir = ".";
$hpo_file = "$Bin/../../db/Genes/CardiomyopathyRelatedHPO/HPOTerms.tsv";

GetOptions(
    "gene-file=s" => \$gene_file,
    "term-file=s" => \$hpo_file,
    "outdir=s"    => \$outdir,
    "help|h"      => \$help,
) or die usage();

if ($help) {
    print usage();
    exit 0;
}

die usage() unless defined $gene_file && -s $gene_file;
make_path($outdir) unless -d $outdir;
$outdir = abs_path($outdir) || $outdir;


############################################################
# 1. Load phenotype_to_genes.txt
############################################################

open(my $GIN, "<", $gene_file)
    or die "Cannot open $gene_file: $!\n";

my %hpo_genes;

while (my $line = <$GIN>) {
    chomp $line;
    $line =~ s/\r$//;

    next if $line =~ /^\s*$/;
    next if $line =~ /^\s*hpo_id\b/i;

    my ($hpo_id, $hpo_name, $gene_id, $gene, $disease_id);

    if ($line =~ /^\s*(HP:\d+)\s+(.+?)\s+(\d+)\s+(\S+)\s+(\S+)\s*$/) {
        ($hpo_id, $hpo_name, $gene_id, $gene, $disease_id) =
            ($1, $2, $3, $4, $5);
    } else {
        warn "Skip malformed gene annotation: $line\n";
        next;
    }

    $hpo_name =~ s/^\s+|\s+$//g;
    $hpo_name =~ s/\s+/ /g;

    push @{$hpo_genes{$hpo_name}}, [$gene, $disease_id];
}

close($GIN);

############################################################
# 2. Read HPOTerms.tsv and output one file per HPO
############################################################

open(my $HIN, "<", $hpo_file)
    or die "Cannot open $hpo_file: $!\n";

my %seen_hpo;

while (my $line = <$HIN>) {
    chomp $line;
    $line =~ s/\r$//;

    next if $line =~ /^\s*$/;
    next if $line =~ /^\s*hpo_id\b/i;

    my ($hpo_id, $hpo_name) =
        $line =~ /^\s*(\S+)\s+(.+?)\s*$/;

    unless (defined $hpo_id && defined $hpo_name) {
        warn "Skip malformed HPO term: $line\n";
        next;
    }

    $hpo_name =~ s/^\s+|\s+$//g;
    $hpo_name =~ s/\s+/ /g;

    next if $seen_hpo{$hpo_name}++;

    # Generate filename
    my $filename = $hpo_name;
    $filename =~ s/[^A-Za-z0-9]+/_/g;
    $filename =~ s/^_+|_+$//g;
    $filename .= ".genes";

    open(my $OUT, ">", "$outdir/$filename")
        or die "Cannot create $filename: $!\n";

    print $OUT "Gene\tEvidence\tNote\n";

    my %seen_gene;

    if (exists $hpo_genes{$hpo_name}) {
        for my $item (@{$hpo_genes{$hpo_name}}) {
            my ($gene, $disease_id) = @$item;

            # Remove duplicate gene-disease associations
            my $key = "$gene\t$disease_id";
            next if $seen_gene{$key}++;

            print $OUT join("\t",
                $gene,
                $hpo_name,
                $disease_id
            ), "\n";
        }
    }

    close($OUT);

    print "[DONE] $filename\n";
}

close($HIN);

print "[DONE] All HPO gene files generated.\n";



sub usage {
    return <<'USAGE';
Usage:
  perl ExtractHPOAssociatedGenes.pl \
    --gene-file /path/to/phenotype_to_genes.txt \
    [--term-file /path/to/HPOTerms.tsv] \
    [--outdir output_directory]

Required:
  --gene-file    HPO phenotype-to-gene association table.

Optional:
  --term-file    HPO term table. Defaults to the bundled resource.
  --outdir       Output directory. Default: current directory.
  --help         Show this help message.
USAGE
}
