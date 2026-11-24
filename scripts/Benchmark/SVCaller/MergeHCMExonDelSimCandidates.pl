#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename qw(basename dirname);
use File::Spec;
use Cwd qw(abs_path);

# ============================================================
# MergeHCMExonDelSimCandidates.pl
#
# Purpose:
#   Merge HCMExonDel simulation candidate TSV files:
#
#   results/
#   ├── HCM_VCL_DEL_50_100_rep1/
#   │   └── S5_Candidates/
#   │       └── merged_candidates.tsv
#   ├── HCM_VCL_DEL_50_100_rep2/
#   │   └── S5_Candidates/
#   │       └── merged_candidates.tsv
#   └── ...
#
# Sample name example:
#   HCM_VCL_DEL_50_100_rep3
#
# Parsed information:
#   Sim_Gene      = VCL
#   Sim_DEL_Min   = 50
#   Sim_DEL_Max   = 100
#   Replicate     = 3
#
# Output:
#   One merged TSV file containing all HCMExonDel candidates.
# ============================================================


# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------

my $input_dir = ".";
my $output    = "HCMExonDel.sim_candidates.merged.tsv";
my $help      = 0;

GetOptions(
    "input-dir|i=s" => \$input_dir,
    "output|o=s"    => \$output,
    "help|h"        => \$help,
) or usage();

usage() if $help;


# ------------------------------------------------------------
# Check input directory
# ------------------------------------------------------------

die "ERROR: input directory does not exist: $input_dir\n"
    unless -d $input_dir;

$input_dir = abs_path($input_dir);


# ------------------------------------------------------------
# Find candidate files
# ------------------------------------------------------------

my @files = glob(
    File::Spec->catfile(
        $input_dir,
        "*",
        "S5_Candidates",
        "merged_candidates.tsv"
    )
);

die "ERROR: no merged_candidates.tsv files found under:\n"
    . "  $input_dir/*/S5_Candidates/merged_candidates.tsv\n"
    unless @files;


# ------------------------------------------------------------
# Parse sample information for sorting
# ------------------------------------------------------------

my @samples;

for my $file (@files) {

    #
    # file:
    # results/HCM_VCL_DEL_50_100_rep3/S5_Candidates/merged_candidates.tsv
    #
    # dirname(file)
    #   -> .../HCM_VCL_DEL_50_100_rep3/S5_Candidates
    #
    # dirname(dirname(file))
    #   -> .../HCM_VCL_DEL_50_100_rep3
    #

    my $sample_dir = dirname(dirname($file));
    my $sample_id  = basename($sample_dir);

    my (
        $sim_gene,
        $del_min,
        $del_max,
        $replicate
    );

    if (
        $sample_id =~
        /^HCM_(.+?)_DEL_(\d+)_(\d+)_rep(\d+)$/
    ) {
        $sim_gene  = $1;
        $del_min   = $2;
        $del_max   = $3;
        $replicate = $4;
    }
    else {
        warn "WARNING: cannot parse sample name: $sample_id\n";

        $sim_gene  = "NA";
        $del_min   = 0;
        $del_max   = 0;
        $replicate = 0;
    }

    push @samples, {
        file       => $file,
        sample_id  => $sample_id,
        sim_gene   => $sim_gene,
        del_min    => $del_min,
        del_max    => $del_max,
        replicate  => $replicate,
    };
}


# ------------------------------------------------------------
# Sort:
#
#   Gene
#   DEL minimum size
#   DEL maximum size
#   replicate
# ------------------------------------------------------------

@samples = sort {

       $a->{sim_gene} cmp $b->{sim_gene}

    || $a->{del_min} <=> $b->{del_min}

    || $a->{del_max} <=> $b->{del_max}

    || $a->{replicate} <=> $b->{replicate}

} @samples;


# ------------------------------------------------------------
# Open output
# ------------------------------------------------------------

open my $OUT, ">", $output
    or die "ERROR: cannot write $output: $!\n";


# ------------------------------------------------------------
# Merge files
# ------------------------------------------------------------

my $master_header;

my $file_count      = 0;
my $candidate_count = 0;
my $empty_count     = 0;


for my $sample (@samples) {

    my $file      = $sample->{file};
    my $sample_id = $sample->{sample_id};

    open my $IN, "<", $file
        or die "ERROR: cannot open $file: $!\n";

    my $header;
    my $row_count = 0;

    while (my $line = <$IN>) {

        chomp $line;
        $line =~ s/\r$//;

        next if $line =~ /^\s*$/;


        # ----------------------------------------------------
        # Header
        # ----------------------------------------------------

        if (!defined $header) {

            $header = $line;

            my @header_fields = split /\t/, $header, -1;

            die "ERROR: empty header in file:\n$file\n"
                unless @header_fields;


            # First file: create output header
            if (!defined $master_header) {

                $master_header = $header;

                #
                # Original:
                #
                # SampleID Gene Chrom Cluster_ID ...
                #
                # Output:
                #
                # SampleID
                # Sim_Gene
                # Sim_DEL_Min
                # Sim_DEL_Max
                # Replicate
                # Gene
                # Chrom
                # ...
                #

                my $first_column = shift @header_fields;

                print $OUT join(
                    "\t",
                    $first_column,
                    "Sim_Gene",
                    "Sim_DEL_Min",
                    "Sim_DEL_Max",
                    "Replicate",
                    @header_fields
                ), "\n";
            }
            else {

                #
                # Make sure all candidate files have
                # identical columns
                #
                if ($header ne $master_header) {
                    die "\nERROR: header mismatch:\n"
                      . "File: $file\n\n"
                      . "Expected:\n$master_header\n\n"
                      . "Observed:\n$header\n";
                }
            }

            next;
        }


        # ----------------------------------------------------
        # Candidate record
        # ----------------------------------------------------

        my @fields = split /\t/, $line, -1;

        next unless @fields;

        my $original_sample_id = shift @fields;


        # Normally original SampleID should equal directory name.
        # Use original value if it exists.
        if (
            !defined $original_sample_id
            || $original_sample_id eq ""
        ) {
            $original_sample_id = $sample_id;
        }


        print $OUT join(
            "\t",

            $original_sample_id,

            $sample->{sim_gene},
            $sample->{del_min},
            $sample->{del_max},
            $sample->{replicate},

            @fields

        ), "\n";

        $row_count++;
        $candidate_count++;
    }

    close $IN;

    $file_count++;

    if ($row_count == 0) {

        $empty_count++;

        warn sprintf(
            "WARNING: no candidate records: %s\n",
            $sample_id
        );
    }
}


close $OUT;


# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

print "\n";
print "============================================================\n";
print "HCMExonDel simulation candidate merge\n";
print "============================================================\n";
print "Input directory     : $input_dir\n";
print "Input files         : $file_count\n";
print "Candidate records   : $candidate_count\n";
print "Empty samples       : $empty_count\n";
print "Output              : $output\n";
print "============================================================\n";


# ============================================================
# Usage
# ============================================================

sub usage {

    print <<"USAGE";

Usage:

  perl MergeHCMExonDelSimCandidates.pl [options]

Options:

  -i, --input-dir DIR
      HCMExonDel results directory.

      Default:
          .

  -o, --output FILE
      Output merged TSV file.

      Default:
          HCMExonDel.sim_candidates.merged.tsv

  -h, --help
      Show this help message.


Example:

  perl MergeHCMExonDelSimCandidates.pl \\
      --input-dir /path/to/SVCallerBenchmark/HCMExonDel/results \\
      --output HCMExonDel.sim_candidates.merged.tsv

USAGE

    exit;
}

