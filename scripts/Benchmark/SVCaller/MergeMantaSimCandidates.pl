#!/usr/bin/env perl
use strict;
use warnings;
use Cwd qw(abs_path);

# ============================================================
# Script:
#   MergeMantaSimCandidates.pl
#
# Description:
#   Merge Manta DEL VCF results from simulated deletion samples
#   into one tab-delimited XLS file.
#
#   Empty VCF samples are recorded separately for downstream
#   false-negative / sensitivity analysis.
#
# Working directory:
#   The caller-results root directory supplied as the first argument.
#
# Input:
#   Manta/*/*/*.manta.DEL.vcf
#
# Output:
#   1. Manta_sim_merged_candidates.xls
#   2. Manta_sim_empty_vcf_samples.xls
#
# Coordinate definition for DEL:
#
#   Manta VCF:
#
#       POS = left anchor base
#       END = INFO/END
#
#   Standardized actual deleted interval:
#
#       Start = POS + 1
#       End   = END
#       Size  = abs(SVLEN)
#
#   Thus Start-End is a 1-based inclusive deletion interval.
#
# Example:
#
#       POS   = 1000
#       END   = 1087
#       SVLEN = -87
#
#       Start = 1001
#       End   = 1087
#       Size  = 87
#
# ============================================================


# ============================================================
# Configuration
# ============================================================

my $sv_caller_dir = shift @ARGV || ".";
$sv_caller_dir = abs_path($sv_caller_dir) || $sv_caller_dir;

my $input_dir =
    "$sv_caller_dir/Manta";

my $output =
    "$sv_caller_dir/Manta_sim_merged_candidates.xls";

my $empty_output =
    "$sv_caller_dir/Manta_sim_empty_vcf_samples.xls";


# ============================================================
# Find Manta VCF files
# ============================================================

my @files = sort glob(
    "$input_dir/*/*/*.manta.DEL.vcf"
);

die "ERROR: No Manta DEL VCF files found:\n"
    . "$input_dir/*/*/*.manta.DEL.vcf\n"
    unless @files;


# ============================================================
# Open output files
# ============================================================

open(my $OUT, ">", $output)
    or die "ERROR: Cannot write $output: $!\n";

open(my $EMPTY, ">", $empty_output)
    or die "ERROR: Cannot write $empty_output: $!\n";


# ============================================================
# Main output header
#
# Keep consistent with real-data Manta merged table
# ============================================================

my @header = (

    # Standardized SV fields
    "SampleID",
    "Chrom",
    "Start",
    "End",
    "Size",
    "SVTYPE",

    # Original VCF coordinates
    "VCF_POS",
    "VCF_END",
    "SVLEN",

    # Basic VCF fields
    "ID",
    "QUAL",
    "FILTER",

    # Manta INFO fields
    "CIPOS",
    "CIEND",
    "CIGAR",
    "HOMLEN",
    "HOMSEQ",

    # FORMAT/sample fields
    "GT",
    "FT",
    "GQ",

    # Paired-read evidence
    "PR",
    "PR_Ref",
    "PR_Alt",

    # Split-read evidence
    "SR",
    "SR_Ref",
    "SR_Alt",

    # Genotype likelihoods
    "PL",

    # Original raw VCF fields
    "Raw_REF",
    "Raw_ALT",
    "Raw_INFO",
    "Raw_FORMAT",
    "Raw_SAMPLE"
);

print $OUT join("\t", @header), "\n";


# ============================================================
# Empty VCF output header
# ============================================================

print $EMPTY join("\t",
    "SampleID",
    "Gene",
    "Length_Bin",
    "Replicate",
    "VCF_File"
), "\n";


# ============================================================
# Statistics
# ============================================================

my $file_count     = 0;
my $record_count   = 0;
my $empty_vcf      = 0;

my $del_count      = 0;
my $non_del_count  = 0;


# ============================================================
# Run information
# ============================================================

print "\n";
print "============================================================\n";
print " Manta simulated-data merger\n";
print "============================================================\n";
print "SV caller dir      : $sv_caller_dir\n";
print "Input directory    : $input_dir\n";
print "Input VCF files    : ", scalar(@files), "\n";
print "Merged output      : $output\n";
print "Empty VCF output   : $empty_output\n";
print "============================================================\n";


# ============================================================
# Process each Manta VCF
# ============================================================

for my $file (@files) {

    $file_count++;

    print "[$file_count/", scalar(@files), "] $file\n";

    open(my $IN, "<", $file)
        or die "ERROR: Cannot open $file: $!\n";


    my $sample_id;
    my $file_record_count = 0;


    # ========================================================
    # Read VCF
    # ========================================================

    while (my $line = <$IN>) {

        chomp $line;
        $line =~ s/\r$//;


        # ----------------------------------------------------
        # Skip VCF metadata
        # ----------------------------------------------------

        next if $line =~ /^##/;


        # ----------------------------------------------------
        # Parse #CHROM header
        # ----------------------------------------------------

        if ($line =~ /^#CHROM/) {

            my @h = split(/\t/, $line, -1);

            if (@h < 10) {

                die "ERROR: No sample column found in:\n"
                    . "$file\n";
            }

            $sample_id = $h[9];

            next;
        }


        # Skip other comments
        next if $line =~ /^#/;

        # Skip blank lines
        next if $line =~ /^\s*$/;


        # ----------------------------------------------------
        # Sample ID must already be known
        # ----------------------------------------------------

        if (!defined $sample_id) {

            die "ERROR: Variant record found before "
                . "#CHROM header:\n$file\n";
        }


        # ====================================================
        # Parse VCF columns
        # ====================================================

        my @f = split(/\t/, $line, -1);

        if (@f < 10) {

            warn "WARNING: Invalid VCF record skipped:\n";
            warn "$file\n";
            warn "$line\n";

            next;
        }


        my (
            $chrom,
            $vcf_pos,
            $id,
            $ref,
            $alt,
            $qual,
            $filter,
            $info_raw,
            $format_raw,
            $sample_raw
        ) = @f[0..9];


        # ====================================================
        # Parse INFO
        # ====================================================

        my %INFO;

        for my $item (split(/;/, $info_raw)) {

            next unless length $item;

            if ($item =~ /^([^=]+)=(.*)$/) {

                $INFO{$1} = $2;

            }
            else {

                # Flag fields
                $INFO{$item} = 1;
            }
        }


        # ====================================================
        # Basic SV information
        # ====================================================

        my $svtype =
            exists $INFO{"SVTYPE"}
            ? $INFO{"SVTYPE"}
            : "NA";

        my $vcf_end =
            exists $INFO{"END"}
            ? $INFO{"END"}
            : "NA";

        my $svlen =
            exists $INFO{"SVLEN"}
            ? $INFO{"SVLEN"}
            : "NA";


        # ====================================================
        # Standardized coordinates
        # ====================================================

        my $start = $vcf_pos;
        my $end   = $vcf_end;
        my $size  = "NA";


        # ----------------------------------------------------
        # DEL
        # ----------------------------------------------------

        if ($svtype eq "DEL") {

            $del_count++;


            # Actual deletion starts after VCF anchor
            if ($vcf_pos =~ /^\d+$/) {

                $start = $vcf_pos + 1;
            }


            # Prefer SVLEN
            if ($svlen =~ /^-?\d+$/) {

                $size = abs($svlen);
            }

            # Fallback:
            # Size = END - POS
            elsif (
                $vcf_pos =~ /^\d+$/
                &&
                $vcf_end =~ /^\d+$/
            ) {

                $size = $vcf_end - $vcf_pos;
            }
        }


        # ----------------------------------------------------
        # Unexpected non-DEL record
        #
        # Current input files are *.manta.DEL.vcf, so normally
        # all records should be DEL. Keep unexpected records
        # rather than silently discarding them.
        # ----------------------------------------------------

        else {

            $non_del_count++;

            $start = $vcf_pos;
            $end   = $vcf_end;


            if ($svlen =~ /^-?\d+$/) {

                $size = abs($svlen);
            }

            elsif (
                $vcf_pos =~ /^\d+$/
                &&
                $vcf_end =~ /^\d+$/
            ) {

                $size = abs(
                    $vcf_end - $vcf_pos
                );
            }


            warn "WARNING: Non-DEL record found in DEL VCF: "
                . "$sample_id $chrom:$vcf_pos "
                . "SVTYPE=$svtype\n";
        }


        # ====================================================
        # Extract Manta INFO fields
        # ====================================================

        my $cipos =
            exists $INFO{"CIPOS"}
            ? $INFO{"CIPOS"}
            : "NA";

        my $ciend =
            exists $INFO{"CIEND"}
            ? $INFO{"CIEND"}
            : "NA";

        my $cigar =
            exists $INFO{"CIGAR"}
            ? $INFO{"CIGAR"}
            : "NA";

        my $homlen =
            exists $INFO{"HOMLEN"}
            ? $INFO{"HOMLEN"}
            : "NA";

        my $homseq =
            exists $INFO{"HOMSEQ"}
            ? $INFO{"HOMSEQ"}
            : "NA";


        # ====================================================
        # Parse FORMAT / SAMPLE
        # ====================================================

        my @format_keys =
            split(/:/, $format_raw, -1);

        my @sample_values =
            split(/:/, $sample_raw, -1);

        my %FORMAT;


        for (my $i = 0; $i < @format_keys; $i++) {

            my $key = $format_keys[$i];

            my $value =
                defined $sample_values[$i]
                ? $sample_values[$i]
                : "NA";

            $FORMAT{$key} = $value;
        }


        # ====================================================
        # Extract FORMAT fields
        # ====================================================

        my $gt =
            exists $FORMAT{"GT"}
            ? $FORMAT{"GT"}
            : "NA";

        my $ft =
            exists $FORMAT{"FT"}
            ? $FORMAT{"FT"}
            : "NA";

        my $gq =
            exists $FORMAT{"GQ"}
            ? $FORMAT{"GQ"}
            : "NA";

        my $pr =
            exists $FORMAT{"PR"}
            ? $FORMAT{"PR"}
            : "NA";

        my $sr =
            exists $FORMAT{"SR"}
            ? $FORMAT{"SR"}
            : "NA";

        my $pl =
            exists $FORMAT{"PL"}
            ? $FORMAT{"PL"}
            : "NA";


        # ====================================================
        # Split PR into reference / alternate supporting reads
        #
        # Example:
        #
        #   PR = 21,13
        #
        #   PR_Ref = 21
        #   PR_Alt = 13
        #
        # ====================================================

        my $pr_ref = "NA";
        my $pr_alt = "NA";

        if ($pr ne "NA") {

            my @tmp = split(/,/, $pr, -1);

            $pr_ref = defined $tmp[0] ? $tmp[0] : "NA";
            $pr_alt = defined $tmp[1] ? $tmp[1] : "NA";
        }


        # ====================================================
        # Split SR into reference / alternate supporting reads
        #
        # Example:
        #
        #   SR = 22,15
        #
        #   SR_Ref = 22
        #   SR_Alt = 15
        #
        # ====================================================

        my $sr_ref = "NA";
        my $sr_alt = "NA";

        if ($sr ne "NA") {

            my @tmp = split(/,/, $sr, -1);

            $sr_ref = defined $tmp[0] ? $tmp[0] : "NA";
            $sr_alt = defined $tmp[1] ? $tmp[1] : "NA";
        }


        # ====================================================
        # Write merged output
        # ====================================================

        print $OUT join("\t",

            # Standardized fields
            $sample_id,
            $chrom,
            $start,
            $end,
            $size,
            $svtype,

            # Original coordinates
            $vcf_pos,
            $vcf_end,
            $svlen,

            # Basic VCF
            $id,
            $qual,
            $filter,

            # INFO
            $cipos,
            $ciend,
            $cigar,
            $homlen,
            $homseq,

            # FORMAT
            $gt,
            $ft,
            $gq,

            # Paired-read evidence
            $pr,
            $pr_ref,
            $pr_alt,

            # Split-read evidence
            $sr,
            $sr_ref,
            $sr_alt,

            # Genotype likelihood
            $pl,

            # Raw VCF
            $ref,
            $alt,
            $info_raw,
            $format_raw,
            $sample_raw

        ), "\n";


        $record_count++;
        $file_record_count++;
    }


    close $IN;


    # ========================================================
    # Empty VCF handling
    # ========================================================

    if ($file_record_count == 0) {

        $empty_vcf++;


        # ----------------------------------------------------
        # Recover SampleID from filename if needed
        #
        # Example:
        #
        # HCM_ACTN2_DEL_1000_2000_rep1.manta.DEL.vcf
        #
        # ----------------------------------------------------

        if (!defined $sample_id || $sample_id eq "") {

            if (
                $file =~
                /\/([^\/]+)\.manta\.DEL\.vcf$/
            ) {

                $sample_id = $1;
            }

            else {

                $sample_id = "NA";
            }
        }


        # ----------------------------------------------------
        # Parse simulated-sample metadata
        #
        # Example:
        #
        # HCM_ACTN2_DEL_1000_2000_rep1
        #
        # Gene       = ACTN2
        # Length_Bin = 1000_2000
        # Replicate  = 1
        #
        # ----------------------------------------------------

        my $gene       = "NA";
        my $length_bin = "NA";
        my $replicate  = "NA";


        if (
            $sample_id =~
            /^HCM_(.+)_DEL_(\d+_\d+)_rep(\d+)$/
        ) {

            $gene       = $1;
            $length_bin = $2;
            $replicate  = $3;
        }


        # ----------------------------------------------------
        # Write empty VCF table
        # ----------------------------------------------------

        print $EMPTY join("\t",

            $sample_id,
            $gene,
            $length_bin,
            $replicate,
            $file

        ), "\n";


        warn "WARNING: Empty VCF: $sample_id\n";
    }
}


# ============================================================
# Close outputs
# ============================================================

close $OUT;
close $EMPTY;


# ============================================================
# Summary
# ============================================================

print "\n";
print "============================================================\n";
print " Merge completed\n";
print "============================================================\n";

print "VCF files processed : $file_count\n";
print "Empty VCF files     : $empty_vcf\n";
print "Non-empty VCF files : ", $file_count - $empty_vcf, "\n";
print "Total SV records    : $record_count\n";

print "------------------------------------------------------------\n";

print "DEL records         : $del_count\n";
print "Non-DEL records     : $non_del_count\n";

print "------------------------------------------------------------\n";

print "Merged output       : $output\n";
print "Empty VCF output    : $empty_output\n";

print "============================================================\n";
print "\n";

