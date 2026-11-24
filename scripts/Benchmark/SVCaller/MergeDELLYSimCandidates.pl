#!/usr/bin/env perl
use strict;
use warnings;
use Cwd qw(abs_path);

# ============================================================
# Script:
#   MergeDELLYSimCandidates.pl
#
# Description:
#   Merge DELLY VCF results from simulated deletion samples
#   into one tab-delimited XLS file.
#
#   In addition, record all empty VCF samples in a separate
#   table for downstream false-negative analysis.
#
# Working directory:
#   The caller-results root directory supplied as the first argument.
#
# Input:
#   DELLY/*/*/*.delly.vcf
#
# Output:
#
#   1. DELLY_sim_merged_candidates.xls
#
#   2. DELLY_sim_empty_vcf_samples.xls
#
# Coordinate definition:
#
#   For SVTYPE=DEL:
#
#       VCF_POS = original VCF POS
#       VCF_END = INFO/END
#
#       Start = VCF_POS + 1
#       End   = VCF_END
#       Size  = abs(SVLEN)
#
#   Start-End therefore represents the actual deleted
#   reference interval using 1-based inclusive coordinates.
#
#   Example:
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
    "$sv_caller_dir/DELLY";

my $output =
    "$sv_caller_dir/DELLY_sim_merged_candidates.xls";

my $empty_output =
    "$sv_caller_dir/DELLY_sim_empty_vcf_samples.xls";


# ============================================================
# Find input DELLY VCF files
# ============================================================

my @files = sort glob(
    "$input_dir/*/*/*.delly.vcf"
);

die "ERROR: No DELLY VCF files found:\n"
    . "$input_dir/*/*/*.delly.vcf\n"
    unless @files;


# ============================================================
# Open main output
# ============================================================

open(my $OUT, ">", $output)
    or die "ERROR: Cannot write $output: $!\n";


# ============================================================
# Open empty-VCF output
# ============================================================

open(my $EMPTY, ">", $empty_output)
    or die "ERROR: Cannot write $empty_output: $!\n";


# ============================================================
# Main output header
#
# Keep consistent with real-data DELLY merged result
# ============================================================

my @header = (

    # Standardized SV fields
    "SampleID",
    "Chrom",
    "Start",
    "End",
    "Size",
    "SVTYPE",

    # Original DELLY coordinates
    "VCF_POS",
    "VCF_END",
    "SVLEN",

    # Basic VCF fields
    "ID",
    "QUAL",
    "FILTER",

    # DELLY INFO fields
    "PRECISE",
    "PE",
    "SR",
    "MAPQ",
    "SRMAPQ",
    "CT",
    "CIPOS",
    "CIEND",
    "INSLEN",
    "HOMLEN",
    "CE",

    # FORMAT / genotype fields
    "GT",
    "GQ",
    "FT",
    "DR",
    "DV",
    "RR",
    "RV",

    # Raw original VCF fields
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
my $dup_count      = 0;
my $inv_count      = 0;
my $ins_count      = 0;
my $bnd_count      = 0;
my $other_count    = 0;


# ============================================================
# Run information
# ============================================================

print "\n";
print "============================================================\n";
print " DELLY simulated-data merger\n";
print "============================================================\n";
print "SV caller dir      : $sv_caller_dir\n";
print "Input directory    : $input_dir\n";
print "Input VCF files    : ", scalar(@files), "\n";
print "Merged output      : $output\n";
print "Empty VCF output   : $empty_output\n";
print "============================================================\n";


# ============================================================
# Process each VCF
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
        # Skip metadata
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


        # Other comments
        next if $line =~ /^#/;

        # Blank lines
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

                # Flag fields:
                # PRECISE / IMPRECISE
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
        # Standardized Start / End / Size
        # ====================================================

        my $start = $vcf_pos;
        my $end   = $vcf_end;
        my $size  = "NA";


        # ----------------------------------------------------
        # DEL
        #
        # DELLY VCF:
        #
        #   POS = left anchor
        #
        # Standardized actual deletion:
        #
        #   Start = POS + 1
        #   End   = END
        #
        # ----------------------------------------------------

        if ($svtype eq "DEL") {

            $del_count++;


            if ($vcf_pos =~ /^\d+$/) {

                $start = $vcf_pos + 1;
            }


            # Prefer SVLEN
            if ($svlen =~ /^-?\d+$/) {

                $size = abs($svlen);
            }

            # Fallback:
            #
            # Size = END - POS
            #
            elsif (
                $vcf_pos =~ /^\d+$/
                &&
                $vcf_end =~ /^\d+$/
            ) {

                $size = $vcf_end - $vcf_pos;
            }
        }


        # ----------------------------------------------------
        # Other SV types
        # ----------------------------------------------------

        else {

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


            # -----------------------------------------------
            # SV-type statistics
            # -----------------------------------------------

            if ($svtype eq "DUP") {

                $dup_count++;
            }

            elsif ($svtype eq "INV") {

                $inv_count++;
            }

            elsif ($svtype eq "INS") {

                $ins_count++;
            }

            elsif (
                $svtype eq "BND"
                ||
                $svtype eq "TRA"
            ) {

                $bnd_count++;
            }

            else {

                $other_count++;
            }
        }


        # ====================================================
        # PRECISE
        # ====================================================

        my $precise =
            exists $INFO{"PRECISE"}
            ? "Yes"
            : "No";


        # ====================================================
        # Extract DELLY INFO fields
        # ====================================================

        my $pe =
            exists $INFO{"PE"}
            ? $INFO{"PE"}
            : "NA";

        my $sr =
            exists $INFO{"SR"}
            ? $INFO{"SR"}
            : "NA";

        my $mapq =
            exists $INFO{"MAPQ"}
            ? $INFO{"MAPQ"}
            : "NA";

        my $srmapq =
            exists $INFO{"SRMAPQ"}
            ? $INFO{"SRMAPQ"}
            : "NA";

        my $ct =
            exists $INFO{"CT"}
            ? $INFO{"CT"}
            : "NA";

        my $cipos =
            exists $INFO{"CIPOS"}
            ? $INFO{"CIPOS"}
            : "NA";

        my $ciend =
            exists $INFO{"CIEND"}
            ? $INFO{"CIEND"}
            : "NA";

        my $inslen =
            exists $INFO{"INSLEN"}
            ? $INFO{"INSLEN"}
            : "NA";

        my $homlen =
            exists $INFO{"HOMLEN"}
            ? $INFO{"HOMLEN"}
            : "NA";

        my $ce =
            exists $INFO{"CE"}
            ? $INFO{"CE"}
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

        my $gq =
            exists $FORMAT{"GQ"}
            ? $FORMAT{"GQ"}
            : "NA";

        my $ft =
            exists $FORMAT{"FT"}
            ? $FORMAT{"FT"}
            : "NA";

        my $dr =
            exists $FORMAT{"DR"}
            ? $FORMAT{"DR"}
            : "NA";

        my $dv =
            exists $FORMAT{"DV"}
            ? $FORMAT{"DV"}
            : "NA";

        my $rr =
            exists $FORMAT{"RR"}
            ? $FORMAT{"RR"}
            : "NA";

        my $rv =
            exists $FORMAT{"RV"}
            ? $FORMAT{"RV"}
            : "NA";


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

            # Basic VCF fields
            $id,
            $qual,
            $filter,

            # DELLY INFO
            $precise,
            $pe,
            $sr,
            $mapq,
            $srmapq,
            $ct,
            $cipos,
            $ciend,
            $inslen,
            $homlen,
            $ce,

            # FORMAT
            $gt,
            $gq,
            $ft,
            $dr,
            $dv,
            $rr,
            $rv,

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
    # Empty VCF
    # ========================================================

    if ($file_record_count == 0) {

        $empty_vcf++;


        # ----------------------------------------------------
        # If SampleID cannot be obtained from VCF header,
        # recover it from filename.
        #
        # Example:
        #
        # HCM_TRIM63_DEL_50_100_rep1.delly.vcf
        #
        # ----------------------------------------------------

        if (!defined $sample_id || $sample_id eq "") {

            if ($file =~ /\/([^\/]+)\.delly\.vcf$/) {

                $sample_id = $1;
            }
            else {

                $sample_id = "NA";
            }
        }


        # ----------------------------------------------------
        # Parse sample metadata
        #
        # Example:
        #
        # HCM_TRIM63_DEL_50_100_rep1
        #
        # Gene       = TRIM63
        # Length_Bin = 50_100
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


        # ----------------------------------------------------
        # Terminal warning
        # ----------------------------------------------------

        warn "WARNING: Empty VCF: "
            . "$sample_id\n";
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
print "DUP records         : $dup_count\n";
print "INV records         : $inv_count\n";
print "INS records         : $ins_count\n";
print "BND/TRA records     : $bnd_count\n";
print "Other records       : $other_count\n";
print "------------------------------------------------------------\n";
print "Merged output       : $output\n";
print "Empty VCF output    : $empty_output\n";
print "============================================================\n";
print "\n";

