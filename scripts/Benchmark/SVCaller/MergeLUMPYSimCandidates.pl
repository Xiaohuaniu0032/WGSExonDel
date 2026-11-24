#!/usr/bin/env perl
use strict;
use warnings;
use Cwd qw(abs_path);

# ============================================================
# Script:
#   MergeLUMPYSimCandidates.pl
#
# Description:
#   Merge LUMPY VCF results from simulated deletion samples
#   into one tab-delimited XLS file.
#
#   Empty VCF samples are recorded separately for downstream
#   false-negative / sensitivity analysis.
#
# Working directory:
#   The caller-results root directory supplied as the first argument.
#
# Input:
#   LUMPY/*/*/*.lumpy.vcf
#
# Output:
#   1. LUMPY_sim_merged_candidates.xls
#   2. LUMPY_sim_empty_vcf_samples.xls
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
#   For non-DEL SV types:
#
#       Start = VCF_POS
#       End   = VCF_END
#       Size  = abs(SVLEN)
#
# ============================================================


# ============================================================
# Configuration
# ============================================================

my $sv_caller_dir = shift @ARGV || ".";
$sv_caller_dir = abs_path($sv_caller_dir) || $sv_caller_dir;

my $input_dir =
    "$sv_caller_dir/LUMPY";

my $output =
    "$sv_caller_dir/LUMPY_sim_merged_candidates.xls";

my $empty_output =
    "$sv_caller_dir/LUMPY_sim_empty_vcf_samples.xls";


# ============================================================
# Find input VCF files
# ============================================================

my @files = sort glob(
    "$input_dir/*/*/*.lumpy.vcf"
);

die "ERROR: No LUMPY VCF files found:\n"
    . "$input_dir/*/*/*.lumpy.vcf\n"
    unless @files;


# ============================================================
# Open outputs
# ============================================================

open(my $OUT, ">", $output)
    or die "ERROR: Cannot write $output: $!\n";

open(my $EMPTY, ">", $empty_output)
    or die "ERROR: Cannot write $empty_output: $!\n";


# ============================================================
# Main output header
#
# Keep consistent with real-data LUMPY merged result
# ============================================================

my @header = (

    # Standardized fields
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

    # LUMPY INFO fields
    "Precision",
    "STRANDS",
    "SU",
    "PE",
    "SR",
    "CIPOS",
    "CIEND",
    "CIPOS95",
    "CIEND95",

    # FORMAT/sample fields
    "GT",
    "FORMAT_SU",
    "FORMAT_PE",
    "FORMAT_SR",

    # Original VCF fields
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

my $file_count   = 0;
my $record_count = 0;
my $empty_vcf    = 0;

my $del_count    = 0;
my $dup_count    = 0;
my $inv_count    = 0;
my $bnd_count    = 0;
my $ins_count    = 0;
my $other_count  = 0;


# ============================================================
# Run information
# ============================================================

print "\n";
print "============================================================\n";
print " LUMPY simulated-data merger\n";
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

        # Skip empty lines
        next if $line =~ /^\s*$/;


        # ----------------------------------------------------
        # Check SampleID
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
                # IMPRECISE / PRECISE
                $INFO{$item} = 1;
            }
        }


        # ====================================================
        # Main SV information
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
        #
        # Standardized deletion interval:
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

            # Fallback to coordinates
            elsif (
                $vcf_pos =~ /^\d+$/
                &&
                $vcf_end =~ /^\d+$/
            ) {

                $size = $vcf_end - $vcf_pos;
            }
        }


        # ----------------------------------------------------
        # Non-DEL records
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
            # Statistics
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
        # Precision
        # ====================================================

        my $precision;

        if (exists $INFO{"IMPRECISE"}) {

            $precision = "Imprecise";
        }

        elsif (exists $INFO{"PRECISE"}) {

            $precision = "Precise";
        }

        else {

            $precision = "NA";
        }


        # ====================================================
        # Extract LUMPY INFO fields
        # ====================================================

        my $strands =
            exists $INFO{"STRANDS"}
            ? $INFO{"STRANDS"}
            : "NA";

        my $su =
            exists $INFO{"SU"}
            ? $INFO{"SU"}
            : "NA";

        my $pe =
            exists $INFO{"PE"}
            ? $INFO{"PE"}
            : "NA";

        my $sr =
            exists $INFO{"SR"}
            ? $INFO{"SR"}
            : "NA";

        my $cipos =
            exists $INFO{"CIPOS"}
            ? $INFO{"CIPOS"}
            : "NA";

        my $ciend =
            exists $INFO{"CIEND"}
            ? $INFO{"CIEND"}
            : "NA";

        my $cipos95 =
            exists $INFO{"CIPOS95"}
            ? $INFO{"CIPOS95"}
            : "NA";

        my $ciend95 =
            exists $INFO{"CIEND95"}
            ? $INFO{"CIEND95"}
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

        my $format_su =
            exists $FORMAT{"SU"}
            ? $FORMAT{"SU"}
            : "NA";

        my $format_pe =
            exists $FORMAT{"PE"}
            ? $FORMAT{"PE"}
            : "NA";

        my $format_sr =
            exists $FORMAT{"SR"}
            ? $FORMAT{"SR"}
            : "NA";


        # ====================================================
        # Output merged record
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

            # LUMPY INFO
            $precision,
            $strands,
            $su,
            $pe,
            $sr,
            $cipos,
            $ciend,
            $cipos95,
            $ciend95,

            # FORMAT
            $gt,
            $format_su,
            $format_pe,
            $format_sr,

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
        # Recover SampleID from filename if necessary
        #
        # HCM_ACTN2_DEL_1000_2000_rep1.lumpy.vcf
        #
        # ----------------------------------------------------

        if (!defined $sample_id || $sample_id eq "") {

            if ($file =~ /\/([^\/]+)\.lumpy\.vcf$/) {

                $sample_id = $1;
            }

            else {

                $sample_id = "NA";
            }
        }


        # ----------------------------------------------------
        # Parse sample metadata
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
        # Write empty VCF sample
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
# Close output files
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

