#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(GetOptions Configure);
Configure(qw(no_ignore_case no_auto_abbrev));

# ============================================================
# CUSUM-based read-depth deletion caller
#
# Input:
# chr pos depth
#
# Example:
# chr20 42778748 53
# chr20 42778749 53
#
# Output:
# CHROM START END LENGTH N_BINS BASELINE MEAN_DEPTH
# MEAN_RATIO MIN_RATIO MAX_CUSUM TYPE
#
# Algorithm:
# 1. Bin depth by --BinSize
# 2. Estimate baseline depth
# 3. Calculate depth ratio = depth / baseline
# 4. Apply one-sided CUSUM for downward shift:
#    S_i = max(0, S_{i-1} + (1 - ratio_i - k))
# 5. If S_i >= h, trigger candidate DEL
# 6. End DEL when ratio recovers for --RecoverBins
#    consecutive bins
# 7. Filter by min length, min bins, and mean ratio
#
# Special handling:
# If the input contains valid depth records but all depth values
# are exactly zero, write a header-only output and exit normally.
# ============================================================

my $in;
my $out;

my $bin_size    = 1;
my $baseline_arg = "auto";

my $k = 0.10;              # allowed normal fluctuation
my $h = 1.00;              # CUSUM alarm threshold

my $edge_ratio    = 0.80;  # used to refine DEL start
my $recover_ratio = 0.80;  # ratio >= this value means recovery
my $recover_bins  = 3;     # consecutive recovered bins to close event

my $del_ratio = 0.65;      # final mean ratio threshold for DEL
my $min_bins  = 3;
my $min_len   = 1;

my $pseudo = 0.1;

GetOptions(
    "In=s"            => \$in,
    "Out=s"           => \$out,
    "BinSize=i"       => \$bin_size,
    "Baseline=s"      => \$baseline_arg,
    "K=f"             => \$k,
    "H=f"             => \$h,
    "EdgeRatio=f"     => \$edge_ratio,
    "RecoverRatio=f"  => \$recover_ratio,
    "RecoverBins=i"   => \$recover_bins,
    "DelRatio=f"      => \$del_ratio,
    "MinBins=i"       => \$min_bins,
    "MinLen=i"        => \$min_len,
) or die usage();

die usage()
    unless defined $in
    && defined $out
    && $in ne ""
    && $out ne "";

die "ERROR: --BinSize must be >= 1\n"
    unless $bin_size >= 1;

die "ERROR: --RecoverBins must be >= 1\n"
    unless $recover_bins >= 1;

die "ERROR: --MinBins must be >= 1\n"
    unless $min_bins >= 1;

# ------------------------------------------------------------
# Step 1. Read depth and bin data
# ------------------------------------------------------------

my %data;
my @chr_order;

open my $IN, "<", $in
    or die "Cannot open $in: $!\n";

my ($cur_chr, $cur_bin_start, $cur_bin_end);
my ($sum_depth, $count_depth) = (0, 0);

# Count original valid depth records.
#
# These counters are used only to determine whether the input
# contains valid records but every depth value is exactly zero.
my $total_depth_records   = 0;
my $nonzero_depth_records = 0;

sub flush_bin {
    return unless defined $cur_chr;

    my $mean_depth = $count_depth > 0
        ? $sum_depth / $count_depth
        : 0;

    if (!exists $data{$cur_chr}) {
        push @chr_order, $cur_chr;
        $data{$cur_chr} = [];
    }

    push @{$data{$cur_chr}}, {
        start => $cur_bin_start,
        end   => $cur_bin_end,
        depth => $mean_depth,
        npos  => $count_depth,
    };
}

while (my $line = <$IN>) {
    chomp $line;

    next if $line =~ /^\s*$/;
    next if $line =~ /^\s*#/;

    my @f = split /\s+/, $line;

    die "ERROR: input line should have at least 3 columns: "
      . "$line\n"
        if @f < 3;

    my ($chr, $pos, $depth) = @f[0, 1, 2];

    die "ERROR: invalid chromosome: $line\n"
        unless defined $chr && $chr ne "";

    die "ERROR: invalid position: $line\n"
        unless $pos =~ /^\d+$/;

    die "ERROR: position must be >= 1: $line\n"
        unless $pos >= 1;

    die "ERROR: invalid depth: $line\n"
        unless $depth =~ /^-?\d+(?:\.\d+)?$/;

    # Record-level statistics used for the all-zero check.
    $total_depth_records++;

    if ($depth != 0) {
        $nonzero_depth_records++;
    }

    my $bin_start =
        int(($pos - 1) / $bin_size) * $bin_size + 1;

    my $bin_end =
        $bin_start + $bin_size - 1;

    if (!defined $cur_chr) {
        (
            $cur_chr,
            $cur_bin_start,
            $cur_bin_end
        ) = (
            $chr,
            $bin_start,
            $bin_end
        );

        ($sum_depth, $count_depth) = ($depth, 1);
    }
    elsif (
        $chr eq $cur_chr
        && $bin_start == $cur_bin_start
    ) {
        $sum_depth += $depth;
        $count_depth++;
    }
    else {
        flush_bin();

        (
            $cur_chr,
            $cur_bin_start,
            $cur_bin_end
        ) = (
            $chr,
            $bin_start,
            $bin_end
        );

        ($sum_depth, $count_depth) = ($depth, 1);
    }
}

flush_bin();

close $IN
    or die "Cannot close $in: $!\n";

# ------------------------------------------------------------
# Handle input files in which all depth values are zero
# ------------------------------------------------------------

if (
    $total_depth_records > 0
    && $nonzero_depth_records == 0
) {
    open my $OUT, ">", $out
        or die "Cannot write $out: $!\n";

    print $OUT join(
        "\t",
        qw/
            CHROM
            START
            END
            LENGTH
            N_BINS
            BASELINE
            MEAN_DEPTH
            MEAN_RATIO
            MIN_RATIO
            MAX_CUSUM
            TYPE
        /
    ), "\n";

    close $OUT
        or die "Cannot close $out: $!\n";

    print STDERR
        "[WARN] All depth values are zero in $in; "
      . "CUSUM analysis is skipped.\n";

    print STDERR
        "[DONE] Header-only output written to $out\n";

    exit 0;
}

# ------------------------------------------------------------
# Step 2. Estimate baseline depth
# ------------------------------------------------------------

my %baseline;

if ($baseline_arg eq "auto") {
    for my $chr (@chr_order) {
        my @depths = map {
            $_->{depth}
        }
        grep {
            $_->{depth} > 0
        }
        @{$data{$chr}};

        die "ERROR: no valid depth for chromosome $chr\n"
            unless @depths;

        $baseline{$chr} = median(\@depths);
    }
}
else {
    die "ERROR: --Baseline should be auto or numeric value\n"
        unless $baseline_arg =~ /^\d+(?:\.\d+)?$/;

    for my $chr (@chr_order) {
        $baseline{$chr} = $baseline_arg;
    }
}

for my $chr (@chr_order) {
    die "ERROR: baseline for $chr must be > 0\n"
        unless defined $baseline{$chr}
        && $baseline{$chr} > 0;
}

# ------------------------------------------------------------
# Step 3. Calculate ratio
# ------------------------------------------------------------

for my $chr (@chr_order) {
    my $base = $baseline{$chr};

    for my $b (@{$data{$chr}}) {
        my $ratio =
            ($b->{depth} + $pseudo)
            /
            ($base + $pseudo);

        $b->{ratio} = $ratio;
    }
}

# ------------------------------------------------------------
# Step 4. CUSUM scan for DEL regions
# ------------------------------------------------------------

open my $OUT, ">", $out
    or die "Cannot write $out: $!\n";

print $OUT join(
    "\t",
    qw/
        CHROM
        START
        END
        LENGTH
        N_BINS
        BASELINE
        MEAN_DEPTH
        MEAN_RATIO
        MIN_RATIO
        MAX_CUSUM
        TYPE
    /
), "\n";

for my $chr (@chr_order) {
    my $bins = $data{$chr};
    my $n    = scalar @{$bins};

    next if $n == 0;

    my $base = $baseline{$chr};

    my $S = 0;

    # Possible DEL start before the alarm threshold is reached.
    my $candidate_start_i;

    # Whether the current position is inside a DEL event.
    my $in_event = 0;

    my $event_start_i;
    my $trigger_i;
    my $max_cusum    = 0;
    my $recover_count = 0;

    for (my $i = 0; $i < $n; $i++) {
        my $ratio = $bins->[$i]{ratio};

        # Downward-shift CUSUM increment.
        #
        # If the depth ratio is substantially lower than 1,
        # the increment becomes positive.
        my $inc = 1 - $ratio - $k;

        # Record a possible event start when CUSUM begins
        # accumulating.
        if (
            !$in_event
            && $S == 0
            && $inc > 0
        ) {
            $candidate_start_i = $i;
        }

        $S += $inc;

        if ($S < 0) {
            $S = 0;
        }

        # Trigger a DEL event when CUSUM reaches the alarm
        # threshold.
        if (
            !$in_event
            && $S >= $h
        ) {
            $in_event = 1;
            $trigger_i = $i;

            # Refine the event start by locating the first bin
            # between the candidate start and trigger position
            # whose ratio is <= edge_ratio.
            $event_start_i = refine_start(
                $bins,
                defined $candidate_start_i
                    ? $candidate_start_i
                    : $i,
                $trigger_i,
                $edge_ratio
            );

            $max_cusum     = $S;
            $recover_count = 0;
        }

        if ($in_event) {
            if ($S > $max_cusum) {
                $max_cusum = $S;
            }

            # Close the event after the depth ratio has recovered
            # for the required number of consecutive bins.
            if ($ratio >= $recover_ratio) {
                $recover_count++;
            }
            else {
                $recover_count = 0;
            }

            if ($recover_count >= $recover_bins) {
                my $event_end_i =
                    $i - $recover_count;

                output_event(
                    $OUT,
                    $chr,
                    $bins,
                    $event_start_i,
                    $event_end_i,
                    $base,
                    $max_cusum,
                    $del_ratio,
                    $min_bins,
                    $min_len
                );

                # Reset status after closing the event.
                $S                 = 0;
                $candidate_start_i = undef;
                $in_event          = 0;
                $event_start_i     = undef;
                $trigger_i         = undef;
                $max_cusum         = 0;
                $recover_count     = 0;
            }
        }

        # If CUSUM returns to zero before an alarm is triggered,
        # clear the candidate start.
        if (
            !$in_event
            && $S == 0
        ) {
            $candidate_start_i = undef;
        }
    }

    # If an event reaches the chromosome or input-file end,
    # evaluate and output the remaining event.
    if ($in_event) {
        my $event_end_i = $n - 1;

        output_event(
            $OUT,
            $chr,
            $bins,
            $event_start_i,
            $event_end_i,
            $base,
            $max_cusum,
            $del_ratio,
            $min_bins,
            $min_len
        );
    }
}

close $OUT
    or die "Cannot close $out: $!\n";

print STDERR "[DONE] Output written to $out\n";

# ============================================================
# Subroutines
# ============================================================

sub output_event {
    my (
        $OUT,
        $chr,
        $bins,
        $start_i,
        $end_i,
        $base,
        $max_cusum,
        $del_ratio,
        $min_bins,
        $min_len
    ) = @_;

    return unless defined $start_i;
    return unless defined $end_i;
    return if $end_i < $start_i;

    my $start = $bins->[$start_i]{start};
    my $end   = $bins->[$end_i]{end};

    my $n_bins =
        $end_i - $start_i + 1;

    my $len =
        $end - $start + 1;

    return if $n_bins < $min_bins;
    return if $len < $min_len;

    my $sum_depth = 0;
    my $sum_ratio = 0;
    my $min_ratio;

    for my $i ($start_i .. $end_i) {
        $sum_depth += $bins->[$i]{depth};
        $sum_ratio += $bins->[$i]{ratio};

        if (
            !defined $min_ratio
            || $bins->[$i]{ratio} < $min_ratio
        ) {
            $min_ratio = $bins->[$i]{ratio};
        }
    }

    my $mean_depth =
        $sum_depth / $n_bins;

    my $mean_ratio =
        $sum_ratio / $n_bins;

    # Final deletion filter.
    return if $mean_ratio > $del_ratio;

    print $OUT join(
        "\t",
        $chr,
        $start,
        $end,
        $len,
        $n_bins,
        sprintf("%.4f", $base),
        sprintf("%.4f", $mean_depth),
        sprintf("%.4f", $mean_ratio),
        sprintf("%.4f", $min_ratio),
        sprintf("%.4f", $max_cusum),
        "DEL"
    ), "\n";
}

sub refine_start {
    my (
        $bins,
        $from_i,
        $to_i,
        $edge_ratio
    ) = @_;

    for my $i ($from_i .. $to_i) {
        if ($bins->[$i]{ratio} <= $edge_ratio) {
            return $i;
        }
    }

    return $from_i;
}

sub median {
    my ($arr) = @_;

    my @v = sort {
        $a <=> $b
    } @{$arr};

    my $n = scalar @v;

    return $v[int($n / 2)]
        if $n % 2 == 1;

    return (
        $v[$n / 2 - 1]
        +
        $v[$n / 2]
    ) / 2;
}

sub usage {
    return <<"USAGE";
Usage:
  perl CUSUMDepthDel.pl --In depth.txt --Out del.cusum.tsv [options]

Required:
  --In              Input depth file, format: chr pos depth
  --Out             Output DEL result file

Options:
  --BinSize        Bin size in bp, default: 1
                    1 means per-base depth.

  --Baseline        auto or numeric value, default: auto
                    auto means median depth per chromosome.
                    For small local regions, numeric baseline
                    is recommended.

  --K               CUSUM tolerance parameter, default: 0.10
                    Larger k makes the algorithm more
                    conservative.

  --H               CUSUM alarm threshold, default: 1.00
                    Larger h requires stronger continuous
                    depth decrease.

  --EdgeRatio      Ratio used to refine DEL start,
                    default: 0.80

  --RecoverRatio   Ratio used to define recovery,
                    default: 0.80

  --RecoverBins    Consecutive recovered bins required to
                    close DEL, default: 3

  --DelRatio       Final mean ratio threshold for DEL,
                    default: 0.65

  --MinBins        Minimum bins required for DEL,
                    default: 3

  --MinLen         Minimum DEL length in bp,
                    default: 1

Examples:
  perl CUSUMDepthDel.pl \\
    --In depth.txt \\
    --Out del.cusum.tsv \\
    --Baseline 50 \\
    --BinSize 1 \\
    --K 0.10 \\
    --H 1.00 \\
    --DelRatio 0.65 \\
    --MinBins 3

  perl CUSUMDepthDel.pl \\
    --In depth.txt \\
    --Out del.cusum.tsv \\
    --Baseline auto \\
    --BinSize 10 \\
    --K 0.10 \\
    --H 1.50 \\
    --RecoverRatio 0.80 \\
    --RecoverBins 2 \\
    --MinLen 50
USAGE
}

