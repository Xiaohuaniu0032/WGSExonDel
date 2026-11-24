#!/usr/bin/env perl
# -*- coding: utf-8 -*-
use strict;
use warnings;
use Getopt::Long qw(GetOptions Configure);
Configure(qw(no_ignore_case no_auto_abbrev));
use File::Basename qw(basename dirname);
use File::Path qw(make_path);
use File::Temp qw(tempdir);
use Cwd qw(abs_path);

# ============================================================
# RunSampleCUSUMDepth.pl
#
# Function:
#   Sample-level wrapper for CUSUM depth deletion calling.
#
# Input:
#   gene_depth/*.depth.tsv
#   gene_mean_depth.tsv (metadata mapping for each per-gene depth file)
#
# For each per-gene depth file:
#   call CUSUMDepthDel.pl
#
# Output:
#   One sample-level CUSUM event table. Per-gene CUSUM files are
#   generated in a temporary directory and removed automatically.
#
# Output format of merged file:
#   SampleID Gene Classification Transcript Chrom Start End Length
#   N_Bins Baseline Mean_Depth Mean_Ratio Min_Ratio Max_CUSUM Type Depth_File
#
# Notes:
#   1. This script does not define algorithm defaults.
#   2. All CUSUM parameters must be explicitly provided in config.
#   3. CUSUM_BASELINE can be "auto" or a numeric value.
#   4. Coordinates are inherited from CUSUMDepthDel.pl.
# ============================================================

my %opt;

GetOptions(
    "Config=s" => \$opt{config},
    "InDir=s"       => \$opt{in_dir},
    "GeneSummary=s" => \$opt{gene_summary},
    "Out=s"         => \$opt{out},
    "Perl=s"   => \$opt{perl},
    "Script=s" => \$opt{script},
    "Help"     => \$opt{help},
) or die usage();

if ($opt{help}) {
    print usage();
    exit 0;
}

for my $k (qw/config in_dir gene_summary out perl script/) {
    die "[ERROR] Missing required option --" . option_name($k) . "\n" . usage()
        unless defined $opt{$k} && $opt{$k} ne "";
}

die "[ERROR] Config file not found: $opt{config}\n"
    unless -s $opt{config};

die "[ERROR] Input gene depth directory not found: $opt{in_dir}\n"
    unless -d $opt{in_dir};

die "[ERROR] Gene-depth summary not found: $opt{gene_summary}\n"
    unless -s $opt{gene_summary};

die "[ERROR] CUSUM script not found: $opt{script}\n"
    unless -s $opt{script};

check_executable($opt{perl}, "perl");

$opt{config} = abs_path($opt{config});
$opt{in_dir} = abs_path($opt{in_dir});
$opt{gene_summary} = abs_path($opt{gene_summary});
$opt{script} = abs_path($opt{script});

my $out_parent = dirname($opt{out});
make_path($out_parent) if defined $out_parent && $out_parent ne "." && !-d $out_parent;

my %conf = read_config($opt{config});

my $baseline     = get_conf_required(\%conf, "CUSUM_BASELINE");
my $bin_size     = get_conf_required(\%conf, "CUSUM_BIN_SIZE");
my $k            = get_conf_required(\%conf, "CUSUM_K");
my $h            = get_conf_required(\%conf, "CUSUM_H");
my $del_ratio    = get_conf_required(\%conf, "CUSUM_DEL_RATIO");
my $min_bins     = get_conf_required(\%conf, "CUSUM_MIN_BINS");
my $min_len      = get_conf_required(\%conf, "CUSUM_MIN_LEN");
my $edge_ratio   = get_conf_required(\%conf, "CUSUM_EDGE_RATIO");
my $recover_ratio = get_conf_required(\%conf, "CUSUM_RECOVER_RATIO");
my $recover_bins = get_conf_required(\%conf, "CUSUM_RECOVER_BINS");

validate_baseline("CUSUM_BASELINE", $baseline);
validate_positive_integer("CUSUM_BIN_SIZE", $bin_size);
validate_nonnegative_number("CUSUM_K", $k);
validate_positive_number("CUSUM_H", $h);
validate_fraction("CUSUM_DEL_RATIO", $del_ratio);
validate_positive_integer("CUSUM_MIN_BINS", $min_bins);
validate_positive_integer("CUSUM_MIN_LEN", $min_len);
validate_fraction("CUSUM_EDGE_RATIO", $edge_ratio);
validate_fraction("CUSUM_RECOVER_RATIO", $recover_ratio);
validate_positive_integer("CUSUM_RECOVER_BINS", $recover_bins);

log_msg("Sample-level CUSUM depth calling started");
log_msg("Config             : $opt{config}");
log_msg("Input depth dir    : $opt{in_dir}");
log_msg("Gene summary       : $opt{gene_summary}");
log_msg("Merged output      : $opt{out}");
log_msg("Perl               : $opt{perl}");
log_msg("CUSUM script       : $opt{script}");
log_msg("CUSUM_BASELINE     : $baseline");
log_msg("CUSUM_BIN_SIZE     : $bin_size");
log_msg("CUSUM_K            : $k");
log_msg("CUSUM_H            : $h");
log_msg("CUSUM_DEL_RATIO    : $del_ratio");
log_msg("CUSUM_MIN_BINS     : $min_bins");
log_msg("CUSUM_MIN_LEN      : $min_len");
log_msg("CUSUM_EDGE_RATIO   : $edge_ratio");
log_msg("CUSUM_RECOVER_RATIO: $recover_ratio");
log_msg("CUSUM_RECOVER_BINS : $recover_bins");

my %gene_metadata = read_gene_summary($opt{gene_summary});

my @depth_files = sort glob("$opt{in_dir}/*.depth.tsv");

die "[ERROR] No *.depth.tsv files found in $opt{in_dir}\n"
    unless @depth_files;

log_msg("Input depth files found: " . scalar(@depth_files));

my $temp_dir = tempdir("HCMExonDel_CUSUM_XXXXXX", TMPDIR => 1, CLEANUP => 1);
log_msg("Temporary CUSUM directory: $temp_dir");

my @cusum_files;

for my $depth_file (@depth_files) {
    my $base = basename($depth_file);

    my $out_base = $base;
    $out_base =~ s/\.depth\.tsv$/.cusum.tsv/;

    if ($out_base eq $base) {
        $out_base .= ".cusum.tsv";
    }

    my $cusum_out = "$temp_dir/$out_base";

    log_msg("Running CUSUM: $base");

    my @cmd = (
        $opt{perl},
        $opt{script},
        "--In",           $depth_file,
        "--Out",          $cusum_out,
        "--BinSize",      $bin_size,
        "--Baseline",     $baseline,
        "--K",            $k,
        "--H",            $h,
        "--EdgeRatio",    $edge_ratio,
        "--RecoverRatio", $recover_ratio,
        "--RecoverBins",  $recover_bins,
        "--DelRatio",     $del_ratio,
        "--MinBins",      $min_bins,
        "--MinLen",       $min_len,
    );

    run_cmd(@cmd);

    die "[ERROR] CUSUM output not generated: $cusum_out\n"
        unless -e $cusum_out;

    push @cusum_files, $cusum_out;
}

merge_cusum_files(\@cusum_files, $opt{out}, \%gene_metadata);

log_msg("Per-gene CUSUM files processed: " . scalar(@cusum_files));
log_msg("Merged CUSUM output written   : $opt{out}");
log_msg("Sample-level CUSUM depth calling finished");

exit 0;

# ============================================================
# Read gene-depth metadata
# ============================================================

sub read_gene_summary {
    my ($file) = @_;

    open my $IN, "<", $file
        or die "[ERROR] Cannot open gene-depth summary $file: $!\n";

    my $header = <$IN>;
    die "[ERROR] Empty gene-depth summary: $file\n"
        unless defined $header;
    chomp $header;
    $header =~ s/\r$//;
    my @columns = split /\t/, $header, -1;
    my %index = map { $columns[$_] => $_ } 0 .. $#columns;

    for my $required (qw/SampleID Gene Classification Transcript Depth_File/) {
        die "[ERROR] Required column '$required' missing from $file\n"
            unless exists $index{$required};
    }

    my %metadata;
    while (my $line = <$IN>) {
        chomp $line;
        $line =~ s/\r$//;
        next if $line eq "";
        my @fields = split /\t/, $line, -1;
        die "[ERROR] Column count mismatch in $file: $line\n"
            unless @fields == @columns;
        my $depth_file = $fields[$index{Depth_File}];
        die "[ERROR] Empty Depth_File in $file\n"
            unless defined $depth_file && $depth_file ne "";
        die "[ERROR] Duplicate Depth_File '$depth_file' in $file\n"
            if exists $metadata{$depth_file};
        $metadata{$depth_file} = {
            SampleID      => $fields[$index{SampleID}],
            Gene          => $fields[$index{Gene}],
            Classification => $fields[$index{Classification}],
            Transcript     => $fields[$index{Transcript}],
        };
    }
    close $IN;

    die "[ERROR] No gene metadata loaded from $file\n"
        unless keys %metadata;

    return %metadata;
}

# ============================================================
# Merge per-gene CUSUM outputs
# ============================================================

sub merge_cusum_files {
    my ($files_ref, $out_file, $metadata_ref) = @_;

    open my $OUT, ">", $out_file
        or die "[ERROR] Cannot write merged CUSUM output $out_file: $!\n";

    print $OUT join(
        "\t",
        qw/
          SampleID
          Gene
          Classification
          Transcript
          Chrom
          Start
          End
          Length
          N_Bins
          Baseline
          Mean_Depth
          Mean_Ratio
          Min_Ratio
          Max_CUSUM
          Type
          Depth_File
        /
    ), "\n";

    my $total_records = 0;

    for my $file (@$files_ref) {
        open my $IN, "<", $file
            or die "[ERROR] Cannot open per-gene CUSUM output $file: $!\n";

        my $header = <$IN>;

        if (!defined $header) {
            close $IN;
            next;
        }

        chomp $header;
        $header =~ s/\r$//;

        my @cols = split /\t/, $header, -1;

        my $expected = join(
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
        );

        my $observed = join("\t", @cols);

        die "[ERROR] Unexpected CUSUM header in $file\n"
          . "[ERROR] Expected: $expected\n"
          . "[ERROR] Observed: $observed\n"
            unless $observed eq $expected;

        my $depth_file_name = basename($file);
        $depth_file_name =~ s/\.cusum\.tsv$/.depth.tsv/;
        my $meta = $metadata_ref->{$depth_file_name};
        die "[ERROR] No gene-summary metadata for depth file: $depth_file_name\n"
            unless $meta;

        while (my $line = <$IN>) {
            chomp $line;
            $line =~ s/\r$//;

            next if $line =~ /^\s*$/;

            my @fields = split /\t/, $line, -1;
            die "[ERROR] Unexpected CUSUM record in $file: $line\n"
                unless @fields == 11;

            print $OUT join("\t",
                $meta->{SampleID},
                $meta->{Gene},
                $meta->{Classification},
                $meta->{Transcript},
                @fields[0 .. 10],
                $depth_file_name,
            ), "\n";
            $total_records++;
        }

        close $IN;
    }

    close $OUT;

    log_msg("Merged CUSUM records: $total_records");
}

# ============================================================
# Config and validation
# ============================================================

sub read_config {
    my ($file) = @_;

    open my $fh, "<", $file
        or die "[ERROR] Cannot open config file $file: $!\n";

    my %conf;

    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;

        next if $line =~ /^\s*$/;
        next if $line =~ /^\s*#/;

        $line =~ s/\s+#.*$//;

        next unless $line =~ /^\s*([^=\s]+)\s*=\s*(.*?)\s*$/;

        my $key = $1;
        my $val = $2;

        $val =~ s/^\s+|\s+$//g;
        $val =~ s/^["']//;
        $val =~ s/["']$//;

        $conf{$key} = $val;
    }

    close $fh;

    return %conf;
}

sub get_conf_required {
    my ($conf_ref, $key) = @_;

    die "[ERROR] Required config key '$key' is missing\n"
        unless exists $conf_ref->{$key}
        && defined $conf_ref->{$key}
        && $conf_ref->{$key} ne "";

    return $conf_ref->{$key};
}

sub validate_baseline {
    my ($name, $v) = @_;

    return 1 if defined $v && $v eq "auto";

    validate_positive_number($name, $v);

    return 1;
}

sub validate_fraction {
    my ($name, $v) = @_;

    die "[ERROR] $name must be numeric: $v\n"
        unless defined $v
        && $v =~ /^(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;

    die "[ERROR] $name must be > 0 and <= 1: $v\n"
        unless $v > 0 && $v <= 1;

    return 1;
}

sub validate_positive_integer {
    my ($name, $v) = @_;

    die "[ERROR] $name must be a positive integer: $v\n"
        unless defined $v && $v =~ /^\d+$/ && $v >= 1;

    return 1;
}

sub validate_positive_number {
    my ($name, $v) = @_;

    die "[ERROR] $name must be numeric: $v\n"
        unless defined $v
        && $v =~ /^(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;

    die "[ERROR] $name must be > 0: $v\n"
        unless $v > 0;

    return 1;
}

sub validate_nonnegative_number {
    my ($name, $v) = @_;

    die "[ERROR] $name must be numeric: $v\n"
        unless defined $v
        && $v =~ /^(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;

    die "[ERROR] $name must be >= 0: $v\n"
        unless $v >= 0;

    return 1;
}

# ============================================================
# General helpers
# ============================================================

sub check_executable {
    my ($cmd, $name) = @_;

    if ($cmd =~ /\//) {
        die "[ERROR] $name executable not found or not executable: $cmd\n"
            unless -x $cmd;
    }
    else {
        my $ret = system("command -v $cmd >/dev/null 2>&1");
        die "[ERROR] $name executable not found in PATH: $cmd\n"
            if $ret != 0;
    }

    return 1;
}

sub run_cmd {
    my (@cmd) = @_;

    print STDERR "[CMD] ", join(" ", map { shell_quote($_) } @cmd), "\n";

    system(@cmd) == 0
        or die "[ERROR] Command failed: " . join(" ", @cmd) . "\n";
}

sub shell_quote {
    my ($s) = @_;

    return "''" if !defined $s || $s eq "";

    if ($s =~ /^[A-Za-z0-9_\.\-\/:=]+$/) {
        return $s;
    }

    $s =~ s/'/'\\''/g;

    return "'$s'";
}

sub option_name {
    my ($key) = @_;

    my %name = (
        config      => "Config",
        in_dir      => "InDir",
        out_dir     => "OutDir",
        out         => "Out",
        perl        => "Perl",
        script      => "Script",
        ratio_all   => "RatioAll",
        summary_all => "SummaryAll",
    );

    return exists $name{$key} ? $name{$key} : $key;
}

sub log_msg {
    my ($msg) = @_;

    print STDERR "[INFO] $msg\n";
}

sub usage {
    return <<'USAGE';
Usage:
  perl bin/RunSampleCUSUMDepth.pl \
    --Config /absolute/path/hcm_exondel.conf \
    --InDir /absolute/path/S2_Depth/gene_depth \
    --GeneSummary /absolute/path/S2_Depth/gene_mean_depth.tsv \
    --Out /absolute/path/S2_Depth/cusum_events.tsv \
    --Perl /usr/bin/perl \
    --Script /absolute/path/bin/CUSUMDepthDel.pl

Required arguments:
  --Config        HCMExonDel configuration file.
  --InDir         Directory containing per-gene base-depth files.
  --GeneSummary   Gene-level summary containing Depth_File metadata.
  --Out           Combined sample-level CUSUM event table.
  --Perl          Perl executable.
  --Script        CUSUMDepthDel.pl path.

Per-gene CUSUM outputs are temporary and are removed after merging.
USAGE
}
