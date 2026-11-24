#!/usr/bin/env perl
# -*- coding: utf-8 -*-

use strict;
use warnings;
use Getopt::Long qw(GetOptions Configure);
Configure(qw(no_ignore_case no_auto_abbrev));
use File::Path qw(make_path);
use File::Basename qw(dirname);
use Cwd qw(abs_path);
use FindBin qw($Bin);

# ============================================================
# HCMExonDel.pl
#
# Main pipeline generator for HCM exon-level deletion detection.
#
# Function:
#   Generate one executable run shell for each BAM sample.
#
# Input:
#   1. Config file
#   2. BAM list, two TAB-delimited columns:
#        SampleID    /absolute/path/to/sample.bam
#
# Output:
#   OUTDIR/SAMPLE/run_SAMPLE.sh
#
# Workflow:
#   Step 00. extract_target_bam
#   Step 01. gene_mean_depth_and_window_ratio
#   Step 02. cusum_depth_evidence
#   Step 03. plot_depth_ratio
#   Step 04. extract_raw_sa_reads
#   Step 05. annotate_split_read_breakpoints
#   Step 06. cluster_split_reads
#   Step 07. extract_valid_tlen
#   Step 08. plot_tlen_distribution
#   Step 09. extract_discordant_pairs
#   Step 10. cluster_discordant_pairs
#   Step 11. merge_evidence
#   Step 12. annotate_candidates
#
# Design notes:
#   1. HCMExonDel only analyzes HCM core genes.
#   2. The original WGS BAM is used only by ExtractTargetBAM.pl.
#   3. All downstream analyses use SAMPLE.target.bam.
#   4. MergeEvidence.pl writes high-quality depth, split and discordant
#      candidate tables, then integrates them using split-read clusters
#      as the candidate event backbone.
#   5. Depth and discordant signals remain supporting evidence; the
#      existing candidate-level evidence rules are unchanged.
#   6. HCMExonDel.pl reads EVIDENCE_OVERLAP_FRACTION from config and
#      passes it to MergeEvidence.pl as --MinOverlap.
# ============================================================

my $config;
my $bam_list;
my $outdir = "results";
my $genome_build = "hg19";
my $help   = 0;

GetOptions(
    "Config=s"   => \$config,
    "BamList=s"  => \$bam_list,
    "OutDir=s"      => \$outdir,
    "GenomeBuild=s" => \$genome_build,
    "Help|h"        => \$help,
) or die usage();

if ($help) {
    print usage();
    exit 0;
}

die usage() unless defined $config && defined $bam_list;

$config = abs_path($config);
die "[ERROR] Config file not found: $config\n"
    unless defined $config && -s $config;

$bam_list = abs_path($bam_list);
die "[ERROR] BAM list file not found: $bam_list\n"
    unless defined $bam_list && -s $bam_list;

my $project_root = abs_path($Bin);
die "[ERROR] Cannot resolve HCMExonDel project root: $Bin\n"
    unless defined $project_root && -d $project_root;

$genome_build = lc($genome_build);
die "[ERROR] GenomeBuild must be hg19 or hg38: $genome_build\n"
    unless $genome_build eq "hg19" || $genome_build eq "hg38";

my %CONF = read_config($config);

# ------------------------------------------------------------
# Runtime executables
# ------------------------------------------------------------
my $perl   = abs_path($^X) || $^X;
my $python = get_conf_required(\%CONF, "PYTHON");
my $samtools = get_conf_required(\%CONF, "SAMTOOLS");

for my $entry (
    [PERL => $perl],
    [PYTHON => $python],
    [SAMTOOLS => $samtools],
) {
    my ($name, $path) = @$entry;
    die "[ERROR] $name must be an absolute path: $path\n"
        unless $path =~ m{^/};
}

check_executable($perl, "PERL");
check_executable($python, "PYTHON");
check_executable($samtools, "SAMTOOLS");

# ------------------------------------------------------------
# Workflow switches
# ------------------------------------------------------------
my $run_tlen_qc = normalize_bool(
    get_conf_required(\%CONF, "RUN_TLEN_QC"),
    "RUN_TLEN_QC"
);

my @required_config_keys = qw(
    PYTHON SAMTOOLS HCM_CORE_GENE_LIST
    TARGET_REGION_FLANK TARGET_BAM_THREADS
    MIN_MAPQ EXCLUDE_DUPLICATES
    WINDOW_SIZE WINDOW_STEP MIN_WINDOW_SIZE MIN_GENE_MEAN_DEPTH
    DEL_DEPTH_RATIO_CUTOFF MIN_CONSECUTIVE_DEL_WINDOWS
    CUSUM_BASELINE CUSUM_BIN_SIZE CUSUM_K CUSUM_H CUSUM_DEL_RATIO
    CUSUM_MIN_BINS CUSUM_MIN_LEN CUSUM_EDGE_RATIO CUSUM_RECOVER_RATIO CUSUM_RECOVER_BINS
    SPLIT_READ_THREADS SA_SPLIT_CLUSTER_WINDOW SA_SPLIT_MIN_SUPPORT_READS
    RUN_TLEN_QC MIN_DISCORDANT_INSERT_SIZE MIN_DISCORDANT_PAIRS
    DISCORDANT_CLUSTER_DISTANCE FILTER_DELETION_ORIENTATION DISCORDANT_READ_THREADS
    EVIDENCE_OVERLAP_FRACTION MIN_EXON_OVERLAP_FRACTION
);
push @required_config_keys, qw/TLEN_MAX_PAIRS TLEN_THREADS TLEN_RANDOM_SEED/
    if $run_tlen_qc;

for my $key (@required_config_keys) {
    get_conf_required(\%CONF, $key);
}

# ------------------------------------------------------------
# Algorithm parameters checked by main pipeline
# ------------------------------------------------------------
my $target_region_flank = get_conf_required(\%CONF, "TARGET_REGION_FLANK");
validate_nonnegative_integer("TARGET_REGION_FLANK", $target_region_flank);

my $target_bam_threads = get_conf_required(\%CONF, "TARGET_BAM_THREADS");
validate_positive_integer("TARGET_BAM_THREADS", $target_bam_threads);

my $split_read_threads = get_conf_required(\%CONF, "SPLIT_READ_THREADS");
validate_positive_integer("SPLIT_READ_THREADS", $split_read_threads);

my $discordant_read_threads = get_conf_required(\%CONF, "DISCORDANT_READ_THREADS");
validate_positive_integer("DISCORDANT_READ_THREADS", $discordant_read_threads);

my $evidence_overlap_fraction = get_conf_required(\%CONF, "EVIDENCE_OVERLAP_FRACTION");
validate_fraction("EVIDENCE_OVERLAP_FRACTION", $evidence_overlap_fraction);

my %CUSUM = (
    baseline      => get_conf_required(\%CONF, "CUSUM_BASELINE"),
    bin_size      => get_conf_required(\%CONF, "CUSUM_BIN_SIZE"),
    k             => get_conf_required(\%CONF, "CUSUM_K"),
    h             => get_conf_required(\%CONF, "CUSUM_H"),
    del_ratio     => get_conf_required(\%CONF, "CUSUM_DEL_RATIO"),
    min_bins      => get_conf_required(\%CONF, "CUSUM_MIN_BINS"),
    min_len       => get_conf_required(\%CONF, "CUSUM_MIN_LEN"),
    edge_ratio    => get_conf_required(\%CONF, "CUSUM_EDGE_RATIO"),
    recover_ratio => get_conf_required(\%CONF, "CUSUM_RECOVER_RATIO"),
    recover_bins  => get_conf_required(\%CONF, "CUSUM_RECOVER_BINS"),
);

validate_cusum_params(\%CUSUM);

my %TLEN;

if ($run_tlen_qc) {
    %TLEN = (
        max_pairs   => get_conf_required(\%CONF, "TLEN_MAX_PAIRS"),
        threads     => get_conf_required(\%CONF, "TLEN_THREADS"),
        random_seed => get_conf_required(\%CONF, "TLEN_RANDOM_SEED"),
    );

    validate_positive_integer("TLEN_MAX_PAIRS",  $TLEN{max_pairs});
    validate_positive_integer("TLEN_THREADS",    $TLEN{threads});
    validate_integer("TLEN_RANDOM_SEED",         $TLEN{random_seed});
}

# ------------------------------------------------------------
# Genome-build-specific MANE annotation
# ------------------------------------------------------------
my $annotation_dir = "$project_root/db/Annotation/$genome_build";
my $exon_txt = abs_path("$annotation_dir/MANERefSeqExons.tsv");
my $gene_txt = abs_path("$annotation_dir/MANERefSeqGenes.tsv");

my $core_gene_list = resolve_config_path(
    get_conf_required(\%CONF, "HCM_CORE_GENE_LIST")
);

die "[ERROR] MANE exon annotation not found for $genome_build: $annotation_dir/MANERefSeqExons.tsv\n"
    unless defined $exon_txt && -s $exon_txt;

die "[ERROR] MANE gene annotation not found for $genome_build: $annotation_dir/MANERefSeqGenes.tsv\n"
    unless defined $gene_txt && -s $gene_txt;

die "[ERROR] HCM_CORE_GENE_LIST file not found: $core_gene_list\n"
    unless -s $core_gene_list;

check_exon_txt_format($exon_txt);
check_gene_txt_format($gene_txt);
check_core_gene_list_format($core_gene_list);

# ------------------------------------------------------------
# Output directory
# ------------------------------------------------------------
$outdir = abs_path($outdir) || $outdir;

make_path($outdir) unless -d $outdir;

$outdir = abs_path($outdir);

# ------------------------------------------------------------
# Pipeline scripts
# ------------------------------------------------------------
my %SCRIPT = (
    target_bam        => "$Bin/bin/ExtractTargetBAM.pl",
    depth             => "$Bin/bin/RunGeneMeanDepth.pl",
    sample_cusum      => "$Bin/bin/RunSampleCUSUMDepth.pl",
    cusum_depth       => "$Bin/bin/CUSUMDepthDel.pl",
    plot_depth        => "$Bin/bin/PlotDepthRatio.py",
    split_extract     => "$Bin/bin/ExtractSASplitReads.pl",
    split_identify    => "$Bin/bin/IdentifyIntragenicDeletions.pl",
    split_cluster     => "$Bin/bin/ClusterSASplitReads.pl",
    tlen_extract      => "$Bin/bin/ExtractValidTLEN.pl",
    tlen_plot         => "$Bin/bin/PlotTLENDistribution.py",
    discordant_extract  => "$Bin/bin/ExtractDiscordantPairs.pl",
    discordant_cluster  => "$Bin/bin/ClusterDiscordantReads.pl",
    merge             => "$Bin/bin/MergeEvidence.pl",
    annotate          => "$Bin/bin/AnnotateCandidates.pl",
);

my @required_scripts = qw/
    target_bam
    depth
    sample_cusum
    cusum_depth
    plot_depth
    split_extract
    split_identify
    split_cluster
    discordant_extract
    discordant_cluster
    merge
    annotate
/;

push @required_scripts, qw/tlen_extract tlen_plot/
    if $run_tlen_qc;

for my $key (@required_scripts) {
    check_script($SCRIPT{$key}, $key);
}

# ------------------------------------------------------------
# BAM list
# ------------------------------------------------------------
my @samples = read_bam_list($bam_list);

die "[ERROR] No valid sample found in BAM list: $bam_list\n"
    unless @samples;

print_summary(
    config                      => $config,
    project_root                => $project_root,
    genome_build                => $genome_build,
    bam_list                    => $bam_list,
    exon_txt                    => $exon_txt,
    gene_txt                    => $gene_txt,
    core_gene_list              => $core_gene_list,
    outdir                      => $outdir,
    perl                        => $perl,
    python                      => $python,
    samtools                    => $samtools,
    run_tlen_qc                 => $run_tlen_qc,
    target_region_flank         => $target_region_flank,
    target_bam_threads          => $target_bam_threads,
    split_read_threads          => $split_read_threads,
    discordant_read_threads     => $discordant_read_threads,
    evidence_overlap_fraction   => $evidence_overlap_fraction,
    cusum_ref                   => \%CUSUM,
    tlen_ref                    => \%TLEN,
    sample_num                  => scalar(@samples),
);

for my $item (@samples) {
    my $sample = $item->{sample};
    my $bam    = $item->{bam};

    my %DIR = prepare_sample_dirs($outdir, $sample);
    my %OUT = build_sample_outputs($sample, \%DIR);

    my $sample_shell = "$DIR{sample}/run_${sample}.sh";

    write_sample_shell(
        shell                       => $sample_shell,
        sample                      => $sample,
        bam                         => $bam,
        config                      => $config,
        genome_build                => $genome_build,
        core_gene_list              => $core_gene_list,
        exon_txt                    => $exon_txt,
        gene_txt                    => $gene_txt,
        perl                        => $perl,
        python                      => $python,
        samtools                    => $samtools,
        conf_ref                    => \%CONF,
        run_info                    => "$DIR{sample}/run_info.tsv",
        run_tlen_qc                 => $run_tlen_qc,
        depth_ratio_cutoff          => get_conf_required(\%CONF, "DEL_DEPTH_RATIO_CUTOFF"),
        evidence_overlap_fraction   => $evidence_overlap_fraction,
        script_ref                  => \%SCRIPT,
        dir_ref                     => \%DIR,
        out_ref                     => \%OUT,
    );

    chmod 0755, $sample_shell
        or die "[ERROR] Cannot chmod sample shell: $sample_shell\n";


    print "[INFO] Sample shell generated: $sample_shell\n";
}

print "\n";
print "[INFO] All sample shells generated successfully\n\n";
print "Run example:\n";
print "  bash $outdir/SAMPLE/run_SAMPLE.sh\n\n";

exit 0;

# ============================================================
# Main shell writer
# ============================================================

sub write_sample_shell {
    my %args = @_;

    my $shell                     = $args{shell};
    my $sample                    = $args{sample};
    my $bam                       = $args{bam};
    my $config                    = $args{config};
    my $genome_build              = $args{genome_build};
    my $core_gene_list            = $args{core_gene_list};
    my $exon_txt                  = $args{exon_txt};
    my $gene_txt                  = $args{gene_txt};
    my $perl                      = $args{perl};
    my $python                    = $args{python};
    my $samtools                  = $args{samtools};
    my $conf_ref                  = $args{conf_ref};
    my $run_info                  = $args{run_info};
    my $run_tlen_qc               = $args{run_tlen_qc};
    my $depth_ratio_cutoff        = $args{depth_ratio_cutoff};
    my $evidence_overlap_fraction = $args{evidence_overlap_fraction};
    my $script_ref                = $args{script_ref};
    my $dir_ref                   = $args{dir_ref};
    my $out_ref                   = $args{out_ref};

    open my $SH, ">", $shell
        or die "[ERROR] Cannot write sample shell: $shell\n";

    write_header(
        fh             => $SH,
        sample         => $sample,
        config         => $config,
        original_bam   => $bam,
        target_bam     => $out_ref->{target_bam},
        genome_build   => $genome_build,
        exon_txt       => $exon_txt,
        gene_txt       => $gene_txt,
        core_gene_list => $core_gene_list,
        outdir         => $dir_ref->{sample},
    );

    my $mkdir_cmd = join(
        " ",
        "mkdir -p",
        map { shell_quote($_) } (
            $dir_ref->{log},
            $dir_ref->{target_bam},
            $dir_ref->{depth},
            $dir_ref->{split},
            $dir_ref->{discordant},
            $dir_ref->{candidate},
            $dir_ref->{report},
        ),
    );

    write_block($SH, "Start sample: $sample");
    print $SH "$mkdir_cmd\n\n";

    write_run_info_shell(
        fh              => $SH,
        file            => $run_info,
        sample          => $sample,
        bam             => $bam,
        config          => $config,
        genome_build    => $genome_build,
        gene_annotation => $gene_txt,
        exon_annotation => $exon_txt,
        core_gene_list  => $core_gene_list,
        perl            => $perl,
        python          => $python,
        samtools        => $samtools,
        conf_ref        => $conf_ref,
    );

    my $cmd = join(
        " ",
        shell_quote($perl),
        shell_quote($script_ref->{target_bam}),
        "--Config",         shell_quote($config),
        "--GeneAnnotation", shell_quote($gene_txt),
        "--Sample",         shell_quote($sample),
        "--Bam",            shell_quote($bam),
        "--OutDir",         shell_quote($dir_ref->{target_bam}),
    );
    write_cmd($SH, "Step 00: extract_target_bam", $cmd,
        "$dir_ref->{log}/00.extract_target_bam.log");

    $cmd = join(
        " ",
        shell_quote($perl),
        shell_quote($script_ref->{depth}),
        "--Config",         shell_quote($config),
        "--GeneAnnotation", shell_quote($gene_txt),
        "--Bam",            shell_quote($out_ref->{target_bam}),
        "--Sample",         shell_quote($sample),
        "--Out",            shell_quote($out_ref->{window_candidates}),
    );
    write_cmd($SH, "Step 01: gene_mean_depth_and_window_ratio", $cmd,
        "$dir_ref->{log}/01.depth_analysis.log");

    $cmd = join(
        " ",
        shell_quote($perl),
        shell_quote($script_ref->{sample_cusum}),
        "--Config",      shell_quote($config),
        "--InDir",       shell_quote($out_ref->{gene_depth_dir}),
        "--GeneSummary", shell_quote($out_ref->{gene_mean_depth}),
        "--Out",         shell_quote($out_ref->{cusum_events}),
        "--Perl",        shell_quote($perl),
        "--Script",      shell_quote($script_ref->{cusum_depth}),
    );
    write_cmd($SH, "Step 02: cusum_depth_evidence", $cmd,
        "$dir_ref->{log}/02.cusum_depth.log");

    $cmd = join(
        " ",
        shell_quote($python),
        shell_quote($script_ref->{plot_depth}),
        "--Input",            shell_quote($out_ref->{window_depth_ratio}),
        "--Output",           shell_quote($out_ref->{depth_plot_pdf}),
        "--DepthRatioCutoff", shell_quote($depth_ratio_cutoff),
    );
    write_cmd($SH, "Step 03: plot_depth_ratio", $cmd,
        "$dir_ref->{log}/03.depth_ratio_plot.log");

    $cmd = join(
        " ",
        shell_quote($perl),
        shell_quote($script_ref->{split_extract}),
        "--Config", shell_quote($config),
        "--Bam",    shell_quote($out_ref->{target_bam}),
        "--Sample", shell_quote($sample),
        "--Out",    shell_quote($out_ref->{split_reads}),
    );
    write_cmd($SH, "Step 04: extract_raw_split_reads", $cmd,
        "$dir_ref->{log}/04.extract_split_reads.log");

    $cmd = join(
        " ",
        shell_quote($perl),
        shell_quote($script_ref->{split_identify}),
        "--Input", shell_quote($out_ref->{split_reads}),
        "--Out",   shell_quote($out_ref->{split_reads}),
    );
    write_cmd($SH, "Step 05: annotate_split_read_breakpoints", $cmd,
        "$dir_ref->{log}/05.annotate_split_read_breakpoints.log");

    $cmd = join(
        " ",
        shell_quote($perl),
        shell_quote($script_ref->{split_cluster}),
        "--Config",      shell_quote($config),
        "--Input",       shell_quote($out_ref->{split_reads}),
        "--ClustersOut", shell_quote($out_ref->{split_clusters}),
    );
    write_cmd($SH, "Step 06: cluster_split_reads", $cmd,
        "$dir_ref->{log}/06.cluster_split_reads.log");

    if ($run_tlen_qc) {
        $cmd = join(
            " ",
            shell_quote($perl),
            shell_quote($script_ref->{tlen_extract}),
            "--Config", shell_quote($config),
            "--Bam",    shell_quote($out_ref->{target_bam}),
            "--Out",    shell_quote($out_ref->{valid_tlen_tmp}),
        );
        write_cmd($SH, "Step 07: extract_tlen", $cmd,
            "$dir_ref->{log}/07.extract_tlen.log");

        $cmd = join(
            " ",
            shell_quote($python),
            shell_quote($script_ref->{tlen_plot}),
            "--Input",   shell_quote($out_ref->{valid_tlen_tmp}),
            "--Output",  shell_quote($out_ref->{tlen_plot_pdf}),
            "--Summary", shell_quote($out_ref->{tlen_summary}),
            "--Sample",  shell_quote($sample),
        );
        write_cmd($SH, "Step 08: plot_tlen_distribution", $cmd,
            "$dir_ref->{log}/08.tlen_plot.log");
        print $SH "rm -f " . shell_quote($out_ref->{valid_tlen_tmp}) . "\n\n";
    }
    else {
        my $tlen_cleanup = join(
            " ",
            "rm -f",
            map { shell_quote($_) } (
                $out_ref->{valid_tlen_tmp},
                $out_ref->{tlen_summary},
                $out_ref->{tlen_plot_pdf},
            ),
        );
        print $SH "$tlen_cleanup\n\n";
        write_cmd($SH, "Step 07: TLEN QC skipped",
            "echo " . shell_quote("[INFO] RUN_TLEN_QC=0; TLEN QC skipped"),
            "$dir_ref->{log}/07.extract_tlen.log");
        write_cmd($SH, "Step 08: TLEN plot skipped",
            "echo " . shell_quote("[INFO] RUN_TLEN_QC=0; TLEN plot skipped"),
            "$dir_ref->{log}/08.tlen_plot.log");
    }

    $cmd = join(
        " ",
        shell_quote($perl),
        shell_quote($script_ref->{discordant_extract}),
        "--Config", shell_quote($config),
        "--Bam",    shell_quote($out_ref->{target_bam}),
        "--Sample", shell_quote($sample),
        "--Out",    shell_quote($out_ref->{discordant_pairs}),
    );
    write_cmd($SH, "Step 09: extract_discordant_pairs", $cmd,
        "$dir_ref->{log}/09.extract_discordant_pairs.log");

    $cmd = join(
        " ",
        shell_quote($perl),
        shell_quote($script_ref->{discordant_cluster}),
        "--Config",         shell_quote($config),
        "--Input",          shell_quote($out_ref->{discordant_pairs}),
        "--AssignmentsOut", shell_quote($out_ref->{discordant_assignments}),
        "--ClustersOut",    shell_quote($out_ref->{discordant_clusters}),
    );
    write_cmd($SH, "Step 10: cluster_discordant_pairs", $cmd,
        "$dir_ref->{log}/10.cluster_discordant_pairs.log");

    $cmd = join(
        " ",
        shell_quote($perl),
        shell_quote($script_ref->{merge}),
        "--Sample",               shell_quote($sample),
        "--Depth",                shell_quote($out_ref->{cusum_events}),
        "--Split",                shell_quote($out_ref->{split_clusters}),
        "--Discordant",           shell_quote($out_ref->{discordant_clusters}),
        "--DepthCandidatesOut",   shell_quote($out_ref->{depth_candidates}),
        "--SplitCandidatesOut",   shell_quote($out_ref->{split_candidates}),
        "--DiscordantCandidatesOut", shell_quote($out_ref->{discordant_candidates}),
        "--MinOverlap",           shell_quote($evidence_overlap_fraction),
        "--Out",                  shell_quote($out_ref->{merged}),
    );
    write_cmd($SH, "Step 11: merge_evidence", $cmd,
        "$dir_ref->{log}/11.merge_evidence.log");

    $cmd = join(
        " ",
        shell_quote($perl),
        shell_quote($script_ref->{annotate}),
        "--Config",         shell_quote($config),
        "--GeneAnnotation", shell_quote($gene_txt),
        "--ExonAnnotation", shell_quote($exon_txt),
        "--Input",          shell_quote($out_ref->{merged}),
        "--Out",            shell_quote($out_ref->{annotated}),
    );
    write_cmd($SH, "Step 12: annotate_candidates", $cmd,
        "$dir_ref->{log}/12.annotate_candidates.log");

    print $SH "echo \"[INFO] Genome build: $genome_build\"\n";
    print $SH "echo \"[INFO] Finished sample: $sample \\$(date)\"\n";
    print $SH "echo \"[INFO] Original BAM: $bam\"\n";
    print $SH "echo \"[INFO] Target BAM: $out_ref->{target_bam}\"\n";
    print $SH "echo \"[INFO] Gene mean depth: $out_ref->{gene_mean_depth}\"\n";
    print $SH "echo \"[INFO] Window depth-ratio table: $out_ref->{window_depth_ratio}\"\n";
    print $SH "echo \"[INFO] Window candidates: $out_ref->{window_candidates}\"\n";
    print $SH "echo \"[INFO] CUSUM events: $out_ref->{cusum_events}\"\n";
    print $SH "echo \"[INFO] Depth-ratio plot: $out_ref->{depth_plot_pdf}\"\n";
    print $SH "echo \"[INFO] Annotated split reads: $out_ref->{split_reads}\"\n";
    print $SH "echo \"[INFO] Split clusters: $out_ref->{split_clusters}\"\n";
    print $SH "echo \"[INFO] Filtered discordant pairs: $out_ref->{discordant_pairs}\"\n";
    print $SH "echo \"[INFO] Discordant cluster assignments: $out_ref->{discordant_assignments}\"\n";
    print $SH "echo \"[INFO] Discordant clusters: $out_ref->{discordant_clusters}\"\n";
    print $SH "echo \"[INFO] Depth candidates: $out_ref->{depth_candidates}\"\n";
    print $SH "echo \"[INFO] Split candidates: $out_ref->{split_candidates}\"\n";
    print $SH "echo \"[INFO] Discordant candidates: $out_ref->{discordant_candidates}\"\n";
    print $SH "echo \"[INFO] Merged candidates: $out_ref->{merged}\"\n";
    print $SH "echo \"[INFO] Final report: $out_ref->{annotated}\"\n";

    close $SH;
}

sub write_header {
    my %args = @_;

    my $FH           = $args{fh};
    my $sample       = $args{sample};
    my $config       = $args{config};
    my $genome_build = $args{genome_build};
    my $core_gene_list = $args{core_gene_list};
    my $original_bam = $args{original_bam};
    my $target_bam   = $args{target_bam};
    my $exon_txt     = $args{exon_txt};
    my $gene_txt     = $args{gene_txt};
    my $outdir       = $args{outdir};

    print $FH <<"HEADER";
#!/usr/bin/env bash
set -euo pipefail

############################################################
# Auto-generated by HCMExonDel.pl
#
# Sample       : $sample
# Config       : $config
# Genome build : $genome_build
# Gene annot.  : $gene_txt
# Exon annot.  : $exon_txt
# Gene list    : $core_gene_list
# Original BAM : $original_bam
# Target BAM   : $target_bam
# Outdir       : $outdir
#
# Workflow:
#   Step 00. extract_target_bam
#   Step 01. gene_mean_depth_and_window_ratio
#   Step 02. cusum_depth_evidence
#   Step 03. plot_depth_ratio
#   Step 04. extract_raw_sa_reads
#   Step 05. annotate_split_read_breakpoints
#   Step 06. cluster_split_reads
#   Step 07. extract_valid_tlen
#   Step 08. plot_tlen_distribution
#   Step 09. extract_discordant_pairs
#   Step 10. cluster_discordant_pairs
#   Step 11. merge_evidence
#   Step 12. annotate_candidates
#
# Evidence merge logic:
#   split cluster is the candidate backbone;
#   depth and discordant evidence are used for validation.
############################################################

HEADER
}

sub write_block {
    my ($fh, $message) = @_;
    print $fh "echo \"[INFO] $message \$(date)\"\n";
}

sub write_cmd {
    my ($fh, $message, $cmd, $log) = @_;

    print $fh "echo \"[INFO] $message \$(date)\"\n";
    print $fh "$cmd > " . shell_quote($log) . " 2>&1\n\n";
}

# ============================================================
# Directory and output path preparation
# ============================================================

sub prepare_sample_dirs {
    my ($outdir, $sample) = @_;

    my %dir;

    $dir{sample}     = "$outdir/$sample";
    $dir{log}        = "$dir{sample}/logs";
    $dir{target_bam} = "$dir{sample}/S1_TargetBam";
    $dir{depth}      = "$dir{sample}/S2_Depth";
    $dir{split}      = "$dir{sample}/S3_SplitReads";
    $dir{discordant} = "$dir{sample}/S4_DiscordantReads";
    $dir{candidate}  = "$dir{sample}/S5_Candidates";
    $dir{report}     = "$dir{sample}/S6_Report";

    # The main script only prepares the sample directory and writes the run shell.
    # Analysis-stage directories are created by the generated shell at execution time.
    make_path($dir{sample}) unless -d $dir{sample};

    return %dir;
}

sub build_sample_outputs {
    my ($sample, $dir_ref) = @_;

    my %out = (
        target_bam             => "$dir_ref->{target_bam}/$sample.target.bam",
        target_bai             => "$dir_ref->{target_bam}/$sample.target.bam.bai",
        target_bed             => "$dir_ref->{target_bam}/$sample.target_regions.bed",
        target_region_tsv      => "$dir_ref->{target_bam}/$sample.target_regions.tsv",
        target_bam_summary     => "$dir_ref->{target_bam}/$sample.target_bam.summary.tsv",

        gene_mean_depth        => "$dir_ref->{depth}/gene_mean_depth.tsv",
        window_depth_ratio     => "$dir_ref->{depth}/window_depth_ratio.tsv",
        window_candidates      => "$dir_ref->{depth}/window_candidates.tsv",
        gene_depth_dir         => "$dir_ref->{depth}/gene_depth",
        cusum_events           => "$dir_ref->{depth}/cusum_events.tsv",
        depth_plot_pdf         => "$dir_ref->{depth}/depth_ratio.pdf",

        split_reads            => "$dir_ref->{split}/split_reads.tsv",
        split_clusters         => "$dir_ref->{split}/split_clusters.tsv",

        valid_tlen_tmp         => "$dir_ref->{discordant}/.valid_tlen.tmp.tsv",
        tlen_plot_pdf          => "$dir_ref->{discordant}/tlen_distribution.pdf",
        tlen_summary           => "$dir_ref->{discordant}/tlen_summary.tsv",
        discordant_pairs       => "$dir_ref->{discordant}/discordant_pairs.tsv",
        discordant_assignments => "$dir_ref->{discordant}/discordant_cluster_assignments.tsv",
        discordant_clusters    => "$dir_ref->{discordant}/discordant_clusters.tsv",

        depth_candidates       => "$dir_ref->{candidate}/depth_candidates.tsv",
        split_candidates       => "$dir_ref->{candidate}/split_candidates.tsv",
        discordant_candidates  => "$dir_ref->{candidate}/discordant_candidates.tsv",
        merged                 => "$dir_ref->{candidate}/merged_candidates.tsv",
        annotated              => "$dir_ref->{report}/annotated_candidates.tsv",
    );

    return %out;
}

# ============================================================
# BAM list
# ============================================================

sub read_bam_list {
    my ($file) = @_;

    my @list;
    my %seen_sample;
    my %seen_bam;

    open my $FH, "<", $file
        or die "[ERROR] Cannot open BAM list: $file\n";

    my $line_no = 0;

    while (my $line = <$FH>) {
        $line_no++;
        chomp $line;
        $line =~ s/\r$//;

        next if $line =~ /^\s*$/;
        next if $line =~ /^\s*#/;

        my @fields = split /\t/, $line, -1;

        die "[ERROR] Invalid BAM list format at line $line_no:\n"
          . "        $line\n"
          . "        Expected two TAB-delimited columns:\n"
          . "        SampleID    /absolute/path/to/sample.bam\n"
            unless @fields == 2;

        my ($sample, $bam) = @fields;

        $sample =~ s/^\s+|\s+$//g;
        $bam    =~ s/^\s+|\s+$//g;

        die "[ERROR] Empty sample name at line $line_no\n"
            unless defined $sample && $sample ne "";

        die "[ERROR] Empty BAM path at line $line_no\n"
            unless defined $bam && $bam ne "";

        die "[ERROR] Invalid sample name at line $line_no: $sample\n"
          . "        Sample name must start with a letter or number and may then contain letters, numbers, dot, underscore and hyphen.\n"
            unless $sample =~ /^[A-Za-z0-9][A-Za-z0-9_.-]*$/;

        die "[ERROR] BAM path must be absolute at line $line_no:\n"
          . "        $bam\n"
            unless $bam =~ m{^/};

        die "[ERROR] BAM file does not end with .bam at line $line_no:\n"
          . "        $bam\n"
            unless $bam =~ /\.bam$/i;

        my $abs_bam = abs_path($bam);

        die "[ERROR] BAM file not found at line $line_no:\n"
          . "        $bam\n"
            unless defined $abs_bam && -s $abs_bam;

        check_bam_index($abs_bam);

        die "[ERROR] Duplicated sample name detected at line $line_no: $sample\n"
            if exists $seen_sample{$sample};

        die "[ERROR] Duplicated BAM path detected at line $line_no:\n"
          . "        $abs_bam\n"
            if exists $seen_bam{$abs_bam};

        $seen_sample{$sample} = 1;
        $seen_bam{$abs_bam}  = 1;

        push @list, {
            sample => $sample,
            bam    => $abs_bam,
        };
    }

    close $FH;

    return @list;
}

sub check_bam_index {
    my ($bam) = @_;

    my $bai1 = "$bam.bai";
    my $bai2 = $bam;
    $bai2 =~ s/\.bam$/.bai/i;

    die "[ERROR] BAM index not found for: $bam\n"
      . "        Expected: $bai1 or $bai2\n"
        unless -s $bai1 || -s $bai2;

    return 1;
}

# ============================================================
# Config and path helpers
# ============================================================

sub read_config {
    my ($file) = @_;

    my %conf;
    my $line_no = 0;
    my %allowed = map { $_ => 1 } qw(
        PYTHON SAMTOOLS HCM_CORE_GENE_LIST TARGET_REGION_FLANK TARGET_BAM_THREADS
        MIN_MAPQ EXCLUDE_DUPLICATES WINDOW_SIZE WINDOW_STEP MIN_WINDOW_SIZE
        MIN_GENE_MEAN_DEPTH DEL_DEPTH_RATIO_CUTOFF MIN_CONSECUTIVE_DEL_WINDOWS
        CUSUM_BASELINE CUSUM_BIN_SIZE CUSUM_K CUSUM_H CUSUM_DEL_RATIO
        CUSUM_MIN_BINS CUSUM_MIN_LEN CUSUM_EDGE_RATIO CUSUM_RECOVER_RATIO
        CUSUM_RECOVER_BINS SPLIT_READ_THREADS SA_SPLIT_CLUSTER_WINDOW
        SA_SPLIT_MIN_SUPPORT_READS RUN_TLEN_QC TLEN_MAX_PAIRS TLEN_THREADS
        TLEN_RANDOM_SEED MIN_DISCORDANT_INSERT_SIZE MIN_DISCORDANT_PAIRS
        DISCORDANT_CLUSTER_DISTANCE FILTER_DELETION_ORIENTATION
        DISCORDANT_READ_THREADS EVIDENCE_OVERLAP_FRACTION MIN_EXON_OVERLAP_FRACTION
    );

    open my $FH, "<", $file
        or die "[ERROR] Cannot open config file: $file\n";

    while (my $line = <$FH>) {
        $line_no++;
        chomp $line;
        $line =~ s/\r$//;

        next if $line =~ /^\s*$/;
        next if $line =~ /^\s*#/;

        # Remove inline comments: KEY=value # comment
        $line =~ s/\s+#.*$//;

        die "[ERROR] Invalid config line $file:$line_no: $line\n"
            unless $line =~ /^\s*([A-Z][A-Z0-9_]*)\s*=\s*(.*?)\s*$/;

        my ($key, $val) = ($1, $2);
        $val =~ s/^\s+|\s+$//g;
        $val =~ s/^['"]//;
        $val =~ s/['"]$//;

        die "[ERROR] Empty value for config parameter $key in $file\n"
            if $val eq "";

        die "[ERROR] Duplicate config parameter '$key' in $file\n"
            if exists $conf{$key};

        die "[ERROR] Unknown config parameter '$key' in $file\n"
            unless $allowed{$key};

        $conf{$key} = $val;
    }

    close $FH;

    return %conf;
}

sub get_conf_required {
    my ($conf_ref, $key) = @_;

    die "[ERROR] Required config parameter missing: $key\n"
        unless exists $conf_ref->{$key}
            && defined $conf_ref->{$key}
            && $conf_ref->{$key} ne "";

    return $conf_ref->{$key};
}

sub normalize_bool {
    my ($v, $name) = @_;

    die "[ERROR] $name is undefined\n"
        unless defined $v;

    $v =~ s/^\s+|\s+$//g;

    return 1 if $v =~ /^(1|yes|true|on)$/i;
    return 0 if $v =~ /^(0|no|false|off)$/i;

    die "[ERROR] $name must be a boolean value: 1/0, yes/no, true/false, on/off. Observed: $v\n";
}

sub resolve_config_path {
    my ($path) = @_;

    die "[ERROR] Empty path provided\n"
        unless defined $path && $path ne "";

    die "[ERROR] Configured paths must be absolute: $path\n"
        unless $path =~ m{^/};

    my $abs = abs_path($path);

    die "[ERROR] Path not found: $path\n"
        unless defined $abs;

    return $abs;
}

# ============================================================
# Validation helpers
# ============================================================

sub validate_cusum_params {
    my ($cusum_ref) = @_;

    validate_baseline("CUSUM_BASELINE",              $cusum_ref->{baseline});
    validate_positive_integer("CUSUM_BIN_SIZE",       $cusum_ref->{bin_size});
    validate_nonnegative_number("CUSUM_K",            $cusum_ref->{k});
    validate_positive_number("CUSUM_H",               $cusum_ref->{h});
    validate_fraction("CUSUM_DEL_RATIO",              $cusum_ref->{del_ratio});
    validate_positive_integer("CUSUM_MIN_BINS",       $cusum_ref->{min_bins});
    validate_positive_integer("CUSUM_MIN_LEN",        $cusum_ref->{min_len});
    validate_fraction("CUSUM_EDGE_RATIO",             $cusum_ref->{edge_ratio});
    validate_fraction("CUSUM_RECOVER_RATIO",          $cusum_ref->{recover_ratio});
    validate_positive_integer("CUSUM_RECOVER_BINS",   $cusum_ref->{recover_bins});

    return 1;
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
        unless defined $v && $v =~ /^(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;

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

sub validate_nonnegative_integer {
    my ($name, $v) = @_;

    die "[ERROR] $name must be a non-negative integer: $v\n"
        unless defined $v && $v =~ /^\d+$/ && $v >= 0;

    return 1;
}

sub validate_integer {
    my ($name, $v) = @_;

    die "[ERROR] $name must be an integer: $v\n"
        unless defined $v && $v =~ /^-?\d+$/;

    return 1;
}

sub validate_positive_number {
    my ($name, $v) = @_;

    die "[ERROR] $name must be numeric: $v\n"
        unless defined $v && $v =~ /^(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;

    die "[ERROR] $name must be > 0: $v\n"
        unless $v > 0;

    return 1;
}

sub validate_nonnegative_number {
    my ($name, $v) = @_;

    die "[ERROR] $name must be numeric: $v\n"
        unless defined $v && $v =~ /^(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;

    die "[ERROR] $name must be >= 0: $v\n"
        unless $v >= 0;

    return 1;
}

sub check_executable {
    my ($cmd, $name) = @_;

    die "[ERROR] Empty executable for $name\n"
        unless defined $cmd && $cmd ne "";

    if ($cmd =~ m{/}) {
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

sub check_script {
    my ($path, $name) = @_;

    die "[ERROR] Pipeline script not found for $name: $path\n"
        unless defined $path && -s $path;

    return 1;
}

sub check_exon_txt_format {
    my ($file) = @_;

    my @required = qw/Gene Transcript Exon Chrom Start End Strand/;

    check_table_header(
        file     => $file,
        desc     => "MANERefSeqExons.tsv",
        required => \@required,
    );
}

sub check_gene_txt_format {
    my ($file) = @_;

    my @required = qw/Gene Transcript Chrom Start End Strand/;

    check_table_header(
        file     => $file,
        desc     => "MANERefSeqGenes.tsv",
        required => \@required,
    );
}

sub check_table_header {
    my %args = @_;

    my $file     = $args{file};
    my $desc     = $args{desc};
    my $required = $args{required};

    open my $FH, "<", $file
        or die "[ERROR] Cannot open $desc file: $file\n";

    my $header = <$FH>;

    close $FH;

    die "[ERROR] Empty $desc file: $file\n"
        unless defined $header;

    chomp $header;
    $header =~ s/\r$//;

    my @cols = split /\t/, $header, -1;
    my %seen = map { $_ => 1 } @cols;

    for my $col (@$required) {
        die "[ERROR] Required column '$col' not found in $desc: $file\n"
            unless exists $seen{$col};
    }

    return 1;
}

sub check_core_gene_list_format {
    my ($file) = @_;

    open my $FH, "<", $file
        or die "[ERROR] Cannot open HCM_CORE_GENE_LIST: $file\n";

    my $line_no      = 0;
    my $data_line_no = 0;
    my $gene_count   = 0;

    while (my $line = <$FH>) {
        chomp $line;
        $line =~ s/\r$//;
        $line_no++;

        next if $line =~ /^\s*$/;
        next if $line =~ /^\s*#/;

        $line =~ s/\s+#.*$//;
        $line =~ s/^\s+|\s+$//g;

        next if $line eq "";

        $data_line_no++;

        my @f = split /\s+/, $line;
        my $gene = $f[0];
        $gene = "" unless defined $gene;

        if ($data_line_no == 1 && $gene =~ /^Gene$/i) {
            next;
        }

        die "[ERROR] Empty gene symbol in HCM_CORE_GENE_LIST at line $line_no\n"
            if $gene eq "";

        die "[ERROR] Invalid gene symbol in HCM_CORE_GENE_LIST at line $line_no: $gene\n"
            unless $gene =~ /^[A-Za-z0-9_.-]+$/;

        $gene_count++;
    }

    close $FH;

    die "[ERROR] No core genes found in HCM_CORE_GENE_LIST: $file\n"
        unless $gene_count > 0;

    return 1;
}

# ============================================================
# Shell and output helpers
# ============================================================

sub shell_quote {
    my ($str) = @_;

    die "[ERROR] Undefined shell argument\n"
        unless defined $str;

    return "''" if $str eq "";

    if ($str =~ /^[A-Za-z0-9_\.\-\/:=,+]+$/) {
        return $str;
    }

    $str =~ s/'/'\\''/g;

    return "'$str'";
}

sub print_summary {
    my %args = @_;

    print "\n";
    print "============================================================\n";
    print "HCMExonDel pipeline setup\n";
    print "============================================================\n";
    print "Config             : $args{config}\n";
    print "Project root       : $args{project_root}\n";
    print "Genome build       : $args{genome_build}\n";
    print "BAM list           : $args{bam_list}\n";
    print "Exon TXT           : $args{exon_txt}\n";
    print "Gene TXT           : $args{gene_txt}\n";
    print "HCM core gene list : $args{core_gene_list}\n";
    print "Output directory   : $args{outdir}\n";
    print "PERL               : $args{perl} ($^V)\n";
    print "PYTHON             : $args{python}\n";
    print "SAMTOOLS           : $args{samtools}\n";
    print "Sample number      : $args{sample_num}\n";

    print "\n";
    print "Workflow switches\n";
    print "------------------------------------------------------------\n";
    print "RUN_TLEN_QC                 : $args{run_tlen_qc}\n";

    print "\n";
    print "Target BAM and read extraction parameters\n";
    print "------------------------------------------------------------\n";
    print "TARGET_REGION_FLANK     : $args{target_region_flank}\n";
    print "TARGET_BAM_THREADS      : $args{target_bam_threads}\n";
    print "SPLIT_READ_THREADS      : $args{split_read_threads}\n";
    print "DISCORDANT_READ_THREADS : $args{discordant_read_threads}\n";

    if ($args{run_tlen_qc}) {
        my $t = $args{tlen_ref};

        print "\n";
        print "TLEN QC parameters\n";
        print "------------------------------------------------------------\n";
        print "TLEN_MAX_PAIRS   : $t->{max_pairs}\n";
        print "TLEN_THREADS     : $t->{threads}\n";
        print "TLEN_RANDOM_SEED : $t->{random_seed}\n";
    }

    print "\n";
    print "Evidence merging\n";
    print "------------------------------------------------------------\n";
    print "Merge backbone            : split-read cluster\n";
    print "Supporting evidence       : CUSUM depth + discordant read pairs\n";
    print "EVIDENCE_OVERLAP_FRACTION : $args{evidence_overlap_fraction}\n";

    print "\n";
    print "CUSUM parameters\n";
    print "------------------------------------------------------------\n";

    my $c = $args{cusum_ref};

    print "CUSUM_BASELINE      : $c->{baseline}\n";
    print "CUSUM_BIN_SIZE      : $c->{bin_size}\n";
    print "CUSUM_K             : $c->{k}\n";
    print "CUSUM_H             : $c->{h}\n";
    print "CUSUM_DEL_RATIO     : $c->{del_ratio}\n";
    print "CUSUM_MIN_BINS      : $c->{min_bins}\n";
    print "CUSUM_MIN_LEN       : $c->{min_len}\n";
    print "CUSUM_EDGE_RATIO    : $c->{edge_ratio}\n";
    print "CUSUM_RECOVER_RATIO : $c->{recover_ratio}\n";
    print "CUSUM_RECOVER_BINS  : $c->{recover_bins}\n";
    print "============================================================\n\n";
}

sub write_run_info_shell {
    my %args = @_;
    my $fh = $args{fh};

    my @metadata = (
        [SampleID        => $args{sample}],
        [OriginalBAM     => $args{bam}],
        [ConfigFile      => $args{config}],
        [GenomeBuild     => $args{genome_build}],
        [GeneAnnotation  => $args{gene_annotation}],
        [ExonAnnotation  => $args{exon_annotation}],
        [HCMCoreGeneList => $args{core_gene_list}],
        [PerlPath        => $args{perl}],
        [PythonPath      => $args{python}],
        [SamtoolsPath    => $args{samtools}],
    );

    print $fh "{\n";
    print $fh "  printf 'Key\\tValue\\n'\n";

    for my $item (@metadata) {
        my ($key, $value) = @$item;
        $value = '' unless defined $value;
        $value =~ s/[\t\r\n]+/ /g;
        print $fh "  printf '%s\\t%s\\n' "
            . shell_quote($key) . " " . shell_quote($value) . "\n";
    }

    for my $key (sort keys %{ $args{conf_ref} }) {
        my $value = $args{conf_ref}{$key};
        $value = '' unless defined $value;
        $value =~ s/[\t\r\n]+/ /g;
        print $fh "  printf '%s\\t%s\\n' "
            . shell_quote("Config_$key") . " " . shell_quote($value) . "\n";
    }

    my $perl_version_cmd = shell_quote($args{perl}) . " -e "
        . shell_quote('printf "%vd\\n", $^V');
    my $python_version_cmd = shell_quote($args{python}) . " -c "
        . shell_quote('import sys; print(sys.version.split()[0])');
    my $samtools_version_cmd = shell_quote($args{samtools})
        . " --version | sed -n '1p'";

    print $fh q{  printf '%s\t%s\n' 'PerlVersion' "} . '$(' . $perl_version_cmd . ')"' . "\n";
    print $fh q{  printf '%s\t%s\n' 'PythonVersion' "} . '$(' . $python_version_cmd . ')"' . "\n";
    print $fh q{  printf '%s\t%s\n' 'SamtoolsVersion' "} . '$(' . $samtools_version_cmd . ')"' . "\n";
    print $fh q{  printf '%s\t%s\n' 'GeneratedAt' "$(date '+%Y-%m-%d %H:%M:%S %z')"} . "\n";
    print $fh "} > " . shell_quote($args{file}) . "\n\n";
}

sub usage {
    return <<'USAGE';
Usage:
  perl HCMExonDel.pl \
    --Config conf/hcm_exondel.conf \
    --BamList bam.list \
    --OutDir results \
    --GenomeBuild hg38

Required:
  --Config      HCMExonDel configuration file.
  --BamList     TAB-delimited BAM list with two columns:
                  SampleID    /absolute/path/to/sample.bam

Optional:
  --OutDir      Output directory. Default: results
  --GenomeBuild Reference genome build: hg19 or hg38. Default: hg19.
  --Help        Show this help message.

Generated output:
  OUTDIR/SAMPLE/run_SAMPLE.sh

Workflow:
  Step 00. Extract target BAM
  Step 01. Calculate gene and window depth
  Step 02. Run CUSUM depth analysis
  Step 03. Generate depth-ratio plot
  Step 04. Extract raw SA split reads
  Step 05. Annotate split-read breakpoints in split_reads.tsv
  Step 06. Cluster valid split-read deletion records
  Step 07. Extract TLEN values (optional)
  Step 08. Plot TLEN distribution (optional)
  Step 09. Extract and filter Discordant read pairs directly from BAM
  Step 10. Cluster Discordant read pairs
  Step 11. Integrate Depth, Split, and Discordant evidence
  Step 12. Annotate candidates

Notes:
  1. HCMExonDel analyzes the genes listed in HCM_CORE_GENE_LIST.
  2. The original WGS BAM is used only to generate SAMPLE.target.bam.
  3. All downstream analyses use SAMPLE.target.bam.
  4. Split-read clusters define candidate events; depth and discordant
     evidence are used for supporting validation.
  5. Genome-specific MANE annotations are selected automatically from db/Annotation/<GenomeBuild>/.
  6. Existing run_SAMPLE.sh files are overwritten when regenerated.

Example:
  perl HCMExonDel.pl \
    --Config conf/hcm_exondel.conf \
    --BamList example/input_bam.list \
    --OutDir test/test_results

Then run:
  bash test/test_results/SAMPLE/run_SAMPLE.sh
USAGE
}

