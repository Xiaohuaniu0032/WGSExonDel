# HCMExonDel

HCMExonDel is a lightweight, configurable pipeline for detecting candidate exon-level deletions in hypertrophic cardiomyopathy (HCM)-related core genes from whole-genome sequencing (WGS) BAM files.

The current workflow is **split-read-centered**: split-read clusters define candidate deletion intervals, while CUSUM depth evidence and discordant read-pair evidence are used to validate and prioritize those intervals. Final candidates are annotated against MANE RefSeq exon coordinates for manual review and experimental validation.

> **Research use only.** HCMExonDel does not perform clinical variant classification, and all reported candidates require independent review and validation.

---

## Table of Contents

- [Overview](#overview)
- [Current Design](#current-design)
- [Workflow](#workflow)
- [Evidence Integration](#evidence-integration)
- [Requirements](#requirements)
- [Repository Structure](#repository-structure)
- [Input Files](#input-files)
- [Configuration](#configuration)
- [Quick Start](#quick-start)
- [Output Structure](#output-structure)
- [Key Output Files](#key-output-files)
- [Candidate Interpretation](#candidate-interpretation)
- [Test Example: 25B09089386](#test-example-25b09089386)
- [Manual Review and Validation](#manual-review-and-validation)
- [Running Individual Modules](#running-individual-modules)
- [Troubleshooting](#troubleshooting)
- [Limitations](#limitations)
- [Citation](#citation)
- [License](#license)
- [Contact](#contact)

---

## Overview

Exon-level deletions and other structural variants may be missed by analysis pipelines focused mainly on SNVs and small indels. HCMExonDel was developed for WGS-based reanalysis of genetically unresolved HCM cases and families, with emphasis on deletions affecting established HCM-related genes.

HCMExonDel integrates three complementary signals:

1. **Read-depth reduction** detected by base-level CUSUM analysis.
2. **Split-read evidence** derived from supplementary alignments and breakpoint-supporting reads.
3. **Discordant read-pair evidence** derived from abnormally large, deletion-like paired-end fragments.

The pipeline does not call deletions across the whole genome. It always analyzes genes listed in `HCM_CORE_GENE_LIST`.

---

## Current Design

The current implementation differs from earlier versions of the pipeline in several important ways.

### Core-gene-only analysis

HCMExonDel always analyzes the HCM core genes listed in:

```text
HCM_CORE_GENE_LIST
```

Core-gene-only analysis is now controlled exclusively by `HCM_CORE_GENE_LIST`. `ANALYZE_CORE_GENES_ONLY` is not a supported configuration parameter.

### Target BAM first

For each sample, the original WGS BAM is used only in Step 00 to generate:

```text
SAMPLE.target.bam
```

The target BAM contains the core-gene regions plus `TARGET_REGION_FLANK`. All downstream depth, split-read, TLEN, and discordant-read analyses use this smaller target BAM.

### Split-read-centered candidate generation

A passed split-read cluster is the candidate event backbone. CUSUM depth intervals and discordant-read clusters are evaluated as supporting evidence for that split interval.

Consequently, a depth-only or discordant-only signal is not included in the final merged candidate file when no qualifying split-read cluster is present.

### Candidate gene BAM extraction removed

Candidate gene BAM extraction is no longer part of the standard workflow. The former parameter `CANDIDATE_BAM_FLANK` and the former `05.gene_bam/` output directory have been removed.

For IGV review, use the original WGS BAM and/or the generated `SAMPLE.target.bam`.

---

## Workflow

The main script, `HCMExonDel.pl`, validates the configuration and input BAM list, then generates one executable shell script for each sample. It does not directly execute the complete analysis.

```text
Original coordinate-sorted WGS BAM
                |
                v
Step 00  extract_target_bam
         Core genes + flanking regions -> SAMPLE.target.bam
                |
                v
Step 01  gene_mean_depth_and_window_ratio
         Gene mean depth, all window depth ratios, window candidates
                |
                v
Step 02  cusum_depth_evidence
         Base-level CUSUM deletion-supporting intervals
                |
                v
Step 03  plot_depth_ratio
         Per-gene depth-ratio PDF
                |
                v
Step 04  extract_raw_sa_reads
         Raw Primary/SA alignment pairs; multi-SA reads are excluded
                |
                v
Step 05  annotate_split_read_breakpoints
         Add deletion coordinates/status to split_reads.tsv in place
                |
                v
Step 06  cluster_valid_split_read_deletions
         PASS deletion records -> split_clusters.tsv only
                |
                +---------------------------------------+
                |                                       |
                v                                       v
Step 07  extract_tlen                         Step 09  extract_discordant_pairs
Step 08  plot_tlen_distribution                        Direct BAM parsing and filtering
         Optional library QC                            |
                                                       v
                                             Step 10  cluster_discordant_pairs
                                                       Cluster summaries only
                \                                       /
                 \                                     /
                  v                                   v
Step 11             merge_evidence
        Split-centered candidates with unchanged evidence-level rules
                                |
                                v
Step 12             annotate_candidates
        MANE gene/transcript/exon annotation and final report
```

For each sample, the main script generates:

```text
OUTDIR/SAMPLE/run_SAMPLE.sh
```

Run this shell script manually or submit it to an HPC scheduler.

---

## Evidence Integration

### 1. Window-level depth analysis

`RunGeneMeanDepth.pl` calculates gene-level and window-level depth statistics using the target BAM. A window is considered deletion-supporting when its depth ratio is below `DEL_DEPTH_RATIO_CUTOFF`.

The window-level results are useful for:

- detecting consecutive low-depth windows;
- visualizing depth changes across each core gene;
- checking whether a candidate interval shows an expected heterozygous depth reduction.

The window candidate file is **not** the depth input used by the current evidence-merging step.

### 2. CUSUM depth evidence

`RunSampleCUSUMDepth.pl` applies `CUSUMDepthDel.pl` to retained per-gene base-depth files. The combined output is:

```text
cusum_events.tsv
```

Per-gene CUSUM intermediate files are written to a temporary directory and
removed after the sample-level event table is generated.

This file is the depth-evidence input passed to `MergeEvidence.pl`.

For each split interval, overlapping CUSUM depth intervals are combined. Depth support is assigned when their covered portion reaches `EVIDENCE_OVERLAP_FRACTION` of the split interval.

### 3. Split-read evidence

Split-read analysis is separated into three steps:

1. `ExtractSASplitReads.pl` writes raw Primary/SA alignment information to `split_reads.tsv`. It keeps only reads with exactly one valid SA entry.
2. `IdentifyIntragenicDeletions.pl` annotates each row in `split_reads.tsv` in place. It adds the inferred deletion interval, junction type, breakpoint status, and breakpoint reason while retaining all raw read records. The input and output paths are intentionally the same, and the script uses a temporary file before replacing the original.
3. `ClusterSASplitReads.pl` reads the annotated `split_reads.tsv` and clusters only rows with `Breakpoint_Status=PASS`, `Junction_Type=Deletion_Gap`, and `Breakpoint_Reason=Valid_Deletion_Gap`. It writes only the cluster summary `split_clusters.tsv`; no event or assignment table is generated.

A cluster passes when it contains at least `SA_SPLIT_MIN_SUPPORT_READS` unique supporting read names. Only `Pass` clusters are used as candidate backbones in S5.

### 4. Discordant read-pair evidence

Discordant analysis is separated into two steps while preserving the existing filtering and clustering rules:

1. `ExtractDiscordantPairs.pl` parses the target BAM directly, applies the existing basic SAM flag, MAPQ, optional duplicate, same-chromosome, positive-TLEN, minimum insert-size, optional FR-orientation, and valid-interval filters, then writes one retained record per read pair to `discordant_pairs.tsv`. It also calculates the existing inner-gap interval. No intermediate raw paired-alignment TSV is written.
2. `ClusterDiscordantReads.pl` applies the existing coordinate-distance clustering and minimum unique-pair support rule. It writes only `discordant_clusters.tsv`; no per-pair assignment table is generated.

This refactor combines BAM extraction and Discordant filtering into one pass. It does not change the existing filtering thresholds, interval calculation, or clustering criteria. S5 continues to compare the resulting PASS discordant-cluster intervals with split-read candidate intervals using `EVIDENCE_OVERLAP_FRACTION`.

### 5. Evidence levels

The current confidence model is:

| Evidence level | Required evidence | Interpretation |
|---|---|---|
| `High` | Split + Depth + Discordant | Breakpoint, depth, and paired-end evidence are all present. |
| `Moderate` | Split + Depth, without qualifying Discordant support | Breakpoint and depth evidence are present. |
| `Low` | Split without qualifying Depth support | A split cluster is present, but depth support does not pass the overlap threshold. Discordant evidence alone does not raise the candidate to `Moderate`. |

The final merged file contains all passed split-centered candidates and records their evidence level. The current workflow does not generate a separate `.all.tsv` candidate file.

---

## Requirements

### Operating system

HCMExonDel is designed for Linux environments, including Linux workstations, servers, and HPC clusters.

### Required software

- Perl 5
- Python 3
- samtools
- Bash

### Perl modules

The pipeline mainly uses standard Perl modules, including:

```text
Getopt::Long
File::Basename
File::Path
Cwd
FindBin
```

### Python packages

The plotting scripts require:

```text
pandas
numpy
matplotlib
```

Install them when necessary:

```bash
python3 -m pip install pandas numpy matplotlib
```

### Recommended review software

- IGV

---

## Repository Structure

```text
HCMExonDel/
├── HCMExonDel.pl
├── bin/                         # Core pipeline scripts
├── conf/                        # Example configuration
├── db/
│   ├── Genes/                   # HCM and phenotype-associated gene lists
│   ├── Annotation/
│   │   ├── hg19/                 # MANE RefSeq gene/exon coordinates
│   │   └── hg38/                 # MANE RefSeq gene/exon coordinates
│   └── MANE/                    # MANE source and transcript tables
├── scripts/
│   ├── BuildReference/
│   ├── Simulation/
│   ├── Benchmark/SVCaller/
│   ├── Plotting/
│   ├── Annotation/
│   ├── HPO/
│   └── Utilities/
├── example/
├── test/
├── docs/
│   ├── Figures/
│   └── refactoring_manifest.tsv
├── README.md
└── .gitignore
```

Perl and Python script filenames use UpperCamelCase. Database directories and files use UpperCamelCase, while result TSV/PDF filenames use lowercase snake_case. Sample result directories are named `S1_TargetBam` through `S6_Report`.

## Input Files

### 1. Coordinate-sorted WGS BAM

Each sample must have a coordinate-sorted BAM file and a BAM index.

Accepted index names are:

```text
sample.bam.bai
sample.bai
```

Additional requirements:

- the BAM path in the BAM list must be absolute;
- the filename must end with `.bam`;
- the BAM reference build and chromosome naming must match the annotation files;
- sample names may contain letters, numbers, dots, underscores, and hyphens.

### 2. BAM list

The BAM list must contain exactly two TAB-delimited columns:

```text
SampleID<TAB>/absolute/path/to/sample.bam
```

Example:

```text
Sample_01<TAB>/absolute/path/to/Sample_01.bam
Sample_02<TAB>/absolute/path/to/Sample_02.bam
```

Rules:

- no header is required;
- blank lines are ignored;
- lines beginning with `#` are ignored;
- duplicate sample IDs are not allowed;
- duplicate BAM paths are not allowed.

### 3. MANE RefSeq exon annotation

Selected automatically according to `--GenomeBuild` from:

```text
db/Annotation/<GenomeBuild>/MANERefSeqExons.tsv
```

Required columns:

```text
Gene  Transcript  Exon  Chrom  Start  End  Strand
```

Coordinates must be **1-based closed**.

### 4. MANE RefSeq gene annotation

Selected automatically according to `--GenomeBuild` from:

```text
db/Annotation/<GenomeBuild>/MANERefSeqGenes.tsv
```

Required columns:

```text
Gene  Transcript  Chrom  Start  End  Strand
```

Coordinates must be **1-based closed**.

The file is used to build target regions and perform depth analysis.

### 5. HCM core-gene list

Configured by:

```text
HCM_CORE_GENE_LIST
```

Example:

```text
Gene  Classification
MYBPC3  Definitive
MYH7  Definitive
FHOD3 Moderate
```

Only the first column is used as the gene symbol. A header named `Gene` is allowed. The classification column is optional.

A core gene absent from the selected `MANERefSeqGenes.tsv` is reported and skipped during target-region generation; it does not terminate the entire sample analysis.

---

## Configuration

Copy and edit the example configuration:

```bash
cp conf/hcm_exondel.example.conf conf/hcm_exondel.conf
```

The example configuration contains placeholder absolute paths. Copy it to `conf/hcm_exondel.conf`, then replace every configured path with an existing absolute path on the target system before running the pipeline. The reference genome build is selected on the command line; MANE annotation paths are not configured in the conf file.

### Runtime and gene-list paths

| Parameter | Example default | Description |
|---|---:|---|
| `PYTHON` | `/absolute/path/to/python` | Python executable used for plotting. |
| `SAMTOOLS` | `/absolute/path/to/samtools` | samtools executable. |
| `HCM_CORE_GENE_LIST` | `/absolute/path/to/HCMExonDel/db/Genes/HCMClassic8Genes.tsv` | HCM core-gene list. |

All path-valued configuration entries must be absolute. `HCM_CORE_GENE_LIST` defaults to `HCMClassic8Genes.tsv`. `PYTHON` and `SAMTOOLS` must also be absolute executable paths. The Perl interpreter is taken from the Perl runtime executing `HCMExonDel.pl` and is not configured separately. The current `HCMClassic8Genes.tsv` contains eight classical sarcomeric HCM genes plus FHOD3 (nine gene records); confirm this is the intended default gene scope before publication.

Select the reference build with `--GenomeBuild hg19` or `--GenomeBuild hg38`. The default is `hg19`. The main script automatically selects `db/Annotation/<GenomeBuild>/MANERefSeqGenes.tsv` and `MANERefSeqExons.tsv`, validates both files, and records their paths in the generated shell header. The sample shell writes these paths and the runtime software versions to `run_info.tsv` when execution begins.

### Target BAM and common read filters

| Parameter | Default | Description |
|---|---:|---|
| `TARGET_REGION_FLANK` | `5000` | Bases added on both sides of each core-gene interval when generating the target BAM. |
| `TARGET_BAM_THREADS` | `4` | Threads used during target BAM extraction/sorting. |
| `MIN_MAPQ` | `20` | Minimum mapping quality. |
| `EXCLUDE_DUPLICATES` | `1` | Exclude duplicate-marked reads when supported by the corresponding module. |

### Window-depth analysis

| Parameter | Default | Description |
|---|---:|---|
| `WINDOW_SIZE` | `500` | Window size in bp. |
| `WINDOW_STEP` | `200` | Step between adjacent windows in bp. |
| `MIN_WINDOW_SIZE` | `200` | Minimum retained window length in bp. |
| `MIN_GENE_MEAN_DEPTH` | `20` | Minimum gene mean depth for reliable ratio calculation. |
| `DEL_DEPTH_RATIO_CUTOFF` | `0.65` | Window depth-ratio threshold for deletion-supporting windows. |
| `MIN_CONSECUTIVE_DEL_WINDOWS` | `3` | Minimum number of consecutive deletion-supporting windows. |

### CUSUM depth evidence

| Parameter | Default | Description |
|---|---:|---|
| `CUSUM_BASELINE` | `auto` | Baseline mode for the CUSUM calculation. |
| `CUSUM_BIN_SIZE` | `1` | CUSUM bin size in bp. |
| `CUSUM_K` | `0.1` | CUSUM reference/drift parameter. |
| `CUSUM_H` | `5` | CUSUM decision threshold. |
| `CUSUM_DEL_RATIO` | `0.65` | Ratio used to define deletion-like depth. |
| `CUSUM_MIN_BINS` | `3` | Minimum number of bins in a CUSUM interval. |
| `CUSUM_MIN_LEN` | `1` | Minimum interval length in bp. |
| `CUSUM_EDGE_RATIO` | `0.80` | Edge-extension ratio threshold. |
| `CUSUM_RECOVER_RATIO` | `0.80` | Recovery ratio threshold. |
| `CUSUM_RECOVER_BINS` | `3` | Consecutive recovery bins required to close an event. |

### Split-read clustering

| Parameter | Default | Description |
|---|---:|---|
| `SPLIT_READ_THREADS` | `4` | Threads used during split-read extraction. |
| `SA_SPLIT_CLUSTER_WINDOW` | `20` | Maximum breakpoint-coordinate distance for clustering. |
| `SA_SPLIT_MIN_SUPPORT_READS` | `5` | Minimum number of unique read names required for a passed cluster. |

### TLEN quality control

| Parameter | Default | Description |
|---|---:|---|
| `RUN_TLEN_QC` | `1` | Enable insert-size/TLEN QC. |
| `TLEN_MAX_PAIRS` | `1000000` | Maximum number of valid read pairs sampled. |
| `TLEN_THREADS` | `4` | Threads used for TLEN extraction. |
| `TLEN_RANDOM_SEED` | `20260701` | Random seed for reproducible sampling. |

The TLEN plot reports the median, MAD, robust standard deviation, and a diagnostic threshold of:

```text
median + 3 × 1.4826 × MAD
```

This diagnostic value does not automatically replace `MIN_DISCORDANT_INSERT_SIZE`. Tune the discordant threshold explicitly when library insert-size characteristics differ from the default setting.

### Discordant read-pair analysis

| Parameter | Default | Description |
|---|---:|---|
| `MIN_DISCORDANT_INSERT_SIZE` | `500` | Minimum absolute insert size for a discordant pair. |
| `MIN_DISCORDANT_PAIRS` | `3` | Minimum supporting read pairs in a passed cluster. |
| `DISCORDANT_CLUSTER_DISTANCE` | `500` | Maximum distance used to group nearby discordant pairs. |
| `FILTER_DELETION_ORIENTATION` | `1` | Require deletion-like paired-end orientation. |
| `DISCORDANT_READ_THREADS` | `4` | Threads used during discordant-read extraction. |

### Evidence integration and annotation

| Parameter | Default | Description |
|---|---:|---|
| `EVIDENCE_OVERLAP_FRACTION` | `0.90` | Minimum fraction of the split interval that must be supported by qualifying depth or discordant evidence. |
| `MIN_EXON_OVERLAP_FRACTION` | `0.20` | Minimum fraction of an exon overlapped by a candidate for that exon to be reported as affected. |

### Removed parameters

Do not include the following parameters in a current configuration file:

```text
ANALYZE_CORE_GENES_ONLY
CANDIDATE_BAM_FLANK
```

The strict configuration parser rejects unknown or duplicate keys. These historical parameters are not supported.

---

## Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/Xiaohuaniu0032/HCMExonDel.git
cd HCMExonDel
```

### 2. Prepare the configuration

```bash
cp conf/hcm_exondel.example.conf conf/hcm_exondel.conf
vim conf/hcm_exondel.conf
```

The example configuration uses placeholder absolute paths. Replace them with valid paths on your system before running the pipeline. Set every path-valued configuration entry to an existing absolute path. At minimum, verify:

```text
PYTHON
SAMTOOLS
HCM_CORE_GENE_LIST
```

The MANE annotation files are selected automatically from `db/Annotation/hg19/` or `db/Annotation/hg38/`.

### 3. Prepare the BAM list

```bash
printf 'Sample_01\t/absolute/path/to/Sample_01.bam\n' > input_bam.list
```

### 4. Generate sample-specific shell scripts

```bash
perl HCMExonDel.pl \
  --Config conf/hcm_exondel.conf \
  --BamList input_bam.list \
  --OutDir results \
  --GenomeBuild hg38
```

Expected output:

```text
results/Sample_01/run_Sample_01.sh
```

### 5. Run the analysis

```bash
bash results/Sample_01/run_Sample_01.sh
```

For multiple samples, execute or submit each generated `run_SAMPLE.sh` independently.

### 6. Regenerating sample shell scripts

The main script only generates sample-level shell scripts. `--GenomeBuild` accepts `hg19` or `hg38` and defaults to `hg19`. If `run_SAMPLE.sh` already exists, it is overwritten. Existing analysis result directories are not deleted or modified during shell generation.

Before executing a regenerated shell, review the existing sample results and remove stale outputs when a clean rerun is required.


---

## Output Structure

```text
OUTDIR/
└── SAMPLE/
    ├── run_SAMPLE.sh
    ├── run_info.tsv
    ├── logs/
    ├── S1_TargetBam/
    ├── S2_Depth/
    ├── S3_SplitReads/
    ├── S4_DiscordantReads/
    ├── S5_Candidates/
    └── S6_Report/
```

### Log files

`run_info.tsv` is created when the generated sample shell starts running. It records the sample, input BAM, selected genome build, MANE annotation paths, HCM gene list, all configuration values, executable paths and runtime software versions. Generating or regenerating the shell does not create or overwrite `run_info.tsv`.

The generated workflow writes one log per step under `logs/`:

```text
00.extract_target_bam.log
01.depth_analysis.log
02.cusum_depth.log
03.depth_ratio_plot.log
04.extract_split_reads.log
05.identify_split_deletions.log
06.cluster_split_reads.log
07.extract_tlen.log
08.tlen_plot.log
09.extract_discordant_pairs.log
10.cluster_discordant_pairs.log
11.merge_evidence.log
12.annotate_candidates.log
```

Skipped optional steps still receive a log containing the skip message.

---

## Key Output Files

### `S1_TargetBam/`

| File | Description |
|---|---|
| `SAMPLE.target.bam` | Sorted target BAM containing core genes plus configured flanks. |
| `SAMPLE.target.bam.bai` | BAM index. |
| `SAMPLE.target_regions.tsv` | Unmerged target-region records generated from the core-gene list and gene annotation. |
| `SAMPLE.target_regions.bed` | Merged BED intervals used for BAM extraction. |
| `SAMPLE.target_bam.summary.tsv` | Counts of requested, found, and skipped genes plus output paths and extraction settings. |

The BED file is generated internally; users do not need to provide a BED file.

The current S2 depth output files and their column definitions are documented in [`docs/output_format.md`](docs/output_format.md).

### `S2_Depth/`

| File or directory | Description | Used in final merge? |
|---|---|---:|
| `gene_mean_depth.tsv` | Gene-level mean-depth summary, including the per-gene depth filename used by CUSUM. | Indirectly |
| `window_depth_ratio.tsv` | One row per window, combining window depth, depth ratio and window status. | No |
| `window_candidates.tsv` | Consecutive deletion-supporting window candidates. | No |
| `gene_depth/` | Per-gene base-depth files used as CUSUM input. | **Yes, indirectly** |
| `cusum_events.tsv` | Combined sample-level CUSUM deletion intervals with gene/transcript metadata. | **Yes** |
| `depth_ratio.pdf` | Multi-page depth-ratio plot, one page per core gene. | No |

The former separate window-depth, all-window-ratio and deletion-window tables have
been consolidated into `window_depth_ratio.tsv`. Per-gene CUSUM result files are
temporary and are not retained. Base-level depth-ratio diagnostics are not part
of the standard workflow; optional utility scripts are located in
`scripts/DepthQC/`.
### `S3_SplitReads/`

| File | Description |
|---|---|
| `split_reads.tsv` | Raw Primary/SA alignment fields plus read-level breakpoint/deletion annotations added in place. All extracted rows are retained. No `Cluster_ID` is written. |
| `split_clusters.tsv` | Cluster summaries from valid deletion-gap records, including both `Pass` and `LowSupport` clusters. |

`split_reads.tsv` starts with the 16 raw alignment fields: SampleID, Read_Name, Flag, Primary_Chrom/Start/End/CIGAR/MAPQ/Strand, and SA_Chrom/Start/End/CIGAR/MAPQ/Strand/NM. Step 05 appends Primary/SA segment roles, breakpoints, the inferred deletion interval, `Junction_Type`, `Breakpoint_Status`, and `Breakpoint_Reason`. Multiple-SA reads are excluded during extraction. All extracted rows are retained during annotation. The file has no `Cluster_ID`; S3 clustering writes only the summary table.

### `S4_DiscordantReads/`

| File | Description |
|---|---|
| `discordant_pairs.tsv` | Filtered read-pair records passing the existing Discordant criteria. Contains the existing inner-gap `Start`, `End`, and `Size_bp`, but no cluster fields. |
| `discordant_clusters.tsv` | Cluster summaries, including `Pass` and `LowSupport` clusters. |
| `tlen_summary.tsv` | TLEN summary statistics when `RUN_TLEN_QC=1`. |
| `tlen_distribution.pdf` | TLEN distribution plot when `RUN_TLEN_QC=1`. |

The two processing stages are `ExtractDiscordantPairs.pl` and `ClusterDiscordantReads.pl`. BAM parsing and the existing Discordant filters are performed in one pass. The filtering thresholds, FR-orientation option, inner-gap interval calculation, cluster-distance rule, and minimum support rule are unchanged. The cluster summary stores the observed coordinate ranges. The plotting script reconstructs pair-to-cluster links from these ranges, so no assignment table is required.

### `S5_Candidates/`

The candidate stage writes three evidence-specific high-quality candidate tables and one split-centered integration table.

| File | Selection rule and description |
|---|---|
| `depth_candidates.tsv` | CUSUM events with `Type=DEL`. CUSUM event-generation thresholds are applied upstream in S2. |
| `split_candidates.tsv` | Split-read clusters with `Cluster_Status=Pass`. LowSupport clusters remain available in S3 but are excluded here. |
| `discordant_candidates.tsv` | Discordant-read clusters with `Cluster_Status=Pass`. LowSupport clusters remain available in S4 but are excluded here. |
| `merged_candidates.tsv` | All passed split-centered candidates with depth and discordant support metrics and an evidence level. |

`Candidate_ID` is assigned independently within each evidence-specific table. Split and discordant source clusters are traceable through `Source_Cluster_ID`; each depth event is identified by its `Candidate_ID`, `Depth_File`, and genomic coordinates. The integrated table receives a separate `CANDIDATE_######` identifier and retains `Split_Cluster_ID` as the source link.

S5 performs evidence filtering, interval overlap assessment, and split-centered candidate integration only. Gene, transcript, and exon annotation are deferred to S6, where the genome-build-specific MANE gene and exon tables are used consistently.

### `S6_Report/`

| File | Description |
|---|---|
| `annotated_candidates.tsv` | Compact final report containing MANE gene/transcript annotation, exon overlap summaries, and the integrated sequencing evidence. High, Moderate and Low split-centered candidates are retained. |

---

## Candidate Interpretation

### Coordinates

Candidate coordinates are reported as **1-based closed intervals**. Candidate size is calculated as:

```text
End - Start + 1
```

The split-read cluster defines the candidate interval. Depth and discordant evidence are used only as supporting validation; they do not replace the split-read coordinates.

### Important merged-candidate columns

| Column | Description |
|---|---|
| `Candidate_ID` | Sequential integrated candidate identifier within the sample. |
| `Split_Cluster_ID` | Source split-read cluster identifier. |
| `Chrom`, `Start`, `End`, `Size_bp` | Split-read-defined candidate interval and its length. |
| `Split_Reads` | Number of unique split-read names supporting the cluster. |
| `Depth_Support` | `Yes` when the union of qualifying CUSUM DEL intervals covers at least `EVIDENCE_OVERLAP_FRACTION` of the split interval. |
| `Depth_Covered_Bases` | Number of split-interval bases covered by the union of overlapping CUSUM intervals. |
| `Depth_Coverage_Fraction` | Covered fraction of the split interval. |
| `Depth_Candidate_Count` | Number of CUSUM DEL intervals overlapping the split interval. |
| `Discordant_Support` | `Yes` when at least one passing discordant cluster individually covers at least `EVIDENCE_OVERLAP_FRACTION` of the split interval. |
| `Discordant_Pairs` | Sum of `Support_Pairs` from qualifying overlapping passing discordant clusters. |
| `Discordant_Overlap_Fraction` | Maximum overlap fraction among discordant clusters. |
| `Discordant_Cluster_Count` | Number of passing discordant clusters meeting the overlap threshold. |
| `Median_Insert_Size` | Unique median insert-size values from qualifying discordant clusters. |
| `Evidence_Types` | Comma-separated supporting evidence types. |
| `Evidence_Level` | `High`, `Moderate`, or `Low`, using the unchanged split-centered rule. |

The evidence-level rule is unchanged: `High` = Split + Depth + Discordant; `Moderate` = Split + Depth; `Low` = Split without sufficient Depth support, even when Discordant support is present. A discordant-only or depth-only event does not independently create a row in `merged_candidates.tsv`.

### Final report annotation columns

`annotated_candidates.tsv` contains 20 columns:

| Column | Description |
|---|---|
| `SampleID`, `Candidate_ID` | Sample and integrated candidate identifiers. |
| `Gene` | MANE gene or genes overlapping the candidate; nearest gene is reported for a flanking event. |
| `Gene_Classification` | Evidence category from `HCM_CORE_GENE_LIST`, when available. |
| `Transcript` | MANE RefSeq transcript associated with the annotated gene. |
| `Region_Type` | `Exonic`, `Intronic`, `Flanking`, or `Intergenic`. |
| `Chrom`, `Start`, `End`, `Size_bp` | Candidate interval and length, using 1-based closed coordinates. |
| `Affected_Exons` | Exons meeting `MIN_EXON_OVERLAP_FRACTION`, formatted as `EX12,EX13,EX14`. |
| `Fully_Covered_Exons` | Affected exons completely contained in the candidate interval. |
| `Split_Reads` | Number of unique split reads supporting the candidate. |
| `Depth_Support`, `Depth_Coverage_Fraction` | Depth evidence status and covered fraction. |
| `Discordant_Support`, `Discordant_Pairs` | Discordant evidence status and supporting read-pair count. |
| `Discordant_Overlap_Fraction` | Maximum overlap fraction among passing discordant clusters. |
| `Median_Insert_Size` | Median insert size from supporting discordant clusters; `NA` without discordant support. |
| `Evidence_Level` | Unchanged `High`, `Moderate`, or `Low` classification from S5. |

`Affected_Exon_Count`, `Partially_Overlapped_Exons`, `Depth_Covered_Bases`, and `Evidence_Types` are omitted from the final report because they are derivable or redundant. `Region_Type=Exonic` means that the candidate overlaps at least one MANE exon; the affected-exon list applies the configured overlap-fraction threshold. `Intronic` means overlap with a gene interval but no exon interval. `Flanking` means outside the gene body but within `TARGET_REGION_FLANK` bp of the nearest gene. Other events are `Intergenic`.

A `High` evidence level describes sequencing evidence for a split-centered event. It does not by itself establish coding-exon involvement or pathogenicity. Interpret `Evidence_Level` together with `Region_Type` and exon-overlap columns.

---

## Test Example: 25B09089386

The repository contains an archived example result directory. This result was generated with the historical 35-gene list now retained as `db/Genes/HCMCoreGenesLegacy.tsv`. The current example configuration defaults to `db/Genes/HCMClassic8Genes.tsv`, so the archived result does not represent a fresh run with the current default gene list. It is retained for historical reference only; its S2, S5 and S6 files use legacy output schemas and must not be treated as current-format examples.

The example result directory is:

```text
test/test_results/25B09089386/
```

### Target-region extraction

According to `S1_TargetBam/25B09089386.target_bam.summary.tsv`:

| Metric | Result |
|---|---:|
| Core genes in list | 35 |
| Core genes found in MANE gene annotation | 34 |
| Core genes skipped | 1 |
| Skipped gene | `MT-TI` |
| Raw gene regions | 34 |
| Merged target BED regions | 34 |
| Target-region flank | 5,000 bp |

The skipped gene is reported in the summary rather than causing the sample workflow to terminate.

### Merged candidates

`S5_Candidates/25B09089386.merged_candidates.tsv` is the archived pre-refactor filename and contains seven historical split-centered candidates. It is not an example of the current S5 column schema:

| Evidence level | Count |
|---|---:|
| `High` | 3 |
| `Moderate` | 2 |
| `Low` | 2 |
| **Total** | **7** |

### Exon annotation

Only one of the seven candidate intervals meets the configured MANE exon-overlap requirement:

| Field | Result |
|---|---|
| Sample | `25B09089386` |
| Gene | `FHOD3` |
| Transcript | `NM_001281740` |
| Region | `chr18:34232240-34241309` |
| Size | 9,070 bp |
| Split-read support | 16 unique reads |
| Depth support | Yes |
| Depth coverage fraction | 0.9994 |
| Discordant support | Yes |
| Evidence level | `High` |
| Affected exons | `EX12`, `EX13`, `EX14` |
| Fully covered exons | `EX12`, `EX13`, `EX14` |

The other six intervals are retained as sequencing-evidence candidates but are reported as having no qualifying exon overlap. This distinction is intentional:

```text
Evidence strength != exon consequence
```

The test dataset demonstrates expected file generation and candidate prioritization. It should not be interpreted as a formal sensitivity, specificity, or clinical-performance validation dataset.

---

## Manual Review and Validation

All candidates should be reviewed before biological or clinical interpretation. See [`docs/integration_audit.md`](docs/integration_audit.md) for the code-level review and [`docs/release_checklist.md`](docs/release_checklist.md) for the pre-release verification items and remaining validation requirements.

Recommended review procedure:

1. Open the original WGS BAM and/or `SAMPLE.target.bam` in IGV.
2. Inspect depth across the candidate interval and adjacent exons.
3. Confirm that the depth decrease is consistent across the proposed deletion.
4. Inspect split reads and supplementary alignments near both breakpoints.
5. Inspect discordant read pairs spanning the interval.
6. Check mapping quality, duplicate status, repetitive sequence, segmental duplications, and nearby alignment artifacts.
7. Confirm that chromosome names and reference coordinates match the annotation build.
8. Compare affected and unaffected family members when available.
9. Prioritize candidates affecting coding exons and segregating with the phenotype.

Possible orthogonal validation methods include:

- breakpoint PCR;
- Sanger sequencing across the breakpoint;
- ddPCR;
- MLPA;
- qPCR;
- RNA-level confirmation when relevant and biologically interpretable;
- family segregation analysis.

---

## Running Individual Modules

The recommended approach is to generate and execute `run_SAMPLE.sh`. The following commands illustrate the current module interfaces. Set `PROJECT_ROOT` to the absolute repository path first. Downstream modules should normally use `SAMPLE.target.bam`, not the original WGS BAM.

```bash
PROJECT_ROOT="$(pwd)"
```

### Target BAM extraction

```bash
perl "$PROJECT_ROOT/bin/ExtractTargetBAM.pl" \
  --Config "$PROJECT_ROOT/conf/hcm_exondel.conf" \
  --GeneAnnotation "$PROJECT_ROOT/db/Annotation/hg19/MANERefSeqGenes.tsv" \
  --Sample SAMPLE \
  --Bam /absolute/path/to/SAMPLE.bam \
  --OutDir "$PROJECT_ROOT/results/SAMPLE/S1_TargetBam"
```

### Window-depth analysis

```bash
perl "$PROJECT_ROOT/bin/RunGeneMeanDepth.pl" \
  --Config "$PROJECT_ROOT/conf/hcm_exondel.conf" \
  --GeneAnnotation "$PROJECT_ROOT/db/Annotation/hg19/MANERefSeqGenes.tsv" \
  --Bam "$PROJECT_ROOT/results/SAMPLE/S1_TargetBam/SAMPLE.target.bam" \
  --Sample SAMPLE \
  --Out "$PROJECT_ROOT/results/SAMPLE/S2_Depth/window_candidates.tsv"
```

### Sample-level CUSUM analysis

```bash
perl "$PROJECT_ROOT/bin/RunSampleCUSUMDepth.pl" \
  --Config "$PROJECT_ROOT/conf/hcm_exondel.conf" \
  --InDir "$PROJECT_ROOT/results/SAMPLE/S2_Depth/gene_depth" \
  --GeneSummary "$PROJECT_ROOT/results/SAMPLE/S2_Depth/gene_mean_depth.tsv" \
  --Out "$PROJECT_ROOT/results/SAMPLE/S2_Depth/cusum_events.tsv" \
  --Perl /usr/bin/perl \
  --Script "$PROJECT_ROOT/bin/CUSUMDepthDel.pl"
```

### Depth-ratio plotting

```bash
python3 "$PROJECT_ROOT/bin/PlotDepthRatio.py" \
  --Input "$PROJECT_ROOT/results/SAMPLE/S2_Depth/window_depth_ratio.tsv" \
  --Output "$PROJECT_ROOT/results/SAMPLE/S2_Depth/depth_ratio.pdf" \
  --DepthRatioCutoff 0.65
```

The value passed to `--DepthRatioCutoff` must match `DEL_DEPTH_RATIO_CUTOFF` in the active configuration.

### Split-read extraction, deletion identification, and clustering

```bash
perl "$PROJECT_ROOT/bin/ExtractSASplitReads.pl" \
  --Config "$PROJECT_ROOT/conf/hcm_exondel.conf" \
  --Bam "$PROJECT_ROOT/results/SAMPLE/S1_TargetBam/SAMPLE.target.bam" \
  --Sample SAMPLE \
  --Out "$PROJECT_ROOT/results/SAMPLE/S3_SplitReads/split_reads.tsv"

perl "$PROJECT_ROOT/bin/IdentifyIntragenicDeletions.pl" \
  --Input "$PROJECT_ROOT/results/SAMPLE/S3_SplitReads/split_reads.tsv" \
  --Out "$PROJECT_ROOT/results/SAMPLE/S3_SplitReads/split_reads.tsv"

perl "$PROJECT_ROOT/bin/ClusterSASplitReads.pl" \
  --Config "$PROJECT_ROOT/conf/hcm_exondel.conf" \
  --Input "$PROJECT_ROOT/results/SAMPLE/S3_SplitReads/split_reads.tsv" \
  --ClustersOut "$PROJECT_ROOT/results/SAMPLE/S3_SplitReads/split_clusters.tsv"
```

### Discordant-pair extraction and clustering

```bash
perl "$PROJECT_ROOT/bin/ExtractDiscordantPairs.pl" \
  --Config "$PROJECT_ROOT/conf/hcm_exondel.conf" \
  --Bam "$PROJECT_ROOT/results/SAMPLE/S1_TargetBam/SAMPLE.target.bam" \
  --Sample SAMPLE \
  --Out "$PROJECT_ROOT/results/SAMPLE/S4_DiscordantReads/discordant_pairs.tsv"

perl "$PROJECT_ROOT/bin/ClusterDiscordantReads.pl" \
  --Config "$PROJECT_ROOT/conf/hcm_exondel.conf" \
  --Input "$PROJECT_ROOT/results/SAMPLE/S4_DiscordantReads/discordant_pairs.tsv" \
  --ClustersOut "$PROJECT_ROOT/results/SAMPLE/S4_DiscordantReads/discordant_clusters.tsv"
```

### Evidence merging

```bash
perl "$PROJECT_ROOT/bin/MergeEvidence.pl" \
  --Sample SAMPLE \
  --Depth "$PROJECT_ROOT/results/SAMPLE/S2_Depth/cusum_events.tsv" \
  --Split "$PROJECT_ROOT/results/SAMPLE/S3_SplitReads/split_clusters.tsv" \
  --Discordant "$PROJECT_ROOT/results/SAMPLE/S4_DiscordantReads/discordant_clusters.tsv" \
  --DepthCandidatesOut "$PROJECT_ROOT/results/SAMPLE/S5_Candidates/depth_candidates.tsv" \
  --SplitCandidatesOut "$PROJECT_ROOT/results/SAMPLE/S5_Candidates/split_candidates.tsv" \
  --DiscordantCandidatesOut "$PROJECT_ROOT/results/SAMPLE/S5_Candidates/discordant_candidates.tsv" \
  --MinOverlap 0.90 \
  --Out "$PROJECT_ROOT/results/SAMPLE/S5_Candidates/merged_candidates.tsv"
```

`MergeEvidence.pl` consumes the S2 CUSUM events, S3 split clusters, and S4 discordant clusters. It does not perform gene annotation; S6 uses the selected MANE reference tables for that purpose.

### Candidate annotation

```bash
perl "$PROJECT_ROOT/bin/AnnotateCandidates.pl" \
  --Config "$PROJECT_ROOT/conf/hcm_exondel.conf" \
  --GeneAnnotation "$PROJECT_ROOT/db/Annotation/hg19/MANERefSeqGenes.tsv" \
  --ExonAnnotation "$PROJECT_ROOT/db/Annotation/hg19/MANERefSeqExons.tsv" \
  --Input "$PROJECT_ROOT/results/SAMPLE/S5_Candidates/merged_candidates.tsv" \
  --Out "$PROJECT_ROOT/results/SAMPLE/S6_Report/annotated_candidates.tsv"
```

---

## Troubleshooting

### BAM index not found

```text
[ERROR] BAM index not found
```

Create an index:

```bash
samtools index sample.bam
```

### BAM index is older than the BAM

A stale index can cause incorrect regional extraction or warnings. Rebuild it:

```bash
rm -f sample.bam.bai sample.bai
samtools index sample.bam
```

### Invalid BAM list format

The BAM list must contain exactly two TAB-delimited columns. A sequence of spaces is not equivalent to a TAB.

Check hidden characters with:

```bash
cat -A input_bam.list
```

### Existing sample results

Generating `run_SAMPLE.sh` overwrites only the shell file. It does not clear existing analysis outputs. Remove stale sample result files before executing a rerun when input data, parameters, or software versions have changed.

### Deprecated parameter error

Remove either of these obsolete settings from the configuration:

```text
ANALYZE_CORE_GENES_ONLY
CANDIDATE_BAM_FLANK
```

### Core gene not found in the gene annotation

The target-BAM summary lists skipped genes. Check:

- the gene symbol spelling;
- whether the symbol is present in the MANE gene annotation;
- whether mitochondrial or non-MANE genes require a separate annotation strategy;
- whether the gene list contains aliases rather than the annotation's approved symbol.

A missing gene is skipped and reported; it does not invalidate successfully resolved genes.

### Chromosome naming mismatch

Examples of incompatible naming include:

```text
chr18 vs 18
chrM  vs MT
```

Ensure that the BAM and annotation files use the same reference build and chromosome naming convention.

### Per-gene depth files

Per-gene base-depth files are retained in `S2_Depth/gene_depth/` because the CUSUM step uses them. This directory is recreated at the start of each depth-analysis run, so old per-gene depth files are not reused.


### No split-read candidates

Because the current merge is split-read-centered, no passed split cluster means no merged candidate, even when depth or discordant signals exist.

Possible causes include:

- a true breakpoint in repetitive or poorly mappable sequence;
- missing or filtered supplementary alignments;
- insufficient breakpoint-spanning reads;
- `SA_SPLIT_MIN_SUPPORT_READS` being too strict;
- low mapping quality;
- breakpoint locations falling outside the extracted target region and flank.

Review depth and discordant diagnostic outputs separately when investigating such cases.

### Candidate remains `Low` despite depth overlap

Depth support requires coverage of at least `EVIDENCE_OVERLAP_FRACTION` of the split interval. Partial overlap below the threshold is recorded but does not qualify as depth support.

Review:

```text
Depth_Covered_Bases
Depth_Coverage_Fraction
Depth_Range
```

### Weak discordant support

Possible causes include:

- the deletion is smaller than `MIN_DISCORDANT_INSERT_SIZE`;
- the library insert-size distribution differs from the configured threshold;
- too few supporting read pairs;
- orientation filtering;
- mapping ambiguity around the breakpoints.

Use the TLEN summary as a diagnostic guide, then tune `MIN_DISCORDANT_INSERT_SIZE` explicitly when justified.

### Empty merged output

An empty merged file is a valid outcome when no split cluster passes filtering. The script should retain the header so downstream annotation can terminate cleanly.

Inspect:

```text
S3_SplitReads/split_reads.tsv
S3_SplitReads/split_clusters.tsv
logs/04.extract_split_reads.log
logs/05.cluster_split_reads.log
```

### High-evidence candidate has no exon annotation

`High` means that split, depth, and discordant evidence support the genomic event. It does not guarantee exon overlap.

Inspect:

```text
Gene
Transcript
Region_Type
Affected_Exons
Fully_Covered_Exons
```

Such a candidate may be intronic, flanking, intergenic, below the exon-overlap threshold, affected by annotation/reference mismatch, or artifactual.

---

## Limitations

- HCMExonDel is designed for WGS BAM data and is not optimized for WES or targeted-panel data.
- The pipeline analyzes only genes listed in `HCM_CORE_GENE_LIST`.
- The final merge is split-read-centered; events without a passed split-read cluster are not reported as merged candidates.
- Breakpoints in repetitive, low-complexity, or poorly mappable regions may lack reliable split-read support.
- Read-depth analysis is affected by sequencing depth, GC bias, mappability, alignment quality, and local coverage variability.
- Small deletions may not produce sufficient depth or discordant-pair evidence.
- Large deletions may extend beyond the target gene plus configured flank, reducing breakpoint evidence in the target BAM.
- A fixed `MIN_DISCORDANT_INSERT_SIZE` may not be optimal for every sequencing library.
- Exon annotation depends on the selected MANE transcript file and `MIN_EXON_OVERLAP_FRACTION`.
- Evidence level is not equivalent to pathogenicity, clinical classification, or exon consequence.
- The included test sample demonstrates workflow behavior but does not establish analytical sensitivity, specificity, precision, or limit of detection.
- All candidates require manual review and independent experimental validation.

---

## Citation

A formal publication citation is not yet available. When using HCMExonDel in a report or manuscript, cite the repository and describe the analysis version and configuration parameters.

Suggested methods text:

```text
Candidate exon-level deletions in hypertrophic cardiomyopathy-related core genes were detected from whole-genome sequencing BAM files using HCMExonDel. The pipeline generated a target BAM spanning predefined core genes and flanking regions, identified split-read clusters as candidate event backbones, evaluated supporting base-level CUSUM depth and discordant read-pair evidence, and annotated candidate intervals against MANE RefSeq exon coordinates. All candidates were manually reviewed and required independent validation.
```

Also report, at minimum:

- repository commit or release;
- reference genome build;
- MANE annotation version;
- HCM core-gene list version;
- principal depth, split-read, discordant-read, and overlap thresholds.

---

## License

No license file is currently included in the repository. Add an explicit license before distributing or reusing the software beyond the intended project context.

---

## Contact

Maintainer: `Xiaohuaniu0032`

Repository:

```text
https://github.com/Xiaohuaniu0032/HCMExonDel
```

---

## Disclaimer

HCMExonDel is provided for research use only. It is not intended for direct clinical diagnosis. Candidate deletion events must be reviewed by qualified personnel and validated using independent methods before clinical interpretation or reporting.


