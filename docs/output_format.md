# HCMExonDel output format

## S2_Depth

The standard depth-analysis output is organized as follows:

```text
S2_Depth/
├── gene_mean_depth.tsv
├── window_depth_ratio.tsv
├── window_candidates.tsv
├── gene_depth/
│   └── *.depth.tsv
├── cusum_events.tsv
└── depth_ratio.pdf
```

All genomic coordinates in the TSV files use 1-based, closed intervals.
The per-gene depth files are the direct `samtools depth` output and are
retained because the CUSUM step consumes them. Per-gene CUSUM intermediate
files are temporary and are removed after `cusum_events.tsv` is written.

### `gene_mean_depth.tsv`

One row per analyzed MANE gene/transcript.

| Column | Description |
|---|---|
| `SampleID` | Sample identifier. |
| `Gene` | HCM core gene symbol. |
| `Classification` | Evidence classification from the configured HCM gene list. |
| `Transcript` | MANE Select RefSeq transcript accession. |
| `Chrom` | Chromosome. |
| `Start`, `End` | Gene interval. |
| `Strand` | Gene strand. |
| `Exon_Count` | Number of exons in the selected MANE transcript. |
| `Gene_Length` | Gene interval length in bp. |
| `Gene_Mean_Depth` | Mean per-base depth across the gene interval. |
| `Window_Count` | Number of windows generated for the gene. |
| `Depth_File` | Basename of the corresponding file in `gene_depth/`. |

### `window_depth_ratio.tsv`

One row per generated window. This table consolidates the former window-depth,
all-window-ratio and deletion-window tables. It includes normal, deletion-like,
low-gene-depth and non-evaluable windows.

| Column | Description |
|---|---|
| `SampleID` | Sample identifier. |
| `Gene` | Gene symbol. |
| `Classification` | HCM gene evidence classification. |
| `Transcript` | MANE Select transcript. |
| `Window_ID` | Window identifier within the gene. |
| `Chrom` | Chromosome. |
| `Start`, `End` | Window interval. |
| `Window_Length` | Window length in bp. |
| `Gene_Mean_Depth` | Mean depth across the gene. |
| `Window_Mean_Depth` | Mean depth across the window. |
| `Depth_Ratio` | `Window_Mean_Depth / Gene_Mean_Depth`; `NA` when not calculable. |
| `Window_Status` | `Normal_window`, `Del_window`, `Low_gene_depth` or `NA`. |

`Is_Del_Window` and free-text `Comment` are omitted because their information
is represented by `Window_Status` or the documented parameter rules.

### `window_candidates.tsv`

Consecutive deletion-supporting windows are merged using the existing
`MIN_CONSECUTIVE_DEL_WINDOWS` rule. This is an independent window-level result
and is not the depth input to the current split-read-centered evidence merge.

| Column | Description |
|---|---|
| `SampleID` | Sample identifier. |
| `Candidate_ID` | Sequential depth-candidate identifier within the sample. |
| `Gene` | Gene symbol. |
| `Classification` | HCM gene evidence classification. |
| `Transcript` | MANE Select transcript. |
| `Chrom` | Chromosome. |
| `Start`, `End` | Merged candidate interval. |
| `Length_bp` | Candidate interval length in bp. |
| `Window_Count` | Number of consecutive deletion-supporting windows. |
| `Gene_Mean_Depth` | Mean depth across the gene. |
| `Candidate_Mean_Depth` | Length-weighted mean depth across candidate windows. |
| `Mean_Depth_Ratio` | Mean of candidate-window depth ratios. |
| `Min_Depth_Ratio`, `Max_Depth_Ratio` | Minimum and maximum candidate-window depth ratios. |

### `cusum_events.tsv`

The merged CUSUM event table is used by `MergeEvidence.pl` as depth evidence.
Per-gene CUSUM output files are temporary and are not retained.

| Column | Description |
|---|---|
| `SampleID` | Sample identifier. |
| `Gene` | Gene symbol mapped from `gene_mean_depth.tsv`. |
| `Classification` | HCM gene evidence classification. |
| `Transcript` | MANE Select transcript. |
| `Chrom` | Chromosome. |
| `Start`, `End` | CUSUM event interval. |
| `Length` | Event length in bp. |
| `N_Bins` | Number of bins supporting the event. |
| `Baseline` | CUSUM baseline depth. |
| `Mean_Depth` | Mean depth across the event. |
| `Mean_Ratio` | Mean event depth ratio relative to baseline. |
| `Min_Ratio` | Minimum depth ratio within the event. |
| `Max_CUSUM` | Maximum CUSUM statistic. |
| `Type` | Event type; the current deletion workflow uses `DEL`. |
| `Depth_File` | Basename of the source per-gene depth file. |

### `depth_ratio.pdf`

A multi-page PDF with one page per analyzed gene. The horizontal threshold line
is supplied from `DEL_DEPTH_RATIO_CUTOFF` in the configuration file; the
plotting script does not maintain an independent hard-coded threshold.


## S3_SplitReads

The split-read workflow has three separate stages: raw SA extraction, deletion-event identification, and event clustering.

```text
S3_SplitReads/
├── split_reads.tsv
├── split_deletion_events.tsv
├── split_cluster_assignments.tsv
└── split_clusters.tsv
```

All genomic coordinates use 1-based, closed intervals.

### `split_reads.tsv`

One row represents one primary alignment with exactly one valid SA alignment. Reads with multiple SA entries are excluded during extraction. This table contains raw alignment information only; it does not contain deletion coordinates, `Event_ID`, `Cluster_ID`, or `Cluster_Status`.

The 16 columns are:

```text
SampleID
Read_Name
Flag
Primary_Chrom
Primary_Start
Primary_End
Primary_CIGAR
Primary_MAPQ
Primary_Strand
SA_Chrom
SA_Start
SA_End
SA_CIGAR
SA_MAPQ
SA_Strand
SA_NM
```

`Primary_End` and `SA_End` are calculated from the corresponding alignment start and the reference-consuming CIGAR length.

### `split_deletion_events.tsv`

One row represents one read-level deletion event retained by `IdentifyIntragenicDeletions.pl`. The script applies these checks in order:

1. `Primary_Chrom` must equal `SA_Chrom`.
2. `Primary_Strand` must equal `SA_Strand`.
3. The Primary CIGAR must have a terminal or leading soft clip. A terminal `S` means the Primary alignment is treated as the left alignment and the SA alignment as the right alignment. A leading `S` means the SA alignment is treated as left and the Primary alignment as right.
4. The positive reference gap between the two alignment intervals is reported as the deletion interval.

For Primary-left events:

```text
Start = Primary_End + 1
End   = SA_Start - 1
```

For SA-left events:

```text
Start = SA_End + 1
End   = Primary_Start - 1
```

Rows with non-positive intervals are excluded. Despite the script name, no HCM target-gene overlap check, gene-list filtering, or transcript annotation is performed at this stage; gene annotation is deferred to S6.

Columns:

| Column | Description |
|---|---|
| `Event_ID` | Unique read-level deletion event identifier (`SR_EVENT_######`). |
| `SampleID`, `Read_Name` | Sample and read identifiers. |
| `Chrom`, `Start`, `End`, `Size_bp` | Inferred deletion interval and length. |
| `Strand` | Shared Primary/SA strand. |
| `Primary_Chrom`, `Primary_Start`, `Primary_End`, `Primary_CIGAR`, `Primary_MAPQ`, `Primary_Strand` | Primary alignment details. |
| `SA_Chrom`, `SA_Start`, `SA_End`, `SA_CIGAR`, `SA_MAPQ`, `SA_Strand`, `SA_NM` | SA alignment details. |

### `split_cluster_assignments.tsv`

One row maps one `Event_ID` to its cluster.

| Column | Description |
|---|---|
| `Event_ID` | Read-level deletion event identifier. |
| `Cluster_ID` | Assigned split-read cluster identifier (`SR_CLUSTER_######`). |
| `Cluster_Status` | `Pass` or `LowSupport`. |

### `split_clusters.tsv`

One row represents one cluster. Both `Pass` and `LowSupport` clusters are retained. S5 uses only `Pass` clusters.

| Column | Description |
|---|---|
| `SampleID` | Sample identifier. |
| `Cluster_ID` | Unique cluster identifier within the sample. |
| `Chrom` | Chromosome. |
| `Start`, `End`, `Size_bp` | Median cluster interval and its length. |
| `Support_Reads` | Number of unique read names supporting the cluster. |
| `Support_Records` | Number of unique read-level deletion events in the cluster. |
| `Start_Range`, `End_Range` | Minimum-to-maximum event boundaries. |
| `Mean_Size_bp` | Mean deletion-event length. |
| `Size_Range_bp` | Minimum-to-maximum event length. |
| `Primary_Strand_Count`, `SA_Strand_Count` | Strand counts among supporting events. |
| `Cluster_Status` | `Pass` when `Support_Reads` meets `SA_SPLIT_MIN_SUPPORT_READS`; otherwise `LowSupport`. |

## S4_DiscordantReads

The Discordant workflow parses the target BAM directly, applies the existing read filters in the same pass, and then performs clustering. The existing read filtering, deletion-like interval calculation, cluster-distance rule, and minimum-support rule are retained.

```text
S4_DiscordantReads/
├── discordant_pairs.tsv
├── discordant_cluster_assignments.tsv
├── discordant_clusters.tsv
├── tlen_summary.tsv              # when RUN_TLEN_QC=1
└── tlen_distribution.pdf         # when RUN_TLEN_QC=1
```

### `discordant_pairs.tsv`

One row represents one read pair retained after the existing Discordant filters. Only the positive-TLEN record is retained for each fragment. The existing inner-gap interval is calculated as:

```text
Left_Breakpoint  = Read_End
Right_Breakpoint = Mate_Start
Start            = Read_End + 1
End              = Mate_Start - 1
Size_bp          = End - Start + 1
```

This interval calculation is intentionally retained from the existing implementation. It is the Discordant interval used by the existing S5 overlap comparison; it is not recalculated from a TLEN baseline model.

Columns:

| Column | Description |
|---|---|
| `Pair_ID` | Unique retained pair identifier (`DR_PAIR_######`). |
| `SampleID`, `Read_Name` | Sample and pair query name. |
| `Flag` | SAM flag of the retained positive-TLEN alignment. |
| `Chrom` | Chromosome shared by read and mate. |
| `Read_Start`, `Read_End` | Retained read alignment interval. |
| `Mate_Chrom`, `Mate_Start` | Mate chromosome and alignment start. |
| `TLEN`, `Insert_Size` | Signed template length and its absolute value. |
| `MAPQ`, `CIGAR` | Retained alignment mapping quality and CIGAR. |
| `Read_Strand`, `Mate_Strand`, `Pair_Orientation` | Read/mate strand and pair orientation. |
| `Left_Breakpoint`, `Right_Breakpoint` | Inner boundaries used by the existing interval calculation. |
| `Start`, `End`, `Size_bp` | Existing inferred Discordant interval and length. |

No `Cluster_ID` or `Cluster_Status` is written into this file.

### `discordant_cluster_assignments.tsv`

One row maps one `Pair_ID` to its cluster.

| Column | Description |
|---|---|
| `Pair_ID` | Filtered read-pair identifier. |
| `Cluster_ID` | Assigned cluster identifier (`DR_CLUSTER_######`). |
| `Cluster_Status` | `Pass` or `LowSupport`. |

### `discordant_clusters.tsv`

One row represents one cluster. Both `Pass` and `LowSupport` clusters are retained. S5 uses only `Pass` clusters.

| Column | Description |
|---|---|
| `SampleID` | Sample identifier. |
| `Cluster_ID` | Cluster identifier. |
| `Chrom` | Chromosome. |
| `Start`, `End`, `Size_bp` | Median cluster interval and length, using the existing interval definition. |
| `Support_Pairs` | Number of unique read-pair names supporting the cluster. |
| `Start_Range`, `End_Range` | Minimum-to-maximum interval boundaries among supporting pairs. |
| `Median_Insert_Size`, `Mean_Insert_Size` | Median and mean absolute TLEN. |
| `Insert_Size_Range` | Minimum-to-maximum absolute TLEN. |
| `Pair_Orientation_Count` | Counts by pair orientation. |
| `Cluster_Status` | `Pass` when `Support_Pairs` meets `MIN_DISCORDANT_PAIRS`; otherwise `LowSupport`. |

### TLEN quality control

When `RUN_TLEN_QC=1`, the temporary positive-TLEN sample table is removed after successful plotting. The retained outputs are `tlen_summary.tsv` and `tlen_distribution.pdf`. The diagnostic robust threshold reported in the summary does not replace `MIN_DISCORDANT_INSERT_SIZE` and is not used to change the existing Discordant filtering logic.

## S5_Candidates

S5 writes three evidence-specific candidate tables and one split-centered integration table.

```text
S5_Candidates/
├── depth_candidates.tsv
├── split_candidates.tsv
├── discordant_candidates.tsv
└── merged_candidates.tsv
```

### Candidate selection

| File | Source and selection rule |
|---|---|
| `depth_candidates.tsv` | `cusum_events.tsv` records with `Type=DEL`. |
| `split_candidates.tsv` | `split_clusters.tsv` records with `Cluster_Status=Pass`. |
| `discordant_candidates.tsv` | `discordant_clusters.tsv` records with `Cluster_Status=Pass`. |
| `merged_candidates.tsv` | Each passing split cluster is a candidate backbone. Depth and Discordant evidence are assessed by the existing interval-overlap rule. |

Gene, transcript, and exon annotation are not performed in S5. They are handled by S6 using the selected MANE annotation files.

### `depth_candidates.tsv`

One row per retained CUSUM DEL event. Columns:

```text
SampleID, Candidate_ID, Chrom, Start, End, Size_bp,
N_Bins, Baseline_Depth, Mean_Depth, Mean_Depth_Ratio,
Min_Depth_Ratio, Max_CUSUM, Depth_File
```

`Candidate_ID` uses the `DEPTH_######` prefix.

### `split_candidates.tsv`

One row per passing split cluster. Columns:

```text
SampleID, Candidate_ID, Source_Cluster_ID, Chrom, Start, End,
Size_bp, Support_Reads
```

`Candidate_ID` uses the `SPLIT_######` prefix. `Source_Cluster_ID` retains the original S3 `Cluster_ID`.

### `discordant_candidates.tsv`

One row per passing Discordant cluster. Columns:

```text
SampleID, Candidate_ID, Source_Cluster_ID, Chrom, Start, End,
Size_bp, Support_Pairs, Median_Insert_Size, Mean_Insert_Size,
Insert_Size_Range, Pair_Orientation_Count
```

`Candidate_ID` uses the `DISCORDANT_######` prefix. `Source_Cluster_ID` retains the original S4 `Cluster_ID`.

### `merged_candidates.tsv`

One row per passing split-centered candidate. The split cluster defines the candidate coordinates. Depth and Discordant evidence only provide support and do not replace the split interval.

Columns:

| Column | Description |
|---|---|
| `SampleID`, `Candidate_ID` | Sample and integrated candidate identifier (`CANDIDATE_######`). |
| `Split_Cluster_ID` | Source S3 split cluster. |
| `Chrom`, `Start`, `End`, `Size_bp` | Split-read-defined interval. |
| `Split_Reads` | Number of unique supporting split reads. |
| `Depth_Support` | `Yes` when the union of overlapping CUSUM DEL intervals covers at least `EVIDENCE_OVERLAP_FRACTION` of the split interval. |
| `Depth_Covered_Bases`, `Depth_Coverage_Fraction` | Number and fraction of split-interval bases covered by depth evidence. |
| `Depth_Candidate_Count` | Number of CUSUM DEL intervals overlapping the split interval. |
| `Discordant_Support` | `Yes` when at least one passing Discordant cluster individually covers at least `EVIDENCE_OVERLAP_FRACTION` of the split interval. |
| `Discordant_Pairs` | Sum of supporting pairs from qualifying overlapping Discordant clusters. |
| `Discordant_Overlap_Fraction` | Maximum overlap fraction among passing Discordant clusters. |
| `Discordant_Cluster_Count` | Number of passing Discordant clusters meeting the overlap threshold. |
| `Median_Insert_Size` | Unique median insert-size values from qualifying Discordant clusters. |
| `Evidence_Types` | Comma-separated evidence types. |
| `Evidence_Level` | Existing `High`, `Moderate`, or `Low` classification. |

**Evidence_Level logic is unchanged:**

- `High`: Split + Depth + Discordant.
- `Moderate`: Split + Depth, without qualifying Discordant support.
- `Low`: Split without sufficient Depth support, whether or not Discordant support is present.

A Depth-only or Discordant-only event does not independently create a row in `merged_candidates.tsv`.

## S6_Report

`annotated_candidates.tsv` is the final report. It combines S5 evidence metrics with MANE gene, transcript, and exon annotation. It contains 20 columns in this order:

```text
SampleID
Candidate_ID
Gene
Gene_Classification
Transcript
Region_Type
Chrom
Start
End
Size_bp
Affected_Exons
Fully_Covered_Exons
Split_Reads
Depth_Support
Depth_Coverage_Fraction
Discordant_Support
Discordant_Pairs
Discordant_Overlap_Fraction
Median_Insert_Size
Evidence_Level
```

| Column | Description |
|---|---|
| `Gene` | MANE gene or genes overlapping the candidate; nearest gene is reported for a flanking event. |
| `Gene_Classification` | Classification from the configured HCM gene list. |
| `Transcript` | MANE RefSeq transcript associated with the annotated gene. |
| `Region_Type` | `Exonic`, `Intronic`, `Flanking`, or `Intergenic`. |
| `Affected_Exons` | Exons meeting `MIN_EXON_OVERLAP_FRACTION`, formatted as `EX12,EX13,EX14`. |
| `Fully_Covered_Exons` | Affected exons fully contained in the candidate interval. |
| `Depth_Support`, `Depth_Coverage_Fraction` | Depth support status and covered fraction. |
| `Discordant_Support`, `Discordant_Pairs` | Discordant support status and supporting pair count. |
| `Discordant_Overlap_Fraction` | Maximum overlap fraction among passing Discordant clusters. |
| `Median_Insert_Size` | Median insert size from qualifying Discordant clusters; `NA` without Discordant support. |
| `Evidence_Level` | The unchanged S5 evidence level. |

`Affected_Exon_Count`, `Partially_Overlapped_Exons`, `Depth_Covered_Bases`, and `Evidence_Types` are omitted from the final report as redundant or derivable. `Region_Type=Exonic` indicates any MANE exon overlap; the affected-exon list applies the configured overlap-fraction threshold. `Intronic` means overlap with a gene interval but no exon interval. `Flanking` means outside the gene body but within `TARGET_REGION_FLANK` bp of the nearest gene. Other events are `Intergenic`.
