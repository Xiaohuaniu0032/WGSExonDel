# Final downstream workflow review: S3-S6

## Fixed design decisions

1. S3 split-read processing is separated into raw SA extraction, deletion-event identification, and clustering.
2. Reads with more than one SA entry are excluded during raw extraction.
3. `split_reads.tsv` contains only the 16 raw Primary/SA alignment fields and has no event or cluster identifiers.
4. `IdentifyIntragenicDeletions.pl` applies same-chromosome and same-strand checks, then uses the Primary CIGAR terminal/leading soft clip to determine the alignment order and calculate the positive reference gap.
5. S4 Discordant processing parses the target BAM directly, applies the existing SAM flag/MAPQ/duplicate filters, same-chromosome check, positive-TLEN selection, minimum insert-size threshold, optional FR orientation, inner-gap interval calculation, and deduplication in one pass, then performs clustering separately. Existing thresholds and clustering rules are retained. `ClusterDiscordantReads.pl` also snapshots records when closing each cluster so assignments remain linked to the correct `Pair_ID` across multiple clusters.
6. S4 cluster identifiers and statuses are written to a separate assignment table, not repeated in the complete read-pair table.
7. S5 remains split-read-centered. Depth and Discordant intervals are compared with each passing Split interval using the existing overlap-fraction rule.
8. The `Evidence_Level` rule is unchanged:
   - High: Split + Depth + Discordant.
   - Moderate: Split + Depth, without qualifying Discordant support.
   - Low: Split without sufficient Depth support, regardless of Discordant support.
9. S5 no longer performs provisional gene/transcript assignment. S6 performs the MANE-based gene, transcript, region, and exon annotation.
10. S6 writes a compact 20-column final report.

## Final S3 outputs

```text
S3_SplitReads/
├── split_reads.tsv
├── split_deletion_events.tsv
├── split_cluster_assignments.tsv
└── split_clusters.tsv
```

## Final S4 outputs

```text
S4_DiscordantReads/
├── discordant_pairs.tsv
├── discordant_cluster_assignments.tsv
├── discordant_clusters.tsv
├── tlen_summary.tsv
└── tlen_distribution.pdf
```

The TLEN files are generated when `RUN_TLEN_QC=1`.

## Final S5 outputs

```text
S5_Candidates/
├── depth_candidates.tsv
├── split_candidates.tsv
├── discordant_candidates.tsv
└── merged_candidates.tsv
```

## Final S6 output

```text
S6_Report/
└── annotated_candidates.tsv
```

## Verification performed

- Perl syntax check for the main script and all Perl scripts.
- Python bytecode compilation for the Python scripts.
- Synthetic SAM test for S3 raw extraction, multi-SA exclusion, deletion-event interval inference, event assignments, and cluster summary.
- Synthetic SAM parity tests confirming that the new one-pass `ExtractDiscordantPairs.pl` produces a `discordant_pairs.tsv` identical to the previous two-step implementation under the default settings, with FR-orientation filtering disabled, and with duplicate exclusion disabled.
- Synthetic SAM test for S4 current filtering criteria, Pair_ID assignment, cluster assignments, and cluster summary, including multiple clusters to verify each Pair_ID maps to its correct Cluster_ID. Cluster summary values were unchanged from the previous clustering implementation; the previously incorrect assignment-table reference mapping was corrected.
- Empty-input test confirming that S4 still writes header-only pair, assignment, and cluster tables.
- Synthetic S5 integration test covering High, Moderate, and Low; the observed categories matched the unchanged Evidence_Level rule.
- Synthetic S6 annotation test covering Exonic, Flanking, and Intergenic events, including concise affected-exon labels.
- Main-script shell generation test using an hg38 annotation selection; the generated sample shell passed `bash -n` and contained the revised S3-S6 commands.

These are code-level and synthetic-data checks. They do not replace a complete run on real WGS BAM data or clinical validation of the candidate intervals.
