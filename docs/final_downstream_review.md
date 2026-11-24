# Final downstream workflow review: S3-S6

## Current interface contract

1. S3 raw SA extraction writes one row per primary alignment with exactly one valid SA entry to `split_reads.tsv`.
2. `IdentifyIntragenicDeletions.pl` annotates `split_reads.tsv` in place. It retains every extracted row and appends segment roles, breakpoints, deletion coordinates, junction type, breakpoint status, and breakpoint reason. It uses a temporary file before replacing the input.
3. `ClusterSASplitReads.pl` reads the annotated `split_reads.tsv` directly. Only rows with `Breakpoint_Status=PASS`, `Junction_Type=Deletion_Gap`, and `Breakpoint_Reason=Valid_Deletion_Gap` enter clustering. It writes only `split_clusters.tsv`; no event or assignment TSV is generated.
4. S4 Discordant processing parses the target BAM directly, applies the existing SAM flag/MAPQ/duplicate filters, same-chromosome check, positive-TLEN selection, minimum insert-size threshold, optional FR orientation, inner-gap interval calculation, and deduplication in one pass. The existing thresholds and filtering rules are retained.
5. `ClusterDiscordantReads.pl` reads `discordant_pairs.tsv` and writes only `discordant_clusters.tsv`. The coordinate-distance clustering algorithm and minimum unique-pair support rule are unchanged. No assignment TSV is generated.
6. The plotting script reconstructs split-read and discordant-pair cluster membership from the observed coordinate ranges in the respective cluster summaries. If ranges overlap, it selects the nearest median interval.
7. S5 remains split-read-centered. Depth and Discordant intervals are compared with each passing Split interval using the existing overlap-fraction rule.
8. The `Evidence_Level` rule is unchanged:
   - High: Split + Depth + Discordant.
   - Moderate: Split + Depth, without qualifying Discordant support.
   - Low: Split without sufficient Depth support, regardless of Discordant support.
9. S5 does not perform provisional gene/transcript assignment. S6 performs the MANE-based gene, transcript, region, and exon annotation.
10. S6 writes the final annotated candidate report.

## Final S3 outputs

```text
S3_SplitReads/
├── split_reads.tsv
└── split_clusters.tsv
```

## Final S4 outputs

```text
S4_DiscordantReads/
├── discordant_pairs.tsv
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

## Validation status

This interface synchronization has been reviewed against the current S3/S4 script arguments and the S5/S6 input schemas. A complete rerun of sample `25B09089386` still requires access to the original BAM and the configured reference/annotation files in the HPC environment.

Before treating the sample as successfully rerun, perform:

- Perl syntax checks for `WGSExonDel.pl` and all Perl modules.
- Python bytecode compilation for `PlotHCMExonDelMultiTrack.py`.
- S3 synthetic-data test covering valid deletion gaps, filtered reasons, clustering, and header-only output.
- S4 synthetic-data test covering multiple clusters, coordinate ranges, minimum support, and header-only output.
- S5/S6 integration checks confirming the existing candidate schema and unchanged `Evidence_Level` rule.
- Full `25B09089386` rerun on HPC, followed by comparison of the FHOD3 candidate coordinates and supporting evidence with the previously reviewed result.

These checks are distinct from clinical validation of candidate intervals.
