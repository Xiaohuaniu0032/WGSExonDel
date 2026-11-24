# HCMExonDel Integration Audit

## Scope

This review uses the Phase 7 repository as its starting point. It checks the
main pipeline generator, configuration handling, module interfaces, result
paths, repository documentation and release prerequisites. No analysis
algorithm or candidate evidence threshold was intentionally changed in this
phase. No automated test framework was added.

## Corrections made

1. Fixed the physical line counter in the main configuration parser. Invalid
   configuration entries now report their actual line number.
2. Fixed the BAM-list line counter, which previously advanced twice per input
   line.
3. Restricted sample identifiers to begin with an alphanumeric character.
   This still accepts sample names beginning with digits and prevents special
   path components such as `.` and `..` from being used as sample directories.
4. Moved `run_info.tsv` creation from `HCMExonDel.pl` into the generated sample
   shell. The main script now creates the sample directory and writes or
   overwrites only `run_<SampleID>.sh`; runtime metadata is written when that
   shell begins execution.
5. Expanded runtime metadata to include the selected genome build, MANE gene
   and exon annotation paths, HCM gene list, all configuration values,
   executable paths and runtime Perl, Python and samtools versions.
6. Corrected the README configuration-copy instructions. The example template
   is copied to the ignored local file `conf/hcm_exondel.conf` rather than
   copied onto itself.
7. Corrected the documented Step 01 and Step 08 log filenames to match the
   generated shell.
8. Removed duplicate S3/S5 headings and repeated text in
   `docs/output_format.md`.
9. Renamed the example shell-generation helper to
   `test/GenerateExampleRunShells.sh` to follow the project script naming
   convention. It is a helper for generating per-sample run shells, not an
   automated test suite.
10. Added a pre-release checklist covering real-data validation, reference
    annotation provenance, the default gene-list scope, licensing and citation
    metadata.

## Static and manual checks

- Perl syntax: 30 scripts passed.
- Python syntax: 4 scripts passed.
- Bash syntax: 2 scripts passed.
- Main-script help output: passed.
- Python plotting-script help output: passed.
- Main generator smoke check with a synthetic BAM placeholder: passed.
- Generated sample shell syntax: passed.
- hg38 annotation path selection in the generated shell: passed.
- `run_info.tsv` absent immediately after shell generation and created when
  the generated shell starts: passed.
- Runtime Perl, Python and samtools version capture: passed.
- Sample IDs beginning with digits: accepted and produce `run_<SampleID>.sh`.
- Unsafe sample ID `..`: rejected.
- Invalid configuration line: reported with the correct physical line number.

The synthetic BAM placeholder was deliberately invalid for alignment
processing. Execution was stopped at the first BAM-dependent step after
verifying that runtime metadata was written. This was not a full pipeline run.

## Interface review

The main script's generated commands were checked against the current command
line option declarations of the core Perl and Python modules. The S1-S6 paths
and expected output filenames are consistent with the current main script and
`docs/output_format.md`.

All active configuration keys in `hcm_exondel.example.conf` were checked for
references in the main pipeline or the corresponding analysis module. No
obviously unused active configuration key was identified.

## Not yet verified

- End-to-end analysis with real coordinate-sorted BAM files.
- Real-data concordance of hg19 and hg38 annotation paths.
- Sensitivity, specificity, precision and reproducibility on a validation set.
- Agreement of detected FHOD3 deletion coordinates and affected exons with
  orthogonal validation data.
- The exact RefSeq GTF release/checksum used to generate both MANE coordinate
  tables.
- Whether the nine-record `HCMClassic8Genes.tsv` default list is the intended
  release scope. It currently contains eight classical sarcomeric genes plus
  FHOD3.
- Repository license and citation author metadata.

The code-level integration review is complete, but the repository should not
be described as analytically validated or publication-ready until the
real-BAM items above are completed.

## Post-audit hotfix: ExtractTargetBAM integer validation

- Added the missing `validate_positive_integer()` helper to `bin/ExtractTargetBAM.pl`.
- `TARGET_BAM_THREADS` is now checked as a positive integer before BAM extraction.
- No extraction, interval, or evidence-selection logic was changed.
