# HCMExonDel Pre-release Checklist

This checklist separates code-level readiness from analytical validation. Static
checks do not replace evaluation with representative sequencing data.

## Code and repository review

- [x] Core pipeline scripts use the current S1-S6 directory layout.
- [x] Core Perl and Python script filenames use UpperCamelCase.
- [x] Core command-line options use UpperCamelCase.
- [x] The main configuration parser rejects unknown and duplicate parameters.
- [x] hg19 and hg38 MANE annotation paths are selected by `--GenomeBuild`.
- [x] Sample shell files use the `run_<SampleID>.sh` naming convention.
- [x] Sample shell generation overwrites the shell file without clearing existing result files.
- [x] Runtime metadata is written by the generated sample shell, not by the main shell generator.
- [x] README and output-format documentation describe the current output schema.
- [x] No automated test framework is included in this release-preparation scope.

## Required before public release

- [ ] Run the complete workflow with at least one representative hg19 BAM.
- [ ] Run the complete workflow with at least one representative hg38 BAM.
- [ ] Confirm the BAM reference build agrees with `--GenomeBuild` and the selected MANE tables.
- [ ] Review all S1-S6 outputs from real BAM runs, including empty-result cases.
- [ ] Recheck the known FHOD3 exon-deletion positive sample and compare the detected interval and exon annotation with the established validation results.
- [ ] Review unaffected family members and available negative controls for false-positive candidates.
- [ ] Confirm the final report fields and evidence levels against manually reviewed IGV evidence.
- [ ] Record the exact samtools, Perl, Python and Python-package versions used for the publication analysis.
- [ ] Record the exact MANE summary release and the RefSeq GTF source/version used to build both genome-build coordinate tables.
- [ ] Replace all placeholder absolute paths in the local configuration before analysis.
- [ ] Confirm the contents and intended scope of `HCMClassic8Genes.tsv`. The current file contains eight classical sarcomeric HCM genes plus FHOD3 (nine gene records); its name should be retained or changed only after confirming that this is the intended default gene set.
- [ ] Add a software license selected by the repository owner. No license is currently included.
- [ ] Add citation metadata (`CITATION.cff`) after confirming the preferred author list, maintainer details and publication/repository citation.
- [ ] Create a versioned Git tag and release notes after the real-data validation is complete.

## Suggested analytical review record

For each representative sample, retain the following information with the
analysis record:

- Sample identifier and BAM checksum.
- Reference build and BAM header contig information.
- Configuration file and checksum.
- MANE annotation table checksums and source versions.
- Software versions and execution date.
- Counts of read-level evidence, clusters, S5 candidates and S6 candidates.
- Manually reviewed candidate coordinates, affected exons and evidence level.
- Orthogonal validation result, when available.

HCMExonDel is a research-use candidate-detection workflow. A passing evidence
level is not a clinical classification and does not replace independent
validation.
