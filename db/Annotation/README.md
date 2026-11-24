# MANE RefSeq Coordinate Tables

The annotation tables are organized by reference genome build.

| Genome build | Gene coordinates | Exon coordinates | Unmatched transcripts |
|---|---|---|---|
| hg19 | `hg19/MANERefSeqGenes.tsv` | `hg19/MANERefSeqExons.tsv` | `hg19/UnmatchedTranscripts.tsv` |
| hg38 | `hg38/MANERefSeqGenes.tsv` | `hg38/MANERefSeqExons.tsv` | `hg38/UnmatchedTranscripts.tsv` |

All coordinate tables are tab-delimited and use 1-based, closed genomic coordinates. Their first-line headers are used by the pipeline for input validation.

To rebuild the tables, use `scripts/BuildReference/BuildMANERefSeqCoordinates.pl` with a RefSeq GTF matching the selected genome build and the MANE Select transcript table in `db/MANE/`. The GTF files are not distributed with this repository; supply the appropriate reference annotation file separately.
