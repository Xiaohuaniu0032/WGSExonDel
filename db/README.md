# Database Resources

This directory contains the gene lists and reference annotation tables used by HCMExonDel.

## Directory layout

- `Genes/`: HCM gene lists and phenotype-associated gene resources.
- `Annotation/hg19/`: MANE Select gene and exon coordinates for hg19.
- `Annotation/hg38/`: MANE Select gene and exon coordinates for hg38.
- `MANE/`: source MANE summary and the derived RefSeq MANE Select transcript table.

All tabular resources use tab-delimited text. Coordinates in the annotation tables are 1-based, closed intervals.

Reference-building scripts are maintained under `scripts/BuildReference/`. The GTF files used to build the coordinate tables are not bundled because the previous repository contained only broken symbolic links to an internal server. Obtain the matching RefSeq GTF for the required genome build before rebuilding these resources.

## Default HCM gene list

`Genes/HCMClassic8Genes.tsv` is the default list selected by the example
configuration. It currently contains eight classical sarcomeric HCM genes plus
FHOD3, for nine gene records in total. Confirm that this is the intended default
analysis scope before creating a public release. `Genes/HCMGenes.tsv` is a
broader HCM-associated resource and is not selected automatically by the main
pipeline.
