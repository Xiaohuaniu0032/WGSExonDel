# Database Resources

This directory contains the gene lists and reference annotation tables used by WGSExonDel.

## Directory layout

- `Genes/`: HCM gene lists and phenotype-associated gene resources, including HPO and ClinGen data.
- `Annotation/hg19/`: MANE Select gene and exon coordinates for hg19.
- `Annotation/hg38/`: MANE Select gene and exon coordinates for hg38.
- `MANE/`: source MANE summary and the derived RefSeq MANE Select transcript table.

Gene lists and annotation tables use tab-delimited text. Coordinates in the annotation tables are 1-based, closed intervals. The original ClinGen export (`tableExport.csv`) is a CSV file.

Reference-building scripts are maintained under `scripts/BuildReference/`. The RefSeq GTF files are not included in the repository. Download the matching hg19 or hg38 `ncbiRefSeq.gtf.gz` file from UCSC before rebuilding the annotation tables; download links and instructions are provided in `Annotation/README.md`.

## Default HCM gene list

`Genes/HCMClassic8Genes.tsv` is the default gene list specified by `HCM_CORE_GENE_LIST` in the example configuration. Despite the filename, it currently contains nine gene records: eight classical sarcomeric HCM genes plus `FHOD3`.

Additional HPO-associated gene lists and ClinGen gene-disease curation resources are provided under `Genes/`. See the corresponding subdirectory `README.md` files for their sources and preparation details.
