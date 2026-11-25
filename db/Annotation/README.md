# MANE RefSeq Coordinate Tables

This directory stores gene-level and exon-level genomic coordinate tables for MANE Select RefSeq transcripts. The annotation tables are organized by reference genome build.

| Genome build | Gene coordinates | Exon coordinates | Unmatched transcripts |
|---|---|---|---|
| hg19 | `hg19/MANERefSeqGenes.tsv` | `hg19/MANERefSeqExons.tsv` | `hg19/UnmatchedTranscripts.tsv` |
| hg38 | `hg38/MANERefSeqGenes.tsv` | `hg38/MANERefSeqExons.tsv` | `hg38/UnmatchedTranscripts.tsv` |

All coordinate tables are tab-delimited and use 1-based, closed genomic coordinates. Their first-line headers are used by the pipeline for input validation.

## RefSeq GTF files

The RefSeq GTF files are available from the UCSC Genome Browser download server:

| Genome build | Download directory | GTF file |
|---|---|---|
| hg19 | [UCSC hg19 genes](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/) | [hg19.ncbiRefSeq.gtf.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.ncbiRefSeq.gtf.gz) |
| hg38 | [UCSC hg38 genes](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/) | [hg38.ncbiRefSeq.gtf.gz](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz) |

Download and decompress the GTF file for the required reference build. The script expects an uncompressed GTF file.

For example, download the hg38 file:

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz \
    -O db/Annotation/hg38/ncbiRefSeq.gtf.gz

gunzip -c db/Annotation/hg38/ncbiRefSeq.gtf.gz \
    > db/Annotation/hg38/ncbiRefSeq.gtf
```

For hg19, use the corresponding hg19 download URL and output directory.

## Build MANE RefSeq coordinate tables

Use `scripts/BuildReference/BuildMANERefSeqCoordinates.pl` to extract exon-level and gene-level coordinates for MANE Select RefSeq transcripts from the UCSC NCBI RefSeq GTF file.

The script also requires the MANE Select transcript table:

```text
db/MANE/MANERefSeqSelect.tsv
```

Run the following commands from the repository root.

### hg19

```bash
perl scripts/BuildReference/BuildMANERefSeqCoordinates.pl \
    --gtf db/Annotation/hg19/ncbiRefSeq.gtf \
    --mane db/MANE/MANERefSeqSelect.tsv \
    --outdir db/Annotation/hg19
```

### hg38

```bash
perl scripts/BuildReference/BuildMANERefSeqCoordinates.pl \
    --gtf db/Annotation/hg38/ncbiRefSeq.gtf \
    --mane db/MANE/MANERefSeqSelect.tsv \
    --outdir db/Annotation/hg38
```

### Output files

The script generates the following files in the selected output directory:

| File | Description |
|---|---|
| `MANERefSeqExons.tsv` | Exon-level coordinates for MANE Select transcripts |
| `MANERefSeqGenes.tsv` | Gene-level genomic coordinates and exon counts |
| `UnmatchedTranscripts.tsv` | MANE Select transcripts without matching exon records in the GTF |

The output coordinate convention is 1-based, closed, matching the input GTF.

### Optional arguments

The script supports the following optional arguments:

- `--exon-out`: Specify the output exon-coordinate file. Default: `<outdir>/MANERefSeqExons.tsv`.
- `--gene-out`: Specify the output gene-coordinate file. Default: `<outdir>/MANERefSeqGenes.tsv`.
- `--unmatched`: Specify the unmatched-transcript report. Default: `<outdir>/UnmatchedTranscripts.tsv`.
- `--keep-version`: Retain RefSeq transcript version suffixes. By default, transcript versions are removed (for example, `NM_001281740.3` becomes `NM_001281740`).
- `--help`: Display the usage information.

To view all available options:

```bash
perl scripts/BuildReference/BuildMANERefSeqCoordinates.pl --help
```
