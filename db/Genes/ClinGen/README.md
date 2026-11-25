# ClinGen Gene Lists

This directory stores source data downloaded from ClinGen Expert Panels and manually prepared, analysis-ready gene lists.

## Hereditary Cardiovascular Disease Gene Curation Expert Panel

**Source:** [ClinGen Hereditary Cardiovascular Disease Gene Curation Expert Panel](https://www.clinicalgenome.org/affiliation/40104/)

`tableExport.csv` is the original table downloaded from the ClinGen affiliation page. It contains gene-disease curation records, including the gene, disease, mode of inheritance (MOI), SOP, contribution, classification, and last evaluation date.

The CSV is a source/reference file. Prepare a separate gene-list TSV from this downloaded table before using the data in downstream analyses.

### Preparing the gene-list file

Use the following three-column, tab-delimited format:

| Column | Description |
| --- | --- |
| `Gene` | Gene symbol only |
| `Evidence` | ClinGen `Classification` value |
| `Note` | Disease name and MONDO ID, MOI, SOP, contribution, and last evaluation date |

Important points when preparing the file:

- In the downloaded CSV, the `Gene` field may contain the HGNC identifier appended directly to the gene symbol, for example, `ABCC9HGNC:60`. Extract the gene symbol and remove the appended HGNC identifier from the `Gene` value.
- Retain the original ClinGen classification in the `Evidence` column. Do not silently convert or combine different evidence classifications.
- Include the disease name and MONDO ID, together with the available curation details, in `Note`.
- Keep separate records when the same gene has different associated diseases, modes of inheritance, or curation classifications. Do not deduplicate records by gene symbol alone.

Example:

```tsv
Gene	Evidence	Note
ABCC9	Disputed	Disease: Brugada syndrome (MONDO:0015263); MOI: AD; SOP: SOP11; Contribution: Primary; Last Eval: 10/28/2025
```

## Preparing gene lists for other diseases

Gene lists for other disease areas can be prepared in the same way:

1. Identify the relevant ClinGen Expert Panel or gene-disease curation page.
2. Download the corresponding gene curation table from the ClinGen website and retain the original downloaded file as source data.
3. Prepare a separate TSV gene list using the standardized `Gene`, `Evidence`, and `Note` columns described above.
4. Record the source ClinGen page in the corresponding documentation.
5. Preserve separate gene-disease curation records when a gene is associated with multiple diseases or inheritance modes.

Use descriptive filenames for additional source exports and prepared gene lists so that the disease area or Expert Panel can be identified. The original ClinGen export should be kept separate from the manually prepared TSV file.
