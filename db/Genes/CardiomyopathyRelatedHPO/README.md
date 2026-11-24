# HPO phenotype_to_genes.txt Download Guide

## 1. Data Source

`phenotype_to_genes.txt` is a gene--phenotype association file released
by the Human Phenotype Ontology (HPO) project. It provides associations
between HPO terms, genes, and related diseases.

-   **Data source:** [Human Phenotype Ontology
    (HPO)](https://github.com/obophenotype/human-phenotype-ontology)
-   **Version:** `v2026-09-01`
-   **Official release:**
    https://github.com/obophenotype/human-phenotype-ontology/releases/tag/v2026-09-01
-   **Official file documentation:**
    https://github.com/obophenotype/human-phenotype-ontology/blob/v2026-09-01/docs/annotations/phenotype_to_genes.md

## 2. Download Instructions

### 2.1 Download Using Linux Command Line

Use `wget` to download the file from the specified HPO release:

``` bash
wget -c \
  -O phenotype_to_genes.txt \
  https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2026-09-01/phenotype_to_genes.txt
```

Alternatively, use `curl`:

``` bash
curl -L \
  -o phenotype_to_genes.txt \
  https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2026-09-01/phenotype_to_genes.txt
```

Options:

-   `-c`: Resume an interrupted download (`wget`).
-   `-L`: Follow HTTP redirects (`curl`).
-   `-O` or `-o`: Specify the output filename.

### 2.2 Download Using a Web Browser

1.  Open the official HPO release page:\
    https://github.com/obophenotype/human-phenotype-ontology/releases/tag/v2026-09-01
2.  Expand the **Assets** section.
3.  Locate `phenotype_to_genes.txt`.
4.  Click the filename to download the file.

**Note:** Download `phenotype_to_genes.txt` directly from the release
assets rather than downloading the entire source code archive.

## 3. File Format

The file is a tab-delimited text file (TSV) containing five columns.

  Column           Description
  ---------------- ------------------------------------------------------
  `hpo_id`         HPO term identifier
  `hpo_name`       HPO term name
  `ncbi_gene_id`   NCBI Gene identifier
  `gene_symbol`    Gene symbol
  `disease_id`     Associated disease identifier, such as OMIM or ORPHA

Example:

``` text
hpo_id  hpo_name    ncbi_gene_id    gene_symbol disease_id
HP:0009830  Peripheral neuropathy   773 CACNA1A OMIM:183086
HP:0009830  Peripheral neuropathy   87178   PNPT1   ORPHA:101111
HP:0009830  Peripheral neuropathy   87178   PNPT1   ORPHA:319514
HP:0009830  Peripheral neuropathy   87178   PNPT1   OMIM:608703
```

The actual file uses tab characters to separate columns. The example
above illustrates the file structure.

## 4. Usage and Important Notes

1.  A single HPO term can be associated with multiple genes.
2.  A gene may be associated with multiple diseases. Consequently, the
    same HPO term and gene symbol may appear in multiple rows.
3.  The file includes associations involving ancestor phenotype terms,
    not only the most specific phenotypes directly annotated to a
    disease.
4.  According to the official documentation, the file includes phenotype
    associations restricted to descendants of `Phenotypic abnormality`
    (`HP:0000118`).
5.  Keep the HPO data version consistent with the version used in
    downstream analyses to ensure reproducibility and traceability.

## 5. Verify the Downloaded File

Check that the file exists and inspect its size:

``` bash
ls -lh phenotype_to_genes.txt
```

Display the first five lines:

``` bash
head -n 5 phenotype_to_genes.txt
```

Check the number of tab-delimited fields in the first five lines:

``` bash
awk -F '\t' 'NR<=5 {print NF}' phenotype_to_genes.txt
```

Each data row should contain five columns. If the file includes a
header, the header should also contain five columns.

## 6. Data Version Management

To maintain traceability across releases, store each HPO version
separately rather than overwriting previously downloaded data.

Recommended directory structure:

``` text
human-phenotype-ontology/
├── v2026-09-01/
│   └── phenotype_to_genes.txt
└── README.md
```

When updating the dataset, check the official HPO release page for the
latest version and replace `v2026-09-01` in the download URL with the
desired release tag.

## 7. References

-   HPO official repository:
    https://github.com/obophenotype/human-phenotype-ontology
-   HPO v2026-09-01 release:
    https://github.com/obophenotype/human-phenotype-ontology/releases/tag/v2026-09-01
-   Official documentation for `phenotype_to_genes.txt`:
    https://github.com/obophenotype/human-phenotype-ontology/blob/v2026-09-01/docs/annotations/phenotype_to_genes.md
