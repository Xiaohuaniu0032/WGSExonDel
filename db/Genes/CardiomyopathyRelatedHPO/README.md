# Cardiomyopathy-Related HPO Gene Lists

This directory contains HPO-associated gene lists for hypertrophic cardiomyopathy (HCM) and dilated cardiomyopathy (DCM).

## Files

| File | Description |
|---|---|
| `HPOTerms.tsv` | Defines the HPO terms to extract. It contains the HPO ID and term name for HCM (`0001639`, Hypertrophic cardiomyopathy) and DCM (`0001644`, Dilated cardiomyopathy). |
| `HypertrophicCardiomyopathy.genes` | Genes associated with the HCM HPO term in the HPO gene-phenotype annotation file. |
| `DilatedCardiomyopathy.genes` | Genes associated with the DCM HPO term in the HPO gene-phenotype annotation file. |

`HPOTerms.tsv` acts as the term list for extraction. The `ExtractHPOAssociatedGenes.pl` script reads the terms in this file and selects the corresponding gene-disease associations from the HPO `phenotype_to_genes.txt` file. In the current script implementation, matching is performed using the HPO term name (`hpo_name`), so the term names in `HPOTerms.tsv` must match those in the HPO annotation file.

## Data Source

The gene-phenotype associations are obtained from the Human Phenotype Ontology (HPO) project.

- **Source:** [Human Phenotype Ontology](https://github.com/obophenotype/human-phenotype-ontology)
- **Version:** `v2026-09-01`
- **Release:** [HPO v2026-09-01](https://github.com/obophenotype/human-phenotype-ontology/releases/tag/v2026-09-01)
- **Input file:** [`phenotype_to_genes.txt`](https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2026-09-01/phenotype_to_genes.txt)

Download the annotation file:

```bash
wget -c -O phenotype_to_genes.txt \
  https://github.com/obophenotype/human-phenotype-ontology/releases/download/v2026-09-01/phenotype_to_genes.txt
```

## Generate the Gene Lists

Run the following command from the repository root. The script uses `HPOTerms.tsv` to determine which HPO terms to extract and writes one gene list for each term.

```bash
mkdir -p tmp/CardiomyopathyRelatedHPO

perl scripts/HPO/ExtractHPOAssociatedGenes.pl \
  --gene-file phenotype_to_genes.txt \
  --term-file db/Genes/CardiomyopathyRelatedHPO/HPOTerms.tsv \
  --outdir tmp/CardiomyopathyRelatedHPO
```

The script generates:

- `Hypertrophic_cardiomyopathy.genes`
- `Dilated_cardiomyopathy.genes`

The filenames in this directory use the repository's standardized CamelCase naming convention. To place the generated files here with those names:

```bash
cp tmp/CardiomyopathyRelatedHPO/Hypertrophic_cardiomyopathy.genes \
  db/Genes/CardiomyopathyRelatedHPO/HypertrophicCardiomyopathy.genes

cp tmp/CardiomyopathyRelatedHPO/Dilated_cardiomyopathy.genes \
  db/Genes/CardiomyopathyRelatedHPO/DilatedCardiomyopathy.genes
```

## Output Format

The generated files are tab-delimited and contain three columns:

| Column | Description |
|---|---|
| `Gene` | Gene symbol |
| `Evidence` | HPO term name associated with the gene |
| `Note` | Associated disease identifier, such as an OMIM or ORPHA ID |

The script removes duplicate gene-disease associations within each output file. A gene may therefore appear in multiple rows when it is associated with multiple diseases.
