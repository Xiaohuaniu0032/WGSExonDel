# MANE Reference Files

- `MANE.GRCh38.v1.4.Summary.txt`: source MANE summary downloaded from the NCBI MANE FTP directory.
- `MANERefSeqSelect.tsv`: derived table containing gene symbol, RefSeq transcript accession, and MANE status.

The derived table is generated with `scripts/BuildReference/ParseMANESummaryFile.pl`. The current parser removes transcript version suffixes to preserve the format used by the existing coordinate-generation workflow.

Source: https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/
