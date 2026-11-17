# Ensembl reference files

To reconstruct peptides, Allo-Affinity reads the following reference files for the selected genome assembly:

Only the GRCh37 RefSeq TSV is included in the repository. At the start of a Nextflow run, missing files are downloaded into `data/Ensembl/GRCh37/` or `data/Ensembl/GRCh38/`. Downloads are checked against the Ensembl `CHECKSUMS` files, and their source, checksum and download time (UTC) are recorded in `.ensembl_downloads.json`. Files downloaded by an earlier run are checked against the recorded checksum before reuse.

For GRCh38, `--ensembl_release` selects the release to download. By default, it is derived from `--vep_version` (for example, `release_116.2` gives `116`). Specify `--ensembl_release` explicitly if the VEP tag has another format. Before changing the release, remove the previously downloaded GRCh38 references so that files from different releases are not mixed. GRCh37 FASTA files come from the Ensembl GRCh37 archive; the provided GRCh37 RefSeq TSV is retained unchanged.

The pipeline warns if an identifiable VEP release differs from the GRCh38 reference release. With `--skip_vep_annotation`, it reads the VEP version from the input VCF headers; otherwise, it uses `--vep_version`.

## GRCh38 reference files

- `Homo_sapiens.<REFERENCE-GENOME>.cdna.all.fa.gz`: [download GRCh38 cDNA FASTA](https://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz) ([directory](https://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/cdna/)).
- `Homo_sapiens.<REFERENCE-GENOME>.pep.all.fa.gz`: [download GRCh38 peptide FASTA](https://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz) ([directory](https://ftp.ensembl.org/pub/release-116/fasta/homo_sapiens/pep/)).
- `Homo_sapiens.<REFERENCE-GENOME>.<VEP-VERSION>.refseq.tsv.gz`: [download GRCh38 RefSeq TSV](https://ftp.ensembl.org/pub/release-116/tsv/homo_sapiens/Homo_sapiens.GRCh38.116.refseq.tsv.gz) ([directory](https://ftp.ensembl.org/pub/release-116/tsv/homo_sapiens/)).

These links point to GRCh38 in Ensembl release 116, the latest release when they were documented. For another release, change the release number in the URL and the RefSeq filename. Select files for the genome assembly used for alignment; GRCh37 sources are described below.

> [!NOTE]
> To keep transcript annotations and reference sequences consistent, use the same Ensembl release as the VEP annotation.

## GRCh37 reference files

- `Homo_sapiens.GRCh37.cdna.all.fa.gz`: [download GRCh37 cDNA FASTA](https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.cdna.all.fa.gz) ([directory](https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/cdna/)).
- `Homo_sapiens.GRCh37.pep.all.fa.gz`: [download GRCh37 peptide FASTA](https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/pep/Homo_sapiens.GRCh37.pep.all.fa.gz) ([directory](https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/pep/)).
- `Homo_sapiens.GRCh37.refseq.tsv.gz`: see the generation process just below.

### Generate `Homo_sapiens.GRCh37.refseq.tsv.gz`

This file was exported from [BioMart](https://grch37.ensembl.org/biomart) in late 2026 using Ensembl Genes 116, Human genes (GRCh37.p13), default parameters and no filters. The downloaded file `mart_export.txt` was renamed.

To generate a new GRCh37 export in BioMart:

1. Under `- CHOOSE DATABASE -`, select `Ensembl Genes 116` (the latest release when this procedure was used).
2. Under `- CHOOSE DATASET -`, select `Human Genes (GRCh37.p13)`.
3. Keep the default attributes, including `Transcript stable ID`, and do not apply filters.
4. Click `Count`, then `Results`.
5. Under `Export all results to`, select `File` and `Compressed File (.gz)`.
6. Click `✅ Go`.

Rename the downloaded file for the pipeline and place it in `data/Ensembl/GRCh37/`:

```bash
mv mart_export.txt.gz data/Ensembl/GRCh37/Homo_sapiens.GRCh37.refseq.tsv.gz
```

GRCh37.p13 is the last published patch of GRCh37. The assembly and Ensembl's GRCh37 gene annotation are fixed, so a BioMart export of the same gene and transcript attributes is not expected to change. For instance, the current export contains the same data as the early 2024 export from Ensembl Genes 111.
