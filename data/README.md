# Data

- `Ensembl/`: reference sequences and transcript mappings, organized by genome assembly (see the [Ensembl README](Ensembl/README.md)).
- `transcripts_lists/`: tissue-specific gene lists from the Human Protein Atlas, used to filter mismatches tables with [`args_gene-transcript_filtering.py`](../scripts/args_gene-transcript_filtering.py) (see [scripts](../scripts/)).
- `consequences.txt`: variant consequences used by the pipeline.
- `distrib_AMS.tsv`: reference AMS distributions for related and unrelated pairs, used to generate `distrib.png` by the pipeline in `<output_dir>/runs/<run_name>/plots/`.
