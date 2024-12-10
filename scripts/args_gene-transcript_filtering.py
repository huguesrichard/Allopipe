import sys
import pandas as pd

# Get the command-line arguments
mismatches_table = sys.argv[1]
filtered_mismatches_table = sys.argv[2]
gene_transcript_txt = sys.argv[3]

# Get gene ID or transcript ID list from .txt
with open(gene_transcript_txt, 'r', encoding='utf-8') as f:
    gene_transcript_list = set(line.strip() for line in f if line.strip())

# load the mismatches table
df = pd.read_csv(mismatches_table, sep='\t')

# Keep lines if the gene ID or the transcript ID is in the list
gene_transcript_MM_table = df[
    df['transcripts_x'].apply(lambda x: str(x) in gene_transcript_list) |
    df['genes_x'].apply(lambda x: str(x) in gene_transcript_list) |
    df['transcripts_y'].apply(lambda x: str(x) in gene_transcript_list) |
    df['genes_y'].apply(lambda x: str(x) in gene_transcript_list)
]

# Save into TSV
gene_transcript_MM_table.to_csv(filtered_mismatches_table, sep='\t', index=False)
