import sys
import pandas as pd

# Get the command-line arguments
mismatches_table = sys.argv[1]
filtered_mismatches_table = sys.argv[2]
rs_txt = sys.argv[3]

# Get rs list from .txt
with open(rs_txt, 'r', encoding='utf-8') as f:
    rs_list = set(line.strip() for line in f if line.strip())

# load the mismatches table
df = pd.read_csv(mismatches_table, sep='\t')

# Keep lines if the rsID is in the rs list
rs_MM_table = df[
    df['ID_x'].apply(lambda x: str(x) in rs_list) |
    df['ID_y'].apply(lambda x: str(x) in rs_list)
]

# Save into TSV
rs_MM_table.to_csv(filtered_mismatches_table, sep='\t', index=False)

