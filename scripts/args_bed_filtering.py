import sys
import pandas as pd

# Get the command-line arguments
mismatches_table = sys.argv[1]
filtered_mismatches_table = sys.argv[2]
bed_file = sys.argv[3]

# load the BED file into a DataFrame
bed_df = pd.read_csv(bed_file, sep='\t', header=None, names=['CHROM', 'START', 'END'])

# remove the 'chr' prefix from the 'CHROM' column
bed_df['CHROM'] = bed_df['CHROM'].str.replace('chr', '', regex=False)

# load the mismatches table
mismatches_df = pd.read_csv(mismatches_table, sep='\t')

# convert columns to numeric types
mismatches_df['POS'] = pd.to_numeric(mismatches_df['POS'], errors='coerce')
bed_df['START'] = pd.to_numeric(bed_df['START'], errors='coerce')
bed_df['END'] = pd.to_numeric(bed_df['END'], errors='coerce')

# drop rows with NaN values that may have been created during conversion
mismatches_df = mismatches_df.dropna(subset=['POS'])
bed_df = bed_df.dropna(subset=['START', 'END'])

# adjust the BED file's 'START' and 'END' columns for 1-based comparison (make them inclusive)
bed_df['START'] = bed_df['START'] + 1
bed_df['END'] = bed_df['END'] + 1


# filter the positions to keep only those within an interval
bed_MM_table = mismatches_df[
    mismatches_df.apply(
        lambda row: ((bed_df['CHROM'] == row['CHROM']) &
                     (bed_df['START'] <= row['POS']) &
                     (bed_df['END'] >= row['POS'])).any(), axis=1
    )
]

# save the filtered result to a new TSV file
bed_MM_table.to_csv(filtered_mismatches_table, sep='\t', index=False)

