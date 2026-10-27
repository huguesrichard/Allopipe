import sys
import pandas as pd

# Get the command-line arguments
mismatches_table = sys.argv[1]
filtered_mismatches_table = sys.argv[2]
bed_file = sys.argv[3]

# load the BED file into a DataFrame
bed_df = pd.read_csv(
    bed_file,
    sep="\t",
    header=None,
    usecols=[0, 1, 2],
    names=["CHROM", "START", "END"],
    dtype={"CHROM": str},
)
# remove the 'chr' prefix from the 'CHROM' column for comparison
bed_chrom = bed_df["CHROM"].str.replace(r"^chr", "", regex=True)

# load the mismatches table
mismatches_df = pd.read_csv(mismatches_table, sep="\t", dtype={"CHROM": str})
mismatch_chrom = mismatches_df["CHROM"].str.replace(r"^chr", "", regex=True)

# convert columns to numeric types
starts = pd.to_numeric(bed_df["START"], errors="coerce")
ends = pd.to_numeric(bed_df["END"], errors="coerce")
positions = pd.to_numeric(mismatches_df["POS"], errors="coerce")

# adjust 1-based mismatch positions for comparison with zero-based BED intervals
# filter the positions to keep only those within an interval
keep = pd.Series(False, index=mismatches_df.index)
for chrom, start, end in zip(bed_chrom, starts, ends):
    keep |= (
        (mismatch_chrom == chrom)
        & (positions - 1 >= start)
        & (positions - 1 < end)
    )

# save the filtered result to a new TSV file
mismatches_df.loc[keep].to_csv(filtered_mismatches_table, sep="\t", index=False)
