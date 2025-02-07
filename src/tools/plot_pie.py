import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np  # For color palette
import os

def pie(run_tables, run_plots, pair, run_name):

    # Function to safely read a TSV file
    def safe_read_tsv(file_path):
        try:
            return pd.read_csv(file_path, delimiter='\t', on_bad_lines='skip')  # Skip problematic lines
        except pd.errors.ParserError as e:
            print(f"Error reading {file_path}: {e}")
            return None

    # Directory containing the files
    directory = run_tables

    # List of TSV files in the directory containing the word "mismatch"
    files = [f for f in os.listdir(directory) if f.endswith('.tsv') and 'mismatches' in f]

    # Custom sort function for chromosomes (1-22, then X and Y)
    def custom_sort_key(chrom):
        try:
            return int(chrom)  # Convert numeric chromosomes to int for sorting
        except ValueError:
            return float('inf')  # Place 'X' and 'Y' at the end

    # Iterate through each file and analyze the 'CHROM' column separately
    for file in files:
        file_path = os.path.join(directory, file)
        df = safe_read_tsv(file_path)

        if df is not None and 'CHROM' in df.columns:
            # Ensure the 'CHROM' column is treated as strings
            df['CHROM'] = df['CHROM'].astype(str).str.strip()

            # Count occurrences of each value in the 'CHROM' column for the current file
            counts = df['CHROM'].value_counts()

            # Sort the counts using the custom sort key
            sorted_counts = counts.sort_index(key=lambda x: x.map(custom_sort_key))

            # Convert to DataFrame for easier manipulation
            counts_df = sorted_counts.reset_index()
            counts_df.columns = ['Chromosome', 'Count']

            # Calculate percentages
            total_count = counts_df['Count'].sum()
            counts_df['Percentage'] = (counts_df['Count'] / total_count) * 100
            counts_df['Percentage'] = counts_df['Percentage'].round(1)  # Round percentages to 1 decimal place

            # Display the table with chromosome details
            # print(f"\nTable for {file}:")
            # print(counts_df)

            # Create a pie chart for the current file
            plt.figure(figsize=(14, 10))
            wedges, texts = plt.pie(
                counts_df['Count'],
                labels=counts_df['Chromosome'],
                autopct=None,  # No percentages shown in the chart
                startangle=90,  # Start at the top (90°)
                colors=plt.cm.Blues(np.linspace(0, 1, len(counts_df))),
                counterclock=False  # Ensure the chart goes clockwise
            )
            plt.title(f'Distribution of AMS across chromosoms {pair}')

            # Position the table to the right of the chart
            ax = plt.gca()
            table = ax.table(
                cellText=counts_df.values,
                colLabels=['Chromosome', 'Count', 'Percentage'],
                colWidths=[0.1, 0.1, 0.1],
                cellLoc='center',
                loc='right',
                bbox=[1.05, 0.1, 0.4, 0.8]  # Adjust the table position to the right
            )
            table.auto_set_font_size(False)
            table.set_fontsize(10)
            table.scale(1.2, 1.2)  # Adjust the table size for better readability

            # Save the pie chart with the table as a PNG image
            output_path = f"{run_plots}/{run_name}_{pair + '_' if pair else ''}pie_chart.png"
            plt.savefig(output_path, format='png', dpi=300, bbox_inches='tight')
            plt.close()
        else:
            print(f"File {file} is missing the 'CHROM' column or could not be read.")


