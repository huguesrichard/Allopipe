import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import os

def hist(ams_exp_path, run_plots):

    # Load the TSV file into a DataFrame
    data = pd.read_csv('../data/distrib_AMS.tsv', sep='\t')

    # Filter the data for "related"
    data_related = data[data['Relatedness'] == 'related']

    # Filter the data for "unrelated"
    data_non_related = data[data['Relatedness'] == 'unrelated']

    # Extract the column of interest (AMS) for each group
    AMS_related = data_related['AMS']  # Extract the 'AMS' column for the "related" group
    AMS_unrelated = data_non_related['AMS']  # Extract the 'AMS' column for the "unrelated" group

    # Plot the distributions for the two groups
    plt.figure(figsize=(10, 6))

    # Plot the distribution for the "related" group
    sns.histplot(AMS_related, kde=True, color='blue', label='Related pairs', stat='density', bins=30)

    # Plot the distribution for the "unrelated" group
    sns.histplot(AMS_unrelated, kde=True, color='red', label='Unrelated pairs', stat='density', bins=30)

    # Define a directory for searching files
    directory = ams_exp_path

    # List all files in the directory, excluding subdirectories
    files = [f for f in os.listdir(directory) if not os.path.isdir(os.path.join(directory, f))]

    # Find a .tsv or .csv file
    tsv_file = next((f for f in files if f.endswith('.tsv')), None)
    csv_file = next((f for f in files if f.endswith('.csv')), None)

    # Check for the priority of the .tsv file if available, otherwise use .csv
    if tsv_file:
        # Load the .tsv file
        df_ams = pd.read_csv(os.path.join(directory, tsv_file), sep='\t')

        # Check that the 'ams_giab' column exists in the .tsv file
        if 'ams_giab' in df_ams.columns:
            special_values = df_ams['ams_giab'].values  # Extract the values from the 'ams_giab' column

            # Add a vertical line for each special value and position the label
            for ams in special_values:
                # Draw the vertical line
                plt.axvline(ams, color='red', linestyle='--')
        else:
            print("No 'ams_giab' column in .tsv file")

    elif csv_file:
        # Load the .csv file
        df_ams = pd.read_csv(os.path.join(directory, csv_file))

        # Check that the 'ams' column exists in the .csv file
        if 'ams' in df_ams.columns:
            special_values = df_ams['ams'].values  # Extract the values from the 'ams' column

            # Add a vertical line for each special value and position the label
            for ams in special_values:
                # Draw the vertical line
                plt.axvline(ams, color='red', linestyle='--', label=f'AMS of the pair : {ams}')
        else:
            print("No 'ams' column in the .csv file")

    else:
        print("No .tsv or .csv file in the directory")

    # Add labels and a legend
    plt.title('AMS Distribution Related vs Unrelated + computed AMS')
    plt.xlabel('AMS values')
    plt.ylabel('Density')
    plt.legend()

    # Place the legend at the top of the plot, aligned to the left
    plt.legend(loc='upper left')

    # Save the plot as a PNG file
    plot_path = run_plots + "/" + "distrib.png"
    plt.savefig(plot_path, format='png')
    plt.close()