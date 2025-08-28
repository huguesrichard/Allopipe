#coding:utf-8
import os
import pandas as pd

import tools.aams_helpers as aams_helpers

def pickle_parsing(str_params, args):
    #############################################################################################
    # TO FIX: regarding how we pick up D/R
    donor_or_recipient = "HG002-VEPannotated"
    #############################################################################################
    
    import os
    # removes base_length from str_params
    str_params_split = "_".join(str_params.split("_")[i] for i in [0, 1, 2, 4])
    pickle_path = os.path.join(
        f"../output/runs/",
        f"{args.run_name}",
        f"run_tables",
        f"{donor_or_recipient}_vep_infos_table_{str_params_split}.pkl"
    )

    pickle_df = pd.read_pickle(pickle_path)

    # remove rows where 'INFO' is NaN or empty
    pickle_df = pickle_df[pickle_df['INFO'].notna() & (pickle_df['INFO'] != '')]
    # explode the 'INFO' column to handle multiple entries
    pickle_df = pickle_df.explode("INFO", ignore_index=True)

    # Define VEP field indices
    # see `vcf_vep_parser()` for proper parsing
    ENSG_INDEX = 4   # gene
    ENST_INDEX = 6   # transcript
    ENST_PROT_POS = 14 # protein position

    # Function to extract ENSG and ENST from INFO string
    def extract_ensg_enst(info_str):
        if pd.isna(info_str):
            return pd.Series([None, None, None])
        fields = info_str.split("|")
        ensg = fields[ENSG_INDEX] if len(fields) > ENSG_INDEX else None
        enst = fields[ENST_INDEX] if len(fields) > ENST_INDEX else None
        prot_pos = fields[ENST_PROT_POS] if len(fields) > ENST_PROT_POS else None
        return pd.Series([ensg, enst, prot_pos])

    # Apply to df_mini and drop INFO
    pickle_df[["Gene_id", "Transcript_id", "Protein_position"]] = pickle_df["INFO"].apply(extract_ensg_enst)
    pickle_df = pickle_df.drop(columns=["INFO"])

    return pickle_df


def mm_intersect(mismatches_df, transcripts_pair):
    # Merge on ENSG
    mismatches_df["CHROM"] = mismatches_df["CHROM"].astype(str)
    mismatches_df["Transcript_id"] = mismatches_df["transcripts_x"]
    mismatches_df["Gene_id"] = mismatches_df["genes_x"]
    merged = pd.merge(
        mismatches_df,
        transcripts_pair,
        on=["CHROM", "Transcript_id", "Gene_id"],
        how="outer"
    )
    # Keep only columns present in transcripts_pair
    merged = merged.reindex(columns=transcripts_pair.columns)
    return merged


def add_pep_seq_chop(transcripts_pair, peptides_ensembl):
    transcripts_pair["CHROM"] = transcripts_pair["CHROM"].astype(str)
    peptides_ensembl["CHROM"] = peptides_ensembl["CHROM"].astype(str)
    transcripts_pair = pd.merge(
        transcripts_pair,
        peptides_ensembl,
        how="inner",
        on=["Gene_id", "Transcript_id", "CHROM"],
    )
    return transcripts_pair


def netchop_table_prep(mismatches_df, transcripts_pair, peptides_ensembl):
    chop_table = mm_intersect(mismatches_df, transcripts_pair)
    # Get the rows with missing peptide_ALT
    missing_mask = chop_table["peptide_ALT"].isna()
    missing_rows = chop_table[missing_mask].copy()
    # Preserve the original index as a column
    missing_rows["orig_index"] = missing_rows.index
    # Merge using your function
    filled_rows = add_pep_seq_chop(missing_rows, peptides_ensembl)

    # Now use 'orig_index' to map back to chop_table
    filled_rows = filled_rows.set_index("orig_index")

    # In-place update of desired columns
    chop_table.loc[filled_rows.index, "peptide_ALT"] = filled_rows["Sequence_aa_y"].values
    chop_table.loc[filled_rows.index, "Peptide_id"] = filled_rows["Peptide_id_y"].values
    
    # remove duplicates and empty rows
    chop_table = chop_table.drop_duplicates(
        subset=["CHROM", "Gene_id", "Transcript_id", "Peptide_id", "peptide_ALT"]
    )
    chop_table = chop_table.dropna(subset=["peptide_ALT"])
    chop_table.to_csv("chop_table.tsv", sep="\t", index=False)

    return chop_table


def run_netchop(chop_table, args):
    # prepare input for netchop
    chop_fasta = "chop_peptides.fasta"
    with open(chop_fasta, "w") as fasta_file:
        for _, row in chop_table.iterrows():
            # NetChop accepts only 10 characters for ID: keep the 10 last characters to have full ENSP ID
            fasta_file.write(f">{row['Peptide_id'][-10:]}\n{row['peptide_ALT']}\n")


    print(f"[{args.pair}] Entering NetChop handler : running NetChop may last a few minutes...")
    
    #######################################################################################################
    # netmhc_run_output = os.path.join(
    # netmhc_dir,
    # (args.pair + "_" if args.pair else "") +
    # args.run_name + "_full_run_information.txt")
    chop_output = "netchop_output.txt"
    #######################################################################################################

    os.system(f"netchop {chop_fasta} -verbose -tdir /tmp > {chop_output}")

    return chop_output


def parse_netchop_output(filepath, min_run_length=1):
    result = {}  # Dict: ident -> list of run dicts
    current_run = []
    prev_id = None
    all_positions = {}  # ident -> list of (pos, aa)

    def store_run(ident, run):
        if len(run) >= min_run_length:
            result.setdefault(ident, []).append({
                'start': run[0],
                'end': run[-1],
                'length': len(run),
                'positions': run
            })

    with open(filepath) as f:
        dash_count = 0
        in_data_block = False

        for line in f:
            line = line.strip()

            if not in_data_block:
                if set(line) <= set("-"):
                    dash_count += 1
                    if dash_count == 2:
                        in_data_block = True
                    continue
                else:
                    continue

            if line.lower().startswith("pos"):
                continue

            parts = line.split()
            if len(parts) < 5:
                continue

            pos, aa, cs, score, ident = parts[:5]
            try:
                pos = int(pos)
            except ValueError:
                continue

            all_positions.setdefault(ident, []).append((pos, aa))

            if cs == ".":
                if ident != prev_id:
                    if prev_id and current_run:
                        store_run(prev_id, current_run)
                    current_run = [pos]
                else:
                    if current_run and pos == current_run[-1] + 1:
                        current_run.append(pos)
                    else:
                        store_run(ident, current_run)
                        current_run = [pos]
            else:
                store_run(prev_id, current_run)
                current_run = []

            prev_id = ident

    store_run(prev_id, current_run)

    return result, all_positions


# Load mapping: short ID -> full Ensembl Peptide ID
def load_peptide_id_map(tsv_path):
    mapping = {}
    with open(tsv_path) as f:
        header = f.readline().strip().split('\t')
        try:
            idx = header.index("Peptide_id")
        except ValueError:
            raise ValueError("Column 'Peptide_id' not found in chop_table.tsv")

        for line in f:
            cols = line.strip().split('\t')
            if len(cols) <= idx:
                continue
            full_id = cols[idx]
            short_id = full_id[-10:]  # take last 10 characters (e.g., "0000356701")
            mapping[short_id] = full_id
    return mapping


def postprocess_netchop(chop_output):
    dot_runs, all_positions = parse_netchop_output(chop_output, min_run_length=9)
    id_map = load_peptide_id_map("chop_table.tsv")

    # Write sequences to CSV
    with open("chop_pep.csv", "w") as out:
        out.write("ID,Peptide\n")
        for ident in dot_runs:
            aa_lookup = dict(all_positions[ident])
            short_id = ident[-10:]
            full_id = id_map.get(short_id, short_id)  # fallback to short_id if not found

            for run in dot_runs[ident]:
                seq = "".join(aa_lookup.get(pos, "X") for pos in run["positions"])
                out.write(f"{full_id},{seq}\n")