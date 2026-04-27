#coding:utf-8
import os
import pandas as pd
import csv
from tools import aams_helpers, parsing_functions


def validate_cleavage_imputation(args):
    """
    Cleavage mode is incompatible with no-imputation mode.
    Raises ValueError when this invalid combination is detected.
    """
    imputation_mode = aams_helpers.read_log_field(args, "Imputation")
    if imputation_mode == "no-imputation":
        raise ValueError(
            f"Cleavage mode cannot be run with 'no-imputation' mode. "
            "Please rerun Allo-Count with 'imputation' mode or disable `--cleavage`."
        )


def _get_vep_indices_from_vcf(vcf_path):
    if vcf_path.endswith(".vcf"):
        _, vep_indices = parsing_functions.vcf_vep_parser(vcf_path)
    elif vcf_path.endswith(".vcf.gz"):
        _, vep_indices = parsing_functions.gzvcf_vep_parser(vcf_path)
    return vep_indices


def pickle_parsing(str_params, args, individual):
    if individual == "donor":
        vcf_path_indiv = aams_helpers.read_log_field(args, "Donor")
    elif individual == "recipient":
        vcf_path_indiv = aams_helpers.read_log_field(args, "Recipient")
    else:
        raise ValueError("individual must be 'donor' or 'recipient'")

    sample = os.path.basename(vcf_path_indiv).split(".")[0]

    # removes base_length from str_params
    str_params_split = "_".join(str_params.split("_")[i] for i in [0, 1, 2, 4])
    pickle_path = os.path.join(
        args.output_dir,
        "runs",
        args.run_name,
        "run_tables",
        f"{sample}_vep_infos_table_{str_params_split}.pkl",
    )

    pickle_df = pd.read_pickle(pickle_path)

    # remove rows where 'INFO' is NaN or empty
    pickle_df = pickle_df[pickle_df["INFO"].notna() & (pickle_df["INFO"] != "")]
    # explode the 'INFO' column to handle multiple entries
    pickle_df = pickle_df.explode("INFO", ignore_index=True)

    # Derive field indices from VCF header to avoid hard-coded positions
    vep_indices = _get_vep_indices_from_vcf(vcf_path_indiv)

    # Function to extract ENSG, ENST and protein position from INFO string
    def extract_ensg_enst(info_str):
        if pd.isna(info_str):
            return pd.Series([None, None, None])
        fields = info_str.split("|")
        ensg = fields[vep_indices.gene] if len(fields) > vep_indices.gene else None
        enst = fields[vep_indices.transcript] if len(fields) > vep_indices.transcript else None
        prot_pos = fields[vep_indices.prot] if len(fields) > vep_indices.prot else None
        return pd.Series([ensg, enst, prot_pos])

    # Apply to pickle_df and drop INFO
    pickle_df[["Gene_id", "Transcript_id", "Protein_position"]] = pickle_df["INFO"].apply(extract_ensg_enst)
    pickle_df = pickle_df.drop(columns=["INFO"])

    return pickle_df


def mm_intersect(mismatches_df, transcripts_pair):
    # Merge on ENSG
    mismatches_df["CHROM"] = mismatches_df["CHROM"].astype(str)
    mismatches_df["POS"] = mismatches_df["POS"].astype(str)
    mismatches_df["Transcript_id"] = mismatches_df["transcripts_x"]
    mismatches_df["Gene_id"] = mismatches_df["genes_x"]
    merged = pd.merge(
        mismatches_df,
        transcripts_pair,
        on=["CHROM", "POS", "Transcript_id", "Gene_id"],
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


def netchop_table_prep(mismatches_df, transcripts_pair, peptides_ensembl, args, netchop_dir, sample_suffix=""):
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
    
    # remove duplicates and empty rows and unused columns
    chop_table = chop_table.drop_duplicates(
        subset=["Gene_id", "Transcript_id", "Peptide_id", "peptide_ALT"]
    )
    chop_table = chop_table.dropna(subset=["peptide_ALT"])
    chop_table = chop_table.drop(columns=["aa_REF", "aa_ALT"])
    
    # save table
    chop_table_path = os.path.join(
        netchop_dir,
        f"{args.pair + '_' if args.pair else ''}{args.run_name}{sample_suffix}_netchop_table.csv"
    )
    chop_table.to_csv(chop_table_path, sep="\t", index=False)

    return chop_table, chop_table_path


def run_netchop(chop_table, args, netchop_dir, sample_suffix=""):
    # prepare input for netchop
    chop_fasta = os.path.join(
        netchop_dir,
        f"{args.pair + '_' if args.pair else ''}{args.run_name}{sample_suffix}_netchop_peptides.fa"
    )
    with open(chop_fasta, "w") as fasta_file:
        for _, row in chop_table.iterrows():
            # NetChop accepts only 10 characters for ID: keep the 10 last characters to have full ENSP ID
            fasta_file.write(f">{row['Peptide_id'][-10:]}\n{row['peptide_ALT']}\n")
    
    chop_output = os.path.join(
        netchop_dir,
        f"{args.pair + '_' if args.pair else ''}{args.run_name}{sample_suffix}_netchop_output.txt"
    )

    # netChop command
    os.system(f"netchop {chop_fasta} -verbose > {chop_output}")

    return chop_output


def parse_netchop_output(filepath, min_run_length):
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
    table = pd.read_csv(tsv_path, sep="\t")
    if "Peptide_id" not in table.columns:
        raise ValueError("Column 'Peptide_id' not found in the table")

    for _, row in table.iterrows():
        peptide_id = row["Peptide_id"]
        short_id = str(peptide_id)[-10:]  # take last 10 characters (e.g., "0000356701")
        chrom = str(row["CHROM"]) if "CHROM" in table.columns and pd.notna(row["CHROM"]) else ""
        pos = str(row["POS"]).strip() if "POS" in table.columns and pd.notna(row["POS"]) else ""
        mapping.setdefault(short_id, []).append((peptide_id, chrom, pos))
    return mapping


def postprocess_netchop(chop_output, chop_table_path, args, netchop_dir, sample_suffix=""):
    dot_runs, all_positions = parse_netchop_output(chop_output, args.length)
    id_map = load_peptide_id_map(chop_table_path)

    # save peptides table
    chop_peptides = os.path.join(
        netchop_dir,
        f"{args.pair + '_' if args.pair else ''}{args.run_name}{sample_suffix}_netchop_peptides.csv"
    )
    aggregated = {}
    with open(chop_peptides, "w", newline="") as out:
        for ident in dot_runs:
            aa_lookup = dict(all_positions[ident])
            short_id = ident[-10:]
            metadata_rows = id_map.get(short_id, [(short_id, "", "")])  # fallback to short_id if not found

            for run in dot_runs[ident]:
                seq = "".join(aa_lookup.get(pos, "X") for pos in run["positions"])
                hla_peptides = aams_helpers.peptide_seg(seq, int(args.length))
                for peptide_id, chrom, pos in metadata_rows:
                    for hla_peptide in hla_peptides:
                        key = (chrom, pos, peptide_id, seq, hla_peptide)
                        if key not in aggregated:
                            aggregated[key] = True

        writer = csv.writer(out)
        writer.writerow(["CHROM", "POS", "Peptide_id", "Peptide", "hla_peptides"])
        for chrom, pos, peptide_id, long_peptide, hla_peptide in aggregated.keys():
            writer.writerow([chrom, pos, peptide_id, long_peptide, hla_peptide])

    return chop_peptides


def load_peptides(path):
    peptides = []
    seen = set()
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            peptide_id = str(row.get("Peptide_id", row.get("ID", ""))).strip()
            peptide = str(row.get("Peptide", row.get("hla_peptides", ""))).strip()
            hla_peptide = str(row.get("hla_peptides", row.get("Peptide", ""))).strip()
            chrom = str(row.get("CHROM", "")).strip()
            pos = str(row.get("POS", "")).strip()
            if peptide_id and hla_peptide and peptide:
                item = (peptide_id, peptide, hla_peptide, chrom, pos)
                if item not in seen:
                    seen.add(item)
                    peptides.append(item)
    return peptides


def deduce_cleaved_peptides(donor_peptides_path, recipient_peptides_path, netchop_dir, args, pair_print):
    donor_peptides = load_peptides(donor_peptides_path)
    recipient_peptides = load_peptides(recipient_peptides_path)
    donor_pairs = {(peptide_id, hla_peptide) for peptide_id, _, hla_peptide, _, _ in donor_peptides}
    recipient_pairs = {(peptide_id, hla_peptide) for peptide_id, _, hla_peptide, _, _ in recipient_peptides}

    orientation = aams_helpers.read_log_field(args, "Orientation")
    if orientation == "dr":
        total_source = len(donor_pairs)
        kept_pairs = donor_pairs - recipient_pairs
        source_rows = donor_peptides
    elif orientation == "rd":
        total_source = len(recipient_pairs)
        kept_pairs = recipient_pairs - donor_pairs
        source_rows = recipient_peptides
    else:
        raise ValueError(f"Invalid orientation: {orientation}")

    removed_count = total_source - len(kept_pairs)
    removed_percentage = (removed_count / total_source * 100) if total_source else 0.0
    print(
        f"{pair_print}Cleavage deduced k-mers: kept {len(kept_pairs)}/{total_source}, "
        f"removed {removed_count} ({removed_percentage:.2f}%)."
    )

    deduced_pep_path = os.path.join(
        netchop_dir,
        f"{args.pair + '_' if args.pair else ''}{args.run_name}_netchop_peptides.csv"
    )
    aggregated = {}
    with open(deduced_pep_path, "w", newline="") as out:
        for peptide_id, peptide, hla_peptide, chrom, pos in source_rows:
            if (peptide_id, hla_peptide) in kept_pairs:
                key = (chrom, pos, peptide_id, peptide, hla_peptide)
                if key not in aggregated:
                    aggregated[key] = True

        writer = csv.writer(out)
        writer.writerow(["CHROM", "POS", "Peptide_id", "Peptide", "hla_peptides"])
        for chrom, pos, peptide_id, long_peptide, hla_peptide in aggregated.keys():
            writer.writerow([chrom, pos, peptide_id, long_peptide, hla_peptide])

    return deduced_pep_path


def prepare_cleavage_netmhcpan_inputs(aams_run_tables, args, pep_indiv_path, deduced_pep_path):
    deduced_df = pd.read_csv(deduced_pep_path)
    peptide_length = int(args.length)
    deduced_df = deduced_df.dropna(subset=["Peptide_id", "hla_peptides"]).copy()
    deduced_df["Peptide_id"] = deduced_df["Peptide_id"].astype(str)
    deduced_df["hla_peptides"] = deduced_df["hla_peptides"].astype(str)
    if "Peptide" in deduced_df.columns:
        deduced_df["Peptide"] = deduced_df["Peptide"].astype(str)
    else:
        deduced_df["Peptide"] = deduced_df["hla_peptides"]
    deduced_df = deduced_df[deduced_df["hla_peptides"].str.len() == peptide_length].drop_duplicates()
    deduced_df = deduced_df.reset_index(drop=True)
    deduced_df["input_order"] = deduced_df.index

    # Expand deduced hla peptides rows, keeping CHROM/POS from cleavage outputs
    deduced_df["CHROM"] = deduced_df["CHROM"].astype(str) if "CHROM" in deduced_df.columns else ""
    deduced_df["POS"] = deduced_df["POS"].astype(str) if "POS" in deduced_df.columns else ""

    expanded_rows = []
    for _, row in deduced_df.iterrows():
        chrom = row.get("CHROM", "")
        pos = row.get("POS", "")
        chrom = None if pd.isna(chrom) or str(chrom).strip() == "" else str(chrom).strip()
        pos = None if pd.isna(pos) or str(pos).strip() == "" else str(pos).strip()

        expanded_rows.append(
            {
                "CHROM": chrom,
                "POS": pos,
                "Peptide_id": row["Peptide_id"],
                "peptide": row["Peptide"],
                "hla_peptides": row["hla_peptides"],
                "input_order": row["input_order"],
            }
        )

    if not expanded_rows:
        raise ValueError(
            f"No hla peptides of length {peptide_length} found in deduced cleavage peptides."
        )
    deduced_expanded = pd.DataFrame(expanded_rows)

    # Recover non-cleavage metadata (e.g., Gene_id/Transcript_id) by Peptide_id
    pep_df_original = pd.read_pickle(pep_indiv_path).copy()
    metadata_cols = [
        col for col in pep_df_original.columns
        if col not in ["CHROM", "POS", "peptide", "peptide_REF", "hla_peptides", "hla_peptides_REF"]
    ]
    meta_map = pep_df_original[metadata_cols].drop_duplicates("Peptide_id")
    deduced_expanded = pd.merge(
        deduced_expanded,
        meta_map,
        how="left",
        on=["Peptide_id"],
    )

    deduced_expanded = deduced_expanded.dropna(subset=["Gene_id"]).copy()
    deduced_expanded = deduced_expanded.sort_values("input_order", kind="stable")
    if deduced_expanded.empty:
        raise ValueError(
            "No deduced peptides could be mapped to Gene_id using the original peptide table."
        )

    # Rebuild peptide-level rows with the same column order as pep_df_original
    deduced_expanded["peptide_REF"] = None
    deduced_expanded["hla_peptides_REF_item"] = ""
    group_cols = [col for col in pep_df_original.columns if col not in ["hla_peptides", "hla_peptides_REF"]]
    for col in group_cols:
        if col not in deduced_expanded.columns:
            deduced_expanded[col] = None
    pep_df = (
        deduced_expanded.groupby(group_cols, dropna=False, as_index=False, sort=False)
        .agg({"hla_peptides": list, "hla_peptides_REF_item": list})
        .rename(columns={"hla_peptides_REF_item": "hla_peptides_REF"})
    )
    pep_df = pep_df[pep_df["hla_peptides"].map(len) > 0].copy()
    pep_df = pep_df.reindex(columns=pep_df_original.columns)

    filtered_pep_path = os.path.join(
        aams_run_tables,
        f"{args.pair + '_' if args.pair else ''}{args.run_name}_pep_df_cleavage_filtered.pkl"
    )
    filtered_pep_tsv_path = os.path.join(
        aams_run_tables,
        f"{args.pair + '_' if args.pair else ''}{args.run_name}_pep_df_cleavage_filtered.tsv"
    )

    filtered_fasta_path = os.path.join(
        aams_run_tables,
        f"{args.pair + '_' if args.pair else ''}{args.run_name}_cleavage_deduced_fasta.fa"
    )
    # Build NetMHCpan FASTA from deduced k-mers with per-gene peptide indexing
    fasta_df = pep_df[["Gene_id", "hla_peptides"]].explode("hla_peptides").dropna(subset=["hla_peptides"])
    fasta_df = fasta_df.drop_duplicates(subset=["Gene_id", "hla_peptides"]).reset_index(drop=True)
    with open(filtered_fasta_path, "w", encoding="utf-8") as handle:
        gene_counts = {}
        for _, row in fasta_df.iterrows():
            gene_id = row["Gene_id"]
            gene_counts[gene_id] = gene_counts.get(gene_id, 0) + 1
            handle.write(f">{gene_id}:{gene_counts[gene_id]}\n")
            handle.write(f"{row['hla_peptides']}\n")

    pep_df.to_pickle(filtered_pep_path)
    pep_df.to_csv(filtered_pep_tsv_path, sep="\t", index=False)
    
    return filtered_fasta_path, filtered_pep_path