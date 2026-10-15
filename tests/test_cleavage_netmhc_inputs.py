"""FASTA/table hand-off from cleavage peptide deduction to NetMHCpan."""
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from tools import cleavage


DEDUCED_COLUMNS = ["CHROM", "POS", "Peptide_id", "Peptide", "hla_peptides"]
KMER_A = "ACDEFGHIK"
KMER_B = "LMNPQRSTV"
KMER_C = "WYACDEFGH"


@pytest.fixture
def original_peptides():
    return pd.DataFrame({
        "CHROM": ["1", "1", "2"],
        "POS": [100, 200, 300],
        "Gene_id": ["ENSG_A", "ENSG_A", "ENSG_B"],
        "Transcript_id": ["ENST_A1", "ENST_A2", "ENST_B"],
        "Peptide_id": ["ENSP_A1", "ENSP_A2", "ENSP_B"],
        "Protein_position": [42, 55, 81],
        "peptide_REF": ["REFERENCE_A1", "REFERENCE_A2", "REFERENCE_B"],
        "peptide": ["ORIGINAL_A1", "ORIGINAL_A2", "ORIGINAL_B"],
        "hla_peptides": [[KMER_A], [KMER_A], [KMER_C]],
        "hla_peptides_REF": [["REFERENCE"], ["REFERENCE"], ["REFERENCE"]],
    })


@pytest.fixture
def prepare(tmp_path, original_peptides):
    def run(rows, length=9, pair="P01", columns=DEDUCED_COLUMNS):
        original_path = tmp_path / "original.pkl"
        original_peptides.to_pickle(original_path)
        deduced_path = tmp_path / "deduced.csv"
        pd.DataFrame(rows, columns=columns).to_csv(deduced_path, index=False)
        before = [original_path.read_bytes(), deduced_path.read_bytes()]
        output_dir = tmp_path / "run_tables"
        output_dir.mkdir()
        args = SimpleNamespace(length=str(length), pair=pair, run_name="cleavage_test")
        fasta, table = cleavage.prepare_cleavage_netmhcpan_inputs(
            str(output_dir), args, str(original_path), str(deduced_path)
        )
        assert [original_path.read_bytes(), deduced_path.read_bytes()] == before
        return Path(fasta), Path(table), pd.read_pickle(table)

    return run


def test_fasta_matches_grouped_table_and_preserves_cleavage_metadata(prepare, original_peptides):
    rows = [
        ["X", "900", "ENSP_B", "MWYACDEFGHK", KMER_C],
        ["3", "500", "ENSP_A1", "MACDEFGHIKLMNPQRSTVW", KMER_A],
        ["4", "600", "ENSP_A2", "MACDEFGHIKL", KMER_A],
        ["3", "500", "ENSP_A1", "MACDEFGHIKLMNPQRSTVW", KMER_B],
        ["3", "500", "ENSP_A1", "MACDEFGHIKLMNPQRSTVW", KMER_A],  # Exact duplicate.
    ]
    fasta, table_path, table = prepare(rows)
    expected = pd.DataFrame({
        "CHROM": ["X", "3", "4"], "POS": ["900", "500", "600"],
        "Gene_id": ["ENSG_B", "ENSG_A", "ENSG_A"],
        "Transcript_id": ["ENST_B", "ENST_A1", "ENST_A2"],
        "Peptide_id": ["ENSP_B", "ENSP_A1", "ENSP_A2"],
        "Protein_position": [81, 42, 55],
        "peptide_REF": [None, None, None],
        "peptide": ["MWYACDEFGHK", "MACDEFGHIKLMNPQRSTVW", "MACDEFGHIKL"],
        "hla_peptides": [[KMER_C], [KMER_A, KMER_B], [KMER_A]],
        "hla_peptides_REF": [[""], ["", ""], [""]],
    })
    pd.testing.assert_frame_equal(table.reset_index(drop=True), expected, check_dtype=False)
    assert list(table.columns) == list(original_peptides.columns)
    # Same k-mer under the same gene is written once despite multiple proteins.
    assert fasta.read_text() == (
        f">ENSG_B:1\n{KMER_C}\n>ENSG_A:1\n{KMER_A}\n>ENSG_A:2\n{KMER_B}\n"
    )
    assert table_path.with_suffix(".tsv").read_text() == expected.to_csv(sep="\t", index=False)
    assert not cleavage.fasta_is_empty(str(fasta))


@pytest.mark.parametrize("length", [8, 9, 10])
def test_keeps_only_requested_kmer_length(prepare, length):
    sequence = "ACDEFGHIKLM"
    selected = sequence[:length]
    rows = [
        ["1", "100", "ENSP_A1", sequence, sequence[:length - 1]],
        ["1", "100", "ENSP_A1", sequence, selected],
        ["1", "100", "ENSP_A1", sequence, sequence[:length + 1]],
    ]
    fasta, _, table = prepare(rows, length=length)
    assert table["hla_peptides"].tolist() == [[selected]]
    assert table["peptide"].tolist() == [sequence]
    assert fasta.read_text() == f">ENSG_A:1\n{selected}\n"


@pytest.mark.parametrize("rows", [
    [],
    [["1", "100", "ENSP_A1", KMER_A, "ACDE"], ["1", "100", "ENSP_A1", KMER_A, "ACDEFGHIKL"]],
    [["1", "100", "ENSP_UNKNOWN", KMER_A, KMER_A]],
    [["1", "100", None, KMER_A, KMER_A], ["1", "100", "ENSP_A1", KMER_A, None]],
], ids=["empty-input", "wrong-lengths", "unknown-protein", "missing-required-values"])
def test_no_usable_peptides_produce_consistent_empty_outputs(prepare, original_peptides, rows):
    fasta, table_path, table = prepare(rows)
    pd.testing.assert_frame_equal(table, original_peptides.iloc[0:0])
    assert fasta.read_bytes() == b""
    assert cleavage.fasta_is_empty(str(fasta))
    assert table_path.with_suffix(".tsv").read_text() == original_peptides.iloc[0:0].to_csv(sep="\t", index=False)


def test_discards_unmapped_proteins_without_losing_valid_peptides(prepare):
    rows = [
        ["1", "100", "ENSP_UNKNOWN", KMER_C, KMER_C],
        ["2", "200", "ENSP_A1", KMER_A, KMER_A],
    ]
    fasta, _, table = prepare(rows)
    assert table["Peptide_id"].tolist() == ["ENSP_A1"]
    assert table["Gene_id"].tolist() == ["ENSG_A"]
    assert table["CHROM"].tolist() == ["2"]
    assert table["POS"].tolist() == ["200"]
    assert fasta.read_text() == f">ENSG_A:1\n{KMER_A}\n"


def test_identical_kmers_in_distinct_genes_remain_distinct_fasta_entries(prepare):
    fasta, _, table = prepare([
        ["1", "100", "ENSP_A1", KMER_A, KMER_A],
        ["2", "200", "ENSP_B", KMER_A, KMER_A],
    ])
    assert table["Gene_id"].tolist() == ["ENSG_A", "ENSG_B"]
    assert table["hla_peptides"].tolist() == [[KMER_A], [KMER_A]]
    assert fasta.read_text() == f">ENSG_A:1\n{KMER_A}\n>ENSG_B:1\n{KMER_A}\n"


def test_shared_gene_kmer_keeps_all_genomic_positions_in_table(prepare):
    fasta, _, table = prepare([
        ["1", "100", "ENSP_A1", KMER_A, KMER_A],
        ["X", "900", "ENSP_A1", KMER_A, KMER_A],
    ])
    assert table["Peptide_id"].tolist() == ["ENSP_A1", "ENSP_A1"]
    assert table["CHROM"].tolist() == ["1", "X"]
    assert table["POS"].tolist() == ["100", "900"]
    assert table["hla_peptides"].tolist() == [[KMER_A], [KMER_A]]
    assert fasta.read_text() == f">ENSG_A:1\n{KMER_A}\n"


def test_missing_long_peptide_column_falls_back_to_kmer(prepare):
    fasta, _, table = prepare(
        [["1", "100", "ENSP_A1", KMER_A]],
        columns=["CHROM", "POS", "Peptide_id", "hla_peptides"],
    )
    assert table["peptide"].tolist() == [KMER_A]
    assert table["hla_peptides"].tolist() == [[KMER_A]]
    assert fasta.read_text() == f">ENSG_A:1\n{KMER_A}\n"


@pytest.mark.parametrize("pair", ["", "P01", "P02"])
def test_output_paths_are_scoped_to_run_and_pair(prepare, tmp_path, pair):
    fasta, table_path, _ = prepare([["1", "100", "ENSP_A1", KMER_A, KMER_A]], pair=pair)
    prefix = f"{pair}_" if pair else ""
    stem = prefix + "cleavage_test"
    output_dir = tmp_path / "run_tables"
    assert fasta == output_dir / f"{stem}_cleavage_deduced_fasta.fa"
    assert table_path == output_dir / f"{stem}_pep_df_cleavage_filtered.pkl"
    assert {path.name for path in output_dir.iterdir()} == {
        fasta.name, table_path.name, f"{stem}_pep_df_cleavage_filtered.tsv",
    }
