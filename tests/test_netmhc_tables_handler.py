#coding:utf-8
"""
Tests for netmhc_tables_handler.py module (NetMHCpan output processing)
"""
import pandas as pd
import pytest

from tools import netmhc_tables_handler


@pytest.fixture
def netmhc_output(tmp_path):
    # One column group per HLA, followed by the two shared summary columns.
    columns = [
        "Pos", "Peptide", "ID",
        "HLA-A02:01", "A_icore", "A_score", "A_rank",
        "HLA-B07:02", "B_icore", "B_score", "B_rank", "Ave", "NB",
    ]
    table = pd.DataFrame([
        ["Pos", "Peptide", "ID", "core", "icore", "EL-score", "EL_Rank",
         "core", "icore", "EL-score", "EL_Rank", "Ave", "NB"],
        [1, "ACDEFGHIK", "GENE1", "ACDEFGHIK", "ACDEFGHIK", 0.8, 0.5,
         "ACDEFGHIK", "ACDEFGHIK", 0.1, 4.0, 0.45, 1],
        [2, "LMNPQRSTV", "GENE2", "LMNPQRSTV", "LMNPQRSTV", 0.2, 3.0,
         "LMNPQRSTV", "LMNPQRSTV", 0.9, 0.2, 0.55, 1],
    ], columns=columns)
    path = tmp_path / "netmhc.out"
    table.to_csv(path, sep="\t", index=False)
    return path, columns


class TestFindSubsets:
    """Tests for find_subsets() - finding HLA subset groups"""
    
    def test_finds_hla_alleles(self, netmhc_output):
        """Test detection of HLA allele groups"""
        netmhc_file, columns = netmhc_output
        args = type("Args", (), {"hla_typing": "HLA-A02:01,HLA-B07:02"})()
        subsets, table = netmhc_tables_handler.find_subsets(str(netmhc_file), args)
        assert subsets == [columns[3:7], columns[7:11]]
        assert table.columns.tolist() == columns
        assert table["Peptide"].tolist() == ["Peptide", "ACDEFGHIK", "LMNPQRSTV"]

    def test_counts_correct_alleles(self, netmhc_output):
        """Test correct counting of alleles"""
        netmhc_file, _ = netmhc_output
        args = type("Args", (), {"hla_typing": "HLA-A02:01,HLA-B07:02"})()
        result = netmhc_tables_handler.find_subsets(str(netmhc_file), args)
        subsets, _ = result
        
        assert len(subsets) == 2
        assert [subset[0] for subset in subsets] == args.hla_typing.split(",")
        assert all(len(subset) == 4 for subset in subsets)
        assert {column for subset in subsets for column in subset}.isdisjoint(
            {"Pos", "Peptide", "ID", "Ave", "NB"}
        )


class TestFormatNetMHCpan:
    """Tests for format_netMHCpan() - formatting output"""
    
    def test_renames_columns(self, netmhc_output):
        """Test that columns are correctly renamed"""
        path, columns = netmhc_output
        netmhc_table = pd.read_csv(path, sep="\t", dtype=str)
        subsets = [columns[3:7], columns[7:11]]
        args = type("Args", (), {"class_type": 1})()
        result = netmhc_tables_handler.format_netMHCpan(netmhc_table, subsets, args)
        
        expected = pd.DataFrame({
            "Pos": ["1", "2", "1", "2"],
            "Peptide": ["ACDEFGHIK", "LMNPQRSTV", "ACDEFGHIK", "LMNPQRSTV"],
            "ID": ["GENE1", "GENE2", "GENE1", "GENE2"],
            "HLA": ["HLA-A02:01", "HLA-A02:01", "HLA-B07:02", "HLA-B07:02"],
            "EL-score": ["0.8", "0.2", "0.1", "0.9"],
            "EL_Rank": ["0.5", "3.0", "4.0", "0.2"],
        })
        pd.testing.assert_frame_equal(result[expected.columns], expected, check_names=False)
