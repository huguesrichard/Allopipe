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
        ["Pos", "Peptide", "ID", "core", "icore", "EL_score", "EL_Rank",
         "core", "icore", "EL_score", "EL_Rank", "Ave", "NB"],
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
        args = type("Args", (), {"hla_typing": "HLA-A02:01,HLA-B07:02", "class_type": 1})()
        subsets, table = netmhc_tables_handler.find_subsets(str(netmhc_file), args)
        assert subsets == [columns[3:7], columns[7:11]]
        assert table.columns.tolist() == columns
        assert table["Peptide"].tolist() == ["Peptide", "ACDEFGHIK", "LMNPQRSTV"]

    def test_counts_correct_alleles(self, netmhc_output):
        """Test correct counting of alleles"""
        netmhc_file, _ = netmhc_output
        args = type("Args", (), {"hla_typing": "HLA-A02:01,HLA-B07:02", "class_type": 1})()
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
            "EL_score": ["0.8", "0.2", "0.1", "0.9"],
            "EL_Rank": ["0.5", "3.0", "4.0", "0.2"],
        })
        pd.testing.assert_frame_equal(result[expected.columns], expected, check_names=False)


def test_handle_netmhcpan_42e_output(tmp_path):
    """Parse and filter the two-header, multi-allele output of NetMHCpan 4.2e."""
    path = tmp_path / "netmhc.out"
    path.write_text(
        "#netMHCpan-4.2 -BA -xls\n"
        "\t\t\tHLA-A02:01\t\t\t\t\t\tHLA-B07:02\t\t\t\t\t\t\t\n"
        "Pos\tPeptide\tID\tcore\ticore\tEL_score\tEL_rank\tBA_score\tBA_rank"
        "\tcore\ticore\tEL_score\tEL_rank\tBA_score\tBA_rank\tAve\tNB\n"
        "1\tACDEFGHIK\tGENE1\tACDEFGHIK\tACDEFGHIK\t0.8\t0.5\t0.3\t1.0"
        "\tACDEFGHIK\tACDEFGHIK\t0.1\t4.0\t0.2\t3.0\t0.45\t1\n"
        "2\tLMNPQRSTV\tGENE2\tLMNPQRSTV\tLMNPQRSTV\t0.2\t3.0\t0.1\t5.0"
        "\tLMNPQRSTV\tLMNPQRSTV\t0.9\t0.2\t0.4\t0.5\t0.55\t1\n"
    )
    args = type("Args", (), {
        "hla_typing": "HLA-A02:01,HLA-B07:02", "class_type": 1, "el_rank": 2.0,
    })()

    result = netmhc_tables_handler.handle_netMHCpan(str(path), args)

    assert list(zip(result["Peptide"], result["HLA"])) == [
        ("ACDEFGHIK", "HLA-A02:01"), ("LMNPQRSTV", "HLA-B07:02"),
    ]
    assert result["EL_score"].tolist() == [0.8, 0.9]
    assert result["EL_Rank"].tolist() == [0.5, 0.2]
    assert (tmp_path / "netmhc_netmhc_ELR_2.0.csv").exists()
    for name in ("netmhc_netmhc_ELR_2.0.csv", "netmhc_netmhc_ELR_nofilter.csv"):
        columns = pd.read_csv(tmp_path / name).columns
        assert "EL_score" in columns
        assert "EL-score" not in columns


def test_handle_netmhciipan_43k_output(tmp_path):
    path = tmp_path / "netmhc.out"
    path.write_text(
        "#NetMHCIIpan-4.3k command line\n"
        "\t\t\t\tDRB1_0401\t\t\t\t\t\t\tHLA-DPA10103-DPB10301\t\t\t\t\t\t\t\t\n"
        "Pos\tPeptide\tID\tTarget\tCore\tInverted\tScore_EL\tRank_EL\tScore_BA\tnM\tRank_BA"
        "\tCore\tInverted\tScore_EL\tRank_EL\tScore_BA\tnM\tRank_BA\tAve\tNB\n"
        "1\tGLEEFLSAMQSARGP\tGENE1\tNA\tLSAMQSARG\t0\t0.35\t3.9\t0.57\t105\t2.2"
        "\tEEFLSAMQS\t0\t0.002\t37.3\t0.34\t1260\t15.7\t0.17\t1\n"
    )
    args = type("Args", (), {
        "hla_typing": "DRB1_0401,HLA-DPA10103-DPB10301", "class_type": 2, "el_rank": 10,
    })()

    result = netmhc_tables_handler.handle_netMHCpan(str(path), args)

    assert result[["Peptide", "ID", "HLA"]].values.tolist() == [
        ["GLEEFLSAMQSARGP", "GENE1", "DRB1_0401"]
    ]
    assert result["Score"].tolist() == [0.35]
    assert result["Rank"].tolist() == [3.9]
    assert result["NB"].tolist() == [1]
    assert "Target" not in result.columns


def test_class_ii_filters_rank_score_and_shared_nb(tmp_path):
    table = pd.DataFrame({
        "Pos": [1, 2, 3, 4, 5],
        "Peptide": ["strong", "weak", "rank8_with_binder", "rank8_no_binder", "zero_score"],
        "HLA": ["DRB1_0401"] * 5,
        "Score": [0.5, 0.4, 0.3, 0.3, 0],
        "Rank": [0.5, 3, 8, 8, 2],
        "NB": [1, 1, 1, 0, 1],
    })
    args = type("Args", (), {"class_type": 2, "el_rank": 10})()
    result = netmhc_tables_handler.filter_netMHC_table(table, args, str(tmp_path / "test.out"))
    assert result["Peptide"].tolist() == ["strong", "weak", "rank8_with_binder"]
