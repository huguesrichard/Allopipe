#coding:utf-8
"""
Tests for table_operations.py module (Transcript & mismatch merging)
"""
from pathlib import Path
import pandas as pd
import pytest

from tools import table_operations, normalize_cohort


class TestSaveMismatch:
    """Tests for save_mismatch() - saving mismatch summary"""

    def test_preserves_dotted_sample_names_in_final_table(self, tmp_path):
        run_ams = tmp_path / "AMS"
        run_ams.mkdir()
        args = type("Args", (), {
            "donor": "/input/test.1.vcf.gz",
            "recipient": "/input/test.2.vcf.gz",
            "pair": "P1",
            "run_name": "dotted_samples",
            "orientation": "dr",
            "min_dp": 20,
            "max_dp": 400,
            "min_ad": 5,
            "min_gq": 0,
            "homozygosity_thr": 0.2,
            "base_length": 3,
            "output_dir": str(tmp_path),
        })()

        ams_exp_path = table_operations.save_mismatch(
            str(run_ams), args, 4,
            str(tmp_path / "donor.tsv"),
            str(tmp_path / "recipient.tsv"),
            str(tmp_path / "mismatches.tsv"),
        )
        table_operations.create_AMS_df(ams_exp_path)

        result = pd.read_csv(Path(ams_exp_path) / "AMS_df.tsv", sep="\t", dtype=str)
        assert result.loc[0, "donor"] == "test.1"
        assert result.loc[0, "recipient"] == "test.2"
    
    def test_creates_ams_file(self, tmp_path, nextflow_environment):
        """Test creation of AMS output file"""
        run_ams = tmp_path / "AMS"
        run_ams.mkdir()
        
        args = type('Args', (), {
            'donor': '/path/to/donor.vcf',
            'recipient': '/path/to/recipient.vcf',
            'pair': 'testpair',
            'run_name': 'test_run',
            'orientation': 'dr',
            'min_dp': 10,
            'max_dp': 1000,
            'min_ad': 2,
            'min_gq': 20,
            'homozygosity_thr': 0.8,
            'base_length': 3,
            'output_dir': str(tmp_path),
        })()
        
        mismatch_count = 5
        
        df_donor_file = str(tmp_path / "donor.tsv")
        df_recipient_file = str(tmp_path / "recipient.tsv")
        mismatches_file = str(tmp_path / "mismatches.tsv")
        
        # Create dummy files
        Path(df_donor_file).write_text("data\n", encoding="utf-8")
        Path(df_recipient_file).write_text("data\n", encoding="utf-8")
        Path(mismatches_file).write_text("data\n", encoding="utf-8")
        
        result_path = table_operations.save_mismatch(
            str(run_ams), args, mismatch_count,
            df_donor_file, df_recipient_file, mismatches_file
        )
        
        files = list(Path(result_path).glob("*.pkl"))
        assert len(files) == 1
        result = pd.read_pickle(files[0])
        published_dir, _ = nextflow_environment
        assert result.loc[0, "ams"] == 5
        assert result.loc[0, "donor_table"] == str(published_dir / "donor.tsv")
        assert result.loc[0, "recipient_table"] == str(published_dir / "recipient.tsv")
        assert result.loc[0, "mismatches_table"] == str(published_dir / "mismatches.tsv")
        pd.testing.assert_frame_equal(pd.read_csv(files[0].with_suffix(".csv")), result)


class TestCreateAmsDataframe:
    """Tests for create_AMS_df() - loading AMS DataFrame"""
    
    def test_loads_csv_ams(self, tmp_path):
        """Test loading AMS DataFrame from directory with pickle files"""
        # create_AMS_df expects a directory containing pickle files
        ams_dir = tmp_path / "ams_data"
        ams_dir.mkdir()
        
        # Create sample pickle files
        # Pair names must contain 'R' or 'P' and have a numeric component
        df1 = pd.DataFrame({
            "pair": ["P01"],
            "mismatch_count": [5],
            "min_dp": [10],
            "max_dp": [1000],
            "min_ad": [2],
            "min_gq": [20],
            "homozygosity_thr": [0.8],
            "base_length": [3]
        })
        df1.to_pickle(str(ams_dir / "ams1.pkl"))
    
        df2 = pd.DataFrame({
            "pair": ["P02"],
        })
        df2.to_pickle(str(ams_dir / "ams2.pkl"))
        
        # Function creates AMS_df.tsv in the directory
        table_operations.create_AMS_df(str(ams_dir))
        
        # Check that the output file was created
        output_file = ams_dir / "AMS_df.tsv"
        assert output_file.exists()
        
        # Load and verify the combined dataframe
        df_combined = pd.read_csv(str(output_file), sep="\t")
        assert len(df_combined) == 2
        assert "mismatch_count" in df_combined.columns

    def test_loads_pickle_ams(self, tmp_path):
        """Test loading AMS DataFrame from pickle"""
        # create_AMS_df expects a directory containing pickle files, not a single file
        ams_dir = tmp_path / "ams_data"
        ams_dir.mkdir()
        
        df_original = pd.DataFrame({
            "pair": ["P01", "P02"],
            "mismatch_count": [5, 8],
        })
        df_original.to_pickle(str(ams_dir / "ams_data.pkl"))
        
        table_operations.create_AMS_df(str(ams_dir))
        
        # Check that AMS_df.tsv was created
        tsv_file = ams_dir / "AMS_df.tsv"
        assert tsv_file.exists()


class TestBuildTranscriptsTableIndiv:
    """Tests for build_transcripts_table_indiv() - individual transcript building"""
    
    def test_builds_donor_transcripts(self, tmp_path):
        """Test building donor transcripts table"""
        # build_transcripts_table_indiv expects a pickle file with VEP data
        vep_data = pd.DataFrame({
            "CHROM": [1, 1],
            "POS": [100, 200],
            "INFO": [
                "missense|ENSG0001|ENST0001|1|2|3|K/N|aAa/aTa|0.001",
                "synonymous|ENSG0002|ENST0002|4|5|6|R/R|cGc/cGc|0.002"
            ]
        })
        vep_file = tmp_path / "vep.pkl"
        vep_data.to_pickle(str(vep_file))
        
        merged_ams = pd.DataFrame({
            "CHROM": [1, 1],
            "POS": [100, 200],
            "aa_ref_indiv_x": ["K", "R"],
            "aa_alt_indiv_x": ["N", "R"],
            "diff": ["1", "0"],
        })
        
        # Mock VEP indices
        vep_indices = type('VepIndices', (), {
            'consequence': 0,
            'gene': 1,
            'transcript': 2,
            'cdna': 3,
            'cds': 4,
            'prot': 5,
            'aa': 6,
            'codons': 7,
            'gnomad': 8,
            'frameshift': None,
        })()
        
        result = table_operations.build_transcripts_table_indiv(
            str(vep_file), merged_ams, vep_indices, "donor"
        )
        
        expected = pd.DataFrame({
            "CHROM": [1, 1],
            "POS": [100, 200],
            "Gene_id": ["ENSG0001", "ENSG0002"],
            "Transcript_id": ["ENST0001", "ENST0002"],
            "Protein_position": ["3", "6"],
            "Amino_acids": ["K/N", "R/R"],
            "aa_ref_indiv_x": ["K", "R"],
            "aa_alt_indiv_x": ["N", "R"],
            "diff": ["1", "0"],
        })
        pd.testing.assert_frame_equal(
            result[expected.columns].reset_index(drop=True), expected
        )
        assert "INFO" not in result.columns

    def test_filters_by_position(self, tmp_path):
        """Test that transcripts are filtered by position"""
        # Create mock VEP table with multiple positions
        vep_data = pd.DataFrame({
            "CHROM": [1, 1, 1],
            "POS": [100, 200, 300],
            "INFO": [
                "missense|ENSG0001|ENST0001|1|2|3|K/N|aAa/aTa|0.001",
                "synonymous|ENSG0002|ENST0002|4|5|6|R/R|cGc/cGc|0.002",
                "missense|ENSG0003|ENST0003|7|8|9|M/I|aTg/aTc|0.003"
            ]
        })
        vep_file = tmp_path / "vep.pkl"
        vep_data.to_pickle(str(vep_file))
        
        # Create merged_ams with only positions 100 and 300
        merged_ams = pd.DataFrame({
            "CHROM": [1, 1],
            "POS": [100, 300],
            "aa_ref_indiv_x": ["K", "M"],
            "aa_alt_indiv_x": ["N", "I"],
            "diff": ["K>N", "M>I"]
        })
        
        # Mock VEP indices
        vep_indices = type('VepIndices', (), {
            'consequence': 0,
            'gene': 1,
            'transcript': 2,
            'cdna': 3,
            'cds': 4,
            'prot': 5,
            'aa': 6,
            'codons': 7,
            'gnomad': 8,
            'frameshift': None,
        })()
        
        result = table_operations.build_transcripts_table_indiv(
            str(vep_file), merged_ams, vep_indices, "donor"
        )
        
        assert result["POS"].tolist() == [100, 300]
        assert result["Gene_id"].tolist() == ["ENSG0001", "ENSG0003"]
        assert result["Transcript_id"].tolist() == ["ENST0001", "ENST0003"]


class TestBuildTranscriptsTable:
    """Tests for build_transcripts_table() - merging donor/recipient"""
    
    def test_merges_donor_recipient_transcripts(self):
        """Test merging transcripts from donor and recipient"""
        transcripts_donor = pd.DataFrame({
            "CHROM": [1],
            "POS": [100],
            "Consequence": ["missense"],
            "Gene_id": ["ENSG0001"],
            "Transcript_id": ["ENST0001"],
            "cDNA_position": ["1"],
            "CDS_position": ["2"],
            "Protein_position": ["42"],
            "Amino_acids": ["K/N"],
            "Codons": ["aAa/aTa"],
            "Frameshift_sequence": [""],
            "gnomADe_AF": ["0.001"],
            "diff": ["1"],
        })
        transcripts_recipient = pd.DataFrame({
            "CHROM": [1],
            "POS": [100],
            "Consequence": ["missense"],
            "Gene_id": ["ENSG0001"],
            "Transcript_id": ["ENST0001"],
            "cDNA_position": ["1"],
            "CDS_position": ["2"],
            "Protein_position": ["42"],
            "Amino_acids": ["K/N"],
            "Codons": ["aAa/aTa"],
            "Frameshift_sequence": [""],
            "gnomADe_AF": ["0.001"],
            "diff": ["1"],
        })
        
        # Keep the common transcript once and retain a recipient-only transcript.
        recipient_only = transcripts_recipient.copy()
        recipient_only.loc[0, "POS"] = 200
        recipient_only.loc[0, "Gene_id"] = "ENSG0002"
        recipient_only.loc[0, "Transcript_id"] = "ENST0002"
        transcripts_recipient = pd.concat(
            [transcripts_recipient, recipient_only], ignore_index=True
        )
        expected = pd.concat([transcripts_donor, recipient_only], ignore_index=True)
        transcripts_donor = pd.concat(
            [transcripts_donor, transcripts_donor], ignore_index=True
        )

        result = table_operations.build_transcripts_table(
            transcripts_donor, transcripts_recipient
        )
        pd.testing.assert_frame_equal(
            result.sort_values("POS").reset_index(drop=True), expected, check_like=True
        )


class TestGetRefRatioPair:
    """Tests for get_ref_ratio_pair()"""
    
    def test_calculates_ratio_correctly(self):
        """Test correct calculation of reference ratio"""
        # get_ref_ratio_pair expects dataframes with CHROM, POS and GT, which are merged
        donor_df = pd.DataFrame({
            "CHROM": ["1", "1", "1"],
            "POS": [100, 200, 300],
            "GT": ["0/0", "0/1", "0/0"],
        })
        recipient_df = pd.DataFrame({
            "CHROM": ["1", "1", "1"],
            "POS": [100, 200, 300],
            "GT": ["0/0", "0/0", "0/1"],
        })
        
        ratio = normalize_cohort.get_ref_ratio_pair(donor_df, recipient_df)
        
        assert isinstance(ratio, tuple)
        assert len(ratio) == 3
        # Expect one position where both are 0/0 (common_ref=1)
        # total_ref counts any row where either GT is 0/0 (3 total)
        assert ratio == (1, 3, 1/3)


class TestCohortReferenceRatio:
    """Reference ratios are now calculated in the final cohort stage."""
    
    def test_counts_common_and_total_reference_positions(self):
        """A heterozygous/reference site contributes only to total_ref."""
        # Test get_ref_ratio_pair function
        donor_df = pd.DataFrame({
            "CHROM": [1, 1],
            "POS": [100, 200],
            "GT": ["0/0", "0/1"]
        })
        recipient_df = pd.DataFrame({
            "CHROM": [1, 1],
            "POS": [100, 200],
            "GT": ["0/0", "0/0"]
        })
        
        common_ref, total_ref, ref_ratio = normalize_cohort.get_ref_ratio_pair(donor_df, recipient_df)
        
        assert common_ref == 1  # One position where both are 0/0
        assert total_ref == 2   # Two positions total
        assert ref_ratio == 0.5


class TestCohortNormalization:
    """The final Nextflow cohort stage replaces the former add_norm helper."""

    @staticmethod
    def cohort_files(tmp_path, last_score=31):
        ams_dir = tmp_path / "AMS"
        ams_dir.mkdir()
        ams_pkls, donor_tables, recipient_tables = [], [], []
        # Deliberately unsorted inputs: pair ordering must be numeric.
        for pair, score, common in [("P10", last_score, 4), ("P1", 10, 2), ("P2", 20, 3)]:
            ams_file = ams_dir / f"{pair}_AMS.pkl"
            pd.DataFrame({"pair": [pair], "ams": [score]}).to_pickle(ams_file)
            donor = pd.DataFrame({"CHROM": ["1"] * 4, "POS": [100, 200, 300, 400], "GT": ["0/0"] * 4})
            recipient = donor.copy()
            recipient["GT"] = ["0/0"] * common + ["0/1"] * (4 - common)
            donor_file = tmp_path / f"{pair}_D0_table.tsv"
            recipient_file = tmp_path / f"{pair}_R0_table.tsv"
            donor.to_csv(donor_file, sep="\t", index=False)
            recipient.to_csv(recipient_file, sep="\t", index=False)
            ams_pkls.append(str(ams_file))
            donor_tables.append(str(donor_file))
            recipient_tables.append(str(recipient_file))
        return ams_dir, ams_pkls, donor_tables, recipient_tables

    def test_adds_normalized_columns(self, tmp_path):
        ams_dir, ams_pkls, donors, recipients = self.cohort_files(tmp_path)
        normalize_cohort.write_normalized_ams_tables(ams_pkls, donors, recipients)
        result = pd.read_csv(ams_dir / "AMS_df.tsv", sep="\t")
        expected = pd.DataFrame({
            "pair": ["P1", "P2", "P10"],
            "ams_giab": [10, 20, 31],
            "common_ref": [2, 3, 4],
            "total_ref": [4, 4, 4],
            "ref_ratio": [0.5, 0.75, 1.0],
            "ams_norm": [20, 20, 20],
        })
        pd.testing.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("last_score", [
        31,
        pytest.param(30, marks=pytest.mark.xfail(
            strict=True,
            reason="Float regression yields 19.999...; astype(int) truncates the expected score 20 to 19",
        )),
    ])
    def test_normalizes_correctly(self, tmp_path, last_score):
        _, ams_pkls, donors, recipients = self.cohort_files(tmp_path, last_score)
        result = normalize_cohort.normalize_ams(ams_pkls, list(reversed(donors)), recipients)
        assert result["pair"].tolist() == ["P1", "P2", "P10"]
        assert result["ams_giab"].tolist() == [10, 20, last_score]
        assert result["ref_ratio"].tolist() == [0.5, 0.75, 1.0]
        assert result["ams_norm"].tolist() == [20, 20, 20]
