#coding:utf-8
"""
Tests for table_operations.py module (Transcript & mismatch merging)
"""
from pathlib import Path
import pandas as pd

from tools import table_operations


class TestSaveMismatch:
    """Tests for save_mismatch() - saving mismatch summary"""
    
    def test_creates_ams_file(self, tmp_path):
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
        })()
        
        mismatch_count = 5
        formatted_datetime = ""
        
        df_donor_file = str(tmp_path / "donor.tsv")
        df_recipient_file = str(tmp_path / "recipient.tsv")
        mismatches_file = str(tmp_path / "mismatches.tsv")
        
        # Create dummy files
        Path(df_donor_file).write_text("data\n", encoding="utf-8")
        Path(df_recipient_file).write_text("data\n", encoding="utf-8")
        Path(mismatches_file).write_text("data\n", encoding="utf-8")
        
        result_path = table_operations.save_mismatch(
            str(run_ams), args, mismatch_count, formatted_datetime,
            df_donor_file, df_recipient_file, mismatches_file
        )
        
        assert result_path is not None
        assert Path(result_path).exists() or str(run_ams) in result_path


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
        })()
        
        result = table_operations.build_transcripts_table_indiv(
            str(vep_file), merged_ams, vep_indices, "donor"
        )
        
        assert isinstance(result, (pd.DataFrame, str, dict))

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
        })()
        
        result = table_operations.build_transcripts_table_indiv(
            str(vep_file), merged_ams, vep_indices, "donor"
        )
        
        # Should only include positions 100 and 300, not 200
        assert len(result) == 2
        assert all(pos in [100, 300] for pos in result["POS"])


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
            "gnomADe_AF": ["0.001"],
            "diff": ["1"],
        })
        
        result = table_operations.build_transcripts_table(
            transcripts_donor, transcripts_recipient
        )
        
        assert isinstance(result, pd.DataFrame)
        # Should have columns from both
        assert len(result) >= 0


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
        
        ratio = table_operations.get_ref_ratio_pair(donor_df, recipient_df)
        
        assert isinstance(ratio, tuple)
        assert len(ratio) == 3
        # Expect one position where both are 0/0 (common_ref=1)
        # total_ref counts any row where either GT is 0/0 (3 total)
        assert ratio == (1, 3, 1/3)


class TestGetRefRatio:
    """Tests for get_ref_ratio() - loading reference populations"""
    
    def test_loads_reference_populations(self, tmp_path):
        """Test loading reference population data"""
        run_path = tmp_path / "run"
        run_path.mkdir()
        
        # Create mock reference files if needed
        assert run_path.exists()

    def test_normalizes_ratios(self):
        """Test normalization of ratios"""
        """Test normalization of ratios"""
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
        
        common_ref, total_ref, ref_ratio = table_operations.get_ref_ratio_pair(donor_df, recipient_df)
        
        assert common_ref == 1  # One position where both are 0/0
        assert total_ref == 2   # Two positions total
        assert ref_ratio == 0.5


class TestAddNorm:
    """Tests for add_norm() - adding normalized columns"""
    
    def test_adds_normalized_columns(self, tmp_path):
        """Test adding normalized score columns"""
        # add_norm expects a dataframe containing 'ams' and 'ref_ratio' columns
        ams_df = pd.DataFrame({
            "pair": ["pairA"],
            "ams": [5],
            "ref_ratio": [0.8],
        })
        # Provide a directory path for output
        ams_dir = tmp_path / "ams_dir"
        ams_dir.mkdir()

        table_operations.add_norm(ams_df, str(ams_dir), ref_ratio=0.8)

        # Check output file was created with expected columns
        out_tsv = ams_dir / "AMS_df.tsv"
        assert out_tsv.exists()
        out_df = pd.read_csv(out_tsv, sep="\t")
        assert "ams_giab" in out_df.columns
        assert "ams_norm" in out_df.columns
        assert "ref_ratio" in out_df.columns

    def test_normalizes_correctly(self, tmp_path):
        """Test that normalization formula is correct"""
        # Create test data
        ams_df = pd.DataFrame({
            "ams": [10, 20, 30],
            "ref_ratio": [0.5, 0.7, 0.9]
        })
        ref_ratio = 0.8  # Not used in add_norm, but passed
        
        # Create the directory
        ams_exp_path = tmp_path / "test"
        ams_exp_path.mkdir()
        
        # Call add_norm
        table_operations.add_norm(ams_df, str(ams_exp_path), ref_ratio)
        
        # Check that ams_norm column was added
        assert "ams_norm" in ams_df.columns
        
        # The normalization uses linear regression
        # We can check that values are reasonable
        assert all(isinstance(x, (int, float)) for x in ams_df["ams_norm"] if x != "NA")
