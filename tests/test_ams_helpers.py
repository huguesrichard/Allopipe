#coding:utf-8
"""
Comprehensive tests for ams_helpers.py module (VCF preparation & mismatch detection)
"""
from types import SimpleNamespace
from pathlib import Path
import pandas as pd
import numpy as np

from tools import ams_helpers


class TestCreateRunDirectory:
    """Tests for create_run_directory()"""
    
    def test_creates_all_subdirectories(self, tmp_path):
        """Test that all required subdirectories are created"""
        output_dir = str(tmp_path / "output")
        run_name = "test_run"
        
        run_path, run_tables, run_plots, run_ams, run_logs = ams_helpers.create_run_directory(
            run_name, output_dir
        )
        
        # All paths should exist
        assert Path(run_path).exists()
        assert Path(run_tables).exists()
        assert Path(run_plots).exists()
        assert Path(run_ams).exists()
        assert Path(run_logs).exists()
        
        # Correct structure
        assert run_name in run_path
        assert "run_tables" in run_tables
        assert "plots" in run_plots

    def test_idempotent_creation(self, tmp_path):
        """Test that calling twice doesn't cause errors"""
        output_dir = str(tmp_path / "output")
        run_name = "test_run"
        
        # First call
        ams_helpers.create_run_directory(run_name, output_dir)
        # Second call should not fail
        ams_helpers.create_run_directory(run_name, output_dir)
        
        assert Path(output_dir).exists()


class TestHandleOverwrite:
    """Tests for handle_overwrite()"""
    
    def test_returns_false_when_file_does_not_exist(self, tmp_path):
        """Test when AMS file doesn't exist"""
        output_dir = tmp_path / "output"
        output_dir.mkdir()
        run_name = "new_run"
        
        args = SimpleNamespace(
            output_dir=str(output_dir),
            run_name=run_name,
            pair="",
            orientation="dr",
            min_dp=10,
            max_dp=1000,
            min_ad=2,
            min_gq=20,
            homozygosity_thr=0.8,
            base_length=3,
        )
        
        # Should return False (no existing file)
        result = ams_helpers.handle_overwrite(args)
        assert result == False

    def test_returns_true_when_ams_file_exists(self, tmp_path):
        """Test when AMS file already exists"""
        output_dir = tmp_path / "output"
        run_name = "test_run"
        run_path, run_tables, run_plots, run_ams, _ = ams_helpers.create_run_directory(
            run_name, str(output_dir)
        )
        
        # Create an existing AMS file in the expected overwrite path
        overwrite_dir = Path(run_ams) / f"{run_name}_AMS_10_1000_2_0.8_20_dr_3"
        overwrite_dir.mkdir(parents=True, exist_ok=True)
        ams_file = overwrite_dir / f"AMS_{run_name}_10_1000_2_0.8_20_dr_3.csv"
        ams_file.write_text("pair,mismatch_count\n", encoding="utf-8")
        
        args = SimpleNamespace(
            output_dir=str(output_dir),
            run_name=run_name,
            pair="",
            orientation="dr",
            min_dp=10,
            max_dp=1000,
            min_ad=2,
            min_gq=20,
            homozygosity_thr=0.8,
            base_length=3,
        )
        
        result = ams_helpers.handle_overwrite(args)
        assert result == True


class TestWriteLog:
    """Tests for write_log()"""
    
    def test_log_file_created_with_fields(self, tmp_path):
        """Test that log file is created with proper fields"""
        output_dir = tmp_path / "output"
        run_name = "test_run"
        _, _, _, _, run_logs = ams_helpers.create_run_directory(run_name, str(output_dir))
        
        args = SimpleNamespace(
            donor="/path/to/donor.vcf",
            recipient="/path/to/recipient.vcf",
            orientation="dr",
            imputation="imputation",
            min_dp=10,
            max_dp=1000,
            min_ad=2,
            min_gq=20,
            homozygosity_thr=0.8,
            base_length=3,
            workers=4,
            norm_score=False,
            pair="",
            run_name=run_name,
        )
        
        ams_helpers.write_log(run_logs, args)
        
        # Check log file exists
        log_file = Path(run_logs) / "run.log"
        assert log_file.exists()
        
        # Check content
        content = log_file.read_text()
        assert "Orientation: dr" in content
        assert "Donor:" in content
        assert "Recipient:" in content


class TestExplodeGtAd:
    """Tests for explode_gt_ad()"""
    
    def test_valid_gt_ad_parsing(self):
        """Test parsing valid GT:AD format"""
        df = pd.DataFrame({
            "CHROM": ["1", "1", "1"],
            "POS": [100, 200, 300],
            "REF": ["A", "A", "A"],
            "ALT": ["T", "G", "C"],
            "GT": ["0/1", "0/0", "1/1"],
            "AD": ["100,50", "200,0", "0,150"],
        })

        result = ams_helpers.explode_gt_ad(df)

        # Should have exploded AD/GT entries and keep an index for the allele
        assert "ind" in result.columns
        assert all(result["GT"] == result["ind"])
        assert result["AD"].dtype == int

    def test_missing_ad_values(self):
        """Test with missing AD values"""
        df = pd.DataFrame({
            "CHROM": ["1", "1"],
            "POS": [100, 200],
            "REF": ["A", "A"],
            "ALT": ["T", "G"],
            "GT": ["0/1", "."],
            "AD": ["100,50", "."],
        })

        result = ams_helpers.explode_gt_ad(df)

        # Should handle missing AD/GT values gracefully by dropping those rows
        assert len(result) == 2
        assert (result["POS"] == 100).all()
        assert set(result["AD"].tolist()) == {100, 50}


class TestFilterOnDepth:
    """Tests for filter_on_depth()"""
    
    def test_filters_within_depth_range(self):
        """Test filtering by depth range"""
        df = pd.DataFrame({
            "CHROM": [1, 1, 1],
            "POS": [100, 200, 300],
            "DP": [50, 150, 500],  # Below, within, above min_dp=100, max_dp=400
        })
        
        df_indiv, subset = ams_helpers.filter_on_depth(df, df.copy(), min_dp=100, max_dp=400)

        assert len(df_indiv) == 1
        assert df_indiv["DP"].values[0] == 150
        assert len(subset) == 1
        assert subset["DP"].values[0] == 150

    def test_filters_removes_below_min(self):
        """Test removal of reads below minimum depth"""
        df = pd.DataFrame({
            "CHROM": [1, 1],
            "POS": [100, 200],
            "DP": [50, 150],
        })
        
        df_indiv, subset = ams_helpers.filter_on_depth(df, df.copy(), min_dp=100, max_dp=1000)

        assert len(df_indiv) == 1
        assert df_indiv["DP"].values[0] == 150
        assert len(subset) == 1
        assert subset["DP"].values[0] == 150

    def test_filters_removes_above_max(self):
        """Test removal of reads above maximum depth"""
        df = pd.DataFrame({
            "CHROM": [1, 1],
            "POS": [100, 200],
            "DP": [150, 500],
        })
        
        df_indiv, subset = ams_helpers.filter_on_depth(df, df.copy(), min_dp=100, max_dp=400)

        assert len(df_indiv) == 1
        assert df_indiv["DP"].values[0] == 150
        assert len(subset) == 1
        assert subset["DP"].values[0] == 150


class TestConvert:
    """Tests for convert() - homozygosity conversion"""
    
    def test_hetero_to_homo_conversion(self):
        """Test conversion of heterozygous to homozygous"""
        # Allelic depth (AD) indicates how much of the reference allele is present.
        # High ratio -> 0/0, low ratio -> 1/1.
        df = pd.DataFrame({
            "CHROM": [1, 1, 1],
            "POS": [100, 200, 300],
            "DP": [100, 100, 100],
            "GT": ["0/1", "0/1", "0/1"],
            "AD": ["90,10", "50,50", "10,90"],
        })

        result = ams_helpers.convert(df, df.copy(), homozygosity_threshold=0.8)

        # High ref:alt ratio should become 0/0; low ratio should become 1/1
        assert (result.loc[result["POS"] == 100, "GT"].iloc[0] == "0/0")
        assert (result.loc[result["POS"] == 300, "GT"].iloc[0] == "1/1")

    def test_threshold_boundary(self):
        """Test exactly at threshold"""
        df = pd.DataFrame({
            "CHROM": [1, 1],
            "POS": [100, 200],
            "DP": [100, 100],
            "GT": ["0/1", "0/1"],
            "AD": ["80,20", "79,21"],
        })

        result = ams_helpers.convert(df, df.copy(), homozygosity_threshold=0.8)

        # 0.8 is not < 0.8 but is > 0.2, so both positions should become 0/0
        assert set(result["GT"]) == {"0/0"}


class TestFilterOnGnomadAf:
    """Tests for filter_on_gnomad_af()"""
    
    def test_filters_by_allele_frequency(self):
        """Test filtering by gnomAD allele frequency"""
        df = pd.DataFrame({
            "CHROM": [1, 1, 1],
            "gnomADe_AF": [0.001, 0.05, 0.5],
        })
        
        result = ams_helpers.filter_on_gnomad_af(df, min_af=0.02)
        
        # Should keep only those >= 0.02
        assert len(result) == 2
        assert 0.001 not in result["gnomADe_AF"].values

    def test_handles_missing_af(self):
        """Test handling of missing AF values"""
        df = pd.DataFrame({
            "CHROM": [1, 1],
            "gnomADe_AF": [0.05, None],
        })
        
        result = ams_helpers.filter_on_gnomad_af(df, min_af=0.02)
        
        # Should handle NaN values gracefully
        assert len(result) >= 1


class TestMergeDfs:
    """Tests for merge_dfs() - merging donor/recipient DataFrames"""
    
    def test_merge_orientation_dr(self):
        """Test merge with dr (donor->recipient) orientation"""
        df_donor = pd.DataFrame({
            "CHROM": [1, 1],
            "POS": [100, 200],
            "REF": ["A", "G"],
            "ALT": ["T", "C"],
            "aa_REF": ["A", "G"],
            "aa_ALT": ["T", "C"],
        })
        df_recipient = pd.DataFrame({
            "CHROM": [1, 1],
            "POS": [100, 200],
            "REF": ["A", "G"],
            "ALT": ["A", "C"],
            "aa_REF": ["A", "G"],
            "aa_ALT": ["A", "C"],
        })
        
        merged, side, opposite = ams_helpers.merge_dfs(df_donor, df_recipient, "dr", "imputation")
        
        assert len(merged) > 0
        assert side == "x"
        assert opposite == "y"

    def test_merge_orientation_rd(self):
        """Test merge with rd (recipient->donor) orientation"""
        df_donor = pd.DataFrame({
            "CHROM": [1],
            "POS": [100],
            "REF": ["A"],
            "ALT": ["T"],
            "aa_REF": ["A"],
            "aa_ALT": ["T"],
        })
        df_recipient = pd.DataFrame({
            "CHROM": [1],
            "POS": [100],
            "REF": ["A"],
            "ALT": ["A"],
            "aa_REF": ["A"],
            "aa_ALT": ["A"],
        })
        
        merged, side, opposite = ams_helpers.merge_dfs(df_donor, df_recipient, "rd", "imputation")
        
        assert side == "y"
        assert opposite == "x"


class TestKeepAlt:
    """Tests for keep_alt() - keeping ALT variants"""
    
    def test_removes_ref_ref_vs_nan(self):
        """Test removal of REF/REF vs NaN positions"""
        merged_df = pd.DataFrame({
            "Ref_donor": ["A", "A"],
            "Alt_donor": ["A", "T"],
            "GT_donor": ["0/0", "0/1"],
            "Ref_recipient": ["A", "A"],
            "Alt_recipient": [np.nan, np.nan],
            "GT_recipient": [np.nan, "0/1"],
        })
        
        result = ams_helpers.keep_alt(merged_df, "donor", "recipient")
        
        # First row (A/A vs NaN) should be removed
        assert len(result) == 1
        assert result["Alt_donor"].values[0] != "A"

    def test_keeps_alt_positions(self):
        """Test that ALT positions are kept"""
        merged_df = pd.DataFrame({
            "Ref_donor": ["A", "A"],
            "Alt_donor": ["T", "G"],
            "GT_donor": ["0/1", "0/1"],
            "Ref_recipient": ["A", "A"],
            "Alt_recipient": ["A", "C"],
            "GT_recipient": ["0/1", "0/1"],
        })
        
        result = ams_helpers.keep_alt(merged_df, "donor", "recipient")
        
        # Both should be kept (have ALT)
        assert len(result) == 2


class TestCountMismatches:
    """Tests for count_mismatches()"""
    
    def test_counts_correct_mismatches(self):
        """Test correct mismatch counting"""
        merged_df = pd.DataFrame({
            "GT_x": ["0/0", "0/0", "0/0"],
            "GT_y": ["0/0", "0/0", "0/0"],
            "aa_indiv_x": ["T", "T", "A"],
            "aa_indiv_y": ["A", "C", "A"],
            "aa_REF": ["A", "A", "A"],
        })

        result, mismatch_count = ams_helpers.count_mismatches(merged_df, "dr")

        # Rows 1 and 2 have mismatches; row 3 does not
        assert mismatch_count == 2
        assert len(result) == 2

    def test_returns_zero_for_no_mismatch(self):
        """Test with no mismatches"""
        merged_df = pd.DataFrame({
            "GT_x": ["0/0"],
            "GT_y": ["0/0"],
            "aa_indiv_x": ["A"],
            "aa_indiv_y": ["A"],
            "aa_REF": ["A"],
        })

        result, mismatch_count = ams_helpers.count_mismatches(merged_df, "dr")

        assert mismatch_count == 0


class TestMismatchType:
    """Tests for mismatch_type() - hetero vs homo detection"""
    
    def test_detects_heterozygous_mismatch(self):
        """Test detection of heterozygous mismatch"""
        df = pd.DataFrame({
            "TYPE_x": ["heterozygous"],
            "TYPE_y": ["homozygous"],
        })
        
        result = ams_helpers.mismatch_type(df)
        
        assert "mismatch_type" in result.columns
        assert result["mismatch_type"].iloc[0] == "heterozygous"

    def test_detects_homozygous_mismatch(self):
        """Test detection of homozygous mismatch"""
        df = pd.DataFrame({
            "TYPE_x": ["homozygous"],
            "TYPE_y": ["homozygous"],
        })
        
        result = ams_helpers.mismatch_type(df)
        
        assert "mismatch_type" in result.columns
        assert result["mismatch_type"].iloc[0] == "homozygous"
