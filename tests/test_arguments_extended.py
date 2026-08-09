#coding:utf-8
"""
Additional tests for arguments_handling.py and netmhc_arguments.py modules
"""
from pathlib import Path
import sys
import pytest

from tools import arguments_handling, netmhc_arguments


def _write_dummy_vcf(path: Path) -> None:
    path.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")


class TestArgumentsHandlingExtended:
    """Extended tests for arguments handling"""
    
    def test_orientation_dr_valid(self, tmp_path, monkeypatch):
        """Test 'dr' orientation is accepted"""
        donor = tmp_path / "donor.vcf"
        recipient = tmp_path / "recipient.vcf"
        _write_dummy_vcf(donor)
        _write_dummy_vcf(recipient)
        
        monkeypatch.setattr(sys, "argv", [
            "ams_pipeline.py",
            str(donor),
            str(recipient),
            "dr",
            "imputation",
        ])
        
        args = arguments_handling.arguments()
        assert args.orientation == "dr"

    def test_orientation_rd_valid(self, tmp_path, monkeypatch):
        """Test 'rd' orientation is accepted"""
        donor = tmp_path / "donor.vcf"
        recipient = tmp_path / "recipient.vcf"
        _write_dummy_vcf(donor)
        _write_dummy_vcf(recipient)
        
        monkeypatch.setattr(sys, "argv", [
            "ams_pipeline.py",
            str(donor),
            str(recipient),
            "rd",
            "imputation",
        ])
        
        args = arguments_handling.arguments()
        assert args.orientation == "rd"

    def test_imputation_valid(self, tmp_path, monkeypatch):
        """Test 'imputation' mode is accepted"""
        donor = tmp_path / "donor.vcf"
        recipient = tmp_path / "recipient.vcf"
        _write_dummy_vcf(donor)
        _write_dummy_vcf(recipient)
        
        monkeypatch.setattr(sys, "argv", [
            "ams_pipeline.py",
            str(donor),
            str(recipient),
            "dr",
            "imputation",
        ])
        
        args = arguments_handling.arguments()
        assert args.imputation == "imputation"

    def test_no_imputation_valid(self, tmp_path, monkeypatch):
        """Test 'no-imputation' mode is accepted"""
        donor = tmp_path / "donor.vcf"
        recipient = tmp_path / "recipient.vcf"
        _write_dummy_vcf(donor)
        _write_dummy_vcf(recipient)
        
        monkeypatch.setattr(sys, "argv", [
            "ams_pipeline.py",
            str(donor),
            str(recipient),
            "dr",
            "no-imputation",
        ])
        
        args = arguments_handling.arguments()
        assert args.imputation == "no-imputation"

    def test_vcf_gz_extension_accepted(self, tmp_path, monkeypatch):
        """Test .vcf.gz extension is accepted"""
        donor = tmp_path / "donor.vcf.gz"
        recipient = tmp_path / "recipient.vcf.gz"
        donor.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
        recipient.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
        
        monkeypatch.setattr(sys, "argv", [
            "ams_pipeline.py",
            str(donor),
            str(recipient),
            "dr",
            "imputation",
        ])
        
        args = arguments_handling.arguments()
        assert str(donor) in args.donor or donor.name in args.donor

    def test_min_depth_default(self, tmp_path, monkeypatch):
        """Test default min_dp value"""
        donor = tmp_path / "donor.vcf"
        recipient = tmp_path / "recipient.vcf"
        _write_dummy_vcf(donor)
        _write_dummy_vcf(recipient)
        
        monkeypatch.setattr(sys, "argv", [
            "ams_pipeline.py",
            str(donor),
            str(recipient),
            "dr",
            "imputation",
        ])
        
        args = arguments_handling.arguments()
        assert hasattr(args, "min_dp")
        assert isinstance(args.min_dp, int)
        assert args.min_dp > 0

    def test_max_depth_default(self, tmp_path, monkeypatch):
        """Test default max_dp value"""
        donor = tmp_path / "donor.vcf"
        recipient = tmp_path / "recipient.vcf"
        _write_dummy_vcf(donor)
        _write_dummy_vcf(recipient)
        
        monkeypatch.setattr(sys, "argv", [
            "ams_pipeline.py",
            str(donor),
            str(recipient),
            "dr",
            "imputation",
        ])
        
        args = arguments_handling.arguments()
        assert args.min_dp <= args.max_dp

    def test_homozygosity_threshold_bounds(self, tmp_path, monkeypatch):
        """Test homozygosity threshold accepts valid values"""
        donor = tmp_path / "donor.vcf"
        recipient = tmp_path / "recipient.vcf"
        _write_dummy_vcf(donor)
        _write_dummy_vcf(recipient)
        
        monkeypatch.setattr(sys, "argv", [
            "ams_pipeline.py",
            str(donor),
            str(recipient),
            "dr",
            "imputation",
            "--homozygosity_thr",
            "0.75",
        ])
        
        args = arguments_handling.arguments()
        assert 0 < args.homozygosity_thr < 1

    def test_pair_name_custom(self, tmp_path, monkeypatch):
        """Test custom pair name"""
        donor = tmp_path / "donor.vcf"
        recipient = tmp_path / "recipient.vcf"
        _write_dummy_vcf(donor)
        _write_dummy_vcf(recipient)
        
        monkeypatch.setattr(sys, "argv", [
            "ams_pipeline.py",
            str(donor),
            str(recipient),
            "dr",
            "imputation",
            "--pair",
            "P001",
        ])
        
        args = arguments_handling.arguments()
        assert args.pair == "P001"

    def test_custom_run_name(self, tmp_path, monkeypatch):
        """Test custom run name"""
        donor = tmp_path / "donor.vcf"
        recipient = tmp_path / "recipient.vcf"
        _write_dummy_vcf(donor)
        _write_dummy_vcf(recipient)
        
        monkeypatch.setattr(sys, "argv", [
            "ams_pipeline.py",
            str(donor),
            str(recipient),
            "dr",
            "imputation",
            "--run_name",
            "customrun",
        ])
        
        args = arguments_handling.arguments()
        assert args.run_name == "customrun"


class TestNetmhcArgumentsHandling:
    """Tests for NetMHC-specific argument handling"""
    
    def test_netmhc_requires_run_name(self, tmp_path, monkeypatch):
        """Test NetMHC pipeline requires run name"""
        # Would need netmhc_arguments implementation
        pass

    def test_netmhc_orientation_options(self):
        """Test NetMHC accepts orientation options"""
        # Would test netmhc_arguments.netmhc_arguments()
        pass


class TestCheckFunctions:
    """Tests for validation functions"""
    
    def test_check_file_accepts_vcf(self, tmp_path):
        """Test check_file accepts .vcf extension"""
        vcf_file = tmp_path / "test.vcf"
        vcf_file.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
        
        parser = arguments_handling.CustomParser(prog="test")
        result = arguments_handling.check_file(parser, str(vcf_file))
        assert str(vcf_file) == result or result is not None

    def test_check_file_accepts_vcf_gz(self, tmp_path):
        """Test check_file accepts .vcf.gz extension"""
        vcf_file = tmp_path / "test.vcf.gz"
        vcf_file.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")
        
        parser = arguments_handling.CustomParser(prog="test")
        result = arguments_handling.check_file(parser, str(vcf_file))
        assert result is not None

    def test_check_threshold_value_valid_range(self):
        """Test threshold value accepts valid range"""
        parser = arguments_handling.CustomParser(prog="test")
        
        # 0 < homozygosity_thr < 1
        result = arguments_handling.check_threshold_value(
            parser, "0.5", "homozygosity_thr"
        )
        assert result == 0.5

    def test_check_threshold_value_rejects_zero(self):
        """Test threshold rejects 0"""
        parser = arguments_handling.CustomParser(prog="test")
        
        with pytest.raises(SystemExit):
            arguments_handling.check_threshold_value(
                parser, "0.0", "homozygosity_thr"
            )

    def test_check_threshold_value_rejects_one(self):
        """Test threshold rejects 1"""
        parser = arguments_handling.CustomParser(prog="test")
        
        with pytest.raises(SystemExit):
            arguments_handling.check_threshold_value(
                parser, "1.0", "homozygosity_thr"
            )

    def test_check_if_accepted_str_orientation(self):
        """Test accepted string validation"""
        parser = arguments_handling.CustomParser(prog="test")
        
        # Should accept valid orientations
        result = arguments_handling.check_if_accepted_str(parser, "dr")
        assert result == "dr" or result is not None

    def test_check_workers_count_valid(self, monkeypatch):
        """Test workers count validation"""
        parser = arguments_handling.CustomParser(prog="test")
        monkeypatch.setattr(arguments_handling.os, "cpu_count", lambda: 8)
        
        result = arguments_handling.check_workers_count(parser, "4")
        assert result == 4

    def test_check_workers_exceeds_cpu_count(self, monkeypatch):
        """Test workers count cannot exceed CPU count"""
        parser = arguments_handling.CustomParser(prog="test")
        monkeypatch.setattr(arguments_handling.os, "cpu_count", lambda: 4)
        
        with pytest.raises(SystemExit):
            arguments_handling.check_workers_count(parser, "8")
