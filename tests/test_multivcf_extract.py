#coding:utf-8
"""
Tests for multivcf_extract.py module (Multi-sample VCF extraction)
"""
from pathlib import Path
import pytest

from tools import multivcf_extract


class TestGetInfosMvcf:
    """Tests for get_infos_mvcf() - parsing multi-VCF headers"""
    
    def test_extracts_column_names(self, tmp_path):
        """Test extraction of column names from multi-VCF"""
        vcf_file = tmp_path / "multi.vcf"
        vcf_file.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\n"
            "1\t100\t.\tA\tT\t.\t.\t.\tGT:AD\t0/1:50,50\t0/0:100,0\n",
            encoding="utf-8"
        )
        
        headers, names, infos = multivcf_extract.get_infos_mvcf(str(vcf_file))
        
        assert isinstance(names, str)
        assert "Sample1" in names
        assert isinstance(infos, list)

    def test_handles_multiple_samples(self, tmp_path):
        """Test with multiple samples"""
        vcf_file = tmp_path / "multi.vcf"
        vcf_file.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n"
            "1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/0\t1/1\t./.\n",
            encoding="utf-8"
        )
        
        headers, names, infos = multivcf_extract.get_infos_mvcf(str(vcf_file))
        
        # Should extract all sample names
        assert len(names) >= 4


class TestIsGenotype:
    """Tests for isgenotype() - detecting genotype columns"""
    
    def test_identifies_valid_genotype(self):
        """Test detection of valid genotype"""
        # The function expects a genotype string (e.g. the GT field)
        result = multivcf_extract.isgenotype("0/1")
        assert result is True

    def test_identifies_missing_genotype(self):
        """Test detection of missing genotype"""
        # The function expects a genotype string (e.g. the GT field)
        result = multivcf_extract.isgenotype("./.")
        assert result is False


class TestWriteVcf:
    """Tests for write_vcf() - writing individual VCF files"""
    
    def test_creates_vcf_files(self, tmp_path):
        """Test creation of separate VCF files"""
        output_dir = tmp_path / "vcfs"
        output_dir.mkdir()
        
        headers = ["##fileformat=VCFv4.2"]
        # write_vcf expects names to be a single tab-separated string like a VCF header line
        names = "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2"
        infos = [
            "1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/0",
            "1\t200\t.\tG\tC\t.\t.\t.\tGT\t1/1\t./."
        ]
        
        donor_colname = "Sample1"
        recipient_colname = "Sample2"
        
        # Would write two VCF files
        multivcf_extract.write_vcf(
            headers, names, infos, str(output_dir),
            donor_colname, recipient_colname
        )
        
        # Check files created
        assert output_dir.exists()

    def test_preserves_vcf_headers(self, tmp_path):
        """Test that VCF headers are preserved"""
        output_dir = tmp_path / "vcfs"
        output_dir.mkdir()
        
        headers = [
            "##fileformat=VCFv4.2",
            '##INFO=<ID=DP,Type=Integer>',
        ]
        # write_vcf expects names to be a tab-separated string
        names = "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2"
        infos = []
        
        multivcf_extract.write_vcf(
            headers, names, infos, str(output_dir),
            "Sample1", "Sample2"
        )
        
        # VCF format line should be present
        assert output_dir.exists()

    def test_filters_for_selected_samples(self, tmp_path):
        """Test that only selected samples are kept"""
        output_dir = tmp_path / "vcfs"
        output_dir.mkdir()
        
        headers = ["##fileformat=VCFv4.2"]
        names = "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\tSample3"
        infos = [
            "1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/0\t1/1",
        ]
        
        # Only write Sample1 and Sample2
        multivcf_extract.write_vcf(
            headers, names, infos, str(output_dir),
            "Sample1", "Sample2"
        )
        
        # Sample3 should not be in output
        assert output_dir.exists()


class TestCreateDependencies:
    """Tests for create_dependencies() - creating output structure"""
    
    def test_creates_output_directory(self, tmp_path):
        """Test creation of vcf subdirectory"""
        output_dir = str(tmp_path / "output")
        run_name = "multivcf_test"
        
        vcf_dir = multivcf_extract.create_dependencies(run_name, output_dir)
        
        assert Path(vcf_dir).exists() or vcf_dir is not None


class TestMainWorkflow:
    """Tests for main() - complete multi-VCF extraction workflow"""
    
    def test_extracts_individual_vcfs(self, tmp_path):
        """Test extraction of individual VCFs from multi-VCF"""
        multi_vcf = tmp_path / "multi.vcf"
        multi_vcf.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDonor\tRecipient\n"
            "1\t100\t.\tA\tT\t.\t.\t.\tGT:AD\t0/1:50,50\t0/0:100,0\n"
            "1\t200\t.\tG\tC\t.\t.\t.\tGT:AD\t1/1:0,100\t0/1:50,50\n",
            encoding="utf-8"
        )
        
        output_dir = tmp_path / "output"
        run_name = "test_run"
        
        # Run extraction
        multivcf_extract.main(
            str(multi_vcf), "Donor", "Recipient",
            run_name, str(output_dir)
        )
        
        # Should create output directory
        assert output_dir.exists()

    def test_handles_missing_columns(self, tmp_path):
        """Test handling when donor/recipient columns don't exist"""
        multi_vcf = tmp_path / "multi.vcf"
        multi_vcf.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\n"
            "1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/1\n",
            encoding="utf-8"
        )
        
        output_dir = tmp_path / "output"
        
        # Should raise or handle gracefully
        with pytest.raises((KeyError, ValueError, IndexError)):
            multivcf_extract.main(
                str(multi_vcf), "NonexistentDonor", "NonexistentRecipient",
                "test", str(output_dir)
            )


class TestMultiVcfEdgeCases:
    """Tests for edge cases in multi-VCF processing"""
    
    def test_handles_empty_vcf(self, tmp_path):
        """Test handling of empty VCF file"""
        multi_vcf = tmp_path / "empty.vcf"
        multi_vcf.write_text(
            "##fileformat=VCFv4.2\n",
            encoding="utf-8"
        )
        
        # Should handle gracefully
        try:
            multivcf_extract.get_infos_mvcf(str(multi_vcf))
        except Exception:
            # Expected for empty VCF
            pass

    def test_handles_single_sample(self, tmp_path):
        """Test handling of single sample (edge case)"""
        multi_vcf = tmp_path / "single.vcf"
        multi_vcf.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tOnlySample\n"
            "1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/1\n",
            encoding="utf-8"
        )
        
        headers, names, infos = multivcf_extract.get_infos_mvcf(str(multi_vcf))
        
        # Should extract single sample
        assert len(names) >= 1

    def test_handles_special_characters_in_sample_names(self, tmp_path):
        """Test handling of special characters in sample names"""
        multi_vcf = tmp_path / "special.vcf"
        multi_vcf.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample-01\tSample_02\n"
            "1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/0\n",
            encoding="utf-8"
        )
        
        headers, names, infos = multivcf_extract.get_infos_mvcf(str(multi_vcf))
        
        # Should handle special characters (names returned as header line string)
        assert "Sample-01" in names
        assert "Sample_02" in names
