#coding:utf-8
"""
Integration tests for Python stages invoked by Nextflow (no external tools)
"""
from types import SimpleNamespace
from pathlib import Path
import gzip
import numpy as np
import pandas as pd
import pytest

from tools import (
    ams_helpers, aams_helpers, parsing_functions,
    cleavage, table_operations
)


class TestAmsIntegrationPipeline:
    """Integration tests for the Allo-Count Python stage."""
    
    def test_ams_workflow_creates_output_structure(self, tmp_path, nextflow_environment):
        """Test that AMS creates complete output directory structure"""
        output_dir = str(tmp_path / "output")
        run_name = "ams_test"
        
        # Create run directory
        run_path, run_tables, run_plots, run_ams, run_logs = ams_helpers.create_run_directory(
            run_name, output_dir
        )
        
        # Write log file
        args = SimpleNamespace(
            donor="/path/donor.vcf",
            recipient="/path/recipient.vcf",
            orientation="dr",
            imputation="imputation",
            min_dp=10, max_dp=1000, min_ad=2, min_gq=20,
            homozygosity_thr=0.8, base_length=3,
            frameshift=False,
            pair="", run_name=run_name, output_dir=output_dir,
        )
        ams_helpers.write_log(run_logs, args)
        
        # Verify structure
        assert Path(run_tables).exists()
        assert Path(run_logs).exists()
        assert (Path(run_logs) / "run.log").exists()

    def test_ams_handles_mismatch_detection(self, tmp_path):
        """Test mismatch detection workflow"""
        # Create sample DataFrames
        df_donor = pd.DataFrame({
            "CHROM": [1, 1],
            "POS": [100, 200],
            "REF": ["A", "G"],
            "ALT": ["T", "C"],
            "GT": ["1/1", "1/1"],
            "aa_REF": ["A", "G"],
            "aa_ALT": ["T", "C"],
        })
        df_recipient = pd.DataFrame({
            "CHROM": [1, 1],
            "POS": [100, 200],
            "REF": ["A", "G"],
            "ALT": ["A", "C"],
            "GT": ["0/0", "1/1"],
            "aa_REF": ["A", "G"],
            "aa_ALT": ["A", "C"],
        })
        
        # Merge
        merged, _, _ = ams_helpers.merge_dfs(df_donor, df_recipient, "dr", "imputation")
        
        # Keep ALT
        filtered = ams_helpers.keep_alt(merged, "x", "y")
        
        # Count mismatches
        count_df, mismatch_count = ams_helpers.count_mismatches(filtered, "dr")
        
        assert mismatch_count == 1
        assert count_df["POS"].tolist() == [100]

    def test_ams_saves_output_files(self, tmp_path, nextflow_environment):
        """Test that AMS saves required output files"""
        output_dir = tmp_path / "output"
        run_name = "test"
        _, _, _, run_ams, _ = ams_helpers.create_run_directory(
            run_name, str(output_dir)
        )
        
        mismatch_count = 5
        
        args = SimpleNamespace(
            donor="donor.vcf", recipient="recipient.vcf",
            pair="P01", run_name=run_name, orientation="dr",
            min_dp=10, max_dp=1000, min_ad=2, min_gq=20,
            homozygosity_thr=0.8, base_length=3, output_dir=str(output_dir),
        )
        
        # Create dummy input files
        df_donor_file = str(output_dir / "donor.tsv")
        df_recipient_file = str(output_dir / "recipient.tsv")
        mismatches_file = str(output_dir / "mismatches.tsv")
        Path(df_donor_file).write_text("data\n", encoding="utf-8")
        Path(df_recipient_file).write_text("data\n", encoding="utf-8")
        Path(mismatches_file).write_text("data\n", encoding="utf-8")
        
        result_path = table_operations.save_mismatch(
            str(run_ams), args, mismatch_count,
            df_donor_file, df_recipient_file, mismatches_file
        )
        
        table_operations.create_AMS_df(result_path)
        result = pd.read_csv(Path(result_path) / "AMS_df.tsv", sep="\t")
        assert result["pair"].tolist() == ["P01"]
        assert result["ams"].tolist() == [5]
        published_dir, _ = nextflow_environment
        assert result["donor_table"].tolist() == [str(published_dir / "donor.tsv")]
        assert result["recipient_table"].tolist() == [str(published_dir / "recipient.tsv")]
        assert result["mismatches_table"].tolist() == [str(published_dir / "mismatches.tsv")]


class TestAamsIntegrationPipeline:
    """Integration tests for the Allo-Affinity Python stage."""
    
    def test_aams_creates_dependencies(self, tmp_path):
        """Test that AAMS creates all necessary directories"""
        output_dir = str(tmp_path / "output")
        run_name = "aams_test"
        
        aams_run_tables, netmhc_dir, aams_path, netchop_dir = aams_helpers.create_aams_dependencies(
            run_name, output_dir
        )
        
        assert Path(aams_run_tables).exists()
        assert Path(netmhc_dir).exists()
        assert Path(aams_path).exists()
        assert Path(netchop_dir).exists()

    def test_aams_netchop_workflow(self, tmp_path):
        """Test NetChop portion of AAMS pipeline"""
        netchop_dir = tmp_path / "netchop"
        netchop_dir.mkdir()
        
        # Create mock data
        mismatches_df = pd.DataFrame({
            "CHROM": ["1"],
            "POS": ["100"],
            "transcripts_x": ["ENST0001"],
            "genes_x": ["ENSG0001"],
        })
        transcripts_pair = pd.DataFrame({
            "CHROM": ["1"],
            "POS": ["100"],
            "Transcript_id": ["ENST0001"],
            "Gene_id": ["ENSG0001"],
            "peptide_ALT": ["MVKKA"],
            "Peptide_id": ["ENSP0001"],
            "Sequence_aa": [np.nan],
            "aa_REF": ["K"],
            "aa_alt_indiv": ["N"],
        })
        peptides_ensembl = pd.DataFrame({
            "CHROM": ["1"],
            "Gene_id": ["ENSG0001"],
            "Transcript_id": ["ENST0001"],
            "Peptide_id": ["ENSP0001"],
            "Sequence_aa": ["MVKKA"],
        })
        
        args = SimpleNamespace(pair="", run_name="test")
        
        # Run netchop table prep
        chop_table, chop_path = cleavage.netchop_table_prep(
            mismatches_df, transcripts_pair, peptides_ensembl, args, str(netchop_dir)
        )
        
        assert chop_table["Peptide_id"].tolist() == ["ENSP0001"]
        assert chop_table["peptide_ALT"].tolist() == ["MVKKA"]
        assert "aa_REF" not in chop_table.columns
        saved = pd.read_csv(chop_path, sep="\t", dtype={"CHROM": str, "POS": str})
        pd.testing.assert_frame_equal(saved, chop_table, check_dtype=False)

    def test_aams_peptide_generation_workflow(self, tmp_path):
        """Test peptide generation in AAMS"""
        output_dir = tmp_path / "output"
        run_name = "peptide_test"
        
        aams_run_tables, _, _, _ = aams_helpers.create_aams_dependencies(run_name, str(output_dir))
        
        # Create sample transcripts
        transcripts_pair = pd.DataFrame({
            "Peptide_id": ["ENSP0001", "ENSP0002"],
            "Sequence_aa": ["MVKKA", "LCCA"],
            "Gene_id": ["ENSG0001", "ENSG0002"],
            "Transcript_id": ["ENST0001", "ENST0002"],
            "CHROM": ["1", "1"],
            "POS": [100, 200],
            "peptide_REF": ["MVKKA", "LCCA"],
            "peptide": ["MVKKA", "LCCA"],
            "hla_peptides": ["MVKKA", "LCCA"],
        })
        
        fasta_path = str(Path(aams_run_tables) / "peptides.fasta")
        
        aams_helpers.write_pep_fasta(fasta_path, transcripts_pair)
        
        assert Path(fasta_path).exists()
        content = Path(fasta_path).read_text()
        assert "ENSP0001" in content
        assert "MVKKA" in content


class TestAmsStageIntegration:
    """VCF parsing and metadata hand-off between Python stages."""
    
    def test_parses_staged_vcf(self, tmp_path):
        """Test VEP parsing of a compressed worker input."""
        # Create minimal VCF
        vcf_file = tmp_path / "test.vcf.gz"
        vcf_text = "\n".join([
            "##fileformat=VCFv4.2",
            '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Allele|Consequence|Gene|Feature|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|gnomADe_AF">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "1\t100\t.\tA\tT\t.\t.\tCSQ=T|missense|ENSG0001|ENST0001|1|2|3|K/N|aAa/aTa|0.001",
        ]) + "\n"
        
        with gzip.open(str(vcf_file), "wt", encoding="utf-8") as f:
            f.write(vcf_text)
        
        # Parse VEP
        df_infos, vep_indices = parsing_functions.gzvcf_vep_parser(str(vcf_file), frameshift_mode=False)
        
        assert df_infos["#CHROM"].tolist() == ["1"]
        assert df_infos["POS"].tolist() == ["100"]
        assert vep_indices.gene == 2
        assert vep_indices.transcript == 3
        assert vep_indices.frameshift is None

    def test_log_metadata_is_readable_by_aams(self, tmp_path, nextflow_environment):
        """The AAMS stage reads provenance written by the AMS stage."""
        output_dir = tmp_path / "output"
        run_name = "complete_test"
        
        # Step 1: Create directories
        run_path, run_tables, _, _, run_logs = ams_helpers.create_run_directory(
            run_name, str(output_dir)
        )
        
        # Step 2: Write log
        args = SimpleNamespace(
            donor="donor.vcf", recipient="recipient.vcf",
            orientation="dr", imputation="imputation",
            min_dp=10, max_dp=1000, min_ad=2, min_gq=20,
            homozygosity_thr=0.8, base_length=3,
            frameshift=False,
            pair="", run_name=run_name, output_dir=str(output_dir),
        )
        ams_helpers.write_log(run_logs, args)
        
        # Step 3: Check structure
        assert Path(run_logs).exists()
        assert (Path(run_logs) / "run.log").exists()
        
        log_content = (Path(run_logs) / "run.log").read_text()
        assert aams_helpers.read_log_field(Path(run_logs) / "run.log", "Orientation") == "dr"
        published_dir, command = nextflow_environment
        assert aams_helpers.read_log_field(Path(run_logs) / "run.log", "Output_dir") == str(published_dir)
        assert f"Nextflow_command: {command}\n" in log_content


class TestAamsNetMhcWorkflow:
    """Test AAMS with NetMHCpan integration"""
    
    def test_prepares_netmhc_fasta_input(self, tmp_path):
        """Test preparation of NetMHCpan FASTA input"""
        aams_dir = tmp_path / "aams"
        aams_dir.mkdir()
        
        # Create peptide table
        peptides = pd.DataFrame({
            "Peptide_id": ["ENSP0001", "ENSP0002"],
            "Sequence_aa": ["MVKKA", "LCCA"],
            "Gene_id": ["ENSG0001", "ENSG0002"],
            "Transcript_id": ["ENST0001", "ENST0002"],
            "CHROM": ["1", "1"],
            "POS": [100, 200],
            "peptide_REF": ["MVKKA", "LCCA"],
            "peptide": ["MVKKA", "LCCA"],
            "hla_peptides": ["MVKKA", "LCCA"],
        })
        
        fasta_file = aams_dir / "netmhc_input.fasta"
        aams_helpers.write_pep_fasta(str(fasta_file), peptides)
        
        # Verify FASTA format
        assert fasta_file.exists()
        lines = fasta_file.read_text().strip().split('\n')
        assert lines[0].startswith('>')
        assert all(line.isalpha() or line.startswith('>') for line in lines)


class TestPipelineErrorHandling:
    """Tests for error handling across pipelines"""
    
    def test_handles_missing_input_file(self, tmp_path):
        """Test graceful handling of missing input files"""
        nonexistent = str(tmp_path / "missing.vcf")
        
        with pytest.raises(FileNotFoundError) as error:
            parsing_functions.vcf_vep_parser(nonexistent, frameshift_mode=False)
        assert error.value.filename == nonexistent

    def test_handles_malformed_vcf(self, tmp_path):
        """Test handling of malformed VCF"""
        vcf_file = tmp_path / "malformed.vcf"
        vcf_file.write_text("not a valid vcf file\n", encoding="utf-8")
        
        with pytest.raises(ValueError, match="does not contain the VEP information"):
            parsing_functions.vcf_vep_parser(str(vcf_file), frameshift_mode=False)

    def test_handles_empty_dataframes(self):
        """Test handling of empty DataFrames"""
        empty_df = pd.DataFrame({
            "GT_x": [],
            "GT_y": [],
        })
        result = ams_helpers.keep_alt(empty_df, "x", "y")
        pd.testing.assert_frame_equal(result, empty_df)


class TestPipelinePermissions:
    """Tests for file permissions and access"""
    
    def test_creates_writable_directories(self, tmp_path):
        """Test that created directories are writable"""
        output_dir = tmp_path / "output"
        run_name = "writetest"
        
        run_path, run_tables, _, _, _ = ams_helpers.create_run_directory(
            run_name, str(output_dir)
        )
        
        # Should be able to write test file
        test_file = Path(run_tables) / "test.txt"
        test_file.write_text("test", encoding="utf-8")
        assert test_file.exists()
