#coding:utf-8
"""
Comprehensive tests for cleavage.py module (NetChop processing)
"""
from types import SimpleNamespace
import gzip
import numpy as np
import pandas as pd
import pytest

from tools import cleavage


class TestGetVepIndicesFromVcf:
    """Tests for _get_vep_indices_from_vcf()"""
    
    def test_vcf_gz_file(self, tmp_path):
        """Test VEP indices extraction from gzipped VCF"""
        vcf_path = str(tmp_path / "test.vcf.gz")
        vcf_text = "\n".join([
            "##fileformat=VCFv4.2",
            '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Allele|Consequence|Gene|Feature|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|gnomADe_AF">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "1\t100\t.\tA\tT\t.\t.\tCSQ=T|missense|ENSG0001|ENST0001|1|2|3|K/N|aAa/aTa|0.001",
        ]) + "\n"
        with gzip.open(vcf_path, "wt", encoding="utf-8") as f:
            f.write(vcf_text)
        
        vep_indices = cleavage._get_vep_indices_from_vcf(vcf_path)
        assert vep_indices.gene == 2  # ENSG index
        assert vep_indices.transcript == 3  # ENST index
        assert vep_indices.prot == 6  # Protein position index

    def test_vcf_uncompressed_file(self, tmp_path):
        """Test VEP indices extraction from uncompressed VCF"""
        vcf_path = str(tmp_path / "test.vcf")
        vcf_text = "\n".join([
            "##fileformat=VCFv4.2",
            '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Allele|Consequence|Gene|Feature|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|gnomADe_AF">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "1\t100\t.\tA\tT\t.\t.\tCSQ=T|missense|ENSG0001|ENST0001|1|2|3|K/N|aAa/aTa|0.001",
        ]) + "\n"
        with open(vcf_path, "w", encoding="utf-8") as f:
            f.write(vcf_text)
        
        vep_indices = cleavage._get_vep_indices_from_vcf(vcf_path)
        assert vep_indices.gene == 2
        assert vep_indices.transcript == 3


class TestMmIntersect:
    """Tests for mm_intersect() - mismatches DataFrame intersection"""
    
    def test_merge_valid_dataframes(self):
        """Test correct merge of mismatches with transcripts"""
        mismatches_df = pd.DataFrame({
            "CHROM": [1, 1],
            "POS": [100, 200],
            "transcripts_x": ["ENST0001", "ENST0002"],
            "genes_x": ["ENSG0001", "ENSG0002"],
        })
        transcripts_pair = pd.DataFrame({
            "CHROM": ["1", "1"],
            "Transcript_id": ["ENST0001", "ENST0002"],
            "Gene_id": ["ENSG0001", "ENSG0002"],
            "peptide_ALT": ["MVKK", "LCCA"],
        })
        
        result = cleavage.mm_intersect(mismatches_df, transcripts_pair)
        assert len(result) > 0
        assert list(result.columns) == list(transcripts_pair.columns)

    def test_merge_empty_mismatches(self):
        """Test with empty mismatches DataFrame"""
        mismatches_df = pd.DataFrame({
            "CHROM": pd.Series([], dtype=int),
            "transcripts_x": pd.Series([], dtype=str),
            "genes_x": pd.Series([], dtype=str),
        })
        transcripts_pair = pd.DataFrame({
            "CHROM": ["1"],
            "Transcript_id": ["ENST0001"],
            "Gene_id": ["ENSG0001"],
        })
        
        result = cleavage.mm_intersect(mismatches_df, transcripts_pair)
        assert len(result) == 1

    def test_merge_no_matching_rows(self):
        """Test with no matching rows between DataFrames"""
        mismatches_df = pd.DataFrame({
            "CHROM": [1],
            "transcripts_x": ["ENST0001"],
            "genes_x": ["ENSG0001"],
        })
        transcripts_pair = pd.DataFrame({
            "CHROM": ["2"],
            "Transcript_id": ["ENST0002"],
            "Gene_id": ["ENSG0002"],
        })
        
        result = cleavage.mm_intersect(mismatches_df, transcripts_pair)
        assert len(result) == 2


class TestAddPepSeqChop:
    """Tests for add_pep_seq_chop() - adding peptide sequences"""
    
    def test_merge_peptides_success(self):
        """Test successful merge of transcripts with peptides"""
        transcripts_pair = pd.DataFrame({
            "CHROM": ["1"],
            "Gene_id": ["ENSG0001"],
            "Transcript_id": ["ENST0001"],
            "data": ["x"],
        })
        peptides_ensembl = pd.DataFrame({
            "CHROM": ["1"],
            "Gene_id": ["ENSG0001"],
            "Transcript_id": ["ENST0001"],
            "Peptide_id": ["ENSP0001"],
            "Sequence_aa": ["MVKKA"],
        })
        
        result = cleavage.add_pep_seq_chop(transcripts_pair, peptides_ensembl)
        assert len(result) > 0
        assert "Peptide_id" in result.columns
        assert "Sequence_aa" in result.columns

    def test_merge_no_matching_peptides(self):
        """Test when no peptides match transcripts"""
        transcripts_pair = pd.DataFrame({
            "CHROM": ["1"],
            "Gene_id": ["ENSG0001"],
            "Transcript_id": ["ENST0001"],
        })
        peptides_ensembl = pd.DataFrame({
            "CHROM": ["2"],
            "Gene_id": ["ENSG0002"],
            "Transcript_id": ["ENST0002"],
            "Peptide_id": ["ENSP0002"],
        })
        
        result = cleavage.add_pep_seq_chop(transcripts_pair, peptides_ensembl)
        assert len(result) == 0


class TestLoadPeptideIdMap:
    """Tests for load_peptide_id_map() - loading peptide ID mappings"""
    
    def test_load_mapping_success(self, tmp_path):
        """Test successful loading of peptide ID mapping"""
        tsv_file = tmp_path / "peptides.tsv"
        tsv_file.write_text(
            "Peptide_id\tSequence\n"
            "ENSP0000000001\tMVKKA\n"
            "ENSP0000000002\tLCCAV\n",
            encoding="utf-8"
        )
        
        mapping = cleavage.load_peptide_id_map(str(tsv_file))
        assert "0000000001" in mapping  # Last 10 chars
        assert mapping["0000000001"] == "ENSP0000000001"
        assert len(mapping) == 2

    def test_missing_peptide_id_column(self, tmp_path):
        """Test with missing Peptide_id column"""
        tsv_file = tmp_path / "peptides.tsv"
        tsv_file.write_text(
            "Sequence\tAnnotation\n"
            "MVKKA\tdata\n",
            encoding="utf-8"
        )
        
        with pytest.raises(ValueError, match="Peptide_id"):
            cleavage.load_peptide_id_map(str(tsv_file))

    def test_incomplete_rows(self, tmp_path):
        """Test with incomplete rows (fewer columns than header)"""
        tsv_file = tmp_path / "peptides.tsv"
        tsv_file.write_text(
            "Peptide_id\tSequence\tAnno\n"
            "ENSP0000000001\n"  # Missing columns
            "ENSP0000000002\tLCCAV\tdata\n",
            encoding="utf-8"
        )
        
        mapping = cleavage.load_peptide_id_map(str(tsv_file))
        # Should gracefully skip incomplete row
        assert len(mapping) == 2
        assert "0000000002" in mapping

    def test_duplicate_short_ids(self, tmp_path):
        """Test handling of duplicate short IDs (last 10 chars)"""
        tsv_file = tmp_path / "peptides.tsv"
        tsv_file.write_text(
            "Peptide_id\tSequence\n"
            "PREFIXABC0000000001\tMVKKA\n"
            "XYZABC0000000001\tLCCAV\n",  # Same last 10
            encoding="utf-8"
        )
        
        mapping = cleavage.load_peptide_id_map(str(tsv_file))
        # Last one wins
        assert mapping["0000000001"] == "XYZABC0000000001"


class TestParseNetchopOutput:
    """Tests for parse_netchop_output() - parsing NetChop results"""
    
    def test_min_run_length_filtering(self, tmp_path):
        """Test filtering by minimum run length"""
        netchop_file = tmp_path / "output.txt"
        netchop_file.write_text("\n".join([
            "header", "-----", "-----",
            "Pos AA C S Ident",
            "1 A . 0.1 ID1",
            "2 B . 0.1 ID1",
            "3 C . 0.1 ID1",
            "4 D S 0.1 ID1",
        ]), encoding="utf-8")
        
        runs, all_pos = cleavage.parse_netchop_output(str(netchop_file), min_run_length=3)
        assert runs["ID1"][0]["length"] == 3
        assert runs["ID1"][0]["positions"] == [1, 2, 3]

    def test_non_consecutive_fragments(self, tmp_path):
        """Test with non-consecutive positions (fragmented runs)"""
        netchop_file = tmp_path / "output.txt"
        netchop_file.write_text("\n".join([
            "header", "-----", "-----",
            "Pos AA C S Ident",
            "1 A . 0.1 ID1",
            "2 B . 0.1 ID1",
            "5 C . 0.1 ID1",  # Gap here
            "6 D . 0.1 ID1",
        ]), encoding="utf-8")
        
        runs, _ = cleavage.parse_netchop_output(str(netchop_file), min_run_length=2)
        assert len(runs["ID1"]) == 2  # Two separate runs
        assert runs["ID1"][0]["positions"] == [1, 2]
        assert runs["ID1"][1]["positions"] == [5, 6]

    def test_identifier_changes(self, tmp_path):
        """Test handling of identifier changes"""
        netchop_file = tmp_path / "output.txt"
        netchop_file.write_text("\n".join([
            "header", "-----", "-----",
            "Pos AA C S Ident",
            "1 A . 0.1 ID1",
            "2 B . 0.1 ID1",
            "1 C . 0.1 ID2",
            "2 D . 0.1 ID2",
        ]), encoding="utf-8")
        
        runs, _ = cleavage.parse_netchop_output(str(netchop_file), min_run_length=2)
        assert "ID1" in runs
        assert "ID2" in runs

    def test_malformed_lines(self, tmp_path):
        """Test handling of malformed lines"""
        netchop_file = tmp_path / "output.txt"
        netchop_file.write_text("\n".join([
            "header", "-----", "-----",
            "Pos AA C S Ident",
            "1 A . 0.1 ID1",
            "not_a_number A . 0.1 ID1",  # Malformed
            "2 B . 0.1 ID1",
        ]), encoding="utf-8")
        
        runs, _ = cleavage.parse_netchop_output(str(netchop_file), min_run_length=2)
        # Should handle gracefully, skipping malformed line
        assert "ID1" in runs

    def test_all_positions_captured(self, tmp_path):
        """Test that all positions are captured even if below min_run_length"""
        netchop_file = tmp_path / "output.txt"
        netchop_file.write_text("\n".join([
            "header", "-----", "-----",
            "Pos AA C S Ident",
            "1 A . 0.1 ID1",
            "2 B . 0.1 ID1",
            "3 C . 0.1 ID1",
        ]), encoding="utf-8")
        
        _, all_pos = cleavage.parse_netchop_output(str(netchop_file), min_run_length=5)
        assert all_pos["ID1"] == [(1, "A"), (2, "B"), (3, "C")]


class TestNetchopTablePrep:
    """Tests for netchop_table_prep() - preparing NetChop input table"""
    
    def test_table_creation_with_file_output(self, tmp_path):
        """Test successful table preparation and file output"""
        netchop_dir = tmp_path / "netchop"
        netchop_dir.mkdir()
        
        mismatches_df = pd.DataFrame({
            "CHROM": [1],
            "transcripts_x": ["ENST0001"],
            "genes_x": ["ENSG0001"],
            "peptide_ALT": ["MVKKA"],
        })
        transcripts_pair = pd.DataFrame({
            "CHROM": ["1"],
            "Transcript_id": ["ENST0001"],
            "Gene_id": ["ENSG0001"],
            "peptide_ALT": ["MVKKA"],
            "Peptide_id": ["ENSP0001"],
        })
        peptides_ensembl = pd.DataFrame({
            "CHROM": ["1"],
            "Gene_id": ["ENSG0001"],
            "Transcript_id": ["ENST0001"],
            "Peptide_id": ["ENSP0001"],
            "Sequence_aa_y": ["LCCA"],
        })
        
        args = SimpleNamespace(pair="", run_name="test")
        
        chop_table, path = cleavage.netchop_table_prep(
            mismatches_df, transcripts_pair, peptides_ensembl, args, str(netchop_dir)
        )
        
        assert len(chop_table) > 0
        assert path.endswith("_netchop_table.csv")
        
    def test_remove_duplicates(self, tmp_path):
        """Test that duplicates are removed"""
        netchop_dir = tmp_path / "netchop"
        netchop_dir.mkdir()
        
        mismatches_df = pd.DataFrame({
            "CHROM": ["1", "1"],
            "transcripts_x": ["ENST0001", "ENST0001"],
            "genes_x": ["ENSG0001", "ENSG0001"],
            "peptide_ALT": ["MVKKA", "MVKKA"],
            "Peptide_id": ["ENSP0001", "ENSP0001"],
        })
        transcripts_pair = pd.DataFrame({
            "CHROM": ["1", "1"],
            "Transcript_id": ["ENST0001", "ENST0001"],
            "Gene_id": ["ENSG0001", "ENSG0001"],
            "peptide_ALT": [np.nan, "MVKKA"],
            "Peptide_id": ["ENSP0001", "ENSP0001"],
            "Sequence_aa": [np.nan, np.nan],
        })
        peptides_ensembl = pd.DataFrame({
            "CHROM": ["1"],
            "Gene_id": ["ENSG0001"],
            "Transcript_id": ["ENST0001"],
            "Peptide_id": ["ENSP0001"],
            "Sequence_aa": ["MVKKA"],
        })
        
        args = SimpleNamespace(pair="", run_name="test")
        
        chop_table, _ = cleavage.netchop_table_prep(
            mismatches_df, transcripts_pair, peptides_ensembl, args, str(netchop_dir)
        )
        
        # Should have removed duplicates
        assert len(chop_table) == 1

    def test_remove_nan_peptides(self, tmp_path):
        """Test that peptides with NaN are removed"""
        netchop_dir = tmp_path / "netchop"
        netchop_dir.mkdir()
        
        mismatches_df = pd.DataFrame({
            "CHROM": [1, 1],
            "transcripts_x": ["ENST0001", "ENST0002"],
            "genes_x": ["ENSG0001", "ENSG0002"],
            "peptide_ALT": ["MVKKA", None],
            "Peptide_id": ["ENSP0001", "ENSP0002"],
        })
        transcripts_pair = pd.DataFrame({
            "CHROM": ["1", "1"],
            "Transcript_id": ["ENST0001", "ENST0002"],
            "Gene_id": ["ENSG0001", "ENSG0002"],
            "peptide_ALT": ["MVKKA", None],
            "Peptide_id": ["ENSP0001", "ENSP0002"],
            "Sequence_aa": [np.nan, np.nan],
        })
        peptides_ensembl = pd.DataFrame({
            "CHROM": ["1", "1"],
            "Gene_id": ["ENSG0001", "ENSG0002"],
            "Transcript_id": ["ENST0001", "ENST0002"],
            "Peptide_id": ["ENSP0001", "ENSP0002"],
            "Sequence_aa": ["MVKKA", "LCCA"],
        })
        
        args = SimpleNamespace(pair="", run_name="test")
        
        chop_table, _ = cleavage.netchop_table_prep(
            mismatches_df, transcripts_pair, peptides_ensembl, args, str(netchop_dir)
        )
        
        # NaN peptides should be removed
        assert pd.isna(chop_table["peptide_ALT"]).sum() == 0


class TestPostprocessNetchop:
    """Tests for postprocess_netchop() - post-processing NetChop results"""
    
    def test_peptide_output_creation(self, tmp_path):
        """Test creation of peptides output file"""
        netchop_dir = tmp_path / "netchop"
        netchop_dir.mkdir()
        
        # Create NetChop output
        netchop_output = netchop_dir / "netchop_out.txt"
        netchop_output.write_text("\n".join([
            "header", "-----", "-----",
            "Pos AA C S Ident",
            "1 M . 0.1 0000356701",
            "2 V . 0.1 0000356701",
            "3 K . 0.1 0000356701",
        ]), encoding="utf-8")
        
        # Create peptide ID map
        chop_table_file = netchop_dir / "chop_table.csv"
        chop_table_file.write_text(
            "Peptide_id\tSequence\n"
            "ENSP0000000356701\tMVKKA\n",
            encoding="utf-8"
        )
        
        args = SimpleNamespace(pair="", run_name="test")
        
        cleavage.postprocess_netchop(str(netchop_output), str(chop_table_file), args, str(netchop_dir))
        
        # Check output file created
        peptides_file = netchop_dir / "test_netchop_peptides.txt"
        assert peptides_file.exists()

    def test_missing_id_fallback(self, tmp_path):
        """Test fallback when peptide ID not found in map"""
        netchop_dir = tmp_path / "netchop"
        netchop_dir.mkdir()
        
        netchop_output = netchop_dir / "netchop_out.txt"
        netchop_output.write_text("\n".join([
            "header", "-----", "-----",
            "Pos AA C S Ident",
            "1 M . 0.1 UNKNOWN0001",  # ID not in map
        ]), encoding="utf-8")
        
        chop_table_file = netchop_dir / "chop_table.csv"
        chop_table_file.write_text(
            "Peptide_id\tSequence\n",
            encoding="utf-8"
        )
        
        args = SimpleNamespace(pair="", run_name="test")
        
        # Should not raise error, use short_id as fallback
        cleavage.postprocess_netchop(str(netchop_output), str(chop_table_file), args, str(netchop_dir))
        
        peptides_file = netchop_dir / "test_netchop_peptides.txt"
        assert peptides_file.exists()
