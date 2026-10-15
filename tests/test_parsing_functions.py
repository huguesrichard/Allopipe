#coding:utf-8
"""
Tests for parsing_functions.py module (VCF & VEP parsing)
"""
import gzip
import pandas as pd
import pytest

from tools import parsing_functions


class TestVepIndicesNamedTuple:
    """Tests for VepIndices namedtuple structure"""
    
    def test_contains_required_fields(self):
        """Test that VepIndices has all required fields"""
        # Create a VepIndices with all fields
        indices = parsing_functions.VepIndices(
            consequence=1, gene=2, transcript=3, cdna=4, cds=5,
            prot=6, aa=7, codons=8, gnomad=9, frameshift=None
        )
        
        assert indices.gene == 2
        assert indices.transcript == 3
        assert indices.prot == 6
        assert indices.consequence == 1
        assert indices.gnomad == 9
        assert indices.frameshift is None


class TestVcfVepParser:
    """Tests for vcf_vep_parser() - uncompressed VCF parsing"""
    
    @pytest.mark.parametrize("frameshift_mode", [False, True])
    def test_parses_valid_vcf(self, tmp_path, frameshift_mode):
        """Test parsing of valid uncompressed VCF"""
        vcf_file = tmp_path / "test.vcf"
        vcf_text = "\n".join([
            "##fileformat=VCFv4.2",
            '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Allele|Consequence|Gene|Feature|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|gnomADe_AF">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "1\t100\t.\tA\tT\t.\t.\tCSQ=T|missense_variant|ENSG0001|ENST0001|123|45|42|K/N|aAa/aTa|0.001",
        ]) + "\n"
        if frameshift_mode:
            vcf_text = vcf_text.replace('gnomADe_AF">', 'gnomADe_AF|FrameshiftSequence">')
            vcf_text = vcf_text.replace("|0.001\n", "|0.001|MVKKA\n")
        vcf_file.write_text(vcf_text, encoding="utf-8")
        
        df_infos, vep_indices = parsing_functions.vcf_vep_parser(str(vcf_file), frameshift_mode=frameshift_mode)
        
        assert isinstance(df_infos, pd.DataFrame)
        assert vep_indices.gene == 2
        assert vep_indices.transcript == 3
        assert vep_indices.frameshift == (10 if frameshift_mode else None)
        assert df_infos["#CHROM"].tolist() == ["1"]
        assert df_infos["POS"].tolist() == ["100"]
        assert df_infos["INFO"].iloc[0].endswith("|MVKKA" if frameshift_mode else "|0.001")

    def test_handles_multiple_variants(self, tmp_path):
        """Test parsing with multiple variants"""
        vcf_file = tmp_path / "test.vcf"
        vcf_text = "\n".join([
            "##fileformat=VCFv4.2",
            '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Allele|Consequence|Gene|Feature|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|gnomADe_AF">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "1\t100\t.\tA\tT\t.\t.\tCSQ=T|missense|ENSG0001|ENST0001|1|2|3|K/N|aAa/aTa|0.001",
            "1\t200\t.\tG\tC\t.\t.\tCSQ=C|frameshift|ENSG0002|ENST0002|1|2|3|G/R|gGg/cCc|0.002",
        ]) + "\n"
        vcf_file.write_text(vcf_text, encoding="utf-8")
        
        df_infos, _ = parsing_functions.vcf_vep_parser(str(vcf_file), frameshift_mode=False)
        
        assert df_infos["POS"].tolist() == ["100", "200"]
        assert df_infos["REF"].tolist() == ["A", "G"]
        assert df_infos["ALT"].tolist() == ["T", "C"]
        assert df_infos["INFO"].tolist() == [
            "CSQ=T|missense|ENSG0001|ENST0001|1|2|3|K/N|aAa/aTa|0.001",
            "CSQ=C|frameshift|ENSG0002|ENST0002|1|2|3|G/R|gGg/cCc|0.002",
        ]

    def test_handles_missing_vep_info(self, tmp_path):
        """Test with missing VEP INFO annotation"""
        vcf_file = tmp_path / "test.vcf"
        vcf_text = "\n".join([
            "##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "1\t100\t.\tA\tT\t.\t.\t.",  # No CSQ info
        ]) + "\n"
        vcf_file.write_text(vcf_text, encoding="utf-8")
        
        with pytest.raises(ValueError, match="does not contain the VEP information"):
            parsing_functions.vcf_vep_parser(str(vcf_file), frameshift_mode=False)


class TestGzvcfVepParser:
    """Tests for gzvcf_vep_parser() - compressed VCF parsing"""
    
    def test_parses_valid_gzipped_vcf(self, tmp_path):
        """Test parsing of valid gzipped VCF"""
        vcf_file = tmp_path / "test.vcf.gz"
        vcf_text = "\n".join([
            "##fileformat=VCFv4.2",
            '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Allele|Consequence|Gene|Feature|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|gnomADe_AF">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "1\t100\t.\tA\tT\t.\t.\tCSQ=T|missense_variant|ENSG0001|ENST0001|123|45|42|K/N|aAa/aTa|0.001",
        ]) + "\n"
        
        with gzip.open(str(vcf_file), "wt", encoding="utf-8") as f:
            f.write(vcf_text)
        
        df_infos, vep_indices = parsing_functions.gzvcf_vep_parser(str(vcf_file), frameshift_mode=False)
        
        assert isinstance(df_infos, pd.DataFrame)
        assert vep_indices.gene == 2

    def test_gzip_same_result_as_uncompressed(self, tmp_path):
        """Test that gzipped and uncompressed VCFs give same results"""
        vcf_text = "\n".join([
            "##fileformat=VCFv4.2",
            '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Allele|Consequence|Gene|Feature|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|gnomADe_AF">',
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "1\t100\t.\tA\tT\t.\t.\tCSQ=T|missense_variant|ENSG0001|ENST0001|1|2|3|K/N|aAa/aTa|0.001",
        ]) + "\n"
        
        # Create uncompressed
        vcf_uncompressed = tmp_path / "test.vcf"
        vcf_uncompressed.write_text(vcf_text, encoding="utf-8")
        
        # Create compressed
        vcf_compressed = tmp_path / "test.vcf.gz"
        with gzip.open(str(vcf_compressed), "wt", encoding="utf-8") as f:
            f.write(vcf_text)
        
        df1, idx1 = parsing_functions.vcf_vep_parser(str(vcf_uncompressed), frameshift_mode=False)
        df2, idx2 = parsing_functions.gzvcf_vep_parser(str(vcf_compressed), frameshift_mode=False)
        
        pd.testing.assert_frame_equal(df1, df2)
        assert vars(idx1) == vars(idx2)


class TestExtractAaFromVep:
    """Tests for extract_aa_from_vep()"""
    
    def test_extracts_genes_transcripts_and_amino_acids(self):
        """VEP extraction excludes synonymous variants and preserves amino acids."""
        df_infos = pd.DataFrame({
            "CHROM": ["1", "1", "1"],
            "POS": [100, 200, 300],
            "INFO": [
                "T|missense_variant|ENSG0001|ENST0001|123|45|42|K/N|aAa/aTa|0.001",
                "C|frameshift_variant|ENSG0002|ENST0002|1|2|3|G/R|gGg/cCc|0.002",
                "A|synonymous_variant|ENSG0003|ENST0003|1|2|3|K/K|aAa/aAa|0.003",
            ]
        })
        vep_indices = parsing_functions.VepIndices(
            consequence=1, gene=2, transcript=3, cdna=4, cds=5,
            prot=6, aa=7, codons=8, gnomad=9, frameshift=None
        )
        
        result = parsing_functions.extract_aa_from_vep(df_infos, vep_indices)
        
        assert result["genes"].tolist() == ["ENSG0001", "ENSG0002"]
        assert result["transcripts"].tolist() == ["ENST0001", "ENST0002"]
        assert result["aa_REF"].tolist() == ["K", "G"]
        assert result["aa_ALT"].tolist() == ["N", "R"]
        assert result["Frameshift_sequence"].tolist() == ["", ""]
        assert result["POS"].tolist() == [100, 200]


class TestReadFasta:
    """Tests for read_fasta()"""
    
    def test_reads_valid_fasta(self, tmp_path):
        """Test reading valid FASTA file"""
        fasta_file = tmp_path / "test.fa"
        fasta_file.write_text(
            ">seq1 description\n"
            "MVKKAMVKKA\n"
            "MVKKV\n"
            ">seq2 description\n"
            "LCCA\n",
            encoding="utf-8"
        )
        
        result = parsing_functions.read_fasta(str(fasta_file))
        
        assert result == {"seq1": "MVKKAMVKKAMVKKV", "seq2": "LCCA"}

    def test_handles_empty_fasta(self, tmp_path):
        """Test with empty FASTA file"""
        fasta_file = tmp_path / "empty.fa"
        fasta_file.write_text("", encoding="utf-8")
        
        # Empty FASTA files cause UnboundLocalError in read_fasta
        with pytest.raises(UnboundLocalError):
            parsing_functions.read_fasta(str(fasta_file))

    def test_multiline_sequences(self, tmp_path):
        """Test FASTA with multi-line sequences"""
        fasta_file = tmp_path / "multiline.fa"
        fasta_file.write_text(
            ">seq1 description\n"
            "MVKKAMVKKAMVKKA\n"
            "LCCA\n"
            ">seq2 description\n"
            "AAAA\n",
            encoding="utf-8"
        )
        
        result = parsing_functions.read_fasta(str(fasta_file))
        
        assert result == {"seq1": "MVKKAMVKKAMVKKALCCA", "seq2": "AAAA"}


class TestReadPepFa:
    """Tests for read_pep_fa()"""
    
    def test_reads_peptide_fasta(self, tmp_path):
        """Test reading peptide FASTA (no description)"""
        pep_file = tmp_path / "peptides.fa"
        # read_pep_fa expects protein_coding in header with specific format
        pep_file.write_text(
            ">ENSP0001.1 protein_coding GRCh38:1:1000:2000:-1 gene:GENE1.1 transcript:ENST0001.1 gene_biotype:protein_coding\n"
            "MVKKA\n"
            ">ENSP0002.1 protein_coding GRCh38:1:3000:4000:1 gene:GENE2.1 transcript:ENST0002.1 gene_biotype:protein_coding\n"
            "LCCA\n",
            encoding="utf-8"
        )
        
        result = parsing_functions.read_pep_fa(str(pep_file))
        
        assert isinstance(result, dict)
        # Keys are peptide IDs (without version)
        assert "ENSP0001" in result
        # Values are lists: [gene, coords, transcript, sequence]
        assert result["ENSP0001"][3] == "MVKKA"
