import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import pandas as pd


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from tools import aams_helpers


class CleavageSourceSelectionTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.output_dir = Path(self.temp_dir.name)
        self.run_name = "test-run"
        self.run_tables_dir = (
            self.output_dir / "runs" / self.run_name / "run_tables"
        )
        self.run_tables_dir.mkdir(parents=True)

        columns = [
            "CHROM",
            "POS",
            "transcripts",
            "genes",
            "aa_REF",
            "aa_ALT",
            "aa_ref_indiv",
            "aa_alt_indiv",
            "aa_indiv",
            "GT",
        ]
        pd.DataFrame(
            [["1", "10", "T1", "G1", "B", "Z", "B", "Z", "B,Z", "0/1"]],
            columns=columns,
        ).to_csv(
            self.run_tables_dir / "PAIR_D0_donor.tsv",
            sep="\t",
            index=False,
        )
        pd.DataFrame(
            [["1", "20", "T1", "G1", "B", "Y", None, "Y", "Y", "1/1"]],
            columns=columns,
        ).to_csv(
            self.run_tables_dir / "PAIR_R0_recipient.tsv",
            sep="\t",
            index=False,
        )

        self.refseq_path = self.output_dir / "refseq.tsv"
        pd.DataFrame({"transcript_stable_id": ["T1"]}).to_csv(
            self.refseq_path,
            sep="\t",
            index=False,
        )
        self.mismatches_path = self.output_dir / "mismatches.tsv"
        pd.DataFrame({"CHROM": ["1"], "POS": ["10"]}).to_csv(
            self.mismatches_path,
            sep="\t",
            index=False,
        )
        self.peptides_ensembl = pd.DataFrame(
            {
                "Gene_id": ["G1"],
                "Transcript_id": ["T1"],
                "CHROM": ["1"],
                "Peptide_id": ["P1"],
                "Sequence_aa": ["ABCD"],
            }
        )
        self.args = SimpleNamespace(
            output_dir=str(self.output_dir),
            run_name=self.run_name,
            length=9,
            pair="PAIR",
        )

    def tearDown(self):
        self.temp_dir.cleanup()

    def build_for(self, individual, pos):
        individual_vep_df = pd.DataFrame(
            {
                "CHROM": ["1"],
                "POS": [pos],
                "Transcript_id": ["T1"],
                "Protein_position": ["2"],
            }
        )
        _, individual_proteins_df, _ = aams_helpers.build_peptides(
            args=self.args,
            mismatches_path=self.mismatches_path,
            individual_vep_df=individual_vep_df,
            cleavage_mode=True,
            ens_transcripts={},
            peptides_ensembl=self.peptides_ensembl.copy(),
            refseq_file=self.refseq_path,
            individual=individual,
        )
        return individual_proteins_df

    def test_donor_uses_d0_and_recipient_uses_r0(self):
        donor_proteins = self.build_for("donor", "10")
        recipient_proteins = self.build_for("recipient", "20")

        self.assertEqual(donor_proteins.iloc[0]["peptide_ALT"], "AZCD")
        self.assertEqual(recipient_proteins.iloc[0]["peptide_ALT"], "AYCD")

    def test_cleavage_requires_an_explicit_individual(self):
        with self.assertRaisesRegex(ValueError, "requires individual"):
            self.build_for(None, "10")

    def test_reference_genotype_does_not_apply_vep_alt(self):
        donor_path = self.run_tables_dir / "PAIR_D0_donor.tsv"
        donor_df = pd.read_csv(donor_path, sep="\t")
        donor_df.loc[0, ["GT", "aa_ref_indiv", "aa_alt_indiv", "aa_indiv"]] = [
            "0/0",
            "B",
            None,
            "B",
        ]
        donor_df.to_csv(donor_path, sep="\t", index=False)

        donor_proteins = self.build_for("donor", "10")

        self.assertTrue(donor_proteins.empty)

    def test_only_individual_alt_is_applied_in_mixed_genotypes(self):
        donor_path = self.run_tables_dir / "PAIR_D0_donor.tsv"
        donor_df = pd.read_csv(donor_path, sep="\t")
        reference_row = donor_df.iloc[0].copy()
        reference_row["POS"] = "11"
        reference_row["GT"] = "0/0"
        reference_row["aa_REF"] = "C"
        reference_row["aa_ALT"] = "X"
        reference_row["aa_ref_indiv"] = "C"
        reference_row["aa_alt_indiv"] = None
        reference_row["aa_indiv"] = "C"
        donor_df = pd.concat(
            [donor_df, reference_row.to_frame().T],
            ignore_index=True,
        )
        donor_df.to_csv(donor_path, sep="\t", index=False)
        individual_vep_df = pd.DataFrame(
            {
                "CHROM": ["1", "1"],
                "POS": ["10", "11"],
                "Transcript_id": ["T1", "T1"],
                "Protein_position": ["2", "3"],
            }
        )

        _, donor_proteins, _ = aams_helpers.build_peptides(
            args=self.args,
            mismatches_path=self.mismatches_path,
            individual_vep_df=individual_vep_df,
            cleavage_mode=True,
            ens_transcripts={},
            peptides_ensembl=self.peptides_ensembl.copy(),
            refseq_file=self.refseq_path,
            individual="donor",
        )

        self.assertEqual(donor_proteins.iloc[0]["peptide_ALT"], "AZCD")
        self.assertNotIn("X", donor_proteins.iloc[0]["peptide_ALT"])

    def test_ambiguous_run_table_is_rejected(self):
        duplicate = self.run_tables_dir / "OTHER_D0_donor.tsv"
        duplicate.touch()

        with self.assertRaisesRegex(RuntimeError, "Expected one run table"):
            aams_helpers.find_single_run_table(self.run_tables_dir, "_D0_")


if __name__ == "__main__":
    unittest.main()
