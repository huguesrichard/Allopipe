"""Peptide selection after NetChop, targeting the Nextflow implementation."""
import csv
from pathlib import Path
from types import SimpleNamespace

import pytest

from tools import cleavage


COLUMNS = ["CHROM", "POS", "Peptide_id", "Peptide", "hla_peptides"]


def peptide(protein, kmer, chrom="1", pos="100", sequence=None):
    return dict(zip(COLUMNS, [chrom, pos, protein, sequence or kmer, kmer]))


COMMON = peptide("ENSP_COMMON", "ACDEFGHIK")
DONOR = peptide("ENSP_DONOR", "LMNPQRSTV", "2", "200", "ALMNPQRSTVW")
RECIPIENT = peptide("ENSP_RECIPIENT", "WYACDEFGH", "3", "300", "MWYACDEFGHK")


@pytest.fixture
def deduce(tmp_path):
    def run(donor_rows, recipient_rows, orientation="dr", pair="P01"):
        donor_file = tmp_path / "donor.csv"
        recipient_file = tmp_path / "recipient.csv"
        for path, rows in [(donor_file, donor_rows), (recipient_file, recipient_rows)]:
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=COLUMNS)
                writer.writeheader()
                writer.writerows(rows)
        log_file = tmp_path / "run.log"
        log_file.write_text(f"Orientation: {orientation}\n", encoding="utf-8")
        args = SimpleNamespace(pair=pair, run_name="cleavage_test")
        output = Path(cleavage.deduce_cleaved_peptides(
            str(donor_file), str(recipient_file), str(tmp_path), args, str(log_file)
        ))
        with output.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            assert reader.fieldnames == COLUMNS
            rows = list(reader)
        return output, rows

    return run


@pytest.mark.parametrize("orientation,expected", [("dr", DONOR), ("rd", RECIPIENT)])
def test_keeps_only_source_specific_peptides(deduce, orientation, expected, capsys):
    # Shared k-mers are excluded even when their positions and long peptides differ.
    recipient_common = peptide("ENSP_COMMON", "ACDEFGHIK", "X", "900", "MACDEFGHIKL")
    _, rows = deduce([COMMON, DONOR], [recipient_common, RECIPIENT], orientation)
    assert rows == [expected]
    assert "kept 1/2, removed 1 (50.00%)." in capsys.readouterr().out


@pytest.mark.parametrize("orientation", ["dr", "rd"])
@pytest.mark.parametrize("scenario", ["identical", "disjoint", "empty_donor", "empty_recipient", "both_empty"])
def test_identical_disjoint_and_empty_inputs(deduce, orientation, scenario, capsys):
    cases = {
        "identical": ([COMMON], [COMMON], [], []),
        "disjoint": ([DONOR], [RECIPIENT], [DONOR], [RECIPIENT]),
        "empty_donor": ([], [RECIPIENT], [], [RECIPIENT]),
        "empty_recipient": ([DONOR], [], [DONOR], []),
        "both_empty": ([], [], [], []),
    }
    donor, recipient, expected_dr, expected_rd = cases[scenario]
    _, rows = deduce(donor, recipient, orientation)
    assert rows == (expected_dr if orientation == "dr" else expected_rd)
    if not (donor if orientation == "dr" else recipient):
        assert "kept 0/0, removed 0 (0.00%)." in capsys.readouterr().out


@pytest.mark.parametrize("orientation", ["dr", "rd"])
def test_comparison_uses_protein_and_kmer_not_sequence_alone(deduce, orientation):
    # Equal sequences from distinct proteins must remain source-specific.
    donor = peptide("ENSP_DONOR", "ACDEFGHIK")
    recipient = peptide("ENSP_RECIPIENT", "ACDEFGHIK")
    # Distinct k-mers within the same protein must also remain source-specific.
    donor_variant = peptide("ENSP_VARIANT", "LMNPQRSTV")
    recipient_variant = peptide("ENSP_VARIANT", "WYACDEFGH")
    _, rows = deduce([donor, donor_variant], [recipient, recipient_variant], orientation)
    assert rows == ([donor, donor_variant] if orientation == "dr" else [recipient, recipient_variant])


@pytest.mark.parametrize("orientation", ["dr", "rd"])
def test_deduplicates_exact_rows_but_preserves_distinct_positions(deduce, orientation, capsys):
    other_position = dict(DONOR, CHROM="X", POS="900")
    source = [DONOR, DONOR, other_position, other_position]
    donor, recipient = (source, []) if orientation == "dr" else ([], source)
    _, rows = deduce(donor, recipient, orientation)
    assert rows == [DONOR, other_position]
    # The summary counts unique (protein, k-mer) pairs, not genomic rows.
    assert "kept 1/1, removed 0 (0.00%)." in capsys.readouterr().out


@pytest.mark.parametrize("pair", ["", "P01", "P02"])
def test_output_filename_is_scoped_to_run_and_pair(deduce, tmp_path, pair):
    output, rows = deduce([DONOR], [], pair=pair)
    prefix = f"{pair}_" if pair else ""
    assert output == tmp_path / f"{prefix}cleavage_test_netchop_peptides.csv"
    assert rows == [DONOR]


def test_rejects_invalid_orientation_without_creating_output(deduce, tmp_path):
    with pytest.raises(ValueError, match="^Invalid orientation: invalid$"):
        deduce([DONOR], [RECIPIENT], orientation="invalid")
    assert not (tmp_path / "P01_cleavage_test_netchop_peptides.csv").exists()
