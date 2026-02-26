from pathlib import Path
from types import SimpleNamespace

import pandas as pd

from tools import aams_helpers, cleavage


def test_get_ams_params_parses_mismatches_file_from_output_dir(tmp_path):
    output_dir = tmp_path / "out"
    run_tables = output_dir / "runs" / "runA" / "run_tables"
    run_tables.mkdir(parents=True)
    mismatch_file = run_tables / "P01_runA_mismatches_20_400_5_gq_0_0.2_bl_3.tsv"
    mismatch_file.write_text("CHROM\tPOS\n1\t10\n", encoding="utf-8")
    str_params, mismatches_path = aams_helpers.get_ams_params("runA", str(output_dir))
    assert str_params == "20_400_5_0_0.2_3"
    assert mismatches_path == str(mismatch_file)


def test_read_log_field_uses_output_dir(tmp_path):
    output_dir = tmp_path / "out"
    logs_dir = output_dir / "runs" / "runB" / "logs"
    logs_dir.mkdir(parents=True)
    log_file = logs_dir / "run.log"
    log_file.write_text("Orientation: dr\nDonor: /tmp/donor.vcf.gz\n", encoding="utf-8")
    args = SimpleNamespace(output_dir=str(output_dir), run_name="runB", pair="")
    assert aams_helpers.read_log_field(args, "Orientation") == "dr"
    assert aams_helpers.read_log_field(args, "Donor") == "/tmp/donor.vcf.gz"


def test_pickle_parsing_reads_from_output_dir_and_extracts_fields(tmp_path):
    output_dir = tmp_path / "out"
    run_name = "runC"
    logs_dir = output_dir / "runs" / run_name / "logs"
    run_tables = output_dir / "runs" / run_name / "run_tables"
    logs_dir.mkdir(parents=True)
    run_tables.mkdir(parents=True)

    donor_path = "/tmp/SAMPLE_DONOR.vcf.gz"
    (logs_dir / "run.log").write_text(
        f"Orientation: dr\nDonor: {donor_path}\nRecipient: /tmp/r.vcf.gz\n",
        encoding="utf-8",
    )

    str_params = "20_400_5_0_0.2_3"
    str_params_split = "20_400_5_0.2"
    pickle_path = run_tables / f"SAMPLE_DONOR_vep_infos_table_{str_params_split}.pkl"
    df = pd.DataFrame(
        {
            "CHROM": ["1"],
            "POS": [100],
            "INFO": ["0|1|2|3|ENSG0001|5|ENST0001|7|8|9|10|11|12|13|42"],
        }
    )
    df.to_pickle(pickle_path)

    args = SimpleNamespace(output_dir=str(output_dir), run_name=run_name, pair="")
    parsed = cleavage.pickle_parsing(str_params, args)
    assert "INFO" not in parsed.columns
    assert parsed.loc[0, "Gene_id"] == "ENSG0001"
    assert parsed.loc[0, "Transcript_id"] == "ENST0001"
    assert parsed.loc[0, "Protein_position"] == "42"


def test_parse_netchop_output_min_run_length(tmp_path):
    netchop_out = tmp_path / "netchop.txt"
    netchop_out.write_text(
        "\n".join(
            [
                "header",
                "-----",
                "-----",
                "Pos AA C S Ident",
                "1 A . 0.1 ID1",
                "2 B . 0.1 ID1",
                "3 C S 0.1 ID1",
            ]
        ),
        encoding="utf-8",
    )
    runs, all_positions = cleavage.parse_netchop_output(str(netchop_out), min_run_length=2)
    assert "ID1" in runs
    assert runs["ID1"][0]["positions"] == [1, 2]
    assert all_positions["ID1"] == [(1, "A"), (2, "B"), (3, "C")]
