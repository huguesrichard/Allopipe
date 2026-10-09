"""Assembly of per-pair worker outputs by the Nextflow cohort stage."""
from pathlib import Path
import subprocess
import sys

import pandas as pd
import pytest

from tools import normalize_cohort


RUN_NAME = "cohort_test"


@pytest.fixture
def pair_runs(tmp_path):
    """Two staged tasks, with two score configurations and distinct pair files."""
    staged = []
    for slot, pair, ams, aams in [("01", "P10", 30, 7), ("02", "P2", 12, 3)]:
        stage_dir = tmp_path / "pair_runs" / slot
        run = stage_dir / "runs" / RUN_NAME
        tables = run / "run_tables"
        tables.mkdir(parents=True)
        donor = pd.DataFrame({"CHROM": ["1", "1"], "POS": [100, 200], "GT": ["0/0", "0/0"]})
        recipient = donor.assign(GT=["0/0", "0/1"])
        donor.to_csv(tables / f"{pair}_D0_sample.tsv", sep="\t", index=False)
        recipient.to_csv(tables / f"{pair}_R0_sample.tsv", sep="\t", index=False)
        # Intermediate peptide pickles must not enter the final AAMS summary.
        pd.DataFrame({"peptide": ["ACDEFGHIK"]}).to_pickle(tables / f"{pair}_pep_df.pkl")
        for profile, offset in [("default", 0), ("strict", 100)]:
            ams_dir = run / "AMS" / profile
            aams_dir = run / "AAMS" / profile
            ams_dir.mkdir(parents=True)
            aams_dir.mkdir(parents=True)
            pd.DataFrame({"pair": [pair], "ams": [ams + offset]}).to_pickle(
                ams_dir / f"AMS_{RUN_NAME}_{pair}_{profile}.pkl"
            )
            pd.DataFrame({"pair": [pair], "aams": [aams + offset]}).to_pickle(
                aams_dir / f"{pair}_{RUN_NAME}_AAMS_df.pkl"
            )
            pd.DataFrame({"not_a_score": [1]}).to_pickle(aams_dir / f"{pair}_intermediate.pkl")
        logs = run / "logs"
        logs.mkdir()
        (logs / f"{pair}_run.log").write_text(f"Pair: {pair}\n", encoding="utf-8")
        staged.append(stage_dir)
    return staged


def file_contents(root):
    return {path.relative_to(root): path.read_bytes() for path in root.rglob("*") if path.is_file()}


@pytest.mark.parametrize("layout", ["named_child", "runs_child", "run_itself", "nested_stage"])
def test_finds_supported_staged_run_layouts(tmp_path, layout):
    root = tmp_path / "stage"
    run = {
        "named_child": root / RUN_NAME,
        "runs_child": root / "runs" / RUN_NAME,
        "run_itself": root / RUN_NAME,
        "nested_stage": root / "01" / "worker" / "runs" / RUN_NAME,
    }[layout]
    run.mkdir(parents=True)
    supplied = run if layout == "run_itself" else root
    assert normalize_cohort.find_run_path(supplied, RUN_NAME) == run


def test_missing_run_reports_name_and_search_directory(tmp_path):
    with pytest.raises(FileNotFoundError) as error:
        normalize_cohort.find_run_path(tmp_path, RUN_NAME)
    assert str(error.value) == f"No such run directory for {RUN_NAME} under {tmp_path}"


def test_collects_exact_score_and_genotype_files_across_tasks(pair_runs):
    collected = normalize_cohort.collect_run_files(pair_runs, RUN_NAME)
    expected = [set(), set(), set(), set()]
    for stage, pair in zip(pair_runs, ["P10", "P2"]):
        run = stage / "runs" / RUN_NAME
        for profile in ["default", "strict"]:
            expected[0].add(str(run / "AMS" / profile / f"AMS_{RUN_NAME}_{pair}_{profile}.pkl"))
            expected[1].add(str(run / "AAMS" / profile / f"{pair}_{RUN_NAME}_AAMS_df.pkl"))
        expected[2].add(str(run / "run_tables" / f"{pair}_D0_sample.tsv"))
        expected[3].add(str(run / "run_tables" / f"{pair}_R0_sample.tsv"))
    assert [set(paths) for paths in collected] == expected
    assert all(len(paths) == len(set(paths)) for paths in collected)


@pytest.mark.parametrize("missing,expected_message", [
    ("ams", "AMS pickle files"),
    ("donor", "donor tables"),
    ("recipient", "recipient tables"),
    ("all", "AMS pickle files, donor tables, recipient tables"),
])
def test_reports_missing_required_file_categories(tmp_path, missing, expected_message):
    run = tmp_path / "runs" / RUN_NAME
    for category, relative in [
        ("ams", "AMS/default/AMS_P2.pkl"),
        ("donor", "run_tables/P2_D0_sample.tsv"),
        ("recipient", "run_tables/P2_R0_sample.tsv"),
    ]:
        path = run / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        if missing not in [category, "all"]:
            path.touch()  # Collection inspects filenames, not file contents.
    with pytest.raises(FileNotFoundError) as error:
        normalize_cohort.collect_run_files([tmp_path], RUN_NAME)
    assert str(error.value) == "Missing expected output files: " + expected_message


def test_merges_all_task_files_and_preserves_sources(pair_runs, tmp_path):
    sources = [stage / "runs" / RUN_NAME for stage in pair_runs]
    before = [file_contents(source) for source in sources]
    expected = {relative: content for contents in before for relative, content in contents.items()}
    output_dir = tmp_path / "published"
    final_run = normalize_cohort.merge_run_dirs(pair_runs, RUN_NAME, output_dir)
    assert final_run == output_dir / "runs" / RUN_NAME
    assert file_contents(final_run) == expected
    assert [file_contents(source) for source in sources] == before
    # Repeated assembly must not remove an existing unrelated output.
    retained = final_run / "retained.txt"
    retained.write_text("keep me", encoding="utf-8")
    assert normalize_cohort.merge_run_dirs(pair_runs, RUN_NAME, output_dir) == final_run
    assert file_contents(final_run) == {
        **expected, Path("retained.txt"): b"keep me",
    }


@pytest.mark.parametrize("with_aams", [True, False])
def test_cli_assembles_and_writes_cohort_tables(pair_runs, tmp_path, with_aams):
    if not with_aams:
        for stage in pair_runs:
            for path in (stage / "runs" / RUN_NAME / "AAMS").rglob("*AAMS_df*.pkl"):
                path.unlink()
    sources = [stage / "runs" / RUN_NAME for stage in pair_runs]
    before = [file_contents(source) for source in sources]
    output_dir = tmp_path / "published"
    worker_cwd = tmp_path / "worker_cwd"
    worker_cwd.mkdir()
    result = subprocess.run([
        sys.executable, str(Path(normalize_cohort.__file__).resolve()),
        "--run-name", RUN_NAME, "--output-dir", str(output_dir),
        "--run-dir", *(str(stage) for stage in pair_runs),
    ], cwd=worker_cwd, capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
    final_run = output_dir / "runs" / RUN_NAME
    for contents in before:
        for relative, content in contents.items():
            assert (final_run / relative).read_bytes() == content
    for profile, offset in [("default", 0), ("strict", 100)]:
        ams = pd.read_csv(final_run / "AMS" / profile / "AMS_df.tsv", sep="\t")
        expected_ams = pd.DataFrame({
            "pair": ["P2", "P10"], "ams_giab": [12 + offset, 30 + offset],
            "common_ref": [1, 1], "total_ref": [2, 2], "ref_ratio": [0.5, 0.5],
            "ams_norm": [12 + offset, 30 + offset],
        })
        pd.testing.assert_frame_equal(ams, expected_ams)
        aams_file = final_run / "AAMS" / profile / "AAMS_df.tsv"
        if with_aams:
            expected_aams = pd.DataFrame({"pair": ["P2", "P10"], "aams": [3 + offset, 7 + offset]})
            pd.testing.assert_frame_equal(pd.read_csv(aams_file, sep="\t"), expected_aams)
        else:
            assert not aams_file.exists()
    assert [file_contents(source) for source in sources] == before
