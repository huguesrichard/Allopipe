from pathlib import Path

import pytest

from tools import aams_helpers, ams_helpers


def test_create_run_directory_under_output_dir(tmp_path):
    output_dir = tmp_path / "out"
    run_path, run_tables, run_plots, run_ams, run_logs = ams_helpers.create_run_directory(
        "runX", str(output_dir)
    )
    assert run_path == str(output_dir / "runs" / "runX")
    assert Path(run_tables).is_dir()
    assert Path(run_plots).is_dir()
    assert Path(run_ams).is_dir()
    assert Path(run_logs).is_dir()


@pytest.mark.parametrize("cleavage", [False, True])
def test_aams_create_dependencies_under_output_dir(tmp_path, cleavage):
    output_dir = tmp_path / "root_out"
    aams_run_tables, netmhc_dir, aams_path, netchop_dir = aams_helpers.create_aams_dependencies(
        "runZ", str(output_dir), cleavage=cleavage
    )
    assert Path(aams_run_tables).is_dir()
    assert Path(netmhc_dir).is_dir()
    assert Path(aams_path).is_dir()
    assert Path(netchop_dir).is_dir() is cleavage
    assert str(output_dir / "runs" / "runZ") in aams_run_tables


def test_disabling_cleavage_preserves_existing_netchop_results(tmp_path):
    _, _, _, netchop_dir = aams_helpers.create_aams_dependencies("runZ", str(tmp_path), cleavage=True)
    result = Path(netchop_dir) / "previous.txt"
    result.write_text("keep")
    aams_helpers.create_aams_dependencies("runZ", str(tmp_path), cleavage=False)
    assert result.read_text() == "keep"
