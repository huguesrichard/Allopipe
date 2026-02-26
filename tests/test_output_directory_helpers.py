from pathlib import Path
from types import SimpleNamespace

from tools import aams_helpers, ams_helpers, multivcf_extract


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


def _build_args_for_overwrite(output_dir: Path):
    return SimpleNamespace(
        output_dir=str(output_dir),
        run_name="runX",
        min_dp=20,
        max_dp=400,
        min_ad=5,
        homozygosity_thr=0.2,
        min_gq=0,
        orientation="dr",
        base_length=3,
        pair="P01",
    )


def test_handle_overwrite_true_when_matching_csv_exists(tmp_path):
    args = _build_args_for_overwrite(tmp_path / "output")
    csv_dir = (
        Path(args.output_dir)
        / "runs"
        / args.run_name
        / "AMS"
        / f"{args.run_name}_AMS_{args.min_dp}_{args.max_dp}_{args.min_ad}_{args.homozygosity_thr}_{args.min_gq}_{args.orientation}_{args.base_length}"
    )
    csv_dir.mkdir(parents=True, exist_ok=True)
    csv_file = (
        csv_dir
        / f"AMS_{args.run_name}_{args.pair}_{args.min_dp}_{args.max_dp}_{args.min_ad}_{args.homozygosity_thr}_{args.min_gq}_{args.orientation}_{args.base_length}.csv"
    )
    csv_file.write_text("pair,ams\nP01,1\n", encoding="utf-8")
    assert ams_helpers.handle_overwrite(args) is True


def test_handle_overwrite_false_when_missing(tmp_path):
    args = _build_args_for_overwrite(tmp_path / "output")
    assert ams_helpers.handle_overwrite(args) is False


def test_multivcf_create_dependencies_uses_output_dir(tmp_path):
    output_dir = tmp_path / "root_out"
    path = multivcf_extract.create_dependencies("runY", str(output_dir))
    assert path == str(output_dir / "runs" / "runY" / "vcf_indiv")
    assert Path(path).is_dir()


def test_aams_create_dependencies_under_output_dir(tmp_path):
    output_dir = tmp_path / "root_out"
    aams_run_tables, netmhc_dir, aams_path, netchop_dir = aams_helpers.create_aams_dependencies(
        "runZ", str(output_dir)
    )
    assert Path(aams_run_tables).is_dir()
    assert Path(netmhc_dir).is_dir()
    assert Path(aams_path).is_dir()
    assert Path(netchop_dir).is_dir()
    assert str(output_dir / "runs" / "runZ") in aams_run_tables
