from pathlib import Path
import sys

import pytest

from tools import arguments_handling


def _write_dummy_vcf(path: Path) -> None:
    path.write_text("##fileformat=VCFv4.2\n", encoding="utf-8")


def test_arguments_output_dir_default_is_absolute(tmp_path, monkeypatch):
    donor = tmp_path / "donor.vcf"
    recipient = tmp_path / "recipient.vcf"
    _write_dummy_vcf(donor)
    _write_dummy_vcf(recipient)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "ams_pipeline.py",
            str(donor),
            str(recipient),
            "dr",
            "imputation",
        ],
    )
    args = arguments_handling.arguments()
    assert Path(args.output_dir).is_absolute()
    assert Path(args.output_dir).name == "output"


def test_arguments_output_dir_custom_is_absolute(tmp_path, monkeypatch):
    donor = tmp_path / "donor.vcf"
    recipient = tmp_path / "recipient.vcf"
    _write_dummy_vcf(donor)
    _write_dummy_vcf(recipient)
    custom_out = tmp_path / "my-output"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "ams_pipeline.py",
            str(donor),
            str(recipient),
            "rd",
            "no-imputation",
            "--output_dir",
            str(custom_out),
        ],
    )
    args = arguments_handling.arguments()
    assert args.output_dir == str(custom_out.resolve())


def test_check_file_rejects_wrong_extension(tmp_path):
    wrong_file = tmp_path / "sample.txt"
    wrong_file.write_text("x", encoding="utf-8")
    parser = arguments_handling.CustomParser(prog="prog")
    with pytest.raises(SystemExit):
        arguments_handling.check_file(parser, str(wrong_file))


def test_check_threshold_value_bounds_homozygosity():
    parser = arguments_handling.CustomParser(prog="prog")
    assert arguments_handling.check_threshold_value(parser, "0.2", "homozygosity_thr") == 0.2
    with pytest.raises(SystemExit):
        arguments_handling.check_threshold_value(parser, "0", "homozygosity_thr")


def test_check_workers_count_bounds(monkeypatch):
    parser = arguments_handling.CustomParser(prog="prog")
    monkeypatch.setattr(arguments_handling.os, "cpu_count", lambda: 4)
    assert arguments_handling.check_workers_count(parser, "2") == 2
    with pytest.raises(SystemExit):
        arguments_handling.check_workers_count(parser, "8")


def test_arguments_raise_when_min_dp_gt_max_dp(tmp_path, monkeypatch):
    donor = tmp_path / "donor.vcf"
    recipient = tmp_path / "recipient.vcf"
    _write_dummy_vcf(donor)
    _write_dummy_vcf(recipient)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "ams_pipeline.py",
            str(donor),
            str(recipient),
            "dr",
            "imputation",
            "--min_dp",
            "30",
            "--max_dp",
            "20",
        ],
    )
    with pytest.raises(ValueError, match="minimal Depth"):
        arguments_handling.arguments()


def test_arguments_raise_when_min_ad_gt_min_dp(tmp_path, monkeypatch):
    donor = tmp_path / "donor.vcf"
    recipient = tmp_path / "recipient.vcf"
    _write_dummy_vcf(donor)
    _write_dummy_vcf(recipient)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "ams_pipeline.py",
            str(donor),
            str(recipient),
            "dr",
            "imputation",
            "--min_dp",
            "20",
            "--min_ad",
            "21",
        ],
    )
    with pytest.raises(ValueError, match="minimal Allelic Depth"):
        arguments_handling.arguments()
