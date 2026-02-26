from pathlib import Path
import sys

import pytest

from tools import netmhc_arguments


def test_netmhc_arguments_output_dir_and_run_name_validation(tmp_path, monkeypatch):
    ensembl_dir = tmp_path / "ensembl"
    ensembl_dir.mkdir()
    output_dir = tmp_path / "custom_output"
    (output_dir / "runs" / "runA").mkdir(parents=True)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "aams_pipeline.py",
            "-d",
            str(ensembl_dir),
            "-n",
            "runA",
            "-a",
            "HLA-A*02:01",
            "--output_dir",
            str(output_dir),
        ],
    )
    args = netmhc_arguments.netmhc_arguments()
    assert args.output_dir == str(output_dir.resolve())
    assert args.run_name == "runA"


def test_netmhc_arguments_missing_run_raises(tmp_path, monkeypatch):
    ensembl_dir = tmp_path / "ensembl"
    ensembl_dir.mkdir()
    output_dir = tmp_path / "custom_output"
    output_dir.mkdir()
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "aams_pipeline.py",
            "-d",
            str(ensembl_dir),
            "-n",
            "unknown_run",
            "-a",
            "HLA-A*02:01",
            "--output_dir",
            str(output_dir),
        ],
    )
    with pytest.raises(SystemExit):
        netmhc_arguments.netmhc_arguments()


def test_netmhc_arguments_length_validation_class1(tmp_path, monkeypatch):
    ensembl_dir = tmp_path / "ensembl"
    ensembl_dir.mkdir()
    output_dir = tmp_path / "out"
    (output_dir / "runs" / "run1").mkdir(parents=True)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "aams_pipeline.py",
            "-d",
            str(ensembl_dir),
            "-n",
            "run1",
            "-a",
            "HLA-A*02:01",
            "--output_dir",
            str(output_dir),
            "--length",
            "7",
        ],
    )
    with pytest.raises(SystemExit):
        netmhc_arguments.netmhc_arguments()
