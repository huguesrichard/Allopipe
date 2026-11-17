"""Validation uses the scoring parsers before generated inputs exist."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tools.validate_arguments import validate


@pytest.fixture
def validation_request():
    return {"run_name": "test", "orientation": "rd", "imputation": "no-imputation",
            "mode": "cohort", "frameshift": False, "ensembl_path": "/not/downloaded",
            "count_opts": "", "affinity_opts": "", "pairs": [
                {"pair_id": "P1", "hla_typing": "HLA-A*02:01"},
                {"pair_id": "P2", "hla_typing": "HLA-B*08:01"},
            ]}


def test_accepts_valid_arguments_before_generated_files_exist(validation_request):
    validate(validation_request)


def test_reports_the_pair_with_invalid_hla(validation_request, capsys):
    validation_request["pairs"][1]["hla_typing"] = "bad"
    with pytest.raises(SystemExit):
        validate(validation_request)
    assert "argument validation failed for P2" in capsys.readouterr().err


def test_rejects_invalid_quoting(validation_request):
    validation_request["affinity_opts"] = "--length '9"
    with pytest.raises(ValueError, match="Invalid quoting"):
        validate(validation_request)


def test_options_are_returned_as_shell_quoted_arguments(validation_request):
    validation_request["affinity_opts"] = '--length "9"'
    assert validate(validation_request)["affinity_opts"] == "--length 9"


def test_length_is_validated_for_the_selected_hla_class(validation_request):
    validation_request["pairs"] = [{"pair_id": "P1", "hla_typing": "DRB1_0101"}]
    validation_request["affinity_opts"] = "--class_type 2 --length 8"
    with pytest.raises(SystemExit):
        validate(validation_request)


@pytest.mark.parametrize("options, hla, expected", [
    ("", "HLA-A*02:01", ["netMHCpan"]),
    ("--class_type 2", "DRB1_0101", ["netMHCIIpan"]),
    ("--dry_run", "HLA-A*02:01", []),
    ("--cleavage", "HLA-A*02:01", ["netMHCpan", "netchop"]),
    ("--dry_run --cleavage", "HLA-A*02:01", ["netchop"]),
])
def test_required_executables_follow_validated_options(validation_request, options, hla, expected):
    validation_request["affinity_opts"] = options
    validation_request["pairs"] = [{"pair_id": "P1", "hla_typing": hla}]
    assert validate(validation_request)["required_executables"] == expected
