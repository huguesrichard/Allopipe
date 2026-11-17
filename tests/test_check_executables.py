"""Exercise preflight checks using deterministic executable substitutes."""

import sys
from pathlib import Path
import subprocess

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from tools.check_executables import check_executables


def test_missing_executable_is_reported(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(tmp_path))
    with pytest.raises(ValueError, match="netMHCpan.*not found in PATH"):
        check_executables(["netMHCpan"])


@pytest.mark.parametrize("message", ["no binaries found for Darwin_arm64", "Bad CPU type in executable"])
def test_wrapper_that_cannot_execute_is_reported(tmp_path, monkeypatch, message):
    executable = tmp_path / "netchop"
    executable.write_text(f"#!/bin/sh\nprintf '%s\\n' '{message}'\n")
    executable.chmod(0o755)
    monkeypatch.setenv("PATH", str(tmp_path))
    with pytest.raises(ValueError, match="netchop.*cannot run"):
        check_executables(["netchop"])


def test_usage_exit_does_not_mean_executable_is_broken(tmp_path, monkeypatch):
    executable = tmp_path / "netMHCpan"
    executable.write_text("#!/bin/sh\necho 'usage: netMHCpan -f input'\nexit 1\n")
    executable.chmod(0o755)
    monkeypatch.setenv("PATH", str(tmp_path))
    check_executables(["netMHCpan"])


def test_wrapper_execution_failure_is_rejected(tmp_path, monkeypatch):
    executable = tmp_path / "netchop"
    executable.write_text("#!/bin/sh\nexit 126\n")
    executable.chmod(0o755)
    monkeypatch.setenv("PATH", str(tmp_path))
    with pytest.raises(ValueError, match="netchop.*cannot run"):
        check_executables(["netchop"])


def test_probe_timeout_is_reported(monkeypatch):
    monkeypatch.setattr("tools.check_executables.shutil.which", lambda name: "/fake/netchop")
    def timeout(*args, **kwargs):
        raise subprocess.TimeoutExpired(args[0], 15)
    monkeypatch.setattr("tools.check_executables.subprocess.run", timeout)
    with pytest.raises(ValueError, match="netchop.*could not be checked"):
        check_executables(["netchop"])
