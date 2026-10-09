import os
import sys
import base64
from pathlib import Path
import pytest


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

# normalize_cohort is also a standalone Nextflow script with sibling imports.
sys.path.insert(0, str(SRC / "tools"))


@pytest.fixture
def nextflow_environment(tmp_path, monkeypatch):
    published_dir = tmp_path / "published"
    command = "nextflow run main.nf --run_name test --output_dir " + str(published_dir)
    monkeypatch.setenv("ALLOPIPE_PUBLISHED_OUTPUT_DIR", str(published_dir))
    monkeypatch.setenv("ALLOPIPE_NEXTFLOW_COMMAND_BASE64", base64.b64encode(command.encode()).decode())
    monkeypatch.setenv("ALLOPIPE_VERSION", "v-test")
    return published_dir, command


def pytest_sessionstart(session):
    # Make path-dependent code deterministic in tests.
    os.chdir(ROOT)
