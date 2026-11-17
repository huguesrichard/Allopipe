"""Check that Nextflow gives its tasks a locale with a decimal point."""

import os
from pathlib import Path
import shutil
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[1]


def test_nextflow_tasks_override_french_parent_locale(tmp_path):
    nextflow = shutil.which("nextflow")
    framework = Path(os.environ.get("NXF_HOME", str(Path.home() / ".nextflow"))) / "framework"
    if not nextflow or not framework.is_dir():
        pytest.skip("Nextflow and a cached runtime are required")

    runtime = tmp_path / "nextflow"
    shutil.copytree(framework, runtime / "framework")
    script = tmp_path / "locale.nf"
    script.write_text(
        """process CHECK_LOCALE {
    output:
    stdout
    script:
    '''
    python -c 'import locale, os; locale.setlocale(locale.LC_ALL, ""); print("TASK_LC_ALL=" + os.environ.get("LC_ALL", "")); print("DECIMAL=" + locale.localeconv()["decimal_point"])'
    '''
}
workflow {
    CHECK_LOCALE()
    CHECK_LOCALE.out.view()
}
""", encoding="utf-8",
    )
    env = dict(os.environ, LC_ALL="fr_FR.UTF-8", LANG="fr_FR.UTF-8", NXF_HOME=str(runtime),
               NXF_OFFLINE="true", NXF_DISABLE_CHECK_LATEST="true", NXF_ANSI_LOG="false")
    result = subprocess.run(
        [nextflow, "-C", str(ROOT / "nextflow.config"), "run", str(script), "-work-dir", str(tmp_path / "work")],
        cwd=tmp_path, env=env, capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "TASK_LC_ALL=C" in result.stdout
    assert "DECIMAL=." in result.stdout
