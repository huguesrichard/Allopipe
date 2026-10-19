"""Run the real Nextflow workflow with deterministic scoring/VEP substitutes."""
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
RUN_NAME = "orchestration_test"


@pytest.fixture(scope="session")
def nextflow_runtime(tmp_path_factory):
    executable = shutil.which("nextflow")
    if not executable:
        pytest.skip("Nextflow integration tests require Nextflow")
    cache = Path(os.environ.get("NXF_HOME", str(Path.home() / ".nextflow")))
    framework = cache / "framework"
    if not framework.is_dir():
        pytest.skip("Nextflow integration tests require an already cached runtime (offline execution)")
    runtime = tmp_path_factory.mktemp("nextflow-runtime")
    shutil.copytree(framework, runtime / "framework")
    bcftools = Path(sys.executable).parent / "bcftools"
    if not bcftools.is_file():
        found = shutil.which("bcftools")
        if not found:
            pytest.skip("Nextflow integration tests require bcftools")
        bcftools = Path(found)
    return executable, runtime, bcftools


@pytest.fixture
def workflow(tmp_path, nextflow_runtime):
    executable, runtime, bcftools = nextflow_runtime
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    multi_vcf = inputs / "cohort.vcf"
    multi_vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=1,length=1000>\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\n"
        "1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/0\t0/0\t0/0\n"
        "1\t200\t.\tG\tC\t.\tPASS\t.\tGT\t0/1\t0/1\t0/1\n",
        encoding="utf-8",
    )
    pairs = inputs / "pairs.csv"
    pairs.write_text(
        'donor,recipient,hla\nA,B,"HLA-A02:01,HLA-B07:02"\n'
        "A,C,HLA-A03:01,HLA-C07:02\nC,A,HLA-B08:01\n", encoding="utf-8",
    )
    # The main workflow prefixes this path with projectDir.
    ensembl = tmp_path / "GRCh38"
    ensembl.mkdir()
    relative_ensembl = os.path.relpath(ensembl, ROOT)
    binaries = tmp_path / "bin"
    binaries.mkdir()
    for name, target in [("python", ""), ("vep", "vep ")]:
        wrapper = binaries / name
        wrapper.write_text(
            '#!/bin/sh\nexec "$ALLOPIPE_TEST_PYTHON" "$ALLOPIPE_TEST_STUB" '
            + target + '"$@"\n', encoding="utf-8",
        )
        wrapper.chmod(0o755)
    output = tmp_path / "published"
    config = tmp_path / "test.config"
    config.write_text(
        "process.executor = 'local'\nprocess.maxForks = 2\n"
        "process.container = null\nprocess.shell = ['/bin/bash', '-euo', 'pipefail']\n"
        "docker.enabled = false\nsingularity.enabled = false\nconda.enabled = false\n"
        "params.vep_version = 'test'\nparams.vep_cache = '/unused'\n"
        "params.frameshift = false\nparams.frameshift_plugin_path = null\n"
        "params.force_overwrite = false\nparams.allo_count_opts = ''\nparams.allo_affinity_opts = ''\n",
        encoding="utf-8",
    )
    env = dict(os.environ)
    env.update({
        "NXF_HOME": str(runtime), "NXF_OFFLINE": "true", "NXF_ANSI_LOG": "false",
        "NXF_DISABLE_CHECK_LATEST": "true", "ALLOPIPE_TEST_PYTHON": sys.executable,
        "ALLOPIPE_TEST_STUB": str(ROOT / "tests" / "fixtures" / "nextflow_task_stub.py"),
        "PATH": str(binaries) + os.pathsep + str(bcftools.parent) + os.pathsep + env.get("PATH", ""),
    })

    def run(*extra, csv_text=None, skip_vep=True, mode="cohort"):
        if csv_text is not None:
            pairs.write_text(csv_text, encoding="utf-8")
        command = [
            executable, "-C", str(config), "run", str(ROOT / "main.nf"),
            "-work-dir", str(tmp_path / "work"), "-with-trace", str(tmp_path / "trace.tsv"),
            "--mode", mode, "--run_name", RUN_NAME, "--orientation", "dr",
            "--imputation", "imputation", "--ensembl_path", relative_ensembl,
            "--output_dir", str(output), "--multi_vcf", str(multi_vcf), "--pairs", str(pairs),
            "--skip_vep_annotation", str(skip_vep).lower(), *extra,
        ]
        if mode == "pair":
            for sample in ["A", "B"]:
                selected = inputs / f"{sample}.vcf"
                extraction = subprocess.run([
                    str(bcftools), "view", "-s", sample, "-O", "v", "-o", str(selected), str(multi_vcf),
                ], capture_output=True, text=True, timeout=30)
                assert extraction.returncode == 0, extraction.stderr
            command += ["--donor", str(inputs / "A.vcf"), "--recipient", str(inputs / "B.vcf")]
        result = subprocess.run(command, cwd=tmp_path, env=env, capture_output=True, text=True, timeout=60)
        (tmp_path / "last-run.stdout").write_text(result.stdout, encoding="utf-8")
        (tmp_path / "last-run.stderr").write_text(result.stderr, encoding="utf-8")
        return result

    return run, output / "runs" / RUN_NAME, tmp_path


@pytest.mark.parametrize("skip_vep", [True, False], ids=["preannotated", "vep-stage"])
def test_cohort_routes_shared_samples_and_per_pair_hla_and_publishes_results(workflow, skip_vep):
    run, published, folder = workflow
    result = run(skip_vep=skip_vep)
    assert result.returncode == 0, result.stdout + result.stderr
    expected_pairs = [
        ("P1", "A", "B", "HLA-A02:01,HLA-B07:02"),
        ("P2", "A", "C", "HLA-A03:01,HLA-C07:02"),
        ("P3", "C", "A", "HLA-B08:01"),
    ]
    for pair, donor, recipient, hla in expected_pairs:
        count = json.loads((published / "logs" / f"{pair}_count.json").read_text())
        affinity = json.loads((published / "logs" / f"{pair}_affinity.json").read_text())
        assert count["pair"] == affinity["pair"] == pair
        assert count["donor"] == donor and count["recipient"] == recipient
        assert affinity == dict(count, hla=hla, ensembl_exists=True)
        assert count["output_dir"] == str(published.parent.parent)
        assert count["orientation"] == "dr" and count["imputation"] == "imputation"
        assert count["frameshift"] is False
        assert count["mode"] == "cohort"
        assert (published / "plots" / f"{pair}_pie_chart.png").read_bytes() == b"pair plot placeholder"
        assert "--run_name" in count["nextflow_command"]
    ams = pd.read_csv(published / "AMS" / "test" / "AMS_df.tsv", sep="\t")
    assert ams["pair"].tolist() == ["P1", "P2", "P3"]
    assert ams["ams_giab"].tolist() == ams["ams_norm"].tolist() == [1, 2, 3]
    assert (published / "plots" / "distrib.png").is_file()
    aams = pd.read_csv(published / "AAMS" / "test" / "AAMS_df.tsv", sep="\t")
    pd.testing.assert_frame_equal(aams, pd.DataFrame({"pair": ["P1", "P2", "P3"], "aams": [10, 20, 30]}))
    assert {path.name for path in (published / "vcf_indiv").iterdir()} == {
        f"{sample}.vcf.gz{suffix}" for sample in ["A", "B", "C"] for suffix in ["", ".tbi"]
    }
    with (folder / "trace.tsv").open() as handle:
        trace = list(csv.DictReader(handle, delimiter="\t"))
    assert sum("EXTRACT_SAMPLE" in row["name"] for row in trace) == 3
    assert sum("ALLO_COUNT" in row["name"] for row in trace) == 3
    assert sum("ALLO_AFFINITY" in row["name"] for row in trace) == 3
    assert sum("NORMALIZE_COHORT" in row["name"] for row in trace) == 1
    assert sum("VEP_ANNOTATION" in row["name"] for row in trace) == (0 if skip_vep else 1)
    assert bool(list((published / "vcf_vep").glob("*.vcf.gz"))) == (not skip_vep)


@pytest.mark.parametrize("csv_text,error", [
    ("donor,recipient,hla\n", "header and at least one donor/recipient row"),
    ("sample,recipient,hla\nA,B,HLA-A02:01\n", "donor and recipient columns"),
    ("donor,recipient\nA,B\n", "an hla column"),
    ("donor,hla,recipient\nA,HLA-A02:01,B\n", "hla as the last column"),
    ("donor,recipient,hla\nA\n", "row 2 is missing donor/recipient values"),
    ("donor,recipient,hla\nA,B,\n", "row 2 has an empty hla value"),
])
def test_invalid_pair_csv_fails_before_launching_tasks(workflow, csv_text, error):
    run, published, folder = workflow
    result = run(csv_text=csv_text)
    assert result.returncode != 0
    assert error in result.stdout + result.stderr
    assert not published.exists()
    assert not list((folder / "work").rglob(".command.sh"))


def test_existing_run_is_refused_without_overwrite(workflow):
    run, published, _ = workflow
    published.mkdir(parents=True)
    sentinel = published / "keep.txt"
    sentinel.write_text("keep", encoding="utf-8")
    result = run()
    assert result.returncode != 0
    assert "output run directory already exists" in result.stdout + result.stderr
    assert sentinel.read_text() == "keep"


def test_resume_reuses_sample_extraction_and_publishes_results(workflow):
    run, published, folder = workflow
    first = run()
    assert first.returncode == 0, first.stdout + first.stderr
    # Remove the previous trace only; the Nextflow task cache remains intact.
    (folder / "trace.tsv").unlink()
    resumed = run("-resume")
    assert resumed.returncode == 0, resumed.stdout + resumed.stderr
    assert "RESUME: keeping existing output run directory" in resumed.stdout + resumed.stderr
    with (folder / "trace.tsv").open() as handle:
        trace = list(csv.DictReader(handle, delimiter="\t"))
    extracted = [row for row in trace if "EXTRACT_SAMPLE" in row["name"]]
    assert len(extracted) == 3 and all(row["status"] == "CACHED" for row in extracted)
    # Count receives workflow.commandLine, which changes when -resume is added;
    # do not require all downstream stages to retain the same cache key.
    assert (published / "AMS" / "test" / "AMS_df.tsv").is_file()


def test_force_overwrite_replaces_only_the_selected_run(workflow):
    run, published, _ = workflow
    published.mkdir(parents=True)
    stale = published / "stale.txt"
    stale.write_text("old", encoding="utf-8")
    other_run = published.parent / "other_run"
    other_run.mkdir()
    sentinel = other_run / "keep.txt"
    sentinel.write_text("keep", encoding="utf-8")
    result = run("--force_overwrite")
    assert result.returncode == 0, result.stdout + result.stderr
    assert not stale.exists()
    assert (published / "AMS" / "test" / "AMS_df.tsv").is_file()
    assert sentinel.read_text() == "keep"


@pytest.mark.parametrize("kind", ["file", "symlink"])
def test_force_overwrite_refuses_non_directory_or_symlink(workflow, kind):
    run, published, folder = workflow
    published.parent.mkdir(parents=True)
    protected_dir = folder / "untouched"
    protected_dir.mkdir()
    target = protected_dir / "keep.txt"
    target.write_text("keep", encoding="utf-8")
    if kind == "file":
        published.write_text("keep", encoding="utf-8")
    else:
        published.symlink_to(protected_dir, target_is_directory=True)
    result = run("--force_overwrite")
    assert result.returncode != 0
    assert "refusing to overwrite" in result.stdout + result.stderr
    assert target.read_text() == "keep"
    if kind == "symlink":
        assert published.is_symlink()
        assert (published / "keep.txt").read_text() == "keep"
    else:
        assert published.is_file() and published.read_text() == "keep"


def test_run_name_cannot_escape_the_output_directory(workflow):
    run, published, _ = workflow
    neighbour = published.parent.parent / "escape"
    neighbour.mkdir(parents=True)
    sentinel = neighbour / "keep.txt"
    sentinel.write_text("keep", encoding="utf-8")
    result = run("--run_name", "../escape", "--force_overwrite")
    assert result.returncode != 0
    assert "must be a single directory name, not a path" in result.stdout + result.stderr
    assert sentinel.read_text() == "keep"
    assert not published.exists()


def test_pair_mode_routes_single_donor_recipient_without_cohort_normalization(workflow):
    run, published, folder = workflow
    result = run("--hla_typing", "HLA-A02:01,HLA-B07:02", mode="pair")
    assert result.returncode == 0, result.stdout + result.stderr
    count = json.loads((published / "logs" / "count.json").read_text())
    affinity = json.loads((published / "logs" / "affinity.json").read_text())
    assert count["pair"] == ""
    assert count["mode"] == "pair"
    assert (published / "plots" / "pie_chart.png").read_bytes() == b"pair plot placeholder"
    assert count["donor"] == "A" and count["recipient"] == "B"
    assert affinity == dict(count, hla="HLA-A02:01,HLA-B07:02", ensembl_exists=True)
    assert (published / "AMS" / "test" / "AMS.pkl").is_file()
    assert (published / "AAMS" / "test" / "AAMS_df.pkl").is_file()
    assert not (published / "AMS" / "test" / "AMS_df.tsv").exists()
    assert (published / "plots" / "distrib.png").is_file()
    with (folder / "trace.tsv").open() as handle:
        trace = list(csv.DictReader(handle, delimiter="\t"))
    assert len(trace) == 2
    assert all("ALLO_COUNT" in row["name"] or "ALLO_AFFINITY" in row["name"] for row in trace)
