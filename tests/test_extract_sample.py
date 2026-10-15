"""Integration tests of the actual EXTRACT_SAMPLE shell script, not Nextflow orchestration."""
import gzip
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

import pytest


SAMPLES = ["sample.1", "sample.2"]
HEADER = "\n".join([
    "##fileformat=VCFv4.2",
    "##contig=<ID=1,length=1000>",
    '##INFO=<ID=CSQ,Number=.,Type=String,Description="VEP annotations">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(SAMPLES),
]) + "\n"
RECORDS = [
    "1\t100\t.\tA\tT\t.\tPASS\tCSQ=T|missense_variant|ENST1\tGT:AD:DP\t0/0:12,0:12\t0|1:8,7:15",
    "1\t200\t.\tG\tC\t.\tPASS\tCSQ=C|missense_variant|ENST2\tGT:AD:DP\t0/1:10,10:20\t./.:.,.:0",
    "1\t300\t.\tT\tA\t.\tPASS\tCSQ=A|missense_variant|ENST3\tGT:AD:DP\t1/1:0,18:18\t0/0:22,0:22",
    "1\t400\t.\tC\tG\t.\tPASS\tCSQ=G|missense_variant|ENST4\tGT:AD:DP\t./.:.,.:0\t1/1:0,16:16",
    "1\t500\t.\tA\tC,G\t.\tPASS\tCSQ=C|missense_variant|ENST5,G|stop_gained|ENST6\tGT:AD:DP\t1/2:0,9,9:18\t0/1:10,10,0:20",
    "1\t600\t.\tG\tT\t.\tPASS\tCSQ=T|missense_variant|ENST7\tGT:AD:DP\t.|.:.,.:0\t0/0:11,0:11",
]
QUERY_FORMAT = "%CHROM\t%POS\t%REF\t%ALT\t%INFO/CSQ[\t%GT\t%AD\t%DP]\n"


@pytest.fixture
def bcftools():
    # Prefer the executable bundled with the Python environment running pytest.
    executable = Path(sys.executable).parent / "bcftools"
    if not executable.is_file():
        found = shutil.which("bcftools")
        if not found:
            pytest.skip("EXTRACT_SAMPLE integration tests require a real bcftools executable")
        executable = Path(found)
    return str(executable)


@pytest.fixture
def extraction_script():
    module = Path(__file__).resolve().parents[1] / "modules" / "extract-sample.nf"
    source = module.read_text(encoding="utf-8")
    match = re.search(r'script:\s*"""(.*?)"""', source, flags=re.DOTALL)
    assert match is not None, "EXTRACT_SAMPLE must expose its shell script"
    return match.group(1)


def run_bcftools(executable, *args):
    result = subprocess.run([executable, *args], capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
    return result.stdout


@pytest.fixture
def extract(tmp_path, bcftools, extraction_script):
    def run(sample, compressed=False, records=RECORDS):
        input_dir = tmp_path / "input with spaces"
        input_dir.mkdir()
        source = input_dir / "cohort.vcf"
        source.write_text(HEADER + "\n".join(records) + "\n", encoding="utf-8")
        if compressed:
            compressed_source = input_dir / "cohort.vcf.gz"
            run_bcftools(bcftools, "view", "-O", "z", "-o", str(compressed_source), str(source))
            source = compressed_source
        original = source.read_bytes()
        worker = tmp_path / "worker"
        worker.mkdir()
        # Nextflow stages the input under its basename in the task directory.
        staged_input = worker / source.name
        shutil.copyfile(source, staged_input)
        script = extraction_script.replace("${sample_id}", sample)
        script = script.replace("${multi_vcf}", staged_input.name)
        assert "${" not in script, "Unexpected module placeholders require explicit test support"
        env = dict(os.environ, PATH=str(Path(bcftools).parent) + os.pathsep + os.environ.get("PATH", ""))
        # Nextflow task scripts run with errexit, nounset and pipefail.
        result = subprocess.run(
            ["bash", "-euo", "pipefail", "-c", script], cwd=worker,
            env=env, capture_output=True, text=True, timeout=30,
        )
        assert source.read_bytes() == original
        assert staged_input.read_bytes() == original
        return result, worker / f"{sample}.vcf.gz"

    return run


def expected_records(sample, positions):
    sample_index = 9 + SAMPLES.index(sample)
    expected = []
    for record in RECORDS:
        columns = record.split("\t")
        if int(columns[1]) in positions:
            expected.append("\t".join([
                columns[0], columns[1], columns[3], columns[4],
                columns[7].removeprefix("CSQ="), *columns[sample_index].split(":"),
            ]))
    return "\n".join(expected) + ("\n" if expected else "")


@pytest.mark.parametrize("compressed", [False, True], ids=["vcf", "vcf-gz"])
@pytest.mark.parametrize("sample,positions", [
    ("sample.1", [100, 200, 300, 500]),
    ("sample.2", [100, 300, 400, 500, 600]),
])
def test_extracts_only_requested_sample_and_preserves_records(extract, bcftools, sample, positions, compressed):
    result, output = extract(sample, compressed=compressed)
    assert result.returncode == 0, result.stdout + result.stderr
    assert run_bcftools(bcftools, "query", "-l", str(output)).splitlines() == [sample]
    # Exact REF/ALT, CSQ and GT/AD/DP: includes 0/0, phased GT and multi-allelic GT.
    observed = run_bcftools(bcftools, "query", "-f", QUERY_FORMAT, str(output))
    assert observed == expected_records(sample, positions)
    index = Path(str(output) + ".tbi")
    assert index.is_file() and index.stat().st_size > 0
    assert run_bcftools(bcftools, "index", "--nrecords", str(output)).strip() == str(len(positions))
    regional = run_bcftools(bcftools, "query", "-r", "1:200-500", "-f", QUERY_FORMAT, str(output))
    assert regional == expected_records(sample, [pos for pos in positions if 200 <= pos <= 500])
    with gzip.open(output, "rt", encoding="utf-8") as handle:
        assert "#CHROM\t" in handle.read()
    input_name = "cohort.vcf.gz" if compressed else "cohort.vcf"
    assert {path.name for path in output.parent.iterdir()} == {input_name, output.name, index.name}


@pytest.mark.parametrize("sample,record_index", [("sample.1", 3), ("sample.2", 1)])
def test_all_missing_genotypes_produce_valid_empty_indexed_vcf(extract, bcftools, sample, record_index):
    result, output = extract(sample, records=[RECORDS[record_index]])
    assert result.returncode == 0, result.stdout + result.stderr
    assert run_bcftools(bcftools, "query", "-l", str(output)).splitlines() == [sample]
    assert run_bcftools(bcftools, "query", "-f", QUERY_FORMAT, str(output)) == ""
    assert Path(str(output) + ".tbi").is_file()
    assert run_bcftools(bcftools, "query", "-r", "1:1-1000", "-f", QUERY_FORMAT, str(output)) == ""


def test_absent_sample_fails_and_does_not_create_an_index(extract):
    result, output = extract("absent_sample")
    assert result.returncode != 0
    assert "absent_sample" in result.stderr
    assert "does not exist" in result.stderr
    assert not Path(str(output) + ".tbi").exists()
