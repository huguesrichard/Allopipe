"""Download and validate the Ensembl references used by Allo-Affinity."""

import argparse
from datetime import datetime, timezone
import gzip
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
from urllib.request import urlopen


ENSEMBL_FTP = "https://ftp.ensembl.org/pub"


def reference_urls(assembly, release, base_url=ENSEMBL_FTP):
    root = f"{base_url}/grch37/current" if assembly == "GRCh37" else f"{base_url}/release-{release}"
    fasta = f"{root}/fasta/homo_sapiens"
    files = {
        f"Homo_sapiens.{assembly}.cdna.all.fa.gz": f"{fasta}/cdna",
        f"Homo_sapiens.{assembly}.pep.all.fa.gz": f"{fasta}/pep",
    }
    if assembly == "GRCh38":
        files[f"Homo_sapiens.GRCh38.{release}.refseq.tsv.gz"] = f"{root}/tsv/homo_sapiens"
    return files


def checksum_entry(directory_url, name):
    with urlopen(f"{directory_url}/CHECKSUMS", timeout=60) as response:
        lines = response.read().decode("utf-8").splitlines()
    for line in lines:
        fields = line.split()
        if len(fields) == 3 and fields[2] == name:
            return int(fields[0]), int(fields[1])
    raise ValueError(f"No checksum for {name} in {directory_url}/CHECKSUMS")


def local_checksum(path):
    # Ensembl CHECKSUMS use the BSD sum algorithm (the default for Unix sum).
    result = subprocess.run(["sum", str(path)], check=True, capture_output=True, text=True)
    fields = result.stdout.split()
    return int(fields[0]), int(fields[1])


def download_reference(directory_url, name, destination):
    expected = checksum_entry(directory_url, name)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=destination.parent, prefix=f".{name}.", suffix=".part", delete=False) as output:
            temporary = Path(output.name)
            with urlopen(f"{directory_url}/{name}", timeout=60) as response:
                while block := response.read(1024 * 1024):
                    output.write(block)
        actual = local_checksum(temporary)
        if actual != expected:
            raise ValueError(f"Checksum mismatch for {name}: expected {expected}, got {actual}")
        os.replace(temporary, destination)
        return expected
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def version_from_tag(tag):
    match = re.fullmatch(r"release_(\d+)(?:\.\d+)?", str(tag))
    return match.group(1) if match else None


def vep_version_from_vcf(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##VEP="):
                match = re.search(r'##VEP=["\']?v?(\d+)', line)
                return match.group(1) if match else None
            if line.startswith("#CHROM"):
                break
    return None


def prepare(assembly, release, directory, vep_version=None, vcfs=(), base_url=ENSEMBL_FTP):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    refseq = f"Homo_sapiens.{assembly}.refseq.tsv.gz"
    if assembly == "GRCh37" and not (directory / refseq).is_file():
        raise FileNotFoundError(f"The repository-provided GRCh37 RefSeq file is missing: {directory / refseq}")

    urls = reference_urls(assembly, release, base_url)

    manifest_path = directory / ".ensembl_downloads.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.is_file() else {}
    if assembly == "GRCh38":
        old_refseqs = [name for name in manifest if name.endswith(".refseq.tsv.gz")
                       and name not in urls and (directory / name).is_file()]
        if old_refseqs:
            raise ValueError(f"References from another Ensembl release are already downloaded in {directory}: {', '.join(old_refseqs)}. Remove the previous downloaded references before changing --ensembl_release.")
    for name, source in urls.items():
        target = directory / name
        previous = manifest.get(name)
        if target.is_file():
            if previous and previous.get("url") != f"{source}/{name}":
                raise ValueError(f"{name} was downloaded from {previous['url']}, which differs from the requested reference source {source}/{name}. Remove the previous downloaded references before changing --ensembl_release.")
            if previous and previous.get("checksum") and local_checksum(target) != tuple(previous["checksum"]):
                raise ValueError(f"Checksum mismatch for existing Ensembl reference: {target}")
            continue
        print(f"Downloading Ensembl reference: {source}/{name}", flush=True)
        checksum = download_reference(source, name, target)
        manifest[name] = {"url": f"{source}/{name}", "release": str(release) if assembly == "GRCh38" else "GRCh37 archive", "checksum": checksum, "downloaded_at": datetime.now(timezone.utc).isoformat()}
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    version_sources = {}
    if vep_version:
        version_sources[version_from_tag(vep_version)] = [f"--vep_version {vep_version}"]
    else:
        for path in vcfs:
            version_sources.setdefault(vep_version_from_vcf(path), []).append(f"##VEP header of {path}")
    for version in sorted(v for v in version_sources if v):
        if assembly == "GRCh38" and version != str(release):
            sources = ", ".join(version_sources[version])
            print(f"WARNING: VEP release {version} ({sources}) differs from Ensembl reference release {release} (selected release for GRCh38 reference downloads).", file=sys.stderr, flush=True)
    if not version_sources or None in version_sources:
        print("WARNING: VEP release could not be determined; reference version compatibility was not checked.", file=sys.stderr, flush=True)
    return directory


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assembly", required=True, choices=("GRCh37", "GRCh38"))
    parser.add_argument("--release", required=True, type=int)
    parser.add_argument("--directory", required=True)
    parser.add_argument("--vep-version")
    parser.add_argument("--vcf", action="append", default=[])
    args = parser.parse_args()
    try:
        prepare(args.assembly, args.release, args.directory, args.vep_version, args.vcf)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.exit(1, f"ERROR: Ensembl reference preparation failed: {error}\n")


if __name__ == "__main__":
    main()
