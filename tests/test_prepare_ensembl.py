import gzip
import json
from pathlib import Path

import pytest

from tools import prepare_ensembl


def make_remote(base, assembly="GRCh38", release=116):
    urls = prepare_ensembl.reference_urls(assembly, release, base.as_uri())
    for name, directory_url in urls.items():
        directory = Path(directory_url.removeprefix("file://"))
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / name
        path.write_bytes(gzip.compress(f"fixture {name}\n".encode(), mtime=0))
        checksum, blocks = prepare_ensembl.local_checksum(path)
        (directory / "CHECKSUMS").write_text(f"{checksum:05} {blocks} {name}\n", encoding="utf-8")
    return urls


def test_downloads_validates_and_reuses_without_network(tmp_path, capsys):
    remote = tmp_path / "remote"
    remote.mkdir()
    urls = make_remote(remote)
    target = tmp_path / "GRCh38"

    prepare_ensembl.prepare("GRCh38", 116, target, vep_version="release_116.2", base_url=remote.as_uri())
    assert all((target / name).is_file() for name in urls)
    manifest = json.loads((target / ".ensembl_downloads.json").read_text())
    assert len(manifest) == 3
    for entry in manifest.values():
        timestamp = prepare_ensembl.datetime.fromisoformat(entry["downloaded_at"])
        assert timestamp.utcoffset() == prepare_ensembl.timezone.utc.utcoffset(timestamp)
    first_run = capsys.readouterr()
    assert first_run.out.count("Downloading Ensembl reference:") == 3
    assert first_run.err == ""
    for path in remote.rglob("CHECKSUMS"):
        path.unlink()

    prepare_ensembl.prepare("GRCh38", 116, target, vep_version="release_116.2", base_url=remote.as_uri())
    second_run = capsys.readouterr()
    assert second_run.out == ""
    assert second_run.err == ""


def test_wrong_checksum_leaves_no_downloaded_file(tmp_path):
    remote = tmp_path / "remote"
    remote.mkdir()
    urls = make_remote(remote)
    name, directory_url = next(iter(urls.items()))
    Path(directory_url.removeprefix("file://"), "CHECKSUMS").write_text(f"1 1 {name}\n")
    target = tmp_path / "GRCh38"

    with pytest.raises(ValueError, match="Checksum mismatch"):
        prepare_ensembl.prepare("GRCh38", 116, target, vep_version="release_116.2", base_url=remote.as_uri())
    assert not (target / name).exists()
    assert not list(target.glob("*.part"))


def test_detects_corruption_of_previously_downloaded_file(tmp_path):
    remote = tmp_path / "remote"
    remote.mkdir()
    urls = make_remote(remote)
    target = tmp_path / "GRCh38"
    prepare_ensembl.prepare("GRCh38", 116, target, vep_version="release_116.2", base_url=remote.as_uri())
    (target / next(iter(urls))).write_bytes(b"corrupted")
    with pytest.raises(ValueError, match="Checksum mismatch for existing"):
        prepare_ensembl.prepare("GRCh38", 116, target, vep_version="release_116.2", base_url=remote.as_uri())


def test_prevents_mixing_downloaded_releases(tmp_path):
    remote = tmp_path / "remote"
    remote.mkdir()
    make_remote(remote, release=111)
    target = tmp_path / "GRCh38"
    prepare_ensembl.prepare("GRCh38", 111, target, vep_version="release_111", base_url=remote.as_uri())
    with pytest.raises(ValueError, match="another Ensembl release"):
        prepare_ensembl.prepare("GRCh38", 116, target, vep_version="release_116", base_url=remote.as_uri())
    assert not (target / "Homo_sapiens.GRCh38.116.refseq.tsv.gz").exists()


def test_grch37_preserves_provided_refseq_and_downloads_only_fastas(tmp_path):
    remote = tmp_path / "remote"
    remote.mkdir()
    urls = make_remote(remote, "GRCh37")
    target = tmp_path / "GRCh37"
    target.mkdir()
    refseq = target / "Homo_sapiens.GRCh37.refseq.tsv.gz"
    refseq.write_bytes(b"provided")
    prepare_ensembl.prepare("GRCh37", 116, target, vep_version="release_116.2", base_url=remote.as_uri())
    assert refseq.read_bytes() == b"provided"
    assert len(urls) == 2
    assert all((target / name).is_file() for name in urls)


def test_grch37_requires_the_repo_provided_refseq(tmp_path):
    with pytest.raises(FileNotFoundError, match="repository-provided GRCh37 RefSeq"):
        prepare_ensembl.prepare("GRCh37", 116, tmp_path / "GRCh37")


def test_warns_for_preannotated_vcf_with_different_vep_release(tmp_path, capsys):
    remote = tmp_path / "remote"
    remote.mkdir()
    make_remote(remote)
    vcf = tmp_path / "test.vcf.gz"
    with gzip.open(vcf, "wt") as handle:
        handle.write('##fileformat=VCFv4.2\n##VEP="v113.0"\n#CHROM\tPOS\n')
    prepare_ensembl.prepare("GRCh38", 116, tmp_path / "GRCh38", vcfs=[vcf], base_url=remote.as_uri())
    warning = capsys.readouterr().err
    assert f"VEP release 113 (##VEP header of {vcf})" in warning
    assert "Ensembl reference release 116 (selected release for GRCh38 reference downloads)" in warning


def test_warning_identifies_vep_version_argument(tmp_path, capsys):
    remote = tmp_path / "remote"
    remote.mkdir()
    make_remote(remote)
    prepare_ensembl.prepare("GRCh38", 116, tmp_path / "GRCh38", vep_version="release_113.0", base_url=remote.as_uri())
    assert "VEP release 113 (--vep_version release_113.0)" in capsys.readouterr().err


def test_warns_if_vep_release_cannot_be_identified(tmp_path, capsys):
    remote = tmp_path / "remote"
    remote.mkdir()
    make_remote(remote)
    vcf = tmp_path / "test.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\n", encoding="utf-8")
    prepare_ensembl.prepare("GRCh38", 116, tmp_path / "GRCh38", vcfs=[vcf], base_url=remote.as_uri())
    assert "VEP release could not be determined" in capsys.readouterr().err


def test_custom_vep_tag_does_not_guess_release_from_unrelated_digits(tmp_path, capsys):
    remote = tmp_path / "remote"
    remote.mkdir()
    make_remote(remote)
    prepare_ensembl.prepare("GRCh38", 116, tmp_path / "GRCh38", vep_version="custom-2024", base_url=remote.as_uri())
    assert "VEP release could not be determined" in capsys.readouterr().err
