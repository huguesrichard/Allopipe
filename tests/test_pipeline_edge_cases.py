"""Boundary cases of active Nextflow-era parsing, cleavage and normalization."""
import gzip

import pandas as pd
import pytest

from tools import cleavage, normalize_cohort, parsing_functions


@pytest.fixture
def annotated_vcf(tmp_path):
    def write(compressed, frameshift_annotation):
        # Reordered annotation fields ensure indices come from the VEP header.
        fields = "Allele|Feature|Gene|Consequence|Protein_position|cDNA_position|CDS_position|Amino_acids|Codons|gnomADe_AF"
        if frameshift_annotation:
            fields += "|FrameshiftSequence"
        annotations = [
            "T|ENST1|ENSG1|missense_variant|42|123|45|K/N|aAa/aTa|0.001",
            "C|ENST2|ENSG2|frameshift_variant|55|165|60|G/R|gGg/cCc|0.002",
            "A|ENST3|ENSG3|synonymous_variant|81|243|90|K/K|aAa/aAa|0.003",
        ]
        if frameshift_annotation:
            annotations = [value + "|MVKKA" for value in annotations]
        expected = pd.DataFrame({
            "#CHROM": ["1", "X", "2"], "POS": ["100", "200", "300"],
            "ID": ["."] * 3, "REF": ["A", "G", "T"], "ALT": ["T", "C", "A"],
            "QUAL": ["."] * 3, "FILTER": ["PASS"] * 3,
            "INFO": ["CSQ=" + value for value in annotations],
            "FORMAT": ["GT:AD"] * 3, "sample.1": ["0/0:12,0", "0|1:8,7", "1/1:0,18"],
        })
        text = (
            "##fileformat=VCFv4.2\n"
            f'##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: {fields}">\n'
            + expected.to_csv(sep="\t", index=False)
        )
        path = tmp_path / ("input.vcf.gz" if compressed else "input.vcf")
        if compressed:
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write(text)
        else:
            path.write_text(text, encoding="utf-8")
        parser = parsing_functions.gzvcf_vep_parser if compressed else parsing_functions.vcf_vep_parser
        return path, parser, expected

    return write


@pytest.mark.parametrize("compressed", [False, True], ids=["vcf", "vcf-gz"])
@pytest.mark.parametrize("frameshift_mode", [False, True])
@pytest.mark.parametrize("chunksize,chunk_lengths", [(1, [1, 1, 1]), (2, [2, 1]), (10, [3])])
def test_chunked_parsing_matches_full_read_and_header_indices(
    annotated_vcf, compressed, frameshift_mode, chunksize, chunk_lengths,
):
    path, parser, expected = annotated_vcf(compressed, frameshift_mode)
    full, full_indices = parser(str(path), frameshift_mode=frameshift_mode)
    reader, indices = parser(str(path), frameshift_mode=frameshift_mode, chunksize=chunksize)
    with reader:
        chunks = list(reader)
    assert [len(chunk) for chunk in chunks] == chunk_lengths
    pd.testing.assert_frame_equal(pd.concat(chunks, ignore_index=True), expected)
    pd.testing.assert_frame_equal(full, expected)
    expected_indices = {
        "consequence": 3, "gene": 2, "transcript": 1, "prot": 4,
        "cdna": 5, "cds": 6, "aa": 7, "codons": 8, "gnomad": 9,
        "frameshift": 10 if frameshift_mode else None,
    }
    assert vars(indices) == vars(full_indices) == expected_indices


@pytest.mark.parametrize("compressed", [False, True], ids=["vcf", "vcf-gz"])
@pytest.mark.parametrize("chunksize", [None, 1])
def test_frameshift_requires_plugin_annotation_even_in_chunked_mode(annotated_vcf, compressed, chunksize):
    path, parser, _ = annotated_vcf(compressed, frameshift_annotation=False)
    with pytest.raises(SystemExit) as error:
        parser(str(path), frameshift_mode=True, chunksize=chunksize)
    assert error.value.code == (
        "Frameshift mode is enabled but 'FrameshiftSequence' annotation "
        "is missing from the VEP input file. "
        "Please re-run VEP with the Frameshift plugin enabled."
    )


@pytest.mark.parametrize("imputation,frameshift,error_text", [
    ("imputation", False, None),
    ("no-imputation", False, "Cleavage mode cannot be run with 'no-imputation' mode."),
    ("imputation", True, "Cleavage mode cannot be run with frameshift mode."),
    ("no-imputation", True, "Cleavage mode cannot be run with 'no-imputation' mode."),
])
def test_cleavage_validates_modes_recorded_in_ams_log(tmp_path, imputation, frameshift, error_text):
    log_file = tmp_path / "P01_run.log"
    content = f"Imputation: {imputation}\nFrameshift: {frameshift}\nOrientation: dr\n"
    log_file.write_text(content, encoding="utf-8")
    if error_text is None:
        assert cleavage.validate_cleavage_compatibility(log_file) is None
    else:
        with pytest.raises(ValueError) as error:
            cleavage.validate_cleavage_compatibility(log_file)
        assert str(error.value).startswith(error_text)
        assert "Please rerun Allo-Count" in str(error.value)
    assert log_file.read_text() == content


@pytest.mark.parametrize("empty", [False, True], ids=["no-reference-genotypes", "empty-tables"])
def test_reference_ratio_without_reference_sites_is_undefined(empty, capsys):
    donor = pd.DataFrame({"CHROM": ["1", "1"], "POS": [100, 200], "GT": ["0/1", "1/1"]})
    recipient = donor.assign(GT=["1/1", "0/1"])
    if empty:
        donor, recipient = donor.iloc[0:0], recipient.iloc[0:0]
    assert normalize_cohort.get_ref_ratio_pair(donor, recipient) == (0, 0, None)
    assert capsys.readouterr().out == "Warning: Normalization cannot be computed - division by 0\n"


@pytest.mark.parametrize("mixed_cohort", [False, True], ids=["only-undefined-pair", "mixed-cohort"])
def test_undefined_ratio_disables_normalization_and_serializes_na(tmp_path, mixed_cohort, capsys):
    ams_dir = tmp_path / "AMS"
    ams_dir.mkdir()
    ams_files, donors, recipients = [], [], []
    pairs = [("P1", 5, False)]
    if mixed_cohort:
        pairs.append(("P2", 9, True))
    for pair, score, has_reference in pairs:
        ams_file = ams_dir / f"{pair}_AMS.pkl"
        pd.DataFrame({"pair": [pair], "ams": [score]}).to_pickle(ams_file)
        donor = pd.DataFrame({
            "CHROM": ["1", "1"], "POS": [100, 200],
            "GT": ["0/0", "0/0"] if has_reference else ["0/1", "1/1"],
        })
        recipient = donor.assign(GT=["0/0", "0/1"] if has_reference else ["1/1", "0/1"])
        donor_file = tmp_path / f"{pair}_D0_table.tsv"
        recipient_file = tmp_path / f"{pair}_R0_table.tsv"
        donor.to_csv(donor_file, sep="\t", index=False)
        recipient.to_csv(recipient_file, sep="\t", index=False)
        ams_files.append(str(ams_file))
        donors.append(str(donor_file))
        recipients.append(str(recipient_file))
    expected = pd.DataFrame({
        "pair": [pair for pair, _, _ in pairs], "ams_giab": [score for _, score, _ in pairs],
        "common_ref": [0, 1] if mixed_cohort else [0],
        "total_ref": [0, 2] if mixed_cohort else [0],
        "ref_ratio": ["NA"] * len(pairs), "ams_norm": ["NA"] * len(pairs),
    })
    result = normalize_cohort.normalize_ams(list(reversed(ams_files)), donors, recipients)
    pd.testing.assert_frame_equal(result.reset_index(drop=True), expected)
    assert capsys.readouterr().out == "Warning: Normalization cannot be computed - division by 0\n"
    normalize_cohort.write_normalized_ams_tables(ams_files, donors, recipients)
    # Disable pandas' NA interpretation to verify literal published values.
    saved = pd.read_csv(ams_dir / "AMS_df.tsv", sep="\t", keep_default_na=False)
    pd.testing.assert_frame_equal(saved, expected)
    assert capsys.readouterr().out == "Warning: Normalization cannot be computed - division by 0\n"
