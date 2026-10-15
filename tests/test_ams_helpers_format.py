import sys
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from tools import ams_helpers


class TestParseFormat(unittest.TestCase):
    def test_returns_required_indices_in_expected_order(self):
        self.assertEqual(
            ams_helpers.parse_format(["AD", "DP", "GT", "GQ"]),
            [2, 3, 0, 1],
        )

    def test_includes_optional_ft(self):
        self.assertEqual(
            ams_helpers.parse_format(["GT", "FT", "DP", "AD", "GQ"]),
            [0, 4, 3, 2, 1],
        )

    def test_reports_all_missing_required_fields(self):
        with self.assertRaisesRegex(
            ValueError,
            r"^Invalid VCF FORMAT: missing required field\(s\): GQ, DP\. "
            r"Observed FORMAT: GT:AD\.$",
        ):
            ams_helpers.parse_format(["GT", "AD"])

    def test_reports_file_sample_and_record_context(self):
        dataframe = pd.DataFrame(
            {
                "#CHROM": ["chr1"],
                "POS": ["12345"],
                "REF": ["A"],
                "ALT": ["T"],
                "FORMAT": ["GT:GQ:DP"],
                "donor": ["0/1:99:30"],
            }
        )

        with self.assertRaisesRegex(
            ValueError,
            r"^Invalid VCF FORMAT in file 'donor\.vcf', sample 'donor', "
            r"record 1:12345: missing required field\(s\): AD\. "
            r"Observed FORMAT: GT:GQ:DP\.$",
        ):
            ams_helpers.get_read_counts(
                dataframe,
                "/input/donor.vcf",
                min_ad=5,
                min_gq=0,
                base_length=3,
            )


if __name__ == "__main__":
    unittest.main()
