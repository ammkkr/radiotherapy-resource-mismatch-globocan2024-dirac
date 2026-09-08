"""Regression checks for overlapping aggregates and absent site records."""
import unittest

from build_radiotherapy_resource_mismatch import build_burden_tables, RT_SITE_RULES, YEARS


def fixture():
    rows = []
    for year in YEARS:
        for measure in ("incidence", "mortality"):
            for code in [r["cancer_code"] for r in RT_SITE_RULES] + ["39", "40"]:
                count = "120" if code == "39" else "100" if code == "40" else "1"
                rows.append(dict(country_iso3="AAA", sex_label="both", year=year,
                                 measure=measure, cancer_code=code, predicted_count=count, pop="1000"))
    return rows


class BurdenIntegrityTests(unittest.TestCase):
    def test_aggregate_order_never_changes_total(self):
        for rows in (fixture(), list(reversed(fixture()))):
            data, _ = build_burden_tables(rows)
            self.assertEqual(data["AAA"]["all_cancer_cases_2024"], 100)
            self.assertEqual(data["AAA"]["selected_site_cases_2024"], 14)

    def test_missing_rectum_is_not_zero(self):
        with self.assertRaisesRegex(ValueError, "Incomplete burden"):
            build_burden_tables([r for r in fixture() if r["cancer_code"] != "9"])

    def test_missing_code40_cannot_fall_back(self):
        with self.assertRaisesRegex(ValueError, "Incomplete burden"):
            build_burden_tables([r for r in fixture() if r["cancer_code"] != "40"])

    def test_duplicate_site_is_rejected(self):
        rows = fixture()
        with self.assertRaisesRegex(ValueError, "Duplicate burden"):
            build_burden_tables(rows + [rows[0]])


if __name__ == "__main__":
    unittest.main()
