"""Build reviewer-facing additional analyses for the Radiology Research Letter.

The analyses isolate the GCO burden-version update while holding the latest
DIRAC resource snapshot fixed, decouple the growth and density denominators,
quantify threshold stability, and partially identify regional proportions when
DIRAC resource observations are missing.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CURRENT_DATASET = PROJECT_ROOT / "data" / "radiotherapy_resource_mismatch_country_2024_2050.csv"
GCO_2022_PROJECTIONS = PROJECT_ROOT / "data" / "source" / "gco_cancer_tomorrow_2022_predictions_long.csv"
EXISTING_ROBUSTNESS = PROJECT_ROOT / "data" / "table_14_radiology_research_letter_robustness.csv"
OUT_DATASET = PROJECT_ROOT / "data" / "radiotherapy_resource_mismatch_country_v6_sensitivity.csv"
TABLE_DIR = PROJECT_ROOT / "data"
REPORT_DIR = PROJECT_ROOT / "docs"
LOG_PATH = PROJECT_ROOT / "data" / "additional_analyses_log.json"

SELECTED_SITE_CODES = {1, 3, 4, 5, 6, 9, 14, 15, 20, 23, 24, 27, 30, 31}
SUBSTANTIVE_SCENARIOS = [
    "All-cancer burden",
    "2040 horizon",
    "RT centres metric",
    "Mortality growth",
    "Excluding small baseline burden",
]
DISPLAY_COUNTRIES = {
    "CIV": "Cote d'Ivoire",
    "COD": "Democratic Republic of the Congo",
    "PRK": "Democratic People's Republic of Korea",
    "TZA": "Tanzania",
}


def flag(series: pd.Series) -> pd.Series:
    return series.astype(str).isin({"1", "1.0", "True", "true"})


def quantile(series: pd.Series, q: float) -> float:
    return float(pd.to_numeric(series, errors="coerce").dropna().quantile(q, interpolation="linear"))


def iso_set(value: Any) -> set[str]:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return set()
    return {item.strip() for item in str(value).split(";") if item.strip()}


def frame_to_markdown(frame: pd.DataFrame, float_digits: int = 3) -> str:
    def render(value: Any) -> str:
        if pd.isna(value):
            return ""
        if isinstance(value, (float, np.floating)):
            return f"{float(value):.{float_digits}f}"
        return str(value).replace("|", "\\|")

    columns = [str(column) for column in frame.columns]
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in frame.itertuples(index=False, name=None):
        lines.append("| " + " | ".join(render(value) for value in row) + " |")
    return "\n".join(lines)


def build_gco_2022_burden() -> pd.DataFrame:
    usecols = [
        "country_iso3",
        "population_label",
        "sex_label",
        "measure",
        "cancer_code",
        "year",
        "predicted_count",
        "pop",
    ]
    gco = pd.read_csv(GCO_2022_PROJECTIONS, usecols=usecols)
    gco["cancer_code"] = pd.to_numeric(gco["cancer_code"], errors="coerce")
    gco["year"] = pd.to_numeric(gco["year"], errors="coerce")
    gco["predicted_count"] = pd.to_numeric(gco["predicted_count"], errors="coerce")
    gco["pop"] = pd.to_numeric(gco["pop"], errors="coerce")

    selected = gco[
        gco["country_iso3"].notna()
        & gco["sex_label"].eq("both")
        & gco["measure"].eq("incidence")
        & gco["cancer_code"].isin(SELECTED_SITE_CODES)
        & gco["year"].isin([2022, 2040, 2050])
    ].copy()
    burden = (
        selected.groupby(["country_iso3", "year"], as_index=False)["predicted_count"]
        .sum()
        .pivot(index="country_iso3", columns="year", values="predicted_count")
        .rename(columns={2022: "gco2022_selected_cases_2022", 2040: "gco2022_selected_cases_2040", 2050: "gco2022_selected_cases_2050"})
        .reset_index()
    )
    burden["gco2022_selected_growth_2022_2050"] = (
        burden["gco2022_selected_cases_2050"] - burden["gco2022_selected_cases_2022"]
    ) / burden["gco2022_selected_cases_2022"]
    return burden


def add_version_comparison(df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    burden_2022 = build_gco_2022_burden()
    out = df.merge(burden_2022, on="country_iso3", how="left", validate="one_to_one")
    matched_common = flag(out["dirac_matched"]) & out["mv_therapy_units"].notna() & out["gco2022_selected_growth_2022_2050"].notna()
    out["gco2022_mv_per_1000_selected_cases_2050"] = np.where(
        matched_common,
        out["mv_therapy_units"] * 1000 / out["gco2022_selected_cases_2050"],
        np.nan,
    )
    growth_q75 = quantile(out.loc[matched_common, "gco2022_selected_growth_2022_2050"], 0.75)
    density_q25 = quantile(out.loc[matched_common, "gco2022_mv_per_1000_selected_cases_2050"], 0.25)
    out["gco2022_high_growth_q75"] = out["gco2022_selected_growth_2022_2050"].ge(growth_q75) & matched_common
    out["gco2022_lower_density_q25"] = out["gco2022_mv_per_1000_selected_cases_2050"].le(density_q25) & matched_common
    out["gco2022_screen_positive"] = out["gco2022_high_growth_q75"] & out["gco2022_lower_density_q25"]
    out["gco2022_growth_percentile"] = np.nan
    out["gco2024_growth_percentile_common"] = np.nan
    out.loc[matched_common, "gco2022_growth_percentile"] = (
        out.loc[matched_common, "gco2022_selected_growth_2022_2050"].rank(method="average", pct=True) * 100
    )
    out.loc[matched_common, "gco2024_growth_percentile_common"] = (
        out.loc[matched_common, "selected_site_relative_case_growth_2050"].rank(method="average", pct=True) * 100
    )

    out["current_screen_positive"] = flag(out["meets_both_thresholds"])
    out["burden_version_transition"] = "Not screen-positive in either version"
    out.loc[matched_common & out["gco2022_screen_positive"] & out["current_screen_positive"], "burden_version_transition"] = "Retained"
    out.loc[matched_common & ~out["gco2022_screen_positive"] & out["current_screen_positive"], "burden_version_transition"] = "New in GCO version 2024"
    out.loc[matched_common & out["gco2022_screen_positive"] & ~out["current_screen_positive"], "burden_version_transition"] = "No longer screen-positive"
    out.loc[~matched_common, "burden_version_transition"] = "Not comparable"

    old_set = set(out.loc[out["gco2022_screen_positive"], "country_iso3"])
    current_common_set = set(out.loc[matched_common & out["current_screen_positive"], "country_iso3"])
    metrics = {
        "common_matched_countries": int(matched_common.sum()),
        "gco2022_growth_q75": growth_q75,
        "gco2022_density_q25": density_q25,
        "gco2022_screen_positive_count": len(old_set),
        "gco2024_screen_positive_count_common": len(current_common_set),
        "retained_count": len(old_set & current_common_set),
        "new_in_2024": sorted(current_common_set - old_set),
        "no_longer_in_2024": sorted(old_set - current_common_set),
        "jaccard": len(old_set & current_common_set) / len(old_set | current_common_set),
    }
    return out, metrics


def add_decoupled_and_threshold_analyses(df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any], pd.DataFrame]:
    out = df.copy()
    matched = flag(out["dirac_matched"]) & out["mv_therapy_units"].notna()
    growth = pd.to_numeric(out["selected_site_relative_case_growth_2050"], errors="coerce")
    cases_2024 = pd.to_numeric(out["selected_site_cases_2024"], errors="coerce")
    density_2050 = pd.to_numeric(out["mv_units_per_1000_selected_site_cases_2050"], errors="coerce")
    out["mv_units_per_1000_selected_cases_2024"] = np.where(matched, out["mv_therapy_units"] * 1000 / cases_2024, np.nan)

    growth_q75 = quantile(growth[matched], 0.75)
    density_2050_q25 = quantile(density_2050[matched], 0.25)
    density_2024_q25 = quantile(out.loc[matched, "mv_units_per_1000_selected_cases_2024"], 0.25)
    out["screen_positive_2024_density"] = (
        matched & growth.ge(growth_q75) & out["mv_units_per_1000_selected_cases_2024"].le(density_2024_q25)
    )

    growth_q80 = quantile(growth[matched], 0.80)
    density_q20 = quantile(density_2050[matched], 0.20)
    out["screen_positive_q80_q20"] = matched & growth.ge(growth_q80) & density_2050.le(density_q20)

    current_set = set(out.loc[out["current_screen_positive"], "country_iso3"])
    baseline_set = set(out.loc[out["screen_positive_2024_density"], "country_iso3"])
    strict_set = set(out.loc[out["screen_positive_q80_q20"], "country_iso3"])
    metrics = {
        "growth_q75": growth_q75,
        "density_2050_q25": density_2050_q25,
        "density_2024_q25": density_2024_q25,
        "baseline_density_screen_positive_count": len(baseline_set),
        "baseline_density_primary_retained": len(current_set & baseline_set),
        "baseline_density_not_retained": sorted(current_set - baseline_set),
        "growth_q80": growth_q80,
        "density_q20": density_q20,
        "strict_screen_positive_count": len(strict_set),
        "strict_primary_retained": len(current_set & strict_set),
        "strict_not_retained": sorted(current_set - strict_set),
    }

    grid_rows: list[dict[str, Any]] = []
    for growth_quantile, density_quantile in [(0.67, 0.33), (0.70, 0.30), (0.75, 0.25), (0.80, 0.20)]:
        growth_threshold = quantile(growth[matched], growth_quantile)
        density_threshold = quantile(density_2050[matched], density_quantile)
        scenario_set = set(
            out.loc[matched & growth.ge(growth_threshold) & density_2050.le(density_threshold), "country_iso3"]
        )
        grid_rows.append(
            {
                "Analysis": f"Q{int(growth_quantile * 100)} growth + Q{int(density_quantile * 100)} density",
                "Growth threshold, %": 100 * growth_threshold,
                "Density threshold": density_threshold,
                "Countries screen-positive, No.": len(scenario_set),
                "Primary countries retained, No./23": f"{len(current_set & scenario_set)}/23",
                "Screen-positive ISO3": ";".join(sorted(scenario_set)),
            }
        )
    return out, metrics, pd.DataFrame(grid_rows)


def add_sensitivity_support(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    robustness = pd.read_csv(EXISTING_ROBUSTNESS)
    scenario_sets = {
        row["Analysis"]: iso_set(row["Flagged ISO3"])
        for _, row in robustness.iterrows()
        if row["Analysis"] in SUBSTANTIVE_SCENARIOS
    }
    out["support_baseline_density"] = out["screen_positive_2024_density"]
    out["support_strict_q80_q20"] = out["screen_positive_q80_q20"]
    for analysis in SUBSTANTIVE_SCENARIOS:
        column = "support_" + analysis.lower().replace("-", "_").replace(" ", "_")
        out[column] = out["country_iso3"].isin(scenario_sets[analysis])
    support_columns = ["support_baseline_density", "support_strict_q80_q20"] + [
        "support_" + analysis.lower().replace("-", "_").replace(" ", "_")
        for analysis in SUBSTANTIVE_SCENARIOS
    ]
    out["sensitivity_support_count_7"] = out[support_columns].sum(axis=1).astype(int)
    return out


def build_missingness_bounds(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []

    def add_row(stratum: str, category: str, sub: pd.DataFrame) -> None:
        n_all = len(sub)
        matched = flag(sub["dirac_matched"])
        observed = int(sub["current_screen_positive"].sum())
        high_growth_missing = int((~matched & flag(sub["high_growth_q75"])).sum())
        rows.append(
            {
                "Stratum": stratum,
                "Category": category,
                "All GCO country records, No.": n_all,
                "DIRAC matched, No.": int(matched.sum()),
                "DIRAC missing, No.": int((~matched).sum()),
                "Observed screen-positive, No.": observed,
                "High-growth resource-unknown, No.": high_growth_missing,
                "Lower-bound proportion, %": 100 * observed / n_all if n_all else np.nan,
                "Upper-bound proportion, %": 100 * (observed + high_growth_missing) / n_all if n_all else np.nan,
                "Matched-only proportion, %": 100 * observed / matched.sum() if matched.sum() else np.nan,
            }
        )

    add_row("Global", "All countries", df)
    for column, stratum in [("who_region", "WHO region"), ("income_label", "Income group"), ("hdi_label", "HDI group")]:
        labels = df[column].fillna("").replace("", "Missing")
        for category in sorted(labels.unique()):
            add_row(stratum, category, df[labels.eq(category)])
    return pd.DataFrame(rows)


def build_construct_diagnostics(df: pd.DataFrame) -> pd.DataFrame:
    matched = df[flag(df["dirac_matched"])].copy()
    matched["population_growth_2024_2050"] = (
        matched["population_2050"] - matched["population_2024"]
    ) / matched["population_2024"]
    common = matched[matched["gco2022_selected_growth_2022_2050"].notna()].copy()

    pairs = [
        (
            "Selected-site growth vs all-cancer growth, GCO version 2024",
            matched["selected_site_relative_case_growth_2050"],
            matched["all_cancer_relative_case_growth_2050"],
        ),
        (
            "Selected-site growth vs population growth, GCO version 2024",
            matched["selected_site_relative_case_growth_2050"],
            matched["population_growth_2024_2050"],
        ),
        (
            "Selected-site growth vs 2050-denominator MV-unit density",
            matched["selected_site_relative_case_growth_2050"],
            matched["mv_units_per_1000_selected_site_cases_2050"],
        ),
        (
            "Selected-site growth vs 2024-denominator MV-unit density",
            matched["selected_site_relative_case_growth_2050"],
            matched["mv_units_per_1000_selected_cases_2024"],
        ),
        (
            "Selected-site growth, GCO version 2022 vs version 2024",
            common["gco2022_selected_growth_2022_2050"],
            common["selected_site_relative_case_growth_2050"],
        ),
        (
            "2050 MV-unit density, GCO version 2022 vs version 2024",
            common["gco2022_mv_per_1000_selected_cases_2050"],
            common["mv_units_per_1000_selected_site_cases_2050"],
        ),
    ]
    rows = []
    for analysis, x, y in pairs:
        valid = pd.to_numeric(x, errors="coerce").notna() & pd.to_numeric(y, errors="coerce").notna()
        rho = float(spearmanr(pd.to_numeric(x[valid]), pd.to_numeric(y[valid])).statistic)
        rows.append({"Diagnostic": analysis, "Countries, No.": int(valid.sum()), "Spearman rho": rho})
    return pd.DataFrame(rows)


def build_version_table(df: pd.DataFrame) -> pd.DataFrame:
    common = df[df["burden_version_transition"].ne("Not comparable")].copy()
    priority = {
        "New in GCO version 2024": 0,
        "No longer screen-positive": 1,
        "Retained": 2,
        "Not screen-positive in either version": 3,
    }
    common["transition_order"] = common["burden_version_transition"].map(priority)
    columns = [
        "country_iso3",
        "country_display",
        "who_region",
        "income_label",
        "gco2022_selected_cases_2022",
        "gco2022_selected_cases_2050",
        "gco2022_selected_growth_2022_2050",
        "gco2022_growth_percentile",
        "gco2022_mv_per_1000_selected_cases_2050",
        "selected_site_cases_2024",
        "selected_site_cases_2050",
        "selected_site_relative_case_growth_2050",
        "gco2024_growth_percentile_common",
        "mv_units_per_1000_selected_site_cases_2050",
        "gco2022_screen_positive",
        "current_screen_positive",
        "burden_version_transition",
    ]
    return common.sort_values(["transition_order", "country_iso3"])[columns].rename(
        columns={"country_display": "country"}
    )


def build_country_stability_table(df: pd.DataFrame) -> pd.DataFrame:
    primary = df[df["current_screen_positive"]].copy()
    primary["Region / income"] = primary["who_region"].fillna("") + " / " + primary["income_label"].fillna("")
    primary["Selected-site cases, 2024 to 2050 (growth)"] = primary.apply(
        lambda row: f"{int(row['selected_site_cases_2024']):,} to {int(row['selected_site_cases_2050']):,} (+{100 * row['selected_site_relative_case_growth_2050']:.0f}%)",
        axis=1,
    )
    primary["MV units (record year)"] = primary.apply(
        lambda row: f"{int(row['mv_therapy_units'])} ({int(row['dirac_last_update_year'])})",
        axis=1,
    )
    primary["2050 cases/MV unit"] = (primary["selected_site_cases_2050"] / primary["mv_therapy_units"]).round().astype(int)
    primary["GCO 2022 burden-version status"] = np.where(
        primary["gco2022_screen_positive"], "Retained", "New in version 2024"
    )
    primary["2024-density sensitivity"] = np.where(
        primary["screen_positive_2024_density"], "Retained", "Not retained"
    )
    primary["Sensitivity support, No./7"] = primary["sensitivity_support_count_7"].astype(str) + "/7"
    primary["sort_pressure"] = primary["selected_site_cases_2050"] / primary["mv_therapy_units"]
    return primary.sort_values("sort_pressure", ascending=False)[
        [
            "country_display",
            "Region / income",
            "Selected-site cases, 2024 to 2050 (growth)",
            "MV units (record year)",
            "2050 cases/MV unit",
            "GCO 2022 burden-version status",
            "2024-density sensitivity",
            "Sensitivity support, No./7",
        ]
    ].rename(columns={"country_display": "Country"})


def write_report(metrics: dict[str, Any], diagnostics: pd.DataFrame, bounds: pd.DataFrame) -> None:
    global_row = bounds[(bounds["Stratum"] == "Global")].iloc[0]
    afro = bounds[(bounds["Stratum"] == "WHO region") & (bounds["Category"] == "AFRO")].iloc[0]
    lines = [
        "# Radiology v6 Additional Analyses",
        "",
        "## Burden-Version Update",
        "",
        f"- Common matched countries: {metrics['version']['common_matched_countries']}.",
        f"- GCO version 2022 screen-positive countries: {metrics['version']['gco2022_screen_positive_count']}.",
        f"- GCO version 2024 screen-positive countries in the common set: {metrics['version']['gco2024_screen_positive_count_common']}.",
        f"- Retained across versions: {metrics['version']['retained_count']}.",
        f"- New in version 2024: {', '.join(metrics['version']['new_in_2024'])}.",
        f"- No longer screen-positive in version 2024: {', '.join(metrics['version']['no_longer_in_2024'])}.",
        "- The same latest-reported DIRAC resource snapshot was used for both burden versions; this comparison isolates burden-data revision rather than resource change.",
        "",
        "## Denominator and Threshold Sensitivity",
        "",
        f"- Using 2024 rather than 2050 selected-site cases in the MV-unit-density denominator retained {metrics['decoupled']['baseline_density_primary_retained']} of 23 primary countries.",
        f"- Countries not retained under the 2024-density denominator: {', '.join(metrics['decoupled']['baseline_density_not_retained'])}.",
        f"- A stricter Q80 growth plus Q20 density definition retained {metrics['decoupled']['strict_primary_retained']} of 23 primary countries.",
        "",
        "## Missing-Resource Bounds",
        "",
        f"- Global country-level proportion: {global_row['Lower-bound proportion, %']:.1f}% to {global_row['Upper-bound proportion, %']:.1f}%.",
        f"- AFRO country-level proportion: {afro['Lower-bound proportion, %']:.1f}% to {afro['Upper-bound proportion, %']:.1f}%.",
        "- Lower bounds treat all DIRAC-missing countries as screen-negative; upper bounds treat high-growth DIRAC-missing countries as meeting the resource criterion.",
        "",
        "## Construct Diagnostics",
        "",
        frame_to_markdown(diagnostics),
        "",
        "No hypothesis tests were used. These analyses describe the complete available country records and assess classification stability.",
    ]
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    (REPORT_DIR / "report_08_radiology_v6_additional_analyses.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_DATASET.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(CURRENT_DATASET)
    df["country_display"] = df["country_iso3"].map(DISPLAY_COUNTRIES).fillna(df["gco_country"])
    numeric_columns = [
        "mv_therapy_units",
        "dirac_last_update_year",
        "selected_site_cases_2024",
        "selected_site_cases_2040",
        "selected_site_cases_2050",
        "selected_site_relative_case_growth_2050",
        "mv_units_per_1000_selected_site_cases_2050",
        "all_cancer_relative_case_growth_2050",
        "population_2024",
        "population_2050",
    ]
    for column in numeric_columns:
        df[column] = pd.to_numeric(df[column], errors="coerce")

    df, version_metrics = add_version_comparison(df)
    df, decoupled_metrics, threshold_grid = add_decoupled_and_threshold_analyses(df)
    df = add_sensitivity_support(df)
    bounds = build_missingness_bounds(df)
    diagnostics = build_construct_diagnostics(df)
    version_table = build_version_table(df)
    country_stability = build_country_stability_table(df)

    df.to_csv(OUT_DATASET, index=False, encoding="utf-8-sig")
    version_table.to_csv(TABLE_DIR / "table_19_gco_2022_2024_burden_version_reclassification.csv", index=False, encoding="utf-8-sig")
    threshold_grid.to_csv(TABLE_DIR / "table_20_decoupled_density_threshold_stability.csv", index=False, encoding="utf-8-sig")
    bounds.to_csv(TABLE_DIR / "table_21_missing_resource_identification_bounds.csv", index=False, encoding="utf-8-sig")
    diagnostics.to_csv(TABLE_DIR / "table_22_construct_diagnostics.csv", index=False, encoding="utf-8-sig")
    country_stability.to_csv(TABLE_DIR / "table_23_radiology_v6_country_stability.csv", index=False, encoding="utf-8-sig")
    (TABLE_DIR / "table_23_radiology_v6_country_stability.md").write_text(
        frame_to_markdown(country_stability, float_digits=0), encoding="utf-8"
    )

    metrics = {"version": version_metrics, "decoupled": decoupled_metrics}
    write_report(metrics, diagnostics, bounds)
    LOG_PATH.write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    print(json.dumps(metrics, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
