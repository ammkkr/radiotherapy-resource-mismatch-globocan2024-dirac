"""Build reproducibility analyses for the associated Research Letter.

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
GCO_2024_PROJECTIONS = PROJECT_ROOT / "data" / "source" / "gco_cancer_tomorrow_2024_predictions_long.csv"
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


def build_gco_selected_site_burden(path: Path, version: str) -> pd.DataFrame:
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
    gco = pd.read_csv(path, usecols=usecols)
    gco["cancer_code"] = pd.to_numeric(gco["cancer_code"], errors="coerce")
    gco["year"] = pd.to_numeric(gco["year"], errors="coerce")
    gco["predicted_count"] = pd.to_numeric(gco["predicted_count"], errors="coerce")
    gco["pop"] = pd.to_numeric(gco["pop"], errors="coerce")

    selected = gco[
        gco["country_iso3"].notna()
        & gco["sex_label"].eq("both")
        & gco["measure"].eq("incidence")
        & gco["cancer_code"].isin(SELECTED_SITE_CODES)
        & gco["year"].isin([2022, 2024, 2025, 2040, 2050])
    ].copy()
    burden = (
        selected.groupby(["country_iso3", "year"], as_index=False)["predicted_count"]
        .sum()
        .pivot(index="country_iso3", columns="year", values="predicted_count")
        .reset_index()
    )
    burden = burden.rename(
        columns={year: f"gco{version}_selected_cases_{int(year)}" for year in burden.columns if year != "country_iso3"}
    )
    baseline_year = 2022 if version == "2022" else 2024
    burden[f"gco{version}_selected_growth_{baseline_year}_2050"] = (
        burden[f"gco{version}_selected_cases_2050"] - burden[f"gco{version}_selected_cases_{baseline_year}"]
    ) / burden[f"gco{version}_selected_cases_{baseline_year}"]
    burden[f"gco{version}_selected_growth_2025_2050"] = (
        burden[f"gco{version}_selected_cases_2050"] - burden[f"gco{version}_selected_cases_2025"]
    ) / burden[f"gco{version}_selected_cases_2025"]
    return burden


def classification_metrics(old_set: set[str], new_set: set[str], common_count: int) -> dict[str, Any]:
    union = old_set | new_set
    return {
        "common_matched_countries": common_count,
        "gco2022_screen_positive_count": len(old_set),
        "gco2024_screen_positive_count_common": len(new_set),
        "retained_count": len(old_set & new_set),
        "new_in_2024": sorted(new_set - old_set),
        "no_longer_in_2024": sorted(old_set - new_set),
        "jaccard": len(old_set & new_set) / len(union) if union else np.nan,
    }


def add_version_comparison(df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    burden_2022 = build_gco_selected_site_burden(GCO_2022_PROJECTIONS, "2022")
    burden_2024 = build_gco_selected_site_burden(GCO_2024_PROJECTIONS, "2024")
    out = df.merge(burden_2022, on="country_iso3", how="left", validate="one_to_one")
    out = out.merge(burden_2024, on="country_iso3", how="left", validate="one_to_one")

    # Legacy operational comparison: each release uses its native baseline year.
    matched_common = (
        flag(out["dirac_matched"])
        & out["mv_therapy_units"].notna()
        & out["gco2022_selected_growth_2022_2050"].notna()
    )
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
    legacy_metrics = classification_metrics(old_set, current_common_set, int(matched_common.sum()))
    legacy_metrics.update(
        {
            "comparison_baseline": "native release years (2022 and 2024)",
            "gco2022_growth_q75": growth_q75,
            "gco2022_density_q25": density_q25,
        }
    )

    # Reviewer-facing comparison: both burden versions use the same 2025 baseline.
    common_year_matched = (
        flag(out["dirac_matched"])
        & out["mv_therapy_units"].notna()
        & out["gco2022_selected_growth_2025_2050"].notna()
        & out["gco2024_selected_growth_2025_2050"].notna()
    )
    out["gco2024_mv_per_1000_selected_cases_2050"] = np.where(
        common_year_matched,
        out["mv_therapy_units"] * 1000 / out["gco2024_selected_cases_2050"],
        np.nan,
    )
    old_growth_q75 = quantile(out.loc[common_year_matched, "gco2022_selected_growth_2025_2050"], 0.75)
    old_density_q25 = quantile(out.loc[common_year_matched, "gco2022_mv_per_1000_selected_cases_2050"], 0.25)
    new_growth_q75 = quantile(out.loc[common_year_matched, "gco2024_selected_growth_2025_2050"], 0.75)
    new_density_q25 = quantile(out.loc[common_year_matched, "gco2024_mv_per_1000_selected_cases_2050"], 0.25)
    out["gco2022_common_year_screen_positive"] = (
        common_year_matched
        & out["gco2022_selected_growth_2025_2050"].ge(old_growth_q75)
        & out["gco2022_mv_per_1000_selected_cases_2050"].le(old_density_q25)
    )
    out["gco2024_common_year_screen_positive"] = (
        common_year_matched
        & out["gco2024_selected_growth_2025_2050"].ge(new_growth_q75)
        & out["gco2024_mv_per_1000_selected_cases_2050"].le(new_density_q25)
    )
    out["gco2022_common_year_growth_percentile"] = np.nan
    out["gco2024_common_year_growth_percentile"] = np.nan
    out.loc[common_year_matched, "gco2022_common_year_growth_percentile"] = (
        out.loc[common_year_matched, "gco2022_selected_growth_2025_2050"].rank(method="average", pct=True) * 100
    )
    out.loc[common_year_matched, "gco2024_common_year_growth_percentile"] = (
        out.loc[common_year_matched, "gco2024_selected_growth_2025_2050"].rank(method="average", pct=True) * 100
    )
    out["common_year_burden_version_transition"] = "Not screen-positive in either version"
    out.loc[
        common_year_matched & out["gco2022_common_year_screen_positive"] & out["gco2024_common_year_screen_positive"],
        "common_year_burden_version_transition",
    ] = "Retained"
    out.loc[
        common_year_matched & ~out["gco2022_common_year_screen_positive"] & out["gco2024_common_year_screen_positive"],
        "common_year_burden_version_transition",
    ] = "New in GCO version 2024"
    out.loc[
        common_year_matched & out["gco2022_common_year_screen_positive"] & ~out["gco2024_common_year_screen_positive"],
        "common_year_burden_version_transition",
    ] = "No longer screen-positive"
    out.loc[~common_year_matched, "common_year_burden_version_transition"] = "Not comparable"

    old_common_year_set = set(out.loc[out["gco2022_common_year_screen_positive"], "country_iso3"])
    new_common_year_set = set(out.loc[out["gco2024_common_year_screen_positive"], "country_iso3"])
    metrics = classification_metrics(old_common_year_set, new_common_year_set, int(common_year_matched.sum()))
    metrics.update(
        {
            "comparison_baseline": "2025 for both burden versions",
            "gco2022_growth_q75": old_growth_q75,
            "gco2022_density_q25": old_density_q25,
            "gco2024_growth_q75": new_growth_q75,
            "gco2024_density_q25": new_density_q25,
            "legacy_native_baseline_comparison": legacy_metrics,
        }
    )
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


def add_resource_recency_analysis(
    df: pd.DataFrame, cutoff_year: int = 2023
) -> tuple[pd.DataFrame, dict[str, Any], pd.DataFrame, pd.DataFrame]:
    out = df.copy()
    matched = flag(out["dirac_matched"]) & out["mv_therapy_units"].notna()
    record_year = pd.to_numeric(out["dirac_last_update_year"], errors="coerce")
    recent = matched & record_year.ge(cutoff_year)
    growth = pd.to_numeric(out["selected_site_relative_case_growth_2050"], errors="coerce")
    density = pd.to_numeric(out["mv_units_per_1000_selected_site_cases_2050"], errors="coerce")

    recent_growth_q75 = quantile(growth[recent], 0.75)
    recent_density_q25 = quantile(density[recent], 0.25)
    out["recent_dirac_record"] = recent
    out["recent_record_screen_positive"] = (
        recent & growth.ge(recent_growth_q75) & density.le(recent_density_q25)
    )

    primary_set = set(out.loc[out["current_screen_positive"], "country_iso3"])
    recent_set = set(out.loc[out["recent_record_screen_positive"], "country_iso3"])
    eligible_primary_set = set(out.loc[out["current_screen_positive"] & recent, "country_iso3"])
    old_primary_set = primary_set - eligible_primary_set
    retained_eligible_set = eligible_primary_set & recent_set

    metrics = {
        "cutoff_year": cutoff_year,
        "primary_matched_countries": int(matched.sum()),
        "recent_record_countries": int(recent.sum()),
        "older_record_countries": int((matched & ~recent).sum()),
        "recent_growth_q75": recent_growth_q75,
        "recent_density_q25": recent_density_q25,
        "recent_screen_positive_count": len(recent_set),
        "primary_recent_record_count": len(eligible_primary_set),
        "primary_older_record_count": len(old_primary_set),
        "primary_older_record_iso3": sorted(old_primary_set),
        "eligible_primary_retained_count": len(retained_eligible_set),
        "eligible_primary_not_retained_iso3": sorted(eligible_primary_set - recent_set),
        "new_recent_screen_positive_iso3": sorted(recent_set - primary_set),
        "recent_screen_positive_iso3": sorted(recent_set),
    }

    summary = pd.DataFrame(
        [
            {
                "Analysis": "Primary, all latest-reported DIRAC years",
                "Countries analyzed, No.": int(matched.sum()),
                "Growth threshold, %": 100 * quantile(growth[matched], 0.75),
                "Density threshold": quantile(density[matched], 0.25),
                "Countries meeting both thresholds, No.": len(primary_set),
                "Primary recent-record countries retained, No.": len(eligible_primary_set),
                "ISO3 meeting both thresholds": ";".join(sorted(primary_set)),
            },
            {
                "Analysis": f"Restricted to DIRAC record year {cutoff_year} or later",
                "Countries analyzed, No.": int(recent.sum()),
                "Growth threshold, %": 100 * recent_growth_q75,
                "Density threshold": recent_density_q25,
                "Countries meeting both thresholds, No.": len(recent_set),
                "Primary recent-record countries retained, No.": len(retained_eligible_set),
                "ISO3 meeting both thresholds": ";".join(sorted(recent_set)),
            },
        ]
    )

    detail = out.loc[matched, [
        "country_iso3",
        "country_display",
        "who_region",
        "income_label",
        "dirac_last_update_year",
        "selected_site_relative_case_growth_2050",
        "mv_units_per_1000_selected_site_cases_2050",
        "current_screen_positive",
        "recent_dirac_record",
        "recent_record_screen_positive",
    ]].copy()
    detail["recency_disposition"] = "Not classified"
    detail.loc[
        detail["current_screen_positive"] & ~detail["recent_dirac_record"], "recency_disposition"
    ] = f"Primary country; record before {cutoff_year}"
    detail.loc[
        detail["current_screen_positive"] & detail["recent_dirac_record"] & detail["recent_record_screen_positive"],
        "recency_disposition",
    ] = "Primary country retained"
    detail.loc[
        detail["current_screen_positive"] & detail["recent_dirac_record"] & ~detail["recent_record_screen_positive"],
        "recency_disposition",
    ] = "Primary country not retained"
    detail.loc[
        ~detail["current_screen_positive"] & detail["recent_record_screen_positive"], "recency_disposition"
    ] = "New in recent-record analysis"
    return out, metrics, summary, detail


def summarize_all_cancer_sensitivity(df: pd.DataFrame) -> dict[str, Any]:
    robustness = pd.read_csv(EXISTING_ROBUSTNESS)
    row = robustness.loc[robustness["Analysis"].eq("All-cancer burden")].iloc[0]
    all_cancer_set = iso_set(row["Flagged ISO3"])
    primary_set = set(df.loc[df["current_screen_positive"], "country_iso3"])
    return {
        "screen_positive_count": len(all_cancer_set),
        "primary_retained_count": len(primary_set & all_cancer_set),
        "primary_not_retained_iso3": sorted(primary_set - all_cancer_set),
        "additional_iso3": sorted(all_cancer_set - primary_set),
        "growth_threshold_percent": float(row["Growth threshold, %"]),
        "density_threshold": float(row["Resource threshold"]),
    }


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


def build_common_year_version_table(df: pd.DataFrame) -> pd.DataFrame:
    common = df[df["common_year_burden_version_transition"].ne("Not comparable")].copy()
    priority = {
        "New in GCO version 2024": 0,
        "No longer screen-positive": 1,
        "Retained": 2,
        "Not screen-positive in either version": 3,
    }
    common["transition_order"] = common["common_year_burden_version_transition"].map(priority)
    columns = [
        "country_iso3",
        "country_display",
        "who_region",
        "income_label",
        "gco2022_selected_cases_2025",
        "gco2022_selected_cases_2050",
        "gco2022_selected_growth_2025_2050",
        "gco2022_common_year_growth_percentile",
        "gco2022_mv_per_1000_selected_cases_2050",
        "gco2024_selected_cases_2025",
        "gco2024_selected_cases_2050",
        "gco2024_selected_growth_2025_2050",
        "gco2024_common_year_growth_percentile",
        "gco2024_mv_per_1000_selected_cases_2050",
        "gco2022_common_year_screen_positive",
        "gco2024_common_year_screen_positive",
        "common_year_burden_version_transition",
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
    version_labels = {
        "Retained": "Retained",
        "New in GCO version 2024": "New in version 2024",
        "No longer screen-positive": "Not in version-2024 common-year screen",
        "Not screen-positive in either version": "Not in common-year screen",
    }
    primary["Common-2025 burden-version status"] = primary[
        "common_year_burden_version_transition"
    ].map(version_labels).fillna("Not comparable")
    primary["Recent-record sensitivity"] = np.select(
        [
            ~primary["recent_dirac_record"],
            primary["recent_record_screen_positive"],
        ],
        [
            "Record before 2023",
            "Retained",
        ],
        default="Not retained",
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
            "Common-2025 burden-version status",
            "Recent-record sensitivity",
            "Sensitivity support, No./7",
        ]
    ].rename(columns={"country_display": "Country"})


def build_jama_health_forum_stability_table(
    df: pd.DataFrame, metrics: dict[str, Any]
) -> pd.DataFrame:
    version = metrics["version"]
    decoupled = metrics["decoupled"]
    recency = metrics["recency"]
    all_cancer = metrics["all_cancer"]
    primary_count = int(df["current_screen_positive"].sum())
    strict_losses = ", ".join(decoupled["strict_not_retained"])
    old_records = ", ".join(recency["primary_older_record_iso3"])
    return pd.DataFrame(
        [
            {
                "Analysis": "Primary selected-site screen",
                "Countries analyzed, No.": 150,
                "Growth / density threshold": f"119.2% / {decoupled['density_2050_q25']:.3f}",
                "Countries meeting both, No.": primary_count,
                "Stability result": "Reference",
                "Differences": "None",
            },
            {
                "Analysis": "Burden versions, common 2025 baseline",
                "Countries analyzed, No.": version["common_matched_countries"],
                "Growth / density threshold": (
                    f"v2022 {100 * version['gco2022_growth_q75']:.1f}% / {version['gco2022_density_q25']:.3f}; "
                    f"v2024 {100 * version['gco2024_growth_q75']:.1f}% / {version['gco2024_density_q25']:.3f}"
                ),
                "Countries meeting both, No.": (
                    f"{version['gco2022_screen_positive_count']} / {version['gco2024_screen_positive_count_common']}"
                ),
                "Stability result": (
                    f"{version['retained_count']}/{version['gco2022_screen_positive_count']} "
                    "version-2022 priorities retained"
                ),
                "Differences": (
                    f"New: {', '.join(version['new_in_2024'])}; "
                    f"no longer: {', '.join(version['no_longer_in_2024'])}"
                ),
            },
            {
                "Analysis": "All-cancer burden",
                "Countries analyzed, No.": 150,
                "Growth / density threshold": (
                    f"{all_cancer['growth_threshold_percent']:.1f}% / {all_cancer['density_threshold']:.3f}"
                ),
                "Countries meeting both, No.": all_cancer["screen_positive_count"],
                "Stability result": f"{all_cancer['primary_retained_count']}/23 primary retained",
                "Differences": f"Added: {', '.join(all_cancer['additional_iso3'])}",
            },
            {
                "Analysis": "2024-case density denominator",
                "Countries analyzed, No.": 150,
                "Growth / density threshold": (
                    f"{100 * decoupled['growth_q75']:.1f}% / {decoupled['density_2024_q25']:.3f}"
                ),
                "Countries meeting both, No.": decoupled["baseline_density_screen_positive_count"],
                "Stability result": f"{decoupled['baseline_density_primary_retained']}/23 primary retained",
                "Differences": f"Not retained: {', '.join(decoupled['baseline_density_not_retained'])}",
            },
            {
                "Analysis": "Q80-growth / Q20-density screen",
                "Countries analyzed, No.": 150,
                "Growth / density threshold": (
                    f"{100 * decoupled['growth_q80']:.1f}% / {decoupled['density_q20']:.3f}"
                ),
                "Countries meeting both, No.": decoupled["strict_screen_positive_count"],
                "Stability result": f"{decoupled['strict_primary_retained']}/23 primary retained",
                "Differences": f"Not retained: {strict_losses}",
            },
            {
                "Analysis": "DIRAC record year 2023 or later",
                "Countries analyzed, No.": recency["recent_record_countries"],
                "Growth / density threshold": (
                    f"{100 * recency['recent_growth_q75']:.1f}% / {recency['recent_density_q25']:.3f}"
                ),
                "Countries meeting both, No.": recency["recent_screen_positive_count"],
                "Stability result": (
                    f"{recency['eligible_primary_retained_count']}/{recency['primary_recent_record_count']} "
                    "eligible primary retained"
                ),
                "Differences": (
                    f"Older primary records: {old_records}; added: "
                    f"{', '.join(recency['new_recent_screen_positive_iso3'])}"
                ),
            },
        ]
    )


def write_report(metrics: dict[str, Any], diagnostics: pd.DataFrame, bounds: pd.DataFrame) -> None:
    global_row = bounds[(bounds["Stratum"] == "Global")].iloc[0]
    afro = bounds[(bounds["Stratum"] == "WHO region") & (bounds["Category"] == "AFRO")].iloc[0]
    version = metrics["version"]
    legacy = version["legacy_native_baseline_comparison"]
    recency = metrics["recency"]
    all_cancer = metrics["all_cancer"]
    old_primary = ", ".join(recency["primary_older_record_iso3"]) or "none"
    recency_losses = ", ".join(recency["eligible_primary_not_retained_iso3"]) or "none"
    recency_additions = ", ".join(recency["new_recent_screen_positive_iso3"]) or "none"
    lines = [
        "# Radiology v6 Additional Analyses",
        "",
        "## Burden-Version Update With a Common Baseline",
        "",
        f"- Both burden versions used 2025 as the baseline and 2050 as the horizon across {version['common_matched_countries']} common matched countries.",
        f"- GCO version 2022 identified {version['gco2022_screen_positive_count']} countries and version 2024 identified {version['gco2024_screen_positive_count_common']}.",
        f"- Retained across versions: {version['retained_count']}.",
        f"- New in version 2024: {', '.join(version['new_in_2024']) or 'none'}.",
        f"- No longer screen-positive in version 2024: {', '.join(version['no_longer_in_2024']) or 'none'}.",
        "- The same latest-reported DIRAC resource snapshot was fixed in both versions, so the common-year comparison isolates burden-data revision from baseline-year and resource changes.",
        f"- The legacy native-baseline comparison retained {legacy['retained_count']} countries and changed {len(legacy['new_in_2024']) + len(legacy['no_longer_in_2024'])} classifications.",
        "",
        "## Denominator and Threshold Sensitivity",
        "",
        f"- Using 2024 rather than 2050 selected-site cases in the MV-unit-density denominator retained {metrics['decoupled']['baseline_density_primary_retained']} of 23 primary countries.",
        f"- Countries not retained under the 2024-density denominator: {', '.join(metrics['decoupled']['baseline_density_not_retained'])}.",
        f"- A stricter Q80 growth plus Q20 density definition retained {metrics['decoupled']['strict_primary_retained']} of 23 primary countries.",
        f"- The all-cancer burden definition retained {all_cancer['primary_retained_count']} of 23 primary countries and identified {all_cancer['screen_positive_count']} countries overall.",
        "",
        "## DIRAC Record-Year Sensitivity",
        "",
        f"- Restricting to DIRAC records from {recency['cutoff_year']} or later retained {recency['recent_record_countries']} of {recency['primary_matched_countries']} matched countries.",
        f"- {recency['primary_older_record_count']} primary countries had older resource records: {old_primary}.",
        f"- Among {recency['primary_recent_record_count']} primary countries eligible for the restricted analysis, {recency['eligible_primary_retained_count']} were retained; not retained: {recency_losses}.",
        f"- Newly identified in the recent-record analysis: {recency_additions}.",
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
    (REPORT_DIR / "additional_analyses.md").write_text("\n".join(lines), encoding="utf-8")


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
    df, recency_metrics, recency_summary, recency_detail = add_resource_recency_analysis(df)
    all_cancer_metrics = summarize_all_cancer_sensitivity(df)
    metrics = {
        "version": version_metrics,
        "decoupled": decoupled_metrics,
        "recency": recency_metrics,
        "all_cancer": all_cancer_metrics,
    }
    bounds = build_missingness_bounds(df)
    diagnostics = build_construct_diagnostics(df)
    version_table = build_version_table(df)
    common_year_version_table = build_common_year_version_table(df)
    country_stability = build_country_stability_table(df)
    jama_health_forum_stability = build_jama_health_forum_stability_table(df, metrics)

    df.to_csv(OUT_DATASET, index=False, encoding="utf-8-sig")
    version_table.to_csv(TABLE_DIR / "table_19_burden_version_reclassification.csv", index=False, encoding="utf-8-sig")
    threshold_grid.to_csv(TABLE_DIR / "table_20_denominator_threshold_stability.csv", index=False, encoding="utf-8-sig")
    bounds.to_csv(TABLE_DIR / "table_21_missing_resource_identification_bounds.csv", index=False, encoding="utf-8-sig")
    diagnostics.to_csv(TABLE_DIR / "table_22_construct_diagnostics.csv", index=False, encoding="utf-8-sig")
    country_stability.to_csv(TABLE_DIR / "table_23_country_stability.csv", index=False, encoding="utf-8-sig")
    common_year_version_table.to_csv(
        TABLE_DIR / "table_24_common_2025_baseline_burden_version_reclassification.csv",
        index=False,
        encoding="utf-8-sig",
    )
    recency_summary.to_csv(
        TABLE_DIR / "table_25_dirac_record_recency_sensitivity.csv", index=False, encoding="utf-8-sig"
    )
    recency_detail.to_csv(
        TABLE_DIR / "table_26_dirac_record_recency_country_detail.csv", index=False, encoding="utf-8-sig"
    )
    country_stability.to_csv(
        TABLE_DIR / "table_27_country_stability_updated.csv", index=False, encoding="utf-8-sig"
    )
    jama_health_forum_stability.to_csv(
        TABLE_DIR / "table_28_stability_summary.csv", index=False, encoding="utf-8-sig"
    )
    (TABLE_DIR / "table_28_stability_summary.md").write_text(
        frame_to_markdown(jama_health_forum_stability), encoding="utf-8"
    )
    (TABLE_DIR / "table_23_radiology_v6_country_stability.md").write_text(
        frame_to_markdown(country_stability, float_digits=0), encoding="utf-8"
    )

    write_report(metrics, diagnostics, bounds)
    LOG_PATH.write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    print(json.dumps(metrics, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
