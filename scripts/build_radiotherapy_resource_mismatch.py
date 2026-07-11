"""Build GCO Cancer Tomorrow / DIRAC radiotherapy resource-density dataset."""

from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_INTERIM = PROJECT_ROOT / "data" / "interim"
DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"
TABLE_DIR = PROJECT_ROOT / "results" / "tables"
REPORT_DIR = PROJECT_ROOT / "results" / "reports"
PROTOCOL_DIR = PROJECT_ROOT / "docs" / "protocol"
LOG_DIR = PROJECT_ROOT / "logs"

GCO_PREDICTIONS = DATA_INTERIM / "gco_cancer_tomorrow_2024_predictions_long.csv"
DIRAC_RESOURCES = DATA_INTERIM / "dirac_country_resources.csv"

RT_SITE_RULES = [
    {"cancer_code": "1", "label": "Lip, oral cavity", "selected_site_relevance": "high", "selected_site_weight": "1", "included_selected_site_set": "1", "rationale": "Head and neck cancer site with common radiotherapy use."},
    {"cancer_code": "3", "label": "Oropharynx", "selected_site_relevance": "high", "selected_site_weight": "1", "included_selected_site_set": "1", "rationale": "Head and neck cancer site with common radiotherapy use."},
    {"cancer_code": "4", "label": "Nasopharynx", "selected_site_relevance": "high", "selected_site_weight": "1", "included_selected_site_set": "1", "rationale": "Head and neck cancer site with common radiotherapy use."},
    {"cancer_code": "5", "label": "Hypopharynx", "selected_site_relevance": "high", "selected_site_weight": "1", "included_selected_site_set": "1", "rationale": "Head and neck cancer site with common radiotherapy use."},
    {"cancer_code": "6", "label": "Oesophagus", "selected_site_relevance": "moderate_high", "selected_site_weight": "1", "included_selected_site_set": "1", "rationale": "Radiotherapy or chemoradiotherapy is commonly part of treatment in selected settings."},
    {"cancer_code": "9", "label": "Rectum", "selected_site_relevance": "moderate_high", "selected_site_weight": "1", "included_selected_site_set": "1", "rationale": "Rectal cancer is separated from grouped colorectum where available."},
    {"cancer_code": "14", "label": "Larynx", "selected_site_relevance": "high", "selected_site_weight": "1", "included_selected_site_set": "1", "rationale": "Head and neck cancer site with common radiotherapy use."},
    {"cancer_code": "15", "label": "Trachea, bronchus and lung", "selected_site_relevance": "moderate", "selected_site_weight": "1", "included_selected_site_set": "1", "rationale": "Radiotherapy is used across curative and palliative lung cancer pathways."},
    {"cancer_code": "20", "label": "Breast", "selected_site_relevance": "high", "selected_site_weight": "1", "included_selected_site_set": "1", "rationale": "Radiotherapy is a core component after breast-conserving surgery and in selected postmastectomy settings."},
    {"cancer_code": "23", "label": "Cervix uteri", "selected_site_relevance": "high", "selected_site_weight": "1", "included_selected_site_set": "1", "rationale": "Radiotherapy and brachytherapy are central to treatment of locally advanced disease."},
    {"cancer_code": "24", "label": "Corpus uteri", "selected_site_relevance": "moderate", "selected_site_weight": "1", "included_selected_site_set": "1", "rationale": "Radiotherapy is used in selected adjuvant and advanced settings."},
    {"cancer_code": "27", "label": "Prostate", "selected_site_relevance": "moderate_high", "selected_site_weight": "1", "included_selected_site_set": "1", "rationale": "External-beam radiotherapy and brachytherapy are common prostate cancer treatments."},
    {"cancer_code": "30", "label": "Bladder", "selected_site_relevance": "moderate", "selected_site_weight": "1", "included_selected_site_set": "1", "rationale": "Radiotherapy is used in bladder-preserving and palliative settings."},
    {"cancer_code": "31", "label": "Brain, central nervous system", "selected_site_relevance": "high", "selected_site_weight": "1", "included_selected_site_set": "1", "rationale": "Radiotherapy is central for many primary CNS tumour pathways."},
]

ALL_CANCER_EXCL_NMSC_CODE = "40"
ALL_CANCER_CODE = "39"
YEARS = ["2024", "2030", "2040", "2050"]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, records: list[dict[str, Any]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = sorted({key for row in records for key in row})
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)


def to_float(value: Any) -> float:
    if value is None or value == "":
        return float("nan")
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def safe_div(numer: float, denom: float) -> float:
    if denom == 0 or math.isnan(denom):
        return float("nan")
    return numer / denom


def quantile(values: list[float], q: float) -> float:
    clean = sorted(v for v in values if not math.isnan(v))
    if not clean:
        return float("nan")
    pos = (len(clean) - 1) * q
    lower = math.floor(pos)
    upper = math.ceil(pos)
    if lower == upper:
        return clean[int(pos)]
    return clean[lower] + (clean[upper] - clean[lower]) * (pos - lower)


def percentile_ranks(values_by_key: dict[str, float], reverse: bool = False) -> dict[str, float]:
    clean = [(key, value) for key, value in values_by_key.items() if not math.isnan(value)]
    if not clean:
        return {}
    clean.sort(key=lambda item: item[1])
    n = len(clean)
    if n == 1:
        return {clean[0][0]: 1.0}
    ranks: dict[str, float] = {}
    for idx, (key, _) in enumerate(clean):
        base = idx / (n - 1)
        ranks[key] = 1 - base if reverse else base
    return ranks


def write_site_rules() -> None:
    path = PROTOCOL_DIR / "selected_cancer_site_rules.csv"
    write_csv(
        path,
        RT_SITE_RULES,
        ["cancer_code", "label", "selected_site_relevance", "selected_site_weight", "included_selected_site_set", "rationale"],
    )


def build_burden_tables(predictions: list[dict[str, str]]) -> tuple[dict[str, dict[str, Any]], dict[str, dict[str, Any]]]:
    rt_codes = {row["cancer_code"] for row in RT_SITE_RULES if row["included_selected_site_set"] == "1"}
    country_meta: dict[str, dict[str, Any]] = {}
    rt_incidence: dict[tuple[str, str], float] = defaultdict(float)
    rt_mortality: dict[tuple[str, str], float] = defaultdict(float)
    all_incidence: dict[tuple[str, str], float] = defaultdict(float)
    all_mortality: dict[tuple[str, str], float] = defaultdict(float)
    population: dict[tuple[str, str], float] = {}

    for row in predictions:
        iso3 = row.get("country_iso3", "").strip().upper()
        if not iso3:
            continue
        if row.get("sex_label") != "both":
            continue
        year = row.get("year")
        if year not in YEARS:
            continue

        country_meta.setdefault(
            iso3,
            {
                "country_iso3": iso3,
                "gco_country": row.get("population_label", ""),
                "country_code": row.get("country_code", ""),
                "who_region": row.get("who_region", ""),
                "hdi_label": row.get("hdi_label", ""),
                "income_label": row.get("income_label", ""),
            },
        )
        cancer_code = row.get("cancer_code", "")
        measure = row.get("measure", "")
        count = to_float(row.get("predicted_count"))
        pop = to_float(row.get("pop"))
        if not math.isnan(pop):
            population[(iso3, year)] = pop

        if cancer_code in rt_codes:
            if measure == "incidence":
                rt_incidence[(iso3, year)] += count
            elif measure == "mortality":
                rt_mortality[(iso3, year)] += count
        if cancer_code == ALL_CANCER_EXCL_NMSC_CODE:
            if measure == "incidence":
                all_incidence[(iso3, year)] += count
            elif measure == "mortality":
                all_mortality[(iso3, year)] += count
        elif cancer_code == ALL_CANCER_CODE:
            # Stored for validation/fallback only when excl-NMSC code is not present.
            key = (iso3, year)
            if measure == "incidence" and all_incidence.get(key, 0) == 0:
                all_incidence[key] += count
            elif measure == "mortality" and all_mortality.get(key, 0) == 0:
                all_mortality[key] += count

    burden: dict[str, dict[str, Any]] = {}
    for iso3, meta in country_meta.items():
        rec = dict(meta)
        for year in YEARS:
            rec[f"rt_relevant_cases_{year}"] = rt_incidence.get((iso3, year), float("nan"))
            rec[f"rt_relevant_deaths_{year}"] = rt_mortality.get((iso3, year), float("nan"))
            rec[f"all_cancer_cases_{year}"] = all_incidence.get((iso3, year), float("nan"))
            rec[f"all_cancer_deaths_{year}"] = all_mortality.get((iso3, year), float("nan"))
            rec[f"population_{year}"] = population.get((iso3, year), float("nan"))

        rec["selected_site_case_increase_2050"] = rec["selected_site_cases_2050"] - rec["selected_site_cases_2024"]
        rec["selected_site_relative_case_growth_2050"] = safe_div(
            rec["selected_site_case_increase_2050"], rec["selected_site_cases_2024"]
        )
        rec["selected_site_death_increase_2050"] = rec["selected_site_deaths_2050"] - rec["selected_site_deaths_2024"]
        rec["selected_site_relative_death_growth_2050"] = safe_div(
            rec["selected_site_death_increase_2050"], rec["selected_site_deaths_2024"]
        )
        rec["all_cancer_absolute_case_increase_2050"] = rec["all_cancer_cases_2050"] - rec["all_cancer_cases_2024"]
        rec["all_cancer_relative_case_growth_2050"] = safe_div(
            rec["all_cancer_absolute_case_increase_2050"], rec["all_cancer_cases_2024"]
        )
        burden[iso3] = rec

    site_summary = build_site_summary(predictions, rt_codes)
    return burden, site_summary


def build_site_summary(predictions: list[dict[str, str]], rt_codes: set[str]) -> dict[str, dict[str, Any]]:
    summary: dict[str, dict[str, Any]] = {}
    for row in predictions:
        if row.get("country_code") != "900":
            continue
        if row.get("sex_label") != "both" or row.get("measure") != "incidence":
            continue
        if row.get("year") not in ["2024", "2050"]:
            continue
        cancer_code = row.get("cancer_code", "")
        if cancer_code not in rt_codes:
            continue
        rec = summary.setdefault(
            cancer_code,
            {
                "cancer_code": cancer_code,
                "cancer_label": row.get("cancer_label", ""),
                "included_included_selected_site_set": "1",
            },
        )
        rec[f"world_cases_{row.get('year')}"] = to_float(row.get("predicted_count"))
    for rec in summary.values():
        rec["world_absolute_case_increase_2050"] = rec.get("world_cases_2050", float("nan")) - rec.get(
            "world_cases_2024", float("nan")
        )
        rec["world_relative_case_growth_2050"] = safe_div(
            rec["world_absolute_case_increase_2050"], rec.get("world_cases_2024", float("nan"))
        )
    return summary


def merge_with_dirac(burden: dict[str, dict[str, Any]], dirac: list[dict[str, str]]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    dirac_by_iso = {row["dirac_iso3"].strip().upper(): row for row in dirac if row.get("dirac_iso3")}
    records: list[dict[str, Any]] = []

    for iso3, burden_row in burden.items():
        rec = dict(burden_row)
        drow = dirac_by_iso.get(iso3)
        rec["dirac_matched"] = "1" if drow else "0"
        if drow:
            rec.update(drow)
        else:
            for field in [
                "dirac_iso3",
                "dirac_country",
                "dirac_region",
                "dirac_last_update_year",
                "rt_centres_total",
                "rt_centres_with_rt",
                "mv_therapy_units",
                "proton_ion_units",
                "kv_therapy_units",
                "brachytherapy_units",
                "ct_units",
                "simulators",
                "treatment_planning_systems",
            ]:
                rec[field] = ""

        mv_units = to_float(rec.get("mv_therapy_units"))
        centres = to_float(rec.get("rt_centres_with_rt"))
        brachy = to_float(rec.get("brachytherapy_units"))
        rt_cases_2050 = to_float(rec.get("selected_site_cases_2050"))
        all_cases_2050 = to_float(rec.get("all_cancer_cases_2050"))
        pop_2024 = to_float(rec.get("population_2024"))

        rec["mv_units_per_1000_selected_site_cases_2050"] = safe_div(mv_units * 1000, rt_cases_2050)
        rec["selected_site_cases_per_mv_unit_2050"] = safe_div(rt_cases_2050, mv_units)
        rec["rt_centres_per_10000_selected_site_cases_2050"] = safe_div(centres * 10000, rt_cases_2050)
        rec["mv_units_per_1000_all_cancer_cases_2050"] = safe_div(mv_units * 1000, all_cases_2050)
        rec["mv_units_per_million_population_2024"] = safe_div(mv_units * 1_000_000, pop_2024)
        rec["brachy_units_per_1000_selected_site_cases_2050"] = safe_div(brachy * 1000, rt_cases_2050)
        records.append(rec)

    dirac_unmatched = sorted(set(dirac_by_iso) - set(burden))
    gco_unmatched = sorted(set(burden) - set(dirac_by_iso))
    merge_summary = {
        "gco_country_records": len(burden),
        "dirac_country_records": len(dirac_by_iso),
        "merged_records": len(records),
        "matched_records": sum(row["dirac_matched"] == "1" for row in records),
        "gco_unmatched_count": len(gco_unmatched),
        "dirac_unmatched_count": len(dirac_unmatched),
        "gco_unmatched_iso3": gco_unmatched,
        "dirac_unmatched_iso3": dirac_unmatched,
    }
    return records, merge_summary


def add_mismatch_metrics(records: list[dict[str, Any]]) -> dict[str, Any]:
    complete = [
        row
        for row in records
        if row["dirac_matched"] == "1"
        and not math.isnan(to_float(row.get("selected_site_relative_case_growth_2050")))
        and not math.isnan(to_float(row.get("selected_site_case_increase_2050")))
        and not math.isnan(to_float(row.get("mv_units_per_1000_selected_site_cases_2050")))
    ]
    growth_values = [to_float(row["selected_site_relative_case_growth_2050"]) for row in complete]
    resource_values = [to_float(row["mv_units_per_1000_selected_site_cases_2050"]) for row in complete]
    growth_q75 = quantile(growth_values, 0.75)
    resource_q25 = quantile(resource_values, 0.25)

    total_increment = sum(to_float(row["selected_site_case_increase_2050"]) for row in records if not math.isnan(to_float(row["selected_site_case_increase_2050"])))

    for row in records:
        rel_growth = to_float(row.get("selected_site_relative_case_growth_2050"))
        resource = to_float(row.get("mv_units_per_1000_selected_site_cases_2050"))
        abs_increase = to_float(row.get("selected_site_case_increase_2050"))
        row["high_growth_q75"] = "1" if not math.isnan(rel_growth) and rel_growth >= growth_q75 else "0"
        row["lower_quartile_unit_density"] = "1" if not math.isnan(resource) and resource <= resource_q25 else "0"
        row["acceleration_risk"] = "1" if row["high_growth_q75"] == "1" and row["lower_quartile_unit_density"] == "1" else "0"
        row["share_global_selected_site_case_increment"] = safe_div(abs_increase, total_increment)

    pressure_ranked = sorted(
        [
            row
            for row in records
            if row["dirac_matched"] == "1"
            and not math.isnan(to_float(row.get("selected_site_cases_per_mv_unit_2050")))
        ],
        key=lambda row: to_float(row["selected_site_cases_per_mv_unit_2050"]),
        reverse=True,
    )
    for idx, row in enumerate(pressure_ranked, start=1):
        row["unit_pressure_rank"] = idx
    for row in records:
        row.setdefault("unit_pressure_rank", "")

    return {
        "complete_case_records": len(complete),
        "growth_q75": growth_q75,
        "resource_q25": resource_q25,
        "global_rt_relevant_case_increment": total_increment,
        "acceleration_risk_count": sum(row["acceleration_risk"] == "1" for row in records),
    }


def region_summary(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in records:
        grouped[row.get("who_region", "") or "Missing"].append(row)
    out: list[dict[str, Any]] = []
    for region, rows in sorted(grouped.items()):
        matched = [row for row in rows if row["dirac_matched"] == "1"]
        mismatch = [row for row in rows if row["acceleration_risk"] == "1"]
        out.append(
            {
                "who_region": region,
                "gco_countries": len(rows),
                "dirac_matched_countries": len(matched),
                "acceleration_risk_countries": len(mismatch),
                "selected_site_cases_2050": sum(to_float(row["selected_site_cases_2050"]) for row in rows if not math.isnan(to_float(row["selected_site_cases_2050"]))),
                "selected_site_case_increment_2050": sum(to_float(row["selected_site_case_increase_2050"]) for row in rows if not math.isnan(to_float(row["selected_site_case_increase_2050"]))),
                "mv_therapy_units": sum(to_float(row["mv_therapy_units"]) for row in matched if not math.isnan(to_float(row["mv_therapy_units"]))),
            }
        )
    return out


def write_report(summary: dict[str, Any], top_rows: list[dict[str, Any]], paths: dict[str, Path]) -> None:
    lines = [
        "# Radiotherapy Resource Mismatch Dataset QC",
        "",
        f"Run finished UTC: {datetime.now(timezone.utc).isoformat()}",
        "",
        "## Inputs",
        "",
        f"- GCO Cancer Tomorrow API data version 2024: `{GCO_PREDICTIONS.relative_to(PROJECT_ROOT)}`",
        f"- DIRAC country resources: `{DIRAC_RESOURCES.relative_to(PROJECT_ROOT)}`",
        "",
        "## Merge Summary",
        "",
        f"- GCO country records: {summary['merge']['gco_country_records']}",
        f"- DIRAC country records: {summary['merge']['dirac_country_records']}",
        f"- Matched records: {summary['merge']['matched_records']}",
        f"- GCO records without DIRAC match: {summary['merge']['gco_unmatched_count']}",
        f"- DIRAC records without GCO match: {summary['merge']['dirac_unmatched_count']}",
        "",
        "## Acceleration-Risk Thresholds",
        "",
        f"- Complete-case records: {summary['metrics']['complete_case_records']}",
        f"- High-growth threshold, relative selected-site case growth q75: {summary['metrics']['growth_q75']:.4f}",
        f"- Lower-quartile unit-density threshold, MV units per 1000 selected-site 2050 cases q25: {summary['metrics']['resource_q25']:.4f}",
        f"- Acceleration-risk countries: {summary['metrics']['acceleration_risk_count']}",
        "",
        "## Interpretation Notes",
        "",
        "- Countries absent from the DIRAC country table are treated as missing resource observations without assigned resource-density values.",
        "- The primary mismatch analysis is therefore a complete-case comparison among GLOBOCAN country records with a DIRAC ISO3 match.",
        "",
        "## Top 10 Countries by Selected-Site Cases per MV Unit",
        "",
        "| Rank | Country | ISO3 | WHO region | Income | Selected-site cases 2050 | Relative growth | MV units | Selected-site cases per MV unit |",
        "|---:|---|---|---|---|---:|---:|---:|---:|",
    ]
    for row in top_rows[:10]:
        lines.append(
                "| {rank} | {country} | {iso3} | {region} | {income} | {cases:.0f} | {growth:.3f} | {units:.0f} | {density:.0f} |".format(
                rank=row.get("unit_pressure_rank", ""),
                country=row.get("gco_country", ""),
                iso3=row.get("country_iso3", ""),
                region=row.get("who_region", ""),
                income=row.get("income_label", ""),
                cases=to_float(row.get("selected_site_cases_2050")),
                growth=to_float(row.get("selected_site_relative_case_growth_2050")),
                units=to_float(row.get("mv_therapy_units")),
                density=to_float(row.get("selected_site_cases_per_mv_unit_2050")),
            )
        )
    lines.extend(
        [
            "",
            "## Output Files",
            "",
        ]
    )
    for label, path in paths.items():
        lines.append(f"- {label}: `{path.relative_to(PROJECT_ROOT)}`")
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    (REPORT_DIR / "report_06_radiotherapy_resource_mismatch_qc.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    for directory in [DATA_PROCESSED, TABLE_DIR, REPORT_DIR, PROTOCOL_DIR, LOG_DIR]:
        directory.mkdir(parents=True, exist_ok=True)
    write_site_rules()

    predictions = read_csv(GCO_PREDICTIONS)
    dirac = read_csv(DIRAC_RESOURCES)
    burden, site_summary = build_burden_tables(predictions)
    records, merge = merge_with_dirac(burden, dirac)
    metrics = add_mismatch_metrics(records)

    records_sorted = sorted(
        [row for row in records if not math.isnan(to_float(row.get("selected_site_cases_per_mv_unit_2050")))],
        key=lambda row: to_float(row.get("selected_site_cases_per_mv_unit_2050")),
        reverse=True,
    )
    top_rows = records_sorted
    region_rows = region_summary(records)
    site_rows = sorted(site_summary.values(), key=lambda row: to_float(row["world_cases_2050"]), reverse=True)

    output_main = DATA_PROCESSED / "radiotherapy_resource_mismatch_country_2024_2050.csv"
    output_top = TABLE_DIR / "table_09_unit_pressure_top_countries.csv"
    output_region = TABLE_DIR / "table_10_region_summary.csv"
    output_site = TABLE_DIR / "table_11_selected_site_burden.csv"
    output_log = LOG_DIR / "radiotherapy_resource_mismatch_build_log.json"

    write_csv(output_main, records)
    top_fields = [
        "unit_pressure_rank",
        "country_iso3",
        "gco_country",
        "who_region",
        "hdi_label",
        "income_label",
        "selected_site_cases_2024",
        "selected_site_cases_2050",
        "selected_site_case_increase_2050",
        "selected_site_relative_case_growth_2050",
        "rt_centres_with_rt",
        "mv_therapy_units",
        "brachytherapy_units",
        "selected_site_cases_per_mv_unit_2050",
        "mv_units_per_1000_selected_site_cases_2050",
        "rt_centres_per_10000_selected_site_cases_2050",
        "acceleration_risk",
        "dirac_last_update_year",
    ]
    write_csv(output_top, top_rows[:30], top_fields)
    write_csv(output_region, region_rows)
    write_csv(output_site, site_rows)

    summary = {
        "script": str(Path(__file__).resolve()),
        "run_finished_utc": datetime.now(timezone.utc).isoformat(),
        "merge": merge,
        "metrics": metrics,
        "rt_site_codes": [row["cancer_code"] for row in RT_SITE_RULES],
        "outputs": {
            "main_country_dataset": str(output_main),
            "top_countries_table": str(output_top),
            "region_summary_table": str(output_region),
            "site_burden_table": str(output_site),
            "site_rules": str(PROTOCOL_DIR / "selected_cancer_site_rules.csv"),
        },
    }
    output_log.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    write_report(
        summary,
        top_rows,
        {
            "country dataset": output_main,
            "top countries": output_top,
            "region summary": output_region,
            "site burden": output_site,
            "site rules": PROTOCOL_DIR / "selected_cancer_site_rules.csv",
            "build log": output_log,
        },
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
