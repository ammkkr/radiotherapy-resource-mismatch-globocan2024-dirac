"""Fetch GCO / Cancer Tomorrow 2024-2050 projections.

Main extraction:
  - type 0 incidence and type 1 mortality
  - sex 0 both sexes
  - all GCO population codes from the downloaded population dictionary
  - main analysis cancer sites plus all-cancer denominator codes 39 and 40
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


PROJECT_ROOT = Path(__file__).resolve().parents[2]
POPULATION_CSV = PROJECT_ROOT / "data" / "interim" / "gco_populations_2024_dictionary.csv"
RULES_CSV = PROJECT_ROOT / "data" / "interim" / "gco_2024_cancer_site_analysis_rules.csv"
RAW_DIR = PROJECT_ROOT / "data" / "raw" / "gco_cancer_tomorrow_2024"
INTERIM_DIR = PROJECT_ROOT / "data" / "interim"
LOG_DIR = PROJECT_ROOT / "logs"

BASE_URL = "https://gco.iarc.who.int/gateway_prod/api/globocan/v3/2024"
GLOBOCAN_VERSION = "2024"
TYPE_MAP = {"0": "incidence", "1": "mortality"}
SEX_MAP = {"0": "both"}
DENOMINATOR_CODES = {"39", "40"}
OUTPUT_CSV = INTERIM_DIR / "gco_cancer_tomorrow_2024_predictions_long.csv"
FETCH_LOG = LOG_DIR / "gco_2024_cancer_tomorrow_predictions_fetch_log.json"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def chunks(values: list[str], size: int) -> list[list[str]]:
    return [values[i : i + size] for i in range(0, len(values), size)]


def fetch_json(url: str) -> tuple[Any, dict[str, Any]]:
    request = Request(
        url,
        headers={
            "Accept": "application/json",
            "User-Agent": f"globocan-{GLOBOCAN_VERSION}-mismatch-data-prep/0.1",
        },
    )
    started = datetime.now(timezone.utc)
    try:
        with urlopen(request, timeout=120) as response:
            raw = response.read()
            status = getattr(response, "status", None)
            content_type = response.headers.get("Content-Type")
    except HTTPError as exc:
        raise RuntimeError(f"HTTP error {exc.code} for {url}") from exc
    except URLError as exc:
        raise RuntimeError(f"URL error for {url}: {exc.reason}") from exc

    payload = json.loads(raw.decode("utf-8"))
    meta = {
        "url": url,
        "status": status,
        "content_type": content_type,
        "bytes": len(raw),
        "started_utc": started.isoformat(),
        "finished_utc": datetime.now(timezone.utc).isoformat(),
    }
    return payload, meta


def endpoint(type_code: str, sex_code: str, population_codes: list[str], cancer_codes: list[str]) -> str:
    pop_part = "_".join(population_codes)
    cancer_part = "_".join(cancer_codes)
    return (
        f"{BASE_URL}/data/prediction/{type_code}/{sex_code}/{pop_part}/{cancer_part}/"
        "?prediction_annual=1&ages_group=0_17"
    )


def flatten_projection_row(
    row: dict[str, Any],
    population_lookup: dict[str, dict[str, str]],
    cancer_lookup: dict[str, dict[str, str]],
) -> dict[str, Any]:
    population_code = str(row.get("id", ""))
    cancer_code = str(row.get("cancer", ""))
    type_code = str(row.get("type", ""))
    sex_code = str(row.get("sex", ""))
    population_meta = population_lookup.get(population_code, {})
    cancer_meta = cancer_lookup.get(cancer_code, {})
    population_component = row.get("population", {}) if isinstance(row.get("population"), dict) else {}
    risk_component = row.get("risk", {}) if isinstance(row.get("risk"), dict) else {}

    return {
        "country_code": population_code,
        "population_label": row.get("id_label", population_meta.get("label", "")),
        "country_iso3": population_meta.get("country_iso3", ""),
        "grouping": population_meta.get("grouping", ""),
        "area_label": population_meta.get("area_label", ""),
        "who_region": population_meta.get("who_region", ""),
        "hdi_label": population_meta.get("hdi_label", ""),
        "income_label": population_meta.get("income_label", ""),
        "sex": sex_code,
        "sex_label": SEX_MAP.get(sex_code, sex_code),
        "type": type_code,
        "measure": TYPE_MAP.get(type_code, type_code),
        "cancer_code": cancer_code,
        "cancer_label": row.get("cancer_label", cancer_meta.get("label", "")),
        "is_main_site": cancer_meta.get("is_main_site", ""),
        "is_summary_denominator": cancer_meta.get("is_summary_denominator", ""),
        "year": row.get("year", ""),
        "predicted_count": row.get("cases_pred", ""),
        "cases_pred": row.get("cases_pred", ""),
        "cases_base": row.get("cases_base", ""),
        "change": row.get("change", ""),
        "percent": row.get("percent", ""),
        "pop": row.get("pop", ""),
        "APC": row.get("APC", ""),
        "population_cases": population_component.get("cases", ""),
        "population_change": population_component.get("change", ""),
        "population_percent": population_component.get("percent", ""),
        "risk_change": risk_component.get("change", ""),
        "risk_percent": risk_component.get("percent", ""),
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> list[str]:
    fieldnames = sorted({key for row in rows for key in row})
    preferred = [
        "country_code",
        "population_label",
        "country_iso3",
        "grouping",
        "area_label",
        "who_region",
        "hdi_label",
        "income_label",
        "sex",
        "sex_label",
        "type",
        "measure",
        "cancer_code",
        "cancer_label",
        "is_main_site",
        "is_summary_denominator",
        "year",
        "predicted_count",
        "cases_base",
        "change",
        "percent",
        "pop",
        "population_cases",
        "population_change",
        "population_percent",
        "risk_change",
        "risk_percent",
        "APC",
    ]
    ordered = [x for x in preferred if x in fieldnames] + [x for x in fieldnames if x not in preferred]
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=ordered)
        writer.writeheader()
        writer.writerows(rows)
    return ordered


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=1,
        help="Population codes per request. Default 1 avoids silent under-return from large multi-population requests.",
    )
    parser.add_argument("--sleep", type=float, default=0.10)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--limit-populations", type=int, default=None)
    args = parser.parse_args()

    RAW_DIR.mkdir(parents=True, exist_ok=True)
    INTERIM_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    populations = read_csv(POPULATION_CSV)
    if args.limit_populations is not None:
        populations = populations[: args.limit_populations]
    rules = read_csv(RULES_CSV)

    population_lookup = {row["country_code"]: row for row in populations}
    cancer_lookup = {row["cancer_code"]: row for row in rules}
    population_codes = [row["country_code"] for row in populations]
    cancer_codes = [
        row["cancer_code"]
        for row in rules
        if row["is_main_site"] == "True" or row["cancer_code"] in DENOMINATOR_CODES
    ]
    cancer_codes = sorted(set(cancer_codes), key=lambda x: int(x))

    all_rows: list[dict[str, Any]] = []
    successes: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []
    pop_chunks = chunks(population_codes, args.chunk_size)

    for type_code in ["0", "1"]:
        for chunk_index, pop_chunk in enumerate(pop_chunks, start=1):
            url = endpoint(type_code, "0", pop_chunk, cancer_codes)
            raw_path = RAW_DIR / f"prediction_type{type_code}_sex0_popchunk{chunk_index:02d}.json"
            try:
                if raw_path.exists() and not args.force:
                    payload = json.loads(raw_path.read_text(encoding="utf-8"))
                    meta = {
                        "url": url,
                        "status": "cached",
                        "bytes": raw_path.stat().st_size,
                        "started_utc": None,
                        "finished_utc": None,
                    }
                else:
                    payload, meta = fetch_json(url)
                    raw_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
                    time.sleep(args.sleep)

                dataset = payload.get("dataset", []) if isinstance(payload, dict) else []
                if not isinstance(dataset, list):
                    raise RuntimeError("Payload did not contain a list-valued dataset.")

                for row in dataset:
                    if isinstance(row, dict):
                        all_rows.append(flatten_projection_row(row, population_lookup, cancer_lookup))

                successes.append(
                    {
                        "type": type_code,
                        "chunk": chunk_index,
                        "population_count": len(pop_chunk),
                        "record_count": len(dataset),
                        "raw_file": str(raw_path.relative_to(PROJECT_ROOT)),
                        "meta": meta,
                    }
                )
                print(
                    f"OK type={type_code} chunk={chunk_index}/{len(pop_chunks)} "
                    f"populations={len(pop_chunk)} rows={len(dataset)}"
                )
            except Exception as exc:  # noqa: BLE001 - log and continue.
                failures.append(
                    {
                        "type": type_code,
                        "chunk": chunk_index,
                        "population_codes": pop_chunk,
                        "url": url,
                        "error": str(exc),
                    }
                )
                print(f"FAIL type={type_code} chunk={chunk_index}: {exc}")

    out_csv = OUTPUT_CSV
    fields = write_csv(out_csv, all_rows)
    expected_rows = len(population_codes) * len(cancer_codes) * 2 * 7

    run_log = {
        "script": str(Path(__file__).resolve()),
        "run_finished_utc": datetime.now(timezone.utc).isoformat(),
        "globocan_version": GLOBOCAN_VERSION,
        "api_base": BASE_URL,
        "parameters": {
            "types": ["0", "1"],
            "sex": "0",
            "population_codes": len(population_codes),
            "cancer_codes": cancer_codes,
            "prediction_annual": 1,
            "ages_group": "0_17",
            "chunk_size": args.chunk_size,
        },
        "chunks_requested": len(pop_chunks) * 2,
        "chunks_successful": len(successes),
        "chunks_failed": len(failures),
        "rows_expected": expected_rows,
        "rows_written": len(all_rows),
        "fields": fields,
        "output_csv": str(out_csv.relative_to(PROJECT_ROOT)),
        "successes": successes,
        "failures": failures,
    }
    log_path = FETCH_LOG
    log_path.write_text(json.dumps(run_log, ensure_ascii=False, indent=2), encoding="utf-8")

    print(
        json.dumps(
            {
                "chunks_requested": run_log["chunks_requested"],
                "chunks_successful": run_log["chunks_successful"],
                "chunks_failed": run_log["chunks_failed"],
                "rows_expected": expected_rows,
                "rows_written": len(all_rows),
                "output_csv": run_log["output_csv"],
            },
            ensure_ascii=False,
            indent=2,
        )
    )
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
