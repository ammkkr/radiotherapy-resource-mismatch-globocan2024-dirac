"""Fetch country-level radiotherapy resource data from IAEA DIRAC.

Outputs:
  - raw JSON under data/source/raw/dirac/
  - parsed country resource CSV under data/source/
  - a fetch log under data/source/logs/
  - a short codebook under docs/
"""

from __future__ import annotations

import csv
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = PROJECT_ROOT / "data" / "source" / "raw" / "dirac"
INTERIM_DIR = PROJECT_ROOT / "data" / "source"
CODEBOOK_DIR = PROJECT_ROOT / "docs"
LOG_DIR = PROJECT_ROOT / "data" / "source" / "logs"

DIRAC_COUNTRIES_URL = "https://dirac.iaea.org/api/DataGridWebApi/GetCountriesAndRegions"

KEEP_FIELDS = [
    "CountryISO3",
    "Country",
    "RegionName",
    "LastUpdate",
    "RTCenters",
    "RTCentersWithRT",
    "HePhotonAndElectronBeamRt",
    "ProtonIonTherapy",
    "XRayGenerator",
    "BrachyTherapy",
    "BrachyTherapyIncEl",
    "CT",
    "Simulators",
    "TPS",
]

RENAME_FIELDS = {
    "CountryISO3": "dirac_iso3",
    "Country": "dirac_country",
    "RegionName": "dirac_region",
    "LastUpdate": "dirac_last_update_year",
    "RTCenters": "rt_centres_total",
    "RTCentersWithRT": "rt_centres_with_rt",
    "HePhotonAndElectronBeamRt": "mv_therapy_units",
    "ProtonIonTherapy": "proton_ion_units",
    "XRayGenerator": "kv_therapy_units",
    "BrachyTherapy": "brachytherapy_flag",
    "BrachyTherapyIncEl": "brachytherapy_units",
    "CT": "ct_units",
    "Simulators": "simulators",
    "TPS": "treatment_planning_systems",
}

NUMERIC_FIELDS = [
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
]


def fetch_json(url: str) -> tuple[dict[str, Any], dict[str, Any]]:
    request = Request(
        url,
        headers={
            "Accept": "application/json",
            "User-Agent": "globocan-rt-resource-mismatch/0.1",
        },
    )
    started = datetime.now(timezone.utc)
    try:
        with urlopen(request, timeout=90) as response:
            raw = response.read()
            status = getattr(response, "status", None)
            content_type = response.headers.get("Content-Type")
    except HTTPError as exc:
        raise RuntimeError(f"HTTP error {exc.code} for {url}") from exc
    except URLError as exc:
        raise RuntimeError(f"URL error for {url}: {exc.reason}") from exc

    payload = json.loads(raw.decode("utf-8"))
    if not isinstance(payload, dict):
        raise RuntimeError(f"Expected JSON object from DIRAC, got {type(payload).__name__}")
    meta = {
        "url": url,
        "status": status,
        "content_type": content_type,
        "bytes": len(raw),
        "started_utc": started.isoformat(),
        "finished_utc": datetime.now(timezone.utc).isoformat(),
    }
    return payload, meta


def clean_record(record: dict[str, Any]) -> dict[str, Any]:
    cleaned: dict[str, Any] = {}
    for field in KEEP_FIELDS:
        out_field = RENAME_FIELDS[field]
        value = record.get(field)
        if out_field == "dirac_iso3" and isinstance(value, str):
            value = value.strip().upper()
        elif isinstance(value, str):
            value = value.strip()
        cleaned[out_field] = value

    for field in NUMERIC_FIELDS:
        value = cleaned.get(field)
        if value is None or value == "":
            cleaned[field] = ""
            continue
        try:
            if isinstance(value, float) and value.is_integer():
                cleaned[field] = int(value)
            else:
                cleaned[field] = int(value)
        except (TypeError, ValueError):
            cleaned[field] = value

    return cleaned


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def write_csv(path: Path, records: list[dict[str, Any]]) -> None:
    fieldnames = [RENAME_FIELDS[field] for field in KEEP_FIELDS]
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)


def summarize(records: list[dict[str, Any]]) -> dict[str, Any]:
    duplicate_iso3 = sorted(
        {
            row["dirac_iso3"]
            for row in records
            if row["dirac_iso3"]
            and sum(1 for x in records if x["dirac_iso3"] == row["dirac_iso3"]) > 1
        }
    )
    summary = {
        "records": len(records),
        "countries_with_iso3": sum(1 for row in records if row.get("dirac_iso3")),
        "duplicate_iso3": duplicate_iso3,
        "latest_update_year": max(
            int(row["dirac_last_update_year"])
            for row in records
            if str(row.get("dirac_last_update_year", "")).isdigit()
        ),
        "earliest_update_year": min(
            int(row["dirac_last_update_year"])
            for row in records
            if str(row.get("dirac_last_update_year", "")).isdigit()
        ),
        "sum_rt_centres_with_rt": sum(int(row["rt_centres_with_rt"] or 0) for row in records),
        "sum_mv_therapy_units": sum(int(row["mv_therapy_units"] or 0) for row in records),
        "sum_brachytherapy_units": sum(int(row["brachytherapy_units"] or 0) for row in records),
        "countries_with_no_reported_mv_units": sum(int(row["mv_therapy_units"] or 0) == 0 for row in records),
    }
    return summary


def write_codebook(path: Path, summary: dict[str, Any], raw_path: Path, csv_path: Path) -> None:
    lines = [
        "# DIRAC Country Radiotherapy Resource Codebook",
        "",
        "Source:",
        "",
        f"- Endpoint: `{DIRAC_COUNTRIES_URL}`",
        f"- Raw JSON: `{raw_path.relative_to(PROJECT_ROOT)}`",
        f"- Parsed CSV: `{csv_path.relative_to(PROJECT_ROOT)}`",
        "- Citation wording suggested by DIRAC: The IAEA Directory of Radiotherapy Centres (DIRAC), https://dirac.iaea.org/ (accessed on mm/yyyy).",
        "",
        "Fetch scope:",
        "",
        f"- Records: {summary['records']}",
        f"- Records with ISO3: {summary['countries_with_iso3']}",
        f"- Update years: {summary['earliest_update_year']} to {summary['latest_update_year']}",
        f"- Sum RT centres with RT: {summary['sum_rt_centres_with_rt']}",
        f"- Sum MV therapy units: {summary['sum_mv_therapy_units']}",
        f"- Sum brachytherapy units: {summary['sum_brachytherapy_units']}",
        f"- Countries with no reported MV therapy units: {summary['countries_with_no_reported_mv_units']}",
        "",
        "Variables:",
        "",
        "| Variable | Source field | Description |",
        "|---|---|---|",
        "| `dirac_iso3` | CountryISO3 | ISO3 country code in DIRAC |",
        "| `dirac_country` | Country | DIRAC country name |",
        "| `dirac_region` | RegionName | DIRAC region name |",
        "| `dirac_last_update_year` | LastUpdate | Latest DIRAC update year shown for country record |",
        "| `rt_centres_total` | RTCenters | Total radiotherapy centres in DIRAC country record |",
        "| `rt_centres_with_rt` | RTCentersWithRT | Centres with radiotherapy |",
        "| `mv_therapy_units` | HePhotonAndElectronBeamRt | Megavoltage photon/electron therapy units shown in country table |",
        "| `proton_ion_units` | ProtonIonTherapy | Proton/light-ion therapy units |",
        "| `kv_therapy_units` | XRayGenerator | kV therapy units |",
        "| `brachytherapy_flag` | BrachyTherapy | DIRAC brachytherapy text flag |",
        "| `brachytherapy_units` | BrachyTherapyIncEl | Brachytherapy units including electronic units |",
        "| `ct_units` | CT | CT units listed in DIRAC country record |",
        "| `simulators` | Simulators | Simulator units listed in DIRAC country record |",
        "| `treatment_planning_systems` | TPS | Treatment planning systems listed in DIRAC country record |",
        "",
        "Notes:",
        "",
        "- `mv_therapy_units` is the preferred primary resource denominator for the first mismatch analysis.",
        "- `rt_centres_with_rt` and `brachytherapy_units` should be used as sensitivity resource variables.",
        "- Zero resource values should not be assumed for countries absent from DIRAC; absence is handled as missing after merging with GLOBOCAN.",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--from-json",
        type=Path,
        default=None,
        help="Process an already downloaded DIRAC JSON file instead of fetching online.",
    )
    args = parser.parse_args()

    for directory in [RAW_DIR, INTERIM_DIR, CODEBOOK_DIR, LOG_DIR]:
        directory.mkdir(parents=True, exist_ok=True)

    raw_path = RAW_DIR / "dirac_countries_and_regions_raw.json"
    csv_path = INTERIM_DIR / "dirac_country_resources.csv"
    log_path = LOG_DIR / "dirac_country_resources_fetch_log.json"
    codebook_path = CODEBOOK_DIR / "dirac_country_resources_codebook.md"

    if args.from_json is not None:
        payload = json.loads(args.from_json.read_text(encoding="utf-8-sig"))
        meta = {
            "url": DIRAC_COUNTRIES_URL,
            "status": "local_file",
            "content_type": "application/json",
            "bytes": args.from_json.stat().st_size,
            "started_utc": datetime.now(timezone.utc).isoformat(),
            "finished_utc": datetime.now(timezone.utc).isoformat(),
            "local_file": args.from_json.name,
        }
        if args.from_json.resolve() != raw_path.resolve():
            write_json(raw_path, payload)
    else:
        payload, meta = fetch_json(DIRAC_COUNTRIES_URL)
        write_json(raw_path, payload)

    records_raw = payload.get("data")
    if not isinstance(records_raw, list):
        raise RuntimeError("DIRAC payload does not contain a list under key 'data'.")

    records = [clean_record(row) for row in records_raw if isinstance(row, dict)]
    write_csv(csv_path, records)

    summary = summarize(records)
    log = {
        "script": "scripts/fetch_dirac_country_resources.py",
        "run_started_utc": meta["started_utc"],
        "run_finished_utc": datetime.now(timezone.utc).isoformat(),
        "source": meta,
        "summary": summary,
        "outputs": {
            "raw_json": raw_path.relative_to(PROJECT_ROOT).as_posix(),
            "csv": csv_path.relative_to(PROJECT_ROOT).as_posix(),
            "codebook": codebook_path.relative_to(PROJECT_ROOT).as_posix(),
        },
    }
    write_json(log_path, log)
    write_codebook(codebook_path, summary, raw_path, csv_path)

    print(json.dumps({"summary": summary, "outputs": log["outputs"]}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
