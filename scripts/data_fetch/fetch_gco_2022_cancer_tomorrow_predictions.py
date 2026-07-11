"""Fetch GCO Cancer Tomorrow 2022 projections for burden-version comparison.

This wrapper reuses the verified 2024 extraction implementation while changing
only version-specific API, metadata, cache, and output paths.
"""

from __future__ import annotations

import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.data_fetch import fetch_gco_2024_cancer_tomorrow_predictions as core  # noqa: E402


core.GLOBOCAN_VERSION = "2022"
core.POPULATION_CSV = PROJECT_ROOT / "data" / "source" / "gco_populations_2022_dictionary.csv"
core.RAW_DIR = PROJECT_ROOT / "data" / "source" / "raw" / "gco_cancer_tomorrow_2022"
core.BASE_URL = "https://gco.iarc.who.int/gateway_prod/api/globocan/v3/2022"
core.OUTPUT_CSV = PROJECT_ROOT / "data" / "source" / "gco_cancer_tomorrow_2022_predictions_long.csv"
core.FETCH_LOG = PROJECT_ROOT / "data" / "source" / "logs" / "gco_2022_cancer_tomorrow_predictions_fetch_log.json"


if __name__ == "__main__":
    raise SystemExit(core.main())
