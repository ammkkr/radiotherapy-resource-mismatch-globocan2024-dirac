"""Create a 3-panel Nature-style figure for selected-site/MV-unit mismatch.

Backend: Python/matplotlib only, following the nature-figure contract.
"""

from __future__ import annotations

import json
import math
import textwrap
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.collections import PatchCollection
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Polygon


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATASET = PROJECT_ROOT / "data" / "radiotherapy_resource_mismatch_country_2024_2050.csv"
BUILD_LOG = PROJECT_ROOT / "data" / "build_log_sanitized.json"
MAP_GEOJSON = PROJECT_ROOT / "data" / "natural_earth" / "ne_110m_admin_0_countries.geojson"
OUT_DIR = PROJECT_ROOT / "figures"
DOC_DIR = PROJECT_ROOT / "docs"
COUNTRY_SOURCE_DATA = PROJECT_ROOT / "data" / "figure_1_country_source_data.csv"
REGION_SOURCE_DATA = PROJECT_ROOT / "data" / "figure_1_region_source_data.csv"


PALETTE = {
    "blue_main": "#0F4D92",
    "blue_mid": "#6F91C7",
    "signal": "#B64342",
    "signal_dark": "#7F2A18",
    "signal_soft": "#F6CFCB",
    "neutral_pale": "#F4F4F4",
    "neutral_light": "#D8D8D8",
    "neutral_mid": "#767676",
    "neutral_dark": "#4D4D4D",
    "black": "#272727",
    "missing": "#E9E9E9",
    "not_included": "#F7F7F7",
}

REGION_COLORS = {
    "AFRO": "#B64342",
    "EMRO": "#C97B48",
    "EURO": "#6F91C7",
    "PAHO": "#42949E",
    "SEARO": "#9A4D8E",
    "WPRO": "#0F4D92",
}

COUNTRY_SHORT = {
    "Congo, Democratic Republic of": "DR Congo",
    "Tanzania, United Republic of": "Tanzania",
    "Bolivia (Plurinational State of)": "Bolivia",
    "Korea, Democratic People's Republic of": "DPR Korea",
    "Democratic Republic of the Congo": "DR Congo",
    "Tanzania, United Republic of": "Tanzania",
    "CÃƒÂ´te d'Ivoire": "Cote d'Ivoire",
    "CÃ´te d'Ivoire": "Cote d'Ivoire",
}

LABEL_OFFSETS = {
    "ETH": (20, 8),
    "COD": (-44, -8),
    "YEM": (30, 8),
    "ZWE": (25, 6),
    "MOZ": (34, -20),
    "PRK": (34, 6),
    "NER": (-32, 10),
    "MWI": (28, -8),
    "AGO": (-40, -12),
    "NGA": (-35, 13),
    "UGA": (22, 2),
    "ZMB": (30, -15),
    "TZA": (28, 6),
}


def apply_publication_style() -> None:
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
    plt.rcParams["svg.fonttype"] = "none"
    plt.rcParams["pdf.fonttype"] = 42
    mpl.rcParams.update(
        {
            "font.size": 6.6,
            "font.weight": "bold",
            "axes.labelsize": 6.7,
            "axes.titlesize": 7.2,
            "axes.labelweight": "bold",
            "axes.titleweight": "bold",
            "xtick.labelsize": 5.8,
            "ytick.labelsize": 5.8,
            "legend.fontsize": 5.7,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.55,
            "xtick.major.width": 0.45,
            "ytick.major.width": 0.45,
            "xtick.major.size": 2.0,
            "ytick.major.size": 2.0,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def to_num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def short_country(name: str) -> str:
    return COUNTRY_SHORT.get(name, name)


def wrap_country_label(rank: int, iso3: str, name: str) -> str:
    label = short_country(name)
    if len(label) > 17:
        label = "\n".join(textwrap.wrap(label, width=17, break_long_words=False))
    return f"{rank}. {label}"


def add_panel_label(ax: plt.Axes, label: str, x: float = -0.035, y: float = 1.02) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8.0,
        fontweight="bold",
        color=PALETTE["black"],
    )


def bolden_figure_text(fig: plt.Figure) -> None:
    for text in fig.findobj(match=mpl.text.Text):
        text.set_fontweight("bold")


def feature_iso3(properties: dict[str, Any]) -> str | None:
    for key in ["ISO_A3", "ADM0_A3", "SOV_A3", "GU_A3", "WB_A3", "ISO_A3_EH"]:
        value = properties.get(key)
        if isinstance(value, str) and value and value != "-99":
            return value.upper()
    return None


def polygon_exteriors(geometry: dict[str, Any]) -> list[list[tuple[float, float]]]:
    geom_type = geometry.get("type")
    coords = geometry.get("coordinates")
    exteriors: list[list[tuple[float, float]]] = []
    if geom_type == "Polygon" and isinstance(coords, list):
        if coords and isinstance(coords[0], list):
            exteriors.append([(float(x), float(y)) for x, y in coords[0]])
    elif geom_type == "MultiPolygon" and isinstance(coords, list):
        for poly in coords:
            if poly and isinstance(poly[0], list):
                exteriors.append([(float(x), float(y)) for x, y in poly[0]])
    return exteriors


def polygon_area_and_centroid(points: list[tuple[float, float]]) -> tuple[float, tuple[float, float]]:
    if len(points) < 3:
        xs = [p[0] for p in points]
        ys = [p[1] for p in points]
        return 0.0, (float(np.mean(xs)), float(np.mean(ys)))
    area2 = 0.0
    cx_num = 0.0
    cy_num = 0.0
    pts = points if points[0] == points[-1] else points + [points[0]]
    for (x0, y0), (x1, y1) in zip(pts[:-1], pts[1:]):
        cross = x0 * y1 - x1 * y0
        area2 += cross
        cx_num += (x0 + x1) * cross
        cy_num += (y0 + y1) * cross
    if abs(area2) < 1e-9:
        xs = [p[0] for p in points]
        ys = [p[1] for p in points]
        return 0.0, (float(np.mean(xs)), float(np.mean(ys)))
    return abs(area2 / 2.0), (cx_num / (3.0 * area2), cy_num / (3.0 * area2))


def load_data() -> tuple[pd.DataFrame, dict[str, Any]]:
    df = pd.read_csv(DATASET)
    numeric_cols = [
        "selected_site_cases_2024",
        "selected_site_cases_2050",
        "selected_site_case_increase_2050",
        "selected_site_relative_case_growth_2050",
        "mv_therapy_units",
        "mv_units_per_1000_selected_site_cases_2050",
    ]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = to_num(df[col])
    df["selected_site_cases_per_mv_unit_2050"] = np.where(
        (df["dirac_matched"].astype(str).eq("1")) & (df["mv_therapy_units"] > 0),
        df["selected_site_cases_2050"] / df["mv_therapy_units"],
        np.nan,
    )
    df["unit_pressure_rank"] = np.nan
    pressure_index = (
        df.dropna(subset=["selected_site_cases_per_mv_unit_2050"])
        .sort_values("selected_site_cases_per_mv_unit_2050", ascending=False)
        .index
    )
    df.loc[pressure_index, "unit_pressure_rank"] = np.arange(1, len(pressure_index) + 1)
    with BUILD_LOG.open("r", encoding="utf-8") as handle:
        log = json.load(handle)
    return df, log


def draw_world_map(
    ax: plt.Axes,
    df: pd.DataFrame,
    cmap: LinearSegmentedColormap,
    norm: LogNorm,
) -> dict[str, tuple[float, float]]:
    if not MAP_GEOJSON.exists():
        raise FileNotFoundError(f"Missing Natural Earth GeoJSON: {MAP_GEOJSON}")

    by_iso = df.set_index("country_iso3", drop=False)
    geo = json.loads(MAP_GEOJSON.read_text(encoding="utf-8"))

    value_patches: list[Polygon] = []
    value_facecolors: list[Any] = []
    missing_patches: list[Polygon] = []
    not_included_patches: list[Polygon] = []
    highlight_patches: list[Polygon] = []
    centroids: dict[str, tuple[float, float]] = {}

    for feature in geo.get("features", []):
        props = feature.get("properties", {})
        iso3 = feature_iso3(props)
        if iso3 == "ATA":
            continue
        in_data = iso3 in by_iso.index if iso3 is not None else False
        row = by_iso.loc[iso3] if in_data else None
        if isinstance(row, pd.DataFrame):
            row = row.iloc[0]

        if row is not None and str(row.get("dirac_matched", "")) == "1":
            density = float(row["mv_units_per_1000_selected_site_cases_2050"])
            facecolor = cmap(norm(max(density, norm.vmin)))
            collection = "value"
        elif row is not None:
            facecolor = PALETTE["missing"]
            collection = "missing"
        else:
            facecolor = PALETTE["not_included"]
            collection = "not_included"

        largest_area = -1.0
        largest_centroid: tuple[float, float] | None = None
        for exterior in polygon_exteriors(feature.get("geometry", {})):
            if len(exterior) < 3:
                continue
            patch = Polygon(exterior, closed=True)
            if collection == "value":
                value_patches.append(patch)
                value_facecolors.append(facecolor)
            elif collection == "missing":
                missing_patches.append(patch)
            else:
                not_included_patches.append(patch)
            if row is not None and str(row.get("acceleration_risk", "")) == "1":
                highlight_patches.append(Polygon(exterior, closed=True))
            area, centroid = polygon_area_and_centroid(exterior)
            if area > largest_area:
                largest_area = area
                largest_centroid = centroid

        if iso3 is not None and largest_centroid is not None:
            centroids[iso3] = largest_centroid

    if not_included_patches:
        ax.add_collection(
            PatchCollection(
                not_included_patches,
                facecolor=PALETTE["not_included"],
                edgecolor="white",
                linewidth=0.25,
                zorder=1,
            )
        )
    if value_patches:
        ax.add_collection(
            PatchCollection(
                value_patches,
                facecolor=value_facecolors,
                edgecolor="white",
                linewidth=0.30,
                zorder=2,
            )
        )
    if missing_patches:
        ax.add_collection(
            PatchCollection(
                missing_patches,
                facecolor=PALETTE["missing"],
                edgecolor="white",
                linewidth=0.25,
                hatch="////",
                zorder=3,
            )
        )
    if highlight_patches:
        ax.add_collection(
            PatchCollection(
                highlight_patches,
                facecolor=(1, 1, 1, 0),
                edgecolor=PALETTE["signal_dark"],
                linewidth=0.82,
                zorder=4,
            )
        )

    ax.set_xlim(-180, 180)
    ax.set_ylim(-58, 85)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")
    return centroids


def annotate_top_map_countries(ax: plt.Axes, df: pd.DataFrame, centroids: dict[str, tuple[float, float]]) -> None:
    top = df.dropna(subset=["unit_pressure_rank"]).sort_values("unit_pressure_rank").head(10)
    for _, row in top.iterrows():
        iso3 = row["country_iso3"]
        if iso3 not in centroids:
            continue
        x, y = centroids[iso3]
        dx, dy = LABEL_OFFSETS.get(iso3, (18, 6))
        ax.annotate(
            f"{int(row['unit_pressure_rank'])} {iso3}",
            xy=(x, y),
            xytext=(x + dx, y + dy),
            ha="left" if dx > 0 else "right",
            va="center",
            fontsize=5.7,
            color=PALETTE["black"],
            arrowprops={
                "arrowstyle": "-",
                "lw": 0.45,
                "color": PALETTE["neutral_mid"],
                "shrinkA": 1.5,
                "shrinkB": 1.5,
            },
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 0.12},
            zorder=8,
        )


def draw_panel_b(ax: plt.Axes, df: pd.DataFrame, resource_q25: float) -> None:
    top = df.dropna(subset=["unit_pressure_rank"]).sort_values("unit_pressure_rank").head(15).copy()
    top = top.sort_values("unit_pressure_rank", ascending=False)
    y = np.arange(len(top))
    density = top["mv_units_per_1000_selected_site_cases_2050"].to_numpy()
    cases = top["selected_site_cases_2050"].to_numpy()
    sizes = 20 + 150 * np.sqrt(cases / np.nanmax(cases))
    colors = np.where(top["acceleration_risk"].astype(str).eq("1"), PALETTE["signal"], PALETTE["neutral_mid"])

    ax.hlines(y, 0.008, density, color="#BFBFBF", lw=0.65, zorder=1)
    ax.scatter(density, y, s=sizes, c=colors, edgecolor="white", linewidth=0.55, zorder=3)
    ax.axvline(resource_q25, color=PALETTE["signal_dark"], lw=0.75, ls="--")
    ax.text(
        resource_q25 * 1.05,
        len(top) - 0.15,
        "country q25\ndensity",
        fontsize=5.2,
        color=PALETTE["signal_dark"],
        ha="left",
        va="top",
        rotation=90,
    )

    labels = [
        wrap_country_label(int(row["unit_pressure_rank"]), row["country_iso3"], row["gco_country"])
        for _, row in top.iterrows()
    ]
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xscale("log")
    ax.set_xlim(0.008, 0.75)
    ax.set_xticks([0.01, 0.03, 0.1, 0.3, 0.75])
    ax.set_xticklabels(["0.01", "0.03", "0.1", "0.3", "0.75"])
    ax.set_xlabel("Latest-reported MV units per 1000 selected-site cases, 2050", labelpad=5.0)
    ax.tick_params(axis="y", length=0, pad=2)

    for yi, (_, row) in zip(y, top.iterrows()):
        growth = 100 * float(row["selected_site_relative_case_growth_2050"])
        mv = int(float(row["mv_therapy_units"]))
        text = f"+{growth:.0f}%; {mv} MV"
        ax.text(0.76, yi, text, ha="left", va="center", fontsize=5.15, color=PALETTE["neutral_dark"])

    ax.text(
        0.01,
        -0.31,
        "Red = high growth + lower-quartile unit density;\ngrey = other high-pressure country. Point area encodes 2050 selected-site cases.",
        transform=ax.transAxes,
        fontsize=5.25,
        color=PALETTE["neutral_mid"],
        ha="left",
        va="top",
        linespacing=1.15,
    )


def make_region_summary(df: pd.DataFrame) -> pd.DataFrame:
    matched = df[df["dirac_matched"].astype(str).eq("1")].copy()
    matched = matched[matched["who_region"].notna() & (matched["who_region"].astype(str) != "")]
    rows = []
    for region, sub in matched.groupby("who_region", sort=False):
        cases_2024 = float(sub["selected_site_cases_2024"].sum())
        cases_2050 = float(sub["selected_site_cases_2050"].sum())
        mv_units = float(sub["mv_therapy_units"].sum())
        rows.append(
            {
                "who_region": region,
                "matched_countries": int(len(sub)),
                "acceleration_risk_countries": int((sub["acceleration_risk"].astype(str) == "1").sum()),
                "acceleration_risk_share": int((sub["acceleration_risk"].astype(str) == "1").sum()) / len(sub) if len(sub) else math.nan,
                "selected_site_cases_2024_matched": cases_2024,
                "selected_site_cases_2050_matched": cases_2050,
                "selected_site_case_increase_2050_matched": cases_2050 - cases_2024,
                "selected_site_relative_case_growth_2050_matched": (cases_2050 - cases_2024) / cases_2024 if cases_2024 else math.nan,
                "mv_therapy_units_matched": mv_units,
                "mv_units_per_1000_selected_site_cases_2050_matched": mv_units * 1000 / cases_2050 if cases_2050 else math.nan,
            }
        )
    return pd.DataFrame(rows)


def draw_panel_c(ax: plt.Axes, region_df: pd.DataFrame, growth_q75: float, resource_q25: float) -> None:
    plot_df = region_df[region_df["who_region"].isin(REGION_COLORS)].copy()
    plot_df = plot_df.sort_values(["acceleration_risk_share", "acceleration_risk_countries"], ascending=True)
    y = np.arange(len(plot_df))
    x = 100 * plot_df["acceleration_risk_share"].to_numpy()
    matched = plot_df["matched_countries"].to_numpy()
    sizes = 55 + 255 * np.sqrt(matched / np.nanmax(matched))

    ax.axvspan(0, 25, color=PALETTE["neutral_pale"], alpha=0.55, zorder=0)
    ax.axvline(25, color=PALETTE["neutral_mid"], lw=0.55, ls=":")
    ax.hlines(y, 0, x, color=PALETTE["neutral_light"], lw=1.4, zorder=1)

    for idx, (_, row) in enumerate(plot_df.iterrows()):
        region = row["who_region"]
        color = REGION_COLORS.get(region, PALETTE["neutral_mid"])
        xi = 100 * float(row["acceleration_risk_share"])
        si = sizes[idx]
        ax.scatter(xi, idx, s=si, color=color, alpha=0.92, edgecolor="white", linewidth=0.75, zorder=3)
        count_label = f"{int(row['acceleration_risk_countries'])}/{int(row['matched_countries'])}"
        ax.text(min(xi + 4.2, 95), idx, count_label, ha="left", va="center", fontsize=5.7, color=PALETTE["black"], zorder=4)

    ax.set_xlim(0, 100)
    ax.set_ylim(-0.55, len(plot_df) - 0.45)
    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["who_region"].tolist())
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xticklabels(["0", "25", "50", "75", "100"])
    ax.set_xlabel("Acceleration-risk countries among matched countries (%)")
    ax.set_ylabel("")
    ax.tick_params(axis="y", length=0)
    ax.text(
        0.98,
        0.03,
        "Labels: acceleration-risk countries /\nmatched countries; point area encodes matched countries",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=5.15,
        color=PALETTE["neutral_mid"],
    )


def make_figure() -> dict[str, Path]:
    apply_publication_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    DOC_DIR.mkdir(parents=True, exist_ok=True)

    df, log = load_data()
    metrics = log["metrics"]
    growth_q75 = float(metrics["growth_q75"])
    resource_q25 = float(metrics["resource_q25"])
    highlow_count = int(metrics["acceleration_risk_count"])
    complete_cases = int(metrics["complete_case_records"])

    cmap = LinearSegmentedColormap.from_list(
        "rt_resource_density",
        [PALETTE["signal"], PALETTE["signal_soft"], PALETTE["neutral_pale"], PALETTE["blue_mid"], PALETTE["blue_main"]],
    )
    norm = LogNorm(vmin=0.01, vmax=5.0, clip=True)

    fig = plt.figure(figsize=(7.65, 6.08))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.90, 1.05], hspace=0.43, wspace=0.42)
    fig.subplots_adjust(left=0.055, right=0.985, top=0.970, bottom=0.105)
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])

    add_panel_label(ax_a, "a", -0.035, 1.026)
    centroids = draw_world_map(ax_a, df, cmap, norm)
    annotate_top_map_countries(ax_a, df, centroids)
    ax_a.text(
        0.030,
        0.160,
        f"{complete_cases} matched countries\n{highlow_count} acceleration-risk countries",
        transform=ax_a.transAxes,
        fontsize=5.8,
        color=PALETTE["black"],
        ha="left",
        va="bottom",
        linespacing=1.15,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 1.2},
    )

    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cax = ax_a.inset_axes([0.335, -0.105, 0.330, 0.034])
    cbar = fig.colorbar(sm, cax=cax, orientation="horizontal", ticks=[0.01, 0.05, 0.1, 0.5, 1, 5])
    cbar.set_label("Latest-reported MV units per 1000 selected-site cases, 2050", fontsize=5.6, labelpad=1.0)
    cbar.ax.set_xticklabels(["0.01", "0.05", "0.1", "0.5", "1", "5"])
    cbar.ax.tick_params(labelsize=5.2, length=1.8, pad=1.0)
    cbar.outline.set_linewidth(0.45)

    legend_handles = [
        Patch(facecolor=(1, 1, 1, 0), edgecolor=PALETTE["signal_dark"], linewidth=0.9, label="High growth + lower-quartile\nunit density"),
        Patch(facecolor=PALETTE["missing"], edgecolor=PALETTE["neutral_mid"], hatch="////", linewidth=0.45, label="DIRAC missing"),
    ]
    ax_a.legend(handles=legend_handles, loc="lower left", bbox_to_anchor=(0.030, 0.020), borderpad=0.2, handlelength=1.2)

    add_panel_label(ax_b, "b", -0.16, 1.015)
    draw_panel_b(ax_b, df, resource_q25)

    region_df = make_region_summary(df)
    add_panel_label(ax_c, "c", -0.13, 1.015)
    draw_panel_c(ax_c, region_df, growth_q75, resource_q25)

    bolden_figure_text(fig)

    out_base = OUT_DIR / "figure_1_radiotherapy_resource_mismatch_3panel"
    paths = {
        "svg": out_base.with_suffix(".svg"),
        "pdf": out_base.with_suffix(".pdf"),
        "png": out_base.with_suffix(".png"),
    }
    fig.savefig(paths["svg"], bbox_inches="tight")
    fig.savefig(paths["pdf"], bbox_inches="tight")
    fig.savefig(paths["png"], dpi=600, bbox_inches="tight")
    plt.close(fig)

    country_source_cols = [
        "country_iso3",
        "gco_country",
        "who_region",
        "income_label",
        "dirac_matched",
        "unit_pressure_rank",
        "selected_site_cases_per_mv_unit_2050",
        "high_growth_q75",
        "lower_quartile_unit_density",
        "acceleration_risk",
        "selected_site_cases_2024",
        "selected_site_cases_2050",
        "selected_site_case_increase_2050",
        "selected_site_relative_case_growth_2050",
        "mv_therapy_units",
        "mv_units_per_1000_selected_site_cases_2050",
    ]
    df[country_source_cols].to_csv(
        COUNTRY_SOURCE_DATA,
        index=False,
        encoding="utf-8-sig",
    )
    region_df.to_csv(
        REGION_SOURCE_DATA,
        index=False,
        encoding="utf-8-sig",
    )
    write_design_note(paths, complete_cases, highlow_count, growth_q75, resource_q25)
    return paths


def write_design_note(
    paths: dict[str, Path],
    complete_cases: int,
    highlow_count: int,
    growth_q75: float,
    resource_q25: float,
) -> None:
    lines = [
        "# Radiotherapy Resource Mismatch Figure Contract",
        "",
        "Core conclusion: Countries with upper-quartile projected selected cancer-site incidence growth are concentrated in settings with sparse latest-reported megavoltage unit density.",
        "Figure archetype: asymmetric mixed-modality figure.",
        "Target journal/output: Research Letter main figure, double-column width, editable SVG/PDF plus high-resolution PNG.",
        "Backend: Python/matplotlib only.",
        "Final size: 7.65 x 6.08 inches before tight bounding-box export.",
        "",
        "## Panel Map",
        "",
        "- a: Hero world map; fill encodes latest-reported MV units per 1000 projected 2050 selected-site cancer cases, red outline marks countries crossing both high-growth and lower-quartile unit-density thresholds, hatch marks DIRAC-missing observations.",
        "- b: Country pressure profile plot; y-axis is selected-site cases per MV-unit pressure rank, x-axis is MV density, point area encodes projected 2050 selected-site cases, right text gives relative growth and latest-reported MV units.",
        "- c: WHO-region acceleration-risk proportion plot; x-axis is acceleration-risk countries as a share of matched countries, labels give acceleration-risk count over matched countries, and point area encodes the number of matched countries.",
        "",
        "## Evidence Hierarchy",
        "",
        "- Hero evidence: geographic co-localisation of projected selected-site incidence growth and sparse latest-reported MV-unit density.",
        "- Validation evidence: ranked country profiles show countries with the highest selected-site cases per latest-reported MV unit.",
        "- Regional synthesis: WHO-region proportions show where acceleration-risk countries are concentrated among matched DIRAC records.",
        "",
        "## Thresholds and n",
        "",
        f"- Complete-case countries: {complete_cases}.",
        f"- High-growth threshold: relative selected-site case growth q75 = {growth_q75:.4f}.",
        f"- Lower-quartile unit-density threshold: MV units per 1000 projected 2050 selected-site cases q25 = {resource_q25:.4f}.",
        f"- Acceleration-risk countries: {highlow_count}.",
        "",
        "## Reviewer-Risk Notes",
        "",
        "- DIRAC-absent countries are treated as missing resource observations without assigned measured resource density.",
        "- The figure visualises selected cancer-site incidence, not modelled radiotherapy demand or utilisation.",
        "- MV units are latest-reported DIRAC country-table counts and are not projected to 2050.",
        "- Panel c reports regional proportions of country-level classifications, avoiding classification of regional aggregate ratios.",
        "",
        "## Exported Files",
        "",
    ]
    for key, path in paths.items():
        lines.append(f"- {key}: `{path.relative_to(PROJECT_ROOT)}`")
    lines.extend(
        [
            f"- country source data: `{COUNTRY_SOURCE_DATA.relative_to(PROJECT_ROOT)}`",
            f"- region source data: `{REGION_SOURCE_DATA.relative_to(PROJECT_ROOT)}`",
        ]
    )
    (DOC_DIR / "figure_contract.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    paths = make_figure()
    print(json.dumps({key: str(path.relative_to(PROJECT_ROOT)) for key, path in paths.items()}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
