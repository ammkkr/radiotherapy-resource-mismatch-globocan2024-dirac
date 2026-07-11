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
DATASET = PROJECT_ROOT / "data" / "radiotherapy_resource_mismatch_country_v6_sensitivity.csv"
BUILD_LOG = PROJECT_ROOT / "data" / "build_log_sanitized.json"
ADDITIONAL_LOG = PROJECT_ROOT / "data" / "additional_analyses_log.json"
MISSING_BOUNDS = PROJECT_ROOT / "data" / "table_21_missing_resource_identification_bounds.csv"
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
    "update_new": "#D18B35",
    "update_exit": "#315F8C",
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
        "gco2022_selected_growth_2022_2050",
        "gco2022_growth_percentile",
        "gco2024_growth_percentile_common",
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
        primary_log = json.load(handle)
    with ADDITIONAL_LOG.open("r", encoding="utf-8") as handle:
        additional_log = json.load(handle)
    return df, {"primary": primary_log, "additional": additional_log}


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
            if row is not None and str(row.get("current_screen_positive", "")).lower() in {"1", "1.0", "true"}:
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


def draw_panel_b(ax: plt.Axes, df: pd.DataFrame) -> None:
    plot_df = df[df["burden_version_transition"].ne("Not comparable")].copy()
    cases = plot_df["selected_site_cases_2050"].to_numpy()
    log_cases = np.log10(np.clip(cases, 1, None))
    sizes = 9 + 54 * (log_cases - log_cases.min()) / (log_cases.max() - log_cases.min())
    plot_df["point_size"] = sizes

    styles = {
        "Not screen-positive in either version": (PALETTE["neutral_light"], 0.52, 0),
        "Retained": (PALETTE["signal"], 0.82, 1),
        "No longer screen-positive": (PALETTE["update_exit"], 0.96, 2),
        "New in GCO version 2024": (PALETTE["update_new"], 1.00, 3),
    }
    for status, (color, alpha, zorder) in styles.items():
        sub = plot_df[plot_df["burden_version_transition"].eq(status)]
        ax.scatter(
            sub["gco2022_growth_percentile"],
            sub["gco2024_growth_percentile_common"],
            s=sub["point_size"],
            color=color,
            alpha=alpha,
            edgecolor="white",
            linewidth=0.45,
            zorder=2 + zorder,
        )

    ax.plot([0, 100], [0, 100], color=PALETTE["neutral_mid"], lw=0.55, ls=":", zorder=0)
    ax.axvline(75, color=PALETTE["signal_dark"], lw=0.65, ls="--", zorder=1)
    ax.axhline(75, color=PALETTE["signal_dark"], lw=0.65, ls="--", zorder=1)
    ax.text(76.5, 3, "q75", fontsize=5.1, color=PALETTE["signal_dark"], ha="left", va="bottom")
    ax.text(3, 76.5, "q75", fontsize=5.1, color=PALETTE["signal_dark"], ha="left", va="bottom")

    offsets = {"GTM": (8, 8), "NGA": (8, -11), "TGO": (8, 8), "ZWE": (-29, -12)}
    for iso3, (dx, dy) in offsets.items():
        row = plot_df[plot_df["country_iso3"].eq(iso3)].iloc[0]
        ax.annotate(
            iso3,
            xy=(row["gco2022_growth_percentile"], row["gco2024_growth_percentile_common"]),
            xytext=(dx, dy),
            textcoords="offset points",
            fontsize=5.4,
            color=PALETTE["update_new"] if iso3 == "GTM" else PALETTE["update_exit"],
            ha="left",
            va="center",
            arrowprops={"arrowstyle": "-", "lw": 0.45, "color": PALETTE["neutral_mid"]},
        )

    ax.set_xlim(0, 101)
    ax.set_ylim(0, 101)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xlabel("Projected growth percentile, GCO version 2022")
    ax.set_ylabel("Projected growth percentile, GCO version 2024")
    handles = [
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=PALETTE["signal"], markeredgecolor="white", markersize=4.8, label="Retained (22)"),
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=PALETTE["update_new"], markeredgecolor="white", markersize=4.8, label="New (1)"),
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=PALETTE["update_exit"], markeredgecolor="white", markersize=4.8, label="No longer (3)"),
    ]
    ax.legend(handles=handles, loc="lower right", borderpad=0.1, handletextpad=0.3, labelspacing=0.25)
    ax.text(
        0.01,
        -0.27,
        "Same latest-reported DIRAC snapshot in both versions.\nPoint area encodes GCO version 2024 selected-site cases in 2050.",
        transform=ax.transAxes,
        fontsize=5.1,
        color=PALETTE["neutral_mid"],
        ha="left",
        va="top",
        linespacing=1.15,
    )


def load_region_bounds() -> pd.DataFrame:
    bounds = pd.read_csv(MISSING_BOUNDS)
    return bounds[
        bounds["Stratum"].eq("WHO region") & bounds["Category"].isin(REGION_COLORS)
    ].copy()


def draw_panel_c(ax: plt.Axes, region_df: pd.DataFrame) -> None:
    plot_df = region_df.sort_values(
        ["Lower-bound proportion, %", "Upper-bound proportion, %"], ascending=True
    ).copy()
    y = np.arange(len(plot_df))
    lower = plot_df["Lower-bound proportion, %"].to_numpy()
    upper = plot_df["Upper-bound proportion, %"].to_numpy()

    ax.hlines(y, lower, upper, color=PALETTE["neutral_light"], lw=2.2, zorder=1)
    ax.scatter(lower, y, s=34, color=PALETTE["signal"], edgecolor="white", linewidth=0.55, zorder=3)
    ax.scatter(upper, y, s=38, facecolor="white", edgecolor=PALETTE["signal_dark"], linewidth=0.85, zorder=4)

    for yi, (_, row) in enumerate(plot_df.iterrows()):
        observed = int(row["Observed screen-positive, No."])
        possible = observed + int(row["High-growth resource-unknown, No."])
        total = int(row["All GCO country records, No."])
        ax.text(
            min(float(row["Upper-bound proportion, %"]) + 2.0, 66.5),
            yi,
            f"{observed}/{total} to {possible}/{total}",
            ha="left",
            va="center",
            fontsize=5.45,
            color=PALETTE["black"],
        )

    ax.set_xlim(0, 68)
    ax.set_ylim(-0.55, len(plot_df) - 0.45)
    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["Category"].tolist())
    ax.set_xticks([0, 20, 40, 60])
    ax.set_xlabel("Countries meeting both thresholds, bounds (%)")
    ax.tick_params(axis="y", length=0)
    handles = [
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=PALETTE["signal"], markeredgecolor="white", markersize=4.8, label="Observed lower bound"),
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor="white", markeredgecolor=PALETTE["signal_dark"], markersize=4.8, label="High-growth missing upper bound"),
    ]
    ax.legend(handles=handles, loc="lower right", borderpad=0.1, handletextpad=0.3, labelspacing=0.25)
    ax.text(
        0.01,
        -0.27,
        "Bounds use all GCO country records. Upper bounds add high-growth\nDIRAC-missing countries; intervals are not confidence intervals.",
        transform=ax.transAxes,
        fontsize=5.1,
        color=PALETTE["neutral_mid"],
        ha="left",
        va="top",
        linespacing=1.15,
    )


def make_figure() -> dict[str, Path]:
    apply_publication_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    DOC_DIR.mkdir(parents=True, exist_ok=True)

    df, log = load_data()
    metrics = log["primary"]["metrics"]
    version_metrics = log["additional"]["version"]
    growth_q75 = float(metrics["growth_q75"])
    resource_q25 = float(metrics["resource_q25"])
    highlow_count = int(metrics["meets_both_thresholds_count"])
    complete_cases = int(metrics["complete_case_records"])
    high_growth_unknown = int((~df["dirac_matched"].astype(str).eq("1") & df["high_growth_q75"].astype(str).eq("1")).sum())

    cmap = LinearSegmentedColormap.from_list(
        "rt_resource_density",
        [PALETTE["signal"], PALETTE["signal_soft"], PALETTE["neutral_pale"], PALETTE["blue_mid"], PALETTE["blue_main"]],
    )
    norm = LogNorm(vmin=0.01, vmax=5.0, clip=True)

    fig = plt.figure(figsize=(7.65, 6.18))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.90, 1.05], hspace=0.43, wspace=0.42)
    fig.subplots_adjust(left=0.060, right=0.985, top=0.970, bottom=0.120)
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])

    add_panel_label(ax_a, "a", -0.035, 1.026)
    centroids = draw_world_map(ax_a, df, cmap, norm)
    ax_a.text(
        0.030,
        0.160,
        f"{complete_cases} matched countries\n{highlow_count} met both thresholds\n{high_growth_unknown} high-growth, resources unknown",
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
        Patch(facecolor=(1, 1, 1, 0), edgecolor=PALETTE["signal_dark"], linewidth=0.9, label="Both thresholds"),
        Patch(facecolor=PALETTE["missing"], edgecolor=PALETTE["neutral_mid"], hatch="////", linewidth=0.45, label="DIRAC missing"),
    ]
    ax_a.legend(handles=legend_handles, loc="lower left", bbox_to_anchor=(0.030, 0.020), borderpad=0.2, handlelength=1.2)

    add_panel_label(ax_b, "b", -0.16, 1.015)
    draw_panel_b(ax_b, df)

    region_df = load_region_bounds()
    add_panel_label(ax_c, "c", -0.13, 1.015)
    draw_panel_c(ax_c, region_df)

    bolden_figure_text(fig)

    out_base = OUT_DIR / "figure_1_radiotherapy_resource_mismatch_3panel"
    paths = {
        "svg": out_base.with_suffix(".svg"),
        "pdf": out_base.with_suffix(".pdf"),
        "png": out_base.with_suffix(".png"),
    }
    fig.savefig(paths["svg"], bbox_inches="tight")
    svg_text = paths["svg"].read_text(encoding="utf-8")
    paths["svg"].write_text(
        "\n".join(line.rstrip() for line in svg_text.splitlines()) + "\n",
        encoding="utf-8",
    )
    fig.savefig(paths["pdf"], bbox_inches="tight")
    fig.savefig(paths["png"], dpi=600, bbox_inches="tight")
    plt.close(fig)

    country_source_cols = [
        "country_iso3",
        "gco_country",
        "country_display",
        "who_region",
        "income_label",
        "dirac_matched",
        "unit_pressure_rank",
        "selected_site_cases_per_mv_unit_2050",
        "high_growth_q75",
        "lower_quartile_unit_density",
        "meets_both_thresholds",
        "selected_site_cases_2024",
        "selected_site_cases_2050",
        "selected_site_case_increase_2050",
        "selected_site_relative_case_growth_2050",
        "mv_therapy_units",
        "mv_units_per_1000_selected_site_cases_2050",
        "gco2022_selected_growth_2022_2050",
        "gco2022_growth_percentile",
        "gco2024_growth_percentile_common",
        "gco2022_screen_positive",
        "current_screen_positive",
        "burden_version_transition",
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
    write_design_note(
        paths,
        complete_cases,
        highlow_count,
        high_growth_unknown,
        growth_q75,
        resource_q25,
        version_metrics,
    )
    return paths


def write_design_note(
    paths: dict[str, Path],
    complete_cases: int,
    highlow_count: int,
    high_growth_unknown: int,
    growth_q75: float,
    resource_q25: float,
    version_metrics: dict[str, Any],
) -> None:
    lines = [
        "# Radiotherapy Resource Mismatch Figure Contract",
        "",
        "Core conclusion: The GCO version 2024 burden update retained a stable country core while changing four classifications, and missing resource observations widened regional uncertainty.",
        "Figure archetype: asymmetric mixed-modality figure.",
        "Target journal/output: Research Letter main figure, double-column width, editable SVG/PDF plus high-resolution PNG.",
        "Backend: Python/matplotlib only.",
        "Final size: 7.65 x 6.18 inches before tight bounding-box export.",
        "",
        "## Panel Map",
        "",
        "- a: Hero world map; fill encodes latest-reported MV units per 1000 projected 2050 selected-site cancer cases, red outline marks countries meeting both thresholds, and hatching marks DIRAC-missing observations.",
        "- b: Burden-version comparison; axes are within-version projected growth percentiles for 149 common matched countries, colors identify retained, new, and no-longer-screen-positive countries, and point area encodes version 2024 selected-site cases in 2050.",
        "- c: WHO-region partial-identification bounds; the lower endpoint is the observed count divided by all regional GCO country records, and the upper endpoint additionally treats high-growth DIRAC-missing countries as meeting the resource criterion.",
        "",
        "## Evidence Hierarchy",
        "",
        "- Hero evidence: geographic co-localisation of projected selected-site incidence growth and sparse latest-reported MV-unit density.",
        "- Update evidence: the burden-version comparison separates the GCO data revision from resource change by holding the same latest-reported DIRAC snapshot fixed.",
        "- Missingness evidence: regional intervals show the range compatible with observed and high-growth resource-unknown country records.",
        "",
        "## Thresholds and n",
        "",
        f"- Complete-case countries: {complete_cases}.",
        f"- High-growth threshold: relative selected-site case growth q75 = {growth_q75:.4f}.",
        f"- Lower-quartile unit-density threshold: MV units per 1000 projected 2050 selected-site cases q25 = {resource_q25:.4f}.",
        f"- Countries meeting both primary thresholds: {highlow_count}.",
        f"- High-growth countries with unknown resource status: {high_growth_unknown}.",
        f"- Countries retained across GCO burden versions: {version_metrics['retained_count']}.",
        f"- New in version 2024: {', '.join(version_metrics['new_in_2024'])}.",
        f"- No longer screen-positive in version 2024: {', '.join(version_metrics['no_longer_in_2024'])}.",
        "",
        "## Reviewer-Risk Notes",
        "",
        "- DIRAC-absent countries are treated as missing resource observations without assigned measured resource density.",
        "- The figure visualises selected cancer-site incidence, not modelled radiotherapy demand or utilisation.",
        "- MV units are latest-reported DIRAC country-table counts and are not projected to 2050.",
        "- Panel b is a burden-version comparison, not a longitudinal resource analysis; the latest-reported DIRAC snapshot is fixed.",
        "- Panel c intervals are deterministic missing-resource bounds, not confidence intervals.",
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
