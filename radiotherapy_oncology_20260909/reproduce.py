"""Reproduce the reviewed analysis from shared country-level aggregates."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def reproduce(data_dir, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(data_dir / "country_data.csv").set_index("country_iso3")
    assert d.index.is_unique and len(d) == 186
    matched = d.units.notna()
    growth = d.cases_2050 / d.cases_2024 - 1
    density = 1000 * d.units / d.cases_2050

    def screen(g, r, eligible=matched, quantiles=(.75, .25), fixed=None):
        valid = eligible & g.notna() & r.notna()
        cut = fixed if fixed is not None else (g[valid].quantile(quantiles[0]), r[valid].quantile(quantiles[1]))
        return valid & g.ge(cut[0]) & r.le(cut[1]), tuple(float(v) for v in cut)

    primary, cut = screen(growth, density)
    assert primary.equals(d.primary)
    assert np.allclose(growth, d.growth)
    assert np.allclose(density[matched], d.density[matched])
    assert np.allclose(density[matched], d.baseline_density[matched] / (1 + growth[matched]))
    scenarios = []

    def record(name, flag, eligible=matched):
        scenarios.append({"Analysis": name, "Countries analyzed": int(eligible.sum()),
            "Countries meeting both": int(flag.sum()), "Eligible primary": int((primary & eligible).sum()),
            "Primary retained": int((primary & flag).sum()),
            "Added": ";".join(d.index[flag & ~primary]),
            "Not retained": ";".join(d.index[primary & eligible & ~flag]),
            "Ineligible primary": ";".join(d.index[primary & ~eligible])})

    record("Primary", primary)
    alternatives = [
        ("All cancers excluding NMSC", d.all40_2050/d.all40_2024-1, 1000*d.units/d.all40_2050, matched, (.75,.25)),
        ("Baseline case denominator", growth, 1000*d.units/d.cases_2024, matched, (.75,.25)),
        ("Stricter percentiles", growth, density, matched, (.8,.2)),
        ("2040 horizon", d.cases_2040/d.cases_2024-1, 1000*d.units/d.cases_2040, matched, (.75,.25)),
        ("Centre density", growth, 10000*d.centres/d.cases_2050, matched, (.75,.25)),
        ("Mortality growth", d.deaths_2050/d.deaths_2024-1, density, matched, (.75,.25)),
        ("Baseline at least 1000 cases", growth, density, matched & d.cases_2024.ge(1000), (.75,.25)),
    ]
    support = pd.Series(0, index=d.index)
    for name, g, r, eligible, q in alternatives:
        selected, _ = screen(g, r, eligible, q)
        support += selected.astype(int)
        record(name, selected, eligible)
    assert support.equals(d.support_7)
    thirteen, _ = screen(d.thirteen_site_cases_2050/d.thirteen_site_cases_2024-1, 1000*d.units/d.thirteen_site_cases_2050)
    record("Thirteen-site extraction without rectum", thirteen)
    for year in (2023, 2024, 2025, 2026):
        eligible = matched & d.last_update.ge(year)
        flag, _ = screen(growth, density, eligible)
        record(f"DIRAC update {year} or later", flag, eligible)
    for amount in (1, 2, 5):
        flag, _ = screen(growth, 1000*(d.units+amount)/d.cases_2050, fixed=cut)
        record(f"Add {amount} unit(s), fixed thresholds", flag)
    for factor in (.8, 1.2, 1.5):
        flag, _ = screen(growth, density*factor, fixed=cut)
        record(f"Inventory x{factor}, fixed thresholds", flag)
    computed = pd.DataFrame(scenarios).set_index("Analysis")
    supplied = pd.read_csv(data_dir/"sensitivity_results.csv").set_index("Analysis").fillna("")
    pd.testing.assert_frame_equal(computed, supplied[computed.columns], check_dtype=False)
    computed.to_csv(output_dir/"sensitivity_results.csv")

    old_g = d.release2022_cases_2050/d.release2022_cases_2025-1
    new_g = d.cases_2050/d.cases_2025-1
    old_d = 1000*d.units/d.release2022_cases_2050
    common = matched & old_g.notna()
    old, old_cut = screen(old_g, old_d, common)
    new, new_cut = screen(new_g, density, common)
    fixed, _ = screen(new_g, density, common, fixed=old_cut)
    supplied_version = pd.read_csv(data_dir/"version_comparison.csv").set_index("country_iso3")
    assert old.equals(supplied_version.old_selected) and new.equals(supplied_version.new_selected)
    version = {"n": int(common.sum()), "old_count": int(old.sum()), "new_count": int(new.sum()),
        "retained": int((old & new).sum()), "new": d.index[new & ~old].tolist(),
        "exited": d.index[old & ~new].tolist(), "jaccard": float((old & new).sum()/(old | new).sum()),
        "common_growth_rho": float(spearmanr(old_g[common], new_g[common]).statistic),
        "thresholds_old": old_cut, "thresholds_new": new_cut,
        "fixed_old_threshold_new_count": int(fixed.sum()),
        "fixed_old_threshold_retained": int((fixed & old).sum()),
        "fixed_old_threshold_additions": d.index[fixed & ~old].tolist(),
        "fixed_old_threshold_exits": d.index[old & ~fixed].tolist()}
    rows = []
    for group, column in [("Global", None), ("WHO region", "who_region"), ("Income group", "income_label"), ("HDI group", "hdi_label")]:
        groups = [("All countries", d)] if column is None else d.groupby(column, dropna=False)
        for category, s in groups:
            n = len(s)
            observed = int(primary.loc[s.index].sum())
            unknown_high = int((s.units.isna() & growth.loc[s.index].ge(cut[0])).sum())
            rows.append({"Group":group, "Category":category, "Countries":n,
                "Matched":int(s.units.notna().sum()), "Unmatched":int(s.units.isna().sum()),
                "Observed selected":observed, "Unmatched high growth":unknown_high,
                "Lower (%)":round(100*observed/n,1), "Upper (%)":round(100*(observed+unknown_high)/n,1)})
    bounds = pd.DataFrame(rows)
    pd.testing.assert_frame_equal(bounds, pd.read_csv(data_dir/"missingness_bounds.csv"), check_dtype=False)
    bounds.to_csv(output_dir/"missingness_bounds.csv", index=False)
    diag = []
    for label, x, y, e in [
        ("Selected-site and all-cancer growth (excluding NMSC)", growth, d.all40_2050/d.all40_2024-1, matched),
        ("Selected-site and population growth", growth, d.population_2050/d.population_2024-1, matched),
        ("Selected-site growth and 2050-case density", growth, density, matched),
        ("Selected-site growth and 2024-case density", growth, 1000*d.units/d.cases_2024, matched),
        ("2022 and 2024 releases, common 2025-2050 growth", old_g, new_g, common),
        ("2022 and 2024 releases, 2050-case density", old_d, density, common)]:
        diag.append({"Comparison":label, "Countries":int(e.sum()), "Spearman rho":round(float(spearmanr(x[e],y[e]).statistic),6)})
    correlations = pd.DataFrame(diag)
    pd.testing.assert_frame_equal(correlations, pd.read_csv(data_dir/"correlation_diagnostics.csv"))
    correlations.to_csv(output_dir/"correlation_diagnostics.csv", index=False)
    summary = {"primary_count":int(primary.sum()), "primary_cutoffs":cut,
        "primary_cases_2024":int(d.loc[primary,"cases_2024"].sum()),
        "primary_cases_2050":int(d.loc[primary,"cases_2050"].sum()),
        "primary_units":int(d.loc[primary,"units"].sum()),
        "primary_support_all7":int(support[primary].eq(7).sum()),
        "primary_support_ge6":int(support[primary].ge(6).sum()), "version":version}
    (output_dir/"summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print("All supplied classifications, sensitivity rows, bounds and correlations reproduced.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    reproduce(args.data, args.output)
