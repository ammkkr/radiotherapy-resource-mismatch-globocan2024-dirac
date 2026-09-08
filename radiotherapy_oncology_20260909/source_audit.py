"""Recompute the Short Communication from archived GCO and DIRAC records."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import os
ROOT = Path(os.environ["RADIOTHERAPY_SOURCE_ARCHIVE"]).resolve()
OUT = ROOT / "results" / "radiotherapy_oncology_review_v2"
SITES = [1, 3, 4, 5, 6, 9, 14, 15, 20, 23, 24, 27, 30, 31]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    checks = []
    def check(name, condition, detail=""):
        checks.append({"check": name, "passed": bool(condition), "detail": detail})
        if not condition:
            raise AssertionError(f"{name}: {detail}")

    inputs = []
    def load_version(version):
        path = ROOT / f"data/interim/gco_cancer_tomorrow_{version}_predictions_long.csv"
        inputs.append(path)
        frame = pd.read_csv(path)
        raw_paths = sorted((ROOT / f"data/raw/gco_cancer_tomorrow_{version}").glob("prediction_type*_sex0_*.json"))
        raw_rows = []
        for file in raw_paths:
            inputs.append(file)
            raw_rows.extend(json.loads(file.read_text(encoding="utf-8-sig"))["dataset"])
        raw = pd.DataFrame(raw_rows)
        keys = ["type", "sex", "id", "cancer", "year"]
        check(f"{version}_raw_unique", not raw.duplicated(keys).any())
        check(f"{version}_zero_APC", raw.APC.eq(0).all())
        long = frame[frame.sex_label.eq("both")].copy()
        check(f"{version}_long_unique", not long.duplicated(["type", "sex", "country_code", "cancer_code", "year"]).any())
        compare = long.merge(raw[keys + ["cases_pred"]], left_on=["type", "sex", "country_code", "cancer_code", "year"], right_on=keys, validate="one_to_one", suffixes=("_long", "_raw"))
        check(f"{version}_raw_long_row_count", len(compare) == len(long) == len(raw), str(len(compare)))
        check(f"{version}_raw_long_values", np.array_equal(compare.predicted_count, compare.cases_pred_raw))
        country = long[long.country_iso3.notna()].copy()
        # The archived main-site extraction omitted rectum (code 9) in both releases.
        extra_paths = sorted((ROOT / "data/raw/radiotherapy_review_rectum_20260909").glob(f"v{version}_type*_chunk*.json"))
        extra_rows = []
        for file in extra_paths:
            inputs.append(file)
            extra_rows.extend(json.loads(file.read_text(encoding="utf-8"))["dataset"])
        extra = pd.DataFrame(extra_rows)
        check(f"{version}_supplemental_unique", not extra.duplicated(keys).any())
        check(f"{version}_supplemental_parameters", extra.APC.eq(0).all() and extra.sex.eq(0).all())
        anchor = extra[extra.cancer.eq(15)].merge(raw[raw.cancer.eq(15)],on=keys,validate="one_to_one",suffixes=("_new","_old"))
        check(f"{version}_lung_anchor_complete", len(anchor) == len(extra[extra.cancer.eq(15)]) > 0)
        check(f"{version}_lung_anchor_unchanged", anchor.cases_pred_new.eq(anchor.cases_pred_old).all())
        country_meta = country[["country_code","country_iso3","population_label","who_region","income_label","hdi_label"]].drop_duplicates()
        rectum = extra[extra.cancer.eq(9)].rename(columns={"id":"country_code","cancer":"cancer_code","cases_pred":"predicted_count"})
        rectum = rectum.merge(country_meta,on="country_code",validate="many_to_one")
        rectum["measure"] = rectum.type.map({0:"incidence",1:"mortality"})
        rectum["sex_label"] = "both"
        check(f"{version}_rectum_complete",len(rectum)==len(country_meta)*14)
        country = pd.concat([country,rectum],ignore_index=True)
        selected = country[country.cancer_code.isin(SITES)]
        counts = selected.groupby(["country_iso3", "measure", "year"]).cancer_code.nunique()
        check(f"{version}_14_sites_every_cell", counts.eq(14).all())
        check(f"{version}_nonnegative_counts", selected.predicted_count.ge(0).all())
        meta = country[["country_iso3", "population_label", "who_region", "income_label", "hdi_label"]].drop_duplicates().set_index("country_iso3")
        check(f"{version}_unique_country_metadata", meta.index.is_unique)
        burden = selected.groupby(["country_iso3", "measure", "year"]).predicted_count.sum().unstack(["measure", "year"])
        all40 = country[country.cancer_code.eq(40)].pivot(index="country_iso3", columns=["measure", "year"], values="predicted_count")
        all39 = country[country.cancer_code.eq(39)].pivot(index="country_iso3", columns=["measure", "year"], values="predicted_count")
        pop = country[country.cancer_code.eq(40) & country.measure.eq("incidence")].pivot(index="country_iso3", columns="year", values="pop")
        return country, meta, burden, all40, all39, pop

    g24, meta, b24, all40, all39, pop = load_version(2024)
    _, _, b22, _, _, _ = load_version(2022)
    dirac_path = ROOT / "data/raw/dirac/dirac_countries_and_regions_raw.json"
    inputs.append(dirac_path)
    resources = pd.DataFrame(json.loads(dirac_path.read_text(encoding="utf-8-sig"))["data"]).set_index("CountryISO3")
    check("dirac_unique_iso", resources.index.is_unique)
    d = meta.copy()
    d["source_population_label"] = d.population_label
    for iso, label in {"CIV":"Cote d'Ivoire", "CUW":"Curacao", "REU":"France, La Reunion", "TUR":"Turkiye"}.items():
        d.loc[iso,"population_label"] = label
    d["units"] = resources.HePhotonAndElectronBeamRt
    d["centres"] = resources.RTCentersWithRT
    d["last_update"] = resources.LastUpdate
    for year in [2024, 2025, 2030, 2040, 2050]:
        d[f"cases_{year}"] = b24[("incidence", year)]
        d[f"deaths_{year}"] = b24[("mortality", year)]
        d[f"all40_{year}"] = all40[("incidence", year)]
        d[f"all39_{year}"] = all39[("incidence", year)]
        d[f"population_{year}"] = pop[year]
    for year in [2025, 2050]:
        d[f"release2022_cases_{year}"] = b22[("incidence", year)]
    no_rectum = g24[g24.cancer_code.isin(set(SITES)-{9}) & g24.measure.eq("incidence")].groupby(["country_iso3", "year"]).predicted_count.sum().unstack()
    for year in [2024, 2050]:
        d[f"thirteen_site_cases_{year}"] = no_rectum[year]
    d["growth"] = d.cases_2050 / d.cases_2024 - 1
    d["density"] = 1000 * d.units / d.cases_2050
    d["baseline_density"] = 1000 * d.units / d.cases_2024
    matched = d.units.notna()
    def screen(growth, density, eligible=matched, gq=.75, dq=.25, fixed=None):
        valid = eligible & growth.notna() & density.notna()
        gt, dt = fixed if fixed is not None else (float(growth[valid].quantile(gq)), float(density[valid].quantile(dq)))
        flag = valid & growth.ge(gt) & density.le(dt)
        return flag, gt, dt
    primary, gt, dt = screen(d.growth, d.density)
    d["primary"] = primary
    check("country_count", len(d) == 186)
    check("matched_count", matched.sum() == 150)
    check("observed_zero_units", (d.units.eq(0)).sum() == 0)

    old_path = ROOT / "data/processed/radiotherapy_resource_mismatch_country_v6_sensitivity.csv"
    inputs.append(old_path)
    old = pd.read_csv(old_path).rename(columns={"selected_site_cases_2050":"rt_relevant_cases_2050"}).set_index("country_iso3")
    check("rectum_addition_increases_cases", d.loc[old.index, "cases_2050"].ge(old.rt_relevant_cases_2050).all())
    check("resource_counts_match_old", np.allclose(d.loc[old.index, "units"], old.mv_therapy_units, equal_nan=True))
    legacy_all = {"old_all_2024_sum": float(old.all_cancer_cases_2024.sum()), "correct_all_excluding_nmsc_2024_sum": int(d.all40_2024.sum()), "correct_all_including_nmsc_2024_sum": int(d.all39_2024.sum())}
    check("legacy_double_count_identified", np.array_equal(old.all_cancer_cases_2024, (d.all40_2024 + d.all39_2024).loc[old.index]))

    scenarios = {}
    for name, growth, density, eligible, gq, dq in [
        ("All cancers excluding NMSC", d.all40_2050/d.all40_2024-1, 1000*d.units/d.all40_2050, matched, .75,.25),
        ("Baseline case denominator", d.growth, d.baseline_density, matched, .75,.25),
        ("Stricter percentiles", d.growth, d.density, matched, .8,.2),
        ("2040 horizon", d.cases_2040/d.cases_2024-1, 1000*d.units/d.cases_2040, matched, .75,.25),
        ("Centre density", d.growth, 10000*d.centres/d.cases_2050, matched, .75,.25),
        ("Mortality growth", d.deaths_2050/d.deaths_2024-1, d.density, matched, .75,.25),
        ("Baseline at least 1000 cases", d.growth, d.density, matched & d.cases_2024.ge(1000), .75,.25),
    ]:
        flag, a, z = screen(growth,density,eligible,gq,dq)
        scenarios[name] = flag
    d["support_7"] = pd.DataFrame(scenarios).sum(axis=1)

    def details(name, flag, eligible=matched, growth_threshold=None, density_threshold=None):
        return {"Analysis":name,"Countries analyzed":int(eligible.sum()),"Countries meeting both":int(flag.sum()),"Eligible primary":int((primary&eligible).sum()),"Primary retained":int((primary&flag).sum()),"Added": ";".join(d.index[flag & ~primary]),"Not retained": ";".join(d.index[primary & eligible & ~flag]),"Ineligible primary": ";".join(d.index[primary & ~eligible]),"Growth threshold":growth_threshold,"Density threshold":density_threshold}
    rows = [details("Primary",primary,growth_threshold=gt,density_threshold=dt)]
    for name, f in scenarios.items():
        rows.append(details(name,f,matched & d.cases_2024.ge(1000) if name.startswith("Baseline at least") else matched))
    thirteen, _, _ = screen(d.thirteen_site_cases_2050/d.thirteen_site_cases_2024-1, 1000*d.units/d.thirteen_site_cases_2050)
    check("thirteen_site_matches_legacy_classification", thirteen.eq(old.current_screen_positive.reindex(d.index)).all())
    rows.append(details("Thirteen-site extraction without rectum", thirteen))

    # Independent baseline alignment and fixed-threshold comparisons.
    old_growth = (b22[("incidence",2050)]/b22[("incidence",2025)]-1).reindex(d.index)
    old_density = (1000*d.units/b22[("incidence",2050)]).reindex(d.index)
    new_growth = d.cases_2050/d.cases_2025-1
    common = matched & old_growth.notna()
    old_flag, old_gt, old_dt = screen(old_growth,old_density,common)
    new_flag,new_gt,new_dt = screen(new_growth,d.density,common)
    fixed_flag,_,_ = screen(new_growth,d.density,common,fixed=(old_gt,old_dt))
    version = {"n":int(common.sum()),"old_count":int(old_flag.sum()),"new_count":int(new_flag.sum()),"retained":int((old_flag&new_flag).sum()),"new":d.index[new_flag&~old_flag].tolist(),"exited":d.index[old_flag&~new_flag].tolist(),"jaccard":float((old_flag&new_flag).sum()/(old_flag|new_flag).sum()),"common_growth_rho":float(spearmanr(old_growth[common],new_growth[common]).statistic),"thresholds_old":[old_gt,old_dt],"thresholds_new":[new_gt,new_dt],"fixed_old_threshold_new_count":int(fixed_flag.sum()),"fixed_old_threshold_retained":int((fixed_flag&old_flag).sum()),"fixed_old_threshold_additions":d.index[fixed_flag&~old_flag].tolist(),"fixed_old_threshold_exits":d.index[old_flag&~fixed_flag].tolist()}
    version_frame = pd.DataFrame({"old_growth":old_growth,"new_growth":new_growth,"old_density":old_density,"new_density":d.density,"old_selected":old_flag,"new_selected":new_flag,"common":common})
    version_frame["old_percentile"] = old_growth[common].rank(method="average",pct=True)*100
    version_frame["new_percentile"] = new_growth[common].rank(method="average",pct=True)*100
    version_frame.index.name="country_iso3"
    version_frame.to_csv(OUT/"version_comparison.csv")
    legacy_primary = old.current_screen_positive.reindex(d.index)
    legacy_all["primary_added_after_rectum"] = d.index[primary & ~legacy_primary].tolist()
    legacy_all["primary_exited_after_rectum"] = d.index[~primary & legacy_primary].tolist()
    legacy_all["13_site_primary_count"] = int(legacy_primary.sum())

    for cutoff in [2023,2024,2025,2026]:
        eligible = matched & d.last_update.ge(cutoff)
        flag,a,z = screen(d.growth,d.density,eligible)
        rows.append(details(f"DIRAC update {cutoff} or later",flag,eligible,a,z))
    for added in [1,2,5]:
        flag,_,_ = screen(d.growth,1000*(d.units+added)/d.cases_2050,fixed=(gt,dt))
        rows.append(details(f"Add {added} unit(s), fixed thresholds",flag))
    for factor in [0.8,1.2,1.5]:
        flag,_,_ = screen(d.growth,d.density*factor,fixed=(gt,dt))
        rows.append(details(f"Inventory x{factor}, fixed thresholds",flag))
    d["additional_units_to_cross_cutoff"] = np.maximum(0,np.floor(dt*d.cases_2050/1000-d.units)+1)
    sensitivity = pd.DataFrame(rows)
    sensitivity.to_csv(OUT/"sensitivity_results.csv",index=False)

    bounds = []
    for group, column in [("Global",None),("WHO region","who_region"),("Income group","income_label"),("HDI group","hdi_label")]:
        groups=[("All countries",d)] if column is None else d.groupby(column,dropna=False)
        for category,s in groups:
            n=len(s); o=int(s.primary.sum()); h=int((s.units.isna() & s.growth.ge(gt)).sum())
            bounds.append({"Group":group,"Category":category,"Countries":n,"Matched":int(s.units.notna().sum()),"Unmatched":int(s.units.isna().sum()),"Observed selected":o,"Unmatched high growth":h,"Lower (%)":round(100*o/n,1),"Upper (%)":round(100*(o+h)/n,1)})
    pd.DataFrame(bounds).to_csv(OUT/"missingness_bounds.csv",index=False)

    diag=[]
    for label,x,y,e in [
        ("Selected-site and all-cancer growth (excluding NMSC)",d.growth,d.all40_2050/d.all40_2024-1,matched),
        ("Selected-site and population growth",d.growth,d.population_2050/d.population_2024-1,matched),
        ("Selected-site growth and 2050-case density",d.growth,d.density,matched),
        ("Selected-site growth and 2024-case density",d.growth,d.baseline_density,matched),
        ("2022 and 2024 releases, common 2025-2050 growth",old_growth,new_growth,common),
        ("2022 and 2024 releases, 2050-case density",old_density,d.density,common),
    ]:
        diag.append({"Comparison":label,"Countries":int(e.sum()),"Spearman rho":round(float(spearmanr(x[e],y[e]).statistic),6)})
    pd.DataFrame(diag).to_csv(OUT/"correlation_diagnostics.csv",index=False)

    check("density_growth_identity",np.allclose(d.loc[matched,"density"],d.loc[matched,"baseline_density"]/(1+d.loc[matched,"growth"])))
    p=d[primary]
    statistics={"country_count":len(d),"matched_count":int(matched.sum()),"primary_count":len(p),"growth_cutoff":gt,"density_cutoff":dt,"selected_2024":int(d.cases_2024.sum()),"selected_2050":int(d.cases_2050.sum()),"primary_cases_2024":int(p.cases_2024.sum()),"primary_cases_2050":int(p.cases_2050.sum()),"primary_units":int(p.units.sum()),"primary_increment":int((p.cases_2050-p.cases_2024).sum()),"matched_increment":int((d.cases_2050-d.cases_2024)[matched].sum()),"primary_low_lower_income":int(p.income_label.isin(["Low income","Lower middle income"]).sum()),"support_all7":int(p.support_7.eq(7).sum()),"support_ge6":int(p.support_7.ge(6).sum()),"support_distribution":{str(k):int(v) for k,v in p.support_7.value_counts().items()},"primary_unit_multiplier_to_cutoff_median":float((dt/p.density).median()),"primary_additional_units_to_cross_cutoff_min":int(p.additional_units_to_cross_cutoff.min()),"primary_additional_units_to_cross_cutoff_median":float(p.additional_units_to_cross_cutoff.median())}
    p.assign(Country=p.population_label)[["Country","who_region","income_label","cases_2024","cases_2050","growth","units","last_update","density","support_7","additional_units_to_cross_cutoff"]].to_csv(OUT/"primary_country_details.csv")
    d.to_csv(OUT/"country_data.csv")
    checks_file={"checks":checks,"passed":sum(c["passed"] for c in checks),"total":len(checks),"legacy_all_cancer_error":legacy_all,"statistics":statistics,"version":version,"software":{"python":__import__('platform').python_version(),"numpy":np.__version__,"pandas":pd.__version__,"scipy":__import__('scipy').__version__}}
    (OUT/"audit_results.json").write_text(json.dumps(checks_file,indent=2),encoding="utf-8")
    provenance=[{"file":str(x.relative_to(ROOT)),"sha256":hashlib.sha256(x.read_bytes()).hexdigest()} for x in dict.fromkeys(inputs)]
    (OUT/"input_checksums.json").write_text(json.dumps(provenance,indent=2),encoding="utf-8")
    print(json.dumps({"checks":len(checks),"statistics":statistics,"version":version,"legacy_all_cancer_error":legacy_all},indent=2))
    print(sensitivity.to_string(index=False))


if __name__ == "__main__":
    main()
