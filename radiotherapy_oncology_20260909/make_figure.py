"""Plot audited 14-site data using the established three-panel composition."""
from pathlib import Path
import json
import sys

import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from PIL import Image

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import make_radiotherapy_resource_mismatch_figure as base

ROOT=Path(__file__).resolve().parent
DATA=ROOT/"data"
OUT=ROOT/"figures"


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    base.apply_publication_style()
    plt.rcParams.update({"font.size":7,"axes.labelsize":7,"xtick.labelsize":6.2,"ytick.labelsize":6.2})
    d=pd.read_csv(DATA/"country_data.csv")
    v=pd.read_csv(DATA/"version_comparison.csv")
    d=d.rename(columns={"population_label":"gco_country","units":"mv_therapy_units","density":"mv_units_per_1000_rt_relevant_cases_2050","cases_2050":"rt_relevant_cases_2050"})
    d["dirac_matched"]=d.mv_therapy_units.notna().astype(int)
    d["high_growth_low_resource"]=d.primary.astype(int)
    d["current_screen_positive"]=d.primary
    d["selected_site_cases_per_mv_unit_2050"]=d.rt_relevant_cases_2050/d.mv_therapy_units
    d["unit_pressure_rank"]=d.selected_site_cases_per_mv_unit_2050.rank(ascending=False,method="first")
    d=d.merge(v,on="country_iso3",validate="one_to_one")
    d["gco2022_common_year_growth_percentile"]=d.old_percentile
    d["gco2024_common_year_growth_percentile"]=d.new_percentile
    d["gco2024_selected_cases_2050"]=d.rt_relevant_cases_2050
    d["common_year_burden_version_transition"]=np.select([~d.common,d.old_selected & d.new_selected,~d.old_selected & d.new_selected,d.old_selected & ~d.new_selected],["Not comparable","Retained","New in GCO version 2024","No longer screen-positive"],default="Not screen-positive in either version")
    palette=base.PALETTE
    cmap=LinearSegmentedColormap.from_list("density",[palette['signal'],palette['signal_soft'],palette['neutral_pale'],palette['blue_mid'],palette['blue_main']])
    norm=LogNorm(.01,5,clip=True)
    fig=plt.figure(figsize=(8.2,6.7))
    a=fig.add_axes([.055,.455,.925,.51])
    b=fig.add_axes([.075,.14,.365,.24])
    c=fig.add_axes([.59,.14,.345,.24])
    d["mv_units_per_1000_selected_site_cases_2050"]=d.mv_units_per_1000_rt_relevant_cases_2050
    base.draw_world_map(a,d,cmap,norm)
    unknown_high = (d.mv_therapy_units.isna() & d.growth.ge(d.loc[d.dirac_matched.eq(1),'growth'].quantile(.75))).sum()
    a.text(.025,.17,f"{d.dirac_matched.sum()} matched countries\n{d.primary.sum()} met both thresholds\n{unknown_high} high-growth, resources unknown",transform=a.transAxes,fontsize=6.3,va="bottom",fontweight="bold",bbox={"facecolor":"white","edgecolor":"none","alpha":.8,"pad":1})
    a.legend(handles=[Patch(facecolor="white",edgecolor=palette['signal_dark'],label="Both thresholds"),Patch(facecolor=palette['missing'],edgecolor=palette['neutral_mid'],hatch="////",label="DIRAC unmatched")],loc="lower left",bbox_to_anchor=(.018,.015),fontsize=6.2)
    cax=a.inset_axes([.335,-.095,.36,.035])
    sm=ScalarMappable(norm=norm,cmap=cmap)
    bar=fig.colorbar(sm,cax=cax,orientation="horizontal",ticks=[.01,.05,.1,.5,1,5])
    bar.ax.set_xticklabels(['0.01','0.05','0.1','0.5','1','5'])
    bar.ax.tick_params(labelsize=6,pad=1,length=2)
    bar.set_label("Recorded MV units per 1000 selected-site cases in 2050",fontsize=6.2,labelpad=2)
    base.draw_panel_b(b,d)
    b.set(xlim=(-2,102),ylim=(-2,102))
    for t in b.texts:
        if t.get_text().startswith("Common 2025"):
            t.set_text("Common 2025 baseline; fixed July 2026 DIRAC extract.\nPoint area scales with log10 projected 2050 cases.")
            t.set_fontsize(5.8)
            t.set_position((0,-.39))
    b.set_xlabel("Projected 2025-2050 growth percentile\nGLOBOCAN 2022 release",labelpad=5)
    b.set_ylabel("Projected 2025-2050 growth percentile\nGLOBOCAN 2024 release",labelpad=5)
    bounds=pd.read_csv(DATA/"missingness_bounds.csv")
    region=bounds[bounds.Group.eq('WHO region') & bounds.Category.notna()].sort_values(['Lower (%)','Upper (%)'])
    y=np.arange(len(region))
    c.hlines(y,region['Lower (%)'],region['Upper (%)'],color=palette['neutral_light'],lw=2)
    c.scatter(region['Lower (%)'],y,s=34,c=palette['signal'],edgecolors="white",zorder=3)
    c.scatter(region['Upper (%)'],y,s=38,facecolors="white",edgecolors=palette['signal_dark'],zorder=4)
    for yi,(_,r) in enumerate(region.iterrows()):
        c.text(r['Upper (%)']+2,yi,f"{r['Observed selected']}/{r.Countries} to {r['Observed selected']+r['Unmatched high growth']}/{r.Countries}",va="center",fontsize=6)
    c.set(xlim=(-2,83),ylim=(-.6,5.6),yticks=y,yticklabels=region.Category,xticks=[0,20,40,60,80])
    c.set_xlabel("Countries meeting both fixed cutoffs (%)",labelpad=5)
    c.tick_params(axis="y",length=0)
    c.legend(handles=[Line2D([],[],marker='o',linestyle='none',color=palette['signal'],label='Observed lower bound'),Line2D([],[],marker='o',linestyle='none',markerfacecolor='white',color=palette['signal_dark'],label='Unmatched high-growth upper bound')],loc="lower right",fontsize=5.6)
    c.text(.0,-.39,"Equal weight per country; recorded counts held fixed.\nBounds address unmatched records, not all uncertainty.",transform=c.transAxes,fontsize=5.8,color=palette['neutral_mid'],va="top")
    for ax,label in [(a,'a'),(b,'b'),(c,'c')]:
        base.add_panel_label(ax,label,0,1.015)
    base.bolden_figure_text(fig)
    stem=OUT/"figure_1_radiotherapy_oncology_v2"
    for ext in ['png','tiff','pdf','svg']:
        kwargs={"dpi":600,"bbox_inches":"tight"}
        if ext=='tiff':kwargs['pil_kwargs']={"compression":"tiff_lzw"}
        if ext=='pdf':kwargs['metadata']={"Creator":None,"Producer":None,"CreationDate":None}
        fig.savefig(stem.with_suffix('.'+ext),**kwargs)
    # Export panels from their actual artists, not fixed pixel cuts.
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    for ax,label in [(a,'a'),(b,'b'),(c,'c')]:
        box=ax.get_tightbbox(renderer).transformed(fig.dpi_scale_trans.inverted()).expanded(1.02,1.02)
        for ext in ['png','tiff']:
            kwargs={"dpi":600,"bbox_inches":box}
            if ext=='tiff':kwargs['pil_kwargs']={"compression":"tiff_lzw"}
            path=OUT/f"figure_1{label}_radiotherapy_oncology_v2.{ext}"
            fig.savefig(path,**kwargs)
            if ext=='tiff':
                with Image.open(path) as im:rgb=im.convert('RGB')
                rgb.save(path,dpi=(600,600),compression='tiff_lzw')
    with Image.open(stem.with_suffix('.tiff')) as im:rgb=im.convert('RGB')
    rgb.save(stem.with_suffix('.tiff'),dpi=(600,600),compression='tiff_lzw')
    plt.close(fig)
    d.to_csv(DATA/"figure_country_source_data.csv",index=False)
    region.to_csv(DATA/"figure_region_source_data.csv",index=False)
    print(stem)


if __name__=="__main__":main()
