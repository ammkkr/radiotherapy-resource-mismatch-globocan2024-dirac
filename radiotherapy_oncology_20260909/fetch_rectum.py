"""Complete the defined site set and verify a contemporaneous anchor."""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import json
from pathlib import Path
from urllib.request import Request, urlopen

import pandas as pd

import os
ROOT = Path(os.environ["RADIOTHERAPY_SOURCE_ARCHIVE"]).resolve()
OUT = ROOT / "data/raw/radiotherapy_review_rectum_20260909"


def fetch(job):
    version, measure, chunk, ids = job
    path = OUT / f"v{version}_type{measure}_chunk{chunk:02d}.json"
    url = f"https://gco.iarc.who.int/gateway_prod/api/globocan/v3/{version}/data/prediction/{measure}/0/{'_'.join(map(str,ids))}/9_15/?prediction_annual=1&ages_group=0_17"
    if path.exists():
        payload = json.loads(path.read_text(encoding="utf-8"))
    else:
        with urlopen(Request(url, headers={"Accept":"application/json","User-Agent":"radiotherapy-resource-review/2.0"}),timeout=45) as response:
            payload = json.load(response)
        path.write_text(json.dumps(payload),encoding="utf-8")
    rows = payload["dataset"]
    if len(rows) != len(ids)*2*7:
        raise ValueError(f"Incomplete API response: {url}: {len(rows)}")
    return {"url":url,"file":path.name,"rows":len(rows)}


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    jobs=[]
    for version in [2022,2024]:
        p=pd.read_csv(ROOT/f"data/interim/gco_cancer_tomorrow_{version}_predictions_long.csv")
        ids=sorted(p[p.country_iso3.notna()].country_code.unique().tolist())
        for measure in [0,1]:
            for chunk,start in enumerate(range(0,len(ids),10),1):
                jobs.append((version,measure,chunk,ids[start:start+10]))
    with ThreadPoolExecutor(max_workers=4) as pool:
        results=list(pool.map(fetch,jobs))
    (OUT/"fetch_log.json").write_text(json.dumps({"accessed_utc":datetime.now(timezone.utc).isoformat(),"queries":results},indent=2),encoding="utf-8")
    print(json.dumps({"queries":len(results),"rows":sum(x['rows'] for x in results)}))


if __name__=="__main__":
    main()
