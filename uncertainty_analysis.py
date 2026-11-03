"""Event-bootstrap intervals with an independent fixed stream for every group."""
from pathlib import Path
import numpy as np,pandas as pd
R=Path(__file__).resolve().parent;A=R/'supplemental_analysis'
SEED=20261101
def ci(v):
    v=np.asarray(v);rng=np.random.default_rng(SEED)
    return np.quantile(v[rng.integers(0,len(v),(2000,len(v)))].mean(axis=1),[.025,.975])
def main():
    ev=pd.read_csv(R/'analysis/main_event_metrics.csv');sel=pd.read_csv(A/'selection_comparison_events.csv')
    v=ev[(ev.interval_min==15)&(ev.method=='linear')].merge(sel[['event_id','season']],on='event_id')
    rows=[]
    for (pol,sc,s),g in v.groupby(['pollutant','scenario','season']):
        rows.append([pol,sc,s,len(g),g.mean_abs_pct.mean(),*ci(g.mean_abs_pct)])
    pd.DataFrame(rows,columns=['pollutant','scenario','season','n','mean','low','high']).to_csv(A/'season_bootstrap_ci.csv',index=False)
    scale=pd.read_csv(A/'interval_correlation_scale.csv');rows=[]
    for (pol,it),g in scale.groupby(['pollutant','interval_min']):
        g=g.replace([np.inf,-np.inf],np.nan).dropna(subset=['acf_1e_minutes','mean_abs_pct'])
        x=g.acf_1e_minutes.to_numpy();y=g.mean_abs_pct.to_numpy()
        rho=pd.Series(x).rank().corr(pd.Series(y).rank());rng=np.random.default_rng(SEED);vals=[]
        for k in range(2000):
            ix=rng.integers(0,len(x),len(x));vals.append(pd.Series(x[ix]).rank().corr(pd.Series(y[ix]).rank()))
        rows.append([pol,it,len(x),rho,*np.nanquantile(vals,[.025,.975])])
    pd.DataFrame(rows,columns=['pollutant','interval_min','n','rho','low','high']).to_csv(A/'spearman_bootstrap_ci.csv',index=False)
if __name__=='__main__':main()
