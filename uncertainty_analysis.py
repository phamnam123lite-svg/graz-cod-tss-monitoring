from pathlib import Path
import numpy as np,pandas as pd
R=Path(__file__).resolve().parent;A=R/'supplemental_analysis'
rng=np.random.default_rng(20261101)
def ci(v):
 v=np.asarray(v);return np.quantile(v[rng.integers(0,len(v),(2000,len(v)))].mean(axis=1),[.025,.975])
ev=pd.read_csv(R/'analysis/main_event_metrics.csv');sel=pd.read_csv(A/'selection_comparison_events.csv')
v=ev[(ev.interval_min==15)&(ev.method=='linear')].merge(sel[['event_id','season']],on='event_id')
rows=[]
for (pol,sc,s),g in v.groupby(['pollutant','scenario','season']):rows.append([pol,sc,s,len(g),g.mean_abs_pct.mean(),*ci(g.mean_abs_pct)])
season=pd.DataFrame(rows,columns=['pollutant','scenario','season','n','mean','low','high']);season.to_csv(A/'season_bootstrap_ci.csv',index=False)
scale=pd.read_csv(A/'interval_correlation_scale.csv');rows=[]
for (pol,it),g in scale.groupby(['pollutant','interval_min']):
 x=g.interval_over_tau.to_numpy();y=g.mean_abs_pct.to_numpy();rho=pd.Series(x).corr(pd.Series(y),method='spearman') if False else pd.Series(x).rank().corr(pd.Series(y).rank())
 vals=[]
 for k in range(2000):
  ix=rng.integers(0,len(x),len(x));vals.append(pd.Series(x[ix]).rank().corr(pd.Series(y[ix]).rank()))
 rows.append([pol,it,len(x),rho,*np.nanquantile(vals,[.025,.975])])
corr=pd.DataFrame(rows,columns=['pollutant','interval_min','n','rho','low','high']);corr.to_csv(A/'spearman_bootstrap_ci.csv',index=False)
