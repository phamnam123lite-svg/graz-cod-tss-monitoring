from pathlib import Path
import sys,json
R=Path(__file__).parent;sys.path.insert(0,str(R/'pylibs'))
import numpy as np,pandas as pd
import analyse_graz as b
O=R/'supplemental_analysis';b.OUT=O
events,_=b.build_events();series={p:b.read_series(path,col) for p,path,col in [('Q','hydraulic data/inflow-analogue.csv','flow_rate'),('COD','pollutant data/calibrated.csv','COD'),('TSS','pollutant data/spectrometer.csv','TSS')]}
rows=[];audit=[]
for e,grid,data in events:
    n=int(e['duration_min']);ix=np.arange(60,60+n+1)
    masks={p:b.at_minutes(s,grid,6)[1] for p,s in series.items()}
    native=ix[np.logical_and.reduce([masks[p][ix] for p in series])]
    if len(native)<3:continue
    # Retain only a complete chain of shared observed knots no more than six minutes apart.
    if np.max(np.diff(native))>6:continue
    match=np.arange(native[0],native[-1]+1)
    rng=np.random.default_rng(b.SEED+int(e['start'].timestamp())%1_000_000);u=rng.random(30)
    qpk=int(np.argmax(data['Q'][ix]))
    audit.append([e['event_id'],len(native),float(native[-1]-native[0]),float(e['duration_min']),int(np.max(np.diff(native)))])
    for pol in ['COD','TSS']:
        ref=b.mass_kg(data['Q'][native],data[pol][native],dt=np.diff(native)*60)
        ref_match=b.mass_kg(data['Q'][match],data[pol][match])
        for sc in ['none','random60','peakQ60']:
            phases=range(15) if sc=='none' else (np.roll(u,7)*15).astype(int)
            for rep,phase in enumerate(phases):
                keep=np.arange(phase,len(grid),15)
                if sc!='none':
                    left=int(np.floor(u[rep]*(n-60+1))) if sc=='random60' else int(np.clip(qpk-30,0,n-60))
                    keep=keep[~((keep>=60+left)&(keep<60+left+60))]
                est=np.interp(native,keep,data[pol][keep]);em=np.interp(match,keep,data[pol][keep])
                err=100*(b.mass_kg(data['Q'][native],est,dt=np.diff(native)*60)-ref)/ref
                err_match=100*(b.mass_kg(data['Q'][match],em)-ref_match)/ref_match
                rows.append([e['event_id'],pol,sc,rep,abs(err),abs(err_match),100*(ref-ref_match)/ref_match])
df=pd.DataFrame(rows,columns=['event_id','pollutant','scenario','replicate','native_knot_abs_error_pct','matched_grid_abs_error_pct','reference_mass_change_pct'])
df.to_csv(O/'native_knot_replicates.csv.gz',index=False,compression='gzip')
ev=df.groupby(['event_id','pollutant','scenario']).mean(numeric_only=True).reset_index();ev.to_csv(O/'native_knot_event_metrics.csv',index=False)
summary=ev.groupby(['pollutant','scenario']).agg(n_events=('event_id','count'),native_knot_abs_error_pct=('native_knot_abs_error_pct','mean'),matched_grid_abs_error_pct=('matched_grid_abs_error_pct','mean')).reset_index()
summary.to_csv(O/'native_knot_summary.csv',index=False)
pd.DataFrame(audit,columns=['event_id','shared_native_points','matched_duration_min','original_duration_min','largest_native_gap_min']).to_csv(O/'native_knot_support.csv',index=False)
print(summary.to_string(index=False));print('native fraction range',pd.read_csv(R/'analysis/events.csv').COD_native_fraction.agg(['min','median','max']).to_dict());print('matched windows',len(audit))
