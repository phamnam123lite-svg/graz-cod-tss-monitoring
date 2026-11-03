"""Source-record audit and disclosed event-exclusion/native-sample sensitivity.

The flagged event is NOT classified as confirmed instrument failure. No
post-hoc concentration threshold is applied to the primary 66-event series.
"""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
import analyse_graz as b
R=Path(__file__).resolve().parent
FLAGGED='20110604T1250'

def mean_ci(v):
    v=np.asarray(v,float);rng=np.random.default_rng(b.SEED)
    draws=v[rng.integers(0,len(v),(2000,len(v)))].mean(axis=1)
    return [v.mean(),*np.quantile(draws,[.025,.975])]

def nearest_native(scheduled,native,tolerance=1):
    """Nearest exact timestamp within tolerance; ties choose the earlier point."""
    j=np.searchsorted(native,scheduled)
    l=native[np.clip(j-1,0,len(native)-1)];r=native[np.clip(j,0,len(native)-1)]
    chosen=np.where(abs(scheduled-l)<=abs(scheduled-r),l,r)
    return np.unique(chosen[abs(chosen-scheduled)<=tolerance])

def main():
    A=R/'analysis';O=R/'supplemental_analysis';O.mkdir(exist_ok=True)
    inventory={p.relative_to(b.RAW).as_posix():pd.read_csv(p,nrows=0).columns.tolist() for p in b.RAW.rglob('*.csv') if 'evaluations' not in p.parts}
    specs=[('Q','hydraulic data/inflow-analogue.csv','flow_rate'),('COD','pollutant data/calibrated.csv','COD'),('TSS','pollutant data/spectrometer.csv','TSS'),('global_COD','pollutant data/spectrometer.csv','COD'),('temperature','pollutant data/spectrometer.csv','temperature')]
    series={label:b.read_series(rel,col).rename(label) for label,rel,col in specs}
    trace=pd.concat([s.loc['2011-06-04 14:30Z':'2011-06-04 16:30Z'] for s in series.values()],axis=1)
    trace.index.name='time_UTC';trace.to_csv(O/'source_qc_20110604.csv')
    labs=pd.read_csv(b.RAW/'lab data/samples.csv');times=pd.to_datetime(labs.time,utc=True)
    labs_same_day=int((times.dt.date==pd.Timestamp('2011-06-04').date()).sum())
    meta={'event_id':FLAGGED,'raw_archive_md5':hashlib.md5((R/'graz/data-set.zip').read_bytes()).hexdigest(),
          'source_columns':inventory,'timestamp_audit':b.timestamp_audit(),
          'source_quality_flag_columns_found':[c for cols in inventory.values() for c in cols if any(k in c.lower() for k in ['flag','quality','status'])],
          'lab_samples_on_2011_06_04':labs_same_day,
          'peak_UTC':'2011-06-04T15:26:00Z','at_peak':trace.loc['2011-06-04 15:26Z'].to_dict(),
          '30_minutes_before_peak_grid_values':{p:float(b.at_minutes(s,pd.DatetimeIndex([pd.Timestamp('2011-06-04 14:56Z')]),6)[0][0]) for p,s in series.items()},
          'before_peak_note':'14:56 UTC is a supported interpolated grid point, not an exact raw-observation timestamp; source_qc_20110604.csv contains only original timestamps.',
          'classification':'unresolved possible probe artefact; source values alone cannot distinguish failure from a real concentration change',
          'policy':'Retain original nonnegative measurements in primary analysis; separately exclude the entire flagged event. No unvalidated concentration threshold is imposed.',
          'source_publication':'https://doi.org/10.3390/data11030045',
          'limitation':'No observation-level quality flags, dated maintenance log, or same-day laboratory samples were found in the supplied archive.'}
    (O/'source_qc_report.json').write_text(json.dumps(meta,indent=2,allow_nan=False))
    ev=pd.read_csv(A/'main_event_metrics.csv');ae=pd.read_csv(O/'additional_event_metrics.csv')
    tables=[];exposure=[];penalty=[];contrasts=[]
    for name,omit in [('all_66',False),('exclude_flagged_65',True)]:
        v=ev[ev.event_id!=FLAGGED] if omit else ev
        v=v[(v.interval_min==15)&(v.method=='linear')]
        for (p,sc),g in v.groupby(['pollutant','scenario']):
            x=g.mean_abs_pct.to_numpy()
            tables.append([name,p,sc,len(g),*mean_ci(x),*np.quantile(x,[.5,.9,.95])])
            penalty.append([name,p,sc,len(g),*mean_ci(g.paired_outage_penalty_pp)])
        a=ae[ae.event_id!=FLAGGED] if omit else ae
        a=a[a.method=='linear']
        for (p,sc),g in a.groupby(['pollutant','scenario']):
            exposure.append([name,p,sc,len(g),*mean_ci(g.mean_abs_pct),g.gap_volume_share_pct.mean(),g.gap_mass_share_pct.mean()])
        for p in ['COD','TSS']:
            t=a[a.pollutant==p].pivot(index='event_id',columns='scenario',values='mean_abs_pct')
            for sc in ['random60','flowweighted60']:
                contrasts.append([name,p,'peak_minus_'+sc,len(t),*mean_ci(t.peakQ60-t[sc])])
    pd.DataFrame(tables,columns=['sample','pollutant','scenario','n','mean','low','high','median','p90','p95']).to_csv(O/'qc_exclusion_table2.csv',index=False)
    pd.DataFrame(exposure,columns=['sample','pollutant','scenario','n','mean','low','high','gap_volume_share_pct','gap_mass_share_pct']).to_csv(O/'qc_exclusion_table3.csv',index=False)
    pd.DataFrame(penalty,columns=['sample','pollutant','scenario','n','mean_penalty_pp','low','high']).to_csv(O/'qc_exclusion_penalties.csv',index=False)
    pd.DataFrame(contrasts,columns=['sample','pollutant','contrast','n','mean_difference_pp','low','high']).to_csv(O/'qc_exclusion_contrasts.csv',index=False)
    b.OUT=O;events,_=b.build_events();rows=[]
    for e,grid,data in events:
        ix=np.arange(60,len(grid)-60)
        for p in ['COD','TSS']:
            native=np.flatnonzero(b.at_minutes(series[p],grid,6)[1])
            ref=b.mass_kg(data['Q'][ix],data[p][ix])
            for interval in [3,6,15]:
                for phase in range(interval):
                    scheduled=np.arange(phase,len(grid),interval)
                    keep=nearest_native(scheduled,native)
                    assert keep[0]<=ix[0] and keep[-1]>=ix[-1]
                    estimate=np.interp(ix,keep,data[p][keep])
                    err=100*(b.mass_kg(data['Q'][ix],estimate)-ref)/ref
                    rows.append([e['event_id'],p,interval,phase,abs(err),len(scheduled),len(keep)])
    df=pd.DataFrame(rows,columns=['event_id','pollutant','interval_min','phase','absolute_error_pct','scheduled_n','unique_native_samples_n'])
    df.to_csv(O/'native_sample_replicates.csv',index=False)
    em=df.groupby(['event_id','pollutant','interval_min']).mean(numeric_only=True).reset_index()
    em.to_csv(O/'native_sample_event_metrics.csv',index=False)
    result=[]
    for (p,it),g in em.groupby(['pollutant','interval_min']):result.append([p,it,len(g),*mean_ci(g.absolute_error_pct)])
    pd.DataFrame(result,columns=['pollutant','interval_min','n','mean','low','high']).to_csv(O/'native_sample_summary.csv',index=False)
    print('Source QC, exclusion tables and native-sample comparison written.',flush=True)

if __name__=='__main__':main()
