"""Reproducible retrospective virtual-monitoring experiment on public Graz data.
No field measurements are generated: sparse schedules/outages are experimental masks.
"""
from pathlib import Path
from importlib.metadata import version
import sys, json, hashlib
ROOT=Path(__file__).resolve().parent
import numpy as np
import pandas as pd

RAW=ROOT/'graz/raw'
OUT=ROOT/'analysis'

SEED=20261101

def mass_kg(q,c,method='linear',dt=60):
    """Exact integral for linear Q and linear or left-step C on each segment."""
    if method=='previous':
        return float(np.sum(dt*(q[:-1]+q[1:])*0.5*c[:-1])/1e6)
    return float(np.sum(dt*(2*q[:-1]*c[:-1]+q[:-1]*c[1:]+q[1:]*c[:-1]+2*q[1:]*c[1:])/6)/1e6)

def read_series(rel,col):
    d=pd.read_csv(RAW/rel,usecols=['time',col])
    t=pd.to_datetime(d.time,utc=True)
    if t.duplicated().any(): raise ValueError('Duplicate timestamps')
    s=pd.Series(d[col].to_numpy(float),index=t,name=col).sort_index()
    return s.where(s>=0)

def at_minutes(s,grid,maxgap):
    """Interpolate ONLY between valid observations separated by <= maxgap minutes.
    Exact observations remain valid; original long gaps never become reference.
    """
    valid=s.loc[grid[0]-pd.Timedelta(minutes=maxgap):grid[-1]+pd.Timedelta(minutes=maxgap)].dropna()
    if valid.empty: return np.full(len(grid),np.nan),np.zeros(len(grid),dtype=bool)
    origin=grid[0].value
    x=(valid.index.as_unit('ns').asi8-origin)/60_000_000_000
    z=(grid.as_unit('ns').asi8-origin)/60_000_000_000
    assert np.all(np.diff(x)>0), 'Minute timestamp conversion must be strictly increasing'
    y=valid.to_numpy()
    j=np.searchsorted(x,z)
    left=np.clip(j-1,0,len(x)-1); right=np.clip(j,0,len(x)-1)
    exact=x[right]==z
    good=exact|((j>0)&(j<len(x))&((x[right]-x[left])<=maxgap))
    vals=np.interp(z,x,y)
    vals[~good]=np.nan
    return vals,exact

def candidate_events(separation=6,tail=3):
    r=read_series('precipitation data/0033901601.csv','0033901601').dropna()
    r=r[r>0]  # Only positive bucket increments define rain groups.
    groups=(r.index.to_series().diff().dt.total_seconds().fillna(1e10)>separation*3600).cumsum().to_numpy()
    events=[]
    for gid in np.unique(groups):
        rain=r.iloc[np.flatnonzero(groups==gid)]
        if rain.sum()<2-1e-8: continue
        start=rain.index[0].floor('min'); end=rain.index[-1].ceil('min')+pd.Timedelta(hours=tail)
        duration=(end-start).total_seconds()/60
        if duration>2880: continue
        events.append(dict(event_id=start.strftime('%Y%m%dT%H%M'),start=start,end=end,
                           rainfall_mm=float(rain.sum()),duration_min=duration,
                           rain_duration_min=(rain.index[-1]-rain.index[0]).total_seconds()/60))
    return events

def build_events(maxgap=6,separation=6,tail=3,saturation=2400):
    OUT.mkdir(parents=True,exist_ok=True)
    q=read_series('hydraulic data/inflow-analogue.csv','flow_rate')
    cod=read_series('pollutant data/calibrated.csv','COD')
    tss=read_series('pollutant data/spectrometer.csv','TSS')
    native=dict(Q=q,COD=cod,TSS=tss)
    entries=[]; retained=[]
    for e in candidate_events(separation,tail):
        a=e['start']-pd.Timedelta(minutes=60); b=e['end']+pd.Timedelta(minutes=60)
        grid=pd.date_range(a,b,freq='min')
        data={}; exact={}
        for label,s in native.items(): data[label],exact[label]=at_minutes(s,grid,maxgap)
        e=e.copy()
        for label in native: e[label+'_support']=float(np.isfinite(data[label]).mean())
        support=all(np.isfinite(v).all() for v in data.values())
        e['flow_saturation_flag']=bool(np.any(data['Q']>=saturation))
        e['retained']=support and not e['flow_saturation_flag']
        e['reason']='retained' if e['retained'] else ';'.join([k+' unsupported' for k,v in data.items() if not np.isfinite(v).all()]+(['flow saturation'] if e['flow_saturation_flag'] else []))
        if e['retained']:
            sl=slice(60,len(grid)-60)
            Q=data['Q'][sl]
            e['Q_peak_Ls']=float(Q.max()); e['Q_min_Ls']=float(Q.min())
            e['volume_m3']=float(np.trapezoid(Q,dx=60)/1000)
            for p in ['COD','TSS']:
                C=data[p][sl]
                e[p+'_mass_kg']=mass_kg(Q,C)
                e[p+'_max_mgL']=float(C.max())
                e[p+'_native_fraction']=float(exact[p][sl].mean())
            e['year']=e['start'].year
            retained.append((e,grid,data))
        entries.append(e)
    pd.DataFrame(entries).to_csv(OUT/f'event_audit_gap{maxgap}_sep{separation}_tail{tail}_sat{saturation}.csv',index=False)
    return retained,entries

def run_event(e,grid,data,quick=False,shift=0):
    """Flow kept intact. Linear is retrospective; carry-forward uses previous sample."""
    n=int(e['duration_min']); ix=np.arange(60,60+n+1); x=np.arange(len(grid))
    Q=data['Q'][ix]
    qpk=int(np.argmax(Q))
    rng=np.random.default_rng(SEED+int(e['start'].timestamp())%1_000_000)
    randoms=rng.random(30)
    scenarios=[('none',0,'none'),('random30',30,'random'),('random60',60,'random'),
               ('random120',120,'random'),('peakQ60',60,'peak'),('early60',60,'early')]
    if quick: scenarios=[scenarios[0],scenarios[2],scenarios[4]]
    rows=[]
    for p in ['COD','TSS']:
        full=data[p]
        if shift: full=np.interp(x+shift,x,full)
        C=full[ix]
        ref=mass_kg(Q,C)
        for delta in [3,6,15,30,60]:
            for scenario,duration,position in scenarios:
                # Every minute phase for no-outage; 30 paired phase/placement replicates otherwise.
                phases=list(range(delta)) if duration==0 else [int(np.floor(v*delta)) for v in np.roll(randoms,7)]
                for rep,phase in enumerate(phases):
                    original_keep=np.arange(phase,len(grid),delta)
                    keep=original_keep.copy()
                    if duration:
                        if position=='random': left=int(np.floor(randoms[rep]*(n-duration+1)))
                        elif position=='peak': left=int(np.clip(qpk-duration//2,0,n-duration))
                        else: left=0
                        keep=keep[~((keep>=60+left)&(keep<60+left+duration))]
                    for method in ['linear','previous']:
                        # Buffers cover one full maximum sampling interval. No forced boundary samples.
                        assert keep[0]<=ix[0] and keep[-1]>=ix[-1]
                        if method=='linear':
                            est=np.interp(ix,keep,full[keep])
                            baseline=np.interp(ix,original_keep,full[original_keep])
                        else:
                            est=full[keep[np.searchsorted(keep,ix,side='right')-1]]
                            baseline=full[original_keep[np.searchsorted(original_keep,ix,side='right')-1]]
                        mass=mass_kg(Q,est,method)
                        err=100*(mass-ref)/ref
                        peakerr=100*((Q*est).max()-(Q*C).max())/(Q*C).max()
                        base_err=100*(mass_kg(Q,baseline,method)-ref)/ref
                        rows.append((e['event_id'],p,delta,scenario,method,rep,err,abs(err),peakerr,abs(base_err),abs(err)-abs(base_err),ref))
    return rows

def summarize(rep,name):
    keys=['event_id','pollutant','interval_min','scenario','method']
    ev=rep.groupby(keys).agg(mean_abs_pct=('absolute_error_pct','mean'),
                             mean_signed_pct=('signed_error_pct','mean'),
                             p90_abs_pct=('absolute_error_pct',lambda x:x.quantile(.9)),
                             probability_over10=('absolute_error_pct',lambda x:float((x>10).mean())),
                             probability_over20=('absolute_error_pct',lambda x:float((x>20).mean())),
                             mean_peak_bias_pct=('peak_error_pct','mean'),
                             paired_outage_penalty_pp=('outage_penalty_pp','mean'),
                             reference_mass_kg=('reference_mass_kg','first')).reset_index()
    ev.to_csv(OUT/f'{name}_event_metrics.csv',index=False)
    summaries=[]
    for key,g in ev.groupby(['pollutant','interval_min','scenario','method']):
        vals=g.mean_abs_pct.to_numpy(); rng=np.random.default_rng(SEED)
        boots=np.mean(vals[rng.integers(0,len(vals),size=(2000,len(vals)))],axis=1)
        summaries.append((*key,len(g),float(vals.mean()),float(np.median(vals)),
            float(g.mean_signed_pct.mean()),float(g.probability_over10.mean()),
            float(g.probability_over20.mean()),float(g.mean_peak_bias_pct.mean()),
            *np.quantile(boots,[.025,.975]).tolist(),float(np.average(g.mean_signed_pct,weights=g.reference_mass_kg)),
            float(np.average(g.mean_abs_pct,weights=g.reference_mass_kg)),float(g.paired_outage_penalty_pp.mean())))
    cols=['pollutant','interval_min','scenario','method','n_events','mean_abs_pct','median_event_mean_abs_pct',
          'mean_signed_pct','probability_over10','probability_over20','mean_peak_bias_pct','ci_low','ci_high',
          'load_weighted_signed_pct','load_weighted_abs_pct','paired_outage_penalty_pp']
    summary=pd.DataFrame(summaries,columns=cols)
    summary.to_csv(OUT/f'{name}_summary.csv',index=False)
    return ev,summary

def main():
    events,audit=build_events()
    print('CANDIDATE',len(audit),'RETAINED',len(events),flush=True)
    if len(events)<10: raise RuntimeError('Too few reference events; revise design transparently')
    pd.DataFrame([e for e,_,_ in events]).to_csv(OUT/'events.csv',index=False)
    rows=[]
    for i,(e,g,d) in enumerate(events):
        rows.extend(run_event(e,g,d))
        if i%10==0: print('event',i,flush=True)
    columns=['event_id','pollutant','interval_min','scenario','method','replicate',
             'signed_error_pct','absolute_error_pct','peak_error_pct','matched_baseline_abs_pct','outage_penalty_pp','reference_mass_kg']
    rep=pd.DataFrame(rows,columns=columns)
    rep.to_csv(OUT/'replicate_metrics.csv.gz',index=False,compression='gzip')
    ev,summary=summarize(rep,'main')
    strict,sa=build_events(maxgap=3)
    strictids=[e['event_id'] for e,_,_ in strict]
    if strictids: summarize(rep[rep.event_id.isin(strictids)],'strict3')
    sat, _=build_events(saturation=2300)
    satids=[e['event_id'] for e,_,_ in sat]
    if satids: summarize(rep[rep.event_id.isin(satids)],'saturation2300')
    # Alternative event-window definitions, with full independent extraction.
    for tag,sep,tail in [('tail2',6,2),('sep8',8,3)]:
        es,au=build_events(separation=sep,tail=tail)
        rs=[]
        for e,g,d in es: rs.extend(run_event(e,g,d,quick=True))
        if rs:
            sr=pd.DataFrame(rs,columns=columns)
            summarize(sr,tag)
            common=sorted(set(sr.event_id)&set(rep.event_id))
            summarize(sr[sr.event_id.isin(common)],tag+'_common')
            summarize(rep[rep.event_id.isin(common)],tag+'_main_common')
    # Sampling-error stability: shifted concentrations define BOTH estimates and reference.
    # This is not clock-offset mass bias against an unshifted reference.
    for lag in [-3,3]:
        rs=[]
        for e,g,d in events: rs.extend(run_event(e,g,d,quick=True,shift=lag))
        summarize(pd.DataFrame(rs,columns=columns),'shift'+str(lag))
    sensitivity=[]
    for e,g,d in events:
        sl=slice(60,len(g)-60)
        q=d['Q'][sl]
        for p in ['COD','TSS']:
            flux=q*d[p][sl]
            full=mass_kg(q,d[p][sl])
            # 3-minute quadrature, including actual final boundary time.
            ids=np.unique(np.r_[np.arange(0,len(flux),3),len(flux)-1])
            coarse=mass_kg(q[ids],d[p][sl][ids],dt=np.diff(ids)*60)
            sensitivity.append((e['event_id'],p,100*(coarse-full)/full))
    pd.DataFrame(sensitivity,columns=['event_id','pollutant','integration_change_pct']).to_csv(OUT/'integration_sensitivity.csv',index=False)
    meta={'seed':SEED,'dataset':'10.5281/zenodo.15783832','archive_md5':hashlib.md5((ROOT/'graz/data-set.zip').read_bytes()).hexdigest(),
          'candidate_events':len(audit),'retained_events':len(events),'strict3_events':len(strictids),
          'replicates':len(rep),'native_maxgap_minutes':6,'event_separation_hours':6,'event_tail_hours':3,
          'min_rainfall_mm':2,'event_duration_hours':[3,48],'buffer_minutes':60,
          'flow_source':'inflow-analogue.csv','flow_saturation_exclusion_Ls':2400,
          'note':'One-minute numerical support grid is interpolated from native 1/3-minute measurements; not independent 1-minute field observations.',
          'timestamp_audit':timestamp_audit(),'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,'matplotlib':version('matplotlib')}
    (OUT/'provenance.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
    print(json.dumps(meta,indent=2),flush=True)
    print(summary[(summary.method=='linear')&(summary.interval_min.isin([3,15,60]))].to_string(index=False),flush=True)

def timestamp_audit():
    result={}
    for label,rel,col in [('Q','hydraulic data/inflow-analogue.csv','flow_rate'),('COD','pollutant data/calibrated.csv','COD'),('TSS','pollutant data/spectrometer.csv','TSS')]:
        s=read_series(rel,col)
        result[label]={'rows':len(s),'subminute_timestamps':int((s.index.as_unit('ns').asi8%60_000_000_000!=0).sum())}
    return result

def flow_predict(q_samples,c_samples,q_target):
    qm=q_samples.mean();qs=q_samples.std()
    if qs<1e-12:return np.full(len(q_target),max(0.,c_samples.mean()))
    coef=np.linalg.lstsq(np.column_stack([np.ones(len(q_samples)),(q_samples-qm)/qs]),c_samples,rcond=None)[0]
    return np.maximum(0,coef[0]+coef[1]*(q_target-qm)/qs)

def verify_calculations():
    idx=pd.date_range('2010-01-01',periods=4,freq='3min',tz='UTC')
    grid=pd.date_range(idx[0],idx[-1],freq='min')
    vals,exact=at_minutes(pd.Series([0.,3.,6.,9.],index=idx),grid,3)
    assert np.allclose(vals,np.arange(10)) and exact.sum()==4
    gapped=pd.Series([0.,9.],index=idx[[0,-1]])
    v,_=at_minutes(gapped,grid,3)
    assert np.isfinite(v).sum()==2
    assert np.isclose(mass_kg(np.full(61,100.),np.full(61,200.)),72.)
    assert np.isclose(mass_kg(np.ones(7),np.array([0.,0.,0.,10.,10.,10.,10.]),'previous',dt=1)*1e6,30.)
    assert np.isclose(mass_kg(np.array([0.,1.]),np.array([0.,1.]),dt=1)*1e6,1/3)
    print('Timestamp, gap and mass-unit checks passed',flush=True)

if __name__=='__main__':
    verify_calculations()
    main()
