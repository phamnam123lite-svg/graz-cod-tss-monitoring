"""Supplemental analyses; reuse the unchanged Graz source and original event definition."""
from pathlib import Path
import sys,json
import numpy as np,pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import analyse_graz as base

def main():
    R=Path(__file__).resolve().parent
    matplotlib.use('Agg')
    O=R/'supplemental_analysis';O.mkdir(exist_ok=True)
    O=R/'supplemental_analysis';O.mkdir(exist_ok=True)
    base.OUT=O
    events,audit=base.build_events()
    E=pd.DataFrame([e for e,_,_ in events])
    REP=pd.read_csv(R/'analysis/replicate_metrics.csv.gz')
    EV=pd.read_csv(R/'analysis/main_event_metrics.csv')
    SEED=base.SEED
    def ci_mean(a):
        a=np.asarray(a,float);rng=np.random.default_rng(SEED)
        b=a[rng.integers(0,len(a),size=(2000,len(a)))].mean(axis=1)
        return [float(a.mean()),*np.quantile(b,[.025,.975]).tolist()]
    rows=[]
    for (pol,sc),g in EV[(EV.interval_min==15)&(EV.method=='linear')].groupby(['pollutant','scenario']):
        a=g.mean_abs_pct.to_numpy()
        rows.append([pol,sc,len(a),a.mean(),*np.quantile(a,[.5,.9,.95]),*ci_mean(g.mean_peak_bias_pct)])
    dist=pd.DataFrame(rows,columns=['pollutant','scenario','n_events','mean_abs_pct','median_event_mean_abs_pct','p90_event_mean_abs_pct','p95_event_mean_abs_pct','mean_peak_bias_pct','peak_ci_low','peak_ci_high'])
    dist.to_csv(O/'distribution_and_peak_ci.csv',index=False)
    paired=[]
    for pol in ['COD','TSS']:
        g=EV[(EV.interval_min==15)&(EV.method=='linear')&(EV.pollutant==pol)].pivot(index='event_id',columns='scenario',values='mean_abs_pct')
        paired.append([pol,len(g),*ci_mean(g.peakQ60-g.random60)])
    pd.DataFrame(paired,columns=['pollutant','n_events','peak_minus_random_pp','ci_low','ci_high']).to_csv(O/'peak_random_paired_contrast.csv',index=False)
    highids=E.loc[(E.COD_native_fraction>.9)&(E.TSS_native_fraction>.9),'event_id']
    high=EV[(EV.event_id.isin(highids))&(EV.interval_min==15)&(EV.method=='linear')&(EV.scenario.isin(['none','random60','peakQ60']))]
    high_summary=high.groupby(['pollutant','scenario']).agg(n_events=('event_id','count'),mean_abs_pct=('mean_abs_pct','mean')).reset_index()
    high_summary.to_csv(O/'high_native_subset.csv',index=False)
    rows=[];acf=[];figure_checks=[]
    rows=[];acf=[];figure_checks=[]
    rows=[];acf=[];figure_checks=[]
    example_ids=[E.iloc[np.argmin(abs(E.Q_peak_Ls-E.Q_peak_Ls.median()))].event_id,E.sort_values('Q_peak_Ls').iloc[-1].event_id]
    for event_no,(e,grid,data) in enumerate(events):
        n=int(e['duration_min']);ix=np.arange(60,60+n+1);Q=data['Q'][ix]
        qpk=int(np.argmax(Q));rng=np.random.default_rng(SEED+int(e['start'].timestamp())%1_000_000)
        u=rng.random(30);phases=(np.roll(u,7)*15).astype(int)
        vol_segments=(Q[:-1]+Q[1:])*.5*60
        cumulative=np.r_[0,np.cumsum(vol_segments)]
        starts=np.arange(n-60+1);vol60=cumulative[starts+60]-cumulative[starts]
        cdf=np.cumsum(vol60)/vol60.sum()
        for pol in ['COD','TSS']:
            full=data[pol];C=full[ix];ref=base.mass_kg(Q,C)
            z=C-C.mean();den=float(z@z);maxlag=max(1,n//4)
            corr=np.array([float(z[:-lag]@z[lag:]/den) if den>0 else np.nan for lag in range(1,maxlag+1)])
            crosses=np.flatnonzero(corr<=np.exp(-1));tau=int(crosses[0]+1) if len(crosses) else np.nan
            acf.append([e['event_id'],pol,tau,maxlag,bool(len(crosses)),float(np.std(C)/np.mean(C))])
            for scenario in ['none','random60','peakQ60','flowweighted60']:
                phase_list=list(range(15)) if scenario=='none' else phases.tolist()
                for rep,phase in enumerate(phase_list):
                    original=np.arange(phase,len(grid),15);keep=original.copy()
                    left=None
                    if scenario!='none':
                        if scenario=='random60':left=int(np.floor(u[rep]*(n-60+1)))
                        elif scenario=='peakQ60':left=int(np.clip(qpk-30,0,n-60))
                        else:left=int(np.searchsorted(cdf,u[rep],side='right'))
                        keep=keep[~((keep>=60+left)&(keep<60+left+60))]
                    volume_share=mass_share=0.
                    if left is not None:
                        sl=slice(left,left+61)
                        volume_share=100*(cumulative[left+60]-cumulative[left])/cumulative[-1]
                        mass_share=100*base.mass_kg(Q[sl],C[sl])/ref
                    for method in ['linear','flow_regression']:
                        if method=='linear':
                            est=np.interp(ix,keep,full[keep]);baseline=np.interp(ix,original,full[original])
                        else:
                            def predict(k):
                                return base.flow_predict(data['Q'][k],full[k],Q)
                            est=predict(keep);baseline=predict(original)
                        error=100*(base.mass_kg(Q,est)-ref)/ref
                        base_error=100*(base.mass_kg(Q,baseline)-ref)/ref
                        rows.append([e['event_id'],pol,scenario,method,rep,phase,left,abs(error),error,abs(error)-abs(base_error),volume_share,mass_share,100*60/n])
            if e['event_id'] in example_ids:
                figure_checks.append({'event_id':e['event_id'],'pollutant':pol,'peak_time':str(grid[60+qpk]),'concentration_at_peak':float(C[qpk]),'concentration_30min_before':float(C[max(0,qpk-30)]),'concentration_30min_after':float(C[min(n,qpk+30)])})
        if event_no%15==0:print('additional event',event_no,flush=True)
    new=pd.DataFrame(rows,columns=['event_id','pollutant','scenario','method','replicate','phase','gap_start_min','absolute_error_pct','signed_error_pct','outage_penalty_pp','gap_volume_share_pct','gap_mass_share_pct','gap_duration_share_pct'])
    new.to_csv(O/'additional_replicates.csv.gz',index=False,compression='gzip')
    check=new[(new.method=='linear')&(new.scenario!='flowweighted60')].sort_values(['event_id','pollutant','scenario','replicate'])
    old=REP[(REP.interval_min==15)&(REP.method=='linear')&(REP.scenario.isin(['none','random60','peakQ60']))].sort_values(['event_id','pollutant','scenario','replicate'])
    assert len(check)==len(old)
    baseline_checks={}
    for col in ['absolute_error_pct','signed_error_pct','outage_penalty_pp']:
        baseline_checks[col]=bool(np.allclose(check[col],old[col],atol=1e-10,rtol=1e-10))
        assert baseline_checks[col],col
    new_ev=new.groupby(['event_id','pollutant','scenario','method']).agg(mean_abs_pct=('absolute_error_pct','mean'),gap_volume_share_pct=('gap_volume_share_pct','mean'),gap_mass_share_pct=('gap_mass_share_pct','mean'),gap_duration_share_pct=('gap_duration_share_pct','first')).reset_index()
    new_ev.to_csv(O/'additional_event_metrics.csv',index=False)
    summ=[]
    for key,g in new_ev.groupby(['pollutant','scenario','method']):
        summ.append([*key,len(g),*ci_mean(g.mean_abs_pct),g.gap_volume_share_pct.mean(),g.gap_mass_share_pct.mean()])
    summary=pd.DataFrame(summ,columns=['pollutant','scenario','method','n_events','mean_abs_pct','ci_low','ci_high','gap_volume_share_pct','gap_mass_share_pct'])
    summary.to_csv(O/'additional_summary.csv',index=False)
    flow_contrast=[]
    for pol in ['COD','TSS']:
        g=new_ev[(new_ev.method=='linear')&(new_ev.pollutant==pol)].pivot(index='event_id',columns='scenario',values='mean_abs_pct')
        flow_contrast.append([pol,*ci_mean(g.peakQ60-g.flowweighted60)])
    pd.DataFrame(flow_contrast,columns=['pollutant','peak_minus_flowweighted_pp','ci_low','ci_high']).to_csv(O/'peak_flowweighted_contrast.csv',index=False)
    q=base.read_series('hydraulic data/inflow-analogue.csv','flow_rate').dropna()
    a=pd.DataFrame(audit)
    a['observed_peak_Ls']=[float(q.loc[e['start']:e['end']].max()) for e in audit]
    a['season']=pd.to_datetime(a.start,utc=True).dt.month.map(lambda m: 'DJF' if m in [12,1,2] else 'MAM' if m in [3,4,5] else 'JJA' if m in [6,7,8] else 'SON')
    a.to_csv(O/'selection_comparison_events.csv',index=False)
    selection=[]
    for retained,g in a.groupby('retained'):
        selection.append({'group':'retained' if retained else 'excluded','n':len(g),'rain_median_mm':float(g.rainfall_mm.median()),'duration_median_min':float(g.duration_min.median()),'observed_peak_median_Ls':float(g.observed_peak_Ls.median()),'observed_peak_available_n':int(g.observed_peak_Ls.notna().sum()),**{s:int((g.season==s).sum()) for s in ['DJF','MAM','JJA','SON']}})
    pd.DataFrame(selection).to_csv(O/'selection_comparison_summary.csv',index=False)
    E['season']=pd.to_datetime(E.start,utc=True).dt.month.map(lambda m:'DJF' if m in [12,1,2] else 'MAM' if m in [3,4,5] else 'JJA' if m in [6,7,8] else 'SON')
    E['rain_group']=np.where(E.rainfall_mm<=E.rainfall_mm.median(),'at_or_below_median','above_median')
    E['flow_group']=np.where(E.Q_peak_Ls<=E.Q_peak_Ls.median(),'at_or_below_median','above_median')
    g=EV[(EV.interval_min==15)&(EV.method=='linear')&(EV.scenario.isin(['none','random60','peakQ60']))].merge(E[['event_id','season','rain_group','flow_group']],on='event_id')
    groups=[]
    for group in ['season','rain_group','flow_group']:
        for key,t in g.groupby([group,'pollutant','scenario']):groups.append([group,*key,len(t),t.mean_abs_pct.mean(),t.mean_abs_pct.median()])
    pd.DataFrame(groups,columns=['grouping','group','pollutant','scenario','n_events','mean_abs_pct','median_event_mean_abs_pct']).to_csv(O/'subgroup_summary.csv',index=False)
    acf=pd.DataFrame(acf,columns=['event_id','pollutant','acf_1e_minutes','max_lag_minutes','crossing_observed','concentration_cv'])
    acf.to_csv(O/'event_correlation_scales.csv',index=False)
    scaled=EV[(EV.scenario=='none')&(EV.method=='linear')].merge(acf,on=['event_id','pollutant'])
    scaled['interval_over_tau']=scaled.interval_min/scaled.acf_1e_minutes
    scaled.to_csv(O/'interval_correlation_scale.csv',index=False)
    correlations=[]
    for (pol,interval),t in scaled.groupby(['pollutant','interval_min']):
        ok=t.dropna(subset=['interval_over_tau'])
        rho=ok.acf_1e_minutes.rank().corr(ok.mean_abs_pct.rank())
        correlations.append([pol,interval,len(ok),rho])
    pd.DataFrame(correlations,columns=['pollutant','interval_min','n_events','spearman_rho']).to_csv(O/'correlation_scale_summary.csv',index=False)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(1,2,figsize=(9,3.5),layout='constrained')
    for ax,pol in zip(axs,['COD','TSS']):
        t=scaled[(scaled.pollutant==pol)&(scaled.interval_min==15)].dropna(subset=['interval_over_tau'])
        ax.scatter(t.acf_1e_minutes,t.mean_abs_pct,s=16,alpha=.6)
        ax.set(xlabel='Concentration correlation time (min)',ylabel='Event mean absolute error (%)',title=pol+' equivalent mass');ax.grid(alpha=.2)
    fig.savefig(O/'figureS1_correlation_scale.png',dpi=300);plt.close(fig)
    fig.savefig(O/'figureS1_correlation_scale.png',dpi=300);plt.close(fig)
    meta={'high_native_event_count':len(highids),'high_native_event_ids':highids.tolist(),'duration_share_60_pct':[float(100*60/E.duration_min.max()),float(100*60/E.duration_min.min())],'duration_share_120_pct':[float(100*120/E.duration_min.max()),float(100*120/E.duration_min.min())],'figure1_checks':figure_checks,'baseline_reproduced':all(baseline_checks.values()),'baseline_checks':baseline_checks,'event_bootstrap_draws':2000,'seed':SEED,'source_version':'0.1','dataset_doi':'10.5281/zenodo.15783832'}
    (O/'review_provenance.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
    print(json.dumps(meta,indent=2));print(summary.to_string(index=False));print(high_summary.to_string(index=False));print(dist.to_string(index=False))
    print(json.dumps(meta,indent=2));print(summary.to_string(index=False));print(high_summary.to_string(index=False));print(dist.to_string(index=False))
    print(json.dumps(meta,indent=2));print(summary.to_string(index=False));print(high_summary.to_string(index=False));print(dist.to_string(index=False))
    print(json.dumps(meta,indent=2));print(summary.to_string(index=False));print(high_summary.to_string(index=False));print(dist.to_string(index=False))

if __name__=='__main__':
    main()
