from pathlib import Path
import sys,json
ROOT=Path(__file__).resolve().parent
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyse_graz import build_events,run_event,mass_kg,SEED

A=ROOT/'analysis'; F=A/'figures'; F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':300,'pdf.fonttype':42,'ps.fonttype':42})
S=pd.read_csv(A/'main_summary.csv'); E=pd.read_csv(A/'events.csv'); V=pd.read_csv(A/'main_event_metrics.csv')
colors={'linear':'#006b86','previous':'#c65d21'}
def save(fig,name):
    fig.savefig(F/(name+'.png'),bbox_inches='tight')
    fig.savefig(F/(name+'.eps'),bbox_inches='tight')
    fig.savefig(F/(name+'.pdf'),bbox_inches='tight')
    plt.close(fig)

fig,axs=plt.subplots(1,2,figsize=(10,3.8),layout='constrained')
for ax,p in zip(axs,['COD','TSS']):
    for method in colors:
        d=S[(S.pollutant==p)&(S.scenario=='none')&(S.method==method)].sort_values('interval_min')
        ax.errorbar(d.interval_min,d.mean_abs_pct,yerr=[d.mean_abs_pct-d.ci_low,d.ci_high-d.mean_abs_pct],marker='o',color=colors[method],capsize=3,label='Linear interpolation' if method=='linear' else 'Previous observation')
    ax.set(xlabel='Observation interval (min)',ylabel='Equal-event mean absolute error (%)',title=p+' equivalent mass')
    ax.set_xticks([3,6,15,30,60]);ax.grid(axis='y',alpha=.2);ax.legend(fontsize=8)
save(fig,'figure2_interval_error')

scenarios=['none','random30','random60','random120','early60','peakQ60']
labels=['No outage','Random 30 min','Random 60 min','Random 120 min','Onset 60 min','Peak-flow 60 min']
fig,axs=plt.subplots(1,2,figsize=(10,4),layout='constrained')
mx=float(S[(S.method=='linear')].mean_abs_pct.max())
for ax,p in zip(axs,['COD','TSS']):
    d=S[(S.pollutant==p)&(S.method=='linear')].pivot(index='scenario',columns='interval_min',values='mean_abs_pct').reindex(scenarios)
    im=ax.imshow(d,aspect='auto',cmap='YlOrRd',vmin=0,vmax=mx)
    ax.set_xticks(range(5),d.columns);ax.set_yticks(range(6),labels);ax.set(xlabel='Observation interval (min)',title=p+' equivalent mass')
    for i in range(6):
        for j in range(5): ax.text(j,i,f'{d.iloc[i,j]:.1f}',ha='center',va='center',color='white' if d.iloc[i,j]>mx*.6 else 'black',fontsize=9)
fig.colorbar(im,ax=axs,shrink=.85,label='Equal-event mean absolute error (%)')
save(fig,'figure3_outage_matrix')

fig,axs=plt.subplots(1,2,figsize=(10,3.8),layout='constrained')
penalty_stats=[]
for ax,p in zip(axs,['COD','TSS']):
    d=V[(V.pollutant==p)&(V.method=='linear')&(V.interval_min==15)]
    for j,sc in enumerate(['random60','early60','peakQ60']):
        vals=d[d.scenario==sc].paired_outage_penalty_pp.to_numpy()
        rng=np.random.default_rng(SEED)
        ax.scatter(j+rng.uniform(-.13,.13,len(vals)),vals,s=13,color=colors['linear'],alpha=.5)
        ax.plot([j-.22,j+.22],[np.median(vals)]*2,color='black',lw=2)
        boot=np.mean(vals[rng.integers(0,len(vals),(2000,len(vals)))],axis=1)
        penalty_stats.append([p,sc,len(vals),vals.mean(),np.median(vals),*np.quantile(boot,[.025,.975])])
    ax.axhline(0,color='gray',lw=.8);ax.set_xticks(range(3),['Random','Onset','Peak flow'])
    ax.set(title=p+' equivalent mass',xlabel='Position of a 60-minute outage',ylabel='Paired increase in absolute error (pp)')
save(fig,'figure4_paired_penalty')
pd.DataFrame(penalty_stats,columns=['pollutant','scenario','n_events','mean_penalty_pp','median_penalty_pp','ci_low','ci_high']).to_csv(A/'paired_penalty_table.csv',index=False)

events,_=build_events()
lookup={e['event_id']:(e,g,d) for e,g,d in events}
ids=[E.iloc[np.argmin(abs(E.Q_peak_Ls-E.Q_peak_Ls.median()))].event_id,E.sort_values('Q_peak_Ls').iloc[-1].event_id]
fig,axs=plt.subplots(2,2,figsize=(10,5.4),layout='constrained')
for j,eid in enumerate(ids):
    e,g,d=lookup[eid];x=np.arange(len(g))/60-1
    axs[0,j].plot(x,d['Q'],color='#404040',lw=1);axs[0,j].set(title=eid[:8]+' UTC',ylabel='Inlet flow (L/s)')
    axs[1,j].plot(x,d['COD'],label='COD equivalent',color='#006b86',lw=1)
    axs[1,j].plot(x,d['TSS'],label='TSS equivalent',color='#c65d21',lw=1)
    for ax in axs[:,j]:
        ax.set_xlim(0,e['duration_min']/60);ax.grid(alpha=.15)
    axs[1,j].set(xlabel='Hours after rainfall onset',ylabel='Concentration (mg/L)');axs[1,j].legend(fontsize=8)
save(fig,'figure1_event_examples')

# All retained-event traces for internal visual review, six events per sheet.
for batch in range(0,len(events),6):
    chunk=events[batch:batch+6]
    fig,axs=plt.subplots(len(chunk),2,figsize=(10,2.0*len(chunk)),squeeze=False,layout='constrained')
    for j,(e,g,d) in enumerate(chunk):
        x=np.arange(len(g))/60-1
        axs[j,0].plot(x,d['Q'],lw=.65,color='black');axs[j,0].set_title(e['event_id']+' | Q',fontsize=9)
        axs[j,1].plot(x,d['COD'],lw=.65,label='COD');axs[j,1].plot(x,d['TSS'],lw=.65,label='TSS')
        for ax in axs[j]:ax.set_xlim(0,e['duration_min']/60)
        axs[j,1].legend(fontsize=7,loc='upper right')
    fig.savefig(ROOT/f'qa_events_{batch//6+1}.png',dpi=130);plt.close(fig)

audit=pd.read_csv(A/'event_audit_gap6_sep6_tail3_sat2400.csv')
stats={'events':len(E),'candidates':len(audit),'saturation_flag_candidates':int(audit.flow_saturation_flag.sum()),
       'years':E.year.value_counts().sort_index().to_dict(),'start':E.start.min(),'end':E.end.max(),
       'ranges':{c:{'min':float(E[c].min()),'median':float(E[c].median()),'max':float(E[c].max())} for c in ['rainfall_mm','duration_min','Q_peak_Ls','volume_m3','COD_mass_kg','TSS_mass_kg','COD_native_fraction','TSS_native_fraction']}}
(A/'descriptive_stats.json').write_text(json.dumps(stats,indent=2),encoding='utf-8')
print(json.dumps(stats,indent=2))
