from pathlib import Path
import pandas as pd
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    import matplotlib;matplotlib.use('Agg')
    A=Path(__file__).resolve().parent/'supplemental_analysis'
    new=pd.read_csv(A/'additional_event_metrics.csv');fig,axs=plt.subplots(1,2,figsize=(9,3.7),layout='constrained')
    new=pd.read_csv(A/'additional_event_metrics.csv');fig,axs=plt.subplots(1,2,figsize=(9,3.7),layout='constrained')
    for ax,pol in zip(axs,['COD','TSS']):
     for sc,label,col in [('random60','Uniform random','#347da3'),('flowweighted60','Flow-weighted random','#a47731'),('peakQ60','Peak flow','#bb3d42')]:
      t=new[(new.pollutant==pol)&(new.method=='linear')&(new.scenario==sc)];ax.scatter(t.gap_mass_share_pct,t.mean_abs_pct,s=13,alpha=.6,color=col,label=label)
     ax.set(xlabel='Reference mass inside the gap (%)',ylabel='Event mean absolute error (%)',title=pol+' equivalent mass');ax.grid(alpha=.2)
    fig.legend(*axs[0].get_legend_handles_labels(),loc='outside upper center',ncol=3,frameon=False,fontsize=9)
    fig.savefig(A/'figure5_gap_exposure.png',dpi=300);plt.close(fig)
    fig.savefig(A/'figure5_gap_exposure.png',dpi=300);plt.close(fig)

if __name__=='__main__':
    main()
