"""Redraw the manuscript's aggregate evidence without participant-level data."""
from pathlib import Path
import shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch

ROOT=Path(__file__).resolve().parents[2]
COLORS={'direct':'#A6D0C6','within_cf_ridge':'#DADDE3','cross_batch_ridge':'#5878A6','CRNP_v0.2.0':'#F4A7AE'}
LABELS={'direct':'Direct','within_cf_ridge':'Within-CF diagnostic','cross_batch_ridge':'Source residual ridge','CRNP_v0.2.0':'CRNP v0.2.0'}

def _panel(ax,label,title):
    ax.set_title(label+'  '+title,loc='left',fontweight='bold',pad=10)
    ax.spines[['top','right']].set_visible(False)
    ax.tick_params(direction='out',length=3)

def _save(fig,out,name):
    fig.savefig(out/(name+'.pdf'),bbox_inches='tight')
    fig.savefig(out/(name+'.png'),dpi=300,bbox_inches='tight');plt.close(fig)

def redraw(output):
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=True)
    tables=ROOT/'results/published'
    plt.rcParams.update({'font.family':'DejaVu Serif','font.size':8,'axes.titlesize':9,'axes.labelsize':8,'pdf.fonttype':42,'ps.fonttype':42})
    fig,ax=plt.subplots(figsize=(7.1,3.5));ax.set_axis_off();ax.set_xlim(0,1);ax.set_ylim(0,1)
    boxes=[(.02,.64,.28,.27,'Frozen representation\n1,222 proteins; 9 modules\nHMACS batches: 722 / 833'),(.36,.64,.28,.27,'Structural coverage masks\nMGH / Sissala protein sets\n9 modules x 4 masks x 2 directions'),(.7,.64,.28,.27,'Score reconstruction\nDirect / source ridge / CRNP\nWithin-CF target diagnostic'),(.02,.14,.28,.29,'Recovery endpoint\nMSE / sample variance of full score\nEqual settings, then equal directions\n2,000 paired participant draws'),(.36,.14,.28,.29,'Computational cost\n24 workloads; 5 repeats\nSingle numeric thread\nInterleaved Python / R methods'),(.7,.14,.28,.29,'Paired-assay evaluation\n35 Emory participants; 6 branches\n52 evaluable concordance cells\n10,000 paired participant draws')]
    for i,(x,y,w,h,text) in enumerate(boxes):
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.008',facecolor=['#D9EEE7','#C7DDEE','#FBD9C9'][i%3],edgecolor='#5878A6',linewidth=.8));ax.text(x+w/2,y+h/2,text,ha='center',va='center',linespacing=1.6)
    for a,b in [((.3,.775),(.36,.775)),((.64,.775),(.7,.775)),((.84,.64),(.84,.43)),((.5,.64),(.5,.43)),((.16,.64),(.16,.43))]:ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=11,color='#5878A6'))
    ax.text(.5,.03,'Frozen modules use joint upstream preprocessing; recovery intervals condition on fixed fits and predictions.',ha='center',fontsize=7)
    _save(fig,out,'Figure1_redrawn')
    overall=pd.read_csv(tables/'HMACS_overall.csv').set_index('method')
    settings=pd.read_csv(tables/'HMACS_72_settings.csv')
    fig,axes=plt.subplots(2,2,figsize=(7.1,5.1),layout='constrained')
    methods=['direct','cross_batch_ridge','CRNP_v0.2.0','within_cf_ridge']
    for ax,metric,title,factor in [(axes[0,0],'risk','Normalized reconstruction error',1),(axes[0,1],'relative_gain','Relative gain vs direct (%)',100)]:
        use=methods if metric=='risk' else methods[1:3]
        for i,m in enumerate(use):
            d=overall.loc[m];v=d[metric]*factor;lo=d[metric+'_ci_low']*factor;hi=d[metric+'_ci_high']*factor
            ax.barh(i,v,color=COLORS[m],edgecolor='#34435D',linewidth=.4,height=.58)
            ax.errorbar(v,i,xerr=[[v-lo],[hi-v]],color='#34435D',capsize=2,linewidth=.8,fmt='none')
            ax.text(v,i,f' {v:.4f}' if factor==1 else f' {v:.2f}%',va='center',fontsize=7)
        ax.set_yticks(range(len(use)),[LABELS[m] for m in use]);ax.invert_yaxis();ax.set_xlabel('NMSE' if factor==1 else 'Percent');ax.margins(x=.24)
        _panel(ax,'a' if metric=='risk' else 'b',title)
    for ax,m,letter in [(axes[1,0],'cross_batch_ridge','c'),(axes[1,1],'CRNP_v0.2.0','d')]:
        settings['column']=settings.platform_label+' / '+settings.direction
        v=settings.pivot(index='module',columns='column',values=m+'__relative_gain')*100
        norm=TwoSlopeNorm(vmin=min(-1,float(v.min().min())),vcenter=0,vmax=max(1,float(v.max().max())))
        im=ax.imshow(v.to_numpy(),aspect='auto',cmap='RdBu',norm=norm)
        ax.set_yticks(range(len(v.index)),v.index);ax.set_xticks(range(len(v.columns)),[c.replace('_',' ') for c in v.columns],rotation=70,ha='right',fontsize=6)
        _panel(ax,letter,LABELS[m]+': cellwise relative gain')
        fig.colorbar(im,ax=ax,shrink=.75,label='Gain (%)')
    _save(fig,out,'Figure2_redrawn')
    timing=pd.read_csv(tables/'Timing_overall.csv').set_index('method');emory=pd.read_csv(tables/'Emory_concordance.csv')
    fig,axes=plt.subplots(1,3,figsize=(7.1,2.8),layout='constrained')
    for ax,metric,letter,title in [(axes[0],'fit_s','a','Fit time (s)'),(axes[1],'prediction_ms','b','Whole-batch score time (ms)')]:
        for i,m in enumerate(['direct','cross_batch_ridge','CRNP_v0.2.0']):
            d=timing.loc[m];v=d[metric+'_median'];lo=d[metric+'_median_across_case_q25'];hi=d[metric+'_median_across_case_q75']
            ax.bar(i,v,color=COLORS[m],edgecolor='#34435D',linewidth=.5)
            ax.errorbar(i,v,yerr=[[v-lo],[hi-v]],fmt='none',color='#34435D',capsize=2,linewidth=.8)
        ax.set_xticks(range(3),['Direct','Source ridge','CRNP'],rotation=35,ha='right');ax.set_ylabel(title);_panel(ax,letter,title)
    ax=axes[2]
    for i,r in enumerate(emory.itertuples()):
        ax.errorbar(r.estimate,i,xerr=[[r.estimate-r.ci_low],[r.ci_high-r.estimate]],fmt='o',color='#5878A6',capsize=3,markersize=4)
    ax.axvline(0,color='#8D96A4',linestyle='--',linewidth=.8);ax.set_yticks(range(3),['All (52 cells)','CSF (27)','Plasma (25)']);ax.invert_yaxis();ax.set_xlabel('Mean Fisher-z difference');_panel(ax,'c','Paired-assay concordance')
    fig.text(.02,-.12,'Timing error bars: IQR across workload medians. Concordance error bars: paired participant bootstrap 95% CI.\nConcordance intervals include zero; superiority and equivalence were not established.',fontsize=7)
    _save(fig,out,'Figure3_redrawn')
    (out/'README.txt').write_text('Redrawn from frozen public aggregate tables. Numerical evidence is preserved; panel layout differs from the submitted PDF snapshots in paper_figures/. No participant matrices or model fits are used. Timing is the original machine-specific summary, not a fresh timing experiment.\n',encoding='utf-8')
    return [str(out/(n+'_redrawn.pdf')) for n in ['Figure1','Figure2','Figure3']]
