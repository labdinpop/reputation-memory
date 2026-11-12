"""Aggregate the approved S3 and draw its single supplementary figure."""
from pathlib import Path
import json
import numpy as np
from run_e1 import read_rows,read_json
from analyze_e1 import save_table
from analyze_phase2 import load_summary,subset,mean_se
from restyle_figures import plt,line,save,COLORS
ROOT=Path(__file__).resolve().parents[1]
def main():
    rows=[]
    for p in sorted((ROOT/'results/S3/points').glob('*/point.json')):
        meta=read_json(p);assert meta['state']=='complete'
        data=read_rows(p.parent/meta['selected_file']);assert len(data)==30
        r={k:float(data[0][k]) for k in ['L','M','b','epsilon','w','K','eta','T0']}
        r.update(kernel='W',R=30,source=meta['source'],still_nonstationary=meta['still_nonstationary'])
        for k in ['f_C','G','all_D','all_C','relative_gap_discrepancy']:
            vals=[float(d[k]) for d in data];v,se,n=mean_se(vals);r.update({k:v,k+'_se':se,k+'_n':n})
        rows.append(r)
    assert len(rows)==7
    save_table(ROOT/'results/S3/point_summary.csv',rows)
    baseline=[r for r in load_summary(ROOT/'results/E1/point_summary.csv') if r['kernel']=='W' and r['b']==1.15 and r['M'] in [1,2,3,5,9,17,33]]
    fig,ax=plt.subplots(figsize=(6.35,3.2),layout='constrained')
    line(ax,baseline,'M',COLORS['W'],'o',label='Mixed selection (Eq. 5)')
    line(ax,rows,'M',COLORS['W'],'s',label='Proportional selection',ls='--',open_=True)
    ax.set_xscale('log',base=2);ax.set_xticks([1,2,3,5,9,17,33],[1,2,3,5,9,17,33])
    ax.set(xlabel='Memory length $M$',ylabel='Mean cooperation $f_C$',ylim=(-.03,1.05))
    ax.legend(loc='lower left',frameon=False);save(fig,'figS3_li_selection')
    print(json.dumps([{k:r[k] for k in ['M','f_C','f_C_se','all_C','all_D','T0','still_nonstationary']} for r in rows],indent=2))
if __name__=='__main__':main()
