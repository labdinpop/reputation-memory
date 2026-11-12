"""Compare zero floor (saved S3) and two approved positive floors."""
from pathlib import Path
import json,numpy as np
from run_e1 import read_rows,read_json
from analyze_e1 import save_table
from analyze_phase2 import mean_se
from restyle_figures import plt,line,save,ordered,AUDIT
ROOT=Path(__file__).resolve().parents[1]
def main():
    rows=[]
    for name in ['S3','S3_floors']:
        for path in sorted((ROOT/f'results/{name}/points').glob('*/point.json')):
            meta=read_json(path);assert meta['state']=='complete'
            data=read_rows(path.parent/meta['selected_file']);assert len(data)==30
            r={k:float(data[0][k]) for k in ['L','M','b','epsilon','w','K','eta','T0']}
            r.update(kernel='W',R=30,source=meta['source'],still_nonstationary=meta['still_nonstationary'])
            for k in ['f_C','G','all_D','all_C','relative_gap_discrepancy']:
                v,se,n=mean_se([float(d[k]) for d in data]);r.update({k:v,k+'_se':se,k+'_n':n})
            rows.append(r)
    assert len(rows)==21
    save_table(ROOT/'results/S3_floors/point_summary.csv',rows)
    fig,ax=plt.subplots(figsize=(6.35,3.8),layout='constrained')
    labels=[r'$\eta=0$',r'$\eta=10^{-6}$ (Li et al. [17])',r'$\eta=0.01$']
    for i,eta in enumerate([0,1e-6,.01]):
        data=[r for r in rows if r['eta']==eta]
        line(ax,data,'M',ordered(3)[i],['o','s','^'][i],label=labels[i],ls=['-','--',':'][i],open_=i==1,offset=[-.035,0,.035][i])
    ax.set_xscale('log',base=2);ax.set_xticks([1,2,3,5,9,17,33],[1,2,3,5,9,17,33])
    ax.set(xlabel='Memory length $M$',ylabel='Mean cooperation $f_C$',ylim=(-.03,1.055))
    ax.legend(loc='lower left',frameon=False);save(fig,'figS3_li_selection')
    (ROOT/'qa/S3_floors/figure_style.json').write_text(json.dumps(AUDIT,indent=2))
    print(json.dumps([{k:r[k] for k in ['eta','M','f_C','f_C_se','all_C','T0','still_nonstationary']} for r in rows],indent=2))
if __name__=='__main__':main()
