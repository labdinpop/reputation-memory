"""Read completed simulations; generate added numbers, drift list, and Figure S3.

No simulation is launched and no previous result file is modified.
"""
from pathlib import Path
import csv, hashlib, json, math, re, statistics
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter

WORK=Path(__file__).resolve().parent
OUT=WORK.parent
ROOT=OUT.parents[1]
# The project edits artigo first; an extracted data release keeps its portable layout.
SRC=ROOT/'artigo' if (ROOT/'artigo').is_dir() else OUT/'sources'
DATA=OUT/'revision_data'
DATA.mkdir(exist_ok=True)

def read(path):
    with path.open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f))

def path_for(value):
    return ROOT/Path(value.replace('\\','/'))

def estimate(rows, field):
    x=[float(r[field]) for r in rows if math.isfinite(float(r[field]))]
    return statistics.mean(x), statistics.stdev(x)/math.sqrt(len(x)) if len(x)>1 else 0

values={}; provenance={}
def put(key,value,source):
    values[key]=str(value);provenance[key]=source

past_cfg=json.loads((ROOT/'config/S3_past.json').read_text())
past=[]
for m in past_cfg['m']:
    folder=ROOT/'results/S3_past/points'/f'PP_a0_m{m:02d}'
    meta=json.loads((folder/'point.json').read_text())
    path=folder/meta['selected_file']
    rows=read(path)
    assert len(rows)==past_cfg['R']==len({r['seed'] for r in rows})
    assert all(float(r['a'])==0 and float(r['w'])==1 and float(r['b'])==past_cfg['b'] for r in rows)
    avg,se=estimate(rows,'f_C')
    rec=dict(m=m,f_C=avg,f_C_se=se,R=len(rows),all_C=sum(int(r['all_C']) for r in rows),all_D=sum(int(r['all_D']) for r in rows),source=str(path.relative_to(ROOT)).replace('\\','/'))
    past.append(rec)
    put(f'r2Past{m}fc',f'{avg:.3f}',dict(file=rec['source'],field='f_C',operation='mean',scale=1))
    put(f'r2Past{m}fcSE',f'{se:.6f}',dict(file=rec['source'],field='f_C',operation='standard_error',scale=1))

put('r2PastW',f"{past_cfg['w']:g}",'config/S3_past.json:w')
full=[r['m'] for r in past if r['all_C']==r['R']]
put('r2PastFullGrid',', '.join(map(str,full)),dict(files=[r['source'] for r in past],operation='all_cooperator_points'))
put('r2PastGrid',', '.join(map(str,past_cfg['m'])),'config/S3_past.json:m')
put('r2PastLongOne',past_cfg['m'][-2],'config/S3_past.json:m[-2]')
put('r2PastLongTwo',past_cfg['m'][-1],'config/S3_past.json:m[-1]')
put('r2PastTotalRuns',sum(r['R'] for r in past),dict(files=[r['source'] for r in past],operation='row_count'))
paper=(SRC/'paper.tex').read_text(encoding='utf-8')
citeorder=list(dict.fromkeys(k for group in re.findall(r'\\cite\{([^}]+)\}',paper) for k in group.split(',')))
put('r2LiReferenceNumber',citeorder.index('li2026')+1,'citation order in sources/paper.tex')

flagged={}
for exp in ['E1','E3','E5','E6','S1','S2']:
    for row in read(ROOT/'results'/exp/'point_summary.csv'):
        if row['still_nonstationary'].lower()!='true':
            continue
        source=row['source'].replace('\\','/')
        if source in flagged:
            flagged[source]['experiments']+=','+exp
            continue
        raw=read(path_for(source)); p=raw[0]
        count=sum(abs(float(r['half_difference']))>.02 for r in raw)
        fraction=count/len(raw)
        assert fraction>.1
        rec={k:p.get(k,'') for k in ['kernel','L','b','M','a','m','epsilon','w','K','eta','T0','T1']}
        rec.update(R=len(raw),drifting_realizations=count,drifting_fraction=fraction,experiments=exp,source=source)
        flagged[source]=rec
drift=list(flagged.values())
assert len(drift)==14
assert sum(int(r['T0'])==20000 for r in drift)==5
with (DATA/'residual_drift.csv').open('w',encoding='utf-8',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=drift[0].keys());writer.writeheader();writer.writerows(drift)
put('r2DriftAt10000',sum(int(r['T0'])==10000 for r in drift),dict(file='revision_data/residual_drift.csv',operation='count_T0_10000'))

base=(WORK/'numbers.tex.before').read_text(encoding='utf-8-sig')
extra='\n% Added values generated from completed simulations by build_revision_assets.py.\n'
extra+='\n'.join(r'\expandafter\def\csname value@'+k+r'\endcsname{'+v+'}' for k,v in values.items())+'\n'
(SRC/'numbers.tex').write_text(base+extra,encoding='utf-8')
(WORK/'revision_numbers_sources.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
(DATA/'past_memory_summary.json').write_text(json.dumps(past,indent=2),encoding='utf-8')

# The table is included only in the reply, not in the manuscript or supplement.
table=[r'\begin{center}',r'\small',r'\begin{tabular}{lrrrrr}',r'\hline',
       r'Memory & $b$ & $\varepsilon$ & $L$ & $T_0$ & Runs with drift \\',r'\hline']
for row in drift:
    if row['kernel']=='PP':
        kernel=r'$\mathrm{PP}('+f"{float(row['a']):g},{int(row['m'])}"+')$'
    else:
        kernel='$'+row['kernel']+'('+str(int(row['M']))+')$'
    table.append(kernel+' & '+f"{float(row['b']):g} & {float(row['epsilon']):g} & {row['L']} & {int(row['T0']):,} & {row['drifting_realizations']}/{row['R']}"+r' \\')
table += [r'\hline',r'\end{tabular}',r'\end{center}']
(OUT/'reply/drift_table.tex').write_text('\n'.join(table)+'\n',encoding='utf-8')

windows=[]
for row in read(ROOT/'results/S3_floors/point_summary.csv'):
    raw=read(path_for(row['source'])); avg,se=estimate(raw,'f_C')
    assert abs(avg-float(row['f_C']))<1e-12
    windows.append(dict(M=int(float(row['M'])),eta=float(row['eta']),f_C=avg,f_C_se=se))
plt.rcParams.update({'font.size':13,'axes.titlesize':13,'legend.fontsize':11,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42})
fig,axes=plt.subplots(1,2,figsize=(10.3,4.5),layout='constrained',sharey=True)
for eta,color,marker,linestyle,offset in [(0,'#0072B2','o','-',-.018),(1e-6,'#0072B2','s','--',0),(.01,'#56B4E9','D',':',.018)]:
    curve=sorted((r for r in windows if r['eta']==eta),key=lambda r:r['M'])
    label=r'$\eta=0$' if eta==0 else (r'$\eta=10^{-6}$' if eta==1e-6 else r'$\eta=0.01$')
    axes[0].errorbar([r['M']*(2**offset) for r in curve],[r['f_C'] for r in curve],yerr=[r['f_C_se'] for r in curve],color=color,marker=marker,ms=6,mfc='white',ls=linestyle,capsize=3,label=label)
basecurve=sorted((r for r in windows if r['eta']==1e-6),key=lambda r:r['M'])
axes[1].errorbar([r['M'] for r in basecurve],[r['f_C'] for r in basecurve],yerr=[r['f_C_se'] for r in basecurve],color='#0072B2',marker='o',mfc='white',ls='--',capsize=3,label=r'W($M$): current action included')
axes[1].errorbar([r['m'] for r in past],[r['f_C'] for r in past],yerr=[r['f_C_se'] for r in past],color='#009E73',marker='^',ls='-',capsize=3,label=r'PP(0,$m$): previous actions')
for ax in axes:
    ax.set_xscale('log',base=2);ax.set_ylim(-.025,1.06);ax.set_xlim(.85,39)
    ax.xaxis.set_major_locator(FixedLocator([1,2,5,10,20,33]));ax.xaxis.set_major_formatter(FuncFormatter(lambda x,pos:f'{x:g}'))
    ax.legend(loc='lower left',frameon=False)
axes[0].set(title='(a) Window memory',xlabel=r'Memory length $M$ (rounds)',ylabel=r'Mean cooperation $f_C$')
axes[1].set(title='(b) Including or excluding the current action',xlabel=r'Number of records ($M$ or $m$)')
for ext in ['pdf','png']:
    fig.savefig(SRC/'figures'/f'figS3_li_selection.{ext}',dpi=200)
plt.close(fig)
report=dict(past_realizations=sum(r['R'] for r in past),residual_drift_combinations=len(drift),after_20000=sum(int(r['T0'])==20000 for r in drift),li_reference_number=values['r2LiReferenceNumber'],new_values=values)
(WORK/'asset_report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print('Generated Figure S3 from saved data, 14 residual-drift rows, and',len(values),'traced values.')
