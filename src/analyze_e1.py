"""Checkpoint 1: aggregate E1, apply the prespecified b rule, draw Figure 2."""
import csv,json,math,shutil
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_e1 import read_json,read_rows,retry_io,json_save,atomic_replace
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/E1'
COLORS={'W':'#2878A5','E':'#D97932','PP':'#8560A9'}

def read_table(path):
    rows=read_rows(path)
    return [{k:v if k in ['kernel','update','seed'] else float(v) for k,v in row.items()} for row in rows]

def save_table(path,rows):
    fields=list(dict.fromkeys(k for row in rows for k in row))
    tmp=path.with_suffix('.tmp')
    def write():
        with tmp.open('w',newline='',encoding='utf-8') as f:
            writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    retry_io(write);atomic_replace(tmp,path)

def select_b(table,cfg):
    ref=sorted([r for r in table if r['kernel']=='W' and r['M']==9],key=lambda r:r['b'])
    selected={};explanation=[]
    for label in ['A','B']:
        rule=cfg['b_selection'][label]
        eligible=[r for r in ref if r['f_C']>=rule['min_fc'] and r['all_D']<=rule['max_all_D']]
        if not eligible:return dict(state='author_decision_required',reason=f'No W(9) point meets the b_{label} thresholds')
        row=eligible[-1];selected[f'b_{label}']=float(row['b'])
        explanation.append(dict(label=label,b=row['b'],f_C=row['f_C'],all_D=row['all_D'],rule=rule))
    if selected['b_B']-selected['b_A']<cfg['b_selection']['minimum_separation']-1e-12:
        selected['b_A']=round(selected['b_A']-.05,3)
        explanation.append(dict(adjustment='b_A reduced by .05 by the prespecified separation rule',b_A=selected['b_A']))
    if selected['b_A'] not in cfg['b']:
        return dict(state='author_decision_required',reason='Separation rule selects b_A outside the E1 grid',**selected,justification=explanation)
    return dict(state='selected',**selected,justification=explanation)

def draw(table,cfg,selection,destination):
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
    fig=plt.figure(figsize=(11,8));gs=fig.add_gridspec(2,2,height_ratios=[1,1.05],hspace=.40,wspace=.24)
    for index,kernel in enumerate(['W','E']):
        ax=fig.add_subplot(gs[0,index]);lookup={(r['M'],r['b']):r['f_C'] for r in table if r['kernel']==kernel}
        data=np.array([[lookup[M,b] for b in cfg['b']] for M in cfg['M']])
        image=ax.imshow(data,origin='lower',aspect='auto',vmin=0,vmax=1,cmap='viridis',extent=[-.5,len(cfg['b'])-.5,-.5,len(cfg['M'])-.5])
        ticks=list(range(0,len(cfg['b']),4));ax.set_xticks(ticks,[f'{cfg["b"][i]:.2f}' for i in ticks])
        ax.set_yticks(range(len(cfg['M'])),cfg['M']);ax.set_xlabel('Temptation b');ax.set_ylabel('Memory length M')
        ax.set_title(f'({chr(97+index)}) {"Window" if kernel=="W" else "Exponential"}')
        if selection['state']=='selected':
            for label,style in [('A','--'),('B',':')]:ax.axvline(cfg['b'].index(selection[f'b_{label}']),color='white',ls=style,lw=1)
        fig.colorbar(image,ax=ax,label='Mean cooperation $f_C$',shrink=.85)
    ax=fig.add_subplot(gs[1,:])
    if selection['state']=='selected':
        for kernel in ['W','E']:
            for label,style in [('A','-'),('B','--')]:
                b=selection[f'b_{label}'];d=sorted([r for r in table if r['kernel']==kernel and abs(r['b']-b)<1e-10],key=lambda r:r['M'])
                x=[r['M'] for r in d];y=[r['f_C'] for r in d]
                ax.errorbar(x,y,yerr=[r['f_C_se'] for r in d],color=COLORS[kernel],ls=style,lw=1.5,capsize=2,label=f'{kernel}, $b_{label}$ = {b:g}')
                ax.scatter(x,y,s=[14+90*r['all_D'] for r in d],color=COLORS[kernel],alpha=.7,zorder=3)
        ax.set_xscale('log',base=2);ax.set_xticks(cfg['M'],[str(x) for x in cfg['M']]);ax.set_ylim(-.03,1.03)
        ax.legend(ncol=2,loc='best');ax.set_ylabel('Mean cooperation $f_C$');ax.set_xlabel('Memory length M')
        ax.set_title('(c) Selected temptations; marker area increases with all-D fraction')
    else:
        ax.axis('off');ax.text(.5,.5,selection['reason'],ha='center',va='center',wrap=True,transform=ax.transAxes)
    for ext in ['png','pdf','svg']:retry_io(lambda:fig.savefig(destination/f'fig02_phase_map.{ext}',bbox_inches='tight'))
    plt.close(fig)

def analyze():
    cfg=read_json(ROOT/'config/E1.json');manifest=read_json(OUT/'manifest.json')
    if manifest['configuration']!=cfg:raise RuntimeError('Analysis configuration differs from run manifest')
    rows=[];tables={};metadata={}
    for folder in sorted((OUT/'points').iterdir()):
        if not folder.is_dir():continue
        m=read_json(folder/'point.json')
        if m['state']!='complete':raise RuntimeError(f'Incomplete point {folder.name}')
        d=read_table(folder/m['selected_file'])
        if len(d)!=cfg['R']:raise RuntimeError('Wrong realization count')
        first=d[0];r={k:first[k] for k in ['kernel','L','M','b','epsilon','w','K','eta','T0']}
        r.update(R=len(d),point=folder.name,source=str((folder/m['selected_file']).relative_to(ROOT)),
            repeated=m['repeated'],still_nonstationary=m['still_nonstationary'])
        for key in ['f_C','all_D','all_C','G','G_time','G_approx','G_from_Cy','relative_gap_discrepancy']:
            v=np.array([row[key] for row in d if math.isfinite(row[key])]);r[key]=float(v.mean()) if len(v) else float('nan');r[f'{key}_se']=float(v.std(ddof=1)/np.sqrt(len(v))) if len(v)>1 else float('nan');r[f'{key}_n']=len(v)
        for key in [k for k in first if k.startswith('C_s_') or k.startswith('C_y_')]:
            v=np.array([row[key] for row in d if math.isfinite(row[key])]);r[key]=float(v.mean()) if len(v) else float('nan');r[f'{key}_se']=float(v.std(ddof=1)/np.sqrt(len(v))) if len(v)>1 else float('nan')
        r['gap_excluded_realizations']=sum(not math.isfinite(row['G']) for row in d)
        r['identity_max_absolute_error']=max((row['identity_absolute_error'] for row in d if math.isfinite(row['identity_absolute_error'])),default=float('nan'))
        rows.append(r);tables[folder.name]=d;metadata[folder.name]=m
    expected=len(cfg['kernels'])*len(cfg['M'])*len(cfg['b'])
    if len(rows)!=expected:raise RuntimeError(f'Expected {expected} points, got {len(rows)}')
    table=rows;save_table(OUT/'point_summary.csv',table)
    selection=select_b(table,cfg);json_save(OUT/'b_selection.json',selection)
    figures=ROOT/'results/figures';figures.mkdir(parents=True,exist_ok=True);draw(table,cfg,selection,figures)
    retained=[]
    if selection['state']=='selected':
        dest=ROOT/'results/E2';dest.mkdir(exist_ok=True)
        for label in ['A','B']:
            for M in [1,9,33]:
                b=selection[f'b_{label}'];name=f'W_M{M:02d}_b{b:.3f}';d=tables[name];meta=metadata[name]
                representative=int(sorted(d,key=lambda r:r['f_C'])[len(d)//2]['realization'])
                for rep in range(cfg['R']):
                    filename=f'{name}_T0{meta["T0"]}_r{rep:02d}.npy';source=OUT/'E2_candidates'/filename;target=dest/filename
                    if source.exists():retry_io(lambda:shutil.copyfile(source,target))
                    elif not target.exists():raise FileNotFoundError(source)
                retained.append(dict(b_label=label,b=b,M=M,T0=meta['T0'],representative_realization=representative,
                    selection_rule='Middle realization in stable sort by f_C; final state shown, no visual selection',
                    file=f'{name}_T0{meta["T0"]}_r{representative:02d}.npy'))
        json_save(dest/'selection.json',retained)
        # Only direct candidate files under this verified project folder are removed.
        candidate_root=(OUT/'E2_candidates').resolve()
        if not candidate_root.is_relative_to(ROOT.resolve()):raise RuntimeError('Candidate folder outside workspace')
        for path in candidate_root.glob('*.npy'):
            if path.resolve().parent!=candidate_root:raise RuntimeError('Unexpected candidate path')
            retry_io(path.unlink)
    budget=read_json(ROOT/'qa/runtime_budget.json');est=budget['estimates'][str(cfg['R'])]
    report=['# Ponto de checagem 1','',f'E1 completo: {expected} pontos, {cfg["R"]} realizações por ponto.',
        f'Pontos repetidos com T0 = 10000: {sum(r["repeated"] for r in table)}; ainda com critério de não estacionariedade: {sum(r["still_nonstationary"] for r in table)}.','',
        'Mapas e cortes: `results/figures/fig02_phase_map.png` (também PDF e SVG).',
        'Médias, erros-padrão, absorções e contagens de G: `results/E1/point_summary.csv`.','',
        '## Escolha de b_A e b_B','',json.dumps(selection,indent=2,ensure_ascii=False),'',
        '## Custo estimado das etapas seguintes','',
        f'Workers: {budget["workers"]}. Estimativa de CPU em horas, sem repetições de estacionariedade:']
    report += [f'- {k}: {v:.2f} h' for k,v in est['worker_hours'].items() if k!='E1']
    report += ['', 'Estes valores são projeções; cada experimento será cronometrado novamente antes da execução.',
        '', '**Parada obrigatória:** os autores precisam aprovar este ponto antes de E3, E5, E6, S1 e S2. Nenhuma etapa seguinte foi iniciada.']
    retry_io(lambda:(ROOT/'checkpoint1.md').write_text('\n'.join(report)+'\n',encoding='utf-8'))
    return selection
if __name__=='__main__':analyze()
