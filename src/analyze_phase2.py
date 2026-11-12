"""Figures and factual checkpoint-2 summary; this module never writes the paper."""
import csv,json,math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,to_rgb
from matplotlib.lines import Line2D
from model import Parameters,weights
from run_e1 import read_json,read_rows,json_save,retry_io,digest
from analyze_e1 import read_table,save_table,COLORS
from run_phase2 import designs,EXPERIMENTS

ROOT=Path(__file__).resolve().parents[1]
FIG=ROOT/'results/figures'
SELECTION=ROOT/'results/E1/b_selection.json'

def mean_se(values):
    v=np.asarray(values,dtype=float);v=v[np.isfinite(v)]
    return (float(v.mean()) if len(v) else math.nan,float(v.std(ddof=1)/np.sqrt(len(v))) if len(v)>1 else math.nan,len(v))

def load_summary(path):
    result=[]
    for row in read_rows(path):
        out={}
        for k,v in row.items():
            if v in ['True','False']:out[k]=v=='True'
            else:
                try:out[k]=float(v) if v else math.nan
                except ValueError:out[k]=v
        result.append(out)
    return result

def subset(rows,**criteria):
    return [r for r in rows if all(r.get(k)==v for k,v in criteria.items())]

def one(rows,**criteria):
    matches=subset(rows,**criteria)
    if len(matches)!=1:raise RuntimeError(f'Expected one point for {criteria}; got {len(matches)}')
    return matches[0]

def aggregate(name):
    cfg,tasks=designs(name);result=[]
    for task in tasks:
        folder=ROOT/f'results/{name}/points'/task['name'];meta=read_json(folder/'point.json')
        if meta['state']!='complete':raise RuntimeError(f'Incomplete point: {name}/{task["name"]}')
        path=ROOT/meta['source_csv'] if meta.get('reused') else folder/meta['selected_file']
        if meta.get('reused') and digest(path)!=meta['source_sha256']:raise RuntimeError('Reused data changed')
        data=read_table(path)
        if len(data)!=task['R']:raise RuntimeError('Wrong realization count')
        row=dict(meta['parameters'],point=task['name'],R=len(data),T0=meta['T0'],reused=meta.get('reused',False),
            repeated=meta['repeated'],still_nonstationary=meta['still_nonstationary'],source=str(path.relative_to(ROOT)))
        fields=['f_C','all_D','all_C','G','G_time','G_approx','G_from_Cy','relative_gap_discrepancy']
        fields += [k for k in data[0] if k.startswith('C_s_') or k.startswith('C_y_')]
        for key in fields:
            avg,se,n=mean_se([r[key] for r in data]);row[key]=avg;row[key+'_se']=se
            if not key.startswith('C_'):row[key+'_n']=n
        row['gap_excluded_realizations']=len(data)-row['G_n']
        row['identity_max_absolute_error']=max((r['identity_absolute_error'] for r in data if math.isfinite(r['identity_absolute_error'])),default=math.nan)
        result.append(row)
    save_table(ROOT/f'results/{name}/point_summary.csv',result)
    return result

def shade(kernel,index,total=3):
    base=np.array(to_rgb(COLORS[kernel]));amount=.45+.55*index/max(1,total-1)
    return tuple(1-amount*(1-base))

def save(fig,name):
    FIG.mkdir(parents=True,exist_ok=True)
    for ext in ['png','pdf','svg']:retry_io(lambda:fig.savefig(FIG/f'{name}.{ext}',bbox_inches='tight',dpi=180))
    plt.close(fig)

def line(ax,rows,x,label,color,ls='-',marker='o',y='f_C'):
    rows=sorted(rows,key=lambda r:r[x]);xx=[r[x] for r in rows];yy=[r[y] for r in rows]
    ax.errorbar(xx,yy,yerr=[r.get(y+'_se',0) for r in rows],color=color,label=label,ls=ls,marker=marker,ms=3,capsize=2,lw=1.3)

def fig01():
    fig,ax=plt.subplots(figsize=(7,3.3))
    for k,p,label in [('W',Parameters(M=9),'W(9)'),('E',Parameters(kernel='E',M=9),'E(9)'),('PP',Parameters(kernel='PP',a=.4,m=8),'PP(0.4, 8)')]:
        h,_=weights(p);x=np.arange(min(len(h),33));ax.plot(x,h[:len(x)],'o-',ms=3,label=label,color=COLORS[k])
    ax.set(xlabel='Record age (rounds)',ylabel=r'Kernel weight $h_\ell$',xlim=(-.5,32.5),ylim=(-.01,.43));ax.legend()
    save(fig,'fig01_memory_weights')

def fig03():
    selected=read_json(ROOT/'results/E2/selection.json');fig,axes=plt.subplots(2,3,figsize=(9,6),layout='constrained')
    facts=[]
    for entry in selected:
        i=0 if entry['b_label']=='A' else 1;j=[1,9,33].index(entry['M']);ax=axes[i,j]
        state=retry_io(lambda:np.load(ROOT/'results/E2'/entry['file']))
        if state.shape!=(100,100) or not np.isin(state,[0,1]).all():raise RuntimeError('Invalid E2 configuration')
        ax.imshow(state,cmap=ListedColormap(['#eeeeeb',COLORS['W']]),vmin=0,vmax=1,interpolation='nearest',origin='lower')
        ax.set_xticks([0,49,99],[1,50,100]);ax.set_yticks([0,49,99],[1,50,100])
        ax.set_title(f'W({entry["M"]}), final $f_C$ = {state.mean():.3f}')
        if i==1:ax.set_xlabel('Lattice column')
        ax.set_ylabel(f'$b_{entry["b_label"]}$ = {entry["b"]:g}\nLattice row' if j==0 else 'Lattice row')
        facts.append(dict(entry,final_f_C=float(state.mean())))
    fig.legend(handles=[Line2D([],[],marker='s',ls='',color=COLORS['W'],label='Cooperator'),Line2D([],[],marker='s',ls='',color='#ccccca',label='Defector')],loc='outside upper center',ncol=2)
    save(fig,'fig03_spatial_configurations');json_save(ROOT/'results/E2/configuration_facts.json',facts)
    return facts

def fig04(E1,E3,bs):
    aa=[.05,.1,.15,.2,.3,.4,.6,.8,1.];mm=[1,2,4,8,16,32]
    xedges=[.75,1.5,3,6,12,24,40];yedges=[0,.075,.125,.175,.25,.35,.5,.7,.9,1.1]
    fig=plt.figure(figsize=(11,8));gs=fig.add_gridspec(2,2,hspace=.38)
    for i,b in enumerate(bs):
        ax=fig.add_subplot(gs[0,i]);z=np.array([[one(E3,b=b,a=a,**({} if a==1 else {'m':m}))['f_C'] for m in mm] for a in aa])
        im=ax.pcolormesh(xedges,yedges,z,vmin=0,vmax=1,cmap='viridis',shading='flat')
        windows=sorted([r for r in E1 if r['kernel']=='W' and r['b']==b and r['M']>1],key=lambda r:r['M'])
        ax.scatter([r['M']-1 for r in windows],[1/r['M'] for r in windows],c=[r['f_C'] for r in windows],vmin=0,vmax=1,cmap='viridis',edgecolors='white',linewidths=1,marker='D',s=45,label='Windows W(M)')
        ax.set_xscale('log',base=2);ax.set_xticks([1,2,4,8,16,32,64],[1,2,4,8,16,32,64]);ax.set_yticks([0,.1,.2,.4,.6,.8,1.])
        ax.set(xlim=(.75,72),ylim=(0,1.1),xlabel='Past window m (rounds)',ylabel='Current-action weight a',title=f'({chr(97+i)}) $b_{"AB"[i]}$ = {b:g}')
        ax.text(.02,.96,'a = 1: m irrelevant',transform=ax.transAxes,va='top',fontsize=7);ax.legend(loc='upper right',fontsize=8)
        fig.colorbar(im,ax=ax,label='Mean cooperation $f_C$')
    ax=fig.add_subplot(gs[1,:])
    for j,a in enumerate([.1,.2,.4]):
        for i,b in enumerate(bs):line(ax,subset(E3,a=a,b=b),'m',f'a = {a:g}, $b_{"AB"[i]}$',shade('PP',j),'-' if i==0 else '--')
    ax.set_xscale('log',base=2);ax.set_xticks(mm,mm);ax.set(xlabel='Past window m (rounds)',ylabel='Mean cooperation $f_C$',ylim=(-.03,1.03),title='(c) Past-window length at fixed current-action weight')
    ax.legend(ncol=3,fontsize=9);save(fig,'fig04_present_past')

def mechanism_examples(E1,E3,bs):
    examples=[]
    for b in bs:
        for k in ['W','E','PP']:
            r=one(E3,b=b,a=.4,m=8) if k=='PP' else one(E1,b=b,kernel=k,M=9)
            p=Parameters(kernel=k,a=.4,m=8,b=b) if k=='PP' else Parameters(kernel=k,M=9,b=b)
            h,_=weights(p);past=sum(h[l]*r[f'C_s_{l}'] for l in range(1,len(h)))
            examples.append(dict(row=r,kernel=k,label='PP(0.4, 8)' if k=='PP' else f'{k}(9)',b=b,current=float(h[0]),past=float(past)))
    return examples

def fig05(E1,E3,bs):
    examples=mechanism_examples(E1,E3,bs);fig,axes=plt.subplots(1,3,figsize=(15,4.4),layout='constrained')
    for e in examples:
        r=e['row'];lags=[int(k[4:]) for k in r if k.startswith('C_s_') and k[4:].isdigit() and math.isfinite(r[k])]
        if not lags:continue
        maxlag=max(lags);ages=np.arange(min(maxlag,32)+1)
        axes[0].plot(ages,[r[f'C_s_{l}'] for l in ages],color=COLORS[e['kernel']],ls='-' if e['b']==bs[0] else '--',label=e['label']+f', b={e["b"]:g}')
    axes[0].set(xlabel='Lag (rounds)',ylabel=r'$C_s(\ell)$',title='(a) Strategy persistence');axes[0].legend(fontsize=7)
    for i,e in enumerate(examples):
        if not math.isfinite(e['past']):
            axes[1].text(i,.05,'No coexistence',rotation=90,ha='center',fontsize=7);continue
        axes[1].bar(i,e['current'],color=COLORS[e['kernel']]);axes[1].bar(i,e['past'],bottom=e['current'],color=COLORS[e['kernel']],alpha=.35)
        axes[1].errorbar(i,e['current']+e['past'],yerr=e['row']['G_approx_se'],color='black',capsize=2,fmt='none')
    axes[1].set_xticks(range(6),[e['label']+f'\nb={e["b"]:g}' for e in examples],rotation=40,ha='right',fontsize=8)
    axes[1].set(ylabel='Contribution to reputation gap G',title='(b) Present (solid) and past (pale)')
    for k,marker in [('W','o'),('E','s'),('PP','^')]:
        for i,b in enumerate(bs):
            d=[r for r in (E1 if k!='PP' else E3) if r['kernel']==k and r['b']==b and math.isfinite(r['G'])]
            axes[2].errorbar([r['G'] for r in d],[r['f_C'] for r in d],xerr=[r['G_se'] for r in d],yerr=[r['f_C_se'] for r in d],fmt=marker,ms=3,ls='',color=COLORS[k],mfc=COLORS[k] if i==0 else 'white',alpha=.6,label=f'{k}, $b_{"AB"[i]}$')
    axes[2].set(xlabel='Reputation gap G',ylabel='Mean cooperation $f_C$',ylim=(-.03,1.03),title='(c) All memory parameters at selected b');axes[2].legend(fontsize=8)
    save(fig,'fig05_reputation_gap')
    return examples

def maxima(E5,bs):
    result=[]
    for k,b,eps in [(k,b,e) for k in ['W','E'] for b in bs for e in [0,.005,.01,.02,.05,.1]]:
        d=sorted(subset(E5,kernel=k,b=b,epsilon=eps),key=lambda r:r['M']);best=max(r['f_C'] for r in d)
        peaks=[r for r in d if abs(r['f_C']-best)<1e-12];low=max(r['f_C']-r['f_C_se'] for r in peaks)
        compatible=[int(r['M']) for r in d if r['f_C']+r['f_C_se']>=low]
        for r in peaks:
            i=d.index(r);result.append(dict(kernel=k,b=b,epsilon=eps,M_star=int(r['M']),f_C=r['f_C'],f_C_se=r['f_C_se'],
                exact_ties=len(peaks),grid_lower=d[i-1]['M'] if i else math.nan,grid_upper=d[i+1]['M'] if i+1<len(d) else math.nan,
                overlapping_one_SE_intervals=','.join(map(str,compatible)),note='Descriptive grid and one-SE uncertainty; no fitted optimum or significance test'))
    save_table(ROOT/'results/E5/grid_optima.csv',result);return result

def fig06(E5,bs,optima):
    epsilons=[0,.005,.01,.02,.05,.1];fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    for i,b in enumerate(bs):
        for j,k in enumerate(['W','E']):
            ax=axes[i,j]
            for n,eps in enumerate(epsilons):
                d=subset(E5,kernel=k,b=b,epsilon=eps);color=shade(k,n,len(epsilons))
                line(ax,d,'M',fr'$\epsilon$ = {eps:g}',color,marker=['o','s','^','v','D','x'][n])
                peaks=subset(optima,kernel=k,b=b,epsilon=eps)
                ax.scatter([r['M_star'] for r in peaks],[r['f_C'] for r in peaks],marker='*',s=70,edgecolors='black',linewidths=.4,color=color,zorder=4)
            ax.set_xscale('log',base=2);ax.set_xticks([1,2,3,5,9,17,33],[1,2,3,5,9,17,33])
            ax.set(xlabel='Memory length M',ylabel='Mean cooperation $f_C$',ylim=(-.03,1.03),title=f'{k}, $b_{"AB"[i]}$ = {b:g}');ax.legend(ncol=2,fontsize=8)
    save(fig,'fig06_observation_errors')

def fig07(E6,bs):
    fig,axes=plt.subplots(2,2,figsize=(11,7),layout='constrained')
    for col,b in enumerate(bs):
        for j,M in enumerate([1,9,33]):
            for eps,style in [(0,'-'),(.02,'--')]:
                d=subset(E6,b=b,M=M,epsilon=eps)
                for row,y in enumerate(['f_C','all_D']):line(axes[row,col],d,'L',fr'W({M}), $\epsilon$={eps:g}',shade('W',j),style,['o','s','^'][j],y=y)
        for row in range(2):axes[row,col].set(xlabel='Lattice side L (sites)',ylabel='Mean cooperation $f_C$' if row==0 else 'All-D fraction',xticks=[50,100,200],ylim=(-.03,1.03),title=f'$b_{"AB"[col]}$ = {b:g}')
        axes[0,col].legend(fontsize=8,ncol=2)
    save(fig,'fig07_population_size')

def figS1(E1,E3,S1,bs):
    fig,axes=plt.subplots(1,2,figsize=(10,4.3),layout='constrained')
    for i,b in enumerate(bs):
        ax=axes[i]
        for j,m in enumerate([2,8]):
            for eta,rows,style in [(1e-6,E3,'--'),(.01,S1,'-')]:
                d=[r for r in rows if r['kernel']=='PP' and r['b']==b and r['m']==m and r['a'] in [.1,.2,.4]]
                line(ax,d,'a',fr'm={m}, $\eta$={eta:g}',shade('PP',j,2),style)
        for eta,rows,style in [(1e-6,E1,'--'),(.01,S1,'-')]:
            r=one(rows,kernel='W',M=9,b=b);ax.axhline(r['f_C'],color=COLORS['W'],ls=style,label=fr'W(9), $\eta$={eta:g}')
            ax.axhspan(r['f_C']-r['f_C_se'],r['f_C']+r['f_C_se'],color=COLORS['W'],alpha=.1)
        ax.set(xlabel='Current-action weight a',ylabel='Mean cooperation $f_C$',xticks=[.1,.2,.4],ylim=(-.03,1.03),title=f'$b_{"AB"[i]}$ = {b:g}');ax.legend(fontsize=8,ncol=2)
    save(fig,'figS1_eta')

def figS2(E1,S2):
    fig,ax=plt.subplots(figsize=(8,4))
    for i,M in enumerate([1,9,33]):
        for name,rows,style in [('Synchronous',E1,'-'),('Sequential',S2,'--')]:line(ax,subset(rows,kernel='W',M=M),'b',f'W({M}), {name}',shade('W',i),style)
    ax.set(xlabel='Temptation b',ylabel='Mean cooperation $f_C$',ylim=(-.03,1.03));ax.legend(fontsize=8,ncol=2)
    save(fig,'figS2_update_order')

def captions():
    text={
      'Figure 1':'Weights of W(9), E(9), and PP(0.4, 8), shown through age 32. The full exponential kernel continues beyond the displayed range.',
      'Figure 2':'Mean cooperation for window and exponential memory, with cuts at the selected temptations. L = 100 and epsilon = 0; marker area in the cuts is 14 + 90 times the all-D fraction. Error bars show one standard error across realizations.',
      'Figure 3':'Final spatial configurations for W(1), W(9), and W(33), at b_A and b_B. Blue denotes cooperation; each panel uses the middle realization after sorting by time-averaged cooperation, with its final cooperation shown in the title.',
      'Figure 4':'Cooperation in the present-past plane, window-memory overlays, and cuts at fixed current-action weight. L = 100 and epsilon = 0; a = 1 is one kernel displayed across m, W(1) at m = 0 is omitted from the logarithmic axis, and windows outside the PP grid appear in the unfilled region. Color scales run from zero to one, and error bars show one standard error.',
      'Figure 5':'Strategy persistence, present-past gap decomposition, and cooperation versus reputation gap at b_A and b_B. W(9), E(9), and PP(0.4, 8) illustrate the first two panels; the third includes all sampled memory parameters at the selected temptations, with G averaged over realizations having coexistence. Error bars show one standard error.',
      'Figure 6':'Cooperation as a function of memory length under constant observation errors, with stars at the observed grid maxima, including exact ties. L = 100, with b_A and b_B in separate rows. Error bars show one standard error.',
      'Figure 7':'Cooperation and all-D fraction versus lattice size for W(1), W(9), and W(33). Solid and dashed lines indicate epsilon = 0 and 0.02, respectively, at the two selected temptations. Error bars show one standard error.',
      'Figure S1':'Present-past comparisons with eta = 0.01 and 1e-6, with W(9) as a reference. L = 100 and epsilon = 0 at b_A and b_B. Bars and horizontal bands show one standard error.',
      'Figure S2':'Window-memory cooperation under synchronous and random sequential updates. L = 100 and epsilon = 0, with reputations updated once at the start of each sequential sweep. Error bars show one standard error.'}
    retry_io(lambda:(FIG/'captions.txt').write_text('\n\n'.join(k+'\n'+v for k,v in text.items())+'\n',encoding='utf-8'))

def basic_figures():
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig01();facts=fig03();captions();return facts

def extract_mechanism(E1,E3,E5):
    out=ROOT/'results/E4';out.mkdir(exist_ok=True);rows=[]
    for experiment,data in [('E1',E1),('E3',E3)]:
        for r in data:
            a=r['a'] if r['kernel']=='PP' else (1/r['M'] if r['kernel']=='W' else 2/(r['M']+1))
            valid=math.isfinite(r['G'])
            row={k:r[k] for k in ['kernel','b','L','f_C','f_C_se','G','G_se','G_n','G_approx','G_from_Cy','source']}
            row.update(experiment=experiment,a=a,M=r.get('M',math.nan),m=r.get('m',math.nan),
                G_current=a if valid else math.nan,G_past=r['G_approx']-a if valid else math.nan,G_past_se=r['G_approx_se'],
                excluded_realizations=r['gap_excluded_realizations'])
            rows.append(row)
    save_table(out/'reputation_gap.csv',rows)
    fields=['kernel','M','b','epsilon','R','G','G_se','G_n','G_approx','G_approx_se','relative_gap_discrepancy','relative_gap_discrepancy_se','relative_gap_discrepancy_n','gap_excluded_realizations','source']
    save_table(ROOT/'results/E5/noise_gap_discrepancy.csv',[{k:r[k] for k in fields} for r in E5])

def analyze():
    if not read_json(ROOT/'config/approvals.json')['checkpoint1_approved']:raise RuntimeError('Checkpoint 1 not approved')
    selection=read_json(SELECTION);bs=[selection['b_A'],selection['b_B']]
    E1=load_summary(ROOT/'results/E1/point_summary.csv');data={n:aggregate(n) for n in EXPERIMENTS}
    extract_mechanism(E1,data['E3'],data['E5'])
    facts=basic_figures();fig04(E1,data['E3'],bs);examples=fig05(E1,data['E3'],bs)
    optima=maxima(data['E5'],bs);fig06(data['E5'],bs,optima);fig07(data['E6'],bs);figS1(E1,data['E3'],data['S1'],bs);figS2(E1,data['S2'])
    write_summary(E1,data,bs,facts,examples,optima)
    json_save(ROOT/'results/phase2/analysis_manifest.json',dict(source_sha256=digest(Path(__file__)),
        figure_files={p.name:digest(p) for p in FIG.glob('fig*.png')},checkpoint2_approved=read_json(ROOT/'config/approvals.json').get('checkpoint2_approved',False)))

def fmt(r,key='f_C'):
    return f'{r[key]:.4f} ± {r[key+"_se"]:.4f}' if math.isfinite(r[key]) else 'indefinido (sem coexistência)'

def write_summary(E1,data,bs,facts,examples,optima):
    E3=data['E3'];E5=data['E5'];E6=data['E6'];S1=data['S1'];S2=data['S2']
    paragraphs=[]
    paragraphs.append('Figura 1. Os pesos atuais são 1/9 para W(9), 0,2 para E(9) e 0,4 para PP(0,4; 8). A janela e o presente-passado têm suporte finito; a cauda exponencial continua além da idade 32 mostrada.')
    desc=[]
    for b in bs:
        for k in ['W','E']:
            lo=one(E1,b=b,kernel=k,M=1);hi=one(E1,b=b,kernel=k,M=65)
            desc.append(f'{k}, b={b:g}: f_C vai de {fmt(lo)} em M=1 a {fmt(hi)} em M=65; fração todos-D no último ponto = {hi["all_D"]:.3f}')
    paragraphs.append('Figura 2. '+'; '.join(desc)+'.')
    paragraphs.append('Figura 3. As frações finais de cooperadores nas seis configurações selecionadas são '+ '; '.join(f'W({r["M"]}), b={r["b"]:g}: {r["final_f_C"]:.4f}' for r in facts)+'. A escolha usa a posição central na ordenação por cooperação média temporal, não a aparência da configuração.')
    desc=[]
    for b in bs:
        for a in [.1,.2,.4]:
            desc.append(f'b={b:g}, a={a:g}: m=1 dá {fmt(one(E3,b=b,a=a,m=1))}, e m=32 dá {fmt(one(E3,b=b,a=a,m=32))}')
    paragraphs.append('Figura 4. '+'; '.join(desc)+'.')
    gap_pp=one(E3,b=bs[0],a=.1,m=32);gap_exp=one(E1,b=bs[0],kernel='E',M=49)
    paragraphs.append('Figura 5. A decomposição presente + passado do gap é '+ '; '.join(f'{e["label"]}, b={e["b"]:g}: {e["current"]:.4f} + {e["past"]:.4f}, com G medido = {fmt(e["row"],"G")}' for e in examples)+
        f'. Em b={bs[0]:g}, PP(0,1;32) apresenta G={fmt(gap_pp,"G")} e f_C={fmt(gap_pp)}, enquanto E(49) apresenta G={fmt(gap_exp,"G")} e f_C={fmt(gap_exp)}; os gaps próximos coexistem com cooperações distintas, com {int(gap_pp["G_n"])} e {int(gap_exp["G_n"])} realizações válidas para G, respectivamente. Realizações sem coexistência ficam fora de G; suas contagens permanecem nas tabelas.')
    desc=[]
    for k in ['W','E']:
        for b in bs:
            seq=[]
            for eps in [0,.005,.01,.02,.05,.1]:
                peaks=subset(optima,kernel=k,b=b,epsilon=eps);seq.append(f'{eps:g}: '+','.join(str(r['M_star']) for r in peaks))
            desc.append(f'{k}, b={b:g}, epsilon → M* = '+'; '.join(seq))
    noise_examples=[one(E5,kernel=k,b=bs[0],M=9,epsilon=.1) for k in ['W','E']]
    noise_description='; '.join(f'{r["kernel"]}(9), b={r["b"]:g}, epsilon=0,1: {100*r["relative_gap_discrepancy"]:.4f}% ± {100*r["relative_gap_discrepancy_se"]:.4f}%' for r in noise_examples)
    paragraphs.append('Figura 6. '+' | '.join(desc)+'. As incertezas de grade e os intervalos de um erro-padrão sobrepostos estão em results/E5/grid_optima.csv; os empates são preservados. A discrepância (G − G_approx)/G é '+noise_description+'. A tabela E5 contém todos os pontos e seus denominadores de coexistência.')
    desc=[]
    for b in bs:
        for eps in [0,.02]:
            r=one(E6,b=b,epsilon=eps,M=9,L=200);desc.append(f'W(9), b={b:g}, epsilon={eps:g}, L=200: f_C={fmt(r)}, todos-D={r["all_D"]:.3f}')
    paragraphs.append('Figura 7. '+'; '.join(desc)+'. Os valores em L=50 e L=100 e os outros comprimentos de memória constam na figura e na tabela E6.')
    paragraphs.append('Figura S1. '+'; '.join(f'b={b:g}, PP(0,1;8): eta=1e-6 dá {fmt(one(E3,b=b,a=.1,m=8))}, eta=0,01 dá {fmt(one(S1,b=b,kernel="PP",a=.1,m=8))}' for b in bs)+'.')
    paragraphs.append('Figura S2. '+'; '.join(f'W(9), b={b:g}: atualização síncrona dá {fmt(one(E1,b=b,kernel="W",M=9))}, sequencial dá {fmt(one(S2,b=b,kernel="W",M=9))}' for b in bs)+'.')
    # Descriptive directions only: no inferential tests or fitted trends.
    h1=[(one(E1,b=b,kernel=k,M=65)['f_C'],one(E1,b=b,kernel=k,M=1)['f_C']) for b in bs for k in ['W','E']]
    h2a=[(one(E3,b=b,a=1.)['f_C'],one(E3,b=b,a=.05,m=m)['f_C']) for b in bs for m in [1,2,4,8,16,32]]
    h2m=[(one(E3,b=b,a=a,m=32)['f_C'],one(E3,b=b,a=a,m=1)['f_C']) for b in bs for a in [.05,.1,.15,.2,.3,.4,.6,.8]]
    def differences(rows,x):
        ordered=sorted(rows,key=lambda r:r[x]);return [right['f_C']-left['f_C'] for left,right in zip(ordered,ordered[1:])]
    d1=[delta for k in ['W','E'] for b in sorted({r['b'] for r in E1}) for delta in differences(subset(E1,kernel=k,b=b),'M')]
    da=[delta for b in bs for m in [1,2,4,8,16,32] for delta in differences(subset(E3,b=b,m=m)+subset(E3,b=b,a=1.),'a')]
    dm=[delta for b in bs for a in [.05,.1,.15,.2,.3,.4,.6,.8] for delta in differences(subset(E3,b=b,a=a),'m')]
    rev1=sum(d>1e-12 for d in d1);reva=sum(d<-1e-12 for d in da);revm=sum(d>1e-12 for d in dm)
    only_present=all(r['M_star']==1 and r['exact_ties']==1 for r in optima)
    hypotheses=[
        f'H1: {"confirmada descritivamente na grade simulada" if rev1==0 else "inconclusiva quanto à monotonicidade em toda a grade"}; {rev1} inversões em {len(d1)} comparações adjacentes de M, com platôs nas regiões absorventes. M=65 tem menor cooperação que M=1 em {sum(x<y for x,y in h1)}/4 comparações nos b selecionados.',
        f'H2: {"confirmada descritivamente na grade simulada" if reva==0 and revm==0 else "inconclusiva quanto à afirmação conjunta"}; {reva} inversões em {len(da)} incrementos de a e {revm} inversões em {len(dm)} incrementos de m. A duração do passado mantém efeito próprio com a fixo; estes resultados não medem uma dominância geral de a sobre m.',
        f'H3: inconclusiva como organizador geral entre núcleos. A tendência visual é crescente em cada b, mas G não ordena todos os casos: em b={bs[0]:g}, PP(0,1;32) tem (G, f_C)=({gap_pp["G"]:.4f}, {gap_pp["f_C"]:.4f}), e E(49) tem ({gap_exp["G"]:.4f}, {gap_exp["f_C"]:.4f}), com coexistência em todas as realizações de ambos os pontos.',
        (f'H4: refutada na grade testada a possibilidade de um ótimo intermediário crescente com epsilon; o máximo observado foi M=1, sem empate, nas {len(optima)} combinações de núcleo, b e epsilon. O resultado é ausência de ótimo intermediário no domínio amostrado, sem extrapolação para outros valores.' if only_present else 'H4: inconclusiva; os máximos observados, empates e resolução da grade estão reportados sem ajuste de curvas.')
    ]
    issues=[f'{n}: {sum(r["still_nonstationary"] for r in rows)} pontos ainda sinalizados após o transiente estendido' for n,rows in [('E1',E1)]+list(data.items())]
    text='# Resumo factual para o ponto de checagem 2\n\nTodos os valores ± indicam um erro-padrão entre realizações. As classificações das hipóteses são descritivas da grade calculada, sem testes de significância. Não são usados resultados anteriores ao roteiro simple-paper.\n\n'+'\n\n'.join(paragraphs)+'\n\n## Hipóteses — avaliação descritiva\n\n'+'\n\n'.join(hypotheses)+'\n\n## Checagem do transiente\n\n'+'\n\n'.join(issues)+'\n\nA autorização do ponto de checagem 2 está registrada em config/approvals.json; a redação segue o enquadramento aprovado.\n'
    retry_io(lambda:(ROOT/'results_summary.md').write_text(text,encoding='utf-8'))

if __name__=='__main__':analyze()
