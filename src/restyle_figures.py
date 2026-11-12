"""Presentation-only redraw from existing summary tables; never runs simulations."""
from pathlib import Path
import json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.text import Text
from PIL import Image
from analyze_phase2 import load_summary,one,subset,mechanism_examples
from model import Parameters,weights
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/figures';QA=ROOT/'qa/form_review'
COLORS={'W':'#0072B2','E':'#D55E00','PP':'#009E73'}
MARKERS=['o','s','^','D','v','P'];KM={'W':'o','E':'s','PP':'^'}
CM=LinearSegmentedColormap.from_list('viridis_dark',plt.colormaps['viridis'](np.linspace(0,.9,256)))
plt.rcParams.update({'font.size':9.5,'axes.labelsize':9.5,'axes.titlesize':9.5,'xtick.labelsize':9,'ytick.labelsize':9,'legend.fontsize':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
AUDIT={}
def ordered(n):return plt.colormaps['viridis'](np.linspace(.05,.85,n))
def save(fig,name,width=16.6):
    fig.canvas.draw();renderer=fig.canvas.get_renderer();box=fig.get_tightbbox(renderer)
    # PDF bbox adds 0.1 in on each side. Audit sizes after scaling in the paper.
    scale=(width/2.54)/(box.width+.2)
    texts=[t for t in fig.findobj(Text) if t.get_visible() and t.get_text()]
    minimum=min(t.get_fontsize()*scale for t in texts)
    if minimum<8:
        for t in texts:t.set_fontsize(t.get_fontsize()*8.15/minimum)
        fig.canvas.draw();box=fig.get_tightbbox(fig.canvas.get_renderer());scale=(width/2.54)/(box.width+.2)
        minimum=min(t.get_fontsize()*scale for t in texts)
    assert minimum>=8,(name,minimum)
    AUDIT[name]={'print_width_cm':width,'minimum_print_font_pt':minimum}
    for ext in ['pdf','png','svg']:fig.savefig(OUT/f'{name}.{ext}',bbox_inches='tight',dpi=220)
    plt.close(fig)
def line(ax,rows,x,color,marker,label='',ls='-',open_=False,y='f_C',offset=0):
    rows=sorted(rows,key=lambda r:r[x]);xx=np.array([r[x] for r in rows]);yy=[r[y] for r in rows]
    # Curves retain exact abscissas; offsets only separate overlapping error bars/markers.
    ax.plot(xx,yy,color=color,ls=ls,lw=1.25,zorder=1)
    ax.errorbar(xx*(1+offset) if x=='M' else xx+offset,yy,yerr=[r[y+'_se'] for r in rows],ls='',marker=marker,ms=4.5,mfc='white' if open_ else color,mec=color,mew=.95,color=color,capsize=2,label=label,zorder=3)
def legend(ax,ncol=2):ax.legend(loc='lower center',bbox_to_anchor=(.5,1.02),ncol=ncol,frameon=False,columnspacing=1,handlelength=2)
def fig01():
    fig,ax=plt.subplots(figsize=(4.8,2.3),layout='constrained')
    for k,p in [('W',Parameters(M=9)),('E',Parameters(kernel='E',M=9)),('PP',Parameters(kernel='PP',a=.4,m=8))]:
        h,_=weights(p);x=np.arange(min(len(h),33));ax.plot(x,h[:len(x)],color=COLORS[k],marker=KM[k],ms=3,label='PP(0.4, 8)' if k=='PP' else k+'(9)')
    ax.set(xlabel='Record age (rounds)',ylabel=r'Weight $h_\ell$',ylim=(-.012,.43));ax.legend(frameon=False)
    save(fig,'fig01_memory_weights',16.6*.77)
def fig02(D,bs):
    cfg=json.loads((ROOT/'config/E1.json').read_text());fig=plt.figure(figsize=(6.35,4.9),layout='constrained');gs=fig.add_gridspec(2,2,height_ratios=[1,1.1])
    for i,k in enumerate(['W','E']):
        ax=fig.add_subplot(gs[0,i]);z=[[one(D,kernel=k,M=M,b=b)['f_C'] for b in cfg['b']] for M in cfg['M']]
        im=ax.imshow(z,origin='lower',aspect='auto',cmap=CM,vmin=0,vmax=1)
        ax.set_xticks(range(0,25,8),[f'{b:.2f}' for b in cfg['b'][::8]]);ax.set_yticks(range(12),cfg['M'])
        for j,b in enumerate(bs):ax.axvline(cfg['b'].index(b),color='white',ls=['-','--'][j],lw=1.2)
        ax.set(xlabel='Temptation $b$',ylabel='Memory length $M$',title=f'({chr(97+i)}) {"Window" if k=="W" else "Exponential"}');fig.colorbar(im,ax=ax,label='$f_C$',ticks=[0,.5,1],fraction=.05)
    ax=fig.add_subplot(gs[1,:])
    for k in ['W','E']:
        for j,b in enumerate(bs):
            d=sorted(subset(D,kernel=k,b=b),key=lambda r:r['M']);x=[r['M'] for r in d];y=[r['f_C'] for r in d]
            ax.errorbar(x,y,yerr=[r['f_C_se'] for r in d],ls=['-','--'][j],color=COLORS[k],capsize=2,lw=1.25)
            ax.scatter(np.asarray(x)*(.95 if k=='W' else 1.05),y,s=[22+58*r['all_D'] for r in d],marker=KM[k],facecolors='white' if j else COLORS[k],edgecolors=COLORS[k],linewidth=.95,zorder=3,label=f'{k}, $b_{"AB"[j]}$')
    ax.set_xscale('log',base=2);ax.set_xticks(cfg['M'],cfg['M']);ax.set(xlabel='Memory length $M$',ylabel='Mean cooperation $f_C$',ylim=(-.04,1.04),title='(c) Cuts; symbol area encodes all-D fraction')
    ax.legend(ncol=2,loc='upper right',frameon=False);save(fig,'fig02_phase_map')
def fig03():
    entries=json.loads((ROOT/'results/E2/selection.json').read_text());fig,axes=plt.subplots(2,3,figsize=(6.35,4.4),layout='constrained')
    for e in entries:
        i='AB'.index(e['b_label']);j=[1,9,33].index(e['M']);ax=axes[i,j];z=np.load(ROOT/'results/E2'/e['file'])
        ax.imshow(z,cmap=ListedColormap(['#eeeeeb',COLORS['W']]),origin='lower',vmin=0,vmax=1,interpolation='nearest')
        ax.set_xticks([0,49,99],[1,50,100]);ax.set_yticks([0,49,99],[1,50,100]);ax.set_title(f'W({e["M"]}), $f_C$={z.mean():.3f}')
        if i==1:ax.set_xlabel('Column')
        if j==0:ax.set_ylabel(f'$b_{e["b_label"]}$={e["b"]:g}\nRow')
    fig.legend(handles=[Line2D([],[],ls='',marker='s',color=COLORS['W'],label='Cooperator'),Line2D([],[],ls='',marker='s',color='#bbbbbb',label='Defector')],loc='outside upper center',ncol=2,frameon=False)
    save(fig,'fig03_spatial_configurations')
def fig04(E1,E3,bs):
    aa=[.05,.1,.15,.2,.3,.4,.6,.8,1];mm=[1,2,4,8,16,32];fig=plt.figure(figsize=(6.35,5.0),layout='constrained');gs=fig.add_gridspec(2,2)
    for i,b in enumerate(bs):
        ax=fig.add_subplot(gs[0,i]);z=[[one(E3,b=b,a=a,**({} if a==1 else {'m':m}))['f_C'] for m in mm] for a in aa]
        im=ax.pcolormesh([.75,1.5,3,6,12,24,40],[0,.075,.125,.175,.25,.35,.5,.7,.9,1.1],z,cmap=CM,vmin=0,vmax=1)
        ww=[r for r in E1 if r['kernel']=='W' and r['b']==b and r['M']>1];xx=[r['M']-1 for r in ww];yy=[1/r['M'] for r in ww]
        ax.scatter(xx,yy,marker='D',s=36,color='white',zorder=3)
        ax.scatter(xx,yy,marker='D',s=19,c=[r['f_C'] for r in ww],cmap=CM,vmin=0,vmax=1,edgecolors=COLORS['W'],linewidth=.9,zorder=4)
        ax.set_xscale('log',base=2);ax.set_xticks([1,4,16,64],[1,4,16,64]);ax.set_yticks([0,.2,.4,.6,.8,1])
        ax.set(xlim=(.75,72),ylim=(0,1.1),xlabel='Past rounds $m$',ylabel='Current weight $a$',title=f'({chr(97+i)}) $b_{"AB"[i]}$ = {b:g}')
        fig.colorbar(im,ax=ax,label='$f_C$',ticks=[0,.5,1],fraction=.05)
    ax=fig.add_subplot(gs[1,:])
    for j,a in enumerate([.1,.2,.4]):
        for i,b in enumerate(bs):line(ax,subset(E3,a=a,b=b),'m',ordered(3)[j],MARKERS[j],f'$a$={a:g}, $b_{"AB"[i]}$',ls=['-','--'][i],open_=bool(i))
    ax.set_xscale('log',base=2);ax.set_xticks(mm,mm);ax.set(xlabel='Number of past rounds $m$',ylabel='Mean cooperation $f_C$',ylim=(-.04,1.04),title='(c) Fixed weight of current action')
    ax.legend(ncol=3,loc='upper right',frameon=False);save(fig,'fig04_present_past')
def fig05(E1,E3,bs):
    ex=mechanism_examples(E1,E3,bs);fig=plt.figure(figsize=(6.35,5.5),layout='constrained');gs=fig.add_gridspec(2,2,height_ratios=[1,1]);axes=[fig.add_subplot(gs[0,0]),fig.add_subplot(gs[0,1]),fig.add_subplot(gs[1,:])]
    for e in ex:
        r=e['row'];ages=[l for l in range(33) if f'C_s_{l}' in r and np.isfinite(r[f'C_s_{l}'])];j=bs.index(e['b']);k=e['kernel']
        axes[0].plot(ages,[r[f'C_s_{l}'] for l in ages],color=COLORS[k],ls=['-','--'][j],marker=KM[k],mfc='white' if j else COLORS[k],ms=3.5,markevery=(["W","E","PP"].index(k)*2+j,6),lw=1.2,label=f'{k}, $b_{"AB"[j]}$')
    axes[0].set(xlabel=r'$\ell$ (rounds)',ylabel=r'$C_s(\ell)$',title='(a) Strategy persistence',ylim=(-.03,1.03));axes[0].legend(ncol=2,loc='upper right',frameon=False,handlelength=1.5,columnspacing=.6)
    for i,e in enumerate(ex):
        axes[1].bar(i,e['current'],color=COLORS[e['kernel']],edgecolor='black',linewidth=.4)
        axes[1].bar(i,e['past'],bottom=e['current'],color=COLORS[e['kernel']],alpha=.3,edgecolor='black',hatch='///',linewidth=.4)
        axes[1].errorbar(i,e['current']+e['past'],yerr=e['row']['G_approx_se'],fmt='none',color='black',capsize=2)
    axes[1].set_xticks(range(6),[e['kernel']+'\n$b_'+('A' if e['b']==bs[0] else 'B')+'$' for e in ex]);axes[1].set(ylabel='Contribution to $G$',title='(b) Present / past (hatched)')
    for k in ['W','E','PP']:
        for i,b in enumerate(bs):
            d=[r for r in (E3 if k=='PP' else E1) if r['kernel']==k and r['b']==b and np.isfinite(r['G'])]
            axes[2].errorbar([r['G'] for r in d],[r['f_C'] for r in d],xerr=[r['G_se'] for r in d],yerr=[r['f_C_se'] for r in d],ls='',marker=KM[k],ms={'W':7,'E':5.5,'PP':4}[k],mfc='white' if i else COLORS[k],mec=COLORS[k],color=COLORS[k],mew=.9,label=f'{k}, $b_{"AB"[i]}$')
    axes[2].set(xlabel='Reputation gap $G$',ylabel='Mean cooperation $f_C$',ylim=(-.04,1.04),title='(c) Cooperation and reputation gap');axes[2].legend(ncol=3,loc='upper left',frameon=False)
    save(fig,'fig05_reputation_gap')
def fig06(E5,bs):
    fig,axes=plt.subplots(2,2,figsize=(6.35,6.0),layout='constrained');eps=[0,.005,.01,.02,.05,.1]
    for i,b in enumerate(bs):
        for j,k in enumerate(['W','E']):
            ax=axes[i,j]
            for n,e in enumerate(eps):
                d=subset(E5,kernel=k,b=b,epsilon=e)
                line(ax,d,'M',ordered(6)[n],MARKERS[n],fr'$\varepsilon$={e:g}',ls=['-','--'][i],open_=bool(i),offset=(n-2.5)*.04)
            ax.set_xscale('log',base=2);ax.set_xticks([1,2,3,5,9,17,33],[1,2,3,5,9,17,33]);ax.set(xlabel='Memory length $M$',ylabel='$f_C$',ylim=(-.04,1.04),title=f'({chr(97+2*i+j)}) {k}, $b_{"AB"[i]}$={b:g}')
            ax.text(.98,.94,'$M^*=1$',transform=ax.transAxes,ha='right',va='top')
    handles=[Line2D([],[],color=c,marker=m,ls='-',label=fr'$\varepsilon$={e:g}') for c,m,e in zip(ordered(6),MARKERS,eps)]
    fig.legend(handles=handles,ncol=3,loc='outside upper center',frameon=False);save(fig,'fig06_observation_errors')
def fig07(E6,bs):
    fig,axes=plt.subplots(2,2,figsize=(6.35,4.5),layout='constrained')
    for col,b in enumerate(bs):
        for j,M in enumerate([1,9,33]):
            for i,e in enumerate([0,.02]):
                for row,y in enumerate(['f_C','all_D']):line(axes[row,col],subset(E6,b=b,M=M,epsilon=e),'L',ordered(3)[j],MARKERS[j],ls=['-','--'][i],open_=bool(i),y=y,offset=(2*j+i-2.5)*3.5)
        for row in range(2):axes[row,col].set(xlabel='Lattice side $L$',ylabel='$f_C$' if row==0 else 'All-D fraction',xticks=[50,100,200],ylim=(-.055,1.055),title=f'({chr(97+2*row+col)}) $b_{"AB"[col]}$={b:g}')
    handles=[Line2D([],[],color=c,marker=marker,ls=ls,mfc=c if e==0 else 'white',mec=c,label=fr'$M={M},\ \varepsilon={e:g}$') for c,marker,M in zip(ordered(3),MARKERS,[1,9,33]) for ls,e in [('-',0),('--',.02)]]
    fig.legend(handles=handles,ncol=3,loc='outside upper center',frameon=False);save(fig,'fig07_population_size')
def supplements(E1,E3,S1,S2,bs):
    fig,axes=plt.subplots(1,2,figsize=(6.35,3.2),layout='constrained')
    for i,b in enumerate(bs):
        for j,m in enumerate([2,8]):
            for n,(eta,rows) in enumerate([(1e-6,E3),(.01,S1)]):
                d=[r for r in subset(rows,kernel='PP',b=b,m=m) if r['a'] in [.1,.2,.4]]
                line(axes[i],d,'a',COLORS['PP'],MARKERS[j],fr'PP: $m$={m}, $\eta$={eta:g}',ls=['-','--'][n],open_=bool(n))
        for n,rows in enumerate([E1,S1]):
            r=one(rows,kernel='W',M=9,b=b);axes[i].axhline(r['f_C'],color=COLORS['W'],ls=['-','--'][n]);axes[i].axhspan(r['f_C']-r['f_C_se'],r['f_C']+r['f_C_se'],color=COLORS['W'],alpha=.15)
        axes[i].set(xlabel='Weight of current action $a$',ylabel='$f_C$',ylim=(-.03,.7),xticks=[.1,.2,.4],title=f'({chr(97+i)}) $b_{"AB"[i]}$ = {b:g}')
    handles=[Line2D([],[],color=COLORS['PP'],marker=MARKERS[j],ls=ls,mfc=COLORS['PP'] if n==0 else 'white',mec=COLORS['PP'],label=fr'PP, $m={m},\ \eta={eta}$') for j,m in enumerate([2,8]) for n,(ls,eta) in enumerate([('-','10^{-6}'),('--','0.01')])]+[Line2D([],[],color=COLORS['W'],ls=ls,label=fr'W(9), $\eta={eta}$') for ls,eta in [('-','10^{-6}'),('--','0.01')]]
    fig.legend(handles=handles,ncol=3,loc='outside upper center',frameon=False);save(fig,'figS1_eta')
    fig,ax=plt.subplots(figsize=(6.35,3),layout='constrained')
    for j,M in enumerate([1,9,33]):
        for i,(label,rows) in enumerate([('Sync.',E1),('Sequential',S2)]):line(ax,subset(rows,kernel='W',M=M),'b',ordered(3)[j],MARKERS[j],f'$M$={M}, {label}',ls=['-','--'][i],open_=bool(i))
    ax.set(xlabel='Temptation $b$',ylabel='Mean cooperation $f_C$',ylim=(-.04,1.04));legend(ax,3);save(fig,'figS2_update_order')
def accessibility():
    # Machado et al. severity 1 deuteranopia matrix, applied in linear RGB.
    matrix=np.array([[.367322,.860646,-.227968],[.280085,.672501,.047413],[-.011820,.042940,.968881]])
    for n in [2,4,5,6,7]:
        p=next(OUT.glob(f'fig{n:02d}_*.png'));a=np.asarray(Image.open(p).convert('RGB'),dtype=float)/255
        lin=np.where(a<=.04045,a/12.92,((a+.055)/1.055)**2.4)
        gray=np.sum(lin*np.array([.2126,.7152,.0722]),axis=-1);gray=np.repeat(gray[:,:,None],3,axis=2)
        for label,data in [('grayscale',gray),('deuteranopia',np.clip(lin@matrix.T,0,1))]:
            rgb=np.where(data<=.0031308,12.92*data,1.055*np.maximum(data,0)**(1/2.4)-.055)
            Image.fromarray(np.uint8(np.clip(rgb,0,1)*255)).save(QA/(p.stem+'_'+label+'.png'))
def main():
    files=[ROOT/f'results/{e}/point_summary.csv' for e in ['E1','E3','E5','E6','S1','S2']]
    before={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files};D={p.parent.name:load_summary(p) for p in files}
    sel=json.loads((ROOT/'results/E1/b_selection.json').read_text());bs=[sel['b_A'],sel['b_B']]
    fig01();fig02(D['E1'],bs);fig03();fig04(D['E1'],D['E3'],bs);fig05(D['E1'],D['E3'],bs);fig06(D['E5'],bs);fig07(D['E6'],bs);supplements(D['E1'],D['E3'],D['S1'],D['S2'],bs);accessibility()
    after={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files};assert before==after
    (QA/'figure_style_audit.json').write_text(json.dumps({'fonts':AUDIT,'unchanged_data_sha256':before,'palette':COLORS,'ordered_viridis_positions':[.05,.85],'deuteranopia_matrix':[[.367322,.860646,-.227968],[.280085,.672501,.047413],[-.011820,.042940,.968881]]},indent=2))
    print(json.dumps(AUDIT,indent=2))
if __name__=='__main__':main()
