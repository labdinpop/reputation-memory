"""Standalone comparison for author review; never edits the article or Fig. S3."""
from pathlib import Path
import json, math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_e1 import read_json,read_rows,digest,json_save
from analyze_e1 import save_table

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/S3_past'


def stats(values):
    x=np.asarray(values,dtype=float);x=x[np.isfinite(x)]
    return (float(x.mean()) if len(x) else math.nan,
            float(x.std(ddof=1)/np.sqrt(len(x))) if len(x)>1 else math.nan,len(x))


def main():
    cfg=read_json(ROOT/'config/S3_past.json')
    state=read_json(OUT/'status.json');assert state['state']=='complete'
    rows=[];realizations=[]
    for m in cfg['m']:
        folder=OUT/'points'/f'PP_a0_m{m:02d}'
        meta=read_json(folder/'point.json');data=read_rows(folder/meta['selected_file'])
        assert len(data)==cfg['R'] and len({d['seed'] for d in data})==cfg['R']
        realizations.extend(data)
        row=dict(kernel='PP',a=0,m=m,L=cfg['L'],b=cfg['b'],eta=cfg['eta'],epsilon=0,
                 R=cfg['R'],T0=meta['T0'],source=meta['source'],
                 repeated=meta['repeated'],still_nonstationary=meta['still_nonstationary'])
        for key in ['f_C','G','all_D','all_C','relative_gap_discrepancy']:
            avg,se,n=stats([d[key] for d in data]);row.update({key:avg,key+'_se':se,key+'_n':n})
        row['all_C_count']=sum(int(d['all_C']) for d in data)
        row['all_D_count']=sum(int(d['all_D']) for d in data)
        row['identity_max_error']=max(float(d['identity_absolute_error']) for d in data)
        rows.append(row)
    save_table(OUT/'point_summary.csv',rows)
    save_table(OUT/'realizations.csv',realizations)
    baseline=[{k:float(v) for k,v in r.items() if k in ['M','f_C','f_C_se','all_C','all_D']}
              for r in read_rows(ROOT/'results/S3_floors/point_summary.csv')
              if float(r['eta'])==cfg['eta'] and float(r['b'])==cfg['b']]
    baseline.sort(key=lambda r:r['M'])
    common=[]
    for r in rows:
        old=next((q for q in baseline if q['M']==r['m']),None)
        if old:common.append(dict(records=r['m'],past_only_fc=r['f_C'],current_included_fc=old['f_C'],difference=r['f_C']-old['f_C']))
    save_table(OUT/'matched_comparison.csv',common)
    top=max(r['f_C'] for r in rows)
    tied=[r['m'] for r in rows if abs(r['f_C']-top)<1e-12]
    finding=('Não apareceu um ótimo intermediário exclusivo: o maior valor também ocorre em m=1.'
             if 1 in tied else 'A maior média ocorre em m='+', '.join(map(str,tied))+'. A incerteza deve ser considerada antes de afirmar um ótimo.')
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(figsize=(7.6,4.8),layout='constrained')
    ax.errorbar([r['M'] for r in baseline],[r['f_C'] for r in baseline],yerr=[r['f_C_se'] for r in baseline],
                color='#0072B2',marker='o',mfc='white',ls='--',capsize=3,label=r'W(M): current action included (saved S3)')
    ax.errorbar([r['m'] for r in rows],[r['f_C'] for r in rows],yerr=[r['f_C_se'] for r in rows],
                color='#009E73',marker='^',ls='-',capsize=3,label=r'PP(0,m): previous actions only (new runs)')
    ax.set_xscale('log',base=2)
    ticks=sorted(set([r['m'] for r in rows]+[r['M'] for r in baseline]))
    ax.set_xticks(ticks,[f'{x:g}' for x in ticks])
    ax.set(xlabel='Number of records (M or m)',ylabel=r'Mean cooperation $f_C$',ylim=(-.025,1.08),
           title=r'Proportional selection: $b=1.15$, $\eta=10^{-6}$, $L=100$, $\varepsilon=0$')
    ax.legend(loc='lower left',frameon=False)
    for ext in ['png','svg']:fig.savefig(OUT/f'comparison.{ext}',dpi=200)
    plt.close(fig)
    text=['# S3: memória apenas das ações anteriores','',finding,'',
          'Resultados separados para avaliação dos autores. Artigo, suplemento e Fig. S3 não foram alterados.',
          '',f'Foram concluídas {len(rows)*cfg["R"]} realizações independentes: 30 por ponto, L=100, b=1.15, eta=1e-6, epsilon=0, seleção proporcional.',
          'A memória PP(0,m) contém exatamente as m ações anteriores e atribui peso zero à ação atual. As demais regras, incluindo K=0.1 com payoff normalizado e inicialização, são as do S3.',
          '', '| m | Cooperação média ± EP | Todos C | Todos D | G (quando definido) | T0 |',
          '|---:|---:|---:|---:|---:|---:|']
    for r in rows:
        gap=f'{r["G"]:.3f}' if np.isfinite(r['G']) else 'indefinido (sem coexistência amostrada)'
        text.append(f'| {r["m"]} | {r["f_C"]:.3f} ± {r["f_C_se"]:.6f} | {r["all_C_count"]}/30 | {r["all_D_count"]}/30 | {gap} | {r["T0"]} |')
    text.extend(['','## Comparação com o S3 já existente','',
                 'A curva W(M) usa o registro atual mais M-1 registros anteriores; PP(0,m) usa m registros anteriores. O eixo compara a quantidade de registros. As grades não são idênticas: comparações numéricas diretas usam apenas os comprimentos comuns. Os pontos antigos foram lidos de results/S3_floors/point_summary.csv, sem novas simulações de W.',
                 '', '| Registros | W, com ação atual | PP, apenas anteriores | Diferença PP − W |',
                 '|---:|---:|---:|---:|'])
    for c in common:text.append(f'| {c["records"]} | {c["current_included_fc"]:.3f} | {c["past_only_fc"]:.3f} | {c["difference"]:.3f} |')
    text.extend(['',f'Maior média observada: {top:.3f}, em m={", ".join(map(str,tied))}.',
                 f'Pontos repetidos com transiente maior: {sum(r["repeated"] for r in rows)}. Pontos ainda com deriva pela regra predefinida: {sum(r["still_nonstationary"] for r in rows)}.',
                 'Barras na figura: ±1 erro-padrão entre realizações. Ausência de coexistência torna G indefinido, não zero.',
                 'Esta comparação isola a exclusão da ação atual no nosso S3; não reproduz integralmente os modelos de Li et al., que também diferem em aprendizagem e payoffs.',
                 '', '![Comparação](comparison.png)',''])
    (OUT/'results_summary.md').write_text('\n'.join(text),encoding='utf-8')
    protected=read_json(OUT/'manifest.json')['protected_sha256']
    assert all(digest(ROOT/k)==v for k,v in protected.items())
    json_save(OUT/'analysis.json',dict(finding=finding,max_mean=top,maximizing_m=tied,
                                     total_realizations=len(realizations),manuscript_unchanged=True))
    print(json.dumps(dict(finding=finding,rows=rows),indent=2),flush=True)


if __name__=='__main__':main()
