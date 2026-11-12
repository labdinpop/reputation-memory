"""Generate manuscript numbers and bibliography from saved results and DOI metadata."""
from pathlib import Path
import csv,json,re,hashlib,html
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'artigo'
def read(name):
    with (ROOT/f'results/{name}/point_summary.csv').open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f))
def main(numbers_only=False):
    OUT.mkdir(exist_ok=True)
    data={n:read(n) for n in ['E1','E3','E5','E6','S1','S2']}
    if (ROOT/'results/S3/point_summary.csv').exists():data['S3']=read('S3')
    if (ROOT/'results/S3_floors/point_summary.csv').exists():data['S3_floors']=read('S3_floors')
    cfg=json.loads((ROOT/'config/E1.json').read_text())
    sel=json.loads((ROOT/'results/E1/b_selection.json').read_text())
    values={};sources={}
    def put(key,value,source):values[key]=str(value);sources[key]=source
    def point(alias,experiment,**criteria):
        rows=[r for r in data[experiment] if all(r[k]==v if isinstance(v,str) else abs(float(r[k])-v)<1e-10 for k,v in criteria.items())]
        assert len(rows)==1,(alias,len(rows))
        r=rows[0]
        for field,label in [('f_C','fc'),('G','gap'),('all_D','abs'),('relative_gap_discrepancy','delta')]:
            scale=100 if field in ['all_D','relative_gap_discrepancy'] else 1
            digits=0 if field=='all_D' else 3
            for suffix in ['', '_se']:
                key=alias+'-'+label+('SE' if suffix else '')
                value=float(r[field+suffix])*scale
                precision=1 if field=='all_D' and abs(value-round(value))>1e-9 else digits
                put(key,f'{value:.{precision}f}',dict(file=r['source'],field=field+suffix,scale=scale))
        # Scaled gaps preserve small standard errors with three decimal places.
        for field,label in [('G','gapMilli'),('G_se','gapMilliSE')]:
            put(alias+'-'+label,f'{1000*float(r[field]):.3f}',dict(file=r['source'],field=field,scale=1000))
        put(alias+'-deltaAbout',f'{100*float(r["relative_gap_discrepancy"]):.1f}',dict(file=r['source'],field='relative_gap_discrepancy',scale=100))
        put(alias+'-gn',int(float(r['G_n'])),dict(file=r['source'],field='G_n'))
        return r
    for suffix,b in [('A',sel['b_A']),('B',sel['b_B'])]:
        put('b'+suffix,f'{b:g}','results/E1/b_selection.json')
        for k in ['W','E']:
            for M in [1,9,33,49,65]:point(f'{k}{M}{suffix}','E1',kernel=k,M=M,b=b)
            point(f'{k}9noise{suffix}','E5',kernel=k,M=9,b=b,epsilon=.1)
            point(f'{k}1noise{suffix}','E5',kernel=k,M=1,b=b,epsilon=.1)
        for a in [.1,.2,.4]:
            for m in [1,8,32]:point(f'PP{a:g}m{m}{suffix}','E3',a=a,m=m,b=b)
        point('eta'+suffix,'S1',kernel='PP',a=.1,m=8,b=b)
        point('seq'+suffix,'S2',kernel='W',M=9,b=b)
        for L in [50,100,200]:
            for M in [1,9,33]:
                for eps in [0,.02]:point(f'L{L}M{M}eps{eps:g}{suffix}','E6',L=L,M=M,b=b,epsilon=eps)
    if 'S3' in data:
        for M in [1,2,3,5,9,17,33]:point(f'S3M{M}','S3',M=M)
    if 'S3_floors' in data:
        for eta,label in [(0,'zero'),(1e-6,'small'),(.01,'large')]:
            for M in [1,2,3,5,9,17,33]:point(f'S3{label}M{M}','S3_floors',M=M,eta=eta)
    point('etaW9B','S1',kernel='W',M=9,b=sel['b_B'])
    for experiment,label in [('E1','sync'),('S2','sequential')]:
        endpoints=[]
        for M in [1,33]:
            r=next(r for r in data[experiment] if r['kernel']=='W' and float(r['M'])==M and float(r['b'])==sel['b_A'])
            endpoints.append(r)
        difference=float(endpoints[0]['f_C'])-float(endpoints[1]['f_C'])
        put('memoryDifference'+label,f'{difference:.3f}',dict(operation='difference_of_means',files=[r['source'] for r in endpoints],field='f_C',scale=1))
    put('liFloor','10^{-6}','refs/1-s2.0-S0303264726001772-main.pdf: page 7, section 5.2 item 3')
    put('li2026S','-0.5','refs/1-s2.0-S0303264726001772-main.pdf: page 7, Eq. 10')
    put('li2025KappaMean','0.1','refs/1-s2.0-S0303264725002400-main.pdf: page 4, section 5.1')
    put('li2026KappaMean','0.5','refs/1-s2.0-S0303264726001772-main.pdf: page 7, section 5.2 item 2')
    coefficient=4*cfg['K']
    put('rawNoiseCoefficient',f'{coefficient:g}',dict(operation='4*K',file='config/E1.json',field='K'))
    for suffix,b in [('A',sel['b_A']),('B',sel['b_B'])]:
        put('rawNoise'+suffix,f'{coefficient*b:g}',dict(operation='4*K*b',files=['config/E1.json','results/E1/b_selection.json'],b_key='b_'+suffix))

    for k in ['L','T0','T1','R','sample_every','w','K']:
        put(k,format(cfg[k],',') if isinstance(cfg[k],int) and cfg[k]>=1000 else f'{cfg[k]:g}' if isinstance(cfg[k],float) else cfg[k],'config/E1.json')
    for key,value in [('eta',cfg['eta']),('tailTolerance',cfg['identity']['EMA_tail_below']),('identityTolerance',cfg['identity']['relative_tolerance'])]:
        exponent=int(round(__import__('math').log10(value)))
        put(key,'10^{'+str(exponent)+'}','config/E1.json')
    put('retryT',format(cfg['stationarity']['retry_T0'],','),'config/E1.json')
    refinement=json.loads((ROOT/'results/transient_refinement/completed.json').read_text())
    put('refineT',format(max(r['meta']['T0'] for r in refinement),','),'results/transient_refinement/completed.json')
    put('stillFlagged',sum(r['meta']['still_nonstationary'] for r in refinement),'results/transient_refinement/completed.json')
    put('refinedPoints',len(refinement),'results/transient_refinement/completed.json')
    flagged={r['source']:r for rows in data.values() for r in rows if str(r['still_nonstationary']).lower() in ['true','1','1.0']}
    put('totalDrift',len(flagged),'qa/form_review/residual_drift.json')
    (ROOT/'qa/form_review/residual_drift.json').write_text(json.dumps([
        {k:r.get(k) for k in ['source','kernel','L','M','b','epsilon','T0','still_nonstationary']} for r in flagged.values()
    ],indent=2),encoding='utf-8')
    constants={0:'binary strategy definition',1:'binary strategy definition',2:'kernel definition',4:'von Neumann degree'}
    def gather(obj,source):
        if isinstance(obj,dict):
            for x in obj.values():gather(x,source)
        elif isinstance(obj,list):
            for x in obj:gather(x,source)
        elif isinstance(obj,(int,float)) and not isinstance(obj,bool):constants[obj]=source
    for path in (ROOT/'config').glob('*.json'):gather(json.loads(path.read_text(encoding='utf-8')),str(path.relative_to(ROOT)))
    for n,source in constants.items():
        put('n'+f'{n:g}',f'{n:g}',source)
    step=min(y-x for x,y in zip(cfg['b'],cfg['b'][1:]))
    put('n'+f'{step:g}',f'{step:g}','config/E1.json: consecutive b differences')
    for name,file,field in [('Mgrid','E1','M'),('Mgrideps','E5','M'),('epsgrid','E5','epsilon'),('agrid','E3','a'),('mgrid','E3','m')]:
        c=json.loads((ROOT/f'config/{file}.json').read_text())
        put(name,', '.join(f'{x:g}' for x in c[field]),f'config/{file}.json')
    opt=list(csv.DictReader((ROOT/'results/E5/grid_optima.csv').open(encoding='utf-8')))
    put('optimaCount',len(opt),'results/E5/grid_optima.csv')
    put('optimum',int(float(opt[0]['M_star'])),'results/E5/grid_optima.csv')
    lines=['% Generated by src/build_paper_assets.py; do not edit numerical values.',r'\newcommand{\val}[1]{\csname value@#1\endcsname}']
    for key,v in values.items():lines.append(r'\expandafter\def\csname value@'+key+r'\endcsname{'+v+'}')
    (OUT/'numbers.tex').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    (OUT/'numbers_sources.json').write_text(json.dumps(sources,indent=2),encoding='utf-8')
    if numbers_only:
        print(f'Generated {len(values)} numerical macros; bibliography preserved.')
        return
    # Retain local author spelling; fill and verify publication fields from DOI records.
    records=json.loads((ROOT/'qa/doi_verification.json').read_text(encoding='utf-8'))
    assert all(r['status']=='verified' for r in records),'Unverified DOI'
    source=(ROOT/'archive/20261002_antes_simple_paper/artigo/referencias.bib').read_text(encoding='utf-8')
    blocks=dict(re.findall(r'@article\{(\w+),(.+?)(?=\n@article|\Z)',source,re.S))
    def tex(s):
        s=html.unescape(re.sub('<[^>]+>','',str(s))).replace('\\&','&')
        for a,b in [('’',"'"),('–','--'),('—','---'),('ő',r'\H{o}'),('ž',r'\v{z}')]:s=s.replace(a,b)
        return s.replace('&',r'\&').replace('_',r'\_').replace('%',r'\%')
    bib=[];bibtex=[]
    paper=(OUT/'paper.tex').read_text(encoding='utf-8') if (OUT/'paper.tex').exists() else ''
    order=[]
    for group in re.findall(r'\\cite\{([^}]+)\}',paper):
        for key in group.split(','):
            if key not in order:order.append(key)
    records.sort(key=lambda r:order.index(r['key']) if r['key'] in order else 999)
    for r in records:
        k=r['key'];m=r['metadata'];fields=dict(re.findall(r'(\w+)\s*=\s*\{(.*?)\},?\s*\n',blocks[k],re.S))
        title=tex(m['title'][0]);authors=tex(fields['author']).replace(' and ', ', ')
        journal=tex(m['container-title'][0]);year=fields['year'];volume=m.get('volume',fields.get('volume',''))
        pages=m.get('page',m.get('article-number',fields.get('pages','')))
        if pages:pages=tex(pages).replace('--','-').replace('-','--')
        citation=authors+'. '+title.rstrip('.')+'. '+r'\textit{'+journal+'} '+r'\textbf{'+volume+'} ('+year+')'+(', '+pages if pages else '')+'. '+r'\href{https://doi.org/'+r['doi']+'}{doi:'+r['doi']+'}.'
        bib.append(r'\bibitem{'+k+'} '+citation)
        outfields={'author':fields['author'],'title':title,'journal':journal,'year':year,'volume':volume,'pages':pages,'doi':r['doi']}
        bibtex.append('@article{'+k+',\n'+',\n'.join('  '+f+' = {'+str(v)+'}' for f,v in outfields.items() if v)+'\n}')
    (OUT/'references.tex').write_text('\\begin{thebibliography}{99}\n'+'\n\n'.join(bib)+'\n\\end{thebibliography}\n',encoding='utf-8')
    (OUT/'referencias.bib').write_text('\n\n'.join(bibtex)+'\n',encoding='utf-8')
    print(f'Generated {len(values)} values and {len(records)} verified references.')
if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--numbers-only',action='store_true',help='Regenerate numerical macros without rebuilding the bibliography.')
    main(numbers_only=parser.parse_args().numbers_only)
