"""Author-authorized, targeted T0=20000 refinement; no new parameter points."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import json,time
import numpy as np
from model import Parameters,simulate
from run_e1 import read_json,json_save,read_rows,save_rows,digest
from analyze_phase2 import load_summary

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/transient_refinement'

def prepare():
    OUT.mkdir(parents=True,exist_ok=True)
    if (OUT/'design.json').exists():return read_json(OUT/'design.json')
    selected=read_json(ROOT/'results/E1/b_selection.json');bs=[selected['b_A'],selected['b_B']];unique={}
    for exp in ['E1','E3','E5','E6']:
        for r in load_summary(ROOT/f'results/{exp}/point_summary.csv'):
            if not r['still_nonstationary'] or (exp=='E1' and r['b'] not in bs):continue
            source=Path(r['source']);key=str(source)
            if key in unique:unique[key]['used_by'].append(exp);continue
            raw=read_rows(ROOT/source);first=raw[0]
            p={k:(first[k] if k in ['kernel','update'] else int(float(first[k])) if k in ['L','M','m'] else float(first[k])) for k in Parameters.__dataclass_fields__ if k in first}
            if p.get('m')==0:p['m']=1
            unique[key]=dict(old_source=key,parameters=p,R=len(raw),used_by=[exp],old_meta=read_json((ROOT/source).parent/'point.json'))
    design=dict(T0=20000,previous_T0=10000,T1=2000,sample_every=10,points=list(unique.values()),authorization='User authorized rerunning only flagged points affecting b_A/b_B or E3-E6; no high-noise extension.')
    json_save(OUT/'design.json',design);return design

def point(task):
    source=ROOT/task['old_source'];folder=source.parent;path=folder/'T0_20000.csv';rows=read_rows(path);p=Parameters(**task['parameters'])
    for rep in range(len(rows),task['R']):
        rows.append(simulate(p,rep=rep,burn=20000)['row']);save_rows(path,rows)
        json_save(OUT/(folder.name+'.json'),dict(completed=len(rows),R=task['R'],source=str(path.relative_to(ROOT))))
    frac=float(np.mean([abs(float(r['half_difference']))>.02 for r in rows]));meta=dict(task['old_meta'])
    meta.update(selected_file=path.name,T0=20000,repeated=True,still_nonstationary=frac>.1,nonstationary_fraction=frac,refinement_authorized=True)
    return dict(task=task,new_source=str(path.relative_to(ROOT)),meta=meta)

def main():
    design=prepare();json_save(OUT/'status.json',dict(state='running',points=len(design['points'])))
    with ThreadPoolExecutor(6) as pool:completed=list(pool.map(point,design['points']))
    for done in completed:json_save((ROOT/done['new_source']).parent/'point.json',done['meta'])
    # Update references to reused realizations, preserving their previous provenance.
    changed={str(Path(d['task']['old_source'])):d for d in completed}
    for exp in ['E5','E6']:
        for mp in (ROOT/f'results/{exp}/points').glob('*/point.json'):
            m=read_json(mp)
            if m.get('reused') and str(Path(m['source_csv'])) in changed:
                d=changed[str(Path(m['source_csv']))];m['previous_source_csv']=m['source_csv'];m['source_csv']=d['new_source'];m['source_sha256']=digest(ROOT/d['new_source'])
                for k in ['T0','repeated','still_nonstationary']:m[k]=d['meta'][k]
                json_save(mp,m)
    json_save(OUT/'completed.json',completed)
    from analyze_e1 import analyze as first
    first()
    # Checkpoint 1 choices stay fixed: its reference W(9) is unaffected.
    from analyze_phase2 import analyze
    analyze()
    json_save(OUT/'status.json',dict(state='complete',points=len(completed),new_realizations=sum(x['task']['R'] for x in completed),still_flagged=sum(x['meta']['still_nonstationary'] for x in completed)))
if __name__=='__main__':main()
