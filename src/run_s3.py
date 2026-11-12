"""Approved S3 selection-only comparison; resumable and isolated from E1--S2."""
from pathlib import Path
import argparse,json,time,os,traceback
from concurrent.futures import ThreadPoolExecutor,as_completed
import numpy as np
from model_s3 import Parameters,simulate,_reference,neighbors,seed_for
from run_e1 import read_json,json_save,read_rows,save_rows,digest
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/S3'
CFG=read_json(ROOT/'config/S3.json')
def parameter(M):return Parameters(M=M,L=100,b=CFG['b'][0],w=1,eta=0)
def profile():
    nn=neighbors(4);r=np.zeros(16)
    assert [_reference(r,nn,0,u,1,0) for u in [.1,.3,.6,.9]]==list(nn[0])
    r[nn[0]]=[0,1,3,0]
    assert _reference(r,nn,0,.24,1,0)==nn[0,1]
    assert _reference(r,nn,0,.26,1,0)==nn[0,2]
    simulate(Parameters(L=4,w=1,eta=0),burn=0,measure=20,every=5)
    row=simulate(parameter(33),rep=99999,early=False)['row']
    seconds=row['elapsed_seconds'];report=dict(seconds_per_realization=seconds,total_realizations=210,worker_hours=210*seconds/3600,estimated_minutes_7_workers=210*seconds/(7*.7*60),source_sha256=digest(ROOT/'src/model_s3.py'),notes='Full L=100 protocol, M=33, excludes JIT; estimate excludes any stationarity retries; timing realization excluded from results.')
    json_save(ROOT/'qa/S3/runtime.json',report);print(json.dumps(report),flush=True)
def point(M):
    p=parameter(M);folder=OUT/'points'/f'W_M{M:02d}';folder.mkdir(parents=True,exist_ok=True)
    mp=folder/'point.json'
    if mp.exists() and read_json(mp)['state']=='complete':return read_json(mp)
    for burn in [CFG['T0'],CFG['stationarity']['retry_T0']]:
        path=folder/f'T0_{burn}.csv';rows=read_rows(path)
        for rep,r in enumerate(rows):assert int(r['realization'])==rep and int(r['seed'])==seed_for(p,rep)
        for rep in range(len(rows),CFG['R']):
            if (OUT/'PAUSE').exists():return dict(state='paused',M=M)
            rows.append(simulate(p,rep=rep,burn=burn,measure=CFG['T1'],every=CFG['sample_every'])['row']);save_rows(path,rows)
            json_save(mp,dict(state='running',M=M,T0=burn,completed=len(rows)))
        fraction=float(np.mean([abs(float(r['half_difference']))>.02 for r in rows]));unstable=fraction>.1
        if unstable and burn==CFG['T0']:
            json_save(folder/'stationarity_retry.json',dict(first_T0=burn,retry_T0=10000,nonstationary_fraction=fraction));continue
        meta=dict(state='complete',M=M,T0=burn,R=CFG['R'],completed=len(rows),selected_file=path.name,source=str(path.relative_to(ROOT)),repeated=burn!=CFG['T0'],still_nonstationary=unstable,nonstationary_fraction=fraction)
        json_save(mp,meta);return meta
def status():
    points=[read_json(p) for p in (OUT/'points').glob('*/point.json')]
    return dict(points=points,completed=sum(p['completed'] for p in points),total=210,complete_points=sum(p['state']=='complete' for p in points))
def run():
    assert CFG['author_approved'];assert read_json(ROOT/'qa/S3/runtime.json')['source_sha256']==digest(ROOT/'src/model_s3.py')
    state=dict(state='running',pid=os.getpid(),workers=7)
    json_save(OUT/'status.json',state)
    json_save(OUT/'manifest.json',dict(source_sha256={p:digest(ROOT/p) for p in ['src/model.py','src/model_s3.py','src/run_s3.py','config/S3.json']}))
    try:
        with ThreadPoolExecutor(max_workers=7) as pool:
            for future in as_completed([pool.submit(point,M) for M in CFG['M']]):
                print(json.dumps(future.result()),flush=True);json_save(OUT/'status.json',state|status())
        json_save(OUT/'status.json',state|status()|dict(state='complete' if status()['complete_points']==7 else 'paused'))
    except BaseException:
        json_save(OUT/'status.json',state|dict(state='failed',error=traceback.format_exc()));raise
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['profile','run','status']);a=parser.parse_args()
    if a.command=='profile':profile()
    elif a.command=='run':run()
    else:print(json.dumps(status(),indent=2))
