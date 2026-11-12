"""Approved 180-run past-only S3 extension. Never writes manuscript artifacts."""
from pathlib import Path
import argparse, json, os, time, traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np
from model_s3_past import Parameters, simulate, weights, seed_for
from run_e1 import read_json, json_save, read_rows, save_rows, digest

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/S3_past'
QA = ROOT / 'qa/S3_past'
CFG = read_json(ROOT / 'config/S3_past.json')


def parameter(m, L=None):
    return Parameters(kernel='PP', a=0., m=m, L=L or CFG['L'], b=CFG['b'],
                      epsilon=0., w=1., eta=CFG['eta'], K=CFG['K'])


def profile():
    QA.mkdir(parents=True, exist_ok=True)
    checks = []
    for m in CFG['m']:
        p = parameter(m, L=6)
        h, tail = weights(p)
        assert h[0] == 0 and np.allclose(h[1:], 1/m) and abs(h.sum()-1) < 1e-14
        # Reconstruct reputations independently from true pre-imitation states.
        result = simulate(p, rep=987, burn=0, measure=60, every=1, early=False, save_trace=True)
        trace = result['trace'].astype(float)
        valid = (trace.sum(axis=1)>0) & (trace.sum(axis=1)<trace.shape[1])
        direct = np.array([trace[[max(0,t-lag) for lag in range(1,m+1)]].mean(axis=0)
                           for t in range(len(trace))])
        s = trace[valid].ravel(); r = direct[valid].ravel()
        assert len(s)>0
        expected = r[s==1].mean()-r[s==0].mean()
        assert abs(expected-result['row']['G']) < 1e-12
        assert np.allclose(result['row']['G_from_Cy'], expected, rtol=1e-9, atol=1e-12)
        checks.append(dict(m=m, passed=True, direct_gap=expected))
    # The adapter leaves existing S3 trajectories and seeds unchanged.
    from model_s3 import Parameters as Base, simulate as base_simulate
    base = base_simulate(Base(L=6,M=9,w=1), rep=321,burn=0,measure=60,every=1,save_trace=True)
    new = simulate(Parameters(L=6,M=9,w=1),rep=321,burn=0,measure=60,every=1,save_trace=True)
    assert np.array_equal(base['trace'],new['trace']) and base['row']['seed']==new['row']['seed']
    json_save(QA/'validation.json',dict(checks=checks, existing_S3_unchanged=True))
    # Full protocol, excluding JIT and timing run from scientific results.
    row=simulate(parameter(max(CFG['m'])),rep=99999,burn=CFG['T0'],
                 measure=CFG['T1'],every=CFG['sample_every'],early=False)['row']
    seconds=row['elapsed_seconds']; total=len(CFG['m'])*CFG['R']
    report=dict(seconds_per_realization=seconds, total_realizations=total,
                worker_hours=total*seconds/3600,
                estimated_minutes=total*seconds/(CFG['workers']*.7*60),
                source_sha256={x:digest(ROOT/x) for x in ['src/model_s3.py','src/model_s3_past.py','src/run_s3_past.py','config/S3_past.json']},
                notes='Full L=100,m=20 protocol; estimate excludes possible T0=10000 retries.')
    json_save(QA/'runtime.json',report)
    print(json.dumps(report,indent=2),flush=True)


def point(m):
    p=parameter(m); folder=OUT/'points'/f'PP_a0_m{m:02d}';folder.mkdir(parents=True,exist_ok=True)
    mp=folder/'point.json'
    if mp.exists() and read_json(mp)['state']=='complete':return read_json(mp)
    for burn in [CFG['T0'],CFG['stationarity']['retry_T0']]:
        path=folder/f'T0_{burn}.csv';rows=read_rows(path)
        for rep,r in enumerate(rows):
            assert int(r['realization'])==rep and int(r['seed'])==seed_for(p,rep)
        for rep in range(len(rows),CFG['R']):
            if (OUT/'PAUSE').exists():return dict(state='paused',m=m,completed=len(rows))
            rows.append(simulate(p,rep=rep,burn=burn,measure=CFG['T1'],every=CFG['sample_every'])['row'])
            save_rows(path,rows)
            json_save(mp,dict(state='running',m=m,T0=burn,completed=len(rows)))
        fraction=float(np.mean([abs(float(r['half_difference']))>CFG['stationarity']['half_difference'] for r in rows]))
        unstable=fraction>CFG['stationarity']['fraction']
        if unstable and burn==CFG['T0']:
            json_save(folder/'stationarity_retry.json',dict(first_T0=burn,retry_T0=CFG['stationarity']['retry_T0'],fraction=fraction))
            continue
        meta=dict(state='complete',m=m,T0=burn,R=CFG['R'],completed=len(rows),
                  selected_file=path.name,source=str(path.relative_to(ROOT)),
                  repeated=burn!=CFG['T0'],still_nonstationary=unstable,nonstationary_fraction=fraction)
        json_save(mp,meta);return meta


def status():
    points=[read_json(p) for p in sorted((OUT/'points').glob('*/point.json'))]
    return dict(points=points,completed=sum(p['completed'] for p in points),
                total=len(CFG['m'])*CFG['R'],complete_points=sum(p['state']=='complete' for p in points))


def run():
    assert CFG['author_approved'] and read_json(QA/'validation.json')['existing_S3_unchanged']
    hashes=read_json(QA/'runtime.json')['source_sha256']
    assert all(digest(ROOT/k)==v for k,v in hashes.items())
    OUT.mkdir(parents=True,exist_ok=True)
    protected={str(p.relative_to(ROOT)):digest(p) for folder in ['artigo','results/figures']
               for p in (ROOT/folder).rglob('*') if p.is_file()}
    json_save(OUT/'manifest.json',dict(source_sha256=hashes,protected_sha256=protected,config=CFG))
    state=dict(state='running',pid=os.getpid(),workers=CFG['workers'])
    json_save(OUT/'status.json',state|status())
    try:
        with ThreadPoolExecutor(max_workers=CFG['workers']) as pool:
            for future in as_completed([pool.submit(point,m) for m in CFG['m']]):
                print(json.dumps(future.result()),flush=True)
                json_save(OUT/'status.json',state|status())
        end=status();end['state']='complete' if end['complete_points']==len(CFG['m']) else 'paused'
        assert all(digest(ROOT/k)==v for k,v in protected.items())
        json_save(OUT/'status.json',state|end|dict(manuscript_unchanged=True))
    except BaseException:
        json_save(OUT/'status.json',state|dict(state='failed',error=traceback.format_exc()));raise


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['profile','run','status'])
    a=parser.parse_args()
    if a.command=='profile':profile()
    elif a.command=='run':run()
    else:print(json.dumps(status(),indent=2))
