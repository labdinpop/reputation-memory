"""Approved experiments E3/E5/E6/S1/S2. No manuscript drafting or power actions."""
import argparse,hashlib,itertools,json,os,time,traceback
from concurrent.futures import ThreadPoolExecutor,wait,FIRST_COMPLETED
from dataclasses import replace
from pathlib import Path
import numpy as np
from model import Parameters,simulate,seed_for
from run_e1 import read_json,json_save,read_rows,save_rows,digest,retry_io

ROOT=Path(__file__).resolve().parents[1]
CONTROL=ROOT/'results/phase2'
STOP=CONTROL/'PAUSE'
EXPERIMENTS=['E3','E5','E6','S1','S2']

def point_id(p):
    payload=json.dumps(p.canonical(),sort_keys=True,separators=(',',':'))
    return p.kernel+'_'+hashlib.sha256(payload.encode()).hexdigest()[:16]

def designs(name):
    c=read_json(ROOT/f'config/{name}.json');selection=read_json(ROOT/'results/E1/b_selection.json')
    bs=[selection.get(b,b) if isinstance(b,str) else b for b in c['b']]
    common={k:c[k] for k in ['w','K','eta','update']}
    ps=[]
    if name=='E3':
        for b,a,m in itertools.product(bs,c['a'],c['m']):
            if a==1 and m!=c['m'][0]:continue
            ps.append(Parameters(kernel='PP',L=c['L'],b=float(b),a=float(a),m=m,**common))
    elif name in ['E5','E6']:
        for b,k,M,eps,L in itertools.product(bs,c.get('kernels',[c.get('kernel')]),c['M'],c['epsilon'],c['L'] if isinstance(c['L'],list) else [c['L']]):
            ps.append(Parameters(kernel=k,M=M,L=L,b=float(b),epsilon=float(eps),**common))
    elif name=='S1':
        for b in bs:
            for a,m in itertools.product(c['PP']['a'],c['PP']['m']):ps.append(Parameters(kernel='PP',a=float(a),m=m,b=float(b),L=c['L'],**common))
            for M in c['window_M']:ps.append(Parameters(kernel='W',M=M,b=float(b),L=c['L'],**common))
    else:
        for b,M in itertools.product(bs,c['M']):ps.append(Parameters(kernel='W',M=M,b=float(b),L=c['L'],**common))
    tasks=[dict(name=point_id(p),parameters=p.canonical(),R=c.get('R_at_L200',c['R']) if p.L==200 else c['R']) for p in ps]
    assert len({t['name'] for t in tasks})==len(tasks)
    return c,tasks

def parameters(task):
    p=dict(task['parameters'])
    if p.get('m')==0:p['m']=1
    return Parameters(**p)

def reuse_folder(name,p):
    if name in ['E5','E6'] and p.L==100 and p.epsilon==0:
        return ROOT/'results/E1/points'/f'{p.kernel}_M{p.M:02d}_b{p.b:.3f}'
    if name=='E6' and p.L==100:
        return ROOT/'results/E5/points'/point_id(p)
    return None

def validate_rows(rows,p,R,c):
    if len(rows)>R:raise RuntimeError('Too many saved realizations')
    for rep,r in enumerate(rows):
        if int(r['realization'])!=rep or int(r['seed'])!=seed_for(p,rep):raise RuntimeError('Saved seed or realization mismatch')
        if int(r['T1'])!=c['T1'] or int(r['sample_every'])!=c['sample_every']:raise RuntimeError('Saved measurement protocol mismatch')

def run_point(experiment,c,task):
    p=parameters(task);R=task['R'];folder=ROOT/f'results/{experiment}/points'/task['name']
    folder.mkdir(parents=True,exist_ok=True);mp=folder/'point.json'
    if mp.exists():
        prior=read_json(mp)
        if prior['state']=='complete':return prior
    source=reuse_folder(experiment,p)
    if source is not None:
        original=read_json(source/'point.json')
        if original['state']!='complete':raise RuntimeError('Reuse source incomplete')
        csv=source/original['selected_file'];rows=read_rows(csv);validate_rows(rows,p,R,c)
        if len(rows)!=R:raise RuntimeError('Reuse realization count mismatch')
        meta=dict(state='complete',point=task['name'],parameters=p.canonical(),R=R,completed=R,reused=True,
            source_csv=str(csv.relative_to(ROOT)),source_sha256=digest(csv),T0=original['T0'],
            repeated=original['repeated'],still_nonstationary=original['still_nonstationary'])
        json_save(mp,meta);return meta
    for burn in [c['T0'],c['stationarity']['retry_T0']]:
        csv=folder/f'T0_{burn}.csv';rows=read_rows(csv);validate_rows(rows,p,R,c)
        for rep in range(len(rows),R):
            if STOP.exists():return dict(state='paused',point=task['name'])
            row=simulate(p,rep=rep,burn=burn,measure=c['T1'],every=c['sample_every'])['row']
            rows.append(row);save_rows(csv,rows)
            json_save(mp,dict(state='running',point=task['name'],parameters=p.canonical(),R=R,completed=len(rows),T0=burn,reused=False))
        fraction=float(np.mean([abs(float(r['half_difference']))>c['stationarity']['half_difference'] for r in rows]))
        unstable=fraction>c['stationarity']['fraction']
        if unstable and burn==c['T0']:
            json_save(folder/'stationarity_retry.json',dict(first_T0=burn,retry_T0=c['stationarity']['retry_T0'],nonstationary_fraction=fraction));continue
        meta=dict(state='complete',point=task['name'],parameters=p.canonical(),R=R,completed=R,reused=False,
            selected_file=csv.name,T0=burn,repeated=burn!=c['T0'],still_nonstationary=unstable,nonstationary_fraction=fraction)
        json_save(mp,meta);return meta

def profile(workers=16):
    CONTROL.mkdir(parents=True,exist_ok=True)
    b=read_json(ROOT/'results/E1/b_selection.json')['b_A']
    probes={'E3':[Parameters(kernel='PP',a=.05,m=32,b=b)],
        'E5':[Parameters(kernel=k,M=33,epsilon=.1,b=b) for k in ['W','E']],
        'E6':[Parameters(M=33,epsilon=.02,b=b)],
        'S1':[Parameters(kernel='PP',a=.1,m=8,eta=.01,b=b),Parameters(M=9,eta=.01,b=b)],
        'S2':[Parameters(M=33,b=b,update='sequential')]}
    simulate(Parameters(L=4),burn=0,measure=20,every=5,early=False)
    report=dict(workers=workers,model_sha256=digest(ROOT/'src/model.py'),experiments={},
        assumptions='Full protocol at L=100; early absorption disabled; cost scales with site count; 70% parallel efficiency. No stationarity retries or I/O included; timing outputs excluded from scientific data.')
    for name,params in probes.items():
        records=[]
        for p in params:
            row=simulate(p,rep=99999,early=False)['row']
            records.append(dict(parameters=p.canonical(),seconds=row['elapsed_seconds'],coexistence_samples=row['coexistence_samples']))
        c,tasks=designs(name);cost=0;new_runs=0;reuse=0
        for t in tasks:
            p=parameters(t)
            if reuse_folder(name,p) is not None:reuse+=t['R'];continue
            seconds=max(r['seconds'] for r in records if r['parameters']['kernel']==p.kernel)
            cost+=seconds*(p.L/100)**2*t['R'];new_runs+=t['R']
        entry=dict(profiles=records,new_realizations=new_runs,reused_realizations=reuse,points=len(tasks),
            worker_hours=cost/3600,estimated_elapsed_minutes=cost/(workers*.7*60))
        report['experiments'][name]=entry
        print(json.dumps(dict(experiment=name,**entry)),flush=True)
    report['total_worker_hours']=sum(e['worker_hours'] for e in report['experiments'].values())
    report['total_estimated_elapsed_minutes']=sum(e['estimated_elapsed_minutes'] for e in report['experiments'].values())
    prior=read_json(ROOT/'qa/runtime_budget.json')['estimates']['30']['worker_hours']['E1']
    report['E1_E6_worker_hours']=prior+sum(report['experiments'][n]['worker_hours'] for n in ['E3','E5','E6'])
    json_save(ROOT/'qa/runtime_phase2.json',report)
    return report

def snapshot():
    state=read_json(CONTROL/'status.json') if (CONTROL/'status.json').exists() else dict(state='not_started')
    experiments={}
    for name in EXPERIMENTS:
        c,tasks=designs(name);entry=dict(total_points=len(tasks),completed_points=0,new_realizations_saved=0,reused_realizations=0,repeated_points=0)
        for mp in (ROOT/f'results/{name}/points').glob('*/point.json'):
            m=read_json(mp)
            if m['state']=='complete':entry['completed_points']+=1
            if m.get('reused'):entry['reused_realizations']+=m['completed']
            else:
                entry['new_realizations_saved']+=m['completed']
                if m['T0']>c['T0']:entry['new_realizations_saved']+=m['R'];entry['repeated_points']+=1
        experiments[name]=entry
    state['experiments']=experiments
    state['new_realizations_saved']=sum(e['new_realizations_saved'] for e in experiments.values())
    return state

def progress(state):
    try:
        snap=snapshot();state['experiments']=snap['experiments'];state['new_realizations_saved']=snap['new_realizations_saved'];state.pop('progress_warning',None)
    except OSError as exc:state['progress_warning']=str(exc)
    state['updated_at']=time.strftime('%Y-%m-%dT%H:%M:%S')
    try:json_save(CONTROL/'status.json',state)
    except OSError as exc:print(f'Progress save delayed: {exc}',flush=True)

def run(workers):
    gates=read_json(ROOT/'config/approvals.json')
    if not gates['checkpoint1_approved'] or not gates['E0_passed']:raise RuntimeError('Author checkpoint 1 approval or E0 missing')
    audit=read_json(ROOT/'qa/E0_report.json')
    for name,sha in audit['source_sha256'].items():
        if digest(ROOT/name)!=sha:raise RuntimeError('E0 source changed')
    budget=read_json(ROOT/'qa/runtime_phase2.json')
    if budget['model_sha256']!=digest(ROOT/'src/model.py'):raise RuntimeError('Timing source changed')
    if budget['E1_E6_worker_hours']>48 and any(designs(n)[0]['R']>20 for n in EXPERIMENTS):raise RuntimeError('Budget exceeds 48h: reduce R to 20 before launch')
    if STOP.exists():raise RuntimeError('Pause requested; use resume explicitly')
    CONTROL.mkdir(parents=True,exist_ok=True)
    lock=CONTROL/'run.lock';fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY);os.write(fd,str(os.getpid()).encode());os.close(fd)
    state=dict(state='running',pid=os.getpid(),workers=workers,stop_after='checkpoint 2',started_at=time.strftime('%Y-%m-%dT%H:%M:%S'))
    try:
        sources=['src/model.py','src/run_e1.py','src/run_phase2.py']+[f'config/{n}.json' for n in EXPERIMENTS]+['results/E1/b_selection.json']
        manifest=dict(source_sha256={n:digest(ROOT/n) for n in sources})
        mp=CONTROL/'manifest.json'
        if mp.exists() and read_json(mp)!=manifest:raise RuntimeError('Phase 2 resume manifest mismatch')
        json_save(mp,manifest);json_save(CONTROL/'status.json',state)
        for name in EXPERIMENTS:
            state['experiment']=name;progress(state)
            print(json.dumps(dict(event='starting_experiment',experiment=name,budget=budget['experiments'][name])),flush=True)
            c,tasks=designs(name);folder=ROOT/f'results/{name}';folder.mkdir(exist_ok=True)
            json_save(folder/'design.json',dict(configuration=c,points=tasks,b_selection=read_json(ROOT/'results/E1/b_selection.json')))
            iterator=iter(tasks)
            with ThreadPoolExecutor(max_workers=workers) as pool:
                pending={pool.submit(run_point,name,c,next(iterator)) for _ in range(min(workers,len(tasks)))}
                while pending:
                    done,pending=wait(pending,timeout=10,return_when=FIRST_COMPLETED)
                    for future in done:
                        print(json.dumps(dict(experiment=name,**future.result())),flush=True)
                        if not STOP.exists():
                            task=next(iterator,None)
                            if task is not None:pending.add(pool.submit(run_point,name,c,task))
                    progress(state)
            if STOP.exists():state['state']='paused';progress(state);return
        state['state']='simulations_complete';state['experiment']=None;progress(state)
        # Analysis has its own source record; it cannot change the simulation design.
        from analyze_phase2 import analyze
        state['state']='analyzing';progress(state);analyze()
        state['state']='awaiting_checkpoint2_review';state['finished_at']=time.strftime('%Y-%m-%dT%H:%M:%S');progress(state)
    except BaseException:
        state['state']='failed';state['error']=traceback.format_exc();json_save(CONTROL/'status.json',state);raise
    finally:retry_io(lambda:lock.unlink(missing_ok=True))

def main():
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['profile','run','resume','status','pause']);parser.add_argument('--workers',type=int,default=16)
    args=parser.parse_args()
    if args.workers<1:raise ValueError('workers must be positive')
    if args.command=='profile':profile(args.workers);return
    if args.command=='status':print(json.dumps(snapshot(),indent=2));return
    if args.command=='pause':CONTROL.mkdir(parents=True,exist_ok=True);STOP.touch();print('Pause requested; current realizations will finish and save.');return
    if args.command=='resume':STOP.unlink(missing_ok=True)
    run(args.workers)
if __name__=='__main__':main()
