"""Resumable E1 only. Completion produces checkpoint 1 and never starts phase 2."""
import argparse,csv,hashlib,json,os,threading,time,traceback
from concurrent.futures import ThreadPoolExecutor,wait,FIRST_COMPLETED
from pathlib import Path
import numpy as np
from model import Parameters,simulate

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/E1'
STOP=OUT/'PAUSE'

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def retry_io(operation,attempts=30):
    """Retry transient Windows access/sharing locks, on reads as well as writes."""
    for attempt in range(attempts):
        try:
            return operation()
        except OSError as exc:
            transient=isinstance(exc,PermissionError) or getattr(exc,'winerror',None) in (5,32,33)
            if not transient or attempt==attempts-1:raise
            time.sleep(min(.1*2**min(attempt,5),2.))

def atomic_replace(tmp,path,attempts=30):
    """Keep the old checkpoint intact while Windows temporarily locks its file."""
    return retry_io(lambda:tmp.replace(path),attempts)

def read_json(path):
    return json.loads(retry_io(lambda:path.read_text(encoding='utf-8')))

def json_save(path,obj):
    tmp=path.with_suffix('.tmp')
    retry_io(lambda:tmp.write_text(json.dumps(obj,indent=2),encoding='utf-8'));atomic_replace(tmp,path)

def read_rows(path):
    if not path.exists():return []
    def read():
        with path.open(newline='',encoding='utf-8') as f:return list(csv.DictReader(f))
    return retry_io(read)

def save_rows(path,rows):
    tmp=path.with_suffix('.tmp')
    def write():
        with tmp.open('w',newline='',encoding='utf-8') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    retry_io(write)
    atomic_replace(tmp,path)

def points(cfg):
    return [(f'{k}_M{M:02d}_b{b:.3f}',Parameters(kernel=k,M=M,b=b,
        **{key:cfg[key] for key in ['L','epsilon','w','K','eta','update']}))
        for M in cfg['M'] for b in cfg['b'] for k in cfg['kernels']]

def summarize(rows):
    vals=np.array([float(r['half_difference']) for r in rows])
    return dict(R=len(rows),nonstationary_fraction=float(np.mean(abs(vals)>.02)))

def run_point(task,cfg):
    name,p=task;folder=OUT/'points'/name;folder.mkdir(parents=True,exist_ok=True)
    meta_path=folder/'point.json'
    if meta_path.exists():
        meta=read_json(meta_path)
        if meta['state']=='complete':return meta
    candidate=p.kernel=='W' and p.M in cfg['temporary_E2_state_candidates']['M']
    for burn in [cfg['T0'],cfg['stationarity']['retry_T0']]:
        path=folder/f'T0_{burn}.csv';rows=read_rows(path)
        if len(rows)>cfg['R']:raise RuntimeError('Unexpected realization count')
        for rep in range(len(rows),cfg['R']):
            if STOP.exists():return dict(state='paused',point=name)
            output=simulate(p,rep=rep,burn=burn,measure=cfg['T1'],every=cfg['sample_every'])
            row=output['row']
            if candidate:
                target=OUT/'E2_candidates'/f'{name}_T0{burn}_r{rep:02d}.npy'
                tmp=target.with_suffix('.tmp')
                def write_state():
                    with tmp.open('wb') as f:np.save(f,output['final'].reshape(p.L,p.L))
                retry_io(write_state)
                atomic_replace(tmp,target)
            rows.append(row);save_rows(path,rows)
            json_save(meta_path,dict(state='running',point=name,T0=burn,completed=len(rows),R=cfg['R']))
        stats=summarize(rows)
        retry=stats['nonstationary_fraction']>cfg['stationarity']['fraction']
        if burn==cfg['T0'] and retry:
            json_save(folder/'stationarity_retry.json',dict(reason='More than 10% of realizations differ by over .02 between measurement halves',
                first_T0=burn,retry_T0=cfg['stationarity']['retry_T0'],**stats))
            continue
        meta=dict(state='complete',point=name,selected_file=path.name,T0=burn,completed=len(rows),
            repeated=burn!=cfg['T0'],still_nonstationary=retry,**stats)
        json_save(meta_path,meta);return meta
    raise RuntimeError('Unreachable point state')

def snapshot():
    status=read_json(OUT/'status.json') if (OUT/'status.json').exists() else dict(state='not_started')
    completed=0;done=0;retry=0;active=[]
    for path in (OUT/'points').glob('*/point.json'):
        m=read_json(path);completed+=m['completed']
        if m['state']=='complete':done+=1
        else:active.append(m['point'])
        if m.get('T0',5000)>5000:retry+=1;completed+=status.get('R',30)
    status.update(completed_realizations=completed,completed_points=done,repeated_points=retry,active_points=active)
    return status

def update_progress(state):
    # Progress is advisory; failure to read it must not abort numerical workers.
    try:
        state.update(snapshot());state.pop('progress_warning',None)
        state['progress_stale']=False
    except OSError as exc:
        state['progress_stale']=True;state['progress_warning']=str(exc)
    state['updated_at']=time.strftime('%Y-%m-%dT%H:%M:%S')
    try:json_save(OUT/'status.json',state)
    except OSError as exc:print(f'Progress save delayed: {exc}',flush=True)

def run(workers):
    cfg=read_json(ROOT/'config/E1.json');gates=read_json(ROOT/'config/approvals.json')
    audit=read_json(ROOT/'qa/E0_report.json')
    if not gates['E0_passed'] or not gates['gap_definition_approved'] or not audit['E0_accepted']:raise RuntimeError('E0 gate closed')
    for rel,sha in audit['source_sha256'].items():
        if digest(ROOT/rel)!=sha:raise RuntimeError(f'E0 source changed: {rel}')
    if STOP.exists():raise RuntimeError('Pause requested. Use resume explicitly.')
    OUT.mkdir(parents=True,exist_ok=True);(OUT/'points').mkdir(exist_ok=True);(OUT/'E2_candidates').mkdir(exist_ok=True)
    lock=OUT/'run.lock'
    try:
        fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY)
    except FileExistsError:raise RuntimeError('Runner lock exists; inspect its PID before removing a stale lock')
    os.write(fd,str(os.getpid()).encode());os.close(fd)
    try:
        sources=['src/model.py','src/run_e1.py','src/analyze_e1.py','config/E1.json']
        manifest=dict(source_sha256={p:digest(ROOT/p) for p in sources},configuration=cfg)
        mp=OUT/'manifest.json'
        if mp.exists() and read_json(mp)!=manifest:raise RuntimeError('Resume manifest mismatch')
        json_save(mp,manifest)
        tasks=points(cfg);state=dict(state='running',pid=os.getpid(),workers=workers,R=cfg['R'],total_points=len(tasks),
            base_realizations=len(tasks)*cfg['R'],started_at=time.strftime('%Y-%m-%dT%H:%M:%S'),stop_after='checkpoint 1')
        json_save(OUT/'status.json',state)
        # Compile before entering the pool. This discarded tiny run is not part of E1.
        simulate(Parameters(L=4),burn=0,measure=20,every=5)
        iterator=iter(tasks)
        with ThreadPoolExecutor(max_workers=workers) as pool:
            pending={pool.submit(run_point,next(iterator),cfg) for _ in range(min(workers,len(tasks)))}
            while pending:
                finished,pending=wait(pending,timeout=10,return_when=FIRST_COMPLETED)
                for future in finished:
                    meta=future.result();print(json.dumps(meta),flush=True)
                    if not STOP.exists():
                        task=next(iterator,None)
                        if task is not None:pending.add(pool.submit(run_point,task,cfg))
                update_progress(state)
        if STOP.exists():state['state']='paused'
        else:
            state['state']='analyzing';json_save(OUT/'status.json',state)
            from analyze_e1 import analyze
            analyze();state['state']='awaiting_checkpoint1_approval'
        state['finished_at']=time.strftime('%Y-%m-%dT%H:%M:%S');json_save(OUT/'status.json',state)
    except BaseException:
        json_save(OUT/'status.json',dict(state='failed',pid=os.getpid(),error=traceback.format_exc()))
        raise
    finally:lock.unlink(missing_ok=True)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['run','resume','status','pause']);parser.add_argument('--workers',type=int,default=16)
    args=parser.parse_args()
    if args.command=='status':print(json.dumps(snapshot(),indent=2));return
    if args.command=='pause':OUT.mkdir(parents=True,exist_ok=True);STOP.touch();print('Pause requested; running realizations will finish and save.');return
    if args.workers<1:raise ValueError('workers must be positive')
    if args.command=='resume':STOP.unlink(missing_ok=True)
    run(args.workers)
if __name__=='__main__':main()
