"""Canonical Szabo–Toke regression, separate from the reputation model.

Four neighbors plus self-interaction, unnormalized payoff, K=.1, sampling focal
sites with replacement. This is not the paper's S2 permutation-sweep rule.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import csv,json,time
import numpy as np
from numba import njit
from model import neighbors

@njit(cache=True,nogil=True)
def canonical(L,b,sweeps,rng):
    n=L*L;s=np.zeros(n,np.uint8);s[:n//2]=1;rng.shuffle(s)
    nn=np.empty((n,4),np.int64)
    for i in range(n):
        y=i//L;x=i%L;nn[i,0]=((y-1)%L)*L+x;nn[i,1]=((y+1)%L)*L+x
        nn[i,2]=y*L+(x-1)%L;nn[i,3]=y*L+(x+1)%L
    c=np.empty(sweeps)
    for t in range(sweeps):
        for _ in range(n):
            i=rng.integers(0,n);j=nn[i,rng.integers(0,4)];ni=0;nj=0
            for k in range(4):ni+=s[nn[i,k]];nj+=s[nn[j,k]]
            fi=ni+1 if s[i] else b*ni;fj=nj+1 if s[j] else b*nj
            if rng.random()<1/(1+np.exp((fi-fj)/.1)):s[i]=s[j]
        c[t]=s.mean()
        if c[t]==0 or c[t]==1:
            c[t:]=c[t];break
    return c

def run_check():
    root=Path(__file__).resolve().parents[1];out=root/'qa';out.mkdir(exist_ok=True)
    canonical(3,1.4,1,np.random.default_rng(2))
    def job(task):
        b,rep=task;seed=202610020000+int(b*100)*10+rep;start=time.perf_counter()
        c=canonical(100,b,12000,np.random.default_rng(seed))
        return dict(L=100,b=b,realization=rep,seed=seed,tail_mean=float(c[6000:].mean()),final=float(c[-1]),seconds=time.perf_counter()-start)
    rows=list(ThreadPoolExecutor(4).map(job,[(b,r) for b in [1.4,1.86] for r in range(4)]))
    with (out/'benchmark_szabo.csv').open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    low=[r['tail_mean'] for r in rows if r['b']==1.4];high=[r['tail_mean'] for r in rows if r['b']==1.86]
    passed=all(x>0 for x in low) and all(r['final']==0 for r in rows if r['b']==1.86) and np.mean(low)>np.mean(high)
    status=dict(passed=bool(passed),test='Canonical-rule coarse-regime regression; not a new estimate of the critical point',
        mean_b14=float(np.mean(low)),mean_b186=float(np.mean(high)),new_realizations=8)
    (out/'benchmark_szabo.json').write_text(json.dumps(status,indent=2));print(json.dumps(status))
    return passed
if __name__=='__main__':
    if not run_check():raise SystemExit('Canonical benchmark acceptance failed; inspect before changing any tolerance.')
