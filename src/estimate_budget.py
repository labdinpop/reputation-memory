"""Timing-only runs authorized by the simulation budget; no scientific results reused."""
from pathlib import Path
import json,os
from model import Parameters,simulate
ROOT=Path(__file__).resolve().parents[1]
def main():
    # Warm compilation separately; timing estimates below exclude compilation.
    simulate(Parameters(L=4,b=1.1),burn=0,measure=20,every=5,early=False)
    profiles=[]
    for kernel in ['W','E']:
        for M in [9,65]:
            out=simulate(Parameters(L=100,b=1.1,kernel=kernel,M=M),early=False)
            profiles.append(dict(kernel=kernel,M=M,L=100,T0=5000,T1=2000,early_absorption_disabled=True,seconds=out['row']['elapsed_seconds']))
    w=max(r['seconds'] for r in profiles if r['kernel']=='W');e=max(r['seconds'] for r in profiles if r['kernel']=='E');workers=min(16,os.cpu_count() or 1)
    estimates={}
    for R in [30,20]:
        # Unique runs: E5 epsilon=0 and E6 L100 are available from E1/E5.
        sec={'E1':300*R*(w+e),'E2':0,'E3':98*R*w,'E4':0,'E5':70*R*1.25*(w+e),
            'E6':12*R*.25*w+12*20*4*w,'S1':14*R*w,'S2':75*R*2*w}
        estimates[str(R)]=dict(worker_hours={k:v/3600 for k,v in sec.items()},E1_E6_worker_hours=sum(sec[k] for k in ['E1','E3','E5','E6'])/3600,
            E1_E6_elapsed_hours_at_70pct_efficiency=sum(sec[k] for k in ['E1','E3','E5','E6'])/(3600*workers*.7),
            all_elapsed_hours_at_70pct_efficiency=sum(sec.values())/(3600*workers*.7))
    report=dict(profiles=profiles,workers=workers,estimates=estimates,
        assumptions='Worst profiled M per kernel; PP uses window upper timing; noise +25%; sequential x2; L scaling by site count. Independent stationarity retries excluded. Recalibrate before each experiment.',
        E1_release='requires successful E0 audit matching current source',
        timing_runs='Discarded from scientific datasets; no parameter-point estimates produced.')
    (ROOT/'qa/runtime_budget.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
if __name__=='__main__':main()
