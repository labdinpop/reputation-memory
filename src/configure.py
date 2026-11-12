"""Create the requested experiment grids without launching simulations."""
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1]
BASE=dict(L=100,T0=5000,T1=2000,sample_every=10,R=30,w=.5,K=.1,eta=1e-6,
          epsilon=0.,update='synchronous',initial='exactly half, uniformly sampled; truthful repeated buffers',
          gap_definition='pooled over all site-time pairs in sampled coexistence rounds',
          identity=dict(observed_records=True,relative_tolerance=1e-9,EMA_tail_below=1e-12),
          stationarity=dict(half_difference=.02,fraction=.10,retry_T0=10000))
def main():
    out=ROOT/'config';out.mkdir(exist_ok=True)
    bs=[round(1+i*.025,3) for i in range(25)]
    designs={
        'E0':dict(identity_L=50,identity_epsilon=[0,.05],required=['weight normalization','truthful records at epsilon=0','W1=E1 pathwise','w=0 neutrality','G=sum(h*C_y), W/E/PP, epsilon=0 and .05, relative tolerance 1e-9 plus EMA tail','Bernoulli flips within three standard errors','Szabo-Toke automated benchmark']),
        'E1':dict(kernels=['W','E'],M=[1,2,3,5,7,9,13,17,25,33,49,65],b=bs,requires='E0 passed and measurement definitions resolved',
            b_selection=dict(reference='W9',A=dict(min_fc=.40,max_all_D=.10),B=dict(min_fc=.10,max_all_D=.10),minimum_separation=.05),
            temporary_E2_state_candidates=dict(kernel='W',M=[1,9,33],retain_selected_only=True)),
        'E2':dict(source='E1',kernel='W',M=[1,9,33],b=['b_A','b_B'],new_runs=0),
        'E3':dict(kernel='PP',a=[.05,.1,.15,.2,.3,.4,.6,.8,1],m=[1,2,4,8,16,32],deduplicate_a_one=True,b=['b_A','b_B']),
        'E4':dict(source=['E1','E3'],new_runs=0,figure_panels=['C_s','present and past contributions to G','f_C versus G']),
        'E5':dict(kernels=['W','E'],M=[1,2,3,5,9,17,33],epsilon=[0,.005,.01,.02,.05,.1],b=['b_A','b_B']),
        'E6':dict(kernel='W',M=[1,9,33],L=[50,100,200],b=['b_A','b_B'],epsilon=[0,.02],R_at_L200=20),
        'S1':dict(eta=.01,PP=dict(a=[.1,.2,.4],m=[2,8]),window_M=[9],b=['b_A','b_B']),
        'S2':dict(kernel='W',M=[1,9,33],b=bs,update='sequential')}
    for name,extra in designs.items():
        cfg=BASE|extra|dict(experiment=name)
        if name not in ['E0','E1']:cfg['requires']='explicit author approval of checkpoint 1'
        path=out/(name+'.json')
        if path.exists():raise FileExistsError(f'Preserve existing configuration: {path}')
        path.write_text(json.dumps(cfg,indent=2),encoding='utf-8')
    state=dict(gap_definition_approved=True,E0_passed=False,checkpoint1_approved=False,checkpoint2_approved=False)
    (out/'approvals.json').write_text(json.dumps(state,indent=2),encoding='utf-8')
    print('Created E0-E6, S1-S2 and approval gates; no simulation launched.')
if __name__=='__main__':main()
