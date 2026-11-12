"""The simple paper model: W, E and present–past memory only.

Numerical loops reuse the project's payoff and reference-selection conventions.
No archived results or simulation modules are imported.
"""
from dataclasses import dataclass,asdict
import hashlib,json,math,time
import numpy as np
from numba import njit

@dataclass(frozen=True)
class Parameters:
    kernel:str='W'
    L:int=100
    b:float=1.3
    M:int=9
    a:float=.2
    m:int=8
    epsilon:float=0.
    w:float=.5
    K:float=.1
    eta:float=1e-6
    update:str='synchronous'
    def __post_init__(self):
        if self.kernel not in ['W','E','PP'] or self.L<3 or self.b<1:raise ValueError('Invalid kernel, lattice or b')
        if self.M<1 or self.m<1 or not 0<self.a<=1:raise ValueError('Invalid memory')
        if not 0<=self.epsilon<=1 or not 0<=self.w<=1 or self.K<=0 or self.eta<=0:raise ValueError('Invalid selection or observation parameter')
        if self.update not in ['synchronous','sequential']:raise ValueError('Unknown update')
    def canonical(self):
        d=asdict(self)
        if self.kernel in ['W','E']:d.pop('a');d.pop('m')
        else:
            d.pop('M')
            if self.a==1:d['m']=0
        return d

def seed_for(p,rep):
    payload=json.dumps(dict(parameters=p.canonical(),realization=int(rep)),sort_keys=True,separators=(',',':'))
    return int.from_bytes(hashlib.blake2b(payload.encode(),digest_size=8,person=b'simple-paper-v1').digest(),'little')

def weights(p):
    if p.kernel=='W':return np.ones(p.M)/p.M,0.
    if p.kernel=='PP':return (np.array([1.]) if p.a==1 else np.r_[p.a,np.full(p.m,(1-p.a)/p.m)]),0.
    alpha=2/(p.M+1)
    H=1 if alpha==1 else math.floor(math.log(1e-12)/math.log1p(-alpha))+1
    h=alpha*(1-alpha)**np.arange(H)
    return h,float((1-alpha)**H)

def neighbors(L):
    grid=np.arange(L*L).reshape(L,L)
    return np.stack([np.roll(grid,1,0),np.roll(grid,-1,0),np.roll(grid,1,1),np.roll(grid,-1,1)],axis=-1).reshape(-1,4)

@njit(cache=True,nogil=True)
def _gain(s,nn,i,b):
    total=0
    for k in range(4):total+=s[nn[i,k]]
    return total*(1. if s[i] else b)/(4*b)

@njit(cache=True,nogil=True)
def _reference(r,nn,i,u,w,eta):
    total=0.
    for k in range(4):total+=r[nn[i,k]]+eta
    cumulative=0.
    for k in range(4):
        j=nn[i,k];cumulative+=(1-w)/4+w*(r[j]+eta)/total
        if u<cumulative:return j
    return nn[i,3]

@njit(cache=True,nogil=True)
def _accept(gi,gj,K):
    z=(gi-gj)/K
    if z>=0:
        ez=math.exp(-z);return ez/(1+ez)
    return 1/(1+math.exp(z))

@njit(cache=True,nogil=True)
def _run(s,nn,h,kernel,alpha,a,b,epsilon,w,K,eta,sequential,burn,measure,every,irng,erng,early,save_trace):
    n=len(s);H=len(h);th=np.empty((H,n),np.uint8);oh=np.empty((H,n),np.uint8)
    for lag in range(H):th[lag]=s;oh[lag]=s
    osum=s.astype(np.float64)*H;r=s.astype(np.float64)
    ns=measure//every;fc=np.empty(ns);gtime=np.full(ns,np.nan);record=0
    lag_sum=np.zeros(H);cross=np.zeros(H);obs_sum=np.zeros(H);obs_cross=np.zeros(H)
    sx=0.;sr=0.;sxr=0.;sample_count=0.;flips=0
    total_steps=burn+measure
    trace=np.empty((total_steps if save_trace else 0,n),np.uint8)
    F=np.empty(n);out=np.empty_like(s);order=np.arange(n)
    for t in range(total_steps):
        nc=int(s.sum())
        if early and epsilon==0 and t<=burn and (nc==0 or nc==n) and not save_trace:
            fc[:]=nc/n
            return s,fc,gtime,lag_sum,cross,obs_sum,obs_cross,0.,0.,0.,0.,0,trace,True
        pos=t%H
        if epsilon>0:draws=erng.random(n)
        for i in range(n):
            y=s[i]
            if epsilon>0 and draws[i]<epsilon:y=1-y;flips+=1
            osum[i]+=float(y)-float(oh[pos,i]);oh[pos,i]=y;th[pos,i]=s[i]
            if kernel==0:r[i]=osum[i]/H
            elif kernel==1:r[i]=(1-alpha)*r[i]+alpha*y
            elif H==1:r[i]=y
            else:r[i]=a*y+(1-a)*(osum[i]-y)/(H-1)
        if t>=burn and (t-burn)%every==0:
            fc[record]=nc/n
        if t>=burn and (t-burn)%every==0 and 0<nc<n:
            rc=0.;rd=0.
            for i in range(n):
                if s[i]:rc+=r[i]
                else:rd+=r[i]
            # Accumulate conditional sums per round to avoid cancellation near absorption.
            sr+=rd;sxr+=rc
            if 0<nc<n:gtime[record]=rc/nc-rd/(n-nc)
            sx+=nc;sample_count+=n
            for lag in range(H):
                old=(pos-lag)%H
                ls=0.;xs=0.;ys=0.;xys=0.
                for i in range(n):
                    ls+=th[old,i];xs+=s[i]*th[old,i];ys+=oh[old,i];xys+=s[i]*oh[old,i]
                lag_sum[lag]+=ls;cross[lag]+=xs;obs_sum[lag]+=ys;obs_cross[lag]+=xys
        if t>=burn and (t-burn)%every==0:record+=1
        if save_trace:trace[t]=s
        if sequential:
            irng.shuffle(order)
            for i in order:
                j=_reference(r,nn,i,irng.random(),w,eta)
                if irng.random()<_accept(_gain(s,nn,i,b),_gain(s,nn,j,b),K):s[i]=s[j]
        else:
            for i in range(n):F[i]=_gain(s,nn,i,b)
            choices=irng.random(n);acceptance=irng.random(n)
            for i in range(n):
                j=_reference(r,nn,i,choices[i],w,eta)
                out[i]=s[j] if acceptance[i]<_accept(F[i],F[j],K) else s[i]
            s[:]=out
    return s,fc,gtime,lag_sum,cross,obs_sum,obs_cross,sx,sr,sxr,sample_count,flips,trace,False

def simulate(p,rep=0,burn=5000,measure=2000,every=10,seed=None,early=True,save_trace=False):
    if burn<0 or measure<=0 or every<1 or measure%every or measure//every%2:raise ValueError('Need even measurement count')
    seed=seed_for(p,rep) if seed is None else int(seed)
    initial_seed,imitation_seed,error_seed=np.random.SeedSequence(seed).spawn(3)
    initial=np.zeros(p.L*p.L,np.uint8);initial[:len(initial)//2]=1;np.random.default_rng(initial_seed).shuffle(initial)
    h,tail=weights(p);start=time.perf_counter()
    result=_run(initial,neighbors(p.L),h,['W','E','PP'].index(p.kernel),2/(p.M+1),p.a,p.b,p.epsilon,p.w,p.K,p.eta,
        p.update=='sequential',burn,measure,every,np.random.default_rng(imitation_seed),np.random.default_rng(error_seed),early,save_trace)
    s,fc,gt,ls,xs,ys,xys,sx,sr,sxr,count,flips,trace,absorbed_early=result
    C=np.full(len(h),np.nan);Co=C.copy();G_pool=float('nan');G_observed=float('nan');G_pred=float('nan')
    if count>0 and 0<sx/count<1:
        prevalence=sx/count;var=prevalence*(1-prevalence)
        C=(xs/count-prevalence*ls/count)/var;Co=(xys/count-prevalence*ys/count)/var
        G_pool=sxr/sx-sr/(count-sx)
        G_pred=float((1-2*p.epsilon)*(h@C));G_observed=float(h@Co)
    valid=np.isfinite(gt);half=len(fc)//2
    row=p.canonical()|dict(realization=int(rep),seed=str(seed),T0=burn,T1=measure,sample_every=every,
        f_C=float(fc.mean()),all_D=int(s.sum()==0),all_C=int(s.sum()==len(s)),G=float(gt[valid].mean()) if valid.any() else float('nan'),
        G_approx=G_pred,G_from_Cy=G_observed,G_time=float(gt[valid].mean()) if valid.any() else float('nan'),
        n_pairs=int(count),coexistence_samples=int(valid.sum()),half_difference=float(fc[half:].mean()-fc[:half].mean()),
        elapsed_seconds=time.perf_counter()-start)
    row['G']=G_pool
    row['relative_gap_discrepancy']=(G_pool-G_pred)/G_pool if np.isfinite(G_pool) and G_pool!=0 else float('nan')
    row['kernel_tail_weight']=tail
    row['identity_absolute_error']=abs(G_pool-G_observed)
    row['identity_tolerance']=1e-9*abs(G_pool)+tail
    if np.isfinite(G_pool) and row['identity_absolute_error']>row['identity_tolerance']:
        raise ArithmeticError(f"Gap identity failed: {row['identity_absolute_error']} > {row['identity_tolerance']}")
    for lag,value in enumerate(C):row[f'C_s_{lag}']=float(value)
    for lag,value in enumerate(Co):row[f'C_y_{lag}']=float(value)
    diagnostics=dict(G_pooled=G_pool,G_time=row['G_time'],G_from_observed_covariance=G_observed,ema_tail=tail,
        noise_flips=flips,observation_count=(burn+measure)*p.L*p.L,early_absorbed=bool(absorbed_early))
    return dict(row=row,diagnostics=diagnostics,final=s,trace=trace)
