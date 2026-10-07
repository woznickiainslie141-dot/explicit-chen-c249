"""Floating exploration of a direct finite Rankin bound for Rosser presieving.

Nothing printed by this script is a theorem or an interval certificate.
The proposed analytic bound must be audited before parameters are accepted.
"""
import json
import math
from pathlib import Path
import numpy as np

G=math.exp(.5772156649015329)
Umin=1.32/G
cb=.3630837292483046
J=math.log(8/3)
ws=[2,3,5,7,11,13,17,19,23,29,31,37,43]
Ts=[1_000_000,3_000_000,10_000_000,30_000_000]
mask=np.ones(max(Ts)+1,dtype=bool)
mask[:2]=False
for p in range(2,math.isqrt(max(Ts))+1):
 if mask[p]:mask[p*p::p]=False
primes=np.flatnonzero(mask)[1:]
del mask
logQs={w:sum(math.log(p) for p in primes[primes<=w]) for w in ws}


def scalar_margin(L,alpha,epsilon,eta,C1=106,C2=108):
 t=.5-alpha
 s=8*t
 k=8*(t-1/3)
 if not 3<s<4 or k<=1:return -1000
 low=8*(2*G*math.log(s-1)/s-C2*epsilon*math.exp(2)*3*math.exp(-s)/s-eta*(2*G/3+C1*epsilon))
 I=G/(4*t)*math.log(6*(3-8*alpha)/(3-18*alpha))
 FM=2*G/k+C1*epsilon
 upper=4*(1+eta)*(1+2/L)*(I+C1*epsilon*J+2*FM/L**2)
 switch=(1+eta)/2*(2*G/t+3*C1*epsilon)*(1+1e-8)*(1+4/L)*(cb+2/L**2+(J+2/L**2)*(6*math.log1p(1e-8)/L+8/L**2))
 return low-upper-switch-.18/(Umin*math.sqrt(L))


best={L:[] for L in [20000,25000,30000,35000,40000,45000]}
for T in Ts:
 ps=primes[primes<=T].astype(float)
 lp=np.log(ps)
 g=1/(ps-1)
 logT=math.log(T)
 b=max(2/(math.sqrt(T)*logT),6.836e-6/math.log(1e12))
 eps=math.expm1(b+2.964e-6/logT+1/(T-1)+.54/T)*1.000001
 if eps<=1e-4:C1,C2=106,108
 elif eps<=1/2000:C1,C2=109,110
 else:C1,C2=113,114
 for lam in np.linspace(1,5,41):
  sigma=lam/logT
  factors=np.log1p(g*np.exp(sigma*lp))-np.log1p(-g)
  reverse=np.cumsum(factors[::-1])[::-1]-factors
  term=np.log(g/(1-g))+3*sigma*lp+reverse
  shift=term.max()
  suffix=np.cumsum(np.exp(term-shift)[::-1])[::-1]
  for w in ws:
   idx=np.searchsorted(ps,w,side='right')
   base=math.log(suffix[idx])+shift
   logQ=logQs[w]
   for u in np.arange(4,8.01,.25):
    eta=math.exp(base-u*lam)
    cost=logQ+u*logT
    for L in best:
     alpha=(cost+3.5*math.log(L)+.0001)/L
     margin=scalar_margin(L,alpha,eps,eta,C1,C2)
     best[L].append({'L':L,'T':T,'w':w,'u':float(u),'lambda':float(lam),
                     'sigma':sigma,'epsilon_bound':eps,'eta_bound':eta,
                     'support_log_cost':cost,'alpha_min':alpha,
                     'C1':C1,'C2':C2,'margin':margin})
 result={str(L):sorted(v,key=lambda x:x['margin'],reverse=True)[:3] for L,v in best.items()}
 print('T completed',T,flush=True)

result={'status':'FLOATING EXPLORATION; Rankin/Rosser analytic interface is not yet certified',
        'results':result}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result,indent=2))
