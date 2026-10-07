"""Floating exploration of parity-specific, minimum-prime-aware Rankin bounds.

Newton identities and subtractions are used only for parameter selection.
The production certificate uses positive Arb recurrences instead.
"""
import json
import math
from fractions import Fraction
from pathlib import Path
import numpy as np
from conditional_chen_rankin_certificate import odd_primes_to,degree_thresholds

G=math.exp(.5772156649015329)
U=1.32/G


def surplus(L,cost,eps,ep,em):
    alpha=math.ceil((cost+3.5*math.log(L))/L*1e6)/1e6
    t=.5-alpha;s=8*t;k=8*(t-1/3)
    lo=2*G*math.log(s-1)/s-108*eps*math.exp(2)*3*math.exp(-s)/s
    up=2*G/3+106*eps
    low=8*((1+ep)*lo-(ep+em)*up)
    I=G/(4*t)*math.log(6*(3-8*alpha)/(3-18*alpha))
    FM=2*G/k+106*eps
    upper=4*(1+ep)*(1+2/L)*(I+106*eps*math.log(8/3)+2*FM/L**2)
    switch=(1+ep)/2*(2*G/t+318*eps)*(1+1e-8)*(1+4/L)*(
        .3630837292483046+2/L**2+(math.log(8/3)+2/L**2)*
        (6*math.log1p(1e-8)/L+8/L**2))
    return low-upper-switch-.18/(U*math.sqrt(L)),alpha


def exclusive(v):
    return np.cumsum(v[::-1])[::-1]-v


def explore():
    T=30_000_000
    ps=np.array(odd_primes_to(T),dtype=float)
    lp=np.log(ps)
    g=1/(ps-1)
    logT=math.log(T)
    invV=np.exp(exclusive(-np.log1p(-g)))
    starts=[30000,32000,34000,35000,36000,37000]
    best={L:[] for L in starts}
    for sn,sd in [(1,6),(2,11),(1,5)]:
        sigma=sn/sd
        h=g*np.exp(sigma*lp)
        powers=[None]+[exclusive(h**j) for j in range(1,8)]
        elementary=[np.ones(len(ps))]
        for k in range(1,8):
            v=np.zeros(len(ps))
            for j in range(1,k+1):
                v+=(-1)**(j-1)*elementary[k-j]*powers[j]
            elementary.append(v/k)
        plus=np.exp(exclusive(np.log1p(h)))
        minus=np.exp(exclusive(np.log1p(-h)))
        even=(plus+minus)/2
        odd=(plus-minus)/2
        base=g/(1-g)*np.exp(3*sigma*lp)*invV
        suffix={}
        for m in range(9):
            if m>=3:
                suffix[m]=(np.append(np.cumsum((base*np.maximum(even,0))[::-1])[::-1],0),
                           np.append(np.cumsum((base*np.maximum(odd,0))[::-1])[::-1],0))
            if m<8:
                if m%2:odd=odd-elementary[m]
                else:even=even-elementary[m]
        for a in range(44,61):
            u=Fraction(a,8)
            thresholds=degree_thresholds(T,u)
            global_degree=math.floor(u-2)
            for w in [3,5,7,11]:
                idxw=np.searchsorted(ps,w,side='right')
                ep=em=0.
                for m in range(global_degree,math.ceil(u)+1):
                    lo=max(idxw,np.searchsorted(ps,thresholds[m]))
                    hi=(len(ps) if m==global_degree
                        else np.searchsorted(ps,thresholds[m-1]))
                    if lo<hi:
                        vp,vm=suffix[m]
                        ep+=max(0,vp[lo]-vp[hi])
                        em+=max(0,vm[lo]-vm[hi])
                normalizer=math.exp(float(u)*sigma*logT)
                ep/=normalizer;em/=normalizer
                cost=sum(math.log(p) for p in ps[ps<=w])+float(u)*logT
                for L in starts:
                    margin,alpha=surplus(L,cost,1/46600,ep,em)
                    best[L].append({'L':L,'T':T,'w':w,'u':str(u),
                                    'sigma':f'{sn}/{sd}','upper_gap':ep,'lower_gap':em,
                                    'log_support_cost':cost,'alpha':alpha,'margin':margin})
        print(f'Adaptive probe completed sigma={sn}/{sd}',flush=True)
    return {'status':'FLOATING ONLY; cancellation and clipped tails are not certificates',
            'results':{str(L):sorted(v,key=lambda r:r['margin'],reverse=True)[:3]
                       for L,v in best.items()}}


if __name__=='__main__':
    result=explore()
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
