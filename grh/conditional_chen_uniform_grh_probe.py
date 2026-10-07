"""FLOAT ONLY: select GRH parameters for the ordered product certificate.

Newton subtractions and clipped tails are exploratory. The selected parameters
must be freshly checked by positive Arb recurrences and a fresh uniform-product
certificate before being adopted in the manuscript.
"""
import json
import math
from functools import lru_cache
from fractions import Fraction
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from flint import arb
from conditional_chen_rankin_certificate import odd_primes_to, degree_thresholds
from conditional_chen_adaptive_probe import exclusive

G = math.exp(.5772156649015329)
U = 1.32/G


@lru_cache(maxsize=None)
def sieve_scalars(alpha):
    t=.5-alpha
    s=8*t
    k=8*(t-1/3)
    F=2*G/s*(1+quad(lambda x:math.log(x-2)/(x-1),3,s)[0])
    split=t-.25
    Jh=quad(lambda b:math.exp(2-8*(t-b))/b,.125,split)[0]
    Jh+=math.log((1/3)/split)
    return t,s,k,F,Jh


def surplus(L,cost,eps,ep,em):
    alpha=math.ceil((cost+3.5*math.log(L))/L*1e6)/1e6
    if not 0 < alpha < 1/16:
        return -1000,alpha
    t,s,k,F,Jh=sieve_scalars(alpha)
    he=3*math.exp(2-s)/s
    low=8*((1+ep)*(2*G*math.log(s-1)/s-108*eps*he)
           -(ep+em)*(F+106*eps*he))
    I=G/(4*t)*math.log(6*(3-8*alpha)/(3-18*alpha))
    FM=2*G/k+106*eps
    upper=4*(1+ep)*(1+2/L)*(I+106*eps*Jh+2*FM/L**2)
    J=math.log(8/3)
    switch=(1+ep)/2*(2*G/t+318*eps)*(1+1e-8)*(1+4/L)*(
        .3630837292483046+2/L**2+(J+2/L**2)*
        (6*math.log1p(1e-8)/L+8/L**2))
    return low-upper-switch-.18/(U*math.sqrt(L)),alpha


def explore():
    envelope=json.loads(Path(__file__).with_name(
        "conditional_chen_uniform_product_certificate.json").read_text())
    rows={r["T"]:r for r in envelope["rows"]}
    Ts=[100000,200000,500000,1000000,2000000,3000000,
        5000000,10000000,20000000,30000000]
    starts=[24000,26000,28000,30000,32000,34000]
    best={L:[] for L in starts}
    all_primes=np.array(odd_primes_to(max(Ts)),dtype=float)
    for T in Ts:
        eps=float(arb(rows[T]["ordered_uniform_epsilon_ball"]))*1.000001
        if eps>1e-4:
            continue
        ps=all_primes[all_primes<=T]
        lp=np.log(ps)
        g=1/(ps-1)
        invV=np.exp(exclusive(-np.log1p(-g)))
        for sn,sd in [(1,6),(2,11),(1,5),(2,9),(1,4)]:
            sigma=sn/sd
            h=g*np.exp(sigma*lp)
            powers=[None]+[exclusive(h**j) for j in range(1,8)]
            elementary=[np.ones(len(ps))]
            for k in range(1,8):
                value=np.zeros(len(ps))
                for j in range(1,k+1):
                    value+=(-1)**(j-1)*elementary[k-j]*powers[j]
                elementary.append(value/k)
            plus=np.exp(exclusive(np.log1p(h)))
            minus=np.exp(exclusive(np.log1p(-h)))
            even,odd=(plus+minus)/2,(plus-minus)/2
            base=g/(1-g)*np.exp(3*sigma*lp)*invV
            suffix={}
            for m in range(9):
                if m>=3:
                    suffix[m]=(
                        np.append(np.cumsum((base*np.maximum(even,0))[::-1])[::-1],0),
                        np.append(np.cumsum((base*np.maximum(odd,0))[::-1])[::-1],0))
                if m<8:
                    if m%2:
                        odd-=elementary[m]
                    else:
                        even-=elementary[m]
            for a in range(44,61):
                u=Fraction(a,8)
                thresholds=degree_thresholds(T,u)
                global_degree=math.floor(u-2)
                for w in [3,5,7,11]:
                    idxw=np.searchsorted(ps,w,side="right")
                    ep=em=0.
                    for m in range(global_degree,math.ceil(u)+1):
                        lo=max(idxw,np.searchsorted(ps,thresholds[m]))
                        hi=(len(ps) if m==global_degree
                            else np.searchsorted(ps,thresholds[m-1]))
                        if lo<hi:
                            vp,vm=suffix[m]
                            ep+=max(0,vp[lo]-vp[hi])
                            em+=max(0,vm[lo]-vm[hi])
                    divisor=math.exp(float(u)*sigma*math.log(T))
                    ep/=divisor
                    em/=divisor
                    cost=sum(math.log(p) for p in ps[ps<=w])+float(u)*math.log(T)
                    for L in starts:
                        margin,alpha=surplus(L,cost,eps,ep,em)
                        best[L].append({
                            "L":L,"T":T,"w":w,"u":str(u),
                            "sigma":f"{sn}/{sd}","epsilon":eps,
                            "upper_gap":ep,"lower_gap":em,
                            "log_support_cost":cost,"alpha":alpha,
                            "margin":margin,
                        })
        print(f"Uniform GRH FLOAT probe completed T={T}",flush=True)
    return {
        "status":"FLOAT ONLY: Newton subtractions and clipped tails are not certificates",
        "results":{str(L):sorted(rows,key=lambda r:r["margin"],reverse=True)[:3]
                   for L,rows in best.items()},
    }


if __name__=="__main__":
    result=explore()
    Path(__file__).with_suffix(".json").write_text(
        json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))
