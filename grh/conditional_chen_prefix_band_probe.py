"""FLOAT ONLY: tilted quartic bounds retaining d<D0 at first failure.

First failure with minimum r has D0/r^2 <= d < D0 when u>3.
For delta=log(d*r^2/D0), R=2 log r and e,b>=0, the function
exp(sigma*delta)*(1+e*delta+b*delta*(R-delta))^2 is >=1 on [0,R].
It is nonnegative everywhere. Floating Newton/log derivatives and clipped
tails below are exploratory and do not certify this polynomial expectation.
"""
import json
import math
import sys
from fractions import Fraction
from pathlib import Path
import numpy as np
from flint import arb
from conditional_chen_rankin_certificate import odd_primes_to, degree_thresholds
from conditional_chen_adaptive_probe import exclusive
from conditional_chen_uniform_grh_probe import surplus


def raw_log_moments(ps, sigma, M, order=4):
    """Newton derivatives; array [log-moment, degree bucket, minimum prime]."""
    lp = np.log(ps)
    h = np.exp(float(sigma)*lp)/(ps-1)
    powers = [[None]*(M) for _ in range(order+1)]
    for j in range(1,M):
        for a in range(order+1):
            powers[a][j] = exclusive(h**j*(j*lp)**a)
    e = np.zeros((order+1,M,len(ps)))
    e[0,0] = 1
    for k in range(1,M):
        for a in range(order+1):
            for j in range(1,k+1):
                for b in range(a+1):
                    e[a,k] += (-1)**(j-1)*math.comb(a,b)*e[b,k-j]*powers[a-b][j]/k
    def product_derivatives(sign):
        v = sign*h
        log0 = exclusive(np.log(np.abs(1+v)))
        p = np.exp(log0)
        if sign < 0:
            p *= np.where(exclusive((h > 1).astype(np.int64)) % 2, -1, 1)
        l1 = exclusive(lp*v/(1+v))
        l2 = exclusive(lp**2*v/(1+v)**2)
        l3 = exclusive(lp**3*v*(1-v)/(1+v)**3)
        l4 = exclusive(lp**4*v*(1-4*v+v*v)/(1+v)**4)
        return np.array([p,p*l1,p*(l2+l1*l1),
                         p*(l3+3*l1*l2+l1**3),
                         p*(l4+4*l1*l3+3*l2*l2+6*l1*l1*l2+l1**4)])
    plus, minus = product_derivatives(1), product_derivatives(-1)
    even = (plus+minus)/2-e[:,0:M:2].sum(axis=1)
    odd = (plus-minus)/2-e[:,1:M:2].sum(axis=1)
    return np.maximum(np.concatenate([e,even[:,None,:],odd[:,None,:]],axis=1),0)


def band_gaps(T, us, sigmas, ws, M=12):
    ps = np.array(odd_primes_to(T),dtype=np.int64)
    lp = np.log(ps)
    g = 1/(ps-1)
    base = g/(1-g)*np.exp(exclusive(-np.log1p(-g)))
    best = {u:np.full((M+2,len(ps)),np.inf) for u in us}
    ordinary = {u:np.full((M+2,len(ps)),np.inf) for u in us}
    R = 2*lp
    for sigma in sigmas:
        raw = raw_log_moments(ps,sigma,M)
        for u in us:
            a = float(u)*math.log(T)-3*lp
            moments = []
            scale = np.exp(-float(sigma)*a)
            for j in range(5):
                value = np.zeros((M+2,len(ps)))
                for k in range(j+1):
                    value += math.comb(j,k)*raw[k]*(-a)**(j-k)
                moments.append(value*scale)
            m0,m1,m2,m3,m4 = moments
            ordinary[u] = np.minimum(ordinary[u],m0)
            best[u] = np.minimum(best[u],m0)
            v0,v1 = m1,R*m1-m2
            A,B,C = m2,R*m2-m3,R*R*m2-2*R*m3+m4
            # The two boundary optima of the nonnegative quadrant.
            for v,denominator in [(v0,A),(v1,C)]:
                allowed = (v < 0) & (denominator > 0)
                subtract = np.zeros_like(m0)
                np.divide(v*v,denominator,out=subtract,where=allowed)
                candidate = np.where(allowed,np.maximum(m0-subtract,0),np.inf)
                best[u] = np.minimum(best[u],candidate)
            determinant = A*C-B*B
            e,b = np.zeros_like(m0),np.zeros_like(m0)
            good_det = determinant > 0
            np.divide(B*v1-C*v0,determinant,out=e,where=good_det)
            np.divide(B*v0-A*v1,determinant,out=b,where=good_det)
            allowed = good_det & (e >= 0) & (b >= 0)
            candidate = np.where(allowed,np.maximum(m0+v0*e+v1*b,0),np.inf)
            best[u] = np.minimum(best[u],candidate)
        print(f"Prefix band FLOAT probe T={T}, sigma={sigma}",flush=True)
    rows = []
    for u in us:
        thresholds = degree_thresholds(T,u)
        minimum_degree = math.floor(u-2)
        for w in ws:
            def parity(values):
                ep = em = 0.
                for k,value in enumerate(values):
                    if k < M:
                        if k < minimum_degree: continue
                        threshold = thresholds[k] if k < len(thresholds) else 1
                        include = (ps > w) & (ps >= threshold)
                        sign = k % 2
                    else:
                        include = ps > w
                        sign = k-M
                    contribution = float((base[include]*value[include]).sum())
                    if sign: em += contribution
                    else: ep += contribution
                return ep,em
            ep,em = parity(best[u])
            op,om = parity(ordinary[u])
            cost = sum(math.log(p) for p in ps[ps <= w])+float(u)*math.log(T)
            rows.append(dict(T=T,w=w,u=str(u),upper_gap=ep,lower_gap=em,
                             ordinary_upper=op,ordinary_lower=om,log_support_cost=cost))
    return rows


def explore(eh):
    uniform = json.loads(Path(__file__).with_name(
        "conditional_chen_uniform_product_certificate.json").read_text())
    uniform_rows = {r["T"]:r for r in uniform["rows"]}
    sigmas = [Fraction(a,50) for a in range(0,19,2)]
    if eh:
        from conditional_chen_local_eh_probe import endpoint,TABLE
        cases=[(.75,20,.01,833),(.9,10,.02,227),(.95,10,.02,179),(.99,10,.02,152)]
        best={str(theta):[] for theta,_,_,_ in cases}
        for T in [10000,22000]:
            rows=band_gaps(T,[Fraction(a,8) for a in range(32,43)],sigmas,[3,5],M=10)
            eps=float(arb(uniform_rows[T]["ordered_uniform_epsilon_ball"]))*1.000001
            C1,C2=next((a,b) for minimum,a,b in TABLE if eps < 1/minimum)
            for theta,budget,claim,high in cases:
                for row in rows:
                    if endpoint(theta,high,row,eps,C1,C2,budget)<=claim: continue
                    low,upper=83,high
                    while upper-low>1:
                        mid=(upper+low)//2
                        if endpoint(theta,mid,row,eps,C1,C2,budget)>claim: upper=mid
                        else: low=mid
                    best[str(theta)].append(dict(row,log_start=upper,
                        margin=endpoint(theta,upper,row,eps,C1,C2,budget)))
        results={k:sorted(v,key=lambda r:(r["log_start"],-r["margin"]))[:3]
                 for k,v in best.items()}
    else:
        starts=[30000,31000,31850]
        best={L:[] for L in starts}
        for T in [3000000,5000000]:
            rows=band_gaps(T,[Fraction(a,8) for a in range(47,52)],sigmas,[3,5,7],M=12)
            eps=float(arb(uniform_rows[T]["ordered_uniform_epsilon_ball"]))*1.000001
            for L in starts:
                for row in rows:
                    margin,alpha=surplus(L,row["log_support_cost"],eps,
                                         row["upper_gap"],row["lower_gap"])
                    best[L].append(dict(row,L=L,margin=margin,alpha=alpha))
        results={str(k):sorted(v,key=lambda r:-r["margin"])[:3] for k,v in best.items()}
    return {"status":"FLOAT ONLY: shifted moment cancellation and clipped tails",
            "sigma_candidates":[str(s) for s in sigmas],"best":results}


if __name__ == "__main__":
    eh=sys.argv[1:] == ["--eh"]
    result=explore(eh)
    output=Path(__file__).with_name("conditional_chen_prefix_band_eh_probe.json") if eh else Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))
