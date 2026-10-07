"""Floating parameter exploration for finite-Rosser distribution reductions.

No output is an effective EH-only bound or a strict arithmetic certificate.
The finite distribution constants C,X remain unspecified inputs.
"""
import json
import math
from pathlib import Path
import numpy as np
from conditional_chen_rankin_certificate import odd_primes_to

G = math.exp(.5772156649015329)
U = 1.32/G


def lower_margin(theta, L, cost, epsilon, eta, C1, C2, error_inverse):
    s = 3*theta-3*cost/L
    if not 2 < s < 3:
        return -1000.0
    factor = 2*G*math.log(s-1)/s-C2*epsilon*math.exp(2-s)
    factor -= eta*(2*G/s+C1*epsilon)
    rho = 1-20.508e-6/L-2*L*math.exp(-L/3)
    deleted = L**3/math.log(2)*(2*math.exp(-L/2)+1.1*theta*L*math.exp(-L))
    return 3*U*rho*factor-1/error_inverse-deleted-L*L*math.exp(-L)


def explore():
    Ts = [10_000,22_000,50_000,100_000,300_000,1_000_000]
    ws = [2,3,5,7,11,13,17,19,23,29,31,37,43]
    all_primes = np.array(odd_primes_to(max(Ts)),dtype=np.int64)
    logQ = {w:sum(math.log(p) for p in all_primes[all_primes<=w]) for w in ws}
    cases = {(.75,20,L):[] for L in [900,1000,1100,1200,1250,1300,1400]}
    cases.update({(.9,10,L):[] for L in [225,250,275,300,325,350]})
    for T in Ts:
        ps = all_primes[all_primes<=T].astype(float)
        lp = np.log(ps)
        g = 1/(ps-1)
        lnT = math.log(T)
        b = max(2/(math.sqrt(T)*lnT),6.836e-6/math.log(1e12))
        epsilon = math.expm1(b+2.964e-6/lnT+1/(T-1)+.54/T)*1.000001
        if epsilon <= 1e-4:
            C1,C2 = 106,108
        elif epsilon <= 1/2000:
            C1,C2 = 109,110
        elif epsilon <= 1/1000:
            C1,C2 = 113,114
        elif epsilon <= 1/700:
            C1,C2 = 116,117
        elif epsilon <= 1/600:
            C1,C2 = 118,119
        else:
            continue
        for lam in np.arange(1.,5.01,.1):
            sigma = lam/lnT
            factors = np.log1p(g*np.exp(sigma*lp))-np.log1p(-g)
            larger = np.cumsum(factors[::-1])[::-1]-factors
            terms = np.log(g/(1-g))+3*sigma*lp+larger
            shift = terms.max()
            suffix = np.cumsum(np.exp(terms-shift)[::-1])[::-1]
            for w in ws:
                idx = np.searchsorted(ps,w,side='right')
                base = math.log(suffix[idx])+shift
                for u in np.arange(3.,7.01,.125):
                    eta = math.exp(base-u*lam)
                    if eta > .5:
                        continue
                    cost = logQ[w]+u*lnT
                    for (theta,budget,L),rows in cases.items():
                        margin = lower_margin(theta,L,cost,epsilon,eta,C1,C2,budget)
                        rows.append({'theta':theta,'L':L,'T':T,'w':w,'u':float(u),
                                     'sigma':sigma,'lambda':float(lam),
                                     'epsilon_bound':epsilon,'eta_bound':eta,
                                     'log_support_cost':cost,'C1':C1,'C2':C2,
                                     'error_inverse':budget,'margin':margin})
        print(f'Distribution search completed T={T}',flush=True)
    return {'status':'FLOATING EXPLORATION ONLY; C,X remain unknown distribution inputs',
            'results':{f'{theta}_{L}':sorted(rows,key=lambda r:r['margin'],reverse=True)[:3]
                       for (theta,budget,L),rows in cases.items()}}


if __name__ == '__main__':
    result = explore()
    Path(__file__).with_suffix('.json').write_text(
        json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
