"""FLOAT ONLY: nonnegative quadratic majorants for first-failing products.

For X>=0 and 0<=b<1, 1_(X>=1) <= (X-b)^2/(1-b)^2.
Moments below are floating exploratory values, including cancellation in
Newton identities and parity tails. They are never a certificate.
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


def quadratic_gaps(T, us, sigmas, ws, M=12):
    ps = np.array(odd_primes_to(T), dtype=np.int64)
    lp = np.log(ps)
    g = 1/(ps-1)
    invV = np.exp(exclusive(-np.log1p(-g)))
    factor = g/(1-g)*invV
    cache = {}

    def moments(sigma):
        if sigma in cache:
            return cache[sigma]
        h = g*np.exp(float(sigma)*lp)
        powers = [None]+[exclusive(h**j) for j in range(1, M)]
        elementary = [np.ones(len(ps))]
        for k in range(1, M):
            value = np.zeros(len(ps))
            for j in range(1, k+1):
                value += (-1)**(j-1)*elementary[k-j]*powers[j]
            elementary.append(value/k)
        plus = np.exp(exclusive(np.log1p(h)))
        # Larger trial exponents can have h>1 at the smallest primes.
        # Preserve the sign of the suffix product; this remains FLOAT ONLY.
        minus = np.exp(exclusive(np.log(np.abs(1-h))))
        minus *= np.where(exclusive((h > 1).astype(np.int64)) % 2, -1, 1)
        even, odd = (plus+minus)/2, (plus-minus)/2
        for k, value in enumerate(elementary):
            if k % 2:
                odd -= value
            else:
                even -= value
        result = np.array([np.maximum(v, 0) for v in elementary]+[
            np.maximum(even, 0), np.maximum(odd, 0)])
        cache[sigma] = result
        return result

    zero = moments(Fraction(0))
    rows = []
    for u in us:
        best = zero.copy()
        ordinary = zero.copy()
        for sigma in sigmas:
            c = np.exp(float(sigma)*(3*lp-float(u)*math.log(T)))
            m1 = moments(sigma)*c
            m2 = moments(2*sigma)*c**2
            ordinary = np.minimum(ordinary, m1)
            best = np.minimum(best, m1)
            allowed = (m2 < m1) & (m1 < zero)
            denominator = zero-2*m1+m2
            allowed &= denominator > 0
            candidate = np.full_like(zero, np.inf)
            np.divide(zero*m2-m1*m1, denominator,
                      out=candidate, where=allowed)
            candidate = np.where(allowed, np.maximum(candidate, 0), np.inf)
            best = np.minimum(best, candidate)
        thresholds = degree_thresholds(T, u)
        minimum_degree = math.floor(u-2)
        for w in ws:
            def parity_sums(values):
                totals = [0., 0.]
                for k, value in enumerate(values):
                    if k < M:
                        if k < minimum_degree:
                            continue
                        threshold = thresholds[k] if k < len(thresholds) else 1
                        include = (ps > w) & (ps >= threshold)
                        sign = k % 2
                    else:
                        include = ps > w
                        sign = k-M
                    totals[sign] += float((factor[include]*value[include]).sum())
                return totals
            ep, em = parity_sums(best)
            op, om = parity_sums(ordinary)
            cost = sum(math.log(p) for p in ps[ps <= w])+float(u)*math.log(T)
            rows.append(dict(T=T, w=w, u=str(u), upper_gap=ep, lower_gap=em,
                             ordinary_upper=op, ordinary_lower=om,
                             log_support_cost=cost))
        print(f"Quadratic FLOAT probe: T={T}, u={u}", flush=True)
    return rows


def explore():
    uniform = json.loads(Path(__file__).with_name(
        "conditional_chen_uniform_product_certificate.json").read_text())
    uniform_rows = {r["T"]: r for r in uniform["rows"]}
    sigmas = [Fraction(a, 100) for a in range(4, 21, 2)]
    starts = [28000, 30000, 31000, 32000]
    best = {L: [] for L in starts}
    for T in [1000000, 3000000, 5000000]:
        eps = float(arb(uniform_rows[T]["ordered_uniform_epsilon_ball"]))*1.000001
        rows = quadratic_gaps(T, [Fraction(a, 8) for a in range(47, 52)],
                              sigmas, [3, 5, 7])
        for row in rows:
            for L in starts:
                margin, alpha = surplus(L, row["log_support_cost"], eps,
                                        row["upper_gap"], row["lower_gap"])
                best[L].append(dict(row, L=L, epsilon=eps, alpha=alpha,
                                    margin=margin))
    return {
        "status": "FLOAT ONLY: quadratic moment cancellation is not a certificate",
        "sigma_candidates": [str(s) for s in sigmas],
        "results": {str(L): sorted(rows, key=lambda r: r["margin"],
                                  reverse=True)[:3]
                    for L, rows in best.items()},
    }


def explore_eh():
    from conditional_chen_local_eh_probe import endpoint, TABLE
    uniform = json.loads(Path(__file__).with_name(
        "conditional_chen_uniform_product_certificate.json").read_text())
    uniform_rows = {r["T"]:r for r in uniform["rows"]}
    cases = [(.75,20,.01,840),(.9,10,.02,230),
             (.95,10,.02,181),(.99,10,.02,154)]
    best = {str(theta):[] for theta,_,_,_ in cases}
    sigmas = [Fraction(a,40) for a in range(2,15)]
    for T in [10000,22000]:
        rows = quadratic_gaps(T,[Fraction(a,8) for a in range(32,43)],
                              sigmas,[3,5],M=10)
        eps = float(arb(uniform_rows[T]["ordered_uniform_epsilon_ball"]))*1.000001
        C1,C2 = next((a,b) for minimum,a,b in TABLE if eps < 1/minimum)
        for theta,budget,claim,limit in cases:
            for row in rows:
                if endpoint(theta,limit,row,eps,C1,C2,budget) <= claim:
                    continue
                low, high = 83, limit
                while high-low > 1:
                    mid = (low+high)//2
                    if endpoint(theta,mid,row,eps,C1,C2,budget) > claim:
                        high = mid
                    else:
                        low = mid
                best[str(theta)].append(dict(row,theta=theta,log_start=high,
                    margin=endpoint(theta,high,row,eps,C1,C2,budget)))
    return {"status":"FLOAT ONLY: not an EH numeric record; C,X unspecified",
            "best":{k:sorted(v,key=lambda r:(r["log_start"],-r["margin"]))[:3]
                    for k,v in best.items()}}


if __name__ == "__main__":
    is_eh = sys.argv[1:] == ["--eh"]
    result = explore_eh() if is_eh else explore()
    output = (Path(__file__).with_name("conditional_chen_quadratic_eh_probe.json")
              if is_eh else Path(__file__).with_suffix(".json"))
    output.write_text(
        json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
