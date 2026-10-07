"""FLOAT ONLY: choose the Rankin exponent separately at each minimum prime.

This probe uses cancellation and clipped floating tails. Only an independent
positive Arb recurrence can certify the selected local-exponent bounds.
"""
import json
import math
from fractions import Fraction
from pathlib import Path
import numpy as np
from flint import arb
from conditional_chen_rankin_certificate import odd_primes_to, degree_thresholds
from conditional_chen_adaptive_probe import exclusive
from conditional_chen_uniform_grh_probe import surplus


def local_gaps(T, us, sigmas, ws):
    ps = np.array(odd_primes_to(T), dtype=np.int64)
    lp = np.log(ps)
    g = 1/(ps-1)
    invV = np.exp(exclusive(-np.log1p(-g)))
    max_degree = math.ceil(max(us))
    best = {u: [np.full(len(ps), np.inf), np.full(len(ps), np.inf)]
            for u in us}
    for sigma in sigmas:
        h = g*np.exp(float(sigma)*lp)
        powers = [None]+[exclusive(h**j) for j in range(1, max_degree+1)]
        elementary = [np.ones(len(ps))]
        for k in range(1, max_degree+1):
            value = np.zeros(len(ps))
            for j in range(1, k+1):
                value += (-1)**(j-1)*elementary[k-j]*powers[j]
            elementary.append(value/k)
        plus = np.exp(exclusive(np.log1p(h)))
        minus = np.exp(exclusive(np.log1p(-h)))
        even, odd = (plus+minus)/2, (plus-minus)/2
        tails = {}
        for k in range(max_degree+1):
            tails[k] = (np.maximum(even, 0), np.maximum(odd, 0))
            if k % 2:
                odd = odd-elementary[k]
            else:
                even = even-elementary[k]
        base = g/(1-g)*np.exp(3*float(sigma)*lp)*invV
        for u in us:
            global_degree = math.floor(u-2)
            thresholds = degree_thresholds(T, u)
            selected = [np.zeros(len(ps)), np.zeros(len(ps))]
            for k in range(global_degree, math.ceil(u)+1):
                lo = np.searchsorted(ps, thresholds[k])
                hi = (len(ps) if k == global_degree else
                      np.searchsorted(ps, thresholds[k-1]))
                if lo < hi:
                    for sign in range(2):
                        selected[sign][lo:hi] = tails[k][sign][lo:hi]
            normalizer = math.exp(float(u*sigma)*math.log(T))
            for sign in range(2):
                best[u][sign] = np.minimum(
                    best[u][sign], base*selected[sign]/normalizer)
    rows = []
    for u in us:
        for w in ws:
            include = ps > w
            ep, em = (float(v[include].sum()) for v in best[u])
            cost = sum(math.log(p) for p in ps[ps <= w])+float(u)*math.log(T)
            rows.append({"T": T, "w": w, "u": str(u),
                         "upper_gap": ep, "lower_gap": em,
                         "log_support_cost": cost})
    return rows


def degreewise_gaps(T, us, sigmas, ws, M=12):
    """Also vary sigma between exact degrees and the two residual tails."""
    ps = np.array(odd_primes_to(T), dtype=np.int64)
    lp = np.log(ps)
    g = 1/(ps-1)
    invV = np.exp(exclusive(-np.log1p(-g)))
    best = {u: [np.full(len(ps), np.inf) for _ in range(M+2)] for u in us}
    for sigma in sigmas:
        h = g*np.exp(float(sigma)*lp)
        powers = [None]+[exclusive(h**j) for j in range(1, M)]
        elementary = [np.ones(len(ps))]
        for k in range(1, M):
            value = np.zeros(len(ps))
            for j in range(1, k+1):
                value += (-1)**(j-1)*elementary[k-j]*powers[j]
            elementary.append(value/k)
        plus = np.exp(exclusive(np.log1p(h)))
        minus = np.exp(exclusive(np.log1p(-h)))
        even, odd = (plus+minus)/2, (plus-minus)/2
        for k, value in enumerate(elementary):
            if k % 2:
                odd = odd-value
            else:
                even = even-value
        moments = [np.maximum(v, 0) for v in elementary]
        moments += [np.maximum(even, 0), np.maximum(odd, 0)]
        base = g/(1-g)*np.exp(3*float(sigma)*lp)*invV
        for u in us:
            normalizer = math.exp(float(u*sigma)*math.log(T))
            for k, moment in enumerate(moments):
                best[u][k] = np.minimum(best[u][k], base*moment/normalizer)
    rows = []
    for u in us:
        thresholds = degree_thresholds(T, u)
        minimum_degree = math.floor(u-2)
        for w in ws:
            ep = em = 0.0
            for k, value in enumerate(best[u]):
                if k < M:
                    if k < minimum_degree:
                        continue
                    threshold = thresholds[k] if k < len(thresholds) else 1
                    include = (ps > w) & (ps >= threshold)
                    sign = k % 2
                else:
                    include = ps > w
                    sign = k-M
                addition = float(value[include].sum())
                if sign:
                    em += addition
                else:
                    ep += addition
            cost = sum(math.log(p) for p in ps[ps <= w])+float(u)*math.log(T)
            rows.append({"T": T, "w": w, "u": str(u),
                         "upper_gap": ep, "lower_gap": em,
                         "log_support_cost": cost})
    return rows


def explore():
    uniform = json.loads(Path(__file__).with_name(
        "conditional_chen_uniform_product_certificate.json").read_text())
    uniform_rows = {r["T"]: r for r in uniform["rows"]}
    sigmas = [Fraction(0)]+[Fraction(n, 100) for n in range(10, 41, 2)]
    starts = [24000, 26000, 28000, 30000, 31000, 32000, 33200]
    best = {L: [] for L in starts}
    for T in [500000, 1000000, 2000000, 3000000, 5000000, 10000000]:
        eps = float(arb(uniform_rows[T]["ordered_uniform_epsilon_ball"]))*1.000001
        rows = degreewise_gaps(T, [Fraction(a, 8) for a in range(46, 55)],
                              sigmas, [3, 5, 7, 11])
        for row in rows:
            for L in starts:
                margin, alpha = surplus(L, row["log_support_cost"], eps,
                                        row["upper_gap"], row["lower_gap"])
                best[L].append(dict(row, L=L, epsilon=eps, alpha=alpha,
                                    margin=margin))
        print(f"Local Rankin FLOAT probe completed T={T}", flush=True)
    return {
        "status": "FLOAT ONLY: cancellation and clipped tails are not certificates",
        "sigma_candidates": [str(s) for s in sigmas],
        "results": {str(L): sorted(rows, key=lambda r: r["margin"],
                                  reverse=True)[:3]
                    for L, rows in best.items()},
    }


if __name__ == "__main__":
    result = explore()
    Path(__file__).with_suffix(".json").write_text(
        json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
