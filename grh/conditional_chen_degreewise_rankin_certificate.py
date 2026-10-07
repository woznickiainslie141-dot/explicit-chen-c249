"""Positive Arb bounds with a Rankin exponent chosen per prime and degree.

For each possible first-failing prefix, Markov's bound can use any sigma>=0.
Its complete-prime subset sums dominate every prime-deleted version. Taking
the minimum of finitely many such bounds at each minimum prime and degree
preserves domination. No floating arithmetic enters the bound.
"""
import json
import math
import time
from fractions import Fraction
from pathlib import Path
from flint import arb, ctx
from conditional_chen_rankin_certificate import odd_primes_to, degree_thresholds


def certify_degreewise(T, w, u, sigmas, M=12, progress=True):
    assert ctx.prec >= 192 and T >= 3 and 2 <= w < T
    u = Fraction(u)
    sigmas = tuple(Fraction(s) for s in sigmas)
    assert u > 3 and M >= math.ceil(u)
    assert sigmas and len(set(sigmas)) == len(sigmas)
    assert all(s >= 0 for s in sigmas)
    started = time.monotonic()
    primes = odd_primes_to(T)
    exact = [p for p in primes if p <= w]
    medium = [p for p in primes if p > w]
    Q = math.prod(exact)
    V_exact = arb(1)
    for p in exact:
        V_exact *= 1-arb(1)/(p-1)
    thresholds = degree_thresholds(T, u)
    minimum_degree = math.floor(u-2)
    states = []
    for sigma in sigmas:
        Ts = arb(T).root(sigma.denominator)**sigma.numerator
        normalizer = Ts.root(u.denominator)**u.numerator
        states.append({
            "sigma": sigma, "normalizer": normalizer,
            "small": [arb(1)]+[arb(0) for _ in range(M-1)],
            "even_tail": arb(0), "odd_tail": arb(0),
            "totals": [arb(0) for _ in range(M+2)],
        })
    invV = arb(1)
    totals = [arb(0) for _ in range(M+2)]
    for count, p in enumerate(reversed(medium), 1):
        g = arb(1)/(p-1)
        allowed = [k for k in range(minimum_degree, M)
                   if k >= len(thresholds) or p >= thresholds[k]]
        allowed.extend([M, M+1])
        best = {k: None for k in allowed}
        for state in states:
            sigma = state["sigma"]
            ps = arb(p).root(sigma.denominator)**sigma.numerator
            base = g/(1-g)*ps**3*invV/state["normalizer"]
            values = state["small"]+[state["even_tail"], state["odd_tail"]]
            for k in allowed:
                candidate = base*values[k]
                assert candidate >= 0
                state["totals"][k] += candidate
                best[k] = candidate if best[k] is None else best[k].min(candidate)
            # All coefficients and tails update using nonnegative additions.
            h = g*ps
            old_even, old_odd = state["even_tail"], state["odd_tail"]
            if M % 2:
                state["even_tail"] = old_even+h*old_odd
                state["odd_tail"] = old_odd+h*(old_even+state["small"][-1])
            else:
                state["even_tail"] = old_even+h*(old_odd+state["small"][-1])
                state["odd_tail"] = old_odd+h*old_even
            for k in range(M-1, 0, -1):
                state["small"][k] += h*state["small"][k-1]
        for k in allowed:
            totals[k] += best[k]
        invV /= 1-g
        if progress and count % 20000 == 0:
            print(f"Degreewise Rankin: {count}/{len(medium)} primes; "
                  f"{time.monotonic()-started:.1f}s", flush=True)

    def parity_sum(values):
        upper = sum((values[k] for k in range(M) if k % 2 == 0), arb(0))
        lower = sum((values[k] for k in range(M) if k % 2), arb(0))
        return upper+values[M], lower+values[M+1]

    upper, lower = parity_sum(totals)
    fixed = []
    for state in states:
        ep, em = parity_sum(state["totals"])
        assert upper <= ep and lower <= em
        fixed.append({"sigma": str(state["sigma"]),
                      "upper_gap_ball": str(ep), "lower_gap_ball": str(em)})
    cost = arb(Q).log()+arb(u.numerator)/u.denominator*arb(T).log()
    return {
        "status": "strict positive finite Arb bound; first-failure lemma remains analytic",
        "precision_bits": ctx.prec,
        "presieve_type": "pure Rosser medium weights, exact smallest-prime weights",
        "T": T, "exact_cutoff": w, "exact_presieve_Q_max": Q,
        "exact_product_ball": str(V_exact), "u": str(u), "D0": f"{T}^({u})",
        "sigma_candidates": [str(s) for s in sigmas],
        "separate_degree_cutoff": M,
        "complete_odd_prime_count": len(primes), "medium_prime_count": len(medium),
        "minimum_larger_subset_degree": minimum_degree,
        "degree_thresholds": thresholds,
        "upper_relative_Rosser_gap_ball": str(upper),
        "lower_relative_Rosser_gap_ball": str(lower),
        "relative_Rosser_gap_ball": str(upper+lower),
        "log_support_cost_ball": str(cost),
        "degree_contributions_ball": [str(v) for v in totals],
        "fixed_sigma_comparisons": fixed,
        "deleting_primes": "minimum of complete-prime dominating bounds remains dominating",
        "arithmetic": "positive elementary and tail recurrences; Arb min enclosures",
        "formula": "sum_r sum_degree min_sigma[D0^(-sigma) g(r)/(1-g(r)) "
                   "r^(3sigma) V_(>r)^(-1) e_degree(g(p)p^sigma)], "
                   "with retained parity degrees and two residual tails",
        "elapsed_seconds": round(time.monotonic()-started, 3),
    }


if __name__ == "__main__":
    ctx.prec = 192
    result = certify_degreewise(5000000, 5, Fraction(25, 4),
                               [Fraction(0), Fraction(1, 10), Fraction(7, 50),
                                Fraction(9, 50), Fraction(11, 50),
                                Fraction(7, 25), Fraction(9, 25)])
    Path(__file__).with_suffix(".json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
