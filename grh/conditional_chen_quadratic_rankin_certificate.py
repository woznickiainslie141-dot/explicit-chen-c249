"""Arb enclosures for nonnegative quadratic first-failure majorants.

Elementary coefficients and parity tails use positive recurrences. The
quadratic combination uses interval subtraction, with nonnegativity proved
termwise and by the finite weighted Cauchy--Schwarz inequality. It dominates
every prime-deleted subset sum; the optimized complete-prime bound does not
assert monotonicity of a quotient under deletion.
"""
import json
import math
import time
from fractions import Fraction
from pathlib import Path
from flint import arb, ctx
from conditional_chen_rankin_certificate import odd_primes_to, degree_thresholds


def quadratic_candidate(m0, m1, m2):
    """Optimized b in (X-b)^2/(1-b)^2, or None if signs are unresolved."""
    if not (m2 < m1 and m1 < m0):
        return None
    denominator = m0-2*m1+m2
    if not denominator > 0:
        return None
    b = (m1-m2)/(m0-m1)
    if not (b > 0 and b < 1):
        return None
    value = (m0*m2-m1**2)/denominator
    # The exact numerator is a finite weighted variance and is nonnegative.
    # Arb max encloses max(0, exact value), equal to that exact value.
    return value.max(arb(0))


def certify_quadratic(T, w, u, ordinary_sigmas, quadratic_sigmas,
                      M=12, progress=True):
    assert ctx.prec >= 192 and T >= 3 and 2 <= w < T
    u = Fraction(u)
    ordinary_sigmas = tuple(Fraction(s) for s in ordinary_sigmas)
    quadratic_sigmas = tuple(Fraction(s) for s in quadratic_sigmas)
    assert u > 3 and M >= math.ceil(u)
    assert ordinary_sigmas and all(s >= 0 for s in ordinary_sigmas)
    assert quadratic_sigmas and all(s > 0 for s in quadratic_sigmas)
    exponents = sorted(set(ordinary_sigmas) | {Fraction(0)} |
                       set(quadratic_sigmas) | {2*s for s in quadratic_sigmas})
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
    states = {}
    for sigma in exponents:
        Ts = arb(T).root(sigma.denominator)**sigma.numerator
        states[sigma] = {
            "normalizer": Ts.root(u.denominator)**u.numerator,
            "small": [arb(1)]+[arb(0) for _ in range(M-1)],
            "even_tail": arb(0), "odd_tail": arb(0)}
    invV = arb(1)
    totals = [arb(0) for _ in range(M+2)]
    ordinary_totals = [arb(0) for _ in range(M+2)]
    quadratic_accepted = 0
    for count, p in enumerate(reversed(medium), 1):
        g = arb(1)/(p-1)
        base = g/(1-g)*invV
        allowed = [k for k in range(minimum_degree, M)
                   if k >= len(thresholds) or p >= thresholds[k]]
        allowed.extend([M, M+1])
        moments, ps_values = {}, {}
        for sigma, state in states.items():
            ps = arb(p).root(sigma.denominator)**sigma.numerator
            ps_values[sigma] = ps
            factor = ps**3/state["normalizer"]
            values = state["small"]+[state["even_tail"], state["odd_tail"]]
            moments[sigma] = {k: factor*values[k] for k in allowed}
        for k in allowed:
            ordinary = moments[ordinary_sigmas[0]][k]
            for sigma in ordinary_sigmas[1:]:
                ordinary = ordinary.min(moments[sigma][k])
            best = ordinary
            for sigma in quadratic_sigmas:
                candidate = quadratic_candidate(
                    moments[Fraction(0)][k], moments[sigma][k],
                    moments[2*sigma][k])
                if candidate is not None:
                    best = best.min(candidate)
                    quadratic_accepted += 1
            # An enclosure of exact zero may have a tiny negative lower
            # endpoint after rounded max/min. The analytic value is >=0.
            assert not best < 0
            totals[k] += base*best
            ordinary_totals[k] += base*ordinary
        for sigma, state in states.items():
            h = g*ps_values[sigma]
            old_even, old_odd = state["even_tail"], state["odd_tail"]
            if M % 2:
                state["even_tail"] = old_even+h*old_odd
                state["odd_tail"] = old_odd+h*(old_even+state["small"][-1])
            else:
                state["even_tail"] = old_even+h*(old_odd+state["small"][-1])
                state["odd_tail"] = old_odd+h*old_even
            for k in range(M-1, 0, -1):
                state["small"][k] += h*state["small"][k-1]
        invV /= 1-g
        if progress and count % 20000 == 0:
            print(f"Quadratic Rankin: {count}/{len(medium)} primes; "
                  f"{time.monotonic()-started:.1f}s", flush=True)

    def parity_sum(values):
        upper = sum((values[k] for k in range(M) if k % 2 == 0), arb(0))
        lower = sum((values[k] for k in range(M) if k % 2), arb(0))
        return upper+values[M], lower+values[M+1]

    upper, lower = parity_sum(totals)
    ordinary_upper, ordinary_lower = parity_sum(ordinary_totals)
    assert upper <= ordinary_upper and lower <= ordinary_lower
    cost = arb(Q).log()+arb(u.numerator)/u.denominator*arb(T).log()
    return {
        "status": "strict Arb finite moment bound; analytic majorant remains proved in manuscript",
        "precision_bits": ctx.prec,
        "presieve_type": "pure Rosser medium weights, exact smallest-prime weights",
        "T": T, "exact_cutoff": w, "exact_presieve_Q_max": Q,
        "exact_product_ball": str(V_exact), "u": str(u), "D0": f"{T}^({u})",
        "sigma_candidates": [str(s) for s in ordinary_sigmas],
        "quadratic_sigma_candidates": [str(s) for s in quadratic_sigmas],
        "moment_exponents": [str(s) for s in exponents],
        "separate_degree_cutoff": M,
        "complete_odd_prime_count": len(primes), "medium_prime_count": len(medium),
        "minimum_larger_subset_degree": minimum_degree,
        "degree_thresholds": thresholds,
        "upper_relative_Rosser_gap_ball": str(upper),
        "lower_relative_Rosser_gap_ball": str(lower),
        "relative_Rosser_gap_ball": str(upper+lower),
        "ordinary_upper_relative_gap_ball": str(ordinary_upper),
        "ordinary_lower_relative_gap_ball": str(ordinary_lower),
        "degree_contributions_ball": [str(v) for v in totals],
        "log_support_cost_ball": str(cost),
        "quadratic_candidates_with_strict_signs": quadratic_accepted,
        "deleting_primes": "complete nonnegative polynomial expectation dominates; optimize only after domination",
        "arithmetic": "positive moment recurrences, Arb subtraction/division/min, proved nonnegative variance",
        "formula": "min(ordinary Rankin moments, (m0*m2-m1^2)/(m0-2*m1+m2)) "
                   "per r and degree when m2<m1<m0 and strict b in (0,1)",
        "elapsed_seconds": round(time.monotonic()-started, 3),
    }


if __name__ == "__main__":
    ctx.prec = 192
    result = certify_quadratic(
        5000000, 5, Fraction(25,4),
        [Fraction(0), Fraction(1,10), Fraction(7,50), Fraction(9,50),
         Fraction(11,50), Fraction(7,25), Fraction(9,25)],
        [Fraction(a,100) for a in range(4,21,2)])
    # The production JSON includes fresh uniform-product metadata.
    # A separate run must not overwrite that production artifact.
    Path(__file__).with_name(Path(__file__).stem+"_standalone.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
