"""Strict finite first-failure bounds retaining the preceding prefix pass.

When D0>T^3 every first-failing prefix satisfies d<D0. Thus
delta=log(d*r^2/D0) lies in [0,2 log r). Nonnegative tilted polynomial
majorants on this interval can dominate complete and prime-deleted sums.
Raw logarithmic moments use positive exponential-generating recurrences.
Shifted Gram calculations use Arb enclosures and strict optimizer signs.
"""
import json
import math
import time
from fractions import Fraction
from pathlib import Path
from flint import arb, arb_poly, ctx
from conditional_chen_rankin_certificate import odd_primes_to, degree_thresholds


def positive_shift(values, coefficients):
    """EGF convolution: add one atom to the logarithm of a subset product."""
    return [sum((coefficients[j-a]*values[a] for a in range(j+1)),arb(0))
            for j in range(len(values))]


def band_candidate(moments, R):
    """Best of zero, two boundary, and certified interior polynomial choices."""
    m0,m1,m2,m3,m4 = moments
    v0,v1 = m1,R*m1-m2
    A,B,C = m2,R*m2-m3,R*R*m2-2*R*m3+m4
    best = m0
    accepted = 0
    for v,denominator in [(v0,A),(v1,C)]:
        if v < 0 and denominator > 0:
            value = (m0-v*v/denominator).max(arb(0))
            best = best.min(value)
            accepted += 1
    determinant = A*C-B*B
    if determinant > 0:
        e = (B*v1-C*v0)/determinant
        b = (B*v0-A*v1)/determinant
        if e >= 0 and b >= 0:
            value = (m0+v0*e+v1*b).max(arb(0))
            best = best.min(value)
            accepted += 1
    return best,accepted


def certify_band(T,w,u,ordinary_sigmas,band_sigmas,M=12,progress=True,
                 atoms=None,log_base=None):
    assert ctx.prec >= 192 and T >= 3 and 2 <= w < T
    u = Fraction(u)
    ordinary_sigmas = tuple(Fraction(s) for s in ordinary_sigmas)
    band_sigmas = tuple(Fraction(s) for s in band_sigmas)
    assert u > 3 and M >= math.ceil(u)
    assert ordinary_sigmas and all(s >= 0 for s in ordinary_sigmas)
    assert band_sigmas and all(s >= 0 for s in band_sigmas)
    assert atoms is not None or log_base is None
    if log_base is not None:
        assert log_base > 1
    primes = odd_primes_to(T) if atoms is None else sorted(atoms)
    assert len(primes) == len(set(primes)) and all(3 <= p <= T for p in primes)
    exact = [p for p in primes if p <= w]
    medium = [p for p in primes if p > w]
    Q = math.prod(exact)
    V_exact = math.prod((1-arb(1)/(p-1) for p in exact),start=arb(1))
    thresholds = degree_thresholds(T,u)
    minimum_degree = math.floor(u-2)
    log_scale = arb(1) if log_base is None else arb(log_base).log()
    logT = arb(T).log()/log_scale
    states = {}
    for sigma in sorted(set(ordinary_sigmas) | set(band_sigmas)):
        length = 5 if sigma in band_sigmas else 1
        Ts = arb(T).root(sigma.denominator)**sigma.numerator
        states[sigma] = {
            "normalizer":Ts.root(u.denominator)**u.numerator,
            "length":length,
            "small":[arb_poly([1] if k == 0 else []) for k in range(M)],
            "even_tail":arb_poly([]),"odd_tail":arb_poly([])}
    totals = [arb(0) for _ in range(M+2)]
    ordinary_totals = [arb(0) for _ in range(M+2)]
    invV = arb(1)
    accepted_total = 0
    started = time.monotonic()
    factorials = [math.factorial(j) for j in range(5)]
    for count,p in enumerate(reversed(medium),1):
        g = arb(1)/(p-1)
        base = g/(1-g)*invV
        lp = arb(p).log()/log_scale
        a = arb(u.numerator)/u.denominator*logT-3*lp
        R = 2*lp
        delta_shift = arb_poly([(-a)**j/factorials[j] for j in range(5)])
        lp_egf = [lp**j/factorials[j] for j in range(5)]
        allowed = [k for k in range(minimum_degree,M)
                   if k >= len(thresholds) or p >= thresholds[k]]
        allowed += [M,M+1]
        best = {k:None for k in allowed}
        ordinary_best = {k:None for k in allowed}
        coefficients = {}
        for sigma,state in states.items():
            ps = arb(p).root(sigma.denominator)**sigma.numerator
            scale = ps**3/state["normalizer"]
            values = state["small"]+[state["even_tail"],state["odd_tail"]]
            coefficients[sigma] = arb_poly([g*ps*v for v in lp_egf[:state["length"]]])
            for k in allowed:
                candidate = scale*values[k][0]
                if sigma in ordinary_sigmas:
                    ordinary_best[k] = (candidate if ordinary_best[k] is None
                                        else ordinary_best[k].min(candidate))
                if sigma in band_sigmas:
                    shifted_poly = (values[k]*delta_shift).truncate(5)*scale
                    shifted = [shifted_poly[j]*factorials[j] for j in range(5)]
                    candidate,accepted = band_candidate(shifted,R)
                    accepted_total += accepted
                best[k] = candidate if best[k] is None else best[k].min(candidate)
        for k in allowed:
            assert not best[k] < 0
            totals[k] += base*best[k]
            ordinary_totals[k] += base*ordinary_best[k]
        for sigma,state in states.items():
            coefficients_sigma = coefficients[sigma]
            length = state["length"]
            old_even,old_odd = state["even_tail"],state["odd_tail"]
            if M % 2:
                even_add = (old_odd*coefficients_sigma).truncate(length)
                odd_add = ((old_even+state["small"][-1])*coefficients_sigma).truncate(length)
            else:
                even_add = ((old_odd+state["small"][-1])*coefficients_sigma).truncate(length)
                odd_add = (old_even*coefficients_sigma).truncate(length)
            state["even_tail"] = old_even+even_add
            state["odd_tail"] = old_odd+odd_add
            for k in range(M-1,0,-1):
                addition = (state["small"][k-1]*coefficients_sigma).truncate(length)
                state["small"][k] += addition
        invV /= 1-g
        if progress and count % 10000 == 0:
            print(f"Prefix band: {count}/{len(medium)}; "
                  f"{time.monotonic()-started:.1f}s",flush=True)
    def parity(values):
        return (sum((v for k,v in enumerate(values[:M]) if k % 2 == 0),arb(0))+values[M],
                sum((v for k,v in enumerate(values[:M]) if k % 2),arb(0))+values[M+1])
    upper,lower = parity(totals)
    ordinary_upper,ordinary_lower = parity(ordinary_totals)
    # The finite minimum proves the exact inequalities. Enclosures of equal
    # exact values can overlap, so reject only a certified reversed inequality.
    assert not upper > ordinary_upper and not lower > ordinary_lower
    return {
        "status":"strict Arb finite prefix-band bound; analytic majorant in manuscript",
        "precision_bits":ctx.prec,
        "prime_generation":"complete integer sieve" if atoms is None else "explicit test atoms",
        "logarithm_base":"e" if log_base is None else str(log_base),
        "presieve_type":"pure Rosser medium weights, exact smallest-prime weights",
        "T":T,"exact_cutoff":w,"exact_presieve_Q_max":Q,
        "exact_product_ball":str(V_exact),"u":str(u),"D0":f"{T}^({u})",
        "sigma_candidates":[str(s) for s in ordinary_sigmas],
        "band_sigma_candidates":[str(s) for s in band_sigmas],
        "separate_degree_cutoff":M,"complete_odd_prime_count":len(primes),
        "medium_prime_count":len(medium),"minimum_larger_subset_degree":minimum_degree,
        "degree_thresholds":thresholds,
        "upper_relative_Rosser_gap_ball":str(upper),
        "lower_relative_Rosser_gap_ball":str(lower),
        "relative_Rosser_gap_ball":str(upper+lower),
        "ordinary_upper_relative_gap_ball":str(ordinary_upper),
        "ordinary_lower_relative_gap_ball":str(ordinary_lower),
        "degree_contributions_ball":[str(v) for v in totals],
        "log_support_cost_ball":str(arb(Q).log()+arb(u.numerator)/u.denominator*arb(T).log()),
        "band_optimizers_with_strict_signs":accepted_total,
        "arithmetic":"positive raw EGF moments; Arb shifted Gram enclosures; proved nonnegative expectations",
        "deleting_primes":"dominate for each fixed nonnegative polynomial before optimizing",
        "formula":"sum_r,degree base * min(ordinary Rankin, exp(sigma*delta) "
                  "expectation of [1+e*delta+b*delta*(2log(r)-delta)]^2), e,b>=0",
        "elapsed_seconds":round(time.monotonic()-started,3)}


if __name__ == "__main__":
    ctx.prec = 192
    result = certify_band(10000,3,Fraction(35,8),
                          [Fraction(0)]+[Fraction(a,20) for a in range(2,13)],
                          [Fraction(3,25),Fraction(6,25),Fraction(9,25)],10)
    output = Path(__file__).with_name(Path(__file__).stem+"_standalone.json")
    output.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))
