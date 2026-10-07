"""Independent exact rational checks of the quadratic first-failure bound."""
import json
import math
from fractions import Fraction
from pathlib import Path
from flint import arb, ctx
from conditional_chen_quadratic_rankin_certificate import certify_quadratic
from conditional_chen_rosser_checks import first_failure_masses


def direct_bound(primes, T, u, ordinary_sigmas, quadratic_sigmas, M):
    assert all(s.denominator == 1 for s in ordinary_sigmas+quadratic_sigmas)
    exponents = set(ordinary_sigmas) | {Fraction(0)} | set(quadratic_sigmas)
    exponents |= {2*s for s in quadratic_sigmas}
    totals = [Fraction() for _ in range(M+2)]
    for r in primes:
        larger = [p for p in primes if p > r]
        base = Fraction(1,r-2)*math.prod(Fraction(p-1,p-2) for p in larger)
        moments = {s:[Fraction() for _ in range(M+2)] for s in exponents}
        for mask in range(1 << len(larger)):
            selected = [p for i,p in enumerate(larger) if mask & (1 << i)]
            k = len(selected)
            if k < math.floor(u-2) or r**3*T**k < T**u:
                continue
            x = Fraction(r**3*math.prod(selected), T**u)
            density = Fraction(1,math.prod(p-1 for p in selected))
            bucket = k if k < M else M+(k % 2)
            for sigma in exponents:
                moments[sigma][bucket] += density*x**sigma.numerator
        for k in range(M+2):
            candidates = [moments[s][k] for s in ordinary_sigmas]
            m0 = moments[Fraction(0)][k]
            for sigma in quadratic_sigmas:
                m1, m2 = moments[sigma][k], moments[2*sigma][k]
                if m2 < m1 < m0:
                    denominator = m0-2*m1+m2
                    b = (m1-m2)/(m0-m1)
                    assert denominator > 0 and 0 < b < 1
                    value = (m0*m2-m1*m1)/denominator
                    assert value >= 0
                    # Independently evaluate the optimized polynomial moment.
                    assert value == (m2-2*b*m1+b*b*m0)/(1-b)**2
                    candidates.append(value)
            totals[k] += base*min(candidates)
    upper = sum((v for k,v in enumerate(totals[:M]) if k % 2 == 0), Fraction())
    lower = sum((v for k,v in enumerate(totals[:M]) if k % 2), Fraction())
    return totals, upper+totals[M], lower+totals[M+1]


def run():
    ctx.prec = 192
    complete = [5,7,11,13,17,19,23,29]
    ordinary = [Fraction(0),Fraction(1),Fraction(2)]
    quadratic = [Fraction(1),Fraction(2)]
    def ball(v):
        return arb(v.numerator)/v.denominator
    count = 0
    configurations = []
    for u in [4,5,6]:
        for M in [u,u+1,12]:
            certified = certify_quadratic(29,3,Fraction(u),ordinary,quadratic,M,False)
            direct, upper, lower = direct_bound(complete,29,u,ordinary,quadratic,M)
            for a,b in zip(certified["degree_contributions_ball"],direct):
                assert abs(arb(a)-ball(b)) < arb("1e-45")
            assert abs(arb(certified["upper_relative_Rosser_gap_ball"])-ball(upper)) < arb("1e-45")
            assert abs(arb(certified["lower_relative_Rosser_gap_ball"])-ball(lower)) < arb("1e-45")
            for mask in range(1 << len(complete)):
                subset = [p for i,p in enumerate(complete) if mask & (1 << i)]
                V = math.prod(Fraction(p-2,p-1) for p in subset)
                even, odd = first_failure_masses(subset,29**u)
                _, subset_upper, subset_lower = direct_bound(
                    subset,29,u,ordinary,quadratic,M)
                assert odd/V <= subset_upper <= upper
                assert even/V <= subset_lower <= lower
                count += 1
            configurations.append({"u":u,"degree_cutoff":M,
                                   "strict_quadratic_candidates":
                                   certified["quadratic_candidates_with_strict_signs"]})
    assert any(r["strict_quadratic_candidates"] > 0 for r in configurations)
    return {"status":"PASS: finite exact rational checks only",
            "precision_bits":ctx.prec,"explicit_coefficient_configurations":len(configurations),
            "first_failure_and_deleted_prime_checks":count,
            "configurations":configurations,
            "not_checked":"external analytic sieve inputs or full Chen bridge"}


if __name__ == "__main__":
    result = run()
    Path(__file__).with_suffix(".json").write_text(
        json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))
