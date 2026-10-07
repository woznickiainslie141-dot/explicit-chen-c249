"""Compare positive degreewise recurrences with direct subset enumeration."""
import json
import math
from fractions import Fraction
from pathlib import Path
from flint import arb, ctx
from conditional_chen_degreewise_rankin_certificate import certify_degreewise
from conditional_chen_rosser_checks import first_failure_masses


def direct_bound(primes, T, u, sigmas, M):
    """Exact rational subset sums for integer exponents, without recurrences."""
    assert all(s.denominator == 1 for s in sigmas)
    totals = [Fraction() for _ in range(M+2)]
    for r in primes:
        larger = [p for p in primes if p > r]
        inverse = math.prod(Fraction(p-1,p-2) for p in larger)
        candidates = []
        for sigma in sigmas:
            values = [Fraction() for _ in range(M+2)]
            for mask in range(1 << len(larger)):
                selected = [p for i,p in enumerate(larger) if mask & (1 << i)]
                k = len(selected)
                if k < math.floor(u-2) or r**3*T**k < T**u:
                    continue
                product = math.prod(selected)
                density = Fraction(1, math.prod(p-1 for p in selected))
                exponent = (r**3*product)**sigma.numerator
                normalizer = (T**u)**sigma.numerator
                value = inverse*Fraction(1,r-2)*density*Fraction(exponent,normalizer)
                bucket = k if k < M else M+(k % 2)
                values[bucket] += value
            candidates.append(values)
        for k in range(M+2):
            local = candidates[0][k]
            for values in candidates[1:]:
                local = min(local,values[k])
            totals[k] += local
    upper = sum((v for k,v in enumerate(totals[:M]) if k % 2 == 0), Fraction())
    lower = sum((v for k,v in enumerate(totals[:M]) if k % 2), Fraction())
    return totals, upper+totals[M], lower+totals[M+1]


def run():
    ctx.prec = 192
    T, w = 29, 3
    complete = [5,7,11,13,17,19,23,29]
    sigmas = [Fraction(0),Fraction(1),Fraction(2)]
    def ball(value):
        return arb(value.numerator)/value.denominator
    count = 0
    configurations = []
    for u in [4,5,6]:
        for M in [u,u+1,12]:
            result = certify_degreewise(T,w,Fraction(u),sigmas,M,False)
            direct, upper, lower = direct_bound(complete,T,u,sigmas,M)
            for a,b in zip(result["degree_contributions_ball"],direct):
                assert abs(arb(a)-ball(b)) < arb("1e-45")
            assert abs(arb(result["upper_relative_Rosser_gap_ball"])-ball(upper)) < arb("1e-45")
            assert abs(arb(result["lower_relative_Rosser_gap_ball"])-ball(lower)) < arb("1e-45")
            for mask in range(1 << len(complete)):
                subset = [p for i,p in enumerate(complete) if mask & (1 << i)]
                V = math.prod(Fraction(p-2,p-1) for p in subset)
                even,odd = first_failure_masses(subset,T**u)
                actual_upper = odd/V
                actual_lower = even/V
                _,subset_upper,subset_lower = direct_bound(subset,T,u,sigmas,M)
                assert actual_upper <= subset_upper <= upper
                assert actual_lower <= subset_lower <= lower
                count += 1
            configurations.append({"u":u,"degree_cutoff":M,"prime_subsets":256})
    return {"status":"PASS: finite degreewise checks only","precision_bits":ctx.prec,
            "explicit_coefficient_configurations":len(configurations),
            "first_failure_and_deleted_prime_checks":count,
            "configurations":configurations,
            "not_checked":"general analytic sieve inputs or full Chen bridge"}


if __name__ == "__main__":
    result = run()
    Path(__file__).with_suffix(".json").write_text(
        json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))
