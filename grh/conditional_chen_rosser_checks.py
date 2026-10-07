"""Exact finite checks of the new Rosser/Rosser composite and gap identity.

All divisor patterns and all integer transition levels in the chosen finite
universes are checked. These checks do not prove the general sieve theorem.
"""
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path
from conditional_chen_weight_checks import divisor_sums, rosser, subset_products


def transitions(primes, minimum):
    points = set()
    for mask in range(1, 1 << len(primes)):
        chosen = sorted((p for i, p in enumerate(primes) if mask & (1 << i)),
                        reverse=True)
        prefix = 1
        for p in chosen:
            prefix *= p
            points.add(prefix*p*p)
    return sorted({minimum, *(x+j for x in points for j in (-1, 0, 1)
                              if x+j >= minimum)})


def mass(coefficients, primes):
    inverse_densities = subset_products([p-1 for p in primes])
    return sum((Fraction(c, d) for c, d in zip(coefficients, inverse_densities)),
               Fraction())


def first_failure_masses(primes, level, failure_test=None):
    totals = [Fraction(), Fraction()]
    for mask in range(1, 1 << len(primes)):
        chosen = sorted((p for i, p in enumerate(primes) if mask & (1 << i)),
                        reverse=True)
        n = len(chosen)
        prefix = 1
        keep = True
        for m, p in enumerate(chosen, 1):
            prefix *= p
            failure = (failure_test(prefix*p*p) if failure_test is not None
                       else prefix*p*p >= level)
            if m < n and m % 2 == n % 2 and failure:
                keep = False
                break
            if m == n and not failure:
                keep = False
        if keep:
            g = Fraction(1, math.prod(p-1 for p in chosen))
            Vbelow = math.prod(Fraction(p-2, p-1) for p in primes if p < chosen[-1])
            totals[n % 2] += g*Vbelow
    return totals  # even, odd


def rankin_sigma_one(primes, level):
    total = Fraction()
    product = Fraction(1)
    for p in sorted(primes, reverse=True):
        g = Fraction(1, p-1)
        total += g/(1-g)*p**3*product
        product *= (1+g*p)/(1-g)
    return total/level


def rankin_parities_sigma_one(primes, level, minimum_degree=0, degree_allowed=None):
    """Independent explicit subset sums, not the production tail recurrence."""
    values = [Fraction(),Fraction()]
    for r in primes:
        larger = [p for p in primes if p>r]
        invV = math.prod(Fraction(p-1,p-2) for p in larger)
        base = Fraction(r**3,r-2)*invV/level
        for mask in range(1 << len(larger)):
            degree = mask.bit_count()
            if degree < minimum_degree:
                continue
            if degree_allowed is not None and not degree_allowed(r,degree):
                continue
            h_product = math.prod(Fraction(p,p-1) for i,p in enumerate(larger)
                                  if mask & (1 << i))
            values[degree % 2] += base*h_product
    return values  # upper (odd prefix), lower (even prefix)


def check_gap_subsets(complete):
    checks = 0
    for level in transitions(complete, max(complete)**2):
        complete_bound = rankin_sigma_one(complete, level)
        complete_upper,complete_lower = rankin_parities_sigma_one(complete,level)
        for mask in range(1 << len(complete)):
            primes = [p for i, p in enumerate(complete) if mask & (1 << i)]
            bp, bm = rosser(primes, level, True), rosser(primes, level, False)
            V = math.prod(Fraction(p-2, p-1) for p in primes)
            Bp, Bm = mass(bp, primes), mass(bm, primes)
            even, odd = first_failure_masses(primes, level)
            assert Bp-V == odd and V-Bm == even
            assert Bp-Bm == even+odd
            subset_bound = rankin_sigma_one(primes, level)
            subset_upper,subset_lower = rankin_parities_sigma_one(primes,level)
            assert subset_upper+subset_lower == subset_bound
            assert Fraction(0) <= odd/V <= subset_upper <= complete_upper
            assert Fraction(0) <= even/V <= subset_lower <= complete_lower
            assert Fraction(0) <= (Bp-Bm)/V <= subset_bound <= complete_bound
            checks += 1
    return checks


def check_composite(medium, large, exact=()):
    Q = math.prod(exact)
    med_products = subset_products(medium)
    large_products = subset_products(large)
    exact_products = subset_products(exact)
    med_levels = transitions(medium, max(medium)**2)
    large_levels = transitions(large, (max(large)+1)**2)
    ns = len(medium)
    ne = len(exact)
    total_bits = ns+len(large)+ne
    patterns = 0
    level_pairs = 0
    for D0 in med_levels:
        bp, bm = rosser(medium, D0, True), rosser(medium, D0, False)
        ap, am = divisor_sums(bp), divisor_sums(bm)
        for mask in range(len(bp)):
            assert am[mask] <= int(mask == 0) <= ap[mask]
            if bp[mask] or bm[mask]:
                assert med_products[mask] < D0
        for D in large_levels:
            lp, lm = rosser(large, D, True), rosser(large, D, False)
            upper = [0] * (1 << total_bits)
            lower = upper[:]
            for mE in range(1 << ne):
                muE = (-1)**mE.bit_count()
                for m0 in range(1 << ns):
                    for m1 in range(1 << len(large)):
                        mask = m0 | (m1 << ns) | (mE << (ns+len(large)))
                        upper[mask] = muE*bp[m0]*lp[m1]
                        lower[mask] = muE*(bp[m0]*lm[m1]-(bp[m0]-bm[m0])*lp[m1])
                        assert abs(upper[mask]) <= 1 and abs(lower[mask]) <= 1
                        if upper[mask] or lower[mask]:
                            d = exact_products[mE]*med_products[m0]*large_products[m1]
                            assert d < Q*D0*D
            up, lo = divisor_sums(upper), divisor_sums(lower)
            for mask in range(len(up)):
                assert lo[mask] <= int(mask == 0) <= up[mask]
                patterns += 1
            level_pairs += 1
    return {'exact_primes':list(exact),'medium_primes':medium,'large_primes':large,
            'medium_transition_levels':len(med_levels),
            'large_transition_levels':len(large_levels),
            'level_pairs':level_pairs,'pointwise_pattern_checks':patterns,
            'gap_identity_and_deleted_subset_checks':check_gap_subsets(medium)}


def check_fractional_power(u_numerator, u_denominator):
    """Compare x with T^(a/b) by the exact integer test x^b >= T^a."""
    medium = [5,7,11,13,17,19,23,29]
    large = [31,37,41]
    exact = [3]
    T = max(medium)
    target = T**u_numerator

    def failed(x):
        return x**u_denominator >= target

    def coefficients(primes, upper):
        values = []
        for mask in range(1 << len(primes)):
            chosen = sorted((p for i,p in enumerate(primes) if mask & (1 << i)),
                            reverse=True)
            prefix = 1
            accepted = True
            for m,p in enumerate(chosen,1):
                prefix *= p
                if m % 2 == int(upper) and failed(prefix*p*p):
                    accepted = False
                    break
            values.append((-1)**len(chosen) if accepted else 0)
        return values

    bound_numerator = rankin_sigma_one(medium,1)
    degree = (u_numerator-2*u_denominator)//u_denominator
    allowed = lambda r,k: r**(3*u_denominator)*T**(k*u_denominator)>=target
    complete_upper,complete_lower = rankin_parities_sigma_one(medium,1,degree,allowed)
    gap_checks = 0
    for mask in range(1 << len(medium)):
        primes = [p for i,p in enumerate(medium) if mask & (1 << i)]
        bp,bm = coefficients(primes,True),coefficients(primes,False)
        V = math.prod(Fraction(p-2,p-1) for p in primes)
        Bp,Bm = mass(bp,primes),mass(bm,primes)
        even,odd = first_failure_masses(primes,None,failed)
        assert Bp-V == odd and V-Bm == even
        gap = (Bp-Bm)/V
        subset_numerator = rankin_sigma_one(primes,1)
        upper_numerator,lower_numerator = rankin_parities_sigma_one(primes,1,degree,allowed)
        assert upper_numerator <= complete_upper and lower_numerator <= complete_lower
        assert (odd/V)**u_denominator*target <= upper_numerator**u_denominator
        assert (even/V)**u_denominator*target <= lower_numerator**u_denominator
        assert subset_numerator <= bound_numerator
        assert gap >= 0 and gap**u_denominator*target <= subset_numerator**u_denominator
        gap_checks += 1

    bp,bm = coefficients(medium,True),coefficients(medium,False)
    products = subset_products(medium)
    for mask in range(len(bp)):
        if bp[mask] or bm[mask]:
            assert products[mask]**u_denominator < target
    Q = math.prod(exact)
    exact_products = subset_products(exact)
    large_products = subset_products(large)
    ns = len(medium)
    nl = len(large)
    patterns = 0
    for D in [(max(large)+1)**2,37**3,41**3+1]:
        lp,lm = rosser(large,D,True),rosser(large,D,False)
        upper = [0]*(1 << (ns+nl+len(exact)))
        lower = upper[:]
        for mE in range(len(exact_products)):
            for m0 in range(1 << ns):
                for m1 in range(1 << nl):
                    mask = m0 | (m1 << ns) | (mE << (ns+nl))
                    sign = (-1)**mE.bit_count()
                    upper[mask] = sign*bp[m0]*lp[m1]
                    lower[mask] = sign*(bp[m0]*lm[m1]-(bp[m0]-bm[m0])*lp[m1])
                    assert abs(upper[mask]) <= 1 and abs(lower[mask]) <= 1
                    if upper[mask] or lower[mask]:
                        d = exact_products[mE]*products[m0]*large_products[m1]
                        assert d**u_denominator < (Q*D)**u_denominator*target
        up,lo = divisor_sums(upper),divisor_sums(lower)
        for mask in range(len(up)):
            assert lo[mask] <= int(mask == 0) <= up[mask]
            patterns += 1
    return {'T':T,'u':f'{u_numerator}/{u_denominator}',
            'medium_primes':medium,'large_primes':large,'exact_primes':exact,
            'pointwise_pattern_checks':patterns,
            'gap_identity_and_deleted_subset_checks':gap_checks,
            'level_comparisons':'integer powers only; no floating roots',
            'nonzero_medium_gap':mass(bp,medium)-mass(bm,medium)>0}


if __name__ == '__main__':
    for Iplus, Iminus, Jplus, Jminus in itertools.product((0,1), repeat=4):
        coefficient = Iplus*Jminus-(Iplus-Iminus)*Jplus
        assert coefficient in (-1,0,1)
    universes = [check_composite([5,7,11],[13,17,19,23],(3,)),
                 check_composite([7,11,13],[17,19,23],(3,5)),
                 check_composite([3,5,7,11],[13,17,19])]
    fractional = [check_fractional_power(5,1),check_fractional_power(9,2),
                  check_fractional_power(35,8),check_fractional_power(25,4),
                  check_fractional_power(17,4),check_fractional_power(49,8),
                  check_fractional_power(39,8),check_fractional_power(33,8)]
    assert all(x['nonzero_medium_gap'] for x in fractional)
    result = {'status':'PASS: finite exact arithmetic checks only',
              'normalized_coefficient_boolean_checks':16,
              'universes':universes,
              'fractional_power_universes':fractional,
              'total_pointwise_pattern_checks':sum(u['pointwise_pattern_checks'] for u in universes+fractional),
              'total_gap_identity_and_deleted_subset_checks':sum(u['gap_identity_and_deleted_subset_checks'] for u in universes+fractional),
              'checks':['first-failure odd/even main-mass identities',
                        'Rankin upper bound with no factor two',
                        'separate upper/lower Rankin bounds by explicit subsets',
                        'rational-power levels exclude impossible low-degree prefixes',
                        'minimum-prime degree restrictions by exact integer powers',
                        'positive bound decreases when primes are deleted',
                        'exact/Rosser/Rosser composite pointwise brackets',
                        'coefficients in {-1,0,1} and support below Q D0 D'],
              'limitation':'Does not certify general analytic main-mass estimates or external inputs.'}
    Path(__file__).with_suffix('.json').write_text(
        json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,ensure_ascii=False,indent=2))
