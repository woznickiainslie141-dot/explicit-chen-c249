"""Finite, exact sanity checks for the manuscript's composite sieve weights.

These checks exercise the pointwise bounds and support across every divisor
pattern in three finite prime universes, including all Rosser transition levels.
They do not prove the general analytic main-mass estimates or the manuscript.
"""
import json
import math
from pathlib import Path


def subset_products(primes):
    products = [1] * (1 << len(primes))
    for mask in range(1, len(products)):
        bit = mask & -mask
        products[mask] = products[mask ^ bit] * primes[bit.bit_length() - 1]
    return products


def divisor_sums(coefficients):
    """Subset zeta transform: one output per possible prime divisor pattern."""
    result = coefficients[:]
    for bit in range(len(coefficients).bit_length() - 1):
        for mask in range(len(result)):
            if mask & (1 << bit):
                result[mask] += result[mask ^ (1 << bit)]
    return result


def rosser(primes, level, upper):
    result = []
    for mask in range(1 << len(primes)):
        chosen = sorted((p for i, p in enumerate(primes) if mask & (1 << i)),
                        reverse=True)
        prefix = 1
        accepted = True
        for m, p in enumerate(chosen, 1):
            prefix *= p
            if m % 2 == int(upper) and prefix * p * p >= level:
                accepted = False
                break
        result.append((-1) ** len(chosen) if accepted else 0)
    return result


def check_universe(small, large, depth, exact=()):
    assert depth % 2 == 0 and depth >= 2
    cutoff = max(small)
    z = max(large) + 1
    assert max(small) < min(large)
    assert not exact or max(exact) < min(small)
    Q = math.prod(exact)
    exact_products = subset_products(exact)
    small_products = subset_products(small)
    large_products = subset_products(large)
    transitions = set()
    for mask in range(1, 1 << len(large)):
        chosen = sorted((p for i, p in enumerate(large) if mask & (1 << i)),
                        reverse=True)
        prefix = 1
        for p in chosen:
            prefix *= p
            transitions.add(prefix * p * p)
    levels = sorted({z, z*z, *(d + j for d in transitions for j in (-1, 0, 1)
                             if d + j >= z)})
    bplus = [(-1)**m.bit_count() if m.bit_count() <= depth else 0
             for m in range(1 << len(small))]
    bminus = [(-1)**m.bit_count() if m.bit_count() < depth else 0
              for m in range(1 << len(small))]
    ap, am = divisor_sums(bplus), divisor_sums(bminus)
    for mask in range(len(ap)):
        indicator = int(mask == 0)
        assert am[mask] <= indicator <= ap[mask]
        assert ap[mask] >= 0
    ns = len(small)
    pattern_checks = 0
    for level in levels:
        lp, lm = rosser(large, level, True), rosser(large, level, False)
        cp, cm = divisor_sums(lp), divisor_sums(lm)
        for mask in range(len(cp)):
            indicator = int(mask == 0)
            assert cm[mask] <= indicator <= cp[mask]
            assert cp[mask] >= 0
            if lp[mask] or lm[mask]:
                assert large_products[mask] < level
        upper = [0] * (1 << (ns + len(large)))
        lower = upper[:]
        for m0 in range(1 << ns):
            for m1 in range(1 << len(large)):
                mask = m0 | (m1 << ns)
                upper[mask] = bplus[m0] * lp[m1]
                lower[mask] = (bplus[m0] * lm[m1]
                               - (bplus[m0] - bminus[m0]) * lp[m1])
                assert abs(upper[mask]) <= 1 and abs(lower[mask]) <= 1
                if upper[mask] or lower[mask]:
                    assert small_products[m0] * large_products[m1] < cutoff**depth * level
        upper_sum, lower_sum = divisor_sums(upper), divisor_sums(lower)
        # Exact inclusion--exclusion is a disjoint convolution. Check it as
        # integer coefficients and divisor sums, including every exact-prime
        # divisor pattern, instead of only checking the final scalar formula.
        full_upper = [0] * (len(upper) * len(exact_products))
        full_lower = full_upper[:]
        shift = ns + len(large)
        for mE in range(len(exact_products)):
            muE = (-1)**mE.bit_count()
            for mask in range(len(upper)):
                full_mask = mask | (mE << shift)
                full_upper[full_mask] = muE * upper[mask]
                full_lower[full_mask] = muE * lower[mask]
                assert abs(full_upper[full_mask]) <= 1 and abs(full_lower[full_mask]) <= 1
                if full_upper[full_mask] or full_lower[full_mask]:
                    m0 = mask & ((1 << ns)-1)
                    m1 = mask >> ns
                    divisor = exact_products[mE]*small_products[m0]*large_products[m1]
                    assert divisor < Q*cutoff**depth*level
        full_up_sum, full_lo_sum = divisor_sums(full_upper), divisor_sums(full_lower)
        for mask in range(len(full_upper)):
            indicator = int(mask == 0)
            assert full_lo_sum[mask] <= indicator <= full_up_sum[mask]
            pattern_checks += 1
    return {'exact_primes':list(exact),'Q':Q,
            'small_primes': small, 'large_primes': large, 'depth': depth,
            'integer_transition_levels': len(levels),
            'complete_prime_divisor_patterns_per_level': 1 << (ns + len(large) + len(exact)),
            'pointwise_pattern_checks': pattern_checks,
            'checks': ['Rosser pointwise brackets and support',
                       'Brun pointwise brackets', 'composite pointwise brackets',
                       'composite coefficients in {-1,0,1}',
                       'exact inclusion-exclusion pointwise brackets',
                       'composite support below Q T^depth D']}


if __name__ == '__main__':
    universes = [check_universe([3, 5, 7, 11], [13, 17, 19, 23, 29, 31], depth)
                 for depth in (2, 4)]
    universes += [check_universe([13, 17, 19, 23], [29, 31, 37, 41, 43, 47], depth)
                  for depth in (2, 4)]
    universes += [check_universe([7, 11, 13, 17], [19, 23, 29, 31, 37, 41], depth,
                                exact=(3,5)) for depth in (2,4)]
    result = {'status': 'PASS: finite exact arithmetic sanity checks only',
              'universes': universes,
              'total_pointwise_pattern_checks': sum(u['pointwise_pattern_checks']
                                                   for u in universes),
              'limitation': 'Does not certify general analytic main masses or the entire proof.'}
    Path(__file__).with_suffix('.json').write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))
