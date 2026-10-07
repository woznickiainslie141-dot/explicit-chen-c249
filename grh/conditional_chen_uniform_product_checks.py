"""Independent finite checks for the ordered-product certificate.

Integer exhaustive pairs check the interval-combination algebra. Direct rational
prime products check the left/right endpoint convention. These finite checks do
not prove the far-range input or the manuscript's general analytic lemmas.
"""
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path
from flint import arb, ctx
from conditional_chen_uniform_product_certificate import merge_ordered_statistics
from conditional_chen_rankin_certificate import odd_primes_to

ctx.prec = 192


def brute_statistics(values):
    return max(values), min(values), max(
        [0]+[values[i]-values[j] for i in range(len(values))
             for j in range(i+1, len(values))])


integer_checks = 0
for sequence in itertools.product(range(4), repeat=6):
    parts = [tuple(arb(x) for x in brute_statistics(sequence[i:i+2]))
             for i in [0, 2, 4]]
    combined = merge_ordered_statistics(
        merge_ordered_statistics(parts[0], parts[1]), parts[2])
    expected = brute_statistics(sequence)
    assert all(value == arb(target) for value, target in zip(combined, expected))
    integer_checks += 1

endpoint_pairs = 0
primes = odd_primes_to(100)
for T in [10, 20, 40]:
    events = [(T, False)]
    events += [(p, side) for p in primes if T < p < 100
               for side in [False, True]]
    events.append((100, False))
    values = []
    for x, right in events:
        product = math.prod((Fraction(p-2, p-1) for p in primes
                             if p < x or (right and p == x)), start=Fraction(1))
        value = arb(x).log().log()
        value += arb(product.numerator).log()-arb(product.denominator).log()
        values.append(value)
    for i, (u, u_right) in enumerate(events):
        for j in range(i+1, len(events)):
            v, v_right = events[j]
            included = [p for p in primes
                        if (p > u or (p == u and not u_right))
                        and (p < v or (p == v and v_right))]
            inverse = math.prod((Fraction(p-1, p-2) for p in included),
                                start=Fraction(1))
            direct = arb(u).log()/arb(v).log()
            direct *= arb(inverse.numerator)/inverse.denominator
            assert abs((values[i]-values[j]).exp()-direct) < arb("1e-45")
            endpoint_pairs += 1
    split = len(values)//2
    def interval_stats(items):
        peak = items[0]
        trough = items[0]
        for item in items[1:]:
            peak = peak.max(item)
            trough = trough.min(item)
        drawdown = arb(0)
        for i in range(len(items)):
            for j in range(i+1, len(items)):
                drawdown = drawdown.max(items[i]-items[j])
        return peak, trough, drawdown
    combined = merge_ordered_statistics(
        interval_stats(values[:split]), interval_stats(values[split:]))
    expected = interval_stats(values)
    assert all(abs(a-b) < arb("1e-45") for a, b in zip(combined, expected))

tail_checks = 0
for Q in [3, 10, 30]:
    for J in [Q+1, Q+17, Q+50]:
        whole = math.prod((1-Fraction(1, (n-1)**2)
                           for n in range(Q+1, J+1)), start=Fraction(1))
        assert whole == Fraction((Q-1)*J, Q*(J-1))
        assert whole > Fraction(Q-1, Q)
        prime_only = math.prod((1-Fraction(1, (p-1)**2)
                               for p in primes if Q < p <= J),
                              start=Fraction(1))
        assert prime_only >= whole
        tail_checks += 1

result = {
    "status": "PASS: finite ordered-product checks only",
    "integer_sequence_checks": integer_checks,
    "direct_prime_product_endpoint_pairs": endpoint_pairs,
    "integer_and_prime_tail_checks": tail_checks,
    "precision_bits": ctx.prec,
    "not_checked": "general analytic lemmas, published far-range envelopes, full Chen bridge",
}
Path(__file__).with_suffix(".json").write_text(
    json.dumps(result, indent=2)+"\n", encoding="utf-8")
print(json.dumps(result, indent=2))
