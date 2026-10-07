"""Positive weighted counts for actual finite lower Rosser support.

The sum bounds coefficients by d**(-sigma), retaining all preceding
Rosser checks. No arithmetic-progression estimate is assumed or proved.
"""
import hashlib
import json
import math
import time
from fractions import Fraction
from pathlib import Path
from flint import arb, arb_poly, ctx
from conditional_chen_rankin_certificate import odd_primes_to
from conditional_chen_prefix_block_certificate import integer_floor, integer_ceil
from conditional_chen_sharp_support_certificate import support_ratio
from conditional_chen_accepted_support_certificate import append, STATES

ROOT = Path(__file__).resolve().parent


def q(value):
    value = Fraction(value)
    return arb(value.numerator) / value.denominator


def certify(T, w, u, sigma=Fraction(1, 2), h=Fraction(1, 100),
            available=None, progress=False, return_states=False):
    assert ctx.prec >= 192
    u, sigma, h = map(Fraction, [u, sigma, h])
    assert u > 3 and 0 <= sigma <= 1 and 0 < h < 1
    started = time.monotonic()
    primes = [p for p in odd_primes_to(T) if p > w]
    allowed = set(primes) if available is None else set(available)
    assert allowed <= set(primes)
    kappa, witness = support_ratio(primes)
    logD = q(u) * arb(T).log()
    logH = logD + q(kappa).log()
    J = integer_ceil(logH / q(h))
    product, K = 1, 0
    for p in primes:
        product *= p
        if product**u.denominator * kappa.denominator**u.denominator < (
                T**u.numerator * kappa.numerator**u.denominator):
            K += 1
        else:
            break
    groups = {}
    for p in reversed(primes):
        groups.setdefault(integer_floor(arb(p).log() / q(h)), []).append(p)
    states = [arb_poly([]) for _ in STATES]
    states[4] = arb_poly([1])
    for number, (label, ps) in enumerate(groups.items(), 1):
        cutoff = integer_ceil((logD - 3 * arb(min(ps)).log()) / q(h))
        assert label > 0 and cutoff > 0
        coeffs = arb_poly([1])
        size = 0
        for p in ps:
            if p in allowed:
                weight = (-q(sigma) * arb(p).log()).exp()
                coeffs = (coeffs * arb_poly([1, weight])).truncate(K + 1)
                size += 1
        powers = states
        updated = [arb_poly([]) for _ in STATES]
        for j in range(min(K, size, (J - 1) // label) + 1):
            for index in range(6):
                updated[index] += coeffs[j] * powers[index]
            powers = append(powers, label, cutoff, J)
            if not any(powers):
                break
        states = updated
        if progress and number % 100 == 0:
            print(f'Weighted support {number}/{len(groups)} blocks; '
                  f'{time.monotonic() - started:.1f}s', flush=True)
    totals = [p(arb(1)) for p in states]
    lower = sum((totals[i] for i in [2, 3, 4, 5]), arb(0))
    union = sum(totals, arb(0))
    # Arb intervals for equal sums can overlap without establishing <=.
    # Positivity and lower-state inclusion follow from the positive sums.
    assert lower > 0 and union > 0
    exact = odd_primes_to(w)
    presieve = math.prod((1 + (-q(sigma) * arb(p).log()).exp() for p in exact), start=arb(1))
    weighted = presieve * lower
    names = [Path(__file__).name, 'conditional_chen_accepted_support_certificate.py',
             'conditional_chen_rankin_certificate.py', 'conditional_chen_prefix_block_certificate.py',
             'conditional_chen_sharp_support_certificate.py']
    report = {
        'status': 'strict weighted lower-support upper bound; positive domination lemma required',
        'precision_bits': ctx.prec, 'T': T, 'exact_cutoff': w, 'u': str(u),
        'sigma': str(sigma), 'h': str(h), 'kappa': str(kappa), 'support_maximizer': witness,
        'maximum_actual_degree': K, 'label_limit': J, 'blocks': len(groups),
        'complete_medium_prime_count': len(primes), 'available_medium_prime_count': len(allowed),
        'lower_medium_weighted_sum_ball': str(lower), 'union_medium_weighted_sum_ball': str(union),
        'exact_presieve_weighted_factor_ball': str(presieve),
        'full_lower_weighted_sum_upper_ball': str(weighted),
        'state_weighted_balls': [{'alive_flags': flags, 'selected_parity': parity, 'weighted_sum_ball': str(value)}
                               for (flags, parity), value in zip(STATES, totals)],
        'complete_cutoffs_retained_after_deletion': True,
        'elapsed_seconds': round(time.monotonic() - started, 3),
        'input_sha256': {n: hashlib.sha256((ROOT / n).read_bytes()).hexdigest() for n in names},
        'not_checked': ['general weighted domination lemma', 'AP estimates', 'full splice'],
    }
    return (report, states) if return_states else report


def run():
    ctx.prec = 192
    result = certify(2000000, 3, Fraction(3001, 1000), progress=True)
    bridge_path = ROOT / 'conditional_chen_full_rosser_bridge_probe.json'
    bridge = json.loads(bridge_path.read_text(encoding='utf-8'))
    for name, expected in bridge['input_sha256'].items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == expected
    B = int(bridge['verified_goldbach_endpoint'])
    L = arb(B).log()
    G = arb(bridge['main_mass_lower_ball'])
    W = arb(result['full_lower_weighted_sum_upper_ball'])
    target = q(bridge['claimed_count_coefficient'])
    exceptional = (L**3 / arb(2).log() + L**2) / B
    sensitivities = []
    for epsilon in map(Fraction, ['1/1000', '1/400', '1/100', '1/50', '1/20']):
        charge_per_unit_A = W * L**2 * (-(q('1/2') - q(epsilon)) * L).exp()
        capacity = (L * G - exceptional - target) / charge_per_unit_A
        sensitivities.append({'epsilon': str(epsilon), 'A_capacity_ball': str(capacity),
                              'charge_per_unit_A_ball': str(charge_per_unit_A)})
    result['local_candidate_interval'] = bridge['candidate_analytic_interval']
    result['centered_pi_Montgomery_type_input'] = (
        '|pi(N;d,N)-pi(N)/phi(d)| <= A*N^(1/2+epsilon)/sqrt(d) '
        'for the actual nonzero lower Rosser moduli, for every even N in the local interval')
    result['A_capacity_sensitivity'] = sensitivities
    result['warning'] = ('This centered pi statement with numerical A and finite onset is not the '
                         'standard asymptotic corrected Montgomery conjecture for psi. '
                         'Deweighting, recentering, prime powers, constants and onset remain to be proved.')
    result['bridge_input_sha256'] = hashlib.sha256(bridge_path.read_bytes()).hexdigest()
    return result


if __name__ == '__main__':
    result = run()
    Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8', newline='\n')
    print(json.dumps({k: result[k] for k in ['full_lower_weighted_sum_upper_ball',
          'lower_medium_weighted_sum_ball', 'maximum_actual_degree', 'label_limit',
          'A_capacity_sensitivity', 'elapsed_seconds']}, indent=2))
