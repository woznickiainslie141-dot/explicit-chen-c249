"""Finite Rosser prefix certificates for u>2, including nonpositive cuts.

The raw first-failure prefixes are <T^u for u>2. The accepted support is
bounded by max(T,kappa*T^u), accounting for lower single-prime weights.
These facts require the general proof in the associated research note.
No arithmetic-progression estimate is proved by this program.
"""
import hashlib
import json
import math
import time
from fractions import Fraction
from pathlib import Path
from flint import arb, arb_poly, ctx
from conditional_chen_rankin_certificate import odd_primes_to
from conditional_chen_prefix_block_certificate import integer_floor, integer_ceil, block_coefficients, action
from conditional_chen_sharp_support_certificate import support_ratio

ROOT = Path(__file__).resolve().parent


def q(value):
    value = Fraction(value)
    return arb(value.numerator) / value.denominator


def source_cut(T, u, p, h):
    # At u=3 and p=T (prime), the logarithmic expression is exactly zero.
    # Computing 3*log(T)-3*log(T) in interval arithmetic loses that equality.
    if u == 3 and p == T:
        return 0
    return integer_ceil((q(u) * arb(T).log() - 3 * arb(p).log()) / q(h))


def grid(T, w, u, h):
    u, h = map(Fraction, [u, h])
    assert T >= 3 and 2 <= w < T and u > 2 and h > 0
    primes = odd_primes_to(T)
    medium = [p for p in primes if p > w]
    assert medium
    logD = q(u) * arb(T).log()
    J = integer_ceil(logD / q(h))
    K, product = 0, 1
    for p in medium:
        product *= p
        if product**u.denominator < T**u.numerator:
            K += 1
        else:
            break
    assert K >= 1
    groups = {}
    for p in reversed(medium):
        groups.setdefault(integer_floor(arb(p).log() / q(h)), []).append(p)
    blocks = []
    for label, ps in groups.items():
        passing = max(0, source_cut(T, u, min(ps), h))
        failure = max(0, source_cut(T, u, max(ps), h) - K + 1)
        assert 0 <= failure <= passing and label > 0
        blocks.append({'label': label, 'primes': ps, 'passing_source_cutoff': passing,
                       'failure_source_cutoff': failure})
    return primes, medium, blocks, K, J


def certify(T, w, u, h=Fraction(1, 200), available=None, progress=False):
    assert ctx.prec >= 192
    u, h = map(Fraction, [u, h])
    started = time.monotonic()
    primes, medium, blocks, K, J = grid(T, w, u, h)
    allowed = set(medium) if available is None else set(available)
    assert allowed <= set(medium)
    states = [arb_poly([1]), arb_poly([]), arb_poly([1]), arb_poly([])]
    upper, lower = arb(0), arb(0)
    for number, block in enumerate(blocks, 1):
        selected = [p for p in block['primes'] if p in allowed]
        d, c = block_coefficients(selected, K)
        powers = states
        updated = [arb_poly([]) for _ in range(4)]
        for j in range(K + 1):
            for index in range(4):
                updated[index] += d[j] * powers[index]
            failure = block['failure_source_cutoff']
            upper += c[j] * powers[0].right_shift(failure)(arb(1))
            lower += c[j] * powers[3].right_shift(failure)(arb(1))
            if j == K or not any(powers):
                break
            powers = action(powers, block['label'], block['passing_source_cutoff'], J)
        states = updated
        if progress and number % 200 == 0:
            print(f'Subcubic T={T}, u={u}: {number}/{len(blocks)} blocks; '
                  f'{time.monotonic()-started:.1f}s', flush=True)
    kappa, witness = support_ratio(medium)
    logD = q(u) * arb(T).log()
    kappa_dominates_T = kappa.numerator**u.denominator * T**u.numerator > (
        kappa.denominator**u.denominator * T**u.denominator)
    logH = logD + q(kappa).log() if kappa_dominates_T else arb(T).log()
    exact = [p for p in primes if p <= w]
    Q = math.prod(exact)
    names = [Path(__file__).name, 'conditional_chen_rankin_certificate.py',
             'conditional_chen_prefix_block_certificate.py', 'conditional_chen_sharp_support_certificate.py']
    return {
        'status': 'strict finite u>2 prefix certificate; AP remainder and general proof review unproved',
        'precision_bits': ctx.prec, 'T': T, 'exact_cutoff': w, 'u': str(u), 'h': str(h),
        'Q': Q, 'kappa': str(kappa), 'support_maximizer': witness,
        'support_formula': 'Q*max(T,kappa*T^u)', 'kappa_D_dominates_T': kappa_dominates_T,
        'log_medium_support_upper_ball': str(logH),
        'log_full_support_upper_ball': str(arb(Q).log() + logH),
        'maximum_first_failure_prefix_degree': K, 'label_limit': J,
        'nonempty_blocks': len(blocks), 'zero_passing_cutoff_blocks': sum(b['passing_source_cutoff']==0 for b in blocks),
        'complete_odd_prime_count': len(primes), 'complete_medium_prime_count': len(medium),
        'available_medium_prime_count': len(allowed),
        'upper_relative_Rosser_gap_ball': str(upper), 'lower_relative_Rosser_gap_ball': str(lower),
        'complete_cutoffs_retained_after_deletion': True,
        'elapsed_seconds': round(time.monotonic()-started, 3),
        'input_sha256': {n: hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names},
    }


def run():
    ctx.prec = 192
    B = 4000000000000000000
    theta = Fraction(99, 100)
    rows = []
    left = B
    families = [(2000000,Fraction(61,20)), (3000000,Fraction(3)),
                (6000000,Fraction(297,100)), (12000000,Fraction(297,100))]
    for T,u in families:
        finite = certify(T, 2, u, progress=True)
        assert finite['kappa_D_dominates_T'] and finite['Q'] == 1 and finite['kappa'] == '3/25'
        # theta=u/3 makes the support comparison independent of scale
        # when the upper sieve cutoff doubles between adjacent intervals.
        power = math.lcm(u.denominator,theta.denominator)
        assert 3**power * T**(u.numerator*(power//u.denominator)) < (
            25**power * left**(theta.numerator*(power//theta.denominator)))
        primes = odd_primes_to(T)
        V = math.prod((1 - arb(1)/(p-1) for p in primes), start=arb(1))
        G = V * (1 - arb(finite['lower_relative_Rosser_gap_ball']))
        assert G > 0
        L = arb(left).log()
        loss = (L**3/arb(2).log() + L**2)/left
        capacity = (L*G - loss - arb(1)/100000)*L**2
        budget = math.floor(float(capacity))
        while not capacity > budget:
            budget -= 1
        rows.append({'left_endpoint': str(left), 'right_endpoint': str(T**3), 'finite_certificate': finite,
                     'complete_density_ball': str(V), 'main_mass_lower_ball': str(G),
                     'signed_C_capacity_ball': str(capacity), 'signed_C_integer_budget': budget,
                     'support_slack_at_left_ball': str(q(theta)*L-arb(finite['log_full_support_upper_ball']))})
        print({k: rows[-1][k] for k in ['left_endpoint','right_endpoint','main_mass_lower_ball','signed_C_integer_budget']}, flush=True)
        left = T**3
    C = min(row['signed_C_integer_budget'] for row in rows)
    for row in rows:
        L = arb(int(row['left_endpoint'])).log()
        G = arb(row['main_mass_lower_ball'])
        margin = L*G-C/L**2-(L**3/arb(2).log()+L**2)/int(row['left_endpoint'])
        assert margin > arb(1)/100000
        row['common_C_normalized_margin_at_left_ball'] = str(margin)
    return {'status': 'PASS: conditional finite coverage; signed AP premise unproved',
            'precision_bits': ctx.prec, 'theta': str(theta), 'tail_u': '297/100', 'rows': rows,
            'common_signed_AP_constant_budget': C, 'claimed_count_coefficient': '1/100000',
            'conditional_analytic_coverage': [str(B),str(left)],
            'support_scaling_exact_test': '3^100*2^297 < 25^100',
            'required_AP_premise': 'For each of the four fixed-T families, R_W(N)>=-C*N/log(N)^4 for every even N in its stated interval.',
            'not_checked': ['signed AP premise', 'independent reproduction', 'independent general proof review',
                            'unbounded-T main-mass control', 'coverage above 12000000^3', 'full splice']}


if __name__ == '__main__':
    result = run()
    result['input_sha256'] = {Path(__file__).name: hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:result[k] for k in ['status','common_signed_AP_constant_budget','conditional_analytic_coverage']},indent=2))
