"""Independent finite checks and exact counterexamples for bridge inputs.

Lucas certificates are verified by integer arithmetic; candidate primes
are not trusted merely because a primality library labels them prime.
"""
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
from flint import arb, ctx
from conditional_chen_weighted_rosser_support_certificate import certify
from conditional_chen_prefix_block_certificate import certify_blocks
from conditional_chen_prefix_block_checks import exact_gaps, exact_relaxed_sequence

ROOT = Path(__file__).resolve().parent
ctx.prec = 192
STATES = [(1, 0), (1, 1), (2, 0), (2, 1), (3, 0), (3, 1)]


def q(value):
    value = Fraction(value)
    return arb(value.numerator) / value.denominator


def close(a, b):
    assert abs(a - b) < arb('1e-43'), (str(a), str(b))


def checked_json(name):
    result = json.loads((ROOT / name).read_text(encoding='utf-8'))
    for path, expected in result.get('input_sha256', {}).items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == expected, path
    return result


def independent_primes(n):
    # Odd-only table, independent of the production full-integer sieve.
    flags = bytearray([1]) * ((n + 1) // 2)
    flags[0] = 0
    for p in range(3, math.isqrt(n) + 1, 2):
        if flags[p // 2]:
            for slot in range(p * p // 2, len(flags), p):
                flags[slot] = 0
    return [2 * i + 1 for i, flag in enumerate(flags) if flag]


def grid(T, ps, u, h):
    labels = {p: int((arb(p).log() / q(h)).floor().unique_fmpz()) for p in ps}
    logD = q(u) * arb(T).log()
    product, K = 1, 0
    for p in ps:
        product *= p
        if product**u.denominator < T**u.numerator:
            K += 1
        else:
            break
    blocks = []
    for label in sorted(set(labels.values()), reverse=True):
        group = sorted([p for p in ps if labels[p] == label], reverse=True)
        passing = int(((logD - 3 * arb(min(group)).log()) / q(h)).ceil().unique_fmpz())
        failure = max(0, int(((logD - 3 * arb(max(group)).log()) / q(h)).ceil().unique_fmpz()) - K + 1)
        blocks.append({'label': label, 'primes': group,
                       'passing_source_cutoff': passing, 'failure_source_cutoff': failure})
    J = int((logD / q(h)).ceil().unique_fmpz())
    return labels, blocks, K, J


def records(T, ps, u, h, K, J):
    labels, blocks, _, _ = grid(T, ps, u, h)
    cutoffs = {block['label']: block['passing_source_cutoff'] for block in blocks}
    rows = []
    for mask in range(1 << len(ps)):
        chosen = sorted([p for i, p in enumerate(ps) if mask >> i & 1], reverse=True)
        product, true, relaxed, B = 1, 3, 3, 0
        degrees = {}
        for m, p in enumerate(chosen, 1):
            flag = 1 if m % 2 else 2
            product *= p
            if (product * p * p)**u.denominator >= T**u.numerator:
                true &= 3 ^ flag
            if B >= cutoffs[labels[p]]:
                relaxed &= 3 ^ flag
            B += labels[p]
            degrees[labels[p]] = degrees.get(labels[p], 0) + 1
        keep = bool(relaxed) and B < J and all(k <= K for k in degrees.values())
        if true:
            assert keep and len(chosen) <= K
            assert true & relaxed == true
        rows.append((mask, product, B, (relaxed, len(chosen) % 2), keep, true))
    return rows


def weighted_hist(rows, available, sigma):
    hist = {s: {} for s in STATES}
    actual_lower = arb(0)
    for mask, product, B, state, keep, true in rows:
        if mask & ~available:
            continue
        weight = q(Fraction(1, product)) if sigma == 1 else arb(1) if sigma == 0 else 1 / arb(product).sqrt()
        if true & 2:
            actual_lower += weight
        if keep:
            hist[state][B] = hist[state].get(B, arb(0)) + weight
    return hist, actual_lower


def compare(polys, hist):
    count = 0
    for state, poly in zip(STATES, polys):
        for B in range(max(len(poly), max(hist[state], default=0) + 1)):
            close(poly[B], hist[state].get(B, arb(0)))
            count += 1
    return count


def small_checks():
    subsets = coefficients = deleted = patterns = 0
    configs = []
    for w in [2, 3]:
        ps = [p for p in independent_primes(29) if p > w]
        for u in [Fraction(3001, 1000), Fraction(61, 20)]:
            for h in [Fraction(1, 2), Fraction(1, 10)]:
                _, blocks, rawK, rawJ = grid(29, ps, u, h)
                gap = certify_blocks(29, w, u, h, False)
                rational = exact_relaxed_sequence(blocks, rawJ, rawK)
                for value, key in zip(rational, ['upper_relative_Rosser_gap_ball', 'lower_relative_Rosser_gap_ball']):
                    close(arb(gap[key]), q(value))
                for sigma in [Fraction(0), Fraction(1, 2), Fraction(1)]:
                    result, production = certify(29, w, u, sigma, h, return_states=True)
                    rows = records(29, ps, u, h, result['maximum_actual_degree'], result['label_limit'])
                    subsets += len(rows)
                    hist, actual = weighted_hist(rows, (1 << len(ps)) - 1, sigma)
                    coefficients += compare(production, hist)
                    assert actual < arb(result['lower_medium_weighted_sum_ball']) + arb('1e-43')
                    if u == Fraction(3001, 1000) and h == Fraction(1, 10) and sigma == 1:
                        for available in range(1 << len(ps)):
                            kept = [p for i, p in enumerate(ps) if available >> i & 1]
                            exact = exact_gaps(kept, 29, u)
                            assert all(0 <= a <= b for a, b in zip(exact, rational))
                            patterns += 1 << len(kept)
                            deleted += 1
                            fresh, polys = certify(29, w, u, sigma, h, kept, return_states=True)
                            reference, true_lower = weighted_hist(rows, available, sigma)
                            coefficients += compare(polys, reference)
                            assert true_lower < arb(fresh['lower_medium_weighted_sum_ball']) + arb('1e-43')
                            assert arb(fresh['lower_medium_weighted_sum_ball']) < arb(result['lower_medium_weighted_sum_ball']) + arb('1e-43')
                    configs.append({'w': w, 'u': str(u), 'h': str(h), 'sigma': str(sigma)})
    return {'configurations': configs, 'enumerated_complete_subsets': subsets,
            'weighted_histogram_coefficients_checked': coefficients,
            'deleted_production_histograms_checked': deleted,
            'pointwise_divisor_patterns_and_first_failure_identity_cases': patterns,
            'rational_weight_exponents': ['0', '1'], 'square_root_weight_exponent': '1/2'}


def prime_certificates(primes):
    import sympy
    nodes = {}
    def create(n):
        if str(n) in nodes:
            return
        if n <= 10000000:
            assert n >= 2 and all(n % r for r in range(2, math.isqrt(n) + 1))
            nodes[str(n)] = {'kind': 'trial division'}
            return
        factors = {int(p): int(e) for p, e in sympy.factorint(n - 1).items()}
        assert math.prod(p**e for p, e in factors.items()) == n - 1
        witnesses = {}
        for p in factors:
            create(p)
            for a in range(2, 1000):
                if pow(a, n - 1, n) == 1 and math.gcd(pow(a, (n - 1) // p, n) - 1, n) == 1:
                    witnesses[str(p)] = a
                    break
            else:
                raise AssertionError(('No Lucas witness', n, p))
        nodes[str(n)] = {'kind': 'full-factor Lucas', 'factors_of_n_minus_1': {str(p): e for p, e in factors.items()},
                         'witnesses': witnesses}
    for n in primes:
        create(n)
    # Independent verification uses only the certificate and integer arithmetic.
    proved = set()
    def verify(n):
        if n in proved:
            return
        node = nodes[str(n)]
        if node['kind'] == 'trial division':
            assert n >= 2 and all(n % r for r in range(2, math.isqrt(n) + 1))
        else:
            fs = {int(p): e for p, e in node['factors_of_n_minus_1'].items()}
            assert math.prod(p**e for p, e in fs.items()) == n - 1
            for p in fs:
                verify(p)
                a = node['witnesses'][str(p)]
                assert pow(a, n - 1, n) == 1
                assert math.gcd(pow(a, (n - 1) // p, n) - 1, n) == 1
        proved.add(n)
    for n in primes:
        verify(n)
    return nodes


def counterexamples(base, weighted, w2):
    cases = [
        {'w': 3, 'u': '3001/1000', 'exact': [3],
         'medium': [140053, 47, 41, 37, 29, 23, 19, 17, 11, 7, 5],
         'N': 4000000000000000048, 'primes': [1515229998243811363],
         'A_budget': weighted['A_capacity_sensitivity'][0]['A_capacity_ball']},
        {'w': 2, 'u': '3001/1000', 'exact': [],
         'medium': [873017, 43, 37, 31, 29, 23, 19, 17, 7, 5, 3],
         'N': 4000000000000000268, 'primes': [3025970232859408883, 1077910698578226113],
         'A_budget': w2['rows'][0]['A_capacity_epsilon_1_1000_ball']},
        {'w': 2, 'u': '61/20', 'exact': [],
         'medium': [1777351, 43, 37, 31, 29, 23, 19, 17, 7, 5, 3],
         'N': 4000000000000000016, 'primes': [2016999920211064361],
         'A_budget': w2['rows'][1]['A_capacity_epsilon_1_1000_ball']},
    ]
    to_prove = set()
    for case in cases:
        to_prove.update(case['exact'] + case['medium'] + case['primes'])
    nodes = prime_certificates(sorted(to_prove))
    for case in cases:
        u = Fraction(case['u'])
        N = case['N']
        factors = case['exact'] + case['medium']
        d = math.prod(factors)
        phi = math.prod(p - 1 for p in factors)
        assert len(set(factors)) == len(factors) and all(3 <= p <= 2000000 for p in factors)
        assert 4000000000000000000 <= N <= 2000000**3 and N % 2 == 0
        assert math.gcd(N, d) == 1
        prefix = 1
        for m, p in enumerate(case['medium'], 1):
            prefix *= p
            if m % 2 == 0:
                assert (prefix * p * p)**u.denominator < 2000000**u.numerator
        for p in case['primes']:
            assert 2 <= p < N and (N - p) % d == 0
        L = arb(N).log()
        mean_upper = q(Fraction(125506, 100000)) * N / (phi * L)
        error_lower = len(case['primes']) - mean_upper
        assert error_lower > 0
        A_min = error_lower * arb(d).sqrt() * (-(arb(1) / 2 + arb(1) / 1000) * L).exp()
        assert A_min > arb(case['A_budget'])
        # The actual modulus is within the proposed level and the literature's
        # q <= x^(1-epsilon) range for epsilon=1/1000.
        assert d**100 < N**99 and d**1000 < N**999
        case.update({'d': str(d), 'phi_d': str(phi), 'N': str(N),
                     'actual_coefficient': (-1)**len(factors),
                     'AP_prime_count_lower': len(case['primes']),
                     'AP_mean_upper_ball': str(mean_upper), 'AP_error_lower_ball': str(error_lower),
                     'necessary_A_epsilon_1_1000_lower_ball': str(A_min),
                     'candidate_A_capacity_ball': case.pop('A_budget'),
                     'status': 'REFUTED: every A within this candidate capacity fails at this actual supported modulus'})
    return cases, nodes


def run():
    base = checked_json('conditional_chen_full_rosser_bridge_probe.json')
    weighted = checked_json('conditional_chen_weighted_rosser_support_certificate.json')
    w2 = checked_json('conditional_chen_full_rosser_w2_probe.json')
    assert weighted['bridge_input_sha256'] == hashlib.sha256((ROOT / 'conditional_chen_full_rosser_bridge_probe.json').read_bytes()).hexdigest()
    primes = independent_primes(2000000)
    assert len(primes) == base['complete_odd_prime_count'] == 148932
    V = math.prod((1 - arb(1) / (p - 1) for p in primes), start=arb(1))
    close(V, arb(base['complete_density_ball']))
    configurations = [(base, 3, Fraction(3001, 1000))] + [(row, 2, Fraction(row['u'])) for row in w2['rows']]
    scalars = []
    for record, w, u in configurations:
        medium = [p for p in primes if p > w]
        kappa = max([Fraction(1, medium[0]**2)] + [Fraction(r, s * s) for r, s in zip(medium, medium[1:])])
        Q = math.prod(p for p in primes if p <= w)
        gap = record.get('finite_rosser_gap_certificate', record.get('gap'))
        G = V * (1 - arb(gap['lower_relative_Rosser_gap_ball']))
        key = 'main_mass_lower_ball' if w == 3 else 'main_mass_ball'
        close(G, arb(record[key]))
        B = 4000000000000000000
        L = arb(B).log()
        power = math.lcm(u.denominator, 100)
        assert (Q * kappa.numerator)**power * 2000000**(u.numerator * (power // u.denominator)) < (
            kappa.denominator**power * B**(99 * (power // 100)))
        capacity = (L * G - (L**3 / arb(2).log() + L**2) / B - arb(1) / 100000) * L**2
        close(capacity, arb(record['signed_AP_constant_capacity_ball']))
        support_record = weighted if w == 3 else record['weights']
        lower_states = sum((arb(row['weighted_sum_ball']) for row in support_record['state_weighted_balls']
                            if row['alive_flags'] & 2), arb(0))
        close(lower_states, arb(support_record['lower_medium_weighted_sum_ball']))
        presieve = math.prod((1 + 1 / arb(p).sqrt() for p in primes if p <= w), start=arb(1))
        W = presieve * lower_states
        close(W, arb(support_record['full_lower_weighted_sum_upper_ball']))
        support_degree, prefix_product = 0, 1
        for p in medium:
            prefix_product *= p
            if (prefix_product * kappa.denominator)**u.denominator < kappa.numerator**u.denominator * 2000000**u.numerator:
                support_degree += 1
            else:
                break
        assert support_degree == support_record['maximum_actual_degree']
        logH = q(u) * arb(2000000).log() + q(kappa).log()
        assert int((logH / q(Fraction(1, 100))).ceil().unique_fmpz()) == support_record['label_limit']
        epsilons = [Fraction(r['epsilon']) for r in weighted['A_capacity_sensitivity']] if w == 3 else [Fraction(1, 1000)]
        for epsilon in epsilons:
            A_capacity = capacity / L**2 / (W * L**2 * (-(arb(1) / 2 - q(epsilon)) * L).exp())
            expected = next(r['A_capacity_ball'] for r in weighted['A_capacity_sensitivity']
                            if Fraction(r['epsilon']) == epsilon) if w == 3 else record['A_capacity_epsilon_1_1000_ball']
            close(A_capacity, arb(expected))
        C = record['signed_AP_constant_budget']
        assert capacity > C and capacity < C + 1
        for N in [B, 5000000000000000000, 6000000000000000000, 8000000000000000000]:
            ln = arb(N).log()
            margin = ln * G - C / ln**2 - (ln**3 / arb(2).log() + ln**2) / N
            assert margin > arb(1) / 100000
        scalars.append({'w': w, 'u': str(u), 'C_budget': C, 'kappa': str(kappa)})
    small = small_checks()
    counter, certificates = counterexamples(base, weighted, w2)
    names = [Path(__file__).name, 'conditional_chen_full_rosser_bridge_probe.py',
             'conditional_chen_full_rosser_bridge_probe.json',
             'conditional_chen_weighted_rosser_support_certificate.py',
             'conditional_chen_weighted_rosser_support_certificate.json',
             'conditional_chen_full_rosser_w2_probe.py', 'conditional_chen_full_rosser_w2_probe.json',
             'conditional_chen_prefix_block_checks.py']
    return {'status': 'PASS: independent finite checks and rigorous numerical-input refutations',
            'precision_bits': ctx.prec, 'complete_odd_prime_count': len(primes),
            'local_scalar_reproductions': scalars, 'small_exact_checks': small,
            'centered_pi_input_counterexamples': counter, 'primality_certificate_nodes': certificates,
            'continuous_interval_proof': 'In the separate research note, not supplied by finite sampling',
            'input_sha256': {n: hashlib.sha256((ROOT / n).read_bytes()).hexdigest() for n in names},
            'not_checked': ['unproved signed AP input', 'independent general proof review',
                            'finite constants implied by standard named conjectures', 'coverage above T^3', 'full splice']}


if __name__ == '__main__':
    result = run()
    Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8', newline='\n')
    print(json.dumps({k: result[k] for k in ['status', 'complete_odd_prime_count', 'local_scalar_reproductions',
                                          'small_exact_checks', 'centered_pi_input_counterexamples']}, indent=2))
