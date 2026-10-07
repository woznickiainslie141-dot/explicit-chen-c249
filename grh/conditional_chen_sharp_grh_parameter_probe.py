"""FLOAT ONLY: compare prefix-band parameters using the sharp support cost.

This search changes no manuscript claim. Newton subtractions, moment
cancellation and clipped parity tails require fresh Arb certification.
"""
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

from flint import arb

from conditional_chen_prefix_band_probe import band_gaps
from conditional_chen_rankin_certificate import odd_primes_to
from conditional_chen_sharp_support_certificate import support_ratio
from conditional_chen_uniform_grh_probe import surplus


ROOT = Path(__file__).resolve().parent


def first_float_endpoint(row, epsilon, claim=1e-5):
    lo, hi = 20000, 40000
    if surplus(hi, row['sharp_log_support_cost'], epsilon,
               row['upper_gap'], row['lower_gap'])[0] <= claim:
        return None
    while hi-lo > 1:
        mid = (lo+hi)//2
        if surplus(mid, row['sharp_log_support_cost'], epsilon,
                   row['upper_gap'], row['lower_gap'])[0] > claim:
            hi = mid
        else:
            lo = mid
    return hi


def run():
    uniform = json.loads((ROOT/'conditional_chen_uniform_product_certificate.json')
                         .read_text(encoding='utf-8'))
    uniform_rows = {r['T']: r for r in uniform['rows']}
    sigmas = [Fraction(a, 50) for a in range(19)]
    us = [Fraction(a, 16) for a in range(96, 101)]
    records = []
    for T in [5000000]:
        rows = band_gaps(T, us, sigmas, [3, 5, 7], M=12)
        primes = odd_primes_to(T)
        epsilon = float(arb(uniform_rows[T]['ordered_uniform_epsilon_ball']))*1.000001
        for row in rows:
            kappa, witness = support_ratio([p for p in primes if p > row['w']])
            row['kappa'] = str(kappa)
            row['kappa_witness'] = witness
            row['sharp_log_support_cost'] = row['log_support_cost']+math.log(kappa)
            row['epsilon_float'] = epsilon
            row['float_log_start'] = first_float_endpoint(row, epsilon)
            row['margin_at_30824_float'], row['alpha_at_30824_float'] = surplus(
                30824, row['sharp_log_support_cost'], epsilon,
                row['upper_gap'], row['lower_gap'])
            records.append(row)
    inputs = ['conditional_chen_uniform_product_certificate.json',
              'conditional_chen_prefix_band_probe.py',
              'conditional_chen_uniform_grh_probe.py',
              'conditional_chen_sharp_support_certificate.py', Path(__file__).name]
    return {'status':'FLOAT ONLY: proposals for a fresh Arb run; no new theorem',
            'sigma_candidates':[str(s) for s in sigmas],
            'u_candidates':[str(u) for u in us],
            'records':sorted(records, key=lambda r:(r['float_log_start'] or 10**9,
                                                   -r['margin_at_30824_float'])),
            'input_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
                            for p in inputs}}


if __name__ == '__main__':
    result = run()
    Path(__file__).with_suffix('.json').write_text(
        json.dumps(result, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps(result['records'][:5], indent=2))
