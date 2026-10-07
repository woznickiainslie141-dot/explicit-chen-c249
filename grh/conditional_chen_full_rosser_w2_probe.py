"""All-prime Rosser probes without the exact-small-prime support cost."""
from pathlib import Path
from fractions import Fraction
from flint import arb, ctx
import json
import hashlib
import math
from conditional_chen_prefix_block_certificate import certify_blocks
from conditional_chen_weighted_rosser_support_certificate import certify, q

ROOT = Path(__file__).resolve().parent


def run():
    ctx.prec = 192
    B = 4000000000000000000
    L = arb(B).log()
    base_path = ROOT / 'conditional_chen_full_rosser_bridge_probe.json'
    base = json.loads(base_path.read_text(encoding='utf-8'))
    for name, expected in base['input_sha256'].items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == expected
    V = arb(base['complete_density_ball'])
    rows = []
    for u in [Fraction(3001, 1000), Fraction(61, 20)]:
        gap = certify_blocks(2000000, 2, u, Fraction(1, 200), True)
        weights = certify(2000000, 2, u, progress=True)
        G = V * (1 - arb(gap['lower_relative_Rosser_gap_ball']))
        W = arb(weights['full_lower_weighted_sum_upper_ball'])
        loss = (L**3 / arb(2).log() + L**2) / B
        charge = W * L**2 * (-(arb(1) / 2 - arb(1) / 1000) * L).exp()
        capacity = (L * G - loss - arb(1) / 100000) / charge
        signed_capacity = (L * G - loss - arb(1) / 100000) * L**2
        C = math.floor(float(signed_capacity))
        while not signed_capacity > C:
            C -= 1
        log_support = q(u) * arb(2000000).log() + q(Fraction(3, 25)).log()
        slack = q(Fraction(99, 100)) * L - log_support
        assert slack > 0 and C > 0 and G > 0
        row = {'u': str(u), 'w': 2, 'T': 2000000, 'Q': 1, 'kappa': '3/25',
               'gap': gap, 'weights': weights, 'main_mass_ball': str(G),
               'A_capacity_epsilon_1_1000_ball': str(capacity),
               'signed_AP_constant_capacity_ball': str(signed_capacity),
               'signed_AP_constant_budget': C, 'log_support_upper_ball': str(log_support),
               'support_slack_ball': str(slack),
               'candidate_analytic_interval': [str(B), str(2000000**3)]}
        rows.append(row)
        print({k: row[k] for k in ['u', 'main_mass_ball', 'A_capacity_epsilon_1_1000_ball',
                                 'signed_AP_constant_budget', 'support_slack_ball']}, flush=True)
    names = [Path(__file__).name, 'conditional_chen_weighted_rosser_support_certificate.py',
             'conditional_chen_prefix_block_certificate.py', 'conditional_chen_full_rosser_bridge_probe.json']
    return {'status': 'finite parameter probes; AP input unproved', 'precision_bits': ctx.prec,
            'rows': rows, 'input_sha256': {n: hashlib.sha256((ROOT / n).read_bytes()).hexdigest() for n in names},
            'not_checked': ['independent reproduction', 'signed AP input', 'standard conjecture implication',
                            'coverage above T^3', 'full splice']}


if __name__ == '__main__':
    report = run()
    Path(__file__).with_suffix('.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8', newline='\n')
