"""Regenerate the selected finite presieve, product and seven-eighths endpoints.

Analytic proofs are supplied separately: a successful run certifies finite
arithmetic, never the correctness of an arbitrary analytic argument.
"""
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from flint import arb, ctx
import qrh_chen_certificate as q
from conditional_chen_prefix_block_certificate import certify_blocks
from conditional_chen_uniform_product_certificate import certify_uniform_products
from qrh_chen_interface_checks import run as check_interfaces

ctx.prec=192
ROOT=Path(__file__).resolve().parent


def run():
    case=q.selected_case()
    print('Checking new interface constants and complete prime-support input',flush=True)
    interfaces=check_interfaces()
    print('Regenerating all 348510 medium-prime masses',flush=True)
    pre=certify_blocks(5000000,5,Fraction(11,2),Fraction(1,200),progress=True)
    assert arb(pre['upper_relative_Rosser_gap_ball'])<arb(case['eta_plus'])
    assert arb(pre['lower_relative_Rosser_gap_ball'])<arb(case['eta_minus'])
    assert arb(pre['upper_relative_Rosser_gap_ball'])<arb('.000309498615')
    assert arb(pre['lower_relative_Rosser_gap_ball'])<arb('.000308430014')
    assert (pre['medium_prime_count'],pre['maximum_actual_prefix_degree'],
            pre['bin_limit'],pre['nonempty_prime_blocks'])==(348510,22,16968,1810)
    print('Regenerating the unconditional ordered-product bound',flush=True)
    uniform=certify_uniform_products(cutoffs=(5000000,),progress=True)
    assert arb(uniform['rows'][0]['ordered_uniform_epsilon_ball'])<arb(case['epsilon'])
    assert arb(uniform['rows'][0]['ordered_uniform_epsilon_ball'])<arb('.000024260536')
    assert uniform['complete_odd_prime_count']==5761454
    kappa=arb(11)/169
    cost=arb(pre['log_support_cost_ball'])+kappa.log()
    assert abs(cost-arb(case['effective_log_support_cost_ball']))<arb('1e-45')
    # Fixed rational parameters for a prospective whole-half-line theorem.
    fixed=q.endpoint(case,51200,Fraction(1,1000000),Fraction(4315,1000000))
    assert arb(fixed['margin'])>arb('.00003')
    assert arb(fixed['support_slack_ball'])>arb('.57')
    positivity=q.endpoint(case,51046,Fraction(1,1000000),Fraction(2158047,500000000))
    assert arb(positivity['margin'])>arb('9e-7')
    assert arb(positivity['support_slack_ball'])>arb('4.0153e-5')
    assert arb(positivity['c4_times_L2'])<arb('2.200019')
    assert interfaces['log_N0']==51046
    for key,upper in [('A_lower_gate','1e-1013'),
                      ('A_upper_gate','.000058771'),
                      ('pi_upper_gate','.000528935'),
                      ('all_small_primitive_conductors_gate','1.328e-938')]:
        assert arb(interfaces[key])<arb(upper),key
    # The coarse decimal intervals printed in Proposition 7.1 are acceptance
    # gates too, not merely rounded display strings.
    paper_bounds={
        'lower':('7.806926624899798362','7.806926624899798363'),
        'half_upper':('6.500392797042414805','6.500392797042414806'),
        'half_switch':('1.306527039983525481','1.306527039983525482'),
        'AP':('3.1474e-10','3.1475e-10'),
        'rectangle':('5.8028784809024e-6','5.8028784809025e-6'),
        'margin':('9.8468e-7','9.8469e-7'),
    }
    for key,(lower,upper) in paper_bounds.items():
        assert arb(lower)<arb(positivity[key])<arb(upper),key
    assert arb(positivity['rectangle_deletions'])+arb(positivity['finite'])<arb('1e-2700')
    # Domain for every relative totient sum in the prime-factor row.
    L=arb(51046)
    w=L-5*L.log()
    assert w/2-10*w.log()-L/3>arb(10**9).log()
    names=['qrh_chen_certificate.py',Path(__file__).name,'qrh_chen_interface_checks.py',
           'qrh_chen_research_20261007.md',
           'conditional_chen_prefix_block_certificate.py','conditional_chen_rankin_certificate.py',
           'conditional_chen_uniform_product_certificate.py','qrh_chen_51046.tex',
           'bjs_source/Explicit_Chen_-_New.tex','qrh_source/paper.tex','qrh_source/Nonvanishing.lean']
    result={'status':'PASS: fresh finite masses, full ordered product, two fixed rational endpoints',
        'finite_presieve':pre,'uniform_product':uniform,'fixed_endpoint':fixed,
        'new_interface_gates':interfaces,
        'positivity_endpoint':positivity,
        'selected_rational_inputs':case,
        'dependency_scope':'No old GRH JSON certificates or old GRH manuscript are read',
        'input_sha256':{n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names},
        'not_checked':['external seven-eighths theorem proof','general new analytic extension',
                       'independent research review','global optimality']}
    (ROOT/'qrh_chen_reproduction.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('FRESH QRH FINITE REPRODUCTION PASSED',flush=True)
    print(json.dumps(fixed,indent=2),flush=True)


if __name__=='__main__':run()
