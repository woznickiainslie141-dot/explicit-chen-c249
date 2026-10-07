"""Fresh finite Arb certificate for a larger prefix-band exponent grid.

The weights, support lemma and external analytic lemmas are unchanged.
Only the majorants used to bound the finite positive first-failure sums
are enlarged. This does not certify those external analytic lemmas.
"""
import hashlib
import json
import math
import sys
from fractions import Fraction
from pathlib import Path

from flint import arb, ctx

from conditional_chen_prefix_band_certificate import certify_band
from conditional_chen_sharp_endpoint_certificate import grh_case, qG_upper, Umin


ROOT = Path(__file__).resolve().parent
FINITE = ROOT/'conditional_chen_dense_band_grh_finite.json'
ctx.prec = 192


def digest(name):
    return hashlib.sha256((ROOT/name).read_bytes()).hexdigest()


def upper_reciprocal(value):
    n = int(float(1/value))
    while not value < arb(1)/n:
        n -= 1
    assert value > arb(1)/(n+1)
    return f'1/{n}'


def endpoint(case, cost, support, L, claim='-1', beta=Fraction(7,2), tight_ap=False):
    beta = Fraction(beta)
    b = arb(beta.numerator)/beta.denominator
    x = arb(L)
    assert arb(7)/2 <= b < 4
    alpha = math.ceil(float((cost+b*x.log())/x)*10**6)
    while not arb(alpha)/10**6*x-b*x.log() > cost:
        alpha += 1
    raw = grh_case(case, cost, support, L, alpha, '1e-230', '1e-1450', '-1')
    assert x/2-b*x.log() > arb(10**9).log()
    assert arb(alpha)/10**6-b/x > 0
    # The lower and half-upper prime-sequence remainders charge 3R_A/2.
    # Counting square-free moduli gives R_A <= .608 * per_modulus * N/L^(beta-1).
    per_modulus = qG_upper(x)+1/(16*arb.pi())+2*(-x/2).exp()/arb(2).log()
    assert arb('.608')*per_modulus < arb('.12')
    leading = arb('.912')*per_modulus if tight_ap else arb('.18')
    leading_ap = leading*((arb(3)-b)*x.log()).exp()/Umin
    old_ap = arb(raw['AP_budget_ball'])
    new_ap = old_ap-arb('.18')/(Umin*x.sqrt())+leading_ap
    margin = arb(raw['margin_ball'])+old_ap-new_ap
    assert new_ap > 0 and margin > arb(claim)
    raw.update({'beta':str(beta),'tight_square_free_AP_coefficient':tight_ap,
                'AP_leading_coefficient_ball':str(leading),
                'AP_leading_expense_ball':str(leading_ap),
                'AP_budget_beta35_comparison_ball':str(old_ap),
                'AP_budget_ball':str(new_ap),'margin_ball':str(margin),
                'claimed_margin':claim,
                'support_slack_ball':str(arb(alpha)/10**6*x-b*x.log()-cost),
                'log_AP_level_at_start_ball':str(x/2-b*x.log()),
                'sieve_interface':'same defining finite weights and sharp support; enlarged majorant grid changes gap bounds; uniform product input reused'})
    return raw


def first_endpoint(case,cost,support,beta=Fraction(7,2),tight_ap=False):
    lo, hi = 30000, 31000
    claim = arb('.00001')
    assert arb(endpoint(case,cost,support,lo,beta=beta,tight_ap=tight_ap)['margin_ball']) < claim
    assert arb(endpoint(case,cost,support,hi,beta=beta,tight_ap=tight_ap)['margin_ball']) > claim
    while hi-lo > 1:
        mid = (lo+hi)//2
        if arb(endpoint(case,cost,support,mid,beta=beta,tight_ap=tight_ap)['margin_ball']) > claim:
            hi = mid
        else:
            lo = mid
    chosen = endpoint(case,cost,support,hi,'.00001',beta,tight_ap)
    preceding = endpoint(case,cost,support,hi-1,beta=beta,tight_ap=tight_ap)
    assert arb(preceding['margin_ball']) < claim
    return chosen,preceding


def direct_budget_check(row):
    """Recompute the modified AP charge directly, without old-budget subtraction."""
    x = arb(row['log_N0'])
    b = arb(row['beta'])
    alpha = arb(row['alpha'])
    epsilon = arb(row['epsilon'])
    eta = arb(row['eta_plus'])
    k = 8*(arb(1)/6-alpha)
    FM = 2*arb.const_euler().exp()/k+106*epsilon
    # The constant .912 is exactly (3/2)*(76/125).
    rho = arb(3)/2*(arb(76)/125)
    assert Fraction(3,2)*Fraction(76,125) == Fraction(114,125)
    assert abs(rho-arb('.912')) < arb('1e-55')
    expected = rho*(qG_upper(x)+1/(16*arb.pi())+2*(-x/2).exp()/arb(2).log())
    expected *= x**(3-b)/Umin
    expected += arb('.055')*x**4*(-x/6).exp()/Umin
    expected += arb('.8')*(1+eta)*(1+arb(1)/10**10)*FM*x**2*(-x/6).exp()
    assert abs(expected-arb(row['AP_budget_ball'])) < arb('1e-45')
    direct_margin = arb(row['lower_ball'])-arb(row['half_upper_ball'])
    direct_margin -= arb(row['half_switch_ball'])+expected
    direct_margin -= arb(row['half_rectangle_error_ball'])+arb(row['finite_error_ball'])
    assert abs(direct_margin-arb(row['margin_ball'])) < arb('1e-45')
    direct_slack = alpha*x-b*x.log()-arb(15).log()
    direct_slack -= arb(49)/8*arb(5000000).log()+arb('11/169').log()
    assert abs(direct_slack-arb(row['support_slack_ball'])) < arb('1e-45')
    assert x > 8500 and alpha-b/x > 0 and arb(1)/2-b/x > 0
    assert b > 3 and x > 240
    assert qG_upper(x) < qG_upper(arb(8500))
    return 'PASS: direct AP, surplus and support formulas agree within 1e-45'


def run(reuse=False):
    inputs = ['conditional_chen_prefix_band_certificate.py',
              'conditional_chen_rankin_certificate.py']
    generation_hashes = {p:digest(p) for p in inputs}
    if reuse:
        saved = json.loads(FINITE.read_text(encoding='utf-8'))
        assert saved['finite_generator_sha256'] == generation_hashes
        finite = saved['presieve']
    else:
        finite = certify_band(
            5000000, 5, Fraction(49,8),
            [Fraction(0),Fraction(1,10),Fraction(7,50),Fraction(9,50),
             Fraction(11,50),Fraction(7,25),Fraction(9,25)],
            [Fraction(a,25) for a in range(2,10)], M=12)
        FINITE.write_text(json.dumps({'presieve':finite,
                                     'finite_generator_sha256':generation_hashes},
                                    indent=2)+'\n',encoding='utf-8',newline='\n')
    sharp = json.loads((ROOT/'conditional_chen_sharp_endpoint_certificate.json')
                      .read_text(encoding='utf-8'))
    for name, value in sharp['input_sha256'].items():
        assert digest(name) == value, name
    old = sharp['GRH_positivity']
    original = old['presieve']
    assert finite['prime_generation'] == 'complete integer sieve'
    assert finite['logarithm_base'] == 'e'
    for key in ['T','exact_cutoff','u','exact_presieve_Q_max','medium_prime_count',
                'sigma_candidates','separate_degree_cutoff','degree_thresholds']:
        assert finite[key] == original[key], key
    assert set(original['band_sigma_candidates']) <= set(finite['band_sigma_candidates'])
    for key in ['upper_relative_Rosser_gap_ball','lower_relative_Rosser_gap_ball']:
        assert arb(finite[key]) < arb(original[key]), key
    uniform = json.loads((ROOT/'conditional_chen_uniform_product_certificate.json')
                        .read_text(encoding='utf-8'))
    row = next(r for r in uniform['rows'] if r['T'] == finite['T'])
    finite['uniform_epsilon_ball'] = row['ordered_uniform_epsilon_ball']
    case = {**old, 'presieve':finite,
            'eta_plus':upper_reciprocal(arb(finite['upper_relative_Rosser_gap_ball'])),
            'eta_minus':upper_reciprocal(arb(finite['lower_relative_Rosser_gap_ball'])),
            'epsilon':upper_reciprocal(arb(finite['uniform_epsilon_ball']))}
    cost = arb(finite['log_support_cost_ball'])+arb(old['sharp_support']['kappa']).log()
    assert abs(cost-arb(old['effective_log_support_cost_ball'])) < arb('1e-50')
    dense_only,_ = first_endpoint(case,cost,old['sharp_support'])
    assert dense_only['log_N0'] == 30779
    beta_grid = [Fraction(n,100) for n in range(350,358)]
    selected = [first_endpoint(case,cost,old['sharp_support'],b,True) for b in beta_grid]
    chosen,preceding = min(selected,key=lambda pair:(pair[0]['log_N0'],
                                                   -float(arb(pair[0]['margin_ball']))))
    at_31000 = max((endpoint(case,cost,old['sharp_support'],31000,beta=b,tight_ap=True)
                   for b in beta_grid),key=lambda r:float(arb(r['margin_ball'])))
    assert arb(at_31000['margin_ball']) > arb('.001')
    at_31000['claimed_margin'] = '.001'
    direct_checks = [direct_budget_check(row) for row in [chosen,at_31000]]
    for row,rect,finite_bound in [(chosen,'1e-248','1e-1659'),
                                  (at_31000,'1e-250','1e-1673')]:
        assert arb(row['half_rectangle_error_ball']) < arb(rect)
        assert arb(row['finite_error_ball']) < arb(finite_bound)
        row['rectangle_display_bound'] = rect
        row['finite_display_bound'] = finite_bound
    # Reproduce the preceding published-in-manuscript endpoint as a control.
    control = endpoint(old,cost,old['sharp_support'],30824,'.00001')
    assert control['alpha'] == old['alpha']
    for key in ['margin_ball','lower_ball','half_upper_ball','half_switch_ball','AP_budget_ball']:
        assert abs(arb(control[key])-arb(old[key])) < arb('1e-45'), key
    used = inputs+['conditional_chen_sharp_endpoint_certificate.py',
                  'conditional_chen_sharp_endpoint_certificate.json',
                  'conditional_chen_uniform_product_certificate.json',
                  FINITE.name,Path(__file__).name]
    result = {'status':'PASS: fresh finite prefix-band sums and strict GRH endpoint',
              'precision_bits':ctx.prec, 'finite_inputs_regenerated_in_this_run':not reuse,
              'uniform_product_inputs_regenerated_in_this_run':False,
              'only_changed_mathematical_input':'larger finite nonnegative majorant grid, tighter rounding, AP coefficient and beta budget',
              'GRH_positivity':chosen, 'GRH_at_31000':at_31000,
              'GRH_dense_only_beta35':dense_only,
              'beta_candidates':[str(b) for b in beta_grid],
              'beta_endpoint_rows':[{'beta':pair[0]['beta'],'log_N0':pair[0]['log_N0'],
                                    'alpha':pair[0]['alpha'],'margin_ball':pair[0]['margin_ball'],
                                    'preceding_margin_ball':pair[1]['margin_ball']}
                                   for pair in selected],
              'preceding_integer_recipe_margin_ball':preceding['margin_ball'],
              'preceding_integer_recipe_alpha':preceding['alpha'],
              'baseline_30824_scalar_reproduction':'PASS: five quantities agree within 1e-45',
              'direct_modified_budget_checks':direct_checks,
              'half_line_checks':'fixed alpha,beta,delta; increasing support slack and AP level; decreasing q_+ and all expenses',
              'input_sha256':{p:digest(p) for p in used},
              'not_checked':['external analytic lemmas','independent mathematical review',
                             'full splice','global optimality','PDF compilation']}
    return result


if __name__ == '__main__':
    result = run(reuse=sys.argv[1:] == ['--reuse-finite'])
    Path(__file__).with_suffix('.json').write_text(
        json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    for key in ['GRH_positivity','GRH_at_31000']:
        row = result[key]
        print(json.dumps({key:{k:row[k] for k in ['log_N0','alpha','beta','epsilon','eta_plus',
                                                'eta_minus','margin_ball','support_slack_ball']}},indent=2))
