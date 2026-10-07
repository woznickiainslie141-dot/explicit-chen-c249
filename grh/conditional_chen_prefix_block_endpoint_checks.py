"""Independent direct formulas and half-line range checks for current endpoints.

Does not call either production scalar evaluator. The finite domination
lemma and external analytic inputs remain mathematical dependencies.
"""
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from flint import arb, acb, ctx

ROOT=Path(__file__).resolve().parent
ctx.prec=192


def q(value):
    value=Fraction(value)
    return arb(value.numerator)/value.denominator


def check_finite(case):
    pre=case['presieve']
    assert pre['prime_generation']=='complete integer sieve'
    assert pre['retains_all_preceding_checked_prefixes'] is True
    assert arb(pre['uniform_epsilon_ball'])<q(case['epsilon'])
    assert arb(pre['upper_relative_Rosser_gap_ball'])<q(case['eta_plus'])
    assert arb(pre['lower_relative_Rosser_gap_ball'])<q(case['eta_minus'])
    support=case['sharp_support']
    cost=arb(support['Q']).log()+q(pre['u'])*arb(pre['T']).log()+q(support['kappa']).log()
    assert abs(cost-arb(case['effective_log_support_cost_ball']))<arb('1e-45')
    for name,expected in pre['input_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected,name
    return cost


def grh_direct(row):
    cost=check_finite(row)
    x=arb(row['log_N0']); a=q(row['alpha']); b=q(row['beta'])
    e=q(row['epsilon']); ep=q(row['eta_plus']); em=q(row['eta_minus'])
    eg=arb.const_euler().exp(); U=arb('1.32')/eg
    t=arb(1)/2-a; s=8*t; k=8*(t-arb(1)/3)
    v=q('1/10000000000'); delta=q('1/100000000'); ell=(1+delta).log()
    assert row['C1']==106 and row['C2']==108 and 3<s<4 and k>1
    assert x>=8500 and a-b/x>0 and arb(1)/2-b/x>0 and b>=3
    assert x/2-b*x.log()>arb(10**9).log()
    assert x/8>arb(row['presieve']['T']).log() and x>240
    slack=a*x-b*x.log()-cost
    assert slack>0
    integral=acb.integral(lambda y,analytic:(y-2).log(analytic=analytic)/(y-1),acb(3),acb(s)).real
    F=2*eg/s*(1+integral)
    f=2*eg/s*(s-1).log(); h=3*(-s).exp()/s
    factor=(1+ep)*f-(ep+em)*F
    factor-=e*arb(2).exp()*h*((1+ep)*108+(ep+em)*106)
    lower=8*(1-v)*factor
    I=eg/(4*t)*(6*(3-8*a)/(3-18*a)).log()
    FM=2*eg/k+106*e; J=q('8/3').log(); b0=t-q('1/4')
    Hh=acb.integral(lambda y,analytic:(2-8*(acb(t)-y)).exp()/y,acb(q('1/8')),acb(b0)).real
    Hh+=(q('1/3')/b0).log()
    upper=4*(1+ep)*(1+2/x)*(1+v)/(1-(-x/8).exp())*(I+106*e*Hh+2*FM/x**2)
    cb=acb.integral(lambda y,analytic:(2-3*y).log(analytic=analytic)/(y*(1-y)),acb(q('1/8')),acb(q('1/3'))).real
    switch=(1+ep)*(1+v)*(eg/t+159*e)*(1+delta)*(1+4/x)
    switch*=cb+2/x**2+(J+2/x**2)*(6*ell/x+8/x**2)
    Q=q('165/1000')+q('12683/1000')/x+q('254980/1000')/x**2
    Q+=q('2607854/1000')/x**3+q('11605056/1000')/x**4
    Q+=(q('1314/1000')*x+q('92/1000')*x.log()+q('60883/1000')
        +q('8250/1000')*x.log()/x+q('939260/1000')/x)*(-x/4).exp()
    multiplier=q('114/125')
    if 'sparse_modulus_count_budget' in row:
        count=row['sparse_modulus_count']; chosen=count['selected']
        fraction=arb(count['exact_divisor_count_max'])/row['sharp_support']['Q']
        fraction*=((q(chosen['sigma'])-1)*arb(count['medium_support_log_ball'])
                   +arb(chosen['log_euler_product_ball'])).exp()
        assert abs(fraction-arb(chosen['count_fraction_ball']))<arb('1e-42')
        assert 0<fraction<q(row['sparse_modulus_count_budget'])<arb('.608')
        assert count['exact_divisor_count_max']==2**count['exact_prime_count']
        multiplier=q('3/2')*q(row['sparse_modulus_count_budget'])
    lead=multiplier*(Q+1/(16*arb.pi())+2*(-x/2).exp()/arb(2).log())
    ap=lead*((3-b)*x.log()).exp()/U+q('11/200')*x**4*(-x/6).exp()/U
    ap+=q('4/5')*(1+ep)*(1+v)*FM*x**2*(-x/6).exp()
    rectangle=(1+delta)*x**5/(2*U*ell)*(q('3/4')*(-x/48).exp()+(-x/8).exp())
    finite=x**2/U*(2*(-x/8).exp()+(-2*x/3).exp()+(-x).exp())
    values={'lower_ball':lower,'half_upper_ball':upper,'half_switch_ball':switch,
            'AP_budget_ball':ap,'half_rectangle_error_ball':rectangle,
            'finite_error_ball':finite,'support_slack_ball':slack,
            'margin_ball':lower-upper-switch-ap-rectangle-finite}
    for key,value in values.items():assert abs(value-arb(row[key]))<arb('1e-43'),key
    assert values['margin_ball']>q(row['claimed_margin'])
    return {'log_N0':row['log_N0'],'alpha':row['alpha'],'beta':row['beta'],
            'quantities_independently_reproduced':len(values),
            'fixed_parameter_half_line_ranges_checked':True}


def distribution_direct(row):
    case=row['case']; cost=check_finite(case)
    x=arb(row['log_start']); theta=q(row['theta']); s=3*theta-3*cost/x
    assert x>3*arb(10**12).log() and x/3>arb(case['presieve']['T']).log() and 2<s<3
    eps=q(case['epsilon']); ep=q(case['eta_plus']); em=q(case['eta_minus'])
    eg=arb.const_euler().exp(); U=arb('1.32')/eg
    main=(q('198/25')/s*((1+ep)*(s-1).log()-ep-em)
          -3*U*eps*(2-s).exp()*((1+ep)*case['C2']+(ep+em)*case['C1']))
    rho=1-q('20508/1000000000')/x-2*x*(-x/3).exp()
    main*=rho
    deleted=(x**3/arb(2).log()+x**2)*(-x).exp()
    margin=main-arb(row['C_upper_budget'])/x**2-deleted
    assert abs(main-arb(row['main_ball']))<arb('1e-45')
    assert abs(margin-arb(row['margin_ball']))<arb('1e-45')
    assert margin>q('1/100000') and main>0 and 0<rho<1
    # l'(s)>0 and u'(s)<0 on 2<s<3; s increases with log N.
    fprime=2*eg*(s/(s-1)-(s-1).log())/s**2
    lprime=fprime+case['C2']*eps*(2-s).exp()
    uprime=-2*eg/s**2-case['C1']*eps*(2-s).exp()
    assert fprime>0 and lprime>0 and uprime<0 and cost>0
    return {'theta':row['theta'],'C_upper_budget':row['C_upper_budget'],
            'log_start':row['log_start'],'independent_main_and_margin':True,
            'half_line_derivative_signs_at_endpoint':True}


def run():
    refined=ROOT/'conditional_chen_prefix_block_refinement.json'
    path=refined if refined.exists() else ROOT/'conditional_chen_prefix_block_endpoint_certificate.json'
    distribution=json.loads(path.read_text(encoding='utf-8'))
    sparse=ROOT/'conditional_chen_sparse_modulus_certificate.json'
    result=json.loads(sparse.read_text(encoding='utf-8')) if sparse.exists() else distribution
    for report in [distribution,result]:
        for name,expected in report['input_sha256'].items():
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected,name
    grh=[grh_direct(result['GRH_positivity'])]
    if 'GRH_at_31000' in result:grh.append(grh_direct(result['GRH_at_31000']))
    direct=[distribution_direct(row) for row in distribution['distribution_endpoints']]
    inputs=[path.name,Path(__file__).name]+([sparse.name] if sparse.exists() else [])
    return {'status':'PASS: independent direct formulas and endpoint domain checks',
            'GRH_checks':grh,'distribution_checks':direct,
            'input_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs},
            'not_checked':['general prefix-block domination and support-count proofs','external analytic inputs',
                           'independent mathematical review','global optimality','full splice']}


if __name__=='__main__':
    r=run()
    Path(__file__).with_suffix('.json').write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(r,indent=2))
