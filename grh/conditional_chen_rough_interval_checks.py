"""Independent scalar formulas, complete integer density and analytic join.

No production endpoint evaluator is called. Finite histogram checks and
general domination proofs are separate from these scalar verifications.
"""
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
from flint import arb,acb,ctx
from conditional_chen_prefix_block_endpoint_checks import check_finite,grh_direct

ROOT=Path(__file__).resolve().parent
ctx.prec=192


def q(value):
    f=Fraction(value)
    return arb(f.numerator)/f.denominator


def primes_direct(limit):
    sieve=bytearray(b'\1')*(limit+1)
    sieve[:2]=b'\0\0'
    for p in range(2,math.isqrt(limit)+1):
        if sieve[p]:
            first=p*p
            sieve[first:limit+1:p]=b'\0'*((limit-first)//p+1)
    return [p for p in range(2,limit+1) if sieve[p]]


def direct(row,L,anchor,V):
    cost=check_finite(row)
    x=arb(L);x0=arb(row['log_N0']);A=arb(anchor)
    alpha=q(row['alpha']);beta=q(row['beta']);t=arb(1)/2-alpha
    eps=q(row['epsilon']);ep=q(row['eta_plus']);em=q(row['eta_minus'])
    eg=arb.const_euler().exp();U=q('33/25')/eg
    s=8*t;k=8*(t-q('1/3'));v=q('1/10000000000');delta=q('1/100000000')
    ell=(1+delta).log()
    assert row['C1']==106 and row['C2']==108 and 3<s<4 and k>1
    assert 0<beta<=3 and x0<=x<=A and x0>=8500
    assert alpha-beta/x0>0 and q('1/2')-beta/x0>0
    assert x0/2-beta*x0.log()>arb(10**9).log()
    assert x0/6-beta*x0.log()>0 and x0/8>arb(row['presieve']['T']).log()
    slack=alpha*x-beta*x.log()-cost
    assert slack>0 and x0>240
    F=2*eg/s*(1+acb.integral(lambda y,an:(y-2).log(an)/(y-1),acb(3),acb(s)).real)
    f=2*eg*(s-1).log()/s;h=3*(-s).exp()/s
    lower=8*(1-v)*((1+ep)*f-(ep+em)*F
        -eps*arb(2).exp()*h*((1+ep)*108+(ep+em)*106))
    I=eg/(4*t)*(6*(3-8*alpha)/(3-18*alpha)).log()
    FM=2*eg/k+106*eps;J=q('8/3').log();b0=t-q('1/4')
    Hh=acb.integral(lambda y,an:(2-8*(acb(t)-y)).exp()/y,acb(q('1/8')),acb(b0)).real
    Hh+=(q('1/3')/b0).log()
    upper=4*(1+ep)*(1+2/x)*(1+v)/(1-(-x/8).exp())*(I+106*eps*Hh+2*FM/x**2)
    cb=acb.integral(lambda y,an:(2-3*y).log(an)/(y*(1-y)),acb(q('1/8')),acb(q('1/3'))).real
    switch=(1+ep)*(1+v)*(eg/t+159*eps)*(1+delta)*(1+4/x)
    switch*=cb+2/x**2+(J+2/x**2)*(6*ell/x+8/x**2)
    def P(y):
        value=q('165/1000')+q('12683/1000')/y+q('254980/1000')/y**2
        value+=q('2607854/1000')/y**3+q('11605056/1000')/y**4
        value+=(q('1314/1000')*y+q('92/1000')*y.log()+q('60883/1000')
            +q('8250/1000')*y.log()/y+q('939260/1000')/y)*(-y/4).exp()
        return value+1/(16*arb.pi())+2*(-y/2).exp()/arb(2).log()
    count=row['binned_support_count'];Q=row['sharp_support']['Q']
    logH=q(row['presieve']['u'])*arb(row['presieve']['T']).log()+q(row['sharp_support']['kappa']).log()
    m=sum(p<=row['presieve']['exact_cutoff'] and p>2 for p in [2,3,5])
    fraction=arb(2**m)/Q*arb(count['count_upper_ball'])*(-logH).exp()
    assert abs(fraction/arb(count['modulus_count_fraction_ball'])-1)<arb('1e-45')
    medium=q(row['binned_modulus_count_budget']);density=q(row['rough_density_budget'])
    assert 0<fraction<medium and V*(1+ep)<density
    floor0=(cost+arb(2).log()-t*x0).exp()
    floorx=(cost+arb(2).log()-t*x).exp()
    leading=q('3/2')*medium*(density+floor0)*P(x0)
    plateau=leading*((3-beta)*A.log()).exp()/U
    pointwise=q('3/2')*medium*(density+floorx)*P(x)*((3-beta)*x.log()).exp()/U
    assert pointwise<plateau
    other=q('11/200')*x**4*(-x/6).exp()/U
    other+=q('4/5')*(1+ep)*(1+v)*FM*x**2*(-x/6).exp()
    ap=plateau+other
    rectangle=(1+delta)*x**5/(2*U*ell)*(q('3/4')*(-x/48).exp()+(-x/8).exp())
    finite=x**2/U*(2*(-x/8).exp()+(-2*x/3).exp()+(-x).exp())
    values={'lower_ball':lower,'half_upper_ball':upper,'half_switch_ball':switch,
        'AP_budget_ball':ap,'half_rectangle_error_ball':rectangle,'finite_error_ball':finite,
        'support_slack_ball':slack,'margin_ball':lower-upper-switch-ap-rectangle-finite,
        'AP_plateau_leading_ball':plateau,'AP_leading_coefficient_ball':leading,
        'rough_floor_error_ratio_ball':floor0,
        'AP_pointwise_leading_at_start_ball':q('3/2')*medium*(density+floor0)*P(x0)*((3-beta)*x0.log()).exp()/U}
    assert values['margin_ball']>q(row['claimed_margin'])
    if L==row['log_N0']:
        for name,value in values.items():
            assert value>0,name
            assert abs(value/arb(row[name])-1)<arb('1e-43'),name
        assert rectangle<arb(row['rectangle_display_bound'])
        assert finite<arb(row['finite_display_bound'])
    return {'log_N':L,'margin_ball':str(values['margin_ball']),
        'pointwise_AP_leading_ball':str(pointwise),'plateau_leading_ball':str(plateau),
        'support_slack_ball':str(slack),'positive_floor_charge_retained':True,
        'independent_endpoint_quantities':len(values) if L==row['log_N0'] else 0}


def run():
    report=json.loads((ROOT/'conditional_chen_rough_interval_certificate.json').read_text(encoding='utf-8'))
    small=json.loads((ROOT/'conditional_chen_binned_support_checks.json').read_text(encoding='utf-8'))
    tail=json.loads((ROOT/'conditional_chen_sparse_modulus_certificate.json').read_text(encoding='utf-8'))
    for r in [report,small,tail]:
        for name,h in r['input_sha256'].items():
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h,name
    row=report['GRH_positivity'];start,anchor=row['valid_log_interval']
    assert start==row['log_N0']<anchor==tail['GRH_positivity']['log_N0']
    assert Fraction(row['claimed_margin'])==Fraction(tail['GRH_positivity']['claimed_margin'])
    assert report['GRH_at_31000']==tail['GRH_at_31000']
    ps=primes_direct(row['presieve']['T'])
    V=math.prod((1-arb(1)/p for p in ps),start=arb(1))
    assert len(ps)==report['complete_integer_rough_product_prime_count']
    assert abs(V/arb(report['integer_rough_product_ball'])-1)<arb('1e-43')
    count=row['binned_support_count'];assert count['sigma']=='0'
    logH=q(row['presieve']['u'])*arb(row['presieve']['T']).log()+q(row['sharp_support']['kappa']).log()
    width=q(count['log_bin_width'])
    assert int((logH/width).ceil().unique_fmpz())==count['label_limit']
    medium=[p for p in ps if p>row['presieve']['exact_cutoff']]
    K=count['maximum_actual_degree']
    assert arb(math.prod(medium[:K])).log()<logH<arb(math.prod(medium[:K+1])).log()
    assert len(medium)==count['medium_prime_count']
    samples=sorted({start,start+(anchor-start)//4,(start+anchor)//2,start+3*(anchor-start)//4,anchor})
    interval=[direct(row,L,anchor,V) for L in samples]
    old_checks=[grh_direct(tail[k]) for k in ['GRH_positivity','GRH_at_31000']]
    names=[Path(__file__).name,'conditional_chen_rough_interval_certificate.json',
        'conditional_chen_binned_support_checks.json','conditional_chen_sparse_modulus_certificate.json',
        'conditional_chen_prefix_block_endpoint_checks.py']
    return {'status':'PASS: separate scalar formulas, complete integer density and analytic interval join',
        'precision_bits':192,'complete_integer_primes_independently_generated':len(ps),
        'complete_medium_prime_count':len(medium),'count_degree_bound':K,
        'binned_histogram_exact_small_subsets':small['binned_count_checks']['enumerated_subsets'],
        'aggregate_prime_factor_modulus_patterns':small['aggregate_prime_factor_checks']['modulus_patterns'],
        'interval_checks':interval,'old_tail_checks':old_checks,
        'analytic_join':{'covered_log_interval':[start,anchor],'existing_half_line_start':anchor,
            'common_claim':row['claimed_margin'],'inclusive_endpoints':True},
        'continuous_interval_dependencies':'fixed alpha; constant AP plateau; increasing support slack and H; decreasing upper masses, other charges and q_+',
        'input_sha256':{n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names},
        'not_checked':['general histogram and density-domination proofs','external analytic inputs',
            'independent research review','global optimality','finite Goldbach splice','PDF compilation']}


if __name__=='__main__':
    r=run()
    Path(__file__).with_suffix('.json').write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:r[k] for k in ['status','complete_integer_primes_independently_generated','analytic_join','interval_checks']},indent=2))
