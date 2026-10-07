"""Independent formulas and prime sieve for accepted-support GRH endpoints.

The production endpoint evaluator is not imported. Continuous coverage
depends on the stated general support and monotonicity arguments.
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
    assert 0<beta<=3 and x0<=x and x0>=8500
    if Fraction(row['beta'])<3:assert x<=A
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
    if Fraction(row['beta'])==3 and L==row['log_N0']:
        assert abs(pointwise/plateau-1)<arb('1e-45')
    else:assert pointwise<plateau
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
    names=['conditional_chen_accepted_interval_certificate.json',
        'conditional_chen_accepted_support_checks.json','conditional_chen_rough_interval_certificate.json',
        'conditional_chen_rough_interval_checks.json','conditional_chen_binned_support_checks.json',
        'conditional_chen_sparse_modulus_certificate.json']
    reports=[]
    for name in names:
        r=json.loads((ROOT/name).read_text(encoding='utf-8'))
        assert r['status'].startswith('PASS')
        for p,h in r['input_sha256'].items():
            assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
        reports.append(r)
    report,small,previous,previous_checks,binned,tail=reports
    weak,strong=report['GRH_positivity'],report['GRH_at_31000']
    start,anchor=weak['valid_log_interval']
    assert start==weak['log_N0']<anchor==previous['GRH_positivity']['log_N0']
    assert previous['GRH_positivity']['valid_log_interval']==[anchor,26820]
    assert tail['GRH_positivity']['log_N0']==26820
    assert report['analytic_join_chain']==[start,anchor,26820,'infinity']
    assert Fraction(weak['claimed_margin'])==Fraction(previous['GRH_positivity']['claimed_margin'])==Fraction(tail['GRH_positivity']['claimed_margin'])
    assert strong['log_N0']==31000 and strong['valid_log_interval']==[31000,None]
    assert Fraction(strong['beta'])==3 and strong['fixed_beta_three_half_line']
    assert Fraction(strong['claimed_margin'])>Fraction(tail['GRH_at_31000']['claimed_margin'])
    ps=primes_direct(5000000)
    V=math.prod((1-arb(1)/p for p in ps),start=arb(1))
    assert len(ps)==report['complete_integer_rough_product_prime_count']==348513
    assert abs(V/arb(report['integer_rough_product_ball'])-1)<arb('1e-43')
    old_counts={r['u']:r for r in previous['binned_counts']}
    count_checks=[]
    for row in [weak,strong]:
        count=row['binned_support_count']
        assert count['sigma']=='0' and count['log_bin_width']=='1/100'
        assert count['union_without_duplication'] and count['retains_all_preceding_checked_prefixes']
        assert count['complete_cutoffs_retained_after_deletion']
        states=count['state_count_balls']
        assert [(r['alive_flags'],r['selected_parity']) for r in states]==[(1,0),(1,1),(2,0),(2,1),(3,0),(3,1)]
        integers=[arb(r['count_ball']).unique_fmpz() for r in states]
        total=arb(count['count_upper_ball']).unique_fmpz()
        assert None not in integers and total is not None and sum(integers)==total
        assert arb(count['count_upper_ball'])<arb(old_counts[row['presieve']['u']]['count_upper_ball'])
        logH=q(row['presieve']['u'])*arb(5000000).log()+q(row['sharp_support']['kappa']).log()
        width=q(count['log_bin_width'])
        assert int((logH/width).ceil().unique_fmpz())==count['label_limit']
        medium=[p for p in ps if p>5];K=count['maximum_actual_degree']
        assert arb(math.prod(medium[:K])).log()<logH<arb(math.prod(medium[:K+1])).log()
        assert len(medium)==count['medium_prime_count']==348510
        count_checks.append({'u':row['presieve']['u'],'degree':K,'label_limit':count['label_limit'],
                             'union_integer_count':str(total),'six_states_sum_exactly':True})
    samples=sorted({start,start+(anchor-start)//4,(start+anchor)//2,start+3*(anchor-start)//4,anchor})
    interval=[direct(weak,L,anchor,V) for L in samples]
    half_line=[direct(strong,L,31000,V) for L in [31000,32000,62000]]
    previous_interval=[direct(previous['GRH_positivity'],L,26820,V) for L in [anchor,26820]]
    old_tails=[grh_direct(tail[k]) for k in ['GRH_positivity','GRH_at_31000']]
    return {'status':'PASS: separate endpoint formulas, six-state count identities and analytic coverage interfaces',
        'precision_bits':192,'complete_integer_primes_independently_generated':len(ps),
        'complete_medium_prime_count':len(medium),'selected_count_checks':count_checks,
        'accepted_exact_small_subsets':small['enumerated_complete_subsets'],
        'accepted_exact_integer_coefficients':small['integer_histogram_coefficients_checked'],
        'accepted_deleted_subset_patterns':small['deleted_subset_patterns'],
        'accepted_deleted_production_histograms':small['deleted_production_histograms_checked'],
        'aggregate_prime_factor_modulus_patterns':binned['aggregate_prime_factor_checks']['modulus_patterns'],
        'interval_checks':interval,'strong_half_line_checks':half_line,
        'previous_interval_checks':previous_interval,'previous_tail_checks':old_tails,
        'analytic_join':{'covered_log_interval':[start,anchor],
            'prior_interval':[anchor,26820],'existing_half_line_start':26820,
            'common_claim':weak['claimed_margin'],'inclusive_endpoints':True},
        'strong_half_line':{'start':31000,'beta':'3','claim':strong['claimed_margin'],
            'dependencies':'fixed alpha; L^(3-beta)=1; decreasing P, rough floor and every other expense'},
        'continuous_interval_dependencies':'general six-state union domination; fixed alpha; constant AP plateau; increasing support slack and H; decreasing upper masses and other charges; prior GRH interval and tail',
        'input_sha256':{n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest()
            for n in names+[Path(__file__).name,'conditional_chen_prefix_block_endpoint_checks.py']},
        'not_checked':['general six-state and density-domination proofs','external analytic inputs',
            'independent research review','global optimality','finite Goldbach splice','PDF compilation']}


if __name__=='__main__':
    r=run()
    Path(__file__).with_suffix('.json').write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:r[k] for k in ['status','analytic_join','strong_half_line','interval_checks','strong_half_line_checks']},indent=2))
