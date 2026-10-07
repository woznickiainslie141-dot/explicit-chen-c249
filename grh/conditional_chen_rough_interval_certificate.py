"""GRH finite-interval endpoints using binned medium counts and rough factors.

The certified interval is joined to the existing GRH half-line at log N=26820.
This is an analytic regime join, not a join to finite Goldbach verification.
The general binned-count and integer-roughness lemmas remain dependencies.
"""
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
from flint import arb,ctx
from conditional_chen_rankin_certificate import odd_primes_to
from conditional_chen_binned_support_certificate import count_bound,ball
from conditional_chen_prefix_block_endpoint_certificate import make_case
from conditional_chen_sharp_endpoint_certificate import grh_case,qG_upper,Umin
from conditional_chen_dense_band_grh_certificate import upper_reciprocal

ROOT=Path(__file__).resolve().parent
CACHE=ROOT/'tmp/binned_support_counts'
ctx.prec=192


def digest(name):
    return hashlib.sha256((ROOT/name).read_bytes()).hexdigest()


def binned(pre,sigma=Fraction(0)):
    CACHE.mkdir(parents=True,exist_ok=True)
    u=Fraction(pre['u']);kappa=Fraction(11,169);h=Fraction(1,100)
    key=f"T{pre['T']}_u{u.numerator}-{u.denominator}_s{sigma.numerator}-{sigma.denominator}.json"
    path=CACHE/key
    names=['conditional_chen_binned_support_certificate.py','conditional_chen_rankin_certificate.py',
           'conditional_chen_prefix_block_certificate.py']
    hashes={n:digest(n) for n in names}
    if path.exists():
        r=json.loads(path.read_text(encoding='utf-8'))
        assert r['input_sha256']==hashes
        return r,False
    r=count_bound(pre['T'],5,u,kappa,sigma,h)
    path.write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8',newline='\n')
    return r,True


def decimal_upper(value,places=12):
    scale=10**places;n=math.ceil(float(value)*scale)
    while not value<arb(n)/scale:n+=1
    return f'{n}/{scale}'


def endpoint(case,L,beta,anchor,claim='-10'):
    x=arb(L);B=arb(anchor);b=ball(beta)
    assert 0<b<=3 and L<=anchor
    cost=arb(case['effective_log_support_cost_ball'])
    a=math.ceil(float((cost+b*x.log())/x)*10**6)
    while not arb(a)/10**6*x-b*x.log()>cost:a+=1
    alpha=arb(a)/10**6;t=arb(1)/2-alpha
    raw=grh_case(case,arb(0),case['sharp_support'],L,a,'1','1','-10')
    assert alpha-b/x>0 and arb(1)/2-b/x>0
    assert x/2-b*x.log()>arb(10**9).log() and x/6-b*x.log()>0
    # The rough-integer upper weight has support below 2 Q H0.
    boundary_ratio=(cost+arb(2).log()-t*x).exp()
    rough=ball(case['rough_density_budget'])+boundary_ratio
    medium=ball(case['binned_modulus_count_budget'])
    P=qG_upper(x)+1/(16*arb.pi())+2*(-x/2).exp()/arb(2).log()
    leading=arb(3)/2*medium*rough*P
    AP_plateau=leading*((3-b)*B.log()).exp()/Umin
    AP_at_start=leading*((3-b)*x.log()).exp()/Umin
    assert AP_at_start<=AP_plateau or abs(AP_at_start-AP_plateau)<arb('1e-50')
    k=8*(t-arb(1)/3);FM=2*arb.const_euler().exp()/k+106*ball(case['epsilon'])
    other=arb('.055')*x**4*(-x/6).exp()/Umin
    other+=arb('.8')*(1+ball(case['eta_plus']))*(1+arb(1)/10**10)*FM*x**2*(-x/6).exp()
    AP=AP_plateau+other
    margin=arb(raw['lower_ball'])-arb(raw['half_upper_ball'])-arb(raw['half_switch_ball'])
    margin-=AP+arb(raw['half_rectangle_error_ball'])+arb(raw['finite_error_ball'])
    assert margin>ball(claim)
    raw.update({'beta':str(beta),'claimed_margin':claim,'AP_budget_ball':str(AP),
        'margin_ball':str(margin),
        'AP_plateau_leading_ball':str(AP_plateau),'AP_pointwise_leading_at_start_ball':str(AP_at_start),
        'AP_leading_coefficient_ball':str(leading),'rough_floor_error_ratio_ball':str(boundary_ratio),
        'effective_log_support_cost_ball':str(cost),'support_slack_ball':str(alpha*x-b*x.log()-cost),
        'valid_log_interval':[L,anchor],
        'sieve_interface':'binned medium support and T-rough large factors; constant AP plateau on a finite interval; existing GRH tail'})
    return raw


def first_endpoint(case,beta,anchor):
    lo,hi=8500,anchor
    assert arb(endpoint(case,lo,beta,anchor)['margin_ball'])<arb('.00001')
    if not arb(endpoint(case,hi,beta,anchor)['margin_ball'])>arb('.00001'):return None
    while hi-lo>1:
        mid=(lo+hi)//2
        if arb(endpoint(case,mid,beta,anchor)['margin_ball'])>arb('.00001'):hi=mid
        else:lo=mid
    selected=endpoint(case,hi,beta,anchor,'.00001')
    before=endpoint(case,hi-1,beta,anchor)
    assert arb(before['margin_ball'])<arb('.00001')
    selected['preceding_integer_recipe_margin_ball']=before['margin_ball']
    return selected


def run():
    tail=json.loads((ROOT/'conditional_chen_sparse_modulus_certificate.json').read_text(encoding='utf-8'))
    base=json.loads((ROOT/'conditional_chen_prefix_block_endpoint_certificate.json').read_text(encoding='utf-8'))
    refined=json.loads((ROOT/'conditional_chen_prefix_block_refinement.json').read_text(encoding='utf-8'))
    checks=json.loads((ROOT/'conditional_chen_binned_support_checks.json').read_text(encoding='utf-8'))
    for result in [tail,base,refined,checks]:
        for n,v in result['input_sha256'].items():assert digest(n)==v,n
    assert checks['status'].startswith('PASS')
    anchor=tail['GRH_positivity']['log_N0'];assert anchor==26820
    uniform=json.loads((ROOT/'conditional_chen_uniform_product_certificate.json').read_text(encoding='utf-8'))
    uniform_rows={r['T']:r for r in uniform['rows']}
    primes=[2]+odd_primes_to(5000000)
    V=math.prod((1-arb(1)/p for p in primes),start=arb(1))
    finite=[pre for result in [base,refined] for pre in result['finite_cases']
            if pre['T']==5000000 and pre['exact_cutoff']==5]
    rows=[];cases=[];counts=[];generation=[]
    for pre in finite:
        count,fresh=binned(pre);counts.append(count)
        generation.append({'u':pre['u'],'h':count['log_bin_width'],'fresh_in_this_run':fresh})
        case=make_case(pre,uniform_rows)
        rough=V*(1+ball(case['eta_plus']))
        case.update({'binned_support_count':count,
            'binned_modulus_count_budget':upper_reciprocal(arb(count['modulus_count_fraction_ball'])),
            'integer_rough_product_ball':str(V),'rough_density_ball':str(rough),
            'rough_density_budget':decimal_upper(rough)})
        cases.append(case)
        candidates=[first_endpoint(case,Fraction(n,10),anchor) for n in range(18,31)]
        candidates=[r for r in candidates if r is not None];rows+=candidates
        print(f'Rough interval u={pre["u"]},gap h={pre["log_bin_width"]}: '+
              str(min((r['log_N0'] for r in candidates),default=None)),flush=True)
    selected=min(rows,key=lambda r:(r['log_N0'],-float(arb(r['margin_ball']))))
    # Refine beta for the two best distinct finite gap cases; finite-grid
    # selection is not an optimality claim.
    best_cases=sorted(cases,key=lambda case:min((r['log_N0'] for r in rows
            if r['presieve']==case['presieve']),default=10**9))[:2]
    for case in best_cases:
        matched=[r for r in rows if r['presieve']==case['presieve']]
        coarse=min(matched,key=lambda r:r['log_N0']);center=Fraction(coarse['beta'])
        for n in range(max(180,int(center*100)-9),min(300,int(center*100)+9)+1):
            row=first_endpoint(case,Fraction(n,100),anchor)
            if row is not None:rows.append(row)
    selected=min(rows,key=lambda r:(r['log_N0'],-float(arb(r['margin_ball']))))
    assert selected['log_N0']<anchor and arb(selected['margin_ball'])>arb('.00001')
    selected['rectangle_display_bound']='1e-194'
    selected['finite_display_bound']='1e-1337'
    assert arb(selected['half_rectangle_error_ball'])<arb(selected['rectangle_display_bound'])
    assert arb(selected['finite_error_ball'])<arb(selected['finite_display_bound'])
    names=[Path(__file__).name,'conditional_chen_binned_support_certificate.py',
        'conditional_chen_binned_support_checks.py','conditional_chen_binned_support_checks.json',
        'conditional_chen_sparse_modulus_certificate.json','conditional_chen_prefix_block_endpoint_certificate.json',
        'conditional_chen_prefix_block_refinement.json','conditional_chen_uniform_product_certificate.json',
        'conditional_chen_prefix_block_endpoint_certificate.py','conditional_chen_sharp_endpoint_certificate.py',
        'conditional_chen_dense_band_grh_certificate.py','conditional_chen_rankin_certificate.py']
    return {'status':'PASS: strict binned support counts, rough-integer bound and finite GRH interval',
        'precision_bits':192,'GRH_positivity':selected,'GRH_at_31000':tail['GRH_at_31000'],
        'analytic_tail_join_log_N':anchor,'analytic_tail_claim':tail['GRH_positivity']['claimed_margin'],
        'complete_integer_rough_product_prime_count':len(primes),'integer_rough_product_ball':str(V),
        'candidate_rows':[{'u':r['presieve']['u'],'gap_h':r['presieve']['log_bin_width'],
            'beta':r['beta'],'log_N0':r['log_N0'],'margin_ball':r['margin_ball']} for r in rows],
        'binned_counts':counts,'binned_regeneration':generation,
        'input_sha256':{n:digest(n) for n in names},
        'not_checked':['general binned-count and integer-roughness lemmas','external analytic inputs',
            'independent research review','global optimality','connection to finite Goldbach verification','PDF compilation']}


if __name__=='__main__':
    r=run()
    Path(__file__).with_suffix('.json').write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:r['GRH_positivity'][k] for k in ['log_N0','alpha','beta','valid_log_interval','margin_ball','AP_budget_ball','binned_modulus_count_budget','rough_density_budget']},indent=2))
