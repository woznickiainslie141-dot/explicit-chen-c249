"""GRH AP charges based on a proved count of the actual coefficient support.

All finite Rankin Euler products are regenerated; saved finite prefix-block
gaps and ordered-product inputs are checked by hash. The separate general
support-count lemma is required. No new distribution hypothesis is used.
"""
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
from flint import arb,ctx
from conditional_chen_rankin_certificate import odd_primes_to
from conditional_chen_prefix_block_endpoint_certificate import make_case
from conditional_chen_sharp_endpoint_certificate import grh_case,qG_upper,Umin
from conditional_chen_dense_band_grh_certificate import upper_reciprocal

ROOT=Path(__file__).resolve().parent
ctx.prec=192


def q(value):
    value=Fraction(value)
    return arb(value.numerator)/value.denominator


def count_products(T,w):
    primes=[p for p in odd_primes_to(T) if p>w]
    sigmas=[Fraction(3,4),Fraction(4,5),Fraction(17,20),Fraction(9,10)]
    rows=[]
    for sigma in sigmas:
        s=q(sigma)
        value=sum(((1+(-s*arb(p).log()).exp()).log() for p in primes),arb(0))
        rows.append({'sigma':str(sigma),'log_euler_product_ball':str(value)})
    print(f'Complete Rankin count products T={T},w={w}, primes={len(primes)}',flush=True)
    return {'T':T,'exact_cutoff':w,'medium_prime_count':len(primes),'rows':rows}


def count_case(case,products):
    pre=case['presieve']; support=case['sharp_support']
    exact=[p for p in odd_primes_to(pre['exact_cutoff'])]
    assert math.prod(exact)==support['Q']
    logH=q(pre['u'])*arb(pre['T']).log()+q(support['kappa']).log()
    rows=[]
    for item in products['rows']:
        delta=arb(2**len(exact))/support['Q']*((q(item['sigma'])-1)*logH
                  +arb(item['log_euler_product_ball'])).exp()
        rows.append({**item,'count_fraction_ball':str(delta)})
    selected=min(rows,key=lambda r:float(arb(r['count_fraction_ball'])))
    assert 0<arb(selected['count_fraction_ball'])<arb('.608')
    budget=upper_reciprocal(arb(selected['count_fraction_ball']))
    return {**case,'sparse_modulus_count_budget':budget,
            'sparse_modulus_count':{'medium_support_log_ball':str(logH),
                'exact_prime_count':len(exact),'exact_divisor_count_max':2**len(exact),
                'selected':selected,'candidate_rows':rows,
                'bound':'# coefficients <= 2^omega(Q) D H0^sigma product_(w<p<=T)(1+p^(-sigma))'}}


def endpoint(case,L,beta=Fraction(3),claim='-10'):
    x=arb(L); b=q(beta); cost=arb(case['effective_log_support_cost_ball'])
    a=math.ceil(float((cost+b*x.log())/x)*10**6)
    while not arb(a)/10**6*x-b*x.log()>cost:a+=1
    alpha=arb(a)/10**6
    # Only the core sieve quantities and analytic ranges are reused. Passing
    # a zero auxiliary support cost avoids its fixed-beta legacy check;
    # the actual support and AP level are checked explicitly below.
    row=grh_case(case,arb(0),case['sharp_support'],L,a,'1','1','-10')
    assert b>=3 and alpha-b/x>0 and arb(1)/2-b/x>0
    assert x/2-b*x.log()>arb(10**9).log()
    delta=q(case['sparse_modulus_count_budget'])
    k=8*(arb(1)/6-alpha); FM=2*arb.const_euler().exp()/k+106*q(case['epsilon'])
    leading=arb(3)/2*delta*(qG_upper(x)+1/(16*arb.pi())+2*(-x/2).exp()/arb(2).log())
    ap=leading*((3-b)*x.log()).exp()/Umin
    ap+=arb('.055')*x**4*(-x/6).exp()/Umin
    ap+=arb('.8')*(1+q(case['eta_plus']))*(1+arb(1)/10**10)*FM*x**2*(-x/6).exp()
    margin=arb(row['lower_ball'])-arb(row['half_upper_ball'])-arb(row['half_switch_ball'])
    margin-=ap+arb(row['half_rectangle_error_ball'])+arb(row['finite_error_ball'])
    assert margin>q(claim)
    row.update({'beta':str(beta),'claimed_margin':claim,'AP_budget_ball':str(ap),
                'AP_leading_coefficient_ball':str(leading),'margin_ball':str(margin),
                'effective_log_support_cost_ball':str(cost),
                'support_slack_ball':str(alpha*x-b*x.log()-cost),
                'sieve_interface':'actual coefficient count from complete medium-prime Rankin product; unchanged pointwise GRH AP input'})
    return row


def first_endpoint(case):
    lo,hi=8500,40000; claim=arb('.00001')
    assert arb(endpoint(case,lo)['margin_ball'])<claim
    if not arb(endpoint(case,hi)['margin_ball'])>claim:return None
    while hi-lo>1:
        mid=(hi+lo)//2
        if arb(endpoint(case,mid)['margin_ball'])>claim:hi=mid
        else:lo=mid
    row=endpoint(case,hi,claim='.00001')
    before=endpoint(case,hi-1)
    assert arb(before['margin_ball'])<claim
    row['preceding_integer_recipe_margin_ball']=before['margin_ball']
    row['rectangle_display_bound']='1e-60'; row['finite_display_bound']='1e-450'
    assert arb(row['half_rectangle_error_ball'])<arb('1e-60')
    assert arb(row['finite_error_ball'])<arb('1e-450')
    return row


def exact_small_count_checks():
    count=0
    for w in [3,5]:
        ps=[p for p in odd_primes_to(29) if p>w]
        for H in [100,1000,10000]:
            exact=sum(math.prod((p for i,p in enumerate(ps) if mask>>i&1))<H
                      for mask in range(1<<len(ps)))
            for sigma in [Fraction(1,2),Fraction(4,5),Fraction(1)]:
                bound=q(H)**q(sigma)*math.prod((1+(-q(sigma)*arb(p).log()).exp() for p in ps),start=arb(1))
                assert arb(exact)<bound
                count+=1
    return count


def run():
    base=json.loads((ROOT/'conditional_chen_prefix_block_endpoint_certificate.json').read_text(encoding='utf-8'))
    refined=json.loads((ROOT/'conditional_chen_prefix_block_refinement.json').read_text(encoding='utf-8'))
    for result in [base,refined]:
        for p,v in result['input_sha256'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==v,p
    uniform=json.loads((ROOT/'conditional_chen_uniform_product_certificate.json').read_text(encoding='utf-8'))
    uniform_rows={r['T']:r for r in uniform['rows']}
    finite=[p for result in [base,refined] for p in result['finite_cases']
            if p['exact_cutoff']==5 and p['T']>=5000000]
    pairs=sorted({(p['T'],p['exact_cutoff']) for p in finite})
    products={(T,w):count_products(T,w) for T,w in pairs}
    rows=[]
    for pre in finite:
        case=count_case(make_case(pre,uniform_rows),products[(pre['T'],pre['exact_cutoff'])])
        row=first_endpoint(case)
        if row is not None:rows.append(row)
    selected=min(rows,key=lambda r:(r['log_N0'],-float(arb(r['margin_ball']))))
    case={k:selected[k] for k in ['presieve','epsilon','eta_plus','eta_minus','C1','C2',
                'effective_log_support_cost_ball','sharp_support','sparse_modulus_count_budget','sparse_modulus_count']}
    strong=endpoint(case,31000)
    n=int((arb(strong['margin_ball'])*10000).floor().unique_fmpz())
    strong['claimed_margin']=f'{n}/10000'
    assert arb(strong['margin_ball'])>q(strong['claimed_margin'])
    names=[Path(__file__).name,'conditional_chen_prefix_block_refinement.json',
            'conditional_chen_prefix_block_endpoint_certificate.json',
            'conditional_chen_prefix_block_endpoint_certificate.py',
            'conditional_chen_sharp_endpoint_certificate.py',
            'conditional_chen_dense_band_grh_certificate.py',
            'conditional_chen_rankin_certificate.py','conditional_chen_uniform_product_certificate.json']
    return {'status':'PASS: strict regenerated support-count products and GRH scalar endpoints',
            'precision_bits':192,'GRH_positivity':selected,'GRH_at_31000':strong,
            'candidate_rows':[{'T':r['presieve']['T'],'u':r['presieve']['u'],'h':r['presieve']['log_bin_width'],
                'log_N0':r['log_N0'],'count_budget':r['sparse_modulus_count_budget'],'margin_ball':r['margin_ball']} for r in rows],
            'complete_Rankin_products':list(products.values()),
            'exact_small_Rankin_count_checks':exact_small_count_checks(),
            'input_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in names},
            'not_checked':['general support-count lemma','external analytic inputs',
                           'independent research review','global optimality','full splice','PDF compilation']}


if __name__=='__main__':
    result=run()
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:{n:result[k][n] for n in ['log_N0','alpha','beta','margin_ball','claimed_margin','sparse_modulus_count_budget']}
                      for k in ['GRH_positivity','GRH_at_31000']},indent=2))
