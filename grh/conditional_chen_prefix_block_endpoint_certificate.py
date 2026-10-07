"""Strict GRH and finite-C endpoint searches with positive prefix blocks.

All selected finite sums are generated here, while the saved uniform
product certificate is reused with a recorded hash. Distribution C,X
remain unspecified inputs; no numerical EH-only record is asserted.
"""
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

from flint import arb,ctx

from conditional_chen_prefix_block_certificate import certify_blocks
from conditional_chen_rankin_certificate import odd_primes_to
from conditional_chen_sharp_support_certificate import support_ratio
from conditional_chen_sharp_endpoint_certificate import grh_case,qG_upper,Umin
from conditional_chen_dense_band_grh_certificate import upper_reciprocal
from conditional_chen_weighted_distribution_certificate import scalar


ROOT=Path(__file__).resolve().parent
ctx.prec=192
TABLE=[(100000,106,107),(10000,106,108),(2000,109,110),(1000,113,114),
       (900,114,115),(800,115,116),(700,116,117),(600,118,119),
       (500,121,122),(400,125,126),(300,133,134),(200,153,154)]


def make_case(pre,uniform_rows):
    pre={**pre,'uniform_epsilon_ball':uniform_rows[pre['T']]['ordered_uniform_epsilon_ball']}
    epsilon=upper_reciprocal(arb(pre['uniform_epsilon_ball']))
    ep=upper_reciprocal(arb(pre['upper_relative_Rosser_gap_ball']))
    em=upper_reciprocal(arb(pre['lower_relative_Rosser_gap_ball']))
    inverse=int(epsilon.split('/')[1])
    C1,C2=next((a,b) for minimum,a,b in TABLE if inverse>=minimum)
    ps=odd_primes_to(pre['T']);medium=[p for p in ps if p>pre['exact_cutoff']]
    kappa,witness=support_ratio(medium)
    assert kappa==(Fraction(11,169) if pre['exact_cutoff']==5 else Fraction(5,49))
    cost=arb(pre['log_support_cost_ball'])+arb(kappa.numerator).log()-arb(kappa.denominator).log()
    support={'T':pre['T'],'exact_cutoff':pre['exact_cutoff'],'u':pre['u'],
             'Q':pre['exact_presieve_Q_max'],'kappa':str(kappa),'maximizer':witness,
             'complete_medium_prime_count':len(medium),
             'original_log_support_cost_ball':pre['log_support_cost_ball'],
             'sharp_log_support_cost_ball':str(cost)}
    return {'presieve':pre,'epsilon':epsilon,'eta_plus':ep,'eta_minus':em,
            'C1':C1,'C2':C2,'effective_log_support_cost_ball':str(cost),'sharp_support':support}


def grh_endpoint(case,L,beta,claim='-10'):
    x=arb(L);b=arb(beta.numerator)/beta.denominator
    cost=arb(case['effective_log_support_cost_ball'])
    assert arb(7)/2<=b<4
    alpha=math.ceil(float((cost+b*x.log())/x)*10**6)
    while not arb(alpha)/10**6*x-b*x.log()>cost:alpha+=1
    r=grh_case(case,cost,case['sharp_support'],L,alpha,'1','1','-10')
    assert x/2-b*x.log()>arb(10**9).log()
    assert arb(alpha)/10**6-b/x>0 and arb(1)/2-b/x>0
    k=8*(arb(1)/6-arb(alpha)/10**6)
    FM=2*arb.const_euler().exp()/k+106*arb(case['epsilon'])
    leading=arb('.912')*(qG_upper(x)+1/(16*arb.pi())+2*(-x/2).exp()/arb(2).log())
    new_ap=leading*x**(3-b)/Umin+arb('.055')*x**4*(-x/6).exp()/Umin
    new_ap+=arb('.8')*(1+arb(case['eta_plus']))*(1+arb(1)/10**10)*FM*x**2*(-x/6).exp()
    margin=arb(r['lower_ball'])-arb(r['half_upper_ball'])-arb(r['half_switch_ball'])
    margin-=new_ap+arb(r['half_rectangle_error_ball'])+arb(r['finite_error_ball'])
    assert margin>arb(claim) and new_ap>0
    r.update({'beta':str(beta),'AP_leading_coefficient_ball':str(leading),
              'AP_budget_ball':str(new_ap),'margin_ball':str(margin),
              'claimed_margin':claim,'support_slack_ball':str(arb(alpha)/10**6*x-b*x.log()-cost),
              'sieve_interface':'positive grouped bounds retain every preceding prefix check; pure Rosser coefficients and sharp support lemma'})
    return r


def first_grh(case,beta,claim='.00001'):
    lo,hi=20000,30747
    assert arb(grh_endpoint(case,lo,beta)['margin_ball'])<arb(claim)
    if not arb(grh_endpoint(case,hi,beta)['margin_ball'])>arb(claim):return None
    while hi-lo>1:
        mid=(hi+lo)//2
        if arb(grh_endpoint(case,mid,beta)['margin_ball'])>arb(claim):hi=mid
        else:lo=mid
    r=grh_endpoint(case,hi,beta,claim)
    before=grh_endpoint(case,hi-1,beta)
    assert arb(before['margin_ball'])<arb(claim)
    r['preceding_integer_recipe_margin_ball']=before['margin_ball']
    return r


def direct_endpoint(case,theta,C):
    lo,hi=83,3000
    if not scalar(case,theta,hi,C)['margin']>arb('.00001'):return None
    while hi-lo>1:
        mid=(hi+lo)//2
        r=scalar(case,theta,mid,C)
        if r is not None and r['margin']>arb('.00001'):hi=mid
        else:lo=mid
    r=scalar(case,theta,hi,C);before=scalar(case,theta,hi-1,C)
    assert r['margin']>arb('.00001')
    assert before is None or before['margin']<arb('.00001')
    return {'theta':str(theta),'C_upper_budget':C,'log_start':hi,
            'threshold_formula':f'max(X, exp({hi}))','case':case,
            **{key+'_ball':str(value) for key,value in r.items()},
            'preceding_margin_ball':str(before['margin']) if before is not None else None}


def run():
    uniform=json.loads((ROOT/'conditional_chen_uniform_product_certificate.json').read_text(encoding='utf-8'))
    uniform_rows={r['T']:r for r in uniform['rows']}
    checks=json.loads((ROOT/'conditional_chen_prefix_block_checks.json').read_text(encoding='utf-8'))
    assert checks['status'].startswith('PASS')
    for p,v in checks['input_sha256'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==v,p
    grh_rows=[];finite=[];direct_rows=[]
    beta_grid=[Fraction(n,100) for n in range(350,358)]
    cached=json.loads((ROOT/'conditional_chen_prefix_block_grh_finite.json').read_text(encoding='utf-8'))
    for n in [41,42,43,44,45,46]:
        u=Fraction(n,8)
        if u==Fraction(11,2):
            pre=cached
            for p,v in pre['input_sha256'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==v,p
        else:pre=certify_blocks(5000000,5,u,Fraction(1,100),False)
        finite.append(pre);case=make_case(pre,uniform_rows)
        for beta in beta_grid:
            r=first_grh(case,beta)
            if r is not None:grh_rows.append(r)
        print(f'Strict prefix-block GRH candidate u={u}',flush=True)
    selected_grh=min(grh_rows,key=lambda r:(r['log_N0'],-float(arb(r['margin_ball']))))
    for T in [10000,22000,70000]:
        for n in range(26,39):
            pre=certify_blocks(T,3,Fraction(n,8),Fraction(1,100),False)
            finite.append(pre);case=make_case(pre,uniform_rows)
            for theta in [Fraction(3,4),Fraction(9,10),Fraction(19,20),Fraction(99,100)]:
                for C in [1,1000,10000,1000000]:
                    r=direct_endpoint(case,theta,C)
                    if r is not None:direct_rows.append(r)
            print(f'Strict prefix-block direct candidate T={T}, u={Fraction(n,8)}',flush=True)
    selected_direct=[]
    for theta in ['3/4','9/10','19/20','99/100']:
        for C in [1,1000,10000,1000000]:
            rows=[r for r in direct_rows if r['theta']==theta and r['C_upper_budget']==C]
            selected_direct.append(min(rows,key=lambda r:(r['log_start'],-float(arb(r['margin_ball'])))))
    inputs=['conditional_chen_prefix_block_certificate.py','conditional_chen_prefix_block_checks.py',
            'conditional_chen_prefix_block_checks.json','conditional_chen_prefix_block_grh_finite.json',
            'conditional_chen_uniform_product_certificate.json','conditional_chen_sharp_support_certificate.py',
            'conditional_chen_sharp_endpoint_certificate.py','conditional_chen_dense_band_grh_certificate.py',
            'conditional_chen_weighted_distribution_certificate.py',Path(__file__).name]
    return {'status':'PASS: strict prefix-block finite sums and scalar endpoint searches; general lemma required',
            'precision_bits':ctx.prec,'GRH_positivity':selected_grh,
            'GRH_candidate_rows':[{'u':r['presieve']['u'],'beta':r['beta'],'log_N0':r['log_N0'],
                                   'margin_ball':r['margin_ball']} for r in grh_rows],
            'distribution_endpoints':selected_direct,'finite_cases':finite,
            'distribution_C_X':'unspecified finite quantitative inputs; no numerical EH-only record',
            'uniform_product_inputs_regenerated':False,
            'input_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs},
            'not_checked':['general domination proof','external analytic inputs','independent mathematical review',
                           'global optimality','full splice','PDF compilation']}


if __name__=='__main__':
    r=run()
    Path(__file__).with_suffix('.json').write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8',newline='\n')
    g=r['GRH_positivity']
    print(json.dumps({'GRH':{p:g[p] for p in ['log_N0','alpha','beta','margin_ball','support_slack_ball']}},indent=2))
    print(json.dumps([{p:x[p] for p in ['theta','C_upper_budget','log_start','margin_ball']} for x in r['distribution_endpoints']],indent=2))
