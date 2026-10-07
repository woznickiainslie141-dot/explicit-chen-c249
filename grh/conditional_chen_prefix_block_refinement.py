"""Bounded strict refinement of the positive prefix-block endpoint grid.

Finite sums are recomputed, with per-case caches tied to generator hashes.
Saved initial-grid and uniform-product inputs are reused by hash. No numeric
distribution constants or verification-range splice are supplied here.
"""
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from flint import arb, ctx
from conditional_chen_prefix_block_certificate import certify_blocks
from conditional_chen_prefix_block_endpoint_certificate import (
    make_case, grh_endpoint, direct_endpoint)

ROOT=Path(__file__).resolve().parent
CACHE=ROOT/'tmp/prefix_block_refinement'
ctx.prec=192


def digest(name):
    return hashlib.sha256((ROOT/name).read_bytes()).hexdigest()


def finite_case(T,w,u,h):
    CACHE.mkdir(parents=True,exist_ok=True)
    key=f'T{T}_w{w}_u{u.numerator}-{u.denominator}_h{h.denominator}.json'
    path=CACHE/key
    hashes={p:digest(p) for p in ['conditional_chen_prefix_block_certificate.py',
                                'conditional_chen_rankin_certificate.py']}
    if path.exists():
        pre=json.loads(path.read_text(encoding='utf-8'))
        assert pre['input_sha256']==hashes
        assert (pre['T'],pre['exact_cutoff'],pre['u'],pre['log_bin_width'])==(T,w,str(u),str(h))
        return pre,False
    pre=certify_blocks(T,w,u,h,False)
    path.write_text(json.dumps(pre,indent=2)+'\n',encoding='utf-8',newline='\n')
    return pre,True


def first_grh(case,beta,claim='.00001'):
    lo,hi=8500,40000
    assert arb(grh_endpoint(case,lo,beta)['margin_ball'])<arb(claim)
    if not arb(grh_endpoint(case,hi,beta)['margin_ball'])>arb(claim):return None
    while hi-lo>1:
        mid=(hi+lo)//2
        if arb(grh_endpoint(case,mid,beta)['margin_ball'])>arb(claim):hi=mid
        else:lo=mid
    row=grh_endpoint(case,hi,beta,claim)
    before=grh_endpoint(case,hi-1,beta)
    assert arb(before['margin_ball'])<arb(claim)
    row['preceding_integer_recipe_margin_ball']=before['margin_ball']
    row['rectangle_display_bound']='1e-60'
    row['finite_display_bound']='1e-450'
    assert arb(row['half_rectangle_error_ball'])<arb(row['rectangle_display_bound'])
    assert arb(row['finite_error_ball'])<arb(row['finite_display_bound'])
    return row


def run():
    base=json.loads((ROOT/'conditional_chen_prefix_block_endpoint_certificate.json').read_text(encoding='utf-8'))
    for p,v in base['input_sha256'].items():assert digest(p)==v,p
    uniform=json.loads((ROOT/'conditional_chen_uniform_product_certificate.json').read_text(encoding='utf-8'))
    uniform_rows={r['T']:r for r in uniform['rows']}
    finite=[]; grh_rows=[base['GRH_positivity']]; direct_rows=list(base['distribution_endpoints'])
    regeneration=[]
    grid=[(5000000,5,Fraction(n,16),Fraction(1,200)) for n in [87,88,89]]
    grid += [(T,5,Fraction(n,8),Fraction(1,100))
             for T in [10000000,20000000] for n in [43,44,45]]
    betas=[Fraction(n,100) for n in range(350,358)]
    for T,w,u,h in grid:
        pre,fresh=finite_case(T,w,u,h); finite.append(pre)
        regeneration.append({'T':T,'w':w,'u':str(u),'h':str(h),'fresh_in_this_run':fresh})
        case=make_case(pre,uniform_rows)
        rows=[first_grh(case,b) for b in betas]
        grh_rows += [r for r in rows if r is not None]
        print(f'GRH refinement T={T},u={u},h={h}: '+str(min(r['log_N0'] for r in rows if r is not None)),flush=True)
    direct_grid={(r['case']['presieve']['T'],3,Fraction(r['case']['presieve']['u']),Fraction(1,200))
                 for r in base['distribution_endpoints']}
    direct_grid |= {(10000,5,Fraction(n,8),Fraction(1,100)) for n in range(26,33)}
    direct_grid |= {(T,3,Fraction(n,4),Fraction(1,100))
                    for T in [100000,200000,500000] for n in [17,18,19]}
    for T,w,u,h in sorted(direct_grid):
        pre,fresh=finite_case(T,w,u,h);finite.append(pre)
        regeneration.append({'T':T,'w':w,'u':str(u),'h':str(h),'fresh_in_this_run':fresh})
        case=make_case(pre,uniform_rows)
        for theta in [Fraction(3,4),Fraction(9,10),Fraction(19,20),Fraction(99,100)]:
            for C in [1,1000,10000,1000000]:
                row=direct_endpoint(case,theta,C)
                if row is not None:direct_rows.append(row)
        print(f'Direct refinement T={T},w={w},u={u},h={h}',flush=True)
    chosen=min(grh_rows,key=lambda r:(r['log_N0'],-float(arb(r['margin_ball']))))
    selected=[]
    for theta in ['3/4','9/10','19/20','99/100']:
        for C in [1,1000,10000,1000000]:
            rows=[r for r in direct_rows if r['theta']==theta and r['C_upper_budget']==C]
            selected.append(min(rows,key=lambda r:(r['log_start'],-float(arb(r['margin_ball'])))))
    # The stronger count is retained at a fixed, independently reproducible N.
    case={k:chosen[k] for k in ['presieve','epsilon','eta_plus','eta_minus','C1','C2',
                                'effective_log_support_cost_ball','sharp_support']}
    stronger=max((grh_endpoint(case,31000,b) for b in betas),key=lambda r:float(arb(r['margin_ball'])))
    claim=arb(stronger['margin_ball']); unit=10000
    n=int((claim*unit).floor().unique_fmpz())
    assert n>10 and arb(n)/unit<claim
    stronger['claimed_margin']=f'{n}/{unit}'
    inputs=['conditional_chen_prefix_block_endpoint_certificate.json',
            'conditional_chen_prefix_block_endpoint_certificate.py',
            'conditional_chen_prefix_block_certificate.py',
            'conditional_chen_rankin_certificate.py',
            'conditional_chen_uniform_product_certificate.json',Path(__file__).name]
    return {'status':'PASS: bounded strict finite prefix-block and endpoint refinement',
            'precision_bits':ctx.prec,'GRH_positivity':chosen,'GRH_at_31000':stronger,
            'GRH_candidate_rows':[{'T':r['presieve']['T'],'u':r['presieve']['u'],
                'h':r['presieve']['log_bin_width'],'beta':r['beta'],'log_N0':r['log_N0'],
                'margin_ball':r['margin_ball']} for r in grh_rows],
            'distribution_endpoints':selected,'finite_cases':finite,'finite_regeneration':regeneration,
            'uniform_product_inputs_regenerated':False,
            'distribution_C_X':'unspecified; each signed-input family must be specified separately',
            'input_sha256':{p:digest(p) for p in inputs},
            'not_checked':['general analytic proof','external inputs','independent research review',
                           'global optimality','verification-range splice','PDF compilation']}


if __name__=='__main__':
    result=run()
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:{n:result[k][n] for n in ['log_N0','alpha','beta','claimed_margin','margin_ball']}
                      for k in ['GRH_positivity','GRH_at_31000']},indent=2))
    print(json.dumps([{k:r[k] for k in ['theta','C_upper_budget','log_start']}
                      for r in result['distribution_endpoints']],indent=2))
