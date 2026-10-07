"""GRH endpoints using the positive union count of actual medium supports.

The weak finite interval joins the existing half-line at log N=24799.
The stronger count is certified separately with beta=3 on a half-line.
"""
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from flint import arb,ctx
from conditional_chen_accepted_support_certificate import certify
from conditional_chen_prefix_block_endpoint_certificate import make_case
from conditional_chen_rough_interval_certificate import endpoint,first_endpoint,decimal_upper
from conditional_chen_dense_band_grh_certificate import upper_reciprocal
from conditional_chen_binned_support_certificate import ball

ROOT=Path(__file__).resolve().parent
CACHE=ROOT/'tmp/accepted_support_counts'
ctx.prec=192


def digest(name):return hashlib.sha256((ROOT/name).read_bytes()).hexdigest()


def checked(name):
    r=json.loads((ROOT/name).read_text(encoding='utf-8'))
    for p,h in r['input_sha256'].items():assert digest(p)==h,p
    return r


def accepted(pre):
    CACHE.mkdir(parents=True,exist_ok=True)
    u=Fraction(pre['u']);key=f"T{pre['T']}_u{u.numerator}-{u.denominator}_h1-100.json"
    p=CACHE/key
    names=['conditional_chen_accepted_support_certificate.py','conditional_chen_rankin_certificate.py',
        'conditional_chen_prefix_block_certificate.py','conditional_chen_sharp_support_certificate.py',
        'conditional_chen_binned_support_certificate.py']
    hashes={n:digest(n) for n in names}
    if p.exists():
        r=json.loads(p.read_text(encoding='utf-8'))
        if r['input_sha256']==hashes:return r,False,'per-case hash-checked cache'
    standalone=ROOT/'conditional_chen_accepted_support_certificate.json'
    if standalone.exists():
        r=json.loads(standalone.read_text(encoding='utf-8'))
        if (r['input_sha256']==hashes and r['T']==pre['T'] and r['exact_cutoff']==5
                and Fraction(r['u'])==u and r['log_bin_width']=='1/100'):
            p.write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8',newline='\n')
            return r,False,'standalone complete count, hash checked'
    r=certify(pre['T'],5,u,progress=True)
    p.write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8',newline='\n')
    return r,True,'fresh complete count'


def run():
    previous=checked('conditional_chen_rough_interval_certificate.json')
    previous_checks=checked('conditional_chen_rough_interval_checks.json')
    small=checked('conditional_chen_accepted_support_checks.json')
    base=checked('conditional_chen_prefix_block_endpoint_certificate.json')
    refined=checked('conditional_chen_prefix_block_refinement.json')
    assert small['status'].startswith('PASS') and previous_checks['status'].startswith('PASS')
    anchor=previous['GRH_positivity']['log_N0'];assert anchor==24799
    # The earlier rough-interval certificate already anchors this whole file.
    # Its older schema has no input_sha256 mapping of its own.
    uniform=json.loads((ROOT/'conditional_chen_uniform_product_certificate.json').read_text(encoding='utf-8'))
    uniform_rows={r['T']:r for r in uniform['rows']}
    old_counts={r['u']:r for r in previous['binned_counts']}
    V=arb(previous['integer_rough_product_ball'])
    finite=[p for report in [base,refined] for p in report['finite_cases']
            if p['T']==5000000 and p['exact_cutoff']==5]
    cases=[];rows=[];counts=[];generation=[]
    for pre in finite:
        count,fresh,source=accepted(pre)
        assert arb(count['count_upper_ball'])<arb(old_counts[pre['u']]['count_upper_ball'])
        case=make_case(pre,uniform_rows);density=V*(1+ball(case['eta_plus']))
        case.update({'binned_support_count':count,
            'binned_modulus_count_budget':upper_reciprocal(arb(count['modulus_count_fraction_ball'])),
            'integer_rough_product_ball':str(V),'rough_density_ball':str(density),
            'rough_density_budget':decimal_upper(density)})
        counts.append(count);cases.append(case)
        generation.append({'u':pre['u'],'fresh_in_this_run':fresh,'source':source})
        candidates=[first_endpoint(case,Fraction(n,10),anchor) for n in range(15,25)]
        candidates=[r for r in candidates if r is not None];rows+=candidates
        print(f'Accepted interval u={pre["u"]},gap h={pre["log_bin_width"]}: '+
              str(min((r['log_N0'] for r in candidates),default=None)),flush=True)
    best_cases=sorted(cases,key=lambda c:min((r['log_N0'] for r in rows if r['presieve']==c['presieve']),default=10**9))[:2]
    for case in best_cases:
        coarse=min((r for r in rows if r['presieve']==case['presieve']),key=lambda r:r['log_N0'])
        center=Fraction(coarse['beta'])
        for n in range(max(150,int(center*100)-9),min(240,int(center*100)+9)+1):
            r=first_endpoint(case,Fraction(n,100),anchor)
            if r is not None:rows.append(r)
    weak=min(rows,key=lambda r:(r['log_N0'],-float(arb(r['margin_ball']))))
    weak['rectangle_display_bound']='1e-188';weak['finite_display_bound']='1e-1300'
    assert arb(weak['half_rectangle_error_ball'])<arb(weak['rectangle_display_bound'])
    assert arb(weak['finite_error_ball'])<arb(weak['finite_display_bound'])
    strong_rows=[endpoint(case,31000,Fraction(3),31000) for case in cases]
    strong=max(strong_rows,key=lambda r:float(arb(r['margin_ball'])))
    claimed=int((arb(strong['margin_ball'])*10000).floor().unique_fmpz())
    strong['claimed_margin']=f'{claimed}/10000';assert arb(strong['margin_ball'])>ball(strong['claimed_margin'])
    strong['valid_log_interval']=[31000,None]
    strong['fixed_beta_three_half_line']=True
    strong['rectangle_display_bound']='1e-250';strong['finite_display_bound']='1e-1673'
    assert arb(strong['half_rectangle_error_ball'])<arb(strong['rectangle_display_bound'])
    assert arb(strong['finite_error_ball'])<arb(strong['finite_display_bound'])
    for r in [weak,strong]:r['sieve_interface']='six-state union of medium Rosser supports and T-rough aggregate large factors'
    names=[Path(__file__).name,'conditional_chen_accepted_support_certificate.py',
        'conditional_chen_accepted_support_checks.py','conditional_chen_accepted_support_checks.json',
        'conditional_chen_rough_interval_certificate.py','conditional_chen_rough_interval_certificate.json',
        'conditional_chen_rough_interval_checks.json','conditional_chen_sparse_modulus_certificate.json',
        'conditional_chen_prefix_block_endpoint_certificate.json','conditional_chen_prefix_block_refinement.json',
        'conditional_chen_uniform_product_certificate.json','conditional_chen_dense_band_grh_certificate.py',
        'conditional_chen_prefix_block_endpoint_certificate.py','conditional_chen_binned_support_certificate.py']
    return {'status':'PASS: strict accepted-support scalar endpoints; separate general domination proof required',
        'precision_bits':192,'GRH_positivity':weak,'GRH_at_31000':strong,
        'analytic_tail_join_log_N':anchor,'analytic_join_chain':[weak['log_N0'],anchor,26820,'infinity'],
        'integer_rough_product_ball':str(V),'complete_integer_rough_product_prime_count':348513,
        'candidate_rows':[{'u':r['presieve']['u'],'gap_h':r['presieve']['log_bin_width'],
            'beta':r['beta'],'log_N0':r['log_N0'],'margin_ball':r['margin_ball']} for r in rows],
        'stronger_candidates':[{'u':r['presieve']['u'],'gap_h':r['presieve']['log_bin_width'],
            'margin_ball':r['margin_ball']} for r in strong_rows],
        'accepted_support_counts':counts,'accepted_regeneration':generation,
        'input_sha256':{n:digest(n) for n in names},
        'not_checked':['general six-state domination proof','external analytic inputs','independent research review',
            'global optimality','finite Goldbach splice','PDF compilation']}


if __name__=='__main__':
    r=run()
    Path(__file__).with_suffix('.json').write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:{n:r[k][n] for n in ['log_N0','alpha','beta','margin_ball','claimed_margin','binned_modulus_count_budget']}
         for k in ['GRH_positivity','GRH_at_31000']},indent=2))
