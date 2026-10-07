"""Strict scalar reconnaissance for the seven-eighths Chen extension.

This is NOT a standalone certification of a Chen theorem.  Finite presieve
inputs are hash-checked from the existing manuscript; the new complex-analysis,
AP, rectangle and half-line arguments are recorded separately for review.
No GRH AP or bilinear estimate is used by this endpoint evaluator.
"""
from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path
from flint import arb, acb, ctx

ROOT = Path(__file__).resolve().parent
ctx.prec = 192
G = arb.const_euler().exp()
UMIN = arb('1.32') / G


def ball(q):
    if isinstance(q, Fraction):
        return arb(q.numerator) / q.denominator
    return arb(q)


def digest(name):
    return hashlib.sha256((ROOT / name).read_bytes()).hexdigest()


def selected_case():
    """Declared rational inputs; the reproduction proves their finite budgets.

    No old GRH certificate, saved prime mass or cached result is read here.
    """
    T, u = 5000000, Fraction(11,2)
    cost = (arb(15)*11/169).log()+ball(u)*arb(T).log()
    return {'presieve':{'T':T,'exact_cutoff':5,'u':str(u),
                       'log_bin_width':'1/200'},
            'epsilon':'1/41219','eta_plus':'1/3231','eta_minus':'1/3242',
            'C1':106,'C2':108,'effective_log_support_cost_ball':str(cost)}


def inputs():
    """Legacy two-case exploratory reader; not used by the new reproduction."""
    r = json.loads((ROOT / 'conditional_chen_accepted_interval_certificate.json').read_text(encoding='utf-8'))
    for name, sha in r['input_sha256'].items():
        if digest(name) != sha:
            raise ArithmeticError('changed presieve input: ' + name)
    # All nine stored finite constructions, independent of their GRH AP costs.
    cases = []
    # Begin with the two best already packaged constructions.  A larger search
    # needs separately regenerated finite masses, rather than guessed budgets.
    for pre in r['GRH_positivity']['presieve'], r['GRH_at_31000']['presieve']:
        if any(c['presieve']['u'] == pre['u'] for c in cases):
            continue
        # Match an already certified finite case, carrying separate masses.
        for source in ['conditional_chen_prefix_block_endpoint_certificate.json',
                       'conditional_chen_prefix_block_refinement.json']:
            data = json.loads((ROOT / source).read_text(encoding='utf-8'))
            for p in data['finite_cases']:
                if p['T'] == pre['T'] and p['exact_cutoff'] == pre['exact_cutoff'] and p['u'] == pre['u'] and p['log_bin_width'] == pre['log_bin_width']:
                    # Use the published budgets attached to these two cases.
                    selected = r['GRH_positivity'] if pre['u'] == r['GRH_positivity']['presieve']['u'] else r['GRH_at_31000']
                    cases.append({k: selected[k] for k in ['presieve', 'epsilon', 'eta_plus', 'eta_minus', 'C1', 'C2', 'effective_log_support_cost_ball']})
                    break
            else:
                continue
            break
        else:
            raise ArithmeticError('missing finite case')
    return cases


def prime_character_v(logX):
    # New reviewed-lemma target: |pi(x,chi)| <= 10^7 x^(61/64)
    # log^2(50rx); r <= log^10 X <= log^10 x and log(50rx)<2log x.
    return arb(40_000_000) * logX**17 * (-3 * logX / 64).exp()


def ap_c4(L):
    # BJS Lemma 30's large-primitive-conductor input; Lemmas 31/43 conversion.
    # The new proof bounds ALL small conductors, including their inductions,
    # by the retained 3.2e-8 allowance; it does not merely delete a zero term.
    c = L.log()
    w = L - 5*c
    E = 4*w**ball('4.5') / (w-5*w.log())**10 + 4/w**ball('5.5')
    E += 18*(-w/12).exp()/w.sqrt() + ball('2.5')*(-w/6).exp()*w**ball('5.5')
    p2 = w**2*(ball('1.1')*10*w.log()*ball('3.2e-8')/w**8 + 27*E
                + (-w/2).exp()/(2*arb(2).log()*w**8) + ball('.4')*w**3*(-w).exp())
    p1 = p2 + (ball('.67')+2*(-w/6).exp())/w**8
    p = p1*(1+1/(L**2*w**3)+1/((1-4/w)*L)) + ball('2.2')/L**2
    return p+ball('.9')*(w/2-L).exp()*L**4/(w**10*c)


def endpoint(case, logN, delta=Fraction(1, 1_000_000), a=None):
    L = arb(logN)
    c = L.log()
    cost = ball(case['effective_log_support_cost_ball'])
    # sqrt(x1)/log^10(x1), reserving sqrt(N)/L^(25/2).
    if a is None:
        scale = 10**9
        n = int(((cost+ball('12.5')*c)/L*scale).ceil().unique_fmpz())
        a = Fraction(n, scale)
    alpha = ball(a)
    t = ball(Fraction(1, 2)) - alpha
    s, k = 8*t, 8*(t-ball(Fraction(1,3)))
    eps, ep, em = [ball(case[n]) for n in ['epsilon', 'eta_plus', 'eta_minus']]
    d, v, ell = ball(delta), ball('1e-8'), (1+ball(delta)).log()
    if not (alpha*L-ball('12.5')*c > cost and 3<s<4 and 1<k<2):
        return None
    assert (case['C1'], case['C2']) == (106,108)
    assert eps < ball('.0001')
    assert L/8 > ball(case['presieve']['T']).log()
    assert c > ball('10.4') and (L-5*c).log() > ball('10.4')
    assert alpha-ball('12.5')/L > 0
    # Large endpoints and Euler tails.  The new gap also proves the reciprocal
    # and product bounds, without any RH prime-count or Mertens input.
    assert L/8 > 4000
    # New gap-based pi-li bound and reciprocal-prime tail, valid uniformly x>=z.
    zlog = L/8
    r = ball(Fraction(3,64))
    Az = zlog + arb(50).log()
    pi_error_z = arb(10_000_000)*(-r*zlog).exp()*Az**2
    reciprocal_error_z = arb(10_000_000)*(-r*zlog).exp()*(Az**2*(1+1/r)+2*Az/r**2+2/r**3)
    product_loss = (reciprocal_error_z+4*(-zlog).exp()).expm1()
    product_loss += 4*L*(-zlog).exp()
    assert product_loss < v
    assert reciprocal_error_z < 1/L**2
    assert pi_error_z * zlog < 1/L**2
    assert 17/zlog < r and 2/Az < r
    correction = acb.integral(lambda x, analytic: (x-2).log(analytic=analytic)/(x-1), acb(3), acb(s)).real
    Fs, fs, hs = 2*G*(1+correction)/s, 2*G*(s-1).log()/s, 3*(-s).exp()/s
    large_lower, large_upper = fs-108*eps*arb(2).exp()*hs, Fs+106*eps*arb(2).exp()*hs
    low = 8*(1-v)*((1+ep)*large_lower-(ep+em)*large_upper)
    I = G/(4*t)*(6*(3-8*alpha)/(3-18*alpha)).log()
    FM, J, b0 = 2*G/k+106*eps, (arb(8)/3).log(), t-arb(1)/4
    Hh = acb.integral(lambda b, analytic: (2-8*(acb(t)-b)).exp()/b, acb(arb(1)/8), acb(b0)).real+((arb(1)/3)/b0).log()
    upper = 4*(1+ep)*(1+2/L)*(1+v)/(1-(-L/8).exp())*(I+106*eps*Hh+2*FM/L**2)
    cb = acb.integral(lambda b, analytic: (2-3*b).log(analytic=analytic)/(b*(1-b)), acb(arb(1)/8), acb(arb(1)/3)).real
    switch = (1+ep)*(1+v)/2*(2*G/t+318*eps)*(1+d)*(1+4/L)*(cb+2/L**2+(J+2/L**2)*(6*ell/L+8/L**2))
    c4 = ap_c4(L)
    assert c4*L**2 < ball('2.21')
    # r(q) is charged against the SAME averaged AP sum, not a GRH pointwise bound.
    # The 2.21/L^2 envelope is the quantity propagated on the half-line.
    AP = (arb('1.5')+arb('.275')*L)*ball('2.21')/(UMIN*L**3)+4*(1+ep)*(1+v)*FM*ball('2.21')/L**4
    # Triple rectangles: X=N/w_j, Y in [N^(1/8),N^(1/3)], XY<= (1+d)N.
    # Fixed cutoff D0=log(X3)^10; X3=N^(1/8)/2.
    x3 = L/8-arb(2).log()
    v3 = prime_character_v(x3)
    assert arb(50).log()+10*x3.log()<x3
    assert 10*x3.log()<L/8
    assert 10*x3.log()>arb(10**9).log()
    assert 17/x3<arb(3)/64
    m = 39*v3+108*x3**16*(-x3).exp()/x3.log()+26*L**5/x3**10
    m += 88*L**5*((-L/3).exp()+(-L/16).exp())+106/x3**6
    boxes = 5*L/(24*ell)+1
    rect = (1+d)*512*m*boxes/(2*UMIN*L)
    # Missing Euler factors in the rectangle main mass; omega(d)<L.
    # Before normalization: per box <= (1+d)*1.1*N*L^2/z.
    # Division by 2*U_N*N/L^2 therefore introduces another L^2.
    rect_deleted = (1+d)*arb('1.1')*L**4*(-L/8).exp()*boxes/(2*UMIN)
    finite = L**2/UMIN*(2*(-L/8).exp()+(-2*L/3).exp()+(-L).exp())
    margin = low-upper-switch-AP-rect-rect_deleted-finite
    return {'log_N0':logN,'log_log_N0_ball':str(c),'alpha':str(a),'delta':str(delta),
        'presieve_u':case['presieve']['u'],'cost_ball':str(cost),'support_slack_ball':str(alpha*L-ball('12.5')*c-cost),
        **{n:str(x) for n,x in [('lower',low),('half_upper',upper),('half_switch',switch),
        ('AP',AP),('c4_times_L2',c4*L**2),('rectangle',rect),('rectangle_deletions',rect_deleted),
        ('finite',finite),('margin',margin),('new_v3',v3),('product_loss',product_loss),
        ('reciprocal_error_at_z',reciprocal_error_z)]}}


def run():
    cases=inputs()
    rows=[]
    for case in cases:
        for delta in [Fraction(1,10**n) for n in [4,5,6,7,8]]:
            lo,hi=40000,80000
            while hi-lo>1:
                mid=(hi+lo)//2
                row=endpoint(case,mid,delta)
                if row is not None and ball(row['margin'])>ball('.00001'): hi=mid
                else: lo=mid
            row=endpoint(case,hi,delta)
            if row is not None and ball(row['margin'])>ball('.00001'):
                before=endpoint(case,hi-1,delta)
                assert before is None or ball(before['margin'])<ball('.00001')
                row['preceding_recipe_margin']=before['margin'] if before else None
                rows.append(row)
    best=min(rows,key=lambda r:r['log_N0'])
    path=Path(__file__).name
    r={'status':'STRICT ENDPOINT ARITHMETIC; new analytic extension requires review',
        'zero_location_external_input':'all Dirichlet nontrivial zeros have Re rho <= 7/8',
        'external_commit':'adc7f1241b42e322a6451854ab7e4b4c146bf78a','precision_bits':ctx.prec,
        'selected':best,'candidate_rows':rows,
        'not_checked':['new uniform explicit-formula lemma','rectangle applicability','half-line analytic propagation',
                       'fresh regeneration of all finite presieve masses','independent mathematical review','global optimality'],
        'input_sha256':{n:digest(n) for n in [path,'conditional_chen_accepted_interval_certificate.json',
                       'qrh_source/paper.tex','qrh_source/Nonvanishing.lean',
                       'bjs_source/Explicit_Chen_-_New.tex','conditional_chen_grh_bridge.tex']}}
    (ROOT/'qrh_chen_certificate.json').write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(best,indent=2))


if __name__=='__main__':
    run()
