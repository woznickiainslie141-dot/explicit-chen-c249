"""2026-10-03 research probe, not a certificate of a Chen theorem.

The product-lemma balls are rigorous (conditional on the cited published
Mertens inequalities). The BS and Wu scans are diagnostic only: they do
not certify the analytic transfer, all remainders, or the entire half-line.
Existing proof files are imported without writing to them.
"""
from __future__ import annotations
import json
import math
from pathlib import Path
import numpy as np
from flint import arb, ctx
from scipy.integrate import quad
import c249_diagnostic as wu

ctx.prec = 192
ROOT = Path(__file__).resolve().parent
T = 1_000_000
G = math.exp(0.5772156649015329)
C1, C2 = 113, 114  # BJS v6 at epsilon <= 1/1000


def primes_to(n):
    sieve = bytearray(b'\x01') * (n + 1)
    sieve[:2] = b'\x00\x00'
    for p in range(2, math.isqrt(n) + 1):
        if sieve[p]:
            sieve[p*p:n+1:p] = b'\x00' * ((n-p*p)//p+1)
    return [p for p in range(2, n+1) if sieve[p]]


def product_bounds(cutoffs):
    ps = primes_to(T)
    tails = {}
    s2, st = arb(0), arb(0)
    j = len(ps) - 1
    for u in sorted(cutoffs, reverse=True):
        while j >= 0 and ps[j] > u:
            p = ps[j]
            s2 += arb(1) / (p*p)
            st += arb(1) / (p*(p-2))
            j -= 1
        r = arb('6.836e-6') / (12 * arb(10).log())
        r += arb('2.964e-6') / arb(u).log()
        r += st + arb(1)/(T-1)
        r += arb('0.54') * (s2 + arb(1)/T)
        eps = r.exp() - 1
        # Outward decimal budget; strict Arb comparison, no float decision.
        q = math.ceil(float(eps.upper()) * 10**12) + 2
        bound = arb(q) / 10**12
        assert eps < bound and bound < arb(1)/1000
        theta = sum((arb(p).log() for p in ps if p <= u), arb(0))
        assert theta < u
        tails[u] = {
            'epsilon_ball': str(eps),
            'epsilon_budget_decimal': f'{q/10**12:.12f}',
            'theta_ball': str(theta),
            'Q_support_log_budget': u,
        }
    return tails


def bs_diagnostic(c, u, eps):
    L = math.exp(c)
    a = (u + 5*c) / L  # log Q <= u, H=sqrt(N)/L^5
    if not 0 < a < 1/24:
        return None
    s, k = 4-8*a, 8*(1/6-a)
    v = .62*L*math.exp(-L/16)
    low = 8*(1-v)*(2*G*math.log(s-1)/s-eps*C2*3*math.exp(2-s)/s)
    upper = 4*(1+2/L)*(1+v)/(1-math.exp(-L/8)) * (
        G/4*((math.log(6)+math.log((3-8*a)/(3-18*a)))/(.5-a)
             +512/(k*L*L))
        +(math.log(8/3)+64/L**2)*eps*C1)
    delta = .0001
    integral = quad(lambda b: math.log(2-3*b)/(b*(1-b)),1/8,1/3)[0]
    # Li upper envelope, replacing the hard-coded very-large-N constants.
    mass = 1+3/L+27/L**2 + L/3*math.exp(-L/6)/(8*math.pi)
    switch = .5*(1+.06*L*math.exp(-L/6))*(2*G/(.5-a)+3*eps*C1)*(1+delta)*mass*(
        integral+64*math.log(26/21)/L**2
        +(math.log(8/3)+64/L**2)*(3*math.log1p(delta)/L+27/L**2))
    Umin=math.exp(-0.5772156649015329)  # U_N=2 e^{-gamma} C_N, C_N>1/2
    ap = .429/L**2/Umin + .429*(.55+1/L)/(2*Umin*L)
    ap += 4*(1+v)*(2*G/k+eps*C1)*.429/L**3
    # m <= 10 is an unproved diagnostic envelope in this script.
    bilinear = .5/Umin*math.exp(-L/48)*(1+delta)*L**5/math.log1p(delta)*(
        .046*10*(1+delta)**(-1/6)+.159*math.exp(-5*L/48))
    finite = (2*math.exp(-L/8)+math.exp(-2*L/3))*L**2/Umin
    return {'c':c,'u':u,'alpha':a,'epsilon':eps,
            'lower':low,'half_upper':upper,'half_switch':switch,
            'diagnostic_ap':ap,'diagnostic_bilinear':bilinear,
            'diagnostic_margin':low-upper-switch-ap-bilinear-finite}


def wu_diagnostic(c, u, eps):
    wu.QBAR=u
    wu.ETA=eps
    wu.EXC_EPS=0
    wu.EXC_DENSITY=1
    wu.Fplus=lambda s: wu.sieve_F(s)+eps*C1*math.e**2*wu.h(s)
    wu.exceptional_lower=lambda s: wu.sieve_f(s)-eps*C2*math.e**2*wu.h(s)
    x=np.array([191/2000,383/4000,6047/20000,121/400])
    r=wu.rows(c,*x)
    return {'c':c,'u':u,'epsilon':eps,'theta':r['theta'],
            'margin_before_new_errors':wu.master(c,x),
            'status':'diagnostic only; no c249 error budgets transferred'}


def main():
    cutoffs=[400,600,800,1000,1500,2000,3000,5000,10000,20000]
    products=product_bounds(cutoffs)
    bs=[]
    for c in [10,11,12,12.5,13,13.5,14]:
        rows=[bs_diagnostic(c,u,float(products[u]['epsilon_budget_decimal'])) for u in cutoffs]
        rows=[r for r in rows if r is not None]
        bs.append(max(rows,key=lambda r:r['diagnostic_margin']))
    ws=[wu_diagnostic(c,400,float(products[400]['epsilon_budget_decimal']))
        for c in [9,9.5,10,10.5,11,12]]
    B=4*10**18
    support={'verification_B':B,'log_B':math.log(B),'c_B':math.log(math.log(B)),
             'B_one_third':B**(1/3),'B_one_twelfth':B**(1/12),
             'Q400_log_ball':products[400]['theta_ball'],
             'sqrt_B_log':math.log(B)/2}
    # Arithmetic for the proposed tighter transfer budget, conditional on
    # verifying the original kernel table for the new analytic argument.
    L11=arb(11).exp()
    d11=arb('0.000003')/(arb(191)/2000*L11)+2*(-arb(191)/2000*L11).exp()
    h11=arb('1.0001').log()/L11
    transfer=arb(17_400_000)*d11+arb(150_000)*h11+arb('.003')
    assert transfer < arb('.0125')
    transfer_certificate={
        'c':11,'ordinary_and_box_d_coefficient':83496,
        'two_switched_rows_d_coefficient':17280000,
        'rounded_d_coefficient':17400000,
        'h_coefficient':150000,'inflation_budget':'.003',
        'normalized_ball':str(transfer),
        'status':'rigorous arithmetic only; new applicability of kernel table unproved'}
    AB=arb(4)*arb(10)**18
    LB=AB.log()
    HB=AB**(arb(3)/4)
    re=1/LB**2
    deletions=arb('1.3841')*LB**3/LB.log()*(2/AB.sqrt()+arb('1.1')*HB.log()/AB)
    unit_exclusion=LB**2/AB  # remove eta=1 from the actual Chen conclusion
    gate=arb('.001')-re-deletions-unit_exclusion
    assert gate > arb('.000425')
    eh_gate={'explicit_extra_assumption':'sum_{d<=N^(3/4)} mu(d)^2 max_a |pi(N;d,a)-pi(N)/phi(d)| <= N/log(N)^4 for every N>=4e18',
             'unproved_arithmetic_module':'pointwise lower weights |lambda|<=1 with sum lambda/phi >= .001/log(N) for every even N>=4e18',
             'AP_bound_at_B':str(re),'deleted_primes_bound_at_B':str(deletions),
             'eta_equals_one_exclusion_at_B':str(unit_exclusion),
             'conditional_gate_margin_at_B':str(gate),
             'status':'WITHDRAWN: the assumed AP bound fails at B by prime-modulus integrality; see conditional_chen_integrality_certificate.py'}
    result={'warning':'Not a new GRH Chen theorem. Product bounds only use rigorous Arb arithmetic.',
            'sources':{'JS':'https://arxiv.org/html/2208.01229#S3',
                       'BS_v1':'https://arxiv.org/html/2211.08844',
                       'BJS_v6':'https://arxiv.org/abs/2207.09452v6'},
            'product_bounds':products,'BS_diagnostics':bs,'Wu_diagnostics':ws,
            'bridge_support':support,'proposed_transfer_arithmetic':transfer_certificate,
            'quantitative_EH_bridge_gate':eh_gate}
    dest=ROOT/'conditional_chen_probe.json'
    dest.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2,ensure_ascii=False))


if __name__=='__main__':
    main()
