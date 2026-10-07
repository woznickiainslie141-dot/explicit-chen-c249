"""Finite gates for the new seven-eighths analytic interfaces.

These gates check constants and starting inequalities, not analytic theorems.
The continuum arguments and precise cited dependencies remain in the proof note.
"""
import hashlib
import json
from fractions import Fraction
from pathlib import Path

from flint import arb, ctx
from conditional_chen_rankin_certificate import odd_primes_to
import qrh_chen_certificate as q

ROOT = Path(__file__).resolve().parent
ctx.prec = 192


def run(logN=51046, verify_prime_set=True):
    L = arb(logN)
    r = arb(3)/64
    w = L-5*L.log()
    log2 = arb(2).log()
    B = arb(50).log()
    # Disk geometry, principal normalization, BC, Cauchy, contour integral.
    assert arb(57)/64 > arb(7)/8
    assert arb(141)/128-arb(140)/128 == arb(1)/128
    assert 282*(B+4)+4 < 600*B
    assert arb(32)/3+arb(32)/61 < 12
    assert 76800+12/B < 80000
    integral_majorant = (arb(100)/81+1)*(B+log2)+1
    assert integral_majorant < 3*B
    assert 80000*3/arb.pi() < 100000
    assert 500000+2+1/(2*arb(100).log()) < 1000000
    partial_summation_factor = 1/log2+1/((arb(61)/64)*log2**2)
    assert partial_summation_factor < 4
    assert (1000000+2)*4+4 < 10000000

    # The li envelopes start at log x=20 and then follow by differentiation.
    e20 = arb(20).exp()
    li20 = arb(20).ei()
    assert e20/20*(1+arb(1)/20) < li20
    assert li20 < e20/20*(1+arb(1)/20+arb(3)/400)

    def relative_prime_error(t):
        return 10000000*(-r*t).exp()*(t+arb(50).log())**2

    # Lower A mass includes deletion of every prime divisor of N.
    A_lower_gate = L**2*relative_prime_error(L)+L**3*(-L).exp()/log2
    A_upper_gate = L**2*relative_prime_error(L)+3/L
    pi_upper_gate = L*(L/3)*relative_prime_error(L/3)+27/L
    assert A_lower_gate < 1 and A_upper_gate < 1 and pi_upper_gate < 1
    assert 3/(L/3) < r  # all three positive error envelopes decrease

    # All small primitive conductors, even when induced to large moduli.
    # After exact conductor decomposition the whole low part is <=
    # 1.1e6 y^(61/64) log(y)^12.  The retained BJS allowance suffices.
    low_conductor_gate = arb('1.1e6')*w**20*(-r*w).exp()
    assert low_conductor_gate < arb('3.2e-8')
    assert 20/w < r
    assert w.log() > arb('10.4')
    assert w/2-10*w.log() > 10*w.log()  # Q1 < H
    assert arb(50).log()+10*w.log() < w

    # Explicit constants in the fixed-cutoff rectangle; the number of dyadic
    # blocks is <= log(D)/log(2) since the cutoff is at least 2.
    constants = {
        'E1_constant': arb('2.42')*16,
        'E2_constant': arb('1.21')*arb('1.3841')*64,
        'E3_constant': arb('1.1')*16/log2,
        'E4_constant': arb('1.1')*16*arb(12).sqrt()/log2,
        'E5_constant': arb('1.1')*4*24,
    }
    for name, budget in [('E1_constant',39),('E2_constant',108),
                         ('E3_constant',26),('E4_constant',88),('E5_constant',106)]:
        assert constants[name] < budget

    # Normalized BJS polynomial terms and fixed-cutoff rectangle derivatives.
    # Dots here are derivatives with respect to log L.
    a = (L-5)/w
    b = (L-5)*(1-5/w)/(w-5*w.log())
    assert 1 < a < arb('1.01') and b > 1
    assert 2-6+arb('1.01')/w.log() < 0
    assert 2+arb('6.5')*arb('1.01')-10 < 0
    assert 2+arb('1.5')*a-(L-5)/12 < 0
    assert 2+arb('7.5')*a-(L-5)/6 < 0
    assert 6+(L-5)/2-L-10*a-1/L.log() < 0
    x3 = L/8-log2
    assert 17/x3 < r and 16/x3 < 1 and x3 > 1
    assert 5-L/16 < 0

    case=q.selected_case()
    fixed=q.endpoint(case,logN,Fraction(1,1000000),Fraction(2158047,500000000))
    assert fixed is not None
    delta=arb(1)/1000000
    boxes=5*L/(24*(1+delta).log())+1
    # Independently reconstruct the normalization of the deleted-prime charge.
    raw_over_N=(1+delta)*arb('1.1')*L**2*(-L/8).exp()
    normalized=raw_over_N*L**2/(2*q.UMIN)*boxes
    assert abs(arb(fixed['rectangle_deletions'])/normalized-1)<arb('1e-45')
    assert arb(fixed['margin'])>arb('9e-7')

    result = {
        'status': 'PASS: finite constants and starting gates only',
        'log_N0': logN,
        'analytic_scope': 'continuum propagation requires the written derivative arguments',
        'BC_bound_at_minimum_M': str(282*(B+4)+4),
        'smoothing_integral_constant': str(80000*3/arb.pi()),
        'partial_summation_factor': str(partial_summation_factor),
        'A_lower_gate': str(A_lower_gate),
        'A_upper_gate': str(A_upper_gate),
        'pi_upper_gate': str(pi_upper_gate),
        'all_small_primitive_conductors_gate': str(low_conductor_gate),
        'rectangle_constants': {k:str(v) for k,v in constants.items()},
        'corrected_rectangle_deletion_ball':str(normalized),
        'fixed_endpoint_margin_ball':fixed['margin'],
    }
    if verify_prime_set:
        T = 5000000
        primes = odd_primes_to(T)
        medium = [p for p in primes if p>5]
        # Exact integer comparisons; consecutive primes maximize r/s^2.
        assert 169 <= 11*medium[0]**2
        assert all(169*p <= 11*q*q for p,q in zip(medium,medium[1:]))
        assert (11,13) in zip(medium,medium[1:])
        C = arb(1)
        for p in primes:
            C *= 1-arb(1)/(p-1)**2
        C_lower = C*(1-arb(1)/T)
        assert C_lower > arb('.66')
        result['fresh_exact_support_and_Euler_product'] = {
            'complete_odd_prime_count':len(primes),
            'medium_prime_count':len(medium),
            'kappa_exact':'11/169',
            'kappa_attaining_pair':[11,13],
            'C_infinity_lower_ball':str(C_lower),
            'tail':'product over all integers n>T equals (T-1)/T',
        }
    names=[Path(__file__).name,'conditional_chen_rankin_certificate.py','qrh_chen_certificate.py',
           'qrh_chen_research_20261007.md']
    result['input_sha256']={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names}
    result['not_checked']=['external seven-eighths theorem proof',
                           'general contour, sieve and distribution arguments',
                           'global threshold optimality']
    return result


if __name__ == '__main__':
    result=run()
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
