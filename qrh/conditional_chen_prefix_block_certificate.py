"""Positive log-bin bounds retaining every preceding Rosser prefix check.

Complete prime blocks share floor(log(p)/h). Passing and failing cutoffs
are relaxed in opposite directions, so an ambiguous prefix can be counted
in both. A positive block polynomial propagates all possible passing
prefixes and bounds their normalized first-failure masses.
This script does not certify external linear-sieve or distribution lemmas.
"""
import hashlib
import json
import math
import time
from fractions import Fraction
from pathlib import Path

from flint import arb, arb_poly, ctx

from conditional_chen_rankin_certificate import odd_primes_to


def integer_floor(value):
    n = value.floor().unique_fmpz()
    assert n is not None, 'Unresolved log-bin integer'
    return int(n)


def integer_ceil(value):
    n = value.ceil().unique_fmpz()
    assert n is not None, 'Unresolved log-bin integer'
    return int(n)


def grid_blocks(T,w,u,h,atoms=None):
    u,h = Fraction(u),Fraction(h)
    assert T >= 3 and 2 <= w < T and u > 3 and h > 0
    primes = odd_primes_to(T) if atoms is None else sorted(atoms)
    assert len(primes) == len(set(primes)) and all(3 <= p <= T for p in primes)
    medium = [p for p in primes if p > w]
    assert medium
    logD = arb(u.numerator)/u.denominator*arb(T).log()
    width = arb(h.numerator)/h.denominator
    limit = integer_ceil(logD/width)
    product,K = 1,0
    for p in medium:
        product *= p
        if product**u.denominator < T**u.numerator:
            K += 1
        else:
            break
    assert K >= 1
    groups = {}
    for p in reversed(medium):
        label = integer_floor(arb(p).log()/width)
        assert label > 0
        groups.setdefault(label,[]).append(p)
    blocks=[]
    for label,ps in groups.items():
        passing = integer_ceil((logD-3*arb(min(ps)).log())/width)
        failure = max(0,integer_ceil((logD-3*arb(max(ps)).log())/width)-K+1)
        assert passing > 0 and failure <= passing
        blocks.append({'label':label,'primes':ps,'passing_source_cutoff':passing,
                       'failure_source_cutoff':failure})
    return primes,medium,blocks,K,limit


def block_coefficients(primes,K):
    """D(A)=prod(1+t_p+t_p A), C(A)=(D(A)-1)/(1+A), t_p=1/(p-2)."""
    d,c = arb_poly([1]),arb_poly([])
    for p in primes:
        t = arb(1)/(p-2)
        factor = arb_poly([1+t,t])
        d = (d*factor).truncate(K+1)
        c = (c*factor).truncate(K+1)+t
    return d,c


def action(states,label,passing,limit):
    up_even,up_odd,lo_even,lo_odd = states
    shift = lambda p:p.left_shift(label).truncate(limit)
    return [shift(up_odd),shift(up_even.truncate(passing)),
            shift(lo_odd.truncate(passing)),shift(lo_even)]


def certify_blocks(T,w,u,h=Fraction(1,100),progress=True,atoms=None):
    assert ctx.prec >= 192
    started=time.monotonic()
    primes,medium,blocks,K,limit=grid_blocks(T,w,u,h,atoms)
    states=[arb_poly([1]),arb_poly([]),arb_poly([1]),arb_poly([])]
    upper,lower=arb(0),arb(0)
    one=arb(1)
    processed=0
    for number,block in enumerate(blocks,1):
        d,c=block_coefficients(block['primes'],K)
        powers=states
        updated=[arb_poly([]) for _ in range(4)]
        for j in range(K+1):
            for i in range(4):
                updated[i] += d[j]*powers[i]
            f=block['failure_source_cutoff']
            upper += c[j]*powers[0].right_shift(f)(one)
            lower += c[j]*powers[3].right_shift(f)(one)
            if j == K or not any(powers):
                break
            powers=action(powers,block['label'],block['passing_source_cutoff'],limit)
        states=updated
        processed+=len(block['primes'])
        if progress and number % 100 == 0:
            print(f'Prefix blocks: {number}/{len(blocks)}; {processed}/{len(medium)} primes; '
                  f'{time.monotonic()-started:.1f}s',flush=True)
    assert processed==len(medium) and upper >= 0 and lower >= 0
    names=[Path(__file__).name,'conditional_chen_rankin_certificate.py']
    root=Path(__file__).resolve().parent
    exact=[p for p in primes if p<=w]
    Q=math.prod(exact)
    V=math.prod((1-arb(1)/(p-1) for p in exact),start=arb(1))
    return {'status':'strict Arb positive prefix-block bounds; analytic domination lemma required',
            'precision_bits':ctx.prec,'T':T,'exact_cutoff':w,'u':str(Fraction(u)),
            'presieve_type':'pure Rosser medium weights, exact smallest-prime weights',
            'exact_presieve_Q_max':Q,'exact_product_ball':str(V),
            'logarithm_base':'e','D0':f'{T}^({Fraction(u)})',
            'log_support_cost_ball':str(arb(Q).log()+arb(Fraction(u).numerator)/Fraction(u).denominator*arb(T).log()),
            'log_bin_width':str(Fraction(h)),'bin_limit':limit,
            'maximum_actual_prefix_degree':K,'nonempty_prime_blocks':len(blocks),
            'prime_generation':'complete integer sieve' if atoms is None else 'explicit test atoms',
            'complete_odd_prime_count':len(primes),'medium_prime_count':len(medium),
            'upper_relative_Rosser_gap_ball':str(upper),
            'lower_relative_Rosser_gap_ball':str(lower),
            'relative_Rosser_gap_ball':str(upper+lower),
            'retains_all_preceding_checked_prefixes':True,
            'relaxation':'floor-log labels; possible-pass and possible-failure cutoffs; positive block coefficients',
            'deleting_primes':'complete-block positive coefficients dominate deleted blocks; actual prefix degree <=K',
            'input_sha256':{p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in names},
            'elapsed_seconds':round(time.monotonic()-started,3),
            'not_checked':['external analytic inputs','global optimality','Chen endpoint','PDF compilation']}


if __name__=='__main__':
    ctx.prec=192
    result=certify_blocks(10000,3,Fraction(7,2),Fraction(1,100))
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',
                                                encoding='utf-8',newline='\n')
    print(json.dumps(result,indent=2))
