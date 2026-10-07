"""Positive six-state counts for the union of the two medium Rosser supports.

Each candidate subset has one relaxed alive-flag state, so the union is
counted once. Complete-block cutoffs are retained after prime deletion.
"""
import hashlib
import json
import math
import time
from fractions import Fraction
from pathlib import Path
from flint import arb,arb_poly,ctx
from conditional_chen_rankin_certificate import odd_primes_to
from conditional_chen_prefix_block_certificate import integer_floor,integer_ceil
from conditional_chen_sharp_support_certificate import support_ratio
from conditional_chen_binned_support_certificate import ball

ROOT=Path(__file__).resolve().parent
STATES=[(1,0),(1,1),(2,0),(2,1),(3,0),(3,1)]


def append(states,label,cutoff,limit):
    ue,uo,le,lo,be,bo=states
    low=lambda p:p.truncate(cutoff)
    high=lambda p:p.right_shift(cutoff).left_shift(cutoff)
    shift=lambda p:p.left_shift(label).truncate(limit)
    return [shift(uo+high(bo)),shift(low(ue)),
            shift(low(lo)),shift(le+high(be)),
            shift(low(bo)),shift(low(be))]


def certify(T,w,u,h=Fraction(1,100),available=None,progress=False,return_states=False):
    assert ctx.prec>=192
    u,h=map(Fraction,[u,h]);assert u>3 and 0<h<1
    started=time.monotonic()
    primes=[p for p in odd_primes_to(T) if p>w]
    assert primes
    if available is None:available=set(primes)
    else:
        available=set(available);assert available<=set(primes)
    kappa,witness=support_ratio(primes)
    logD=ball(u)*arb(T).log();logH=logD+ball(kappa).log()
    J=integer_ceil(logH/ball(h));assert J>0
    product,K=1,0
    for p in primes:
        product*=p
        if product**u.denominator*kappa.denominator**u.denominator < T**u.numerator*kappa.numerator**u.denominator:
            K+=1
        else:break
    groups={}
    for p in reversed(primes):groups.setdefault(integer_floor(arb(p).log()/ball(h)),[]).append(p)
    states=[arb_poly([]) for _ in STATES];states[4]=arb_poly([1])
    for number,(label,ps) in enumerate(groups.items(),1):
        assert label>0
        cutoff=integer_ceil((logD-3*arb(min(ps)).log())/ball(h));assert cutoff>0
        size=sum(p in available for p in ps)
        powers=states;updated=[arb_poly([]) for _ in STATES]
        for j in range(min(K,size,(J-1)//label)+1):
            coefficient=math.comb(size,j)
            for index in range(6):updated[index]+=coefficient*powers[index]
            powers=append(powers,label,cutoff,J)
            if not any(powers):break
        states=updated
        if progress and number%100==0:
            print(f'Accepted support blocks {number}/{len(groups)}; {time.monotonic()-started:.1f}s',flush=True)
    totals=[p(arb(1)) for p in states]
    count=sum(totals,arb(0));assert count>=1
    exact=odd_primes_to(w);Q=math.prod(exact)
    fraction=arb(2**len(exact))/Q*count*(-logH).exp()
    names=[Path(__file__).name,'conditional_chen_rankin_certificate.py',
        'conditional_chen_prefix_block_certificate.py','conditional_chen_sharp_support_certificate.py',
        'conditional_chen_binned_support_certificate.py']
    report={'status':'strict positive union-support count; general six-state domination lemma required',
        'precision_bits':ctx.prec,'T':T,'exact_cutoff':w,'u':str(u),'kappa':str(kappa),
        'support_maximizer':witness,'sigma':'0','log_bin_width':str(h),
        'log_medium_support_ball':str(logH),'maximum_actual_degree':K,'label_limit':J,
        'medium_prime_count':len(primes),'available_medium_prime_count':len(available),
        'nonempty_blocks':len(groups),'count_upper_ball':str(count),
        'modulus_count_fraction_ball':str(fraction),
        'state_count_balls':[{'alive_flags':flags,'selected_parity':parity,'count_ball':str(value)}
            for (flags,parity),value in zip(STATES,totals)],
        'state_convention':'1 upper alive; 2 lower alive; 3 both alive; parity before append',
        'prime_generation':'complete integer sieve','retains_all_preceding_checked_prefixes':True,
        'union_without_duplication':True,'complete_cutoffs_retained_after_deletion':True,
        'elapsed_seconds':round(time.monotonic()-started,3),
        'input_sha256':{n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names},
        'not_checked':['general domination proof','external analytic inputs','Chen endpoint','full splice']}
    return (report,states) if return_states else report


if __name__=='__main__':
    ctx.prec=192
    r=certify(5000000,5,Fraction(87,16),progress=True)
    Path(__file__).with_suffix('.json').write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:r[k] for k in ['T','u','count_upper_ball','modulus_count_fraction_ball','state_count_balls','elapsed_seconds']},indent=2))
