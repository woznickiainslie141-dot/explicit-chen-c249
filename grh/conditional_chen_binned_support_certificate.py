"""Positive logarithmic-bin upper bounds for square-free smooth supports.

Counts all medium-prime subsets below H0; no Rosser acceptance condition
is dropped from the actual sieve weights. Complete prime blocks and
nonnegative arithmetic bound the count directly, without cancellation.
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

ROOT=Path(__file__).resolve().parent


def ball(value):
    value=Fraction(value)
    return arb(value.numerator)/value.denominator


def count_bound(T,w,u,kappa,sigma,h,return_polynomial=False):
    assert ctx.prec>=192
    u,kappa,sigma,h=map(Fraction,[u,kappa,sigma,h])
    assert u>3 and 0<kappa<1 and sigma>=0 and 0<h<1
    started=time.monotonic()
    primes=[p for p in odd_primes_to(T) if p>w]
    logH=ball(u)*arb(T).log()+ball(kappa).log()
    assert logH>0
    width=ball(h);limit=integer_ceil(logH/width)
    product,K=1,0
    for p in primes:
        product*=p
        if product**u.denominator*kappa.denominator**u.denominator < T**u.numerator*kappa.numerator**u.denominator:
            K+=1
        else:break
    assert K>=1
    groups={}
    for p in reversed(primes):
        label=integer_floor(arb(p).log()/width)
        assert label>0
        groups.setdefault(label,[]).append(p)
    state=arb_poly([1]);log_euler=arb(0)
    for label,ps in groups.items():
        if sigma==0:
            coefficients=arb_poly([math.comb(len(ps),j) for j in range(min(K,len(ps))+1)])
            log_euler+=len(ps)*arb(2).log()
        else:
            coefficients=arb_poly([1])
            for p in ps:
                weight=(-ball(sigma)*arb(p).log()).exp()
                coefficients=(coefficients*arb_poly([1,weight])).truncate(K+1)
                log_euler+=(1+weight).log()
        updated=arb_poly([])
        for j in range(min(K,(limit-1)//label)+1):
            updated+=coefficients[j]*state.left_shift(j*label).truncate(limit)
        state=updated
    weighted_count=arb(0)
    for label in range(len(state)):
        ceiling=width*(label+K)
        # Either ceiling is valid for an actual subset. An unresolved
        # comparison safely chooses logH, rather than omitting a term.
        cap=ceiling if ceiling<logH else logH
        weighted_count+=state[label]*(ball(sigma)*cap).exp()
    old=arb(2)**len(primes) if sigma==0 else (ball(sigma)*logH+log_euler).exp()
    assert 1<=weighted_count and weighted_count<=old
    Q=math.prod(p for p in odd_primes_to(w))
    exact=2**len(odd_primes_to(w))
    ratio=arb(exact)/Q*weighted_count*(-logH).exp()
    names=[Path(__file__).name,'conditional_chen_rankin_certificate.py',
           'conditional_chen_prefix_block_certificate.py']
    result={'status':'strict positive log-bin support count; general count lemma required',
        'precision_bits':ctx.prec,'T':T,'exact_cutoff':w,'u':str(u),'kappa':str(kappa),
        'sigma':str(sigma),'log_bin_width':str(h),'log_medium_support_ball':str(logH),
        'maximum_actual_degree':K,'label_limit':limit,'medium_prime_count':len(primes),
        'nonempty_blocks':len(groups),'count_upper_ball':str(weighted_count),
        'ordinary_Rankin_count_upper_ball':str(old),
        'ordinary_Rankin_to_bin_bound_ratio_ball':str(old/weighted_count),
        'log_euler_product_ball':str(log_euler),'modulus_count_fraction_ball':str(ratio),
        'prime_generation':'complete integer sieve','elapsed_seconds':round(time.monotonic()-started,3),
        'input_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in names},
        'not_checked':['general count lemma','external analytic inputs','Chen endpoint','full splice']}
    return (result,state) if return_polynomial else result


if __name__=='__main__':
    ctx.prec=192
    r=count_bound(5000000,5,Fraction(87,16),Fraction(11,169),Fraction(4,5),Fraction(1,100))
    Path(__file__).with_suffix('.json').write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(r,indent=2))
