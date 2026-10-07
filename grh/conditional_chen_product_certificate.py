"""Arb certificate for a fixed-upper-endpoint sieve product.

This is a product lemma, NOT a certificate of a new Chen threshold.
The distinction between fixed and arbitrary upper endpoints is intentional.
"""
import json
import math
from pathlib import Path
from flint import arb, ctx

def primes_to(n):
    sieve = bytearray(b'\x01') * (n + 1)
    sieve[:2] = b'\x00\x00'
    for p in range(2, math.isqrt(n) + 1):
        if sieve[p]:
            sieve[p*p:n+1:p] = b'\x00' * ((n-p*p)//p+1)
    return [p for p in range(2, n+1) if sieve[p]]

ctx.prec = 192
T = 1_000_000
ps = primes_to(T)
gamma = arb.const_euler()
ct = arb(1)
for p in ps[1:]:
    ct *= 1 - arb(1)/(p-1)**2
dinf_lo = 2*(-gamma).exp()*ct*(-arb(1)/(T-1)).exp()
dinf_hi = 2*(-gamma).exp()*ct
prod = arb(1)
peak = arb(0)
peak_prime = None
checked = 0
for p in ps[1:]:
    dp = arb(p).log()*prod
    if p > 11:
        assert dp < dinf_lo
        checked += 1
        if dp > peak:
            peak, peak_prime = dp, p
    prod *= 1-arb(1)/(p-1)
assert arb(T).log()*prod < dinf_lo
r = arb('0.000002964')/arb(T).log()+arb(1)/(T-1)+arb('0.54')/T
delta = 3*arb(40)*(-arb(20)).exp()/(8*arb.pi())
epsilon = r.exp()/(1-delta)-1
assert epsilon < arb('0.000002')
# The strict-all-endpoints variant is different: here an explicit failure.
# u=17, v=19: include p=17 but exclude p=19.
local_ratio = (arb(17).log()/arb(19).log())*(arb(16)/15)
assert local_ratio > 1+arb('0.000002')
data = {
    'status': 'product lemma only; no new Chen threshold certified',
    'T':T, 'precision_bits':ctx.prec,
    'prime_peaks_checked':checked, 'maximum_peak_prime':peak_prime,
    'maximum_peak_ball':str(peak),
    'D_infinity_lower_ball':str(dinf_lo),
    'D_infinity_upper_ball':str(dinf_hi),
    'tail_log_ratio_budget_ball':str(r),
    'epsilon_ball':str(epsilon),
    'epsilon_budget':'0.000002',
    'presieve_Q':1155, 'log_Q_ball':str(arb(1155).log()),
    'range':'1<u<z, z>=exp(40), remove primes <=11',
    'assumption':'RH for the lower bound at z; unconditional finite scan and tail at u',
    'all_upper_endpoints_counterexample_u':17,
    'all_upper_endpoints_counterexample_v':19,
    'counterexample_normalized_ratio_ball':str(local_ratio),
}
dest=Path(__file__).with_suffix('.json')
dest.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps(data,ensure_ascii=False,indent=2))
