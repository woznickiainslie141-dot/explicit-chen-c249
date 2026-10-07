"""A verified support saving for finite Rosser medium-prime weights.

For D0=T^u, u>3, the support is < kappa D0, where kappa is the
maximum of 1/p_min^2 and r/s^2 over r<s in the complete medium set.
This saves log(1/kappa) without changing weights or their masses.
No updated Chen threshold is asserted by this support-only certificate.
"""
import json
import math
from fractions import Fraction
from pathlib import Path
from flint import arb,ctx
from conditional_chen_rankin_certificate import odd_primes_to


def support_ratio(primes):
    assert primes and primes == sorted(set(primes)) and primes[0] > 2
    ratio = Fraction(1,primes[0]**2)
    winner = {"kind":"checked final prefix","minimum_prime":primes[0]}
    # For fixed s, the largest possible r is its immediate predecessor.
    for r,s in zip(primes,primes[1:]):
        value = Fraction(r,s*s)
        if value > ratio:
            ratio = value
            winner = {"kind":"unchecked final prime","r":r,"s":s}
    return ratio,winner


def finite_checks():
    checks = 0
    nonzero = 0
    summaries = []
    complete = [5,7,11,13,17,19,23,29]
    for w in [3,5]:
        atoms = [p for p in complete if p > w]
        kappa,_ = support_ratio(atoms)
        for u in [Fraction(7,2),Fraction(4),Fraction(49,8)]:
            for deleted in range(1 << len(atoms)):
                subset = [p for i,p in enumerate(atoms) if deleted & (1 << i)]
                for mask in range(1 << len(subset)):
                    chosen = sorted((p for i,p in enumerate(subset) if mask & (1 << i)),
                                    reverse=True)
                    d = math.prod(chosen)
                    for upper in [False,True]:
                        prefix = 1
                        accepted = True
                        for j,p in enumerate(chosen,1):
                            prefix *= p
                            if j % 2 == int(upper):
                                accepted &= (prefix*p*p)**u.denominator < 29**u.numerator
                        if accepted:
                            assert (d*kappa.denominator)**u.denominator < (
                                kappa.numerator**u.denominator*29**u.numerator)
                            nonzero += 1
                        checks += 1
            summaries.append({"exact_cutoff":w,"u":str(u),"kappa":str(kappa)})
    return {"status":"PASS: exact integer-power support checks",
            "coefficient_and_deleted_subset_cases":checks,
            "accepted_coefficients_checked":nonzero,"configurations":summaries}


def run():
    ctx.prec = 192
    checked = finite_checks()
    values = []
    for T,w,u in [(5000000,5,Fraction(49,8)),(22000,3,Fraction(39,8)),
                  (10000,3,Fraction(35,8)),(10000,3,Fraction(17,4))]:
        primes = odd_primes_to(T)
        medium = [p for p in primes if p > w]
        kappa,winner = support_ratio(medium)
        Q = math.prod(p for p in primes if p <= w)
        assert kappa < 1 and u > 3
        assert kappa >= Fraction(1,T*T)
        assert kappa.numerator**u.denominator*T**u.numerator > (
            kappa.denominator**u.denominator*T**u.denominator)
        original = arb(Q).log()+arb(u.numerator)/u.denominator*arb(T).log()
        saving = arb(kappa.denominator).log()-arb(kappa.numerator).log()
        sharp = original-saving
        assert sharp < original and saving > 0
        assert kappa == (Fraction(11,169) if w == 5 else Fraction(5,49))
        values.append({"T":T,"exact_cutoff":w,"u":str(u),"Q":Q,
                       "complete_medium_prime_count":len(medium),"kappa":str(kappa),
                       "maximizer":winner,"log_support_saving_ball":str(saving),
                       "original_log_support_cost_ball":str(original),
                       "sharp_log_support_cost_ball":str(sharp)})
    return {"status":"PASS: support-only certificate; not an updated Chen endpoint",
            "precision_bits":ctx.prec,"finite_checks":checked,"cases":values,
            "not_checked":"recharging the full GRH and distribution endpoint budgets"}


if __name__ == "__main__":
    result = run()
    Path(__file__).with_suffix(".json").write_text(json.dumps(result,indent=2)+"\n",
                                                encoding="utf-8")
    print(json.dumps(result,indent=2))
