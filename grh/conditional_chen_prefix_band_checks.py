"""Independent rational checks using formal atoms of sizes 2^j.

Their base-two logarithms, densities, product ratios, shifted moments and
polynomial optimizers are exact rationals. These atoms check the general
finite recurrences and first-failure logic; they are not declared primes.
"""
import json
import math
from fractions import Fraction
from pathlib import Path
from flint import arb,ctx
from conditional_chen_prefix_band_certificate import certify_band
from conditional_chen_rosser_checks import first_failure_masses


def rational_candidate(m,R):
    m0,m1,m2,m3,m4 = m
    v0,v1 = m1,R*m1-m2
    A,B,C = m2,R*m2-m3,R*R*m2-2*R*m3+m4
    values = [m0]
    for v,den in [(v0,A),(v1,C)]:
        if v < 0 and den > 0:
            c = -v/den
            value = m0-v*v/den
            assert value == m0+2*c*v+c*c*den >= 0
            values.append(value)
    det = A*C-B*B
    assert det >= 0
    if det > 0:
        e,b = (B*v1-C*v0)/det,(B*v0-A*v1)/det
        if e >= 0 and b >= 0:
            value = m0+v0*e+v1*b
            assert value == m0+2*e*v0+2*b*v1+e*e*A+2*e*b*B+b*b*C >= 0
            values.append(value)
    return min(values)


def direct_bound(atoms,T,u,ordinary,bands,M):
    logT = T.bit_length()-1
    logD = Fraction(u)*logT
    assert logD.denominator == 1
    D0 = 2**logD.numerator
    totals = [Fraction() for _ in range(M+2)]
    for r in atoms:
        larger = [p for p in atoms if p > r]
        lr = r.bit_length()-1
        base = Fraction(1,r-2)*math.prod(Fraction(p-1,p-2) for p in larger)
        moments = {sigma:[[Fraction() for _ in range(5)] for _ in range(M+2)]
                   for sigma in set(ordinary)|set(bands)}
        for mask in range(1 << len(larger)):
            chosen = [p for i,p in enumerate(larger) if mask & (1 << i)]
            k = len(chosen)
            if k < math.floor(u-2) or r**3*T**k < D0:
                continue
            bucket = k if k < M else M+(k % 2)
            density = Fraction(1,math.prod(p-1 for p in chosen))
            product = math.prod(chosen)
            delta = 3*lr+sum(p.bit_length()-1 for p in chosen)-u*logT
            for sigma,groups in moments.items():
                tilt = Fraction(r**3*product,D0)**sigma.numerator
                for j in range(5):
                    groups[bucket][j] += density*tilt*delta**j
        for k in range(M+2):
            candidates = [moments[s][k][0] for s in ordinary]
            candidates += [rational_candidate(moments[s][k],2*lr) for s in bands]
            totals[k] += base*min(candidates)
    upper = sum((v for k,v in enumerate(totals[:M]) if k % 2 == 0),Fraction())+totals[M]
    lower = sum((v for k,v in enumerate(totals[:M]) if k % 2),Fraction())+totals[M+1]
    return totals,upper,lower


def run():
    ctx.prec = 192
    complete_atoms = [2**j for j in range(2,10)]
    ordinary = [Fraction(0),Fraction(1),Fraction(2)]
    bands = [Fraction(0),Fraction(1)]
    configs = []
    deletion_count = 0
    band_count = 0
    two_prime_count = 0
    for T,u in [(512,Fraction(4)),(512,Fraction(5)),(512,Fraction(6)),
                (256,Fraction(7,2))]:
        atoms = [p for p in complete_atoms if p <= T]
        logT = T.bit_length()-1
        logD = u*logT
        assert logD.denominator == 1
        D0 = 2**logD.numerator
        for M in [math.ceil(u),math.ceil(u)+1,12]:
            certified = certify_band(T,3,u,ordinary,bands,M,False,
                                     atoms=atoms,log_base=2)
            totals,upper,lower = direct_bound(atoms,T,u,ordinary,bands,M)
            for a,b in zip(certified["degree_contributions_ball"],totals):
                exact = arb(b.numerator)/b.denominator
                assert abs(arb(a)-exact) < arb("1e-42")
            for mask in range(1 << len(atoms)):
                subset = [p for i,p in enumerate(atoms) if mask & (1 << i)]
                V = math.prod(Fraction(p-2,p-1) for p in subset)
                even,odd = first_failure_masses(subset,D0)
                _,sub_upper,sub_lower = direct_bound(subset,T,u,ordinary,bands,M)
                assert odd/V <= sub_upper <= upper
                assert even/V <= sub_lower <= lower
                deletion_count += 1
            configs.append({"T":T,"u":str(u),"degree_cutoff":M,
                            "strict_optimizers":certified["band_optimizers_with_strict_signs"]})
        for mask in range(1,1 << len(atoms)):
            selected = sorted((p for i,p in enumerate(atoms) if mask & (1 << i)),
                              reverse=True)
            d = math.prod(selected)
            r = selected[-1]
            prefix = 1
            previous = True
            for n,p in enumerate(selected,1):
                prefix *= p
                if n < len(selected) and n % 2 == len(selected) % 2:
                    previous &= prefix*p*p < D0
            if previous and d*r*r >= D0:
                assert d < D0
                delta = 3*(r.bit_length()-1)+sum(p.bit_length()-1 for p in selected[:-1])-logD
                assert 0 <= delta < 2*(r.bit_length()-1)
                band_count += 1
                two_prime_count += int(len(selected) == 2)
    assert all(c["strict_optimizers"] > 0 for c in configs) and two_prime_count > 0
    return {"status":"PASS: exact formal-atom checks only",
            "precision_bits":ctx.prec,"test_atoms":complete_atoms,"atoms_are_primes":False,
            "explicit_coefficient_configurations":len(configs),
            "first_failure_and_deleted_subset_checks":deletion_count,
            "first_failure_product_band_checks":band_count,
            "two_prime_first_failure_checks":two_prime_count,
            "configs":configs,
            "not_checked":"external analytic sieve inputs and full Chen bridge"}


if __name__ == "__main__":
    result = run()
    Path(__file__).with_suffix(".json").write_text(json.dumps(result,indent=2)+"\n",
                                                encoding="utf-8")
    print(json.dumps(result,indent=2))
