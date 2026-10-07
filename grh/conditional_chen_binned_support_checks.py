"""Independent subset checks for binned support counts and integer roughness."""
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
from flint import arb,ctx
from conditional_chen_rankin_certificate import odd_primes_to
from conditional_chen_binned_support_certificate import count_bound,ball
from conditional_chen_prefix_block_certificate import certify_blocks

ROOT=Path(__file__).resolve().parent


def small_bins():
    configurations=[];subsets=coefficient_checks=0
    for w,kappa in [(3,Fraction(5,49)),(5,Fraction(11,169))]:
        ps=[p for p in odd_primes_to(29) if p>w]
        for u in [Fraction(7,2),Fraction(4),Fraction(49,8)]:
            for sigma in [Fraction(0),Fraction(1)]:
                for h in [Fraction(9,10),Fraction(1,2),Fraction(1,10),Fraction(1,100)]:
                    result,production=count_bound(29,w,u,kappa,sigma,h,True)
                    K=result['maximum_actual_degree'];J=result['label_limit']
                    labels=[int((arb(p).log()/ball(h)).floor().unique_fmpz()) for p in ps]
                    reference={};actual=maximum=0
                    for mask in range(1<<len(ps)):
                        chosen=[p for i,p in enumerate(ps) if mask>>i&1]
                        label=sum(b for i,b in enumerate(labels) if mask>>i&1)
                        n=math.prod(chosen)
                        inside=n**u.denominator*kappa.denominator**u.denominator<29**u.numerator*kappa.numerator**u.denominator
                        if inside:actual+=1;maximum=max(maximum,len(chosen))
                        groups={b:sum(bool(mask>>i&1) for i,c in enumerate(labels) if c==b) for b in set(labels)}
                        if label<J and all(degree<=K for degree in groups.values()):
                            reference[label]=reference.get(label,Fraction())+(Fraction(1,n) if sigma else Fraction(1))
                        subsets+=1
                    assert K==maximum
                    for label in range(max(len(production),max(reference,default=0)+1)):
                        exact=reference.get(label,Fraction())
                        assert abs(production[label]-ball(exact))<arb('1e-43')
                        coefficient_checks+=1
                    logH=ball(u)*arb(29).log()+ball(kappa).log()
                    independent=arb(0)
                    for label,value in reference.items():
                        cap=ball(h)*(label+K)
                        if not cap<logH:cap=logH
                        independent+=ball(value)*(ball(sigma)*cap).exp()
                    assert abs(independent-arb(result['count_upper_ball']))<arb('1e-40')
                    assert arb(actual)<=arb(result['count_upper_ball'])
                    configurations.append({'w':w,'u':str(u),'sigma':str(sigma),'h':str(h),'exact_count':actual,'K':K})
    return {'configurations':configurations,'enumerated_subsets':subsets,
            'exact_rational_histogram_coefficients_checked':coefficient_checks}


def integer_roughness():
    configurations=[];floor_checks=patterns=0
    for T in [29,41]:
        for w,kappa in [(3,Fraction(5,49)),(5,Fraction(11,169))]:
            all_primes=[2]+odd_primes_to(T)
            exact=[p for p in all_primes if p<=w];medium=[p for p in all_primes if p>w]
            level=T**4;coefs={}
            old=certify_blocks(T,w,Fraction(4),Fraction(1,10),False)
            B=Fraction();V=math.prod((Fraction(p-1,p) for p in medium),start=Fraction(1))
            for mask in range(1<<len(medium)):
                chosen=sorted((p for i,p in enumerate(medium) if mask>>i&1),reverse=True)
                prefix=1;accepted=True
                for m,p in enumerate(chosen,1):
                    prefix*=p
                    if m%2 and prefix*p*p>=level:accepted=False;break
                if not accepted:continue
                d=math.prod(chosen);mu=(-1)**len(chosen)
                B+=Fraction(mu,d)
                for exact_mask in range(1<<len(exact)):
                    e=math.prod(p for i,p in enumerate(exact) if exact_mask>>i&1)
                    coefs[d*e]=mu*(-1)**exact_mask.bit_count()
            normalized=(B-V)/V
            assert 0<=normalized and ball(normalized)<=arb(old['upper_relative_Rosser_gap_ball'])
            H=math.prod(exact)*kappa*level
            assert all(d<H for d in coefs) and len(coefs)<H
            product=math.prod((Fraction(p-1,p) for p in all_primes),start=Fraction(1))
            mass=sum((Fraction(c,d) for d,c in coefs.items()),Fraction())
            assert mass==product*(1+normalized)
            for D in [Fraction(31),Fraction(101,2),Fraction(1000),Fraction(10001,2)]:
                count=sum(math.gcd(n,math.prod(all_primes))==1 for n in range(1,math.ceil(D)))
                upper=sum(c*(math.ceil(D/d)-1) for d,c in coefs.items())
                assert count<=upper<=D*mass+len(coefs)
                assert ball(D*mass+len(coefs)) < ball(D)*ball(product)*(1+arb(old['upper_relative_Rosser_gap_ball']))+ball(H)
                floor_checks+=len(coefs);patterns+=math.ceil(D)-1
            configurations.append({'T':T,'w':w,'integer_upper_coefficients':len(coefs),
                                   'exact_normalized_gap':str(normalized)})
    return {'configurations':configurations,'exact_divisor_floor_terms_checked':floor_checks,
            'integer_roughness_patterns_checked':patterns}


def aggregate_prime_factors():
    """Check the entire (q,d) image, not a bound for each q separately."""
    T,z,y=11,31,100
    exact=[3,5];medium=[7,11]
    large=[p for p in odd_primes_to(z-1) if p>T]
    qs=[p for p in odd_primes_to(y-1) if p>=z]
    factors=[math.prod(p for i,p in enumerate(large) if mask>>i&1)
             for mask in range(1<<len(large))]
    exact_factors=[math.prod(p for i,p in enumerate(exact) if mask>>i&1)
                  for mask in range(1<<len(exact))]
    medium_factors=[math.prod(p for i,p in enumerate(medium) if mask>>i&1)
                   for mask in range(1<<len(medium))]
    modulus_patterns=0;rows=[]
    for D in [101,1000,10000,1000000]:
        pairs=[(q,d1) for q in qs for d1 in factors if q*d1<D]
        rough_images={q*d1 for q,d1 in pairs}
        assert len(rough_images)==len(pairs)
        assert all(math.gcd(n,math.prod([2,3,5,7,11]))==1 for n in rough_images)
        images={}
        for q,d1 in pairs:
            for e in exact_factors:
                for d0 in medium_factors:
                    d=e*d0*d1;n=q*d
                    assert math.gcd(q,d)==1
                    assert n not in images,(q,d,n)
                    images[n]=(q,d)
                    modulus_patterns+=1
        rough_count=sum(math.gcd(n,2310)==1 for n in range(1,D))
        assert len(images)==len(exact_factors)*len(medium_factors)*len(rough_images)
        assert len(images)<=len(exact_factors)*len(medium_factors)*rough_count
        rows.append({'D':D,'q_d_pairs':len(images),'distinct_large_factors':len(rough_images),
                     'exact_integer_rough_count':rough_count})
    return {'T':T,'z':z,'q_upper_endpoint':y,'modulus_patterns':modulus_patterns,'rows':rows}


if __name__=='__main__':
    ctx.prec=192
    bins=small_bins();rough=integer_roughness();aggregate=aggregate_prime_factors()
    names=[Path(__file__).name,'conditional_chen_binned_support_certificate.py',
           'conditional_chen_prefix_block_certificate.py','conditional_chen_rankin_certificate.py']
    result={'status':'PASS: independent exact subset histograms, integer roughness and aggregate prime factors',
            'binned_count_checks':bins,'integer_roughness_checks':rough,
            'aggregate_prime_factor_checks':aggregate,
            'input_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in names},
            'not_checked':['general count and density-domination lemmas','external analytic inputs','full splice']}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:{p:v for p,v in value.items() if p!='configurations'}
                     for k,value in [('bins',bins),('rough',rough)]},indent=2))
