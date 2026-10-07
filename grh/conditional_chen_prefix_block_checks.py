"""Independent exact subset and sequential-prime checks for prefix blocks."""
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

from flint import arb,ctx

from conditional_chen_prefix_block_certificate import grid_blocks,certify_blocks,block_coefficients


def exact_gaps(atoms,T,u):
    target=T**u.numerator
    V=math.prod((Fraction(p-2,p-1) for p in atoms),start=Fraction(1))
    masses=[Fraction(),Fraction()]
    failures=[Fraction(),Fraction()]
    weights=[[],[]]
    for mask in range(1<<len(atoms)):
        chosen=sorted((p for i,p in enumerate(atoms) if mask>>i&1),reverse=True)
        density=math.prod((Fraction(1,p-1) for p in chosen),start=Fraction(1))
        mu=(-1)**len(chosen)
        for sign,parity in enumerate([1,0]):
            passed=True;prefix=1
            for m,p in enumerate(chosen,1):
                prefix*=p
                if m%2==parity and (prefix*p*p)**u.denominator>=target:
                    passed=False
                    break
            weights[sign].append(mu if passed else 0)
            masses[sign]+=mu*density*int(passed)
            if chosen and len(chosen)%2==parity:
                prefix=1;previous=True
                for m,p in enumerate(chosen,1):
                    prefix*=p
                    if m<len(chosen) and m%2==parity:
                        previous &= (prefix*p*p)**u.denominator<target
                r=chosen[-1]
                if previous and (prefix*r*r)**u.denominator>=target:
                    below=math.prod((Fraction(p-2,p-1) for p in atoms if p<r),start=Fraction(1))
                    failures[sign]+=density*below
    assert masses[0]-V==failures[0] and V-masses[1]==failures[1]
    for bit in range(len(atoms)):
        for mask in range(1<<len(atoms)):
            if mask>>bit&1:
                for values in weights: values[mask]+=values[mask^(1<<bit)]
    for mask,(up,lo) in enumerate(zip(*weights)):
        assert lo<=int(mask==0)<=up
    return [(masses[0]-V)/V,(V-masses[1])/V]


def exact_relaxed_sequence(blocks,limit,K):
    """One prime at a time; exact rational sparse vectors, no block formula."""
    states=[{0:Fraction(1)},{},{0:Fraction(1)},{}]
    totals=[Fraction(),Fraction()]
    for block in blocks:
        b=block['label'];passing=block['passing_source_cutoff'];failure=block['failure_source_cutoff']
        degrees=[states]+[[{} for _ in range(4)] for _ in range(K)]
        for p in block['primes']:
            t=Fraction(1,p-2)
            for values in degrees:
                totals[0]+=t*sum((v for i,v in values[0].items() if i>=failure),Fraction())
                totals[1]+=t*sum((v for i,v in values[3].items() if i>=failure),Fraction())
            updated=[[{i:(1+t)*v for i,v in s.items()} for s in values] for values in degrees]
            for j in range(1,K+1):
                for dest,source,checked in [(0,1,False),(1,0,True),(2,3,True),(3,2,False)]:
                    for i,v in degrees[j-1][source].items():
                        if i+b<limit and (not checked or i<passing):
                            updated[j][dest][i+b]=updated[j][dest].get(i+b,Fraction())+t*v
            degrees=updated
        states=[{} for _ in range(4)]
        for values in degrees:
            for i,s in enumerate(values):
                for label,value in s.items():states[i][label]=states[i].get(label,Fraction())+value
    return totals


def check_coefficients(primes,K):
    elementary=[Fraction(1)]+[Fraction() for _ in primes]
    for p in primes:
        t=Fraction(1,p-2)
        for k in range(len(primes),0,-1):elementary[k]+=t*elementary[k-1]
    d,c=block_coefficients(primes,K)
    for j in range(K+1):
        exact_d=sum((math.comb(k,j)*elementary[k] for k in range(j,len(primes)+1)),Fraction())
        exact_c=sum((math.comb(k-1,j)*elementary[k] for k in range(j+1,len(primes)+1)),Fraction())
        assert abs(d[j]-arb(exact_d.numerator)/exact_d.denominator)<arb('1e-45')
        assert abs(c[j]-arb(exact_c.numerator)/exact_c.denominator)<arb('1e-45')


def run():
    ctx.prec=192
    configs=[];deleted=identities=patterns=coeff_groups=0
    for w in [3,5]:
        for u in [Fraction(7,2),Fraction(4),Fraction(49,8)]:
            for h in [Fraction(1),Fraction(1,2),Fraction(1,10),Fraction(1,100)]:
                T=29
                _,atoms,blocks,K,limit=grid_blocks(T,w,u,h)
                production=certify_blocks(T,w,u,h,False)
                reference=exact_relaxed_sequence(blocks,limit,K)
                for block in blocks:
                    check_coefficients(block['primes'],K);coeff_groups+=1
                for a,key in zip(reference,['upper_relative_Rosser_gap_ball','lower_relative_Rosser_gap_ball']):
                    assert abs(arb(production[key])-arb(a.numerator)/a.denominator)<arb('1e-42')
                for mask in range(1<<len(atoms)):
                    subset=[p for i,p in enumerate(atoms) if mask>>i&1]
                    exact=exact_gaps(subset,T,u)
                    assert all(0<=a<=b for a,b in zip(exact,reference))
                    deleted+=1;identities+=2;patterns+=1<<len(subset)
                configs.append({'T':T,'w':w,'u':str(u),'h':str(h),'K':K,
                                'complete_upper_bound':str(reference[0]),
                                'complete_lower_bound':str(reference[1])})
    root=Path(__file__).resolve().parent
    names=[Path(__file__).name,'conditional_chen_prefix_block_certificate.py',
           'conditional_chen_rankin_certificate.py']
    return {'status':'PASS: independent exact subset and sequential-prime comparisons',
            'precision_bits':ctx.prec,'configurations':configs,
            'deleted_prime_subsets':deleted,'exact_first_failure_identities':identities,
            'pointwise_divisor_patterns':patterns,'exact_block_coefficient_groups':coeff_groups,
            'input_sha256':{p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in names},
            'not_checked':['general domination proof','external analytic inputs','full Chen splice']}


if __name__=='__main__':
    r=run()
    Path(__file__).with_suffix('.json').write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:v for k,v in r.items() if k!='configurations'},indent=2))
