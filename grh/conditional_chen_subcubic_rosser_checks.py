"""Exact small-model checks for u>2 and scalar finite-coverage checks."""
from pathlib import Path
from fractions import Fraction
from bisect import bisect_right
import math
import json
import hashlib
from flint import arb, ctx
from conditional_chen_subcubic_rosser_certificate import certify
from conditional_chen_prefix_block_checks import exact_gaps, exact_relaxed_sequence

ROOT = Path(__file__).resolve().parent
ctx.prec = 192


def q(value):
    value = Fraction(value)
    return arb(value.numerator)/value.denominator


def close(a,b):
    assert abs(a-b)<arb('1e-43'),(str(a),str(b))


def small_primes(T,w):
    return [p for p in range(max(3,w+1),T+1) if all(p%d for d in range(2,math.isqrt(p)+1))]


def independent_grid(T,ps,u,h):
    logD=q(u)*arb(T).log()
    product,K=1,0
    for p in ps:
        product*=p
        if product**u.denominator<T**u.numerator:K+=1
        else:break
    J=int((logD/q(h)).ceil().unique_fmpz())
    groups={}
    for p in sorted(ps,reverse=True):
        label=int((arb(p).log()/q(h)).floor().unique_fmpz())
        groups.setdefault(label,[]).append(p)
    def cutoff(p):
        if T**u.numerator==p**(3*u.denominator):return 0
        return int(((logD-3*arb(p).log())/q(h)).ceil().unique_fmpz())
    blocks=[]
    for label,group in groups.items():
        blocks.append({'label':label,'primes':group,
            'passing_source_cutoff':max(0,cutoff(min(group))),
            'failure_source_cutoff':max(0,cutoff(max(group))-K+1)})
    return blocks,K,J


def exact_prefix_support(T,ps,u):
    kappa=max([Fraction(1,ps[0]**2)]+[Fraction(r,s*s) for r,s in zip(ps,ps[1:])])
    accepted=failures=singleton_outside_sharp=0
    for mask in range(1<<len(ps)):
        chosen=sorted([p for i,p in enumerate(ps) if mask>>i&1],reverse=True)
        for parity in [0,1]:
            product=1;alive=True
            for m,p in enumerate(chosen,1):
                product*=p
                if m%2==parity and (product*p*p)**u.denominator>=T**u.numerator:
                    assert product**u.denominator<T**u.numerator
                    failures+=1;alive=False;break
            if alive:
                accepted+=1
                if len(chosen)<=1:
                    assert product<=T
                    if (product*kappa.denominator)**u.denominator>=kappa.numerator**u.denominator*T**u.numerator:
                        singleton_outside_sharp+=1
                else:
                    assert (product*kappa.denominator)**u.denominator<kappa.numerator**u.denominator*T**u.numerator
    return accepted,failures,singleton_outside_sharp


def finite_models():
    configs=[(29,w,u,h) for w in [2,3,5]
             for u in map(Fraction,['201/100','11/5','5/2','297/100','3','3001/1000'])
             for h in map(Fraction,['1/2','1/10'])]
    configs += [(T,w,u,Fraction(1,10)) for T,w in [(3,2),(7,3),(7,5)]
                for u in map(Fraction,['201/100','5/2','3'])]
    subsets=accepted=failures=deleted=patterns=fresh_count=zero_blocks=outside=0
    witnesses=[]
    for T,w,u,h in configs:
        ps=small_primes(T,w)
        blocks,K,J=independent_grid(T,ps,u,h)
        result=certify(T,w,u,h)
        assert result['maximum_first_failure_prefix_degree']==K and result['label_limit']==J
        zero_blocks+=sum(b['passing_source_cutoff']==0 for b in blocks)
        reference=exact_relaxed_sequence(blocks,J,K)
        for value,key in zip(reference,['upper_relative_Rosser_gap_ball','lower_relative_Rosser_gap_ball']):
            close(arb(result[key]),q(value))
        exact=exact_gaps(ps,T,u)
        assert all(0<=a<=b for a,b in zip(exact,reference))
        a,f,s=exact_prefix_support(T,ps,u)
        subsets+=1<<len(ps);accepted+=a;failures+=f;outside+=s
        if h==Fraction(1,10) and (T<29 or u in [Fraction(11,5),Fraction(297,100),Fraction(3)]):
            for available in range(1<<len(ps)):
                kept=[p for i,p in enumerate(ps) if available>>i&1]
                true=exact_gaps(kept,T,u)
                assert all(0<=a<=b for a,b in zip(true,reference))
                deleted+=1;patterns+=1<<len(kept)
                if u==Fraction(297,100) or T<29:
                    fresh=certify(T,w,u,h,kept)
                    deleted_blocks=[{**block,'primes':[p for p in block['primes'] if p in kept]} for block in blocks]
                    exact_deleted=exact_relaxed_sequence(deleted_blocks,J,K)
                    for value,key,a,b in zip(exact_deleted,['upper_relative_Rosser_gap_ball','lower_relative_Rosser_gap_ball'],true,reference):
                        close(arb(fresh[key]),q(value));assert a<=value<=b
                    fresh_count+=1
        witnesses.append({'T':T,'w':w,'u':str(u),'h':str(h),
                          'kappa_D_dominates_T':result['kappa_D_dominates_T'],
                          'zero_passing_blocks':result['zero_passing_cutoff_blocks']})
    assert zero_blocks>0 and outside>0
    return {'configurations':witnesses,'enumerated_complete_subsets':subsets,
            'accepted_support_cases':accepted,'first_failure_prefix_cases':failures,
            'singletons_outside_kappa_D_cases':outside,'zero_passing_block_cases':zero_blocks,
            'deleted_prime_sets':deleted,'pointwise_divisor_patterns':patterns,
            'fresh_deleted_gap_certificates_checked':fresh_count}


def independent_primes(n):
    flags=bytearray([1])*((n+1)//2);flags[0]=0
    for p in range(3,math.isqrt(n)+1,2):
        if flags[p//2]:
            for i in range(p*p//2,len(flags),p):flags[i]=0
    return [2*i+1 for i,flag in enumerate(flags) if flag]


def verify_input_hashes(record):
    for n,h in record['input_sha256'].items():
        assert hashlib.sha256((ROOT/n).read_bytes()).hexdigest()==h,n


def coverage(report):
    verify_input_hashes(report)
    assert report['conditional_analytic_coverage']==['4000000000000000000','1728000000000000000000']
    assert report['common_signed_AP_constant_budget']==3489
    assert 3**100*2**297<25**100
    all_primes=independent_primes(12000000)
    C=report['common_signed_AP_constant_budget'];theta=Fraction(report['theta'])
    rows=[];previous=None
    for row in report['rows']:
        finite=row['finite_certificate'];verify_input_hashes(finite)
        T=finite['T'];u=Fraction(finite['u']);left=int(row['left_endpoint']);right=int(row['right_endpoint'])
        assert right==T**3 and (previous is None or left==previous)
        previous=right
        ps=all_primes[:bisect_right(all_primes,T)]
        assert len(ps)==finite['complete_odd_prime_count']
        V=math.prod((1-arb(1)/(p-1) for p in ps),start=arb(1))
        close(V,arb(row['complete_density_ball']))
        G=V*(1-arb(finite['lower_relative_Rosser_gap_ball']));close(G,arb(row['main_mass_lower_ball']))
        assert G>0
        kappa=max([Fraction(1,ps[0]**2)]+[Fraction(r,s*s) for r,s in zip(ps,ps[1:])])
        assert kappa==Fraction(3,25)
        power=math.lcm(u.denominator,theta.denominator)
        assert 3**power*T**(u.numerator*(power//u.denominator))<25**power*left**(theta.numerator*(power//theta.denominator))
        product,K=1,0
        for p in ps:
            product*=p
            if product**u.denominator<T**u.numerator:K+=1
            else:break
        assert K==finite['maximum_first_failure_prefix_degree']
        J=int((q(u)*arb(T).log()/q(Fraction(finite['h']))).ceil().unique_fmpz())
        assert J==finite['label_limit']
        L=arb(left).log()
        capacity=(L*G-(L**3/arb(2).log()+L**2)/left-arb(1)/100000)*L**2
        close(capacity,arb(row['signed_C_capacity_ball']))
        assert capacity>row['signed_C_integer_budget'] and capacity<row['signed_C_integer_budget']+1
        samples=[]
        for N in [left,(left+right)//2,right]:
            ln=arb(N).log()
            margin=ln*G-C/ln**2-(ln**3/arb(2).log()+ln**2)/N
            assert margin>arb(1)/100000
            if N==left:close(margin,arb(row['common_C_normalized_margin_at_left_ball']))
            samples.append({'N':str(N),'normalized_margin_ball':str(margin)})
        rows.append({'T':T,'u':str(u),'complete_odd_prime_count':len(ps),
                     'capacity_integer':row['signed_C_integer_budget'],'samples':samples})
    old=json.loads((ROOT/'conditional_chen_full_rosser_w2_probe.json').read_text(encoding='utf-8'))
    verify_input_hashes(old)
    first=report['rows'][0]['finite_certificate']
    for key in ['upper_relative_Rosser_gap_ball','lower_relative_Rosser_gap_ball']:
        close(arb(first[key]),arb(old['rows'][1]['gap'][key]))
    return {'complete_odd_primes_generated':len(all_primes),'interval_rows':rows,
            'continuous_proof_dependency':'fixed weights and monotonic normalized bound on each closed interval; shared endpoints',
            'general_doubling_support_test':'3^100*2^297<25^100',
            'old_u_61_20_finite_gap_reproduced':True}


def run():
    report=json.loads((ROOT/'conditional_chen_subcubic_rosser_certificate.json').read_text(encoding='utf-8'))
    large=coverage(report)
    small=finite_models()
    names=[Path(__file__).name,'conditional_chen_subcubic_rosser_certificate.py',
           'conditional_chen_subcubic_rosser_certificate.json','conditional_chen_prefix_block_checks.py',
           'conditional_chen_full_rosser_w2_probe.json']
    return {'status':'PASS: exact u>2 models, deleted-prime domination and scalar coverage checks',
            'precision_bits':192,'conditional_coverage_checks':large,'finite_models':small,
            'input_sha256':{n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names},
            'not_checked':['signed AP input','independent general proof review','unbounded-T main-mass bound','full splice']}


if __name__=='__main__':
    result=run()
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':result['status'],'large_prime_count':result['conditional_coverage_checks']['complete_odd_primes_generated'],
                     'small':{k:v for k,v in result['finite_models'].items() if k!='configurations'}},indent=2))
