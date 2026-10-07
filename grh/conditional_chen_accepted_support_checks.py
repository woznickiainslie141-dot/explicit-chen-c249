"""Exact subset and deleted-set references for six-state support counts."""
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
from flint import arb,ctx
from conditional_chen_rankin_certificate import odd_primes_to
from conditional_chen_accepted_support_certificate import certify,STATES

ROOT=Path(__file__).resolve().parent
ctx.prec=192


def q(f):
    f=Fraction(f);return arb(f.numerator)/f.denominator


def reference(T,ps,u,h,kappa,K,J):
    labels={p:int((arb(p).log()/q(h)).floor().unique_fmpz()) for p in ps}
    block_min={label:min(p for p in ps if labels[p]==label) for label in set(labels.values())}
    cutoffs={label:int(((q(u)*arb(T).log()-3*arb(p).log())/q(h)).ceil().unique_fmpz())
             for label,p in block_min.items()}
    records=[]
    for mask in range(1<<len(ps)):
        chosen=sorted([p for i,p in enumerate(ps) if mask>>i&1],reverse=True)
        prefix=1;true_flags=relaxed=3;B=0;degrees={}
        for m,p in enumerate(chosen,1):
            checked=1 if m%2 else 2
            prefix*=p
            if (prefix*p*p)**u.denominator>=T**u.numerator:true_flags&=3^checked
            if B>=cutoffs[labels[p]]:relaxed&=3^checked
            B+=labels[p];degrees[labels[p]]=degrees.get(labels[p],0)+1
        if true_flags:
            assert prefix**u.denominator*kappa.denominator**u.denominator<T**u.numerator*kappa.numerator**u.denominator
            assert len(chosen)<=K and B<J and relaxed
        keep=bool(relaxed) and B<J and all(d<=K for d in degrees.values())
        records.append({'mask':mask,'actual':bool(true_flags),'label':B,
                        'state':(relaxed,len(chosen)%2),'keep':keep})
    return records


def histogram(records,available):
    hist={state:{} for state in STATES};actual=0
    for record in records:
        if record['mask']&~available:continue
        actual+=record['actual']
        if record['keep']:
            state=hist[record['state']];B=record['label'];state[B]=state.get(B,0)+1
    return hist,actual


def compare(production,hist):
    checked=0
    for state,p in zip(STATES,production):
        r=hist[state]
        for B in range(max(len(p),max(r,default=0)+1)):
            assert abs(p[B]-r.get(B,0))<arb('1e-43'),(state,B)
            checked+=1
    return checked


def run():
    configurations=[];subsets=coefchecks=deletionpatterns=deleted_productions=0
    for w in [3,5]:
        ps=[p for p in odd_primes_to(29) if p>w]
        kappa=Fraction(5,49) if w==3 else Fraction(11,169)
        for u in [Fraction(7,2),Fraction(4),Fraction(49,8)]:
            for h in [Fraction(9,10),Fraction(1,2),Fraction(1,10),Fraction(1,100)]:
                result,production=certify(29,w,u,h,return_states=True)
                assert Fraction(result['kappa'])==kappa
                K,J=result['maximum_actual_degree'],result['label_limit']
                records=reference(29,ps,u,h,kappa,K,J);subsets+=len(records)
                full,actual=histogram(records,(1<<len(ps))-1)
                coefchecks+=compare(production,full)
                assert sum(sum(r.values()) for r in full.values())>=actual
                for available in range(1<<len(ps)):
                    hist,count=histogram(records,available)
                    assert count<=arb(result['count_upper_ball'])
                    deletionpatterns+=sum(not record['mask']&~available for record in records)
                    if u==4 and h==Fraction(1,10):
                        available_primes=[p for i,p in enumerate(ps) if available>>i&1]
                        deleted,polys=certify(29,w,u,h,available_primes,return_states=True)
                        coefchecks+=compare(polys,hist);deleted_productions+=1
                        assert count<=arb(deleted['count_upper_ball'])<=arb(result['count_upper_ball'])
                        for p,full_p in zip(polys,production):
                            for B in range(len(p)):assert p[B]<=full_p[B]
                configurations.append({'w':w,'u':str(u),'h':str(h),'exact_union_count':actual,
                    'relaxed_union_count':str(result['count_upper_ball'])})
    names=[Path(__file__).name,'conditional_chen_accepted_support_certificate.py',
           'conditional_chen_rankin_certificate.py','conditional_chen_sharp_support_certificate.py',
           'conditional_chen_prefix_block_certificate.py','conditional_chen_binned_support_certificate.py']
    return {'status':'PASS: exact union histograms, actual acceptance and deleted-prime domination',
        'configurations':configurations,'enumerated_complete_subsets':subsets,
        'integer_histogram_coefficients_checked':coefchecks,'deleted_subset_patterns':deletionpatterns,
        'deleted_production_histograms_checked':deleted_productions,
        'input_sha256':{n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names},
        'not_checked':['general six-state domination lemma','external analytic inputs','Chen endpoints','full splice']}


if __name__=='__main__':
    result=run()
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['configurations','input_sha256']},indent=2))
