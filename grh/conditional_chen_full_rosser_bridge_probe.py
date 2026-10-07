"""Finite all-prime Rosser interface immediately above Goldbach verification.

This does not establish the required signed AP remainder, nor a full-range
Chen theorem. It tests a different family from the two-stage direct sieve.
"""
import hashlib,json,math
from fractions import Fraction
from pathlib import Path
from flint import arb,ctx
from conditional_chen_prefix_block_certificate import certify_blocks
from conditional_chen_rankin_certificate import odd_primes_to
from conditional_chen_sharp_support_certificate import support_ratio

ROOT=Path(__file__).resolve().parent
ctx.prec=192


def q(value):
    f=Fraction(value);return arb(f.numerator)/f.denominator


def run():
    B=4000000000000000000;T=2000000;w=3
    u=Fraction(3001,1000);theta=Fraction(99,100);width=Fraction(1,200)
    finite=certify_blocks(T,w,u,width,progress=True)
    primes=odd_primes_to(T);medium=[p for p in primes if p>w]
    kappa,witness=support_ratio(medium);Q=3
    assert kappa==Fraction(5,49)
    # Check Q*kappa*T^u < B^theta by exact integer powers.
    power=math.lcm(u.denominator,theta.denominator)
    assert (Q*kappa.numerator)**power*T**(u.numerator*(power//u.denominator)) < (
        kappa.denominator**power*B**(theta.numerator*(power//theta.denominator)))
    logB=arb(B).log();log_support=arb(Q).log()+q(kappa).log()+q(u)*arb(T).log()
    support_slack=q(theta)*logB-log_support;assert support_slack>0
    V=math.prod((1-arb(1)/(p-1) for p in primes),start=arb(1))
    gap=arb(finite['lower_relative_Rosser_gap_ball']);assert 0<=gap<1
    main=V*(1-gap);assert main>0
    loss=(logB**3/arb(2).log()+logB**2)/B
    target=q('1/100000');capacity=(logB*main-loss-target)*logB**2
    C=math.floor(float(capacity))
    while not capacity>C:C-=1
    assert C>0
    samples=[]
    for N in [B,5000000000000000000,6000000000000000000,T**3]:
        L=arb(N).log()
        margin=L*main-arb(C)/L**2-(L**3/arb(2).log()+L**2)/N
        assert margin>target
        samples.append({'N':str(N),'normalized_margin_ball':str(margin)})
    names=[Path(__file__).name,'conditional_chen_prefix_block_certificate.py',
        'conditional_chen_rankin_certificate.py','conditional_chen_sharp_support_certificate.py']
    return {'status':'PASS: finite-support and positive-main-mass interface; AP hypothesis unproved',
        'precision_bits':192,'verified_goldbach_endpoint':str(B),
        'candidate_analytic_interval':[str(B),str(T**3)],'T':T,'exact_cutoff':w,
        'u':str(u),'theta':str(theta),'Q':Q,'kappa':str(kappa),'support_maximizer':witness,
        'finite_rosser_gap_certificate':finite,'complete_odd_prime_count':len(primes),
        'complete_medium_prime_count':len(medium),'complete_density_ball':str(V),
        'main_mass_lower_ball':str(main),'log_B_ball':str(logB),
        'log_actual_support_upper_ball':str(log_support),'support_slack_ball':str(support_slack),
        'signed_AP_constant_capacity_ball':str(capacity),'signed_AP_constant_budget':C,
        'claimed_count_coefficient':'1/100000','exceptional_prime_and_unit_loss_at_B_ball':str(loss),
        'interval_samples':samples,
        'required_hypothesis':'For every even N in [B,T^3], the actual deleted-prime finite Rosser sum R_W(N)=sum_d lambda_N^-(d)*(pi(N;d,N)-pi(N)/phi(d)) is >= -C*N/log(N)^4, with the stated C.',
        'continuous_interval_dependencies':'same fixed T,D0 and complete gap bound; deletion domination; prime-count lower bound pi(N)>N/log(N); increasing L*G and decreasing C/L^2 and exceptional losses for L>=log(B)',
        'prime_count_input':{'source':'Rosser-Schoenfeld (1962), Corollary 1','validity':'x>=17',
            'statement':'pi(x)>x/log(x)','doi':'https://doi.org/10.1215/ijm/1255631807'},
        'input_sha256':{n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in names},
        'not_checked':['signed AP hypothesis','independent full-formula reproduction',
            'general Rosser and deletion lemmas','external prime-count and verification sources',
            'coverage above T^3','independent research review','full finite-verification splice']}


if __name__=='__main__':
    report=run()
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['finite_rosser_gap_certificate','input_sha256']},indent=2))
