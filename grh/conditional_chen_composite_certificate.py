"""Arb arithmetic for both composite presieves in the manuscript.

The analytic weight inequalities are proved in the manuscript. This program
certifies the scalar inputs; it is not a mechanical proof of the manuscript
or of the cited external analytic lemmas. No fixed-endpoint-only product
condition is used by the new GRH threshold.
"""
import json
import math
from pathlib import Path
from fractions import Fraction
from flint import arb, acb, ctx
from conditional_chen_rankin_certificate import certify_rankin
from conditional_chen_uniform_product_certificate import certify_uniform_products
from conditional_chen_prefix_band_certificate import certify_band
from conditional_chen_uniform_grh_probe import surplus as propose_surplus

ctx.prec = 192
gamma = arb.const_euler()
G = gamma.exp()

def primes_to(n):
    sieve = bytearray(b'\x01') * (n + 1)
    sieve[:2] = b'\x00\x00'
    for p in range(2, math.isqrt(n) + 1):
        if sieve[p]:
            sieve[p*p:n+1:p] = b'\x00' * ((n-p*p)//p+1)
    return [p for p in range(2, n+1) if sieve[p]]

ps = primes_to(10_000)
ct = arb(1)
for p in ps[1:]:
    ct *= 1 - arb(1)/(p-1)**2
c_infty_lower = ct*(-arb(1)/9999).exp()
assert c_infty_lower > arb('.66')
Umin = arb('1.32')*(-gamma).exp()

# M = gamma + sum_p (log(1-1/p)+1/p); omitted terms are negative.
m_upper = gamma + sum(((1-arb(1)/p).log()+arb(1)/p
                      for p in ps if p <= 1000), arb(0))
assert m_upper < arb('.262')
g_correction = sum((arb(1)/(p*(p-1))
                    for p in ps[1:] if p <= 1000), arb(0))+arb(1)/1000

def f(s):
    return 2*G*(s-1).log()/s

def F_3_to_4(s):
    assert 3 <= s <= 4
    correction = acb.integral(
        lambda x, analytic: (x-2).log(analytic=analytic)/(x-1),
        acb(3), acb(s)).real
    return 2*G*(1+correction)/s

def h(s):
    return 3*(-s).exp()/s if s > 3 else ((-s).exp() if s > 2 else (-arb(2)).exp())

def interface(T_integer, depth, exact_cutoff=2):
    T = arb(T_integer)
    assert 10_000 <= T_integer <= 10**12
    assert depth > 0 and depth % 2 == 0
    assert 2 <= exact_cutoff <= 10_000 and exact_cutoff < T_integer
    exact_primes = [p for p in ps[1:] if p <= exact_cutoff]
    Q = math.prod(exact_primes)
    VQ = arb(1)
    WQ = arb(0)
    for p in exact_primes:
        VQ *= 1-arb(1)/(p-1)
        WQ += arb(1)/(p-1)
    # BJS Lemma 16 up to 10^12; JS Lemma 3.2 above 10^12.
    b = 2/(T.sqrt()*T.log())
    beyond = arb('.000006836')/arb(10**12).log()
    assert b > beyond
    r = arb('.000002964')/T.log()+1/(T-1)+arb('.54')/T
    epsilon = (b+r).exp()-1
    W = T.log().log()+arb('.262')-arb(1)/2+g_correction+b-WQ
    assert W > 0
    # The medium-prime product is at least the complete odd-prime product
    # through T divided by VQ. Removing primes dividing N increases it.
    eta = VQ*W**depth/arb(math.factorial(depth))*T.log()*b.exp()/Umin
    support = arb(Q).log()+depth*T.log()
    return {
        'T':T_integer, 'depth':depth, 'uniform_epsilon_ball':str(epsilon),
        'exact_cutoff':exact_cutoff,'exact_presieve_Q_max':Q,
        'exact_product_ball':str(VQ),'exact_density_sum_ball':str(WQ),
        'reciprocal_upper_error_ball':str(b), 'W_ball':str(W),
        'relative_Brun_gap_ball':str(eta),
        'log_support_cost_ball':str(support)
    }, epsilon, eta, support

def qG_upper(w):
    return (arb('.165')+arb('12.683')/w+arb('254.980')/w**2
      +arb('2607.854')/w**3+arb('11605.056')/w**4
      +(arb('1.314')*w+arb('.092')*w.log()+arb('60.883')
        +arb('8.250')*w.log()/w+arb('939.260')/w)*(-w/4).exp())

uniform = certify_uniform_products()
uniform_rows = {r['T']:r for r in uniform['rows']}
GRH_SIGMAS = [Fraction(0),Fraction(1,10),Fraction(7,50),Fraction(9,50),
              Fraction(11,50),Fraction(7,25),Fraction(9,25)]
EH_SIGMAS = [Fraction(0)]+[Fraction(a,20) for a in range(2,13)]
GRH_BANDS = [Fraction(a,25) for a in range(2,7)]
EH_BANDS = [Fraction(0),Fraction(3,25),Fraction(6,25),Fraction(9,25)]


def band_presieve(T,w,u,sigmas,band_sigmas,M,progress):
    value = certify_band(T,w,u,sigmas,band_sigmas,M,progress)
    row = uniform_rows[T]
    value['uniform_epsilon_ball'] = row['ordered_uniform_epsilon_ball']
    value['coarse_uniform_epsilon_ball'] = row['previous_uniform_epsilon_ball']
    value['uniform_product_input'] = 'ordered finite certificate supplied by caller'
    return value


pre = band_presieve(5_000_000,5,Fraction(49,8),
                   GRH_SIGMAS,GRH_BANDS,12,True)
ep_bound = arb(pre['uniform_epsilon_ball'])
eta_bound = arb(pre['relative_Rosser_gap_ball'])
support = arb(pre['log_support_cost_ball'])
epsilon = arb(1)/41200
def rational_gap_inverse(value):
    proposed = math.floor(1/float(value))
    while not value < arb(1)/proposed:
        proposed -= 1
    return proposed


eta_inverse = rational_gap_inverse(arb(pre['upper_relative_Rosser_gap_ball']))
eta_minus_inverse = rational_gap_inverse(arb(pre['lower_relative_Rosser_gap_ball']))
eta = arb(1)/eta_inverse
eta_minus = arb(1)/eta_minus_inverse
C1, C2 = 106, 108
assert ep_bound < epsilon < arb(1)/10000
assert arb(pre['upper_relative_Rosser_gap_ball']) < eta
assert arb(pre['lower_relative_Rosser_gap_ball']) < eta_minus
assert eta_bound < eta+eta_minus
# Floating arithmetic proposes a fixed endpoint; all claims and conditions
# below are then tested by Arb using freshly generated proof inputs.
for proposed_L in [31550,31600,31650,31700,31750,31800,31850]:
    proposed_margin,proposed_alpha = propose_surplus(
        proposed_L,float(support),float(epsilon),float(eta),float(eta_minus))
    if proposed_margin > .0007:
        log_N0 = proposed_L
        alpha_numerator = round(proposed_alpha*10**6)
        break
else:
    raise AssertionError("No proposed lower GRH endpoint has sufficient margin")
L = arb(log_N0)
alpha = arb(alpha_numerator)/10**6
beta = arb(7)/2
t = arb(1)/2-alpha
s = 8*t
k = 8*(t-arb(1)/3)
delta = arb(1)/10**8
ell = (1+delta).log()
v = arb(1)/10**10
assert alpha*L-beta*L.log() > support
assert alpha-beta/L > 0
assert L/2-beta*L.log() > arb(10**9).log()
assert L/8 > arb(5_000_000).log()
assert s > 3 and s < 4 and k > 1
assert L*(-L/16).exp() < v
assert 3*(L/8)/(8*arb.pi())*(-L/16).exp() < 1/L**2

E20 = arb(20).exp()
li20 = arb(20).ei()
assert E20/20*(1+arb(1)/20+arb(3)/400) > li20
assert li20 > E20/20*(1+arb(1)/20)
assert 3/L**2+L**2/(8*arb.pi())*(-L/2).exp()+2*L**2/(arb(2).log()*L.exp()) < 1/L
assert L**2/(72*arb.pi())*(-L/6).exp() < 1/L-27/L**2

assert qG_upper(L) < arb('.17')
assert qG_upper(arb(8500)) < arb('.166496')
assert qG_upper(L/8) < arb('.64')
assert arb('.165')-arb('237.934')/(L/8)*(-L/16).exp() > arb('.16')
per_modulus = arb('.17')+1/(16*arb.pi())+2*(-L/2).exp()/arb(2).log()
assert arb('.608')*per_modulus < arb('.12')
assert per_modulus < arb('.2')

cb = acb.integral(lambda x,analytic:(2-3*x).log(analytic=analytic)/(x*(1-x)),
                  acb(arb(1)/8),acb(arb(1)/3)).real
assert cb < arb('.363084')
Fs = F_3_to_4(s)
Flo = Fs+C1*epsilon*arb(2).exp()*h(s)
large_lower = f(s)-C2*epsilon*arb(2).exp()*h(s)
assert 0 < large_lower <= Flo
assert f(s) <= Fs < 2*G/3
composite_lower = (1+eta)*large_lower-(eta+eta_minus)*Flo
assert composite_lower > 0
low = 8*(1-v)*composite_lower
I = G/(4*t)*(6*(3-8*alpha)/(3-18*alpha)).log()
FM = 2*G/k+C1*epsilon
J = (arb(8)/3).log()
b0 = t-arb(1)/4
assert arb(1)/8 < b0 < arb(1)/3
Hh = acb.integral(
    lambda b, analytic: (2-8*(acb(t)-b)).exp()/b,
    acb(arb(1)/8), acb(b0)).real + ((arb(1)/3)/b0).log()
assert 0 < Hh < J
upper = 4*(1+eta)*(1+2/L)*(1+v)/(1-(-L/8).exp())*(I+C1*epsilon*Hh+2*FM/L**2)
switch = (1+eta)*(1+v)/2*(2*G/t+3*C1*epsilon)*(1+delta)*(1+4/L)*(cb+2/L**2+(J+2/L**2)*(6*ell/L+8/L**2))
ap = arb('.18')/(Umin*L.sqrt())+arb('.055')*L**4*(-L/6).exp()/Umin+arb('.8')*(1+eta)*(1+v)*FM*L**2*(-L/6).exp()
rect = (1+delta)*L**5/(2*Umin*ell)*(arb('.75')*(-L/48).exp()+(-L/8).exp())
finite = L**2/Umin*(2*(-L/8).exp()+(-2*L/3).exp()+(-L).exp())
margin = low-upper-switch-ap-rect-finite
assert margin > arb('.0006')
assert rect < arb('1e-254')
assert finite < arb('1e-1700')
assert arb('2.2')/(4*arb(2).log()) < 1
mb = arb('2.81')*arb('.64')**(arb(1)/3)
mb += arb('2.809')/(arb('.16')**(arb(2)/3)*(L/8))*(-L/16).exp()
mb += 13*(-L/12).exp()+9*((-5*L/18).exp()+(-L/24).exp())+31*(-13*L/48).exp()
assert mb < 3
assert 13*(1+delta)**(arb(1)/3)*(-L/12).exp() < 1

# Distribution-only reductions. Neither RH nor GRH is used here.
def eh_case(theta, log_start, inverse_error_budget, claimed_margin,
            T_integer,u,exact_cutoff,
            epsilon_inverse,eta_plus_inverse,eta_minus_inverse,C1,C2):
    eh_pre = band_presieve(T_integer,exact_cutoff,u,
                          EH_SIGMAS,EH_BANDS,10,False)
    eh_eps_bound = arb(eh_pre['uniform_epsilon_ball'])
    eh_eta_bound = arb(eh_pre['relative_Rosser_gap_ball'])
    eh_support = arb(eh_pre['log_support_cost_ball'])
    eh_epsilon = arb(1)/epsilon_inverse
    eh_eta_plus = arb(1)/eta_plus_inverse
    eh_eta_minus = arb(1)/eta_minus_inverse
    assert eh_eps_bound < eh_epsilon
    assert arb(eh_pre['upper_relative_Rosser_gap_ball']) < eh_eta_plus
    assert arb(eh_pre['lower_relative_Rosser_gap_ball']) < eh_eta_minus
    assert eh_eta_bound < eh_eta_plus+eh_eta_minus
    assert ((epsilon_inverse >= 2000 and (C1,C2) == (109,110))
            or (epsilon_inverse >= 1000 and (C1,C2) == (113,114))
            or (epsilon_inverse >= 700 and (C1,C2) == (116,117)))
    LE = arb(log_start)
    se = 3*theta-3*eh_support/LE
    assert se > 2 and se < 3
    assert LE/3 > arb(10**12).log()
    assert LE/3 > arb(T_integer).log()
    assert theta*LE > arb(10**9).log()
    rho = 1-arb('.000020508')/LE-2*LE*(-LE/3).exp()
    factor = f(se)-C2*eh_epsilon*(-(se-2)).exp()
    large_upper = 2*G/se+C1*eh_epsilon*(-(se-2)).exp()
    assert 0 < factor <= large_upper
    factor = (1+eh_eta_plus)*factor-(eh_eta_plus+eh_eta_minus)*large_upper
    assert factor > 0
    main = 3*Umin*rho*factor
    deleted = LE**3/arb(2).log()*(2*(-LE/2).exp()+arb('1.1')*theta*LE*(-LE).exp())
    surplus = main-arb(1)/inverse_error_budget-deleted-LE**2*(-LE).exp()
    assert surplus > arb(claimed_margin)
    return {
      'theta':str(theta), 'log_start':log_start,
      'presieve':eh_pre,'epsilon':f'1/{epsilon_inverse}',
      'eta_plus':f'1/{eta_plus_inverse}','eta_minus':f'1/{eta_minus_inverse}',
      'C1':C1,'C2':C2,
      'input':'E_theta(x) <= C x/log(x)^4 for all x>=X; numerical C,X are not supplied',
      'threshold_formula':f'max(X, exp(max({log_start}, sqrt({inverse_error_budget} C))))',
      'assumptions':'the stated finite distribution input only; no RH or GRH',
      's_at_start_ball':str(se),'sieve_factor_ball':str(factor),
      'main_ball':str(main),'surplus_ball':str(surplus),
      'claimed_margin':claimed_margin,'margin_units':'N/log(N)^2',
      'not_an_effective_EH_only_numeric_record':True
    }

eh75 = eh_case(arb(3)/4,823,20,'.01',22_000,Fraction(39,8),3,1717,141,143,113,114)
eh90 = eh_case(arb(9)/10,225,10,'.02',10_000,Fraction(35,8),3,1012,46,45,113,114)
eh95 = eh_case(arb(19)/20,176,10,'.02',10_000,Fraction(17,4),3,1012,34,33,113,114)
eh99 = eh_case(arb(99)/100,150,10,'.02',10_000,Fraction(17,4),3,1012,34,33,113,114)

result = {
 'status':'strict scalar arithmetic for the composite-sieve manuscript; external analytic lemmas remain cited dependencies',
 'precision_bits':ctx.prec, 'C_infinity_lower_ball':str(c_infty_lower),
 'U_min_ball':str(Umin),'Meissel_Mertens_upper_ball':str(m_upper),
 'small_g_correction_upper_ball':str(g_correction),
 'GRH':{'log_N0':log_N0,'log_log_N0_ball':str(L.log()),
   'alpha':f'{alpha_numerator}/1000000','beta':'7/2','epsilon':'1/41200',
   'eta_plus':f'1/{eta_inverse}','eta_minus':f'1/{eta_minus_inverse}',
   'delta':'1/100000000','C1':C1,'C2':C2,'presieve':pre,
   'support_slack_ball':str(alpha*L-beta*L.log()-support),
   'F_s_ball':str(Fs),'h_integral_ball':str(Hh),
   'lower_composite_factor_ball':str(composite_lower),
   'retained_sieve_functions':'F(s)=2 exp(gamma)/s * (1+integral_3^s log(x-2)/(x-1) dx); the h(s_q) correction is integrated piecewise with the original weighted-sum endpoint charge',
   'claimed_margin':'.0006','margin_units':'U_N N/log(N)^2',
   'lower_ball':str(low),'half_upper_ball':str(upper),'half_switch_ball':str(switch),
   'AP_budget_ball':str(ap),'half_rectangle_error_ball':str(rect),
   'finite_error_ball':str(finite),'margin_ball':str(margin),
   'sieve_interface':'uniform at every upper endpoint; exact/Rosser/Rosser composite weights have coefficients of absolute value <=1 and support <Q D0 D; the medium gap is certified by a complete finite Rankin sum'},
 'EH_3_4':eh75,'EH_9_10':eh90,'EH_19_20':eh95,'EH_99_100':eh99,
 'uniform_product_certificate':{
   'S':uniform['S'],'M_cutoff':uniform['M_cutoff'],
   'complete_odd_prime_count':uniform['complete_odd_prime_count'],
   'checked_prime_jumps':uniform['checked_prime_jumps'],
   'regenerated_in_this_run':True,
   'certificate_file':'conditional_chen_uniform_product_certificate.json'},
 'historical_claims':'exp(8500), exp(15000) remain unaudited fixed-endpoint candidates; finite QEH_B bridge is refuted'
}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
Path(__file__).with_name('conditional_chen_prefix_band_certificate.json').write_text(
    json.dumps(pre,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
Path(__file__).with_name('conditional_chen_uniform_product_certificate.json').write_text(
    json.dumps(uniform,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result,ensure_ascii=False,indent=2))
