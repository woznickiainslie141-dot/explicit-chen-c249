"""Numerical certificates for the accompanying conditional Chen manuscript.

This certifies arithmetic sufficiency, not a GRH-only theorem. The small-Q
application of the explicit linear sieve has an unresolved endpoint audit;
see the manuscript and conditional_chen_audit_20261003.md.
Run conditional_chen_product_certificate.py separately for the prime scan.
The bridge cost calculation is retained for audit only: its distribution
hypothesis is REFUTED by conditional_chen_integrality_certificate.py.
"""
import json
from pathlib import Path
from flint import arb, acb, ctx
ctx.prec=192
ga=arb.const_euler(); G=ga.exp(); Umin=(-ga).exp()
L=arb(8500); a=arb(4561)/1000000; eps=arb(1)/500000; v=arb(1)/10**10
s=4-8*a; k=8*(arb(1)/6-a); t=arb(1)/2-a
def f(x): return 2*G*(x-1).log()/x
def h(x): return 3*(-x).exp()/x if x>3 else (-x).exp()
cb=acb.integral(lambda x,analytic:(2-3*x).log(analytic=analytic)/(x*(1-x)),acb(arb(1)/8),acb(arb(1)/3)).real
assert cb < arb('.363084')
beta=arb(7)/2
assert arb(1155).log()<arb('7.1')
assert a*L-beta*L.log()>arb('7.1')
assert a-beta/L>0
assert L/2-beta*L.log()>arb(10**9).log()
assert L/8>1000
assert L*(-L/16).exp() < v
assert 3*(L/8)/(8*arb.pi())*(-L/16).exp() < 1/L**2
# Elementary li envelopes, established at log x=20 and extended by derivatives.
E20=arb(20).exp(); li20=arb(20).ei()
assert E20/20*(1+arb(1)/20+arb(3)/400)>li20
assert li20>E20/20*(1+arb(1)/20)
assert 3/L**2+L**2/(8*arb.pi())*(-L/2).exp()+2*L**2/(arb(2).log()*L.exp()) < 1/L
assert L**2/(72*arb.pi())*(-L/6).exp() < 1/L-27/L**2
low=8*(1-v)*(f(s)-eps*107*arb(1).exp()**2*h(s))
I=G/(4*t)*(6*(3-8*a)/(3-18*a)).log()
FM=2*G/k+eps*106
upper=4*(1+2/L)*(1+v)/(1-(-L/8).exp())*(I+eps*106*(arb(8)/3).log()+2*FM/L**2)
delta=arb(1)/10**8; ell=(1+delta).log()
switch=(1+v)/2*(2*G/t+3*eps*106)*(1+delta)*(1+4/L)*(cb+2/L**2+((arb(8)/3).log()+2/L**2)*(6*ell/L+8/L**2))
# Drop the only negative term of q_G to obtain a decreasing upper envelope.
def qG_upper(log_x):
    w=log_x
    return (arb('.165')+arb('12.683')/w+arb('254.980')/w**2
      +arb('2607.854')/w**3+arb('11605.056')/w**4
      +(arb('1.314')*w+arb('.092')*w.log()+arb('60.883')
        +arb('8.250')*w.log()/w+arb('939.260')/w)*(-w/4).exp())
qG_N=qG_upper(L)
assert qG_N < arb('.17')
assert qG_upper(L/8)<arb('.64')
assert arb('.165')-arb('237.934')/(L/8)*(-L/16).exp()>arb('.16')
# d=1 has r(1)=0; d=2 is excluded since N is even. For d>=3, phi(d)>=2.
per_modulus=arb('.17')+1/(16*arb.pi())+2*(-L/2).exp()/arb(2).log()
assert per_modulus < arb('.195')
assert arb('.608')*per_modulus < arb('.12')
assert per_modulus < arb('.2')
# R_A <= .12 N/L^(5/2). The sum over prime q<y is counted separately:
# R_prime <= .2 N^(5/6)L, avoiding the factor .55 L on the whole R_A.
ap=arb('.18')/(Umin*L.sqrt())+arb('.055')*L**4*(-L/6).exp()/Umin+arb('.8')*(1+v)*FM*L**2*(-L/6).exp()
rect=(1+delta)*L**5/(2*Umin*ell)*(arb('.75')*(-L/48).exp()+(-L/8).exp())
finite=L**2/Umin*(2*(-L/8).exp()+(-2*L/3).exp()+(-L).exp())
margin=low-upper-switch-ap-rect-finite
assert margin > arb('.0018')
assert rect < arb('1e-49')
assert finite < arb('1e-450')
# Both factors in (np,d)=1 must be counted when deleting the restriction.
assert arb('2.2')/(4*arb(2).log())<1
# Uniform elementary envelope for the published bilinear constant.
mb=arb('2.81')*arb('.64')**(arb(1)/3)
mb+=arb('2.809')/(arb('.16')**(arb(2)/3)*(L/8))*(-L/16).exp()
mb+=13*(-L/12).exp()+9*((-5*L/18).exp()+(-L/24).exp())+31*(-13*L/48).exp()
assert mb<3
assert 13*(1+delta)**(arb(1)/3)*(-L/12).exp()<1

# A parameterized EH-style input. Its constants C and starting point X
# remain unknown; this is not a new numerical threshold under ordinary EH.
LE=arb(120); se=arb(9)/4-24/LE
rho_e=1-LE/(8*arb.pi())*(-LE/6).exp()-2*LE*(-LE/3).exp()
main_e=3*Umin*rho_e*(f(se)-107*eps*(-(se-2)).exp())
delete_e=LE**3/arb(2).log()*(2*(-LE/2).exp()+arb('.825')*LE*(-LE).exp())
margin_e=main_e-arb('.05')-delete_e-LE**2*(-LE).exp()
assert se>2 and se<3
assert margin_e>arb('.092')

B=arb(4)*arb(10)**18; LB=B.log(); z=B**(arb(1)/3)
assert B/LB**2 > B.sqrt()*LB/(8*arb.pi())+LB/arb(2).log()
rh_delta=3*z.log()/(8*arb.pi()*z.sqrt())
r=arb('0.000002964')/arb(10**6).log()+arb(1)/(10**6-1)+arb('.54')/10**6
product_eps=r.exp()/(1-rh_delta)-1
assert z>10**6
assert product_eps < arb(1)/700
sb=arb('2.7')-24/LB
fb=f(sb)-arb(117)/700*(-(sb-2)).exp()
assert sb>2 and sb<3 and fb>arb('.072')
rho=1-rh_delta-2*LB/z
mainb=3*Umin*rho*fb
deletions=arb('1.3841')*LB**3/LB.log()*(2/B.sqrt()+arb('0.99')*LB/B)
# log H=.9 L, and 1.1 log H=.99 L.
bridge=mainb-arb('.1')-deletions-LB**2/B
assert bridge > arb('.021')
result={
 'status':'arithmetic certificate only; GRH-only sieve interface pending; former bridge refuted',
 'precision_bits':ctx.prec,
 'GRH_candidate':{'status':'PENDING analytic sieve endpoint audit; not an established GRH-only threshold','log_N0':8500,'log_log_N0':str(L.log()),'alpha':'4561/1000000','beta':'7/2','delta':'1/100000000','epsilon':'1/500000','qG_N_upper_ball':str(qG_N),'per_modulus_upper_ball':str(per_modulus),'squarefree_AP_coefficient_upper_ball':str(arb('.608')*per_modulus),'cbar_ball':str(cb),'lower_ball':str(low),'half_upper_ball':str(upper),'half_switch_ball':str(switch),'AP_budget_ball':str(ap),'half_rectangle_error_ball':str(rect),'finite_error_ball':str(finite),'margin_ball':str(margin),'bilinear_m_upper_ball':str(mb)},
 'EH_parameterized_input':{'status':'PENDING same sieve audit; effective C and X unspecified','input':'for all x>=X, E_(3/4)(x)<= C x/log(x)^4','candidate_threshold':'N >= max(X, exp(max(120,sqrt(20C))))','margin_units':'N/log(N)^2','margin_ball_at_log_N_120':str(margin_e),'main_ball':str(main_e),'deletion_ball':str(delete_e),'note':'ordinary EH does not provide numerical C or X'},
 'bridge':{'B':'4000000000000000000','additional_distribution_hypothesis':'for all x>=B, sum_{d<=x^(9/10)} mu(d)^2 max_{(a,d)=1} |pi(x;d,a)-pi(x)/phi(d)| <= x/(10 log(x)^2)', 'assumptions':'REFUTED hypothesis; cost arithmetic only','product_epsilon_ball':str(product_eps),'s_at_B_ball':str(sb),'sieve_factor_at_B_ball':str(fb),'density_factor_at_B_ball':str(rho),'main_coefficient_ball':str(mainb),'deletion_budget_ball':str(deletions),'margin_ball':str(bridge),'status':'WITHDRAWN: prime-modulus integer counts contradict the assumed upper bound at B; see integrality certificate'}
}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result,ensure_ascii=False,indent=2))
