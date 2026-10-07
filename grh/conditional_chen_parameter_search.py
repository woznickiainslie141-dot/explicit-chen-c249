"""Floating exploration of exact presieving plus truncated Brun weights.

No search output is a theorem or an interval certificate. Accepted manuscript
parameters must be fixed and checked independently with Arb.
"""
import json
import math
from pathlib import Path
from scipy.optimize import minimize_scalar
from functools import lru_cache

G = math.exp(0.5772156649015329)
Umin = 1.32 / G
correction = 0.2740296119016512
cb = 0.36308372924830461
J = math.log(8 / 3)
table = [(74,631,630),(76,559,559),(78,504,504),(80,461,461),
         (85,386,386),(90,336,337),(95,302,302),(100,276,277),
         (120,218,219),(140,189,190),(160,172,173),(180,161,162),
         (200,153,154),(300,133,134),(400,125,126),(500,121,122),
         (600,118,119),(700,116,117),(800,115,116),(900,114,115),
         (1000,113,114),(2000,109,110),(10000,106,108),(100000,106,107)]
small_primes = [p for p in range(3,10001,2)
                if all(p % q for q in range(2,math.isqrt(p)+1))]
finite_density_sum = sum(1/(p-1) for p in small_primes)


@lru_cache(maxsize=None)
def finite_elementary(w, depth):
    values=[1.]+[0.]*depth
    for p in small_primes:
        if p<=w:
            continue
        g=1/(p-1)
        for j in range(depth,0,-1):
            values[j]+=g*values[j-1]
    return tuple(values)


def truncation_gap(w,depth,total_density_upper):
    # Positive coefficient domination for the tail after 10000.
    # This function is exploratory; the accepted parameters use Arb.
    tail=total_density_upper-finite_density_sum
    if tail<0:
        return float('inf')
    return sum(v*tail**(depth-j)/math.factorial(depth-j)
               for j,v in enumerate(finite_elementary(w,depth)))


def presieve(w):
    ps = [p for p in small_primes if p <= w]
    return (sum(math.log(p) for p in ps),
            math.prod(1 - 1/(p-1) for p in ps),
            sum(1/(p-1) for p in ps))


def parameters(L, w, depth, beta, logT):
    T = math.exp(logT)
    logQ, VQ, WQ = presieve(w)
    b = max(2/(math.sqrt(T)*logT), 6.836e-6/math.log(1e12))
    ep = math.expm1(b+2.964e-6/logT+1/(T-1)+.54/T)*1.000001
    eligible = [row for row in table if ep <= 1/row[0]]
    if not eligible:
        return {'margin':-1000}
    row = eligible[-1]
    C1, C2 = row[1:]
    W = math.log(logT)+.262-.5+correction+b-WQ
    gap=truncation_gap(w,depth,W+WQ)
    eta = VQ*logT*math.exp(b)/Umin*gap*1.000001
    cost = logQ+depth*logT
    alpha = (cost+beta*math.log(L)+.0001)/L
    t = .5-alpha
    s = 8*t
    k = 8*(t-1/3)
    if not 3 < s < 4 or k <= 1:
        return {'margin':-1000}
    low = 8*(2*G*math.log(s-1)/s-C2*ep*math.exp(2)*3*math.exp(-s)/s
             -eta*(2*G/3+C1*ep))
    I = G/(4*t)*math.log(6*(3-8*alpha)/(3-18*alpha))
    FM = 2*G/k+C1*ep
    upper = 4*(1+eta)*(1+2/L)*(I+C1*ep*J+2*FM/L**2)
    switch = (1+eta)/2*(2*G/t+3*C1*ep)*(1+1e-8)*(1+4/L)*(cb+2/L**2+(J+2/L**2)*(6*math.log1p(1e-8)/L+8/L**2))
    ap = .18/Umin*L**(3-beta)
    # This explorer only considers L large enough that omitted exponentials
    # are tiny. The Arb acceptance step charges every exponential explicitly.
    return {'margin':low-upper-switch-ap,'L':L,'w':w,'depth':depth,
            'beta':beta,'T':T,'epsilon_bound':ep,'eta_bound':eta,
            'alpha_min':alpha,'support_log_cost':cost,
            'C1':C1,'C2':C2,'table_row':row[0],
            'lower':low,'half_upper':upper,'half_switch':switch,'AP':ap}


def search(L):
    best = []
    for w in [2,3,5,7,11,13,17,19,23,29,31,37,43,47,59,73,97,127,199]:
        for depth in range(4,19,2):
            for beta in [3.5,3.75,4,4.25,4.5]:
                f = lambda u:-parameters(L,w,depth,beta,u)['margin']
                opt = minimize_scalar(f,bounds=(math.log(10000),math.log(1e10)),method='bounded')
                best.append(parameters(L,w,depth,beta,opt.x))
    return sorted(best,key=lambda x:x['margin'],reverse=True)[:5]


def eh_parameters(theta,L,w,depth,logT,error_budget):
    T = math.exp(logT)
    logQ,VQ,WQ = presieve(w)
    b = max(2/(math.sqrt(T)*logT),6.836e-6/math.log(1e12))
    ep = math.expm1(b+2.964e-6/logT+1/(T-1)+.54/T)*1.000001
    eligible = [row for row in table if ep <= 1/row[0]]
    if not eligible:
        return {'margin':-1000}
    C1,C2 = eligible[-1][1:]
    W = math.log(logT)+.262-.5+correction+b-WQ
    gap=truncation_gap(w,depth,W+WQ)
    eta = VQ*logT*math.exp(b)/Umin*gap*1.000001
    cost = logQ+depth*logT
    s = 3*theta-3*cost/L
    if not 2<s<3:
        return {'margin':-1000}
    factor = 2*G*math.log(s-1)/s-C2*ep*math.exp(2-s)-eta*(2*G/s+C1*ep)
    main = 3*Umin*(1-20.508e-6/L-2*L*math.exp(-L/3))*factor
    deleted = L**3/math.log(2)*(2*math.exp(-L/2)+1.1*theta*L*math.exp(-L))
    return {'margin':main-error_budget-deleted-L**2*math.exp(-L),
            'theta':theta,'L':L,'w':w,'depth':depth,'T':T,
            'epsilon_bound':ep,'eta_bound':eta,'C1':C1,'C2':C2,
            's':s,'support_log_cost':cost,'main':main}


def search_eh(theta,L,error_budget):
    best=[]
    for w in [2,3,5,7,11,13,17,19,23,29,31,37,43,47,59,73,97]:
        for depth in range(4,15,2):
            lower=math.log(10000)
            upper=min(math.log(1e9),(L*(theta-2/3)-presieve(w)[0])/depth-.000001)
            if upper<=lower:
                continue
            f=lambda u:-eh_parameters(theta,L,w,depth,u,error_budget)['margin']
            grid=[lower+(upper-lower)*i/24 for i in range(25)]
            idx=min(range(len(grid)),key=lambda i:f(grid[i]))
            opt=minimize_scalar(f,bounds=(grid[max(0,idx-1)],grid[min(24,idx+1)]),method='bounded')
            for u in [opt.x,grid[idx]]:
                best.append(eh_parameters(theta,L,w,depth,u,error_budget))
    return sorted(best,key=lambda x:x['margin'],reverse=True)[:3]


if __name__ == '__main__':
    result = {'status':'FLOATING EXPLORATION ONLY',
              'exact_presieve_composite':{str(L):search(L) for L in (50000,53000,54000)},
              'EH_exploration':{str(theta)+','+str(L):search_eh(theta,L,budget)
                                for theta,L,budget in [(0.75,1000,.05),(0.75,1500,.05),
                                                       (0.9,400,.1),(0.9,500,.1)]}}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
