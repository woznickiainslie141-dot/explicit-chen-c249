"""FLOAT ONLY: degree-specific Rankin exploration for distribution inputs."""
import json
import math
from fractions import Fraction
from pathlib import Path
from flint import arb
from conditional_chen_local_rankin_probe import degreewise_gaps

G = math.exp(.5772156649015329)
U = 1.32/G
TABLE = [(100000,106,107),(10000,106,108),(2000,109,110),
         (1000,113,114),(900,114,115),(800,115,116),(700,116,117),
         (600,118,119),(500,121,122),(400,125,126),(300,133,134),
         (200,153,154)]


def endpoint(theta, L, row, epsilon, C1, C2, budget):
    s = 3*theta-3*row["log_support_cost"]/L
    if not 2 < s < 3:
        return -1e6
    low = 2*G*math.log(s-1)/s-C2*epsilon*math.exp(2-s)
    up = 2*G/s+C1*epsilon*math.exp(2-s)
    factor = (1+row["upper_gap"])*low-(
        row["upper_gap"]+row["lower_gap"])*up
    rho = 1-20.508e-6/L-2*L*math.exp(-L/3)
    deleted = L**3/math.log(2)*(2*math.exp(-L/2)+1.1*theta*L*math.exp(-L))
    return 3*U*rho*factor-1/budget-deleted-L**2*math.exp(-L)


def explore():
    uniform = json.loads(Path(__file__).with_name(
        "conditional_chen_uniform_product_certificate.json").read_text())
    uniform_rows = {r["T"]:r for r in uniform["rows"]}
    cases = [(Fraction(3,4),20,.01,915),
             (Fraction(9,10),10,.02,247),
             (Fraction(19,20),10,.02,247),
             (Fraction(99,100),10,.02,247)]
    sigmas = [Fraction(0)]+[Fraction(a,20) for a in range(2,13)]
    best = {str(theta):[] for theta,_,_,_ in cases}
    count = 0
    for T in [10000,22000,50000,70000,100000]:
        epsilon = float(arb(uniform_rows[T]["ordered_uniform_epsilon_ball"]))*1.000001
        C1,C2 = next((a,b) for minimum,a,b in TABLE if epsilon < 1/minimum)
        rows = degreewise_gaps(T,[Fraction(a,8) for a in range(30,47)],
                              sigmas,[3,5,7],M=10)
        for theta,budget,claim,high_limit in cases:
            for row in rows:
                high = high_limit
                if endpoint(float(theta),high,row,epsilon,C1,C2,budget) <= claim:
                    continue
                low = 83
                while high-low > 1:
                    mid = (high+low)//2
                    if endpoint(float(theta),mid,row,epsilon,C1,C2,budget) > claim:
                        high = mid
                    else:
                        low = mid
                value = endpoint(float(theta),high,row,epsilon,C1,C2,budget)
                best[str(theta)].append(dict(row,theta=str(theta),log_start=high,
                    epsilon=epsilon,C1=C1,C2=C2,margin=value,error_inverse=budget))
                count += 1
        print(f"Degreewise EH FLOAT probe completed T={T}",flush=True)
    return {"status":"FLOAT ONLY: no specified C,X; not a numerical EH record",
            "sigma_candidates":[str(s) for s in sigmas],"successful_rows":count,
            "best":{k:sorted(v,key=lambda r:(r["log_start"],-r["margin"]))[:3]
                    for k,v in best.items()}}


if __name__ == "__main__":
    result = explore()
    Path(__file__).with_suffix(".json").write_text(
        json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))
