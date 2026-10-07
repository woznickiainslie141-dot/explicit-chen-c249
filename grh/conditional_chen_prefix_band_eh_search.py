"""Strict scalar search using the preceding-pass product band."""
import json
import math
from fractions import Fraction
from pathlib import Path
from flint import arb,ctx
from conditional_chen_prefix_band_certificate import certify_band
from conditional_chen_refined_eh_search import first_endpoint,TABLE

ORDINARY = [Fraction(0)]+[Fraction(a,20) for a in range(2,13)]
BANDS = [Fraction(0),Fraction(3,25),Fraction(6,25),Fraction(9,25)]


def explore():
    ctx.prec = 192
    uniform = json.loads(Path(__file__).with_name(
        "conditional_chen_uniform_product_certificate.json").read_text())
    uniform_rows = {r["T"]:r for r in uniform["rows"]}
    cases = [("EH_3_4",arb(3)/4,20,".01",833),
             ("EH_9_10",arb(9)/10,10,".02",227),
             ("EH_19_20",arb(19)/20,10,".02",179),
             ("EH_99_100",arb(99)/100,10,".02",152)]
    best = {key:[] for key,_,_,_,_ in cases}
    tested = 0
    for T in [10000,22000]:
        for w in [3,5]:
            for numerator in [33,34,35,36,39,40]:
                pre = certify_band(T,w,Fraction(numerator,8),ORDINARY,BANDS,10,False)
                bounds = [arb(uniform_rows[T]["ordered_uniform_epsilon_ball"]),
                          arb(pre["upper_relative_Rosser_gap_ball"]),
                          arb(pre["lower_relative_Rosser_gap_ball"])]
                inverses = [math.floor(1/float(v)) for v in bounds]
                epsilon,ep,em = [arb(1)/n for n in inverses]
                assert all(v < arb(1)/n for v,n in zip(bounds,inverses))
                C1,C2 = next((a,b) for minimum,a,b in TABLE if inverses[0] >= minimum)
                cost = arb(pre["log_support_cost_ball"])
                for key,theta,budget,claim,high in cases:
                    candidate = first_endpoint(theta,cost,epsilon,ep,em,C1,C2,budget,claim,high)
                    tested += 1
                    if candidate is None: continue
                    start,margin = candidate
                    best[key].append({
                        "T":T,"w":w,"u":str(Fraction(numerator,8)),
                        "sigma_candidates":[str(s) for s in ORDINARY],
                        "band_sigma_candidates":[str(s) for s in BANDS],
                        "degree_cutoff":10,"log_start":start,
                        "epsilon_inverse":inverses[0],"eta_plus_inverse":inverses[1],
                        "eta_minus_inverse":inverses[2],"C1":C1,"C2":C2,
                        "error_inverse":budget,"claimed_margin":claim,
                        "margin_ball":str(margin),"log_support_cost_ball":str(cost),
                        "presieve":pre})
                print(f"Prefix-band EH strict search T={T}, w={w}, u={Fraction(numerator,8)}",
                      flush=True)
    return {"status":"strict scalar search; production must regenerate selected inputs",
            "tested_candidates":tested,"not_an_effective_EH_only_numeric_record":True,
            "distribution_constants_C_X":"unspecified",
            "best":{k:sorted(v,key=lambda r:(r["log_start"],-float(arb(r["margin_ball"]))))[:3]
                    for k,v in best.items()}}


if __name__ == "__main__":
    result = explore()
    Path(__file__).with_suffix(".json").write_text(json.dumps(result,indent=2)+"\n",
                                                encoding="utf-8")
    print(json.dumps({k:[{a:b for a,b in r.items() if a != "presieve"} for r in rows]
                      for k,rows in result["best"].items()},indent=2))
