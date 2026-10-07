"""Strict finite scalar exploration of the refined EH-style interface.

Every candidate still requires its supplied distribution constants C,X.
Search output is not an effective threshold under ordinary EH alone, nor an
independent proof of the underlying analytic sieve lemmas. Selected parameters
must be rerun in the manuscript's 192-bit production certificate.
"""
import json
import math
from fractions import Fraction
from pathlib import Path
from flint import arb, ctx
from conditional_chen_rankin_certificate import certify_rankin

ctx.prec = 192
G = arb.const_euler().exp()
U = arb("1.32") / G
TABLE = [(100000, 106, 107), (10000, 106, 108), (2000, 109, 110),
         (1000, 113, 114), (700, 116, 117), (600, 118, 119),
         (500, 121, 122)]


def endpoint(theta, L, cost, epsilon, eta_plus, eta_minus, C1, C2, budget):
    L = arb(L)
    s = 3 * theta - 3 * cost / L
    if not 2 < s < 3:
        return None
    correction = (2-s).exp()
    low = 2*G*(s-1).log()/s-C2*epsilon*correction
    upper = 2*G/s+C1*epsilon*correction
    assert upper > low
    factor = (1+eta_plus)*low-(eta_plus+eta_minus)*upper
    rho = 1-arb(".000020508")/L-2*L*(-L/3).exp()
    deleted = L**3/arb(2).log()*(2*(-L/2).exp()
              +arb("1.1")*theta*L*(-L).exp())
    return 3*U*rho*factor-arb(1)/budget-deleted-L**2*(-L).exp()


def first_endpoint(theta, cost, epsilon, eta_plus, eta_minus,
                   C1, C2, budget, claim, high):
    low = max(84, math.floor(3*float(cost)/(3*float(theta)-2)))
    value = endpoint(theta, high, cost, epsilon, eta_plus, eta_minus,
                     C1, C2, budget)
    if value is None or not value > arb(claim):
        return None
    # The proof's half-line monotonicity permits binary endpoint selection.
    while high-low > 1:
        mid = (low+high)//2
        value = endpoint(theta, mid, cost, epsilon, eta_plus, eta_minus,
                         C1, C2, budget)
        if value is not None and value > arb(claim):
            high = mid
        else:
            low = mid
    value = endpoint(theta, high, cost, epsilon, eta_plus, eta_minus,
                     C1, C2, budget)
    assert value > arb(claim)
    assert arb(high)/3 > arb(10**12).log()
    assert 2 < 3*theta-3*cost/high < 3
    return high, value


def explore():
    grids = [
        ("EH_9_10", arb(9)/10, 10, ".02", 293,
         [10000, 22000, 30000, 50000, 70000],
         range(36, 46), [(1, 4), (27, 100), (3, 10)]),
        ("EH_3_4", arb(3)/4, 20, ".01", 1054,
         [50000, 70000, 100000, 130000, 200000],
         range(40, 49), [(23, 100), (1, 4), (27, 100)]),
    ]
    best = {}
    count = 0
    uniform = json.loads(Path(__file__).with_name(
        "conditional_chen_uniform_product_certificate.json").read_text())
    uniform_rows = {r["T"]:r for r in uniform["rows"]}
    for key, theta, budget, claim, high, Ts, us, sigmas in grids:
        rows = []
        for T in Ts:
            for w in [3, 5]:
                for numerator in us:
                    u = Fraction(numerator, 8)
                    for sn, sd in sigmas:
                        pre, eps_bound, _, cost = certify_rankin(
                            T, w, u, sn, sd, progress=False,
                            uniform_epsilon=arb(
                                uniform_rows[T]["ordered_uniform_epsilon_ball"]))
                        # Floats propose simple rational bounds; Arb validates
                        # them before any candidate enters the strict search.
                        eps_inverse = math.floor(1/float(eps_bound))
                        ep_inverse = math.floor(
                            1/float(arb(pre["upper_relative_Rosser_gap_ball"])))
                        em_inverse = math.floor(
                            1/float(arb(pre["lower_relative_Rosser_gap_ball"])))
                        epsilon = arb(1)/eps_inverse
                        ep, em = arb(1)/ep_inverse, arb(1)/em_inverse
                        assert eps_bound < epsilon
                        assert arb(pre["upper_relative_Rosser_gap_ball"]) < ep
                        assert arb(pre["lower_relative_Rosser_gap_ball"]) < em
                        choice = next(((a, b) for minimum, a, b in TABLE
                                       if eps_inverse >= minimum), None)
                        if choice is None:
                            continue
                        C1, C2 = choice
                        candidate = first_endpoint(theta, cost, epsilon, ep, em,
                                                   C1, C2, budget, claim, high)
                        count += 1
                        if candidate is not None:
                            start, margin = candidate
                            rows.append({
                                "T": T, "w": w, "u": str(u),
                                "sigma": f"{sn}/{sd}", "log_start": start,
                                "epsilon_inverse": eps_inverse,
                                "eta_plus_inverse": ep_inverse,
                                "eta_minus_inverse": em_inverse,
                                "C1": C1, "C2": C2,
                                "log_support_cost_ball": str(cost),
                                "margin_ball": str(margin),
                                "claimed_margin": claim,
                                "error_inverse": budget,
                                "presieve": pre,
                            })
            print(f"{key}: completed T={T}; tested {count} candidates", flush=True)
        best[key] = sorted(rows, key=lambda r: (r["log_start"],
                          -float(arb(r["margin_ball"]))))[:3]
    return {
        "status": "strict 192-bit scalar search; selected cases require a fresh production certificate",
        "tested_candidates": count,
        "distribution_constants_C_X": "unspecified inputs",
        "uniform_product_input": "saved ordered finite certificate for parameter search; production must regenerate it",
        "not_an_effective_EH_only_numeric_record": True,
        "best": best,
    }


if __name__ == "__main__":
    result = explore()
    Path(__file__).with_suffix(".json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({k: [{a: b for a, b in row.items() if a != "presieve"}
                        for row in rows]
                      for k, rows in result["best"].items()}, indent=2))
