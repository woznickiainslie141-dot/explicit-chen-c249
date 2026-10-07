"""Recharge fully generated finite inputs using the proved Rosser support bound.

Run conditional_chen_composite_certificate.py first to regenerate the finite
prime/product inputs. This script reuses those saved inputs with recorded
SHA-256 hashes, scans the complete medium-prime sets again for kappa, and checks
all scalar conditions and expenses at the new endpoints. It does not certify
the external analytic lemmas or supply numerical EH constants C,X.
"""
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from flint import arb, acb, ctx
from conditional_chen_sharp_support_certificate import run as support_certificate
from conditional_chen_refined_eh_search import first_endpoint

ROOT = Path(__file__).resolve().parent
ctx.prec = 192
gamma = arb.const_euler()
G = gamma.exp()
Umin = arb("1.32")*(-gamma).exp()


def qG_upper(w):
    return (arb(".165")+arb("12.683")/w+arb("254.980")/w**2
            +arb("2607.854")/w**3+arb("11605.056")/w**4
            +(arb("1.314")*w+arb(".092")*w.log()+arb("60.883")
              +arb("8.250")*w.log()/w+arb("939.260")/w)*(-w/4).exp())


def support_input(case, rows):
    pre = case["presieve"]
    row = rows[(pre["T"], pre["exact_cutoff"], pre["u"])]
    assert pre["prime_generation"] == "complete integer sieve"
    assert pre["logarithm_base"] == "e"
    assert pre["band_optimizers_with_strict_signs"] > 0
    assert row["complete_medium_prime_count"] == pre["medium_prime_count"]
    assert row["Q"] == pre["exact_presieve_Q_max"]
    coarse = arb(pre["log_support_cost_ball"])
    assert abs(coarse-arb(row["original_log_support_cost_ball"])) < arb("1e-50")
    kappa = arb(row["kappa"])
    sharp = coarse+kappa.log()
    assert 0 < sharp < coarse
    assert abs(sharp-arb(row["sharp_log_support_cost_ball"])) < arb("1e-50")
    assert 0 < arb(pre["uniform_epsilon_ball"]) < arb(case["epsilon"])
    assert arb(pre["upper_relative_Rosser_gap_ball"]) < arb(case["eta_plus"])
    assert arb(pre["lower_relative_Rosser_gap_ball"]) < arb(case["eta_minus"])
    return sharp, row


def grh_case(case, cost, support_row, log_N0, alpha_numerator,
             rectangle_bound, finite_bound, claim=".0006"):
    L = arb(log_N0)
    alpha = arb(alpha_numerator)/10**6
    beta = arb(7)/2
    epsilon = arb(case["epsilon"])
    eta, eta_minus = arb(case["eta_plus"]), arb(case["eta_minus"])
    assert (case["C1"], case["C2"]) == (106, 108)
    C1, C2 = 106, 108
    assert epsilon < arb(1)/10000
    t = arb(1)/2-alpha
    s, k = 8*t, 8*(t-arb(1)/3)
    delta, v = arb(1)/10**8, arb(1)/10**10
    ell = (1+delta).log()
    assert alpha*L-beta*L.log() > cost
    assert alpha-beta/L > 0
    assert L/2-beta*L.log() > arb(10**9).log()
    assert L/8 > arb(case["presieve"]["T"]).log()
    assert 3 < s < 4 and k > 1
    assert L*(-L/16).exp() < v
    assert 3*(L/8)/(8*arb.pi())*(-L/16).exp() < 1/L**2
    E20, li20 = arb(20).exp(), arb(20).ei()
    assert E20/20*(1+arb(1)/20+arb(3)/400) > li20
    assert li20 > E20/20*(1+arb(1)/20)
    assert 3/L**2+L**2/(8*arb.pi())*(-L/2).exp()+2*L**2/(arb(2).log()*L.exp()) < 1/L
    assert L**2/(72*arb.pi())*(-L/6).exp() < 1/L-27/L**2
    assert qG_upper(L) < arb(".17")
    assert qG_upper(arb(8500)) < arb(".166496")
    assert qG_upper(L/8) < arb(".64")
    assert arb(".165")-arb("237.934")/(L/8)*(-L/16).exp() > arb(".16")
    per_modulus = arb(".17")+1/(16*arb.pi())+2*(-L/2).exp()/arb(2).log()
    assert arb(".608")*per_modulus < arb(".12") and per_modulus < arb(".2")
    cb = acb.integral(
        lambda x, analytic: (2-3*x).log(analytic=analytic)/(x*(1-x)),
        acb(arb(1)/8), acb(arb(1)/3)).real
    assert cb < arb(".363084")
    correction = acb.integral(
        lambda x, analytic: (x-2).log(analytic=analytic)/(x-1),
        acb(3), acb(s)).real
    Fs = 2*G*(1+correction)/s
    fs = 2*G*(s-1).log()/s
    h = 3*(-s).exp()/s
    Flo, large_lower = Fs+C1*epsilon*arb(2).exp()*h, fs-C2*epsilon*arb(2).exp()*h
    assert 0 < large_lower <= Flo
    assert fs <= Fs < 2*G/3
    composite_lower = (1+eta)*large_lower-(eta+eta_minus)*Flo
    assert composite_lower > 0
    low = 8*(1-v)*composite_lower
    I = G/(4*t)*(6*(3-8*alpha)/(3-18*alpha)).log()
    FM, J = 2*G/k+C1*epsilon, (arb(8)/3).log()
    b0 = t-arb(1)/4
    assert arb(1)/8 < b0 < arb(1)/3
    Hh = acb.integral(
        lambda b, analytic: (2-8*(acb(t)-b)).exp()/b,
        acb(arb(1)/8), acb(b0)).real+((arb(1)/3)/b0).log()
    assert 0 < Hh < J
    upper = 4*(1+eta)*(1+2/L)*(1+v)/(1-(-L/8).exp())*(I+C1*epsilon*Hh+2*FM/L**2)
    switch = (1+eta)*(1+v)/2*(2*G/t+3*C1*epsilon)*(1+delta)*(1+4/L)*(cb+2/L**2+(J+2/L**2)*(6*ell/L+8/L**2))
    ap = arb(".18")/(Umin*L.sqrt())+arb(".055")*L**4*(-L/6).exp()/Umin+arb(".8")*(1+eta)*(1+v)*FM*L**2*(-L/6).exp()
    rect = (1+delta)*L**5/(2*Umin*ell)*(arb(".75")*(-L/48).exp()+(-L/8).exp())
    finite = L**2/Umin*(2*(-L/8).exp()+(-2*L/3).exp()+(-L).exp())
    margin = low-upper-switch-ap-rect-finite
    assert margin > arb(claim)
    assert rect < arb(rectangle_bound) and finite < arb(finite_bound)
    assert arb("2.2")/(4*arb(2).log()) < 1
    mb = arb("2.81")*arb(".64")**(arb(1)/3)
    mb += arb("2.809")/(arb(".16")**(arb(2)/3)*(L/8))*(-L/16).exp()
    mb += 13*(-L/12).exp()+9*((-5*L/18).exp()+(-L/24).exp())+31*(-13*L/48).exp()
    assert mb < 3
    assert 13*(1+delta)**(arb(1)/3)*(-L/12).exp() < 1
    return {
        **case, "log_N0": log_N0, "alpha": f"{alpha_numerator}/1000000",
        "claimed_margin": claim,
        "log_log_N0_ball": str(L.log()), "s_ball": str(s),
        "support_slack_ball": str(alpha*L-beta*L.log()-cost),
        "effective_log_support_cost_ball": str(cost), "sharp_support": support_row,
        "F_s_ball": str(Fs), "h_integral_ball": str(Hh),
        "lower_composite_factor_ball": str(composite_lower),
        "lower_ball": str(low), "half_upper_ball": str(upper),
        "half_switch_ball": str(switch), "AP_budget_ball": str(ap),
        "half_rectangle_error_ball": str(rect), "finite_error_ball": str(finite),
        "margin_ball": str(margin), "rectangle_display_bound": rectangle_bound,
        "finite_display_bound": finite_bound,
        "sieve_interface": "support <kappa Q D0 D; unchanged finite weights, gaps and uniform product inputs",
    }


def eh_case(case, cost, support_row, theta, budget):
    epsilon, ep, em = (arb(case[x]) for x in ["epsilon", "eta_plus", "eta_minus"])
    C1, C2 = case["C1"], case["C2"]
    assert (C1, C2) == (113, 114) and epsilon < arb(1)/1000
    candidate = first_endpoint(theta, cost, epsilon, ep, em, C1, C2,
                               budget, case["claimed_margin"], case["log_start"])
    assert candidate is not None
    start, searched_margin = candidate
    L = arb(start)
    s = 3*theta-3*cost/L
    assert 2 < s < 3
    assert L/3 > arb(10**12).log()
    assert L/3 > arb(case["presieve"]["T"]).log()
    assert theta*L > arb(10**9).log()
    low = 2*G*(s-1).log()/s-C2*epsilon*(2-s).exp()
    upper = 2*G/s+C1*epsilon*(2-s).exp()
    assert 0 < low <= upper
    factor = (1+ep)*low-(ep+em)*upper
    assert factor > 0
    rho = 1-arb(".000020508")/L-2*L*(-L/3).exp()
    main = 3*Umin*rho*factor
    deleted = L**3/arb(2).log()*(2*(-L/2).exp()+arb("1.1")*theta*L*(-L).exp())
    margin = main-arb(1)/budget-deleted-L**2*(-L).exp()
    assert margin > arb(case["claimed_margin"])
    assert abs(margin-searched_margin) < arb("1e-50")
    return {
        **case, "log_start": start,
        "threshold_formula": f"max(X, exp(max({start}, sqrt({budget} C))))",
        "effective_log_support_cost_ball": str(cost), "sharp_support": support_row,
        "s_at_start_ball": str(s), "sieve_factor_ball": str(factor),
        "main_ball": str(main), "surplus_ball": str(margin),
    }


def run():
    if not __debug__:
        raise RuntimeError("Run without Python -O.")
    base = json.loads((ROOT/"conditional_chen_composite_certificate.json").read_text())
    uniform = json.loads((ROOT/"conditional_chen_uniform_product_certificate.json").read_text())
    standalone = json.loads((ROOT/"conditional_chen_prefix_band_certificate.json").read_text())
    assert base["precision_bits"] == 192
    assert base["uniform_product_certificate"]["regenerated_in_this_run"] is True
    assert base["GRH"]["presieve"] == standalone
    assert uniform["complete_odd_prime_count"] == 5761454
    assert uniform["checked_prime_jumps"] == 5760226
    uniform_rows = {row["T"]: row for row in uniform["rows"]}
    for key in ["GRH", "EH_3_4", "EH_9_10", "EH_19_20", "EH_99_100"]:
        pre = base[key]["presieve"]
        assert pre["uniform_epsilon_ball"] == uniform_rows[pre["T"]]["ordered_uniform_epsilon_ball"]
    support = support_certificate()
    rows = {(row["T"], row["exact_cutoff"], row["u"]): row for row in support["cases"]}
    # Reproduce all baseline GRH arithmetic before charging the new support.
    original = base["GRH"]
    coarse = arb(original["presieve"]["log_support_cost_ball"])
    _, sr = support_input(original, rows)
    scaled_alpha = Fraction(original["alpha"])*10**6
    assert scaled_alpha.denominator == 1
    repeated = grh_case(original, coarse, sr, original["log_N0"],
                        scaled_alpha.numerator, "1e-254", "1e-1700")
    for key in ["margin_ball", "lower_ball", "half_upper_ball", "half_switch_ball",
                "AP_budget_ball", "F_s_ball", "h_integral_ball", "support_slack_ball"]:
        assert abs(arb(repeated[key])-arb(original[key])) < arb("1e-45")
    cost, sr = support_input(original, rows)
    result = {
        "status": "PASS: strict endpoint recharging; saved finite inputs with recorded provenance",
        "precision_bits": ctx.prec,
        "GRH": grh_case(original, cost, sr, 31000, 4215, "1e-248", "1e-1670"),
        "baseline_GRH_arithmetic_reproduction": "PASS: eight quantities agree within 1e-45",
        "finite_inputs_regenerated_in_this_run": False,
        "support_maxima_regenerated_in_this_run": True,
        "support_certificate": support,
    }
    for key, theta, budget in [("EH_3_4", arb(3)/4, 20), ("EH_9_10", arb(9)/10, 10),
                               ("EH_19_20", arb(19)/20, 10), ("EH_99_100", arb(99)/100, 10)]:
        cost, sr = support_input(base[key], rows)
        result[key] = eh_case(base[key], cost, sr, theta, budget)
    cost, sr = support_input(original, rows)
    result["GRH_positivity"] = grh_case(
        original, cost, sr, 30824, 4238, "1e-248", "1e-1660", ".00001")
    result["distribution_positivity"] = {}
    for key, theta, budget in [("EH_3_4", arb(3)/4, 20), ("EH_9_10", arb(9)/10, 10),
                               ("EH_19_20", arb(19)/20, 10), ("EH_99_100", arb(99)/100, 10)]:
        cost, sr = support_input(base[key], rows)
        weaker = {**base[key], "claimed_margin": ".00001"}
        result["distribution_positivity"][key] = eh_case(weaker, cost, sr, theta, budget)
    inputs = ["conditional_chen_composite_certificate.py", "conditional_chen_composite_certificate.json",
              "conditional_chen_prefix_band_certificate.py", "conditional_chen_prefix_band_certificate.json",
              "conditional_chen_uniform_product_certificate.py", "conditional_chen_uniform_product_certificate.json",
              "conditional_chen_rankin_certificate.py",
              "conditional_chen_uniform_grh_probe.py", "conditional_chen_adaptive_probe.py",
              "conditional_chen_sharp_support_certificate.py", "conditional_chen_refined_eh_search.py",
              Path(__file__).name]
    result["input_sha256"] = {
        name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in inputs}
    result["not_checked"] = ["external analytic lemmas", "independent research review",
                             "effective EH constants C,X", "full splice", "PDF compilation"]
    return result


if __name__ == "__main__":
    result = run()
    Path(__file__).with_suffix(".json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({key: {name: value for name, value in row.items()
                           if name in ["log_N0", "log_start", "alpha", "margin_ball", "surplus_ball",
                                       "support_slack_ball", "effective_log_support_cost_ball"]}
                      for key, row in result.items()
                      if key.startswith("GRH") or key.startswith("EH_")}, indent=2))
    print(json.dumps({key: {"log_start": row["log_start"],
                           "surplus_ball": row["surplus_ball"]}
                      for key, row in result["distribution_positivity"].items()}, indent=2))
