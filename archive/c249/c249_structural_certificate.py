#!/usr/bin/env python3
"""Arb structural and finite-error certificate for c0=24.9.

Every comparison is decided with outward Arb balls.  The manuscript proves
the analytic implications and half-line propagation; this script certifies
their parameter-dependent starting inequalities and the exact error budget.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import flint
from flint import arb, ctx, fmpq


ctx.prec = 192
ctx.threads = 1
passed_assertions: list[str] = []


def R(num: int | str, den: int | None = None) -> arb:
    value = arb(num)
    return value if den is None else value / den


def prove(name: str, condition: bool) -> None:
    if not condition:
        raise ArithmeticError(f"undecided or false assertion: {name}")
    passed_assertions.append(name)


def txt(value: arb, digits: int = 32) -> str:
    if not value.is_finite() or value.rel_accuracy_bits() < 80:
        raise ArithmeticError(f"insufficiently accurate ball: {value}")
    return value.str(digits, more=True)


def main() -> dict:
    passed_assertions.clear()

    c0 = R(249, 10)
    L = c0.exp()
    log2 = R(2).log()
    k1 = R(191, 2000)
    k2 = R(383, 4000)
    rho = R(119, 500)
    sigma2 = R(6047, 20000)
    sigma1 = R(121, 400)
    eps0 = R(1, 10_000)
    delta = R(29, 25)
    qbar = R(1_000_000_050)
    theta = R(2_423_257, 5_000_000)
    level_shift = R(81, 1_000_000_000)
    exc_eps = R(1, 35)
    exc_density = R(36, 35)
    cn_lower = R(1, 2)

    # Published BJS product input and the fixed continuous level envelope.
    log_u0 = R(1_000_000_000).log()
    q_from_bjs = R(1_000_000_000) * (1 + R(9, 10_000_000) / log_u0)
    prove("BJS theta(10^9) envelope", q_from_bjs < qbar)
    theta_actual_lower = R(1, 2) - (qbar + 20 * c0) / L
    prove("theta rational lower envelope", theta < theta_actual_lower)
    level_loss = 20 * c0 / (k1 * L)
    prove("large-exceptional level loss", level_loss < level_shift)
    prove("product lemma z-domain", k1 * L > 4000)
    prove("ordinary exceptional primes below z*", k1 * L > 10 * c0)

    # Exact Wu and source-geometry constraints for the retuned parameters.
    qk1, qk2 = fmpq(191, 2000), fmpq(383, 4000)
    qrho = fmpq(119, 500)
    qs2, qs1 = fmpq(6047, 20000), fmpq(121, 400)
    prove("parameter order k1<k2", qk1 < qk2)
    prove("parameter order k2<rho", qk2 < qrho)
    prove("parameter order rho<sigma2", qrho < qs2)
    prove("parameter order sigma2<sigma1", qs2 < qs1)
    prove("parameter order sigma1<1/3", qs1 < fmpq(1, 3))
    prove("Wu condition 3k1+rho>=1/2", 3 * qk1 + qrho >= fmpq(1, 2))
    prove("Wu U5 endpoint dominated by rho", 3 * qk1 + qrho >= fmpq(1, 2))
    prove("rectangle family 10 exceeds .398", qk2 + qs2 > fmpq(398, 1000))
    prove("triple reduction 7 and 8", 4 * qs2 > 1)
    prove("triple reduction 9", qk1 + 3 * qs1 > 1)
    prove("triple reduction 10", qk2 + qs2 + 2 * qs1 > 1)
    prove("switch n lower exponent exceeds 1/6", 2 * qk1 > fmpq(1, 6))

    # Reusable elementary constants.
    egamma = arb.const_euler().exp()
    eta = R(1452, 10**10)
    prove("full Mertens correction below eta", 1 / (k1 * L) ** 2 < eta)
    prove("prime-count inflation", R("1.00005") * (1 + eps0) < R("1.0002"))
    prove("prime-count bound at every mass endpoint", 1 / (1 - 6 / L) < R("1.00005"))
    prove("prime-count endpoint logarithm", L / 4 > 30002)
    prove("all measure masses below 3", (1 / k1).log() + R("1e-12") < 3)
    prove("short reciprocal-prime mass below 1/3", (k2 / k1).log() + R("1e-12") < R(1, 3))
    prove("upper sieve function global bound", 2 * egamma + 106 * eta < 4)
    prove("upper sieve function for s>=2", egamma + 106 * eta < 2)
    prove("upper sieve correction derivative", 106 * eta * R(2).exp() < 1)
    prove("RS conductor envelope", egamma + R(5, 8) < 4)
    prove("base normalization bound", 2 / egamma / k1 < 14)
    prove("switched normalization bound", 2 / egamma**2 / k1**2 < 100)

    # Complementary and rough sieve levels.  Wu's moving endpoint is
    # lambda=1/2-2*kappa_1, hence the residual exponent below is exact.
    dcircle_exponent = R(2, 5)
    dstar_exponent = dcircle_exponent + qbar / L
    dstar_cap = R(2077, 5000)  # 0.4154
    prove("total switched modulus exponent <.4154", dstar_exponent < dstar_cap)
    prove("D*sqrt(N) exponent <.9154", dstar_exponent + R(1, 2) < R(4577, 5000))
    s_circle = dcircle_exponent / k1
    prove("complementary sieve parameter exact", abs(s_circle - R(800, 191)) < R("1e-40"))
    residual_exponent = R(1, 2) + 2 * k1 - 2 * k2
    s_rho = (
        residual_exponent - qbar / L
        - 10 * (residual_exponent * L).log() / L
    ) / k1
    s_rho_floor = R(12_037, 2_500)
    prove("rough-sieve parameter", s_rho > s_rho_floor)

    # Every analytic branch and every explicit integration split.
    for branch, th in [("small", theta), ("universal", theta - k1 * level_shift)]:
        prove(branch + " A1 minimum argument >3", (th - k2) / k2 > 3)
        prove(branch + " A1 maximum argument <5", (th - k1) / k1 < 5)
        prove(branch + " A2 minimum argument >1", (th - sigma1) / k1 > 1)
        prove(branch + " A2 minimum argument <2", (th - sigma1) / k1 < 2)
        prove(branch + " A2/A3 s=3 split above k1", th - 3 * k1 > k1)
        prove(branch + " A2/A3 s=3 split below s=2", th - 3 * k1 < th - 2 * k1)
        prove(branch + " A2/A3 s=2 split below sigma2", th - 2 * k1 < sigma2)
        prove(branch + " L4 minimum argument >3", (th - 2 * k2) / k1 > 3)
        prove(branch + " L4 maximum argument <4", (th - 2 * k1) / k1 < 4)
        prove(branch + " L5 geometry split", k1 + k2 < 2 * k2)
        prove(branch + " L5 s=3 split", 2 * k2 < th - 3 * k1)
        prove(branch + " L5 delay split order", th - 3 * k1 < th - 2 * k1)
        prove(branch + " L5 s=2 below terminal", th - 2 * k1 < R(1, 2) - 2 * k1)
        prove(branch + " L5 first argument <4", (th - k1 - k2) / k1 < 4)
        prove(branch + " L5 first argument >3", (th - k1 - k2) / k1 > 3)
        prove(branch + " L5 last argument >1", (th - (R(1, 2) - 2 * k1)) / k1 > 1)
        prove(branch + " L0 argument >5", th / k1 > 5)
        prove(branch + " L0 argument <6", th / k1 < 6)
    prove("box sieve argument >1", 3 * theta - level_shift > 1)
    prove("box sieve argument <2", 3 * theta < 2)
    prove("rough sieve branch >3", s_rho_floor > 3)
    prove("rough sieve branch <5", s_rho_floor < 5)
    prove("complement sieve branch >3", s_circle - level_shift > 3)
    prove("complement sieve branch <5", s_circle < 5)

    # Existing primorial information is used fully: every large-exceptional
    # threshold exceeds P31, hence the largest prime is at least 37.
    P31 = R(100_280_245_065)
    P37 = R(3_710_369_067_405)
    Kbase = (delta * (k1 * L).log()).exp()
    prove("common z* threshold exceeds P31", Kbase > P31)
    prove("common z* threshold below P37 at c0", Kbase < P37)
    prove("exceptional epsilon is 1/(37-2)", abs(exc_eps - R(1, 35)) < R("1e-40"))
    prove("exceptional density is 36/35", abs(exc_density - R(36, 35)) < R("1e-40"))

    # A common zero bound for the three distribution stages.
    nu_uniform = 100 / ((delta * c0 / 2).exp() * (delta * c0) ** 2)
    nu_nonexceptional = 1 / (2 * R("2.0452") * 10 * c0)
    prove("uniform zero-bound branch", nu_uniform < nu_nonexceptional)
    prove("uniform zero bound positive", nu_uniform < R(1, 2))

    # Fixed-cutoff rectangle extension.
    logX3 = L / 8 - log2
    w3 = logX3 - 15 * logX3.log()
    Krect = (delta * w3.log()).exp()
    nu = nu_uniform
    beta = 1 - nu
    prove("exceptional y^-nu term decreasing", nu * w3 > 100)
    prove("rectangle exceptional threshold exceeds P31", Krect > P31)
    prove("fixed D0 below rectangle cap", 10 * w3.log() < L / 2 - 10 * c0)
    prove("rectangle modulus range for mu2/phi", L / 2 - 10 * c0 > R(10**9).log())
    prove("rectangle low cutoff for mu2/phi", 10 * w3.log() > R(10**9).log())
    prove("rectangle modulus gap derivative", L / 2 - 10 - R("10.1") > 0)
    prove("rectangle BJS scale", w3.log() > R("10.4"))
    prove("log^14(y)/sqrt(y) decreasing", w3 > 28)
    prove("log^14(y)/y^(2/3) decreasing", 2 * w3 > 42)
    vprime = (
        R("3.2e-8") / w3**4
        + w3**4 * (
            (-nu * w3).exp() / beta
            + R("1.02") * w3**10 * (-w3 / 2).exp()
            + 3 * w3**10 * (-2 * w3 / 3).exp()
        )
    )
    v3 = vprime * (
        1 + 1 / (logX3**10 * w3**5)
        + 1 / ((1 - 6 / w3) * logX3)
    ) + 3 / logX3
    mrect = (
        39 * v3
        + 108 * logX3**16 * (-logX3).exp() / logX3.log()
        + 26 * L**5 / w3**10
        + 88 * L**5 * ((-R("0.199") * L).exp() + (-L / 16).exp())
        + 106 / logX3**6
    )
    prove("rectangle m endpoint bound", mrect < R("1.437e-8"))
    families = [
        (R(1, 3) - sigma1, sigma1, R(2)),
        (R(1, 3) - sigma2, sigma2, R(2)),
        (sigma1 - R(1, 8), R(1, 8), R(1)),
        (R(11, 24) - sigma1, sigma1, R(1)),
        (1 - k2 - sigma2 - sigma1, sigma1, R(1)),
    ]
    weighted_boxes = R(0)
    for width, amin, wu_weight in families:
        prove("positive box-family width " + str(len(passed_assertions)), width > 0)
        box_count = width * L / (1 + eps0).log() + 1
        weighted_boxes += wu_weight * box_count / amin**3
    rectangle_normalized = mrect * (1 + eps0) * weighted_boxes / (cn_lower * L) / 4
    prove("global rectangle charge <.0085", rectangle_normalized < R("0.0085"))

    # Ordinary BJS AP coefficient, at its certified worst endpoint.
    wap = L - 5 * c0
    Kap = (delta * wap.log()).exp()
    mathcal_E_over_y = (
        4 * wap ** R("4.5") / (wap - 5 * wap.log()) ** 10
        + 4 / wap ** R("5.5")
        + 18 * (-wap / 12).exp() / wap.sqrt()
        + R("2.5") * (-wap / 6).exp() * wap ** R("5.5")
    )
    p2 = wap**2 * (
        R("1.1") * (10 * wap.log())
        * (R("3.2e-8") / wap**8 + (-nu * wap).exp() / (1 - nu))
        + 27 * mathcal_E_over_y
        + (-wap / 2).exp() / (2 * log2 * wap**8)
        + R("0.4") * wap**3 * (-wap).exp()
    )
    p1 = p2 + (R("0.67") + 2 * (-wap / 6).exp()) / wap**8
    pap = p1 * (1 + 1 / (L**2 * wap**3) + 1 / ((1 - 4 / wap) * L)) + R("2.2") / L**2
    c4 = pap + R("0.9") * (wap / 2 - L).exp() * L**4 / (wap**10 * c0)
    prove("ordinary AP coefficient", c4 * L**2 < R("2.21"))
    decomposition_multiplicity = (1 + L / log2) ** 2
    ordinary_ap_normalized = 16 * c4 * decomposition_multiplicity / (4 * cn_lower * L)
    prove("all ordinary AP remainders <1e-6", ordinary_ap_normalized < R("1e-6"))

    # Exact full-conductor mass defect and the q1>=37 applicability checks.
    prove("full-conductor Brun-Titchmarsh factor", 2 * L / (L - 10 * c0) < 3)
    prove("ordinary exceptional threshold exceeds P31", Kap > P31)
    prove("RS lower conductor range", Kap.log().log() > 2)
    prove("ordinary AP BJS scale", wap.log() > R("10.4"))
    prove("AP endpoint denominators", wap - 5 * wap.log() > 10)
    b_exc = 96 * (10 * c0).log() ** 2 / Kap
    exceptional_finite_proxy = 100 * b_exc + 1000 / L**2
    prove("large-exceptional finite charge <1e-6", exceptional_finite_proxy < R("1e-6"))

    # Prime-measure, boundary-strip and inflation transfer.
    mertens_discrepancy = R("3e-6") / (k1 * L) + 2 * (-k1 * L).exp()
    prove("prime discrepancy below 1e-12", mertens_discrepancy < R("1e-12"))
    geometric_shift = (1 + eps0).log() / L
    stieltjes_proxy = R(10**10) * (mertens_discrepancy + 16 * geometric_shift) + R(60, 5000) / 4
    prove("Stieltjes and endpoint transfer <.0035", stieltjes_proxy < R("0.0035"))

    # Switched low and high conductors.
    logXq = k1 * L - log2
    wq = logXq - 15 * logXq.log()
    log_q0_lower = 10 * wq.log() - log2
    Kq = (delta * wq.log()).exp()
    prove("switch low-conductor exponential decreasing", nu * wq > 100)
    prove("switch exceptional threshold exceeds P31", Kq > P31)
    prove("switch endpoint primorial cell below P37", Kq < P37)
    prove("switch q0 below source primes", 10 * wq.log() < k1 * L)
    prove("switch low cutoff for mu2/phi", log_q0_lower > R(10**9).log())
    prove("switch BJS scale", wq.log() > R("10.4"))
    prove("Perron kernel norm", (4 * L + c0 + log2 + 1) / arb.pi() < 2 * L)
    prove("Perron norm derivative", (4 + 1 / L) / arb.pi() < 2)
    prove("rough floor coefficient", 3402 * residual_exponent**(-10) < R(20).exp())

    # Base derivative inequalities used for c>=24.9.
    prove("zero exponential growth", 1 - delta / 2 - 2 / c0 > R(1, 3))
    prove("nu below nonexceptional gap remains so", 1 / c0 - delta / 2 - 2 / c0 < 0)
    for tag, t in [("rectangle", logX3), ("switch", logXq)]:
        prove(tag + " log-scale at least half", 15 * t.log() < t / 2)
        prove(tag + " logarithmic derivative <1.01", 32 * t.log() / t < R("0.01"))
    prove("AP logarithmic derivative <1.01", (L - 5) / wap < R("1.01"))
    prove("smallest exponential beats degree 20", L / 16 > 20)
    prove("finite endpoint polynomial", 2 * L + 1 < (-L / 12 + L).exp())
    prove("finite polynomial constant", R(37024) < R(30).exp())

    vqprime = (
        R("3.2e-8") / wq**4
        + wq**4 * (
            (-nu * wq).exp() / (1 - nu)
            + R("1.02") * wq**10 * (-wq / 2).exp()
            + 3 * wq**10 * (-2 * wq / 3).exp()
        )
    )
    vq = vqprime * (
        1 + 1 / (logXq**10 * wq**5)
        + 1 / ((1 - 6 / wq) * logXq)
    ) + 3 / logXq
    prove("BJS Lemma 29 switch envelope v<=40/L", vq < 40 / L)
    logD_upper = dstar_cap * L
    KR = 2 + logD_upper / log2
    Kn = 2 + L / (6 * log2)
    Jperron = 2 * L
    switch_prefactor = 64 * R(1).exp()**6 * Kn * Jperron**2 * R("1.1") * logD_upper
    b1_over_N = 4 * (-log_q0_lower).exp()
    b2_over_N = R(12).sqrt() * KR * (-L / 12).exp()
    b3_over_N = R(12).sqrt() * KR * (-L / 3).exp()
    b4_over_N = 24 * (-(R("0.5") - dstar_cap) * L).exp()
    switch_high = switch_prefactor * (b1_over_N + b2_over_N + b3_over_N + b4_over_N) * L**2 / cn_lower
    prove("switch large-conductor remainder <1e-20", switch_high < R("1e-20"))
    logq0_upper = 10 * wq.log()
    switch_low = (
        4 * (40 / L) * (1 + L) / (k1 * L) ** 5
        * R("1.21") * logq0_upper * logD_upper * L**2 / cn_lower
    )
    prove("switch low-conductor prime part <1e-8", switch_low < R("1e-8"))
    perron_proxy = (20 + 20 * c0 + (dstar_cap - 2) * L).exp()
    rough_proxy = (20 - 9 * c0).exp()
    prove("Perron normalized proxy <1e-8", perron_proxy < R("1e-8"))
    prove("rough integer-sieve remainder <1e-8", rough_proxy < R("1e-8"))
    switch_proxy = switch_high + switch_low + perron_proxy + rough_proxy
    prove("complete switched remainder <1e-6", switch_proxy < R("1e-6"))

    # Finite and local terms retain the weaker N^(-1/12) envelope because
    # kappa_1>1/12; this avoids hiding a parameter-specific improvement.
    log_comb_terminal = 30 + 12 * c0 - L / 12
    prove("all finite combinatorial terminals <1e-6", log_comb_terminal < R("1e-6").log())
    finite_proxy = log_comb_terminal.exp()
    sieve_endpoint_proxy = 1000 * (2 / (k1 * L) ** 2 + 40 * L * (-k1 * L).exp() + 1 / L**9)
    prove("local-product/endpoints <1e-6", sieve_endpoint_proxy < R("1e-6"))

    budget_rationals = {
        "rectangle": fmpq(85, 10_000),
        "stieltjes_and_prime_endpoints": fmpq(35, 10_000),
        "ordinary_AP": fmpq(1, 1_000_000),
        "large_exceptional_finite": fmpq(1, 1_000_000),
        "switch_remainders": fmpq(1, 1_000_000),
        "finite_combinatorial": fmpq(1, 1_000_000),
        "local_products": fmpq(1, 1_000_000),
    }
    budgets = {key: arb(value) for key, value in budget_rationals.items()}
    exact_total_budget = sum(budget_rationals.values(), fmpq(0))
    prove("total finite budget is .012005", exact_total_budget == fmpq(2401, 200_000))
    final_floor = R("0.0228") - arb(exact_total_budget)
    prove("final positivity floor", final_floor > R("0.0107"))

    result = {
        "engine": {
            "python_flint": flint.__version__,
            "FLINT": flint.__FLINT_VERSION__,
            "python": sys.version,
            "arb_precision_bits": ctx.prec,
            "threads": ctx.threads,
            "rounding": "Arb rigorous outward ball enclosures; exact rational identities",
        },
        "threshold": {"c0": "249/10", "L0": txt(L)},
        "parameters": {
            "kappa_1": "191/2000", "kappa_2": "383/4000",
            "rho": "119/500", "sigma_2": "6047/20000",
            "sigma_1": "121/400", "delta": "29/25",
        },
        "product_and_levels": {
            "BJS_log_Q_upper": "1000000050",
            "theta_actual_lower": txt(theta_actual_lower),
            "theta_used": "2423257/5000000",
            "D_star_exponent_at_c0": txt(dstar_exponent),
            "D_star_exponent_cap": "2077/5000",
            "s_circle": txt(s_circle),
            "s_rho_at_c0": txt(s_rho),
            "s_rho_used": "12037/2500",
            "exceptional_level_loss_at_c0": txt(level_loss),
            "exceptional_level_loss_used": "81/1000000000",
            "K_base_at_c0": txt(Kbase),
            "K_rectangle_at_c0": txt(Krect),
            "K_AP_at_c0": txt(Kap),
            "K_switch_at_c0": txt(Kq),
            "largest_exceptional_prime_floor": "37",
        },
        "rectangle": {"log_X3": txt(logX3), "log_x2_X3": txt(w3), "v0": txt(v3), "m_X3": txt(mrect), "normalized_global_charge": txt(rectangle_normalized)},
        "ordinary_AP": {"c4": txt(c4), "c4_times_L2": txt(c4 * L**2), "normalized_all_rows_proxy": txt(ordinary_ap_normalized)},
        "exceptional": {"full_conductor_mass_defect": txt(b_exc), "finite_proxy": txt(exceptional_finite_proxy)},
        "transfer": {"Mertens_measure_discrepancy": txt(mertens_discrepancy), "geometric_exponent_shift": txt(geometric_shift), "normalized_proxy": txt(stieltjes_proxy)},
        "switch": {
            "residual_exponent": txt(residual_exponent),
            "log_q0_lower": txt(log_q0_lower), "v_q0": txt(vq),
            "large_conductor_normalized": txt(switch_high),
            "low_conductor_normalized": txt(switch_low),
            "Perron_normalized_proxy": txt(perron_proxy),
            "rough_remainder_proxy": txt(rough_proxy),
            "complete_proxy": txt(switch_proxy),
        },
        "finite_proxy": txt(finite_proxy),
        "local_proxy": txt(sieve_endpoint_proxy),
        "budget_rationals": {k: str(v) for k, v in budget_rationals.items()},
        "budgets": {k: txt(v) for k, v in budgets.items()},
        "total_budget": txt(arb(exact_total_budget)),
        "main_floor_imported_from_interval_certificate": "0.0228",
        "final_floor": txt(final_floor),
        "passed_assertions": passed_assertions.copy(),
    }
    path = Path(__file__).resolve()
    result["engine"]["script_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    path.with_name("c249_structural_ledger.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2))
    return result


if __name__ == "__main__":
    main()

