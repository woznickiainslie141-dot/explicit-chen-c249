"""Certify a signed sieve-remainder interface and its finite scalar budgets.

The new remainder estimate is an input, not a theorem about the primes.
Saved finite sieve inputs are checked by hash, not regenerated here.
Exact small cases exercise the signed identity, deleted primes and P2 transfer.
The manuscript contains the general arguments; these tests do not replace them.
"""

from pathlib import Path
from fractions import Fraction
import hashlib
import json
import math
from flint import arb, ctx


ROOT = Path(__file__).resolve().parent
ctx.prec = 192
GAMMA_EXP = arb.const_euler().exp()
U_STAR = arb("1.32") / GAMMA_EXP
CLAIM = arb("0.00001")
LEVELS = [
    ("EH_3_4", Fraction(3, 4)),
    ("EH_9_10", Fraction(9, 10)),
    ("EH_19_20", Fraction(19, 20)),
    ("EH_99_100", Fraction(99, 100)),
]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def ball(q):
    if isinstance(q, Fraction):
        return arb(q.numerator) / q.denominator
    return arb(q)


def scalar(case, theta, log_n, constant):
    """Two algebraically independent expressions, with strict domain checks."""
    L = arb(log_n)
    cost = arb(case["effective_log_support_cost_ball"])
    s = 3 * ball(theta) - 3 * cost / L
    if not (L > 3 * arb(10**12).log() and 2 < s < 3):
        return None
    if not L / 3 > arb(case["presieve"]["T"]).log():
        return None
    eps = arb(case["epsilon"])
    ep, em = arb(case["eta_plus"]), arb(case["eta_minus"])
    correction = (2 - s).exp()
    low = 2 * GAMMA_EXP * (s - 1).log() / s - case["C2"] * eps * correction
    upper = 2 * GAMMA_EXP / s + case["C1"] * eps * correction
    factor = (1 + ep) * low - (ep + em) * upper
    rho = 1 - arb(".000020508") / L - 2 * L * (-L / 3).exp()
    assert 0 < rho <= 1 and low <= upper
    main = 3 * U_STAR * rho * factor
    independently = rho * (
        arb("7.92") / s * ((1 + ep) * (s - 1).log() - ep - em)
        - 3 * U_STAR * eps * correction
        * ((1 + ep) * case["C2"] + (ep + em) * case["C1"])
    )
    assert abs(main - independently) < arb("1e-45")
    deleted = L**3 / arb(2).log() * (-L).exp()
    survivor_one = L**2 * (-L).exp()
    margin = main - arb(constant) / L**2 - deleted - survivor_one
    return {
        "s": s, "large_lower": low, "factor": factor, "rho": rho,
        "main": main, "deleted": deleted, "survivor_one": survivor_one,
        "margin": margin,
        "admissible_C": L**2 * (main - deleted - survivor_one - CLAIM),
    }


def endpoint_rows(saved):
    rows = []
    for key, theta in LEVELS:
        case = saved[key]
        # Reproduce the saved main mass and old absolute-remainder surplus.
        old_L = case["log_start"]
        old = scalar(case, theta, old_L, 0)
        assert abs(old["main"] - arb(case["main_ball"])) < arb("1e-45")
        L = arb(old_L)
        old_deleted = L**3 / arb(2).log() * (
            2 * (-L / 2).exp() + arb("1.1") * ball(theta) * L * (-L).exp()
        )
        budget = 20 if key == "EH_3_4" else 10
        old_margin = old["main"] - arb(1) / budget - old_deleted - old["survivor_one"]
        assert abs(old_margin - arb(case["surplus_ball"])) < arb("1e-45")
        assert 0 < old["deleted"] < old_deleted
        selected = []
        for C in [0, 1, 1000, 10000, 1000000]:
            for start in range(1, 10001):
                result = scalar(case, theta, start, C)
                if result is not None and result["margin"] > CLAIM:
                    assert 0 < result["large_lower"] and result["factor"] > 0
                    before = scalar(case, theta, start - 1, C)
                    assert before is not None and before["margin"] < CLAIM
                    if C == 0:
                        assert before["margin"] < 0
                        assert before["factor"] < 0
                    selected.append({
                        "C_upper_budget": C, "log_start": start,
                        "threshold_formula": f"max(X, exp({start}))",
                        **{name + "_ball": str(value) for name, value in result.items()},
                        "preceding_integer_margin_ball": str(before["margin"]),
                    })
                    break
            else:
                raise AssertionError("No certified endpoint in search domain")
        base_start = selected[1]["log_start"]
        capacity = scalar(case, theta, base_start, 0)["admissible_C"]
        cap_int = int(capacity.floor().unique_fmpz())
        assert arb(cap_int) < capacity < arb(cap_int + 1)
        rows.append({
            "key": key, "theta": str(theta), "fixed_parameters": {
                name: case[name] for name in ("epsilon", "eta_plus", "eta_minus", "C1", "C2")
            },
            "presieve_T": case["presieve"]["T"],
            "support_cost_ball": case["effective_log_support_cost_ball"],
            "reproduced_saved_main_and_surplus": True,
            "endpoint_cases": selected,
            "C_cap_at_first_small_budget_endpoint": cap_int,
            "C_cap_ball": str(capacity),
        })
    return rows


def prime_sieve(limit):
    flags = bytearray(b"\1") * (limit + 1)
    flags[:2] = b"\0\0"
    for p in range(2, math.isqrt(limit) + 1):
        if flags[p]:
            flags[p*p:limit+1:p] = b"\0" * ((limit - p*p) // p + 1)
    return [n for n in range(2, limit + 1) if flags[n]]


def rosser_indicator(primes, level, checked_parity):
    prefix = 1
    for m, p in enumerate(sorted(primes, reverse=True), 1):
        prefix *= p
        if m % 2 == checked_parity and prefix * p * p >= level:
            return 0
    return 1


def coefficient(chosen, w, T, D0, D):
    medium = [p for p in chosen if w < p <= T]
    large = [p for p in chosen if p > T]
    ip = rosser_indicator(medium, D0, 1)
    im = rosser_indicator(medium, D0, 0)
    jp = rosser_indicator(large, D, 1)
    jm = rosser_indicator(large, D, 0)
    mu = (-1) ** len(chosen)
    return mu * (ip * jm - (ip - im) * jp), mu * ip * jp


def exact_checks():
    primes = prime_sieve(4096)
    prime_set = set(primes)
    omega = [0] * 4097
    for p in primes:
        power = p
        while power <= 4096:
            for n in range(power, 4097, power):
                omega[n] += 1
            power *= p
    ns = list(range(32, 2050, 2)) + list(range(2050, 4097, 14))
    configs = [(1, 3, 4), (1, 5, 4), (1, 7, 5), (3, 5, 4), (3, 7, 4)]
    checked = patterns = identities = 0
    min_deleted_slack = 0
    for N in ns:
        z = 1
        while z**3 < N:
            z += 1
        pn = [p for p in primes if p < z and N % p]
        pp = math.prod(pn)
        full_p = [p for p in primes if p <= N]
        complements = [N - p for p in full_p]
        nu = sum(N % p == 0 for p in full_p)
        assert nu <= omega[N]
        A = [N-p for p in full_p if N % p]
        S_A = sum(math.gcd(n, pp) == 1 for n in A)
        survivor_one = int(N - 1 in prime_set)
        pi2 = sum(2 <= n and omega[n] <= 2 for n in complements)
        assert S_A - survivor_one <= pi2
        subsets = []
        for mask in range(1 << len(pn)):
            chosen = [p for i, p in enumerate(pn) if mask >> i & 1]
            d = math.prod(chosen)
            phi = math.prod(p - 1 for p in chosen)
            count = sum(n % d == 0 for n in complements)
            subsets.append((chosen, d, phi, count))
        for w, T, u in configs:
            if T >= z:
                continue
            for multiplier in (1, 3):
                D, D0 = multiplier * z**2, T**u
                coeffs = [coefficient(chosen, w, T, D0, D) for chosen, _, _, _ in subsets]
                assert all(abs(lo) <= 1 and abs(hi) <= 1 for lo, hi in coeffs)
                values = [lo for lo, _ in coeffs]
                uppers = [hi for _, hi in coeffs]
                for bit in range(len(pn)):
                    for mask in range(len(values)):
                        if mask >> bit & 1:
                            values[mask] += values[mask ^ (1 << bit)]
                            uppers[mask] += uppers[mask ^ (1 << bit)]
                for mask, (lo, hi) in enumerate(zip(values, uppers)):
                    indicator = int(mask == 0)
                    assert lo <= indicator <= hi and hi >= 0
                    patterns += 1
                G = sum((Fraction(coeffs[i][0], phi) for i, (_, _, phi, _) in enumerate(subsets)), Fraction())
                weighted_ap = sum(coeffs[i][0] * count for i, (_, _, _, count) in enumerate(subsets))
                R = sum((coeffs[i][0] * (Fraction(count) - Fraction(len(full_p), phi))
                         for i, (_, _, phi, count) in enumerate(subsets)), Fraction())
                assert Fraction(weighted_ap) == len(full_p) * G + R
                weighted_full = weighted_deleted = 0
                for p in full_p:
                    n = N - p
                    mask = sum(1 << i for i, q in enumerate(pn) if n % q == 0)
                    weighted_full += values[mask]
                    if N % p == 0:
                        weighted_deleted += values[mask]
                assert weighted_ap == weighted_full
                assert weighted_deleted <= nu
                assert Fraction(S_A) >= len(full_p) * G + R - nu
                min_deleted_slack = min(min_deleted_slack, weighted_deleted - nu)
                identities += 2
                checked += 1
    return {
        "even_N_count": len(ns), "N_max": max(ns), "composite_cases": checked,
        "pointwise_prime_divisor_patterns": patterns,
        "exact_signed_identity_checks": identities,
        "arithmetic": "integers and Fraction; complete prime sieve",
        "minimum_weighted_deleted_minus_nu": min_deleted_slack,
        "P2_transfer_checked_for_all_selected_N": True,
    }


def factorability_checks(saved):
    p, z, T, D0, D = 999983, 1000000, 5, 625, 10**14
    assert all(p % n for n in range(2, math.isqrt(p) + 1))
    N = z**3
    assert N % p and T < p < z and D >= z**2
    coefficient_p = coefficient([p], 3, T, D0, D)[0]
    H = 3 * D0 * D
    assert coefficient_p == -1 and p**3 > H
    families = []
    for key, theta in LEVELS:
        # For N=2^(3000k), H^(1/3)=2^(1000 theta k).
        # Bertrand supplies a prime in (H^(1/3),2 H^(1/3)).
        exponent = 1000 * theta
        assert exponent.denominator == 1
        e = exponent.numerator
        assert 2**e > saved[key]["presieve"]["T"] and e + 1 < 1000
        slack = 3000 * arb(2).log() * (ball(theta) - ball(Fraction(2, 3)))
        slack -= arb(saved[key]["effective_log_support_cost_ball"])
        assert slack > 0
        families.append({
            "theta": str(theta), "N_family": "2^(3000k), integers k>=1",
            "balanced_factor_exponent_at_k_1": e,
            "Bertrand_prime_interval_exponents_at_k_1": [e, e + 1],
            "z_exponent_at_k_1": 1000,
            "log_D_over_z_squared_at_k_1_ball": str(slack),
            "prime_lower_coefficient": -1,
            "triply_well_factorable_at_N_to_theta": False,
        })
    mpz_delta = Fraction(1, 3)
    assert 180 * mpz_delta == 60 > 7
    assert Fraction(1, 2) + 2 * Fraction(7, 600) == Fraction(157, 300) < Fraction(2, 3)
    return {
        "exact_prime_witness": {"N": N, "z": z, "p": p, "D": D, "D0": D0,
                                "level_H": H, "p_cubed": p**3, "coefficient": coefficient_p},
        "current_parameter_families": families,
        "MPZ_published_region": "600 varpi + 180 delta < 7",
        "delta_to_cover_all_current_prime_coefficients_on_a_half_line": "at least 1/3",
        "published_MPZ_level_supremum": "157/300 < 2/3",
        "two_factor_well_factorability_not_decided": True,
    }


def run():
    if not __debug__:
        raise RuntimeError("Run without Python -O.")
    saved = json.loads((ROOT / "conditional_chen_sharp_endpoint_certificate.json").read_text(encoding="utf-8"))
    assert saved["precision_bits"] == 192
    assert saved["finite_inputs_regenerated_in_this_run"] is False
    assert saved["support_maxima_regenerated_in_this_run"] is True
    for name, expected in saved["input_sha256"].items():
        assert digest(ROOT / name) == expected, "Stale input: " + name
    result = {
        "status": "PASS: exact finite interfaces, factorability witnesses and strict scalar budgets",
        "precision_bits": ctx.prec, "claimed_normalized_lower_count": "0.00001",
        "signed_input": "sum_d lambda_N^-(d) [pi(N;d,N)-pi(N)/phi(d)] >= -C N/log^4(N) for every even N>=X in the sieve domain",
        "signed_input_domain": "L>3 log(10^12), T<N^(1/3), 2<s<3; coefficients used only in this domain",
        "C_and_X_are_unspecified_distribution_inputs": True,
        "finite_prime_product_and_gap_inputs_regenerated": False,
        "endpoint_rows": endpoint_rows(saved),
        "exact_small_cases": exact_checks(),
        "support_class_checks": factorability_checks(saved),
        "input_sha256": {
            **saved["input_sha256"],
            "conditional_chen_sharp_endpoint_certificate.json": digest(ROOT / "conditional_chen_sharp_endpoint_certificate.json"),
            Path(__file__).name: digest(Path(__file__)),
        },
        "not_checked": [
            "truth or effective constants/onset of the new signed distribution input",
            "general analytic sieve lemmas and independent mathematical review",
            "applicability of other well-factorable linear-sieve variants",
            "global theoretical optimality or a full splice with finite verification",
            "PDF compilation and layout",
        ],
    }
    return result


if __name__ == "__main__":
    result = run()
    Path(__file__).with_suffix(".json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({
        "status": result["status"], "exact_small_cases": result["exact_small_cases"],
        "endpoint_table_C_1_1000_10000_1000000": [
            [row["theta"], *[case["log_start"] for case in row["endpoint_cases"] if case["C_upper_budget"]]]
            for row in result["endpoint_rows"]
        ],
    }, ensure_ascii=False))
