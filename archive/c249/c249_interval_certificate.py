#!/usr/bin/env python3
"""Arb interval certificate for the c0=24.9 positivity ledger.

The elementary delay-sieve kernels are reused from the independently
audited c=26 certificate.  This file replaces every parameter-dependent
row, adds the new s=2 splits, and uses the q1>=37 exceptional envelope.
All integrations are ACB ball integrations at 192-bit precision.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import flint
from flint import acb, arb, ctx

import c26_interval_certificate as kernel


ctx.prec = 192
ctx.threads = 1
A = kernel.A

eta = kernel.eta
C1 = kernel.C1
C2 = kernel.C2
egamma = kernel.egamma
emgamma = kernel.emgamma
e2 = kernel.e2

k1 = A(191, 2000)
k2 = A(383, 4000)
sigma2 = A(6047, 20000)
sigma1 = A(121, 400)

# Rational envelopes proved by c249_structural_certificate.py.
theta = A(2_423_257, 5_000_000)
level_shift = A(81, 1_000_000_000)
s_rho = A(12_037, 2_500)
s_circle = A(800, 191)
exc_eps = A(1, 35)
exc_density = A(36, 35)


def integral(name: str, func, lo: acb, hi: acb) -> acb:
    return kernel.integral(name, func, lo, hi)


def fexceptional_le2(s: acb) -> acb:
    # This branch is used only through the trivial inequality S>=0.
    del s
    corr = eta * C2 * e2 * kernel.h_le2(A(1))
    return -4 * exc_eps - (1 + 4 * exc_eps) * corr


def fexceptional_2_3(s: acb, analytic: bool = False) -> acb:
    fv = kernel.f_2_4(s, analytic)
    Fv = kernel.F_le3(s)
    corr2 = eta * C2 * e2 * kernel.h_2_3(s)
    corr1 = eta * C1 * e2 * kernel.h_2_3(s)
    fminus = fv - corr2
    Fplus = Fv + corr1
    return fminus - 3 * exc_eps * (Fplus - fminus)


def fexceptional_3_4(s: acb, analytic: bool = False) -> acb:
    fv = kernel.f_2_4(s, analytic)
    Fv = kernel.F_3_5(s, analytic)
    corr2 = eta * C2 * e2 * kernel.h_ge3(s)
    corr1 = eta * C1 * e2 * kernel.h_ge3(s)
    fminus = fv - corr2
    Fplus = Fv + corr1
    return fminus - 3 * exc_eps * (Fplus - fminus)


def fexceptional_scalar(name: str, s: acb) -> acb:
    u = s - level_shift
    fv = kernel.f_4_6_scalar(name + ".f", u)
    Fv = kernel.F_5_7_scalar(name + ".F", u)
    corr2 = eta * C2 * e2 * kernel.h_ge3(u)
    corr1 = eta * C1 * e2 * kernel.h_ge3(u)
    fminus = fv - corr2
    Fplus = Fv + corr1
    return fminus - 3 * exc_eps * (Fplus - fminus)


def ordinary_rows(exceptional: bool) -> dict[str, acb]:
    th = theta - (k1 * level_shift if exceptional else 0)
    if exceptional:
        F0 = A(2) * fexceptional_scalar("F0.exc", theta / k1) / (k1 * egamma)
    else:
        s0 = theta / k1
        F0 = A(2) * (
            kernel.f_4_6_scalar("F0.small", s0)
            - eta * C2 * e2 * kernel.h_ge3(s0)
        ) / (k1 * egamma)

    F1 = A(2) / egamma * integral(
        "F1.exc" if exceptional else "F1.small",
        lambda t, analytic: kernel.Fplus_3_5((th - t) / t, analytic) / t**2,
        k1,
        k2,
    )

    split3 = th - 3 * k1
    split2 = th - 2 * k1

    def upper_fixed(name: str, upper: acb) -> acb:
        pieces = [
            integral(
                name + ".s_ge_3",
                lambda t, analytic: kernel.Fplus_3_5((th - t) / k1, analytic) / t,
                k1,
                split3,
            ),
            integral(
                name + ".s_2_3",
                lambda t, analytic: kernel.Fplus_2_3((th - t) / k1) / t,
                split3,
                split2,
            ),
            integral(
                name + ".s_1_2",
                lambda t, analytic: kernel.Fplus_le3((th - t) / k1) / t,
                split2,
                upper,
            ),
        ]
        return A(2) * sum(pieces, A(0)) / (k1 * egamma)

    F2 = upper_fixed("F2.exc" if exceptional else "F2.small", sigma1)
    F3 = upper_fixed("F3.exc" if exceptional else "F3.small", sigma2)

    def lower_F4(name: str) -> acb:
        def lower(s: acb, analytic: bool) -> acb:
            arg = (theta - s) / k1
            if exceptional:
                return fexceptional_3_4(arg - level_shift, analytic)
            return kernel.fminus_3_4(arg, analytic)

        left = integral(
            name + ".sum_left",
            lambda s, analytic: lower(s, analytic)
            * ((s - k1) / k1).log(analytic=analytic) / s,
            2 * k1,
            k1 + k2,
        )
        right = integral(
            name + ".sum_right",
            lambda s, analytic: lower(s, analytic)
            * (k2 / (s - k2)).log(analytic=analytic) / s,
            k1 + k2,
            2 * k2,
        )
        return A(2) * (left + right) / (k1 * egamma)

    def lower_F5(name: str) -> acb:
        ctop = A(1, 2) - 2 * k1
        lower_theta = theta - (k1 * level_shift if exceptional else 0)
        s3 = lower_theta - 3 * k1
        s2 = lower_theta - 2 * k1

        def weight_left(s: acb, analytic: bool) -> acb:
            return ((s - k2) * (s - k1) / (k1 * k2)).log(analytic=analytic) / s

        def weight_right(s: acb, analytic: bool) -> acb:
            return (k2 * (s - k1) / (k1 * (s - k2))).log(analytic=analytic) / s

        def lower3(s: acb, analytic: bool) -> acb:
            arg = (theta - s) / k1
            return (
                fexceptional_3_4(arg - level_shift, analytic)
                if exceptional else kernel.fminus_3_4(arg, analytic)
            )

        def lower2(s: acb, analytic: bool) -> acb:
            arg = (theta - s) / k1
            return (
                fexceptional_2_3(arg - level_shift, analytic)
                if exceptional else kernel.fminus_2_3(arg, analytic)
            )

        def lower1(s: acb, analytic: bool) -> acb:
            del analytic
            arg = (theta - s) / k1
            return (
                fexceptional_le2(arg - level_shift)
                if exceptional else kernel.fminus_le2(arg)
            )

        pieces = [
            integral(name + ".sum_left", lambda s, a: lower3(s, a) * weight_left(s, a), k1 + k2, 2 * k2),
            integral(name + ".sum_3", lambda s, a: lower3(s, a) * weight_right(s, a), 2 * k2, s3),
            integral(name + ".sum_2", lambda s, a: lower2(s, a) * weight_right(s, a), s3, s2),
            integral(name + ".sum_1", lambda s, a: lower1(s, a) * weight_right(s, a), s2, ctop),
        ]
        return A(2) * sum(pieces, A(0)) / (k1 * egamma)

    return {
        "F0": F0,
        "F1": F1,
        "F2": F2,
        "F3": F3,
        "F4": lower_F4("F4.exc" if exceptional else "F4.small"),
        "F5": lower_F5("F5.exc" if exceptional else "F5.small"),
    }


def box_rows(exceptional: bool, js: dict[str, acb]) -> dict[str, acb]:
    if exceptional:
        factor = exc_density * 6 * emgamma * kernel.Fplus_le3(3 * theta - level_shift)
    else:
        factor = 6 * emgamma * kernel.Fplus_le3(3 * theta)
    return {"F" + key[1:]: factor * value for key, value in js.items()}


def I_rows() -> dict[str, acb]:
    ratio_log = (k2 / k1).log()
    I13 = ratio_log**4 / 24
    lam = A(1, 2) - 2 * k1
    I14 = integral(
        "I14",
        lambda t, analytic: (
            (t / k1).log(analytic=analytic) ** 2 / (2 * t)
            * ((lam - t) / k2).log(analytic=analytic)
        ),
        k1,
        k2,
    )
    return {"I13": I13, "I14": I14}


def switch_row(exceptional: bool, Is: dict[str, acb]) -> acb:
    rough = emgamma * kernel.Fplus_3_5(s_rho)
    circle = kernel.Fplus_3_5(s_circle - (level_shift if exceptional else 0))
    density = exc_density if exceptional else A(1)
    return (
        density * (1 + eta) * 2 * emgamma * rough * circle
        / k1**2 * (Is["I13"] + Is["I14"])
    )


def ball_text(value: acb | arb, digits: int = 32) -> str:
    return kernel.ball_text(value, digits)


def main() -> dict:
    kernel.integral_records.clear()
    # The shared J/I routines read their parameters from the kernel module.
    old = (kernel.k1, kernel.k2, kernel.sigma2, kernel.sigma1)
    kernel.k1, kernel.k2, kernel.sigma2, kernel.sigma1 = k1, k2, sigma2, sigma1
    try:
        small = ordinary_rows(False)
        large = ordinary_rows(True)
        js = kernel.J_rows()
        Is = I_rows()
    finally:
        kernel.k1, kernel.k2, kernel.sigma2, kernel.sigma1 = old

    small_box = box_rows(False, js)
    large_box = box_rows(True, js)
    small_switch = switch_row(False, Is)
    large_switch = switch_row(True, Is)

    small_master = (
        4 * small["F0"] - small["F1"] - small["F2"] - small["F3"]
        + small["F4"] + small["F5"] - 2 * small_box["F7"]
        - 2 * small_box["F8"] - small_box["F9"] - small_box["F10"]
        - small_switch
    ) / 4
    large_master = (
        4 * large["F0"] - exc_density * (large["F1"] + large["F2"] + large["F3"])
        + large["F4"] + large["F5"] - 2 * large_box["F7"]
        - 2 * large_box["F8"] - large_box["F9"] - large_box["F10"]
        - large_switch
    ) / 4
    negative_main = (
        exc_density * (large["F1"] + large["F2"] + large["F3"])
        + 2 * large_box["F7"] + 2 * large_box["F8"]
        + large_box["F9"] + large_box["F10"] + large_switch
    )
    small_real = kernel.require_real_finite("small master", small_master)
    large_real = kernel.require_real_finite("large master", large_master)
    if not small_real > arb("0.35"):
        raise ArithmeticError(f"small master too small: {small_real}")
    if not large_real > arb("0.0228"):
        raise ArithmeticError(f"large master too small: {large_real}")
    if not kernel.require_real_finite("negative main", negative_main) < arb(60):
        raise ArithmeticError("negative main exceeds transfer envelope")

    result = {
        "engine": {
            "python_flint": flint.__version__,
            "FLINT": flint.__FLINT_VERSION__,
            "python": sys.version,
            "arb_precision_bits": ctx.prec,
            "threads": ctx.threads,
            "rounding": "Arb/ACB rigorous outward ball enclosures",
            "shared_kernel": "c26_interval_certificate.py",
        },
        "rational_parameters": {
            "kappa_1": "191/2000",
            "kappa_2": "383/4000",
            "sigma_2": "6047/20000",
            "sigma_1": "121/400",
            "theta": "2423257/5000000",
            "exceptional_level_shift": "81/1000000000",
            "exceptional_epsilon": "1/35",
            "exceptional_density": "36/35",
            "s_rho": "12037/2500",
            "s_circle": "800/191",
        },
        "ordinary_small": {k: ball_text(v) for k, v in small.items()},
        "ordinary_large_lower": {k: ball_text(v) for k, v in large.items() if k in {"F0", "F4", "F5"}},
        "ordinary_large_upper_base": {k: ball_text(v) for k, v in large.items() if k in {"F1", "F2", "F3"}},
        "J": {k: ball_text(v) for k, v in js.items()},
        "box_small": {k: ball_text(v) for k, v in small_box.items()},
        "box_large": {k: ball_text(v) for k, v in large_box.items()},
        "I": {k: ball_text(v) for k, v in Is.items()},
        "switch_small": ball_text(small_switch),
        "switch_large": ball_text(large_switch),
        "master_before_finite_charges": {
            "small": ball_text(small_master),
            "large": ball_text(large_master),
        },
        "assertions": {"small_master_gt": "0.35", "large_master_gt": "0.0228", "negative_main_lt": "60"},
        "integrals": kernel.integral_records.copy(),
    }
    path = Path(__file__).resolve()
    result["engine"]["script_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    result["engine"]["shared_kernel_sha256"] = hashlib.sha256(
        path.with_name("c26_interval_certificate.py").read_bytes()
    ).hexdigest()
    path.with_name("c249_interval_ledger.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2))
    return result


if __name__ == "__main__":
    main()

