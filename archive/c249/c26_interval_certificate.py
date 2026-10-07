#!/usr/bin/env python3
"""Arb certificate for the c=26 positivity ledger.

All transcendental evaluations and integrations are ball computations made
by Arb through python-flint.  The script deliberately evaluates the main
terms at rational worst-case envelopes rather than at rounded values of
log log N.  It aborts if an integral is non-real, non-finite, or has less
than 90 bits of relative accuracy.

This certificate covers the finite main integrals.  Elementary support,
rectangle, large-sieve, and terminal inequalities are checked in the
companion Arb certificate c26_structural_certificate.py.  An arithmetic
certificate does not by itself discharge the analytic proof obligations.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import flint
from flint import acb, arb, ctx


ctx.prec = 192
ctx.threads = 1
integral_records: list[dict] = []


def A(num: int | str, den: int | None = None) -> acb:
    if den is None:
        return acb(num)
    return acb(num) / den


def require_real_finite(name: str, value: acb, bits: int = 90) -> arb:
    if not value.is_finite() or not value.imag.contains(0):
        raise ArithmeticError(f"{name} is not a finite real ball: {value}")
    out = value.real
    if out.rel_accuracy_bits() < bits:
        raise ArithmeticError(
            f"{name} has only {out.rel_accuracy_bits()} accurate bits: {out}"
        )
    return out


def integral(name: str, func, lo: acb, hi: acb) -> acb:
    value = acb.integral(func, lo, hi, eval_limit=200000, depth_limit=60)
    real = require_real_finite(name, value)
    integral_records.append({
        "name": name,
        "lower_endpoint": lo.real.str(32, more=True),
        "upper_endpoint": hi.real.str(32, more=True),
        "value": real.str(32, more=True),
        "relative_accuracy_bits": real.rel_accuracy_bits(),
    })
    return value


gamma = acb(arb.const_euler())
egamma = gamma.exp()
emgamma = (-gamma).exp()
e2 = A(2).exp()
pi = acb(arb.pi())

eta = A(1452, 10_000_000_000)
C1 = A(106)
C2 = A(107)
Cbar = A(107)

k1 = A(1, 12)
k2 = A(29, 250)
sigma2 = A(141, 500)
sigma1 = A(41, 125)

# Exact rational envelopes proved by c26_structural_certificate.py.
theta = A(4_948_909, 10_000_000)       # theta(N) >= 0.4948909
level_shift = A(1, 30_000_000)         # all exceptional level shifts < this
s_rho = A(25_773, 5_000)               # residual rough sieve parameter
s_circle_inf = A(24, 5)                # D_circle=N^(2/5), for every N


def h_le2(s: acb) -> acb:
    del s
    return A(-2).exp()


def h_2_3(s: acb) -> acb:
    return (-s).exp()


def h_ge3(s: acb) -> acb:
    return A(3) * (-s).exp() / s


def F_le3(s: acb) -> acb:
    return A(2) * egamma / s


def f_2_4(s: acb, analytic: bool = False) -> acb:
    return A(2) * egamma * (s - 1).log(analytic=analytic) / s


def F_3_5(s: acb, analytic: bool = False) -> acb:
    """Closed form obtained from (sF(s))'=f(s-1), valid on 3<s<5."""
    x = s - 2
    primitive = (
        x.log(analytic=analytic) * (1 + x).log(analytic=analytic)
        + (-x).polylog(2)
        + pi**2 / 12
    )
    return A(2) * egamma * (1 + primitive) / s


def f_4_6_scalar(name: str, s: acb) -> acb:
    """Rigorous scalar evaluation of f on 4<s<6."""
    upper = s - 1
    tail = integral(
        name + ".tail",
        lambda u, analytic: F_3_5(u, analytic),
        A(3),
        upper,
    )
    return (A(2) * egamma * A(3).log() + tail) / s


def F_5_7_scalar(name: str, s: acb) -> acb:
    """Rigorous scalar evaluation of F on 5<s<7, reduced to one integral."""
    a = s - 1
    base = A(5) * F_3_5(A(5))
    c0 = A(2) * egamma * A(3).log()
    first = c0 * (a / 4).log()
    second = integral(
        name + ".tail",
        lambda u, analytic: F_3_5(u, analytic)
        * (a / (u + 1)).log(analytic=analytic),
        A(3),
        a - 1,
    )
    return (base + first + second) / s


def Fplus_le3(s: acb) -> acb:
    # All calls to this branch in the ledger have 1<s<2.
    return F_le3(s) + eta * C1 * e2 * h_le2(s)


def Fplus_2_3(s: acb) -> acb:
    return F_le3(s) + eta * C1 * e2 * h_2_3(s)


def Fplus_3_5(s: acb, analytic: bool = False) -> acb:
    return F_3_5(s, analytic) + eta * C1 * e2 * h_ge3(s)


def Fplus_5_7_scalar(name: str, s: acb) -> acb:
    return F_5_7_scalar(name, s) + eta * C1 * e2 * h_ge3(s)


def fminus_3_4(s: acb, analytic: bool = False) -> acb:
    return f_2_4(s, analytic) - eta * C2 * e2 * h_ge3(s)


def fminus_2_3(s: acb, analytic: bool = False) -> acb:
    return f_2_4(s, analytic) - eta * C2 * e2 * h_2_3(s)


def fminus_le2(s: acb) -> acb:
    return -eta * C2 * e2 * h_le2(s)


def fexceptional_le2(s: acb) -> acb:
    """Large-exceptional lower envelope where f=0 and mbar<=1."""
    eps1 = A(1, 29)
    corr2 = eta * C2 * e2 * h_le2(s)
    corrb = eta * Cbar * e2 * h_le2(s)
    return -eps1 - (1 + eps1) * corr2 - 3 * eps1 * (1 + corrb)


def fexceptional_2_3(s: acb, analytic: bool = False) -> acb:
    eps1 = A(1, 29)
    fv = f_2_4(s, analytic)
    Fv = F_le3(s)
    # max(1-f,F-1) <= (1-f)+(F-1)=F-f on this range.
    mbar = Fv - fv
    corr2 = eta * C2 * e2 * h_2_3(s)
    corrb = eta * Cbar * e2 * h_2_3(s)
    return fv - eps1 * (1 - fv) - (1 + eps1) * corr2 \
        - 3 * eps1 * (mbar + corrb)


def fexceptional_3_4(s: acb, analytic: bool = False) -> acb:
    eps1 = A(1, 29)
    fv = f_2_4(s, analytic)
    Fv = F_3_5(s, analytic)
    mbar = Fv - fv
    corr2 = eta * C2 * e2 * h_ge3(s)
    corrb = eta * Cbar * e2 * h_ge3(s)
    return fv - eps1 * (1 - fv) - (1 + eps1) * corr2 \
        - 3 * eps1 * (mbar + corrb)


def fexceptional_scalar(name: str, s: acb) -> acb:
    """Large-exceptional envelope for the scalar F0 argument (5<s<6)."""
    eps1 = A(1, 29)
    u = s - level_shift
    fv = f_4_6_scalar(name + ".f", u)
    Fv = F_5_7_scalar(name + ".F", u)
    mbar = Fv - fv
    corr2 = eta * C2 * e2 * h_ge3(u)
    corrb = eta * Cbar * e2 * h_ge3(u)
    return fv - eps1 * (1 - fv) - (1 + eps1) * corr2 \
        - 3 * eps1 * (mbar + corrb)


def ordinary_rows(exceptional: bool) -> dict[str, acb]:
    th = theta - (k1 * level_shift if exceptional else 0)

    if exceptional:
        F0 = A(2) * fexceptional_scalar("F0.exc", theta / k1) / (k1 * egamma)
    else:
        s0 = theta / k1
        F0 = A(2) * (
            f_4_6_scalar("F0.small", s0) - eta * C2 * e2 * h_ge3(s0)
        ) / (k1 * egamma)

    # F1 remains wholly in 3<s<5.
    F1 = A(2) / egamma * integral(
        "F1.exc" if exceptional else "F1.small",
        lambda t, analytic: Fplus_3_5((th - t) / t, analytic) / t**2,
        k1,
        k2,
    )

    # F2 and F3 cross s=3 at t=theta-3*k1.
    split = th - 3 * k1

    def upper_fixed(name: str, upper: acb) -> acb:
        first = integral(
            name + ".s_ge_3",
            lambda t, analytic: Fplus_3_5((th - t) / k1, analytic) / t,
            k1,
            split,
        )
        second = integral(
            name + ".s_le_3",
            lambda t, analytic: Fplus_2_3((th - t) / k1) / t,
            split,
            upper,
        )
        return A(2) * (first + second) / (k1 * egamma)

    F2 = upper_fixed("F2.exc" if exceptional else "F2.small", sigma1)
    F3 = upper_fixed("F3.exc" if exceptional else "F3.small", sigma2)

    def lower_F4(name: str) -> acb:
        # The original triangle k1 <= t <= u <= k2 has symmetric
        # integrand g(t+u)/(tu).  Reflect it to the square and integrate
        # first along t+u=s.  This removes nested adaptive quadrature.
        def lower(s: acb, analytic: bool) -> acb:
            arg = (theta - s) / k1
            if exceptional:
                return fexceptional_3_4(arg - level_shift, analytic)
            return fminus_3_4(arg, analytic)

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
        # Here k1<=t<=k2, k2<=u and t+u<=1/2-2k1.  On a
        # line t+u=s the t-integral is elementary.  Split both where the
        # geometric weight changes (s=2k2) and where the sieve argument
        # crosses 3 and 2.
        ctop = A(1, 2) - 2 * k1
        lower_theta = theta - (k1 * level_shift if exceptional else 0)
        s3 = lower_theta - 3 * k1
        s2 = lower_theta - 2 * k1

        def weight_left(s: acb, analytic: bool) -> acb:
            return (
                ((s - k2) * (s - k1) / (k1 * k2)).log(analytic=analytic)
                / s
            )

        def weight_right(s: acb, analytic: bool) -> acb:
            return (
                (k2 * (s - k1) / (k1 * (s - k2))).log(analytic=analytic)
                / s
            )

        def lower3(s: acb, analytic: bool) -> acb:
            arg = (theta - s) / k1
            if exceptional:
                return fexceptional_3_4(arg - level_shift, analytic)
            return fminus_3_4(arg, analytic)

        def lower2(s: acb, analytic: bool) -> acb:
            arg = (theta - s) / k1
            if exceptional:
                return fexceptional_2_3(arg - level_shift, analytic)
            return fminus_2_3(arg, analytic)

        def lower1(s: acb, analytic: bool) -> acb:
            del analytic
            arg = (theta - s) / k1
            if exceptional:
                return fexceptional_le2(arg - level_shift)
            return fminus_le2(arg)

        pieces = [
            integral(
                name + ".sum_left",
                lambda s, analytic: lower3(s, analytic) * weight_left(s, analytic),
                k1 + k2,
                2 * k2,
            ),
            integral(
                name + ".sum_3",
                lambda s, analytic: lower3(s, analytic) * weight_right(s, analytic),
                2 * k2,
                s3,
            ),
            integral(
                name + ".sum_2",
                lambda s, analytic: lower2(s, analytic) * weight_right(s, analytic),
                s3,
                s2,
            ),
            integral(
                name + ".sum_1",
                lambda s, analytic: lower1(s, analytic) * weight_right(s, analytic),
                s2,
                ctop,
            ),
        ]
        return A(2) * sum(pieces, A(0)) / (k1 * egamma)

    F4 = lower_F4("F4.exc" if exceptional else "F4.small")
    F5 = lower_F5("F5.exc" if exceptional else "F5.small")
    return {"F0": F0, "F1": F1, "F2": F2, "F3": F3, "F4": F4, "F5": F5}


def J_rows() -> dict[str, acb]:
    def j78(name: str, sig: acb) -> acb:
        return integral(
            name,
            lambda x, analytic: (x - 1).log(analytic=analytic) / x,
            A(2),
            1 / sig - 1,
        )

    def j910(name: str, lo: acb, sig: acb) -> acb:
        return integral(
            name,
            lambda t, analytic: (
                (1 / sig - 1 - t / sig).log(analytic=analytic)
                / (t * (1 - t))
            ),
            lo,
            sig,
        )

    return {
        "J7": j78("J7", sigma1),
        "J8": j78("J8", sigma2),
        "J9": j910("J9", k1, sigma1),
        "J10": j910("J10", k2, sigma2),
    }


def I_rows() -> dict[str, acb]:
    R = (k2 / k1).log()
    I13 = R**4 / 24
    I14 = integral(
        "I14",
        lambda t, analytic: (
            (t / k1).log(analytic=analytic) ** 2
            / (2 * t)
            * ((A(1, 3) - t) / k2).log(analytic=analytic)
        ),
        k1,
        k2,
    )
    return {"I13": I13, "I14": I14}


def box_rows(exceptional: bool, js: dict[str, acb]) -> dict[str, acb]:
    if exceptional:
        s = 3 * theta - level_shift
        factor = A(30, 29) * 6 * emgamma * Fplus_le3(s)
    else:
        factor = 6 * emgamma * Fplus_le3(3 * theta)
    return {"F" + key[1:]: factor * value for key, value in js.items()}


def switch_row(exceptional: bool, Is: dict[str, acb]) -> acb:
    rough = emgamma * Fplus_5_7_scalar("switch.rough", s_rho)
    sc = s_circle_inf - (level_shift if exceptional else 0)
    circle = Fplus_3_5(sc)
    density = A(30, 29) if exceptional else A(1)
    return (
        density * (1 + eta) * 2 * emgamma * rough * circle
        / k1**2 * (Is["I13"] + Is["I14"])
    )


def ball_text(value: acb | arb, digits: int = 32) -> str:
    if isinstance(value, acb):
        value = require_real_finite("output", value)
    return value.str(digits, more=True)


def main() -> dict:
    integral_records.clear()
    small = ordinary_rows(False)
    large = ordinary_rows(True)
    js = J_rows()
    small_box = box_rows(False, js)
    large_box = box_rows(True, js)
    Is = I_rows()
    small_switch = switch_row(False, Is)
    large_switch = switch_row(True, Is)

    small_master = (
        4 * small["F0"] - small["F1"] - small["F2"] - small["F3"]
        + small["F4"] + small["F5"]
        - 2 * small_box["F7"] - 2 * small_box["F8"]
        - small_box["F9"] - small_box["F10"] - small_switch
    ) / 4
    large_master = (
        4 * large["F0"]
        - A(30, 29) * (large["F1"] + large["F2"] + large["F3"])
        + large["F4"] + large["F5"]
        - 2 * large_box["F7"] - 2 * large_box["F8"]
        - large_box["F9"] - large_box["F10"] - large_switch
    ) / 4

    negative_main = (
        A(30, 29) * (large["F1"] + large["F2"] + large["F3"])
        + 2 * large_box["F7"] + 2 * large_box["F8"]
        + large_box["F9"] + large_box["F10"] + large_switch
    )
    if not require_real_finite("negative-main inflation base", negative_main) < arb(60):
        raise ArithmeticError("the prime-count inflation base exceeds 60")

    small_master_real = require_real_finite("small master", small_master)
    large_master_real = require_real_finite("large master", large_master)
    if not small_master_real > arb("0.6435"):
        raise ArithmeticError(f"small master does not exceed 0.6435: {small_master_real}")
    if not large_master_real > arb("0.1105"):
        raise ArithmeticError(f"large master does not exceed 0.1105: {large_master_real}")

    result = {
        "engine": {
            "python_flint": flint.__version__,
            "FLINT": flint.__FLINT_VERSION__,
            "python": sys.version,
            "arb_precision_bits": ctx.prec,
            "threads": ctx.threads,
            "rounding": "Arb/ACB rigorous outward ball enclosures",
        },
        "rational_envelopes": {
            "theta": "4948909/10000000",
            "exceptional_level_shift": "1/30000000",
            "s_rho": "25773/5000",
            "s_circle_inf": "24/5",
        },
        "ordinary_small": {k: ball_text(v) for k, v in small.items()},
        "ordinary_large_lower": {
            k: ball_text(v) for k, v in large.items() if k in {"F0", "F4", "F5"}
        },
        "ordinary_large_upper_base": {
            k: ball_text(v) for k, v in large.items() if k in {"F1", "F2", "F3"}
        },
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
        "assertions": {
            "small_master_gt": "0.6435",
            "large_master_gt": "0.1105",
            "negative_main_lt": "60",
        },
        "integrals": integral_records.copy(),
    }

    script_path = Path(__file__).resolve()
    result["engine"]["script_sha256"] = hashlib.sha256(
        script_path.read_bytes()
    ).hexdigest()
    output_path = script_path.with_name("c26_interval_ledger.json")
    output_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return result


if __name__ == "__main__":
    main()

