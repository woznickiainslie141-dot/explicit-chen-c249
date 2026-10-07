#!/usr/bin/env python3
"""Non-rigorous parameter reconnaissance for the c=24.9 repair.

This script is intentionally *not* a proof certificate.  It uses SciPy
binary floating point only to choose rational parameters before the Arb
certificates are rebuilt.  Unlike the earlier sensitivity scratch work it
uses Wu's parameter-dependent switch geometry

    R >= N**(1/2 + 2*kappa_1 - 2*kappa_2)

and keeps Wu's exact lambda=1/2-2*kappa_1 endpoint.
"""

from __future__ import annotations

import argparse
import math

import numpy as np
from scipy.integrate import quad
from scipy.optimize import differential_evolution
from scipy.special import spence


GAMMA = 0.577215664901532860606512090082402431
EGAMMA = math.exp(GAMMA)
ETA = 1.452e-7
QBAR = 1_000_000_050.0
EXC_EPS = 1.0 / 35.0
EXC_DENSITY = 36.0 / 35.0


def h(s: float) -> float:
    if s <= 2:
        return math.exp(-2)
    if s <= 3:
        return math.exp(-s)
    return 3 * math.exp(-s) / s


def sieve_F(s: float) -> float:
    if s <= 3:
        return 2 * EGAMMA / s
    if s <= 5:
        x = s - 2
        primitive = math.log(x) * math.log1p(x) + spence(1 + x) + math.pi**2 / 12
        return 2 * EGAMMA * (1 + primitive) / s
    if s <= 7:
        base = 5 * sieve_F(5.0)
        tail = quad(sieve_f, 4.0, s - 1.0, epsabs=2e-11, epsrel=2e-11)[0]
        return (base + tail) / s
    raise ValueError(f"F branch exceeded: {s}")


def sieve_f(s: float) -> float:
    if s <= 2:
        return 0.0
    if s <= 4:
        return 2 * EGAMMA * math.log(s - 1) / s
    if s <= 6:
        tail = quad(sieve_F, 3.0, s - 1.0, epsabs=2e-11, epsrel=2e-11)[0]
        return (2 * EGAMMA * math.log(3) + tail) / s
    raise ValueError(f"f branch exceeded: {s}")


def Fplus(s: float) -> float:
    return sieve_F(s) + ETA * 106 * math.e**2 * h(s)


def exceptional_lower(s: float) -> float:
    fv = sieve_f(s)
    Fv = sieve_F(s)
    corr2 = ETA * 107 * math.e**2 * h(s)
    corr1 = ETA * 106 * math.e**2 * h(s)
    fminus = fv - corr2
    Fplus = Fv + corr1
    return fminus - 3 * EXC_EPS * (Fplus - fminus)


def integ(fn, lo: float, hi: float) -> float:
    if hi <= lo:
        return 0.0
    return quad(fn, lo, hi, epsabs=3e-10, epsrel=3e-10, limit=300)[0]


def rows(c: float, k1: float, k2: float, s2: float, s1: float) -> dict[str, float]:
    L = math.exp(c)
    theta = 0.5 - (QBAR + 20 * c) / L
    shift = 20 * c / (k1 * L)
    th = theta - k1 * shift

    out: dict[str, float] = {}
    out["L0"] = 2 * exceptional_lower(theta / k1 - shift) / (k1 * EGAMMA)
    out["A1"] = 2 / EGAMMA * integ(
        lambda t: Fplus((th - t) / t) / t**2, k1, k2
    )
    out["A2"] = 2 / (k1 * EGAMMA) * integ(
        lambda t: Fplus((th - t) / k1) / t, k1, s1
    )
    out["A3"] = 2 / (k1 * EGAMMA) * integ(
        lambda t: Fplus((th - t) / k1) / t, k1, s2
    )
    out["L4"] = 2 / (k1 * EGAMMA) * integ(
        lambda t: integ(
            lambda u: exceptional_lower((theta - t - u) / k1 - shift) / u,
            t,
            k2,
        ) / t,
        k1,
        k2,
    )
    out["L5"] = 2 / (k1 * EGAMMA) * integ(
        lambda t: integ(
            lambda u: exceptional_lower((theta - t - u) / k1 - shift) / u,
            k2,
            0.5 - 2 * k1 - t,
        ) / t,
        k1,
        k2,
    )

    factor = EXC_DENSITY * 6 * math.exp(-GAMMA) * Fplus(3 * theta - shift)
    out["B7"] = factor * integ(lambda v: math.log(v - 1) / v, 2, 1 / s1 - 1)
    out["B8"] = factor * integ(lambda v: math.log(v - 1) / v, 2, 1 / s2 - 1)
    out["B9"] = factor * integ(
        lambda t: math.log((1 - t - s1) / s1) / (t * (1 - t)), k1, s1
    )
    out["B10"] = factor * integ(
        lambda t: math.log((1 - t - s2) / s2) / (t * (1 - t)), k2, s2
    )

    residual = 0.5 + 2 * k1 - 2 * k2
    srough = (residual - QBAR / L - 10 * math.log(residual * L) / L) / k1
    scircle = 0.4 / k1 - shift
    I13 = math.log(k2 / k1) ** 4 / 24
    I14 = integ(
        lambda t: math.log(t / k1) ** 2 / (2 * t)
        * math.log((0.5 - 2 * k1 - t) / k2),
        k1,
        k2,
    )
    out["Bsw"] = (
        EXC_DENSITY * (1 + ETA) * 2 * math.exp(-2 * GAMMA) / k1**2
        * Fplus(srough) * Fplus(scircle) * (I13 + I14)
    )
    out["theta"] = theta
    out["shift"] = shift
    out["srough"] = srough
    out["scircle"] = scircle
    return out


def feasible(x: np.ndarray, gap: float) -> bool:
    k1, k2, s2, s1 = x
    rho = 1 / 3 - k1
    return all([
        k1 >= 1 / 12,
        k1 + gap <= k2,
        k2 + gap <= rho,
        rho + gap <= s2,
        s2 + gap <= s1,
        s1 < 1 / 3,
        k1 < 1 / 8,
        k2 + s2 >= 0.398,
        4 * s2 > 1,
        k1 + 3 * s1 > 1,
        k2 + s2 + 2 * s1 > 1,
        s1 < 11 / 24,
        s1 < 1 - k2 - s2,
    ])


def master(c: float, x: np.ndarray) -> float:
    k1, k2, s2, s1 = x
    r = rows(c, k1, k2, s2, s1)
    return (
        4 * r["L0"]
        - EXC_DENSITY * (r["A1"] + r["A2"] + r["A3"])
        + r["L4"] + r["L5"]
        - 2 * r["B7"] - 2 * r["B8"] - r["B9"] - r["B10"] - r["Bsw"]
    ) / 4


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=float, default=24.9)
    parser.add_argument("--gap", type=float, default=1e-4)
    parser.add_argument("--optimize", action="store_true")
    parser.add_argument("--params", nargs=4, type=float)
    args = parser.parse_args()

    if args.params:
        x = np.array(args.params)
    elif args.optimize:
        def objective(y: np.ndarray) -> float:
            if not feasible(y, args.gap):
                return 1e3 + 1e3 * sum(abs(y))
            try:
                return -master(args.c, y)
            except (ValueError, ZeroDivisionError):
                return 1e3

        result = differential_evolution(
            objective,
            bounds=[(1 / 12, 0.11), (0.0834, 0.14), (0.26, 0.325), (0.27, 0.3325)],
            tol=2e-8,
            popsize=20,
            maxiter=500,
            polish=True,
            seed=249,
            workers=1,
        )
        x = result.x
        print(result)
    else:
        x = np.array([0.0958, 0.0959, 0.3037, 0.3052])

    print("feasible", feasible(x, args.gap))
    print("params", *(f"{v:.15g}" for v in x), "rho", f"{1/3-x[0]:.15g}")
    r = rows(args.c, *x)
    for key, value in r.items():
        print(f"{key:>8} {value:.15g}")
    print("master", f"{master(args.c, x):.15g}")


if __name__ == "__main__":
    main()

