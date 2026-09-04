#!/usr/bin/env python3
"""Single forward verification entry point for the c0=24.9 proof."""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import io
import json
from pathlib import Path

from flint import arb, ctx, fmpq

import c249_interval_certificate as integrals
import c249_structural_certificate as structural


ROOT = Path(__file__).resolve().parent


def exact_ball(q: str) -> arb:
    return arb(fmpq(q))


def endpoints(value: arb, places: int = 12) -> tuple[str, str]:
    scale = 10**places
    lo = int((value.lower() * scale).floor().fmpz())
    hi = int((value.upper() * scale).ceil().fmpz())

    def decimal(k: int) -> str:
        sign = "-" if k < 0 else ""
        k = abs(k)
        return f"{sign}{k // scale}.{k % scale:0{places}d}"

    return decimal(lo), decimal(hi)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-tex", action="store_true")
    args = parser.parse_args()
    with contextlib.redirect_stdout(io.StringIO()):
        icert = integrals.main()
        scert = structural.main()

    schema = [
        ("L0", "F0", "ordinary_large_lower", "4", "lower"),
        ("A1", "F1", "ordinary_large_upper_base", "-36/35", "upper"),
        ("A2", "F2", "ordinary_large_upper_base", "-36/35", "upper"),
        ("A3", "F3", "ordinary_large_upper_base", "-36/35", "upper"),
        ("L4", "F4", "ordinary_large_lower", "1", "lower"),
        ("L5", "F5", "ordinary_large_lower", "1", "lower"),
        ("B7", "F7", "box_large", "-2", "upper"),
        ("B8", "F8", "box_large", "-2", "upper"),
        ("B9", "F9", "box_large", "-1", "upper"),
        ("B10", "F10", "box_large", "-1", "upper"),
        ("Bsw", "", "switch_large", "-1", "upper"),
    ]
    rows = []
    raw = arb(0)
    printed_lower = arb(0)
    for name, key, section, sign, direction in schema:
        ball = arb(icert[section][key] if key else icert[section])
        if not ball.is_finite() or ball.rel_accuracy_bits() < 90:
            raise ArithmeticError(f"unusable row {name}: {ball}")
        coefficient = exact_ball(sign)
        contribution = coefficient * ball / 4
        raw += coefficient * ball
        lo, hi = endpoints(ball)
        printed_lower += coefficient * arb(lo if direction == "lower" else hi) / 4
        rows.append({
            "row": name,
            "source_key": f"{section}.{key}" if key else section,
            "raw_ball_in_Theta_units": ball.str(32, more=True),
            "raw_outward_interval_12dp": [lo, hi],
            "signed_Wu_weight": sign,
            "direction": direction,
            "normalization_divisor": "4",
            "normalized_contribution": contribution.str(32, more=True),
        })

    main_ball = raw / 4
    prior = arb(icert["master_before_finite_charges"]["large"])
    if not main_ball.overlaps(prior):
        raise ArithmeticError("component and forward master expressions disagree")
    if not (main_ball > arb("0.0228") and printed_lower > arb("0.0228")):
        raise ArithmeticError("master or printed table floor is too small")

    budgets = scert["budgets"]
    budget_refs = {
        "rectangle": "eq:box-error",
        "stieltjes_and_prime_endpoints": "eq:transfer-cost",
        "ordinary_AP": "eq:ordinary-AP-cost",
        "large_exceptional_finite": "eq:exc-cost",
        "switch_remainders": "eq:switch-cost",
        "finite_combinatorial": "eq:finite-cost",
        "local_products": "eq:local-cost",
    }
    if set(budgets) != set(budget_refs):
        raise ArithmeticError("an unassigned or additional budget row exists")
    budget_rationals = {key: fmpq(value) for key, value in scert["budget_rationals"].items()}
    if set(budget_rationals) != set(budget_refs):
        raise ArithmeticError("budget rational names disagree")
    for key, value in budget_rationals.items():
        enclosure = arb(budgets[key])
        if not enclosure.is_finite() or enclosure.rel_accuracy_bits() < 90 or not enclosure.contains(arb(value)):
            raise ArithmeticError(f"budget enclosure mismatch: {key}")
    exact_total_budget = sum(budget_rationals.values(), fmpq(0))
    if exact_total_budget != fmpq(2401, 200_000):
        raise ArithmeticError("finite budget total is not 0.012005")
    total_budget = arb(exact_total_budget)
    net = main_ball - total_budget
    if not net > arb("0.0108"):
        raise ArithmeticError("final positivity assertion fails")

    tex_path = ROOT / "explicit_chen_c249_final.tex"
    if args.check_tex:
        tex = tex_path.read_text(encoding="utf-8")
        names = {
            "L0": "L_0", "L4": "L_4", "L5": "L_5",
            "A1": "A_1", "A2": "A_2", "A3": "A_3",
            "B7": "B_7", "B8": "B_8", "B9": "B_9", "B10": "B_{10}",
            "Bsw": r"B_{\rm sw}",
        }
        for row in rows:
            lo, hi = row["raw_outward_interval_12dp"]
            marker = f"% ledger:{row['row']}:{lo}:{hi}:{row['signed_Wu_weight']}"
            if tex.count(marker) != 1:
                raise ArithmeticError(f"missing/duplicate/mismatched TeX row: {marker}")
            visible = ("$" + names[row["row"]] + "$ & $" + lo + "$ & $" + hi
                       + "$ & $" + row["signed_Wu_weight"] + r"$\\")
            if tex.count(visible) != 1:
                raise ArithmeticError(f"visible table differs for {row['row']}")
        for key, ref in budget_refs.items():
            q = budget_rationals[key]
            scaled = q * 1_000_000
            integer = int(scaled)
            if scaled != fmpq(integer) or integer < 0:
                raise ArithmeticError(f"budget is not an exact nonnegative millionth: {key}")
            decimal = f"{integer // 1_000_000}.{integer % 1_000_000:06d}"
            visible = (r"\eqref{" + ref + r"} & $" + decimal + r"$\\")
            if tex.count(visible) != 1:
                raise ArithmeticError(f"visible budget table differs for {key}")
        if tex.count(r"Total & & $0.012005$\\") != 1:
            raise ArithmeticError("visible budget total differs")
        for label in ["eq:master", "sec:uniform", "sec:APformula", "thm:main"]:
            if rf"\label{{{label}}}" not in tex:
                raise ArithmeticError(f"missing TeX anchor {label}")
        if r"\end{document}" not in tex:
            raise ArithmeticError("TeX is incomplete")

    tracked = [
        "c249_verify.py", "c249_interval_certificate.py", "c249_structural_certificate.py",
        "c26_interval_certificate.py", "explicit_chen_c249_final.tex",
        "bjs_source/Explicit_Chen_-_New.tex", "wu_source/ChenDoubleSieve1_Paper.tex",
    ]
    hashes = {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in tracked}
    result = {
        "purpose": "numerical certificate; analytic implications are in the manuscript",
        "target": "all even N with log log N >= 24.9, positivity only",
        "branch_policy": "one universal q1>=37 envelope; exceptional conductors may differ",
        "precision_bits": ctx.prec,
        "rows": rows,
        "main_before_charges": main_ball.str(32, more=True),
        "lower_bound_from_printed_endpoints": printed_lower.str(32, more=True),
        "normalized_finite_budgets": budgets,
        "exact_budget_rationals": {key: str(value) for key, value in budget_rationals.items()},
        "sum_of_finite_budgets": total_budget.str(32, more=True),
        "net_arithmetic_enclosure": net.str(32, more=True),
        "rigorous_positive_floor": "0.0108",
        "tex_table_checked": args.check_tex,
        "main_integral_certificate": icert,
        "structural_certificate": scert,
        "sha256": hashes,
    }
    (ROOT / "c249_final_ledger.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    for row in rows:
        lo, hi = row["raw_outward_interval_12dp"]
        print(f"% ledger:{row['row']}:{lo}:{hi}:{row['signed_Wu_weight']}")
    print("Main:", result["main_before_charges"])
    print("Printed endpoint lower bound:", result["lower_bound_from_printed_endpoints"])
    print("Budget:", result["sum_of_finite_budgets"])
    print("Net:", result["net_arithmetic_enclosure"])
    print("TeX table checked:", args.check_tex)
    print("Written c249_final_ledger.json")


if __name__ == "__main__":
    main()

