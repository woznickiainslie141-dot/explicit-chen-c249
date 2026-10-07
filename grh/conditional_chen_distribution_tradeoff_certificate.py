"""Necessary C,X constraints for all-modulus quantitative distribution inputs.

The general proof is in the manuscript. Arb certifies the displayed witnesses;
exact residue-class histograms provide a separate finite arithmetic check.
No distribution upper bound or sufficient EH onset is proved here.
"""

from fractions import Fraction
from pathlib import Path
import hashlib
import json
from math import isqrt
from flint import arb, ctx


ROOT = Path(__file__).resolve().parent
ctx.prec = 192
B = arb(4) * 10**18


def obstruction(theta, y):
    """x^(-1) log(x)^4 times ((x^theta-1)/(2 log(x^theta))-sqrt(x))."""
    return (y**3 / (2*theta) * ((-(1-theta)*y).exp() - (-y).exp())
            - y**4 * (-y/2).exp())


def direct_obstruction(theta, y):
    x, h = y.exp(), (theta*y).exp()
    assert h >= 72 and h < x
    return ((h-1)/(2*h.log()) - x.sqrt()) * y**4/x


def finite_histogram_checks():
    # Compute actual prime counts in every reduced class, rather than invoking
    # the P mod (q-1) proof used to obtain the general lower bound.
    primes = [n for n in range(2, 201)
              if all(n % d for d in range(2, isqrt(n)+1))]
    cases = 0
    for x in range(3, 201):
        px = [p for p in primes if p <= x]
        total = Fraction(0)
        count_h = 0
        for h in range(2, x):
            if h in primes:
                histogram = [0]*h
                for p in px:
                    histogram[p % h] += 1
                target = Fraction(len(px), h-1)
                total += max(abs(Fraction(histogram[a])-target)
                             for a in range(1, h))
                count_h += 1
            lower = arb(count_h)/2-arb(x).sqrt()
            exact = arb(total.numerator)/total.denominator
            assert exact >= lower
            cases += 1
    return cases


def main():
    if not __debug__:
        raise RuntimeError("Run without Python -O.")
    rows = []
    inputs = [
        (3, 4, 20, "1.1697", 43, 22, None),
        (9, 10, 10, "602.35", 143, 57, 18),
        (19, 20, 10, "4858.17", 336, 130, 57),
        (99, 100, 10, "25860.94", 2246, 834, 492),
    ]
    for p, q, budget, rounded, c1_witness, gate_witness, main_witness in inputs:
        theta = arb(p)/q
        at_b = obstruction(theta, B.log())
        assert at_b > arb(rounded)
        assert abs(at_b-direct_obstruction(theta, B.log())) < arb("1e-45")
        y = arb(c1_witness)
        c1_lower = obstruction(theta, y)
        assert theta*y >= arb(72).log() and c1_lower > 1
        assert abs(c1_lower-direct_obstruction(theta, y)) < arb("1e-45")
        m = arb(gate_witness)
        gate_lower, gate_cap = obstruction(theta, m), m*m/budget
        assert theta*m >= arb(72).log() and gate_lower > gate_cap
        # 3 U_* f(3 theta), with U_*=1.32 exp(-gamma), is an upper
        # bound on the current direct linear-sieve main coefficient.
        main_cap = arb("7.92")*(3*theta-1).log()/(3*theta)
        row = {
            "theta": f"{p}/{q}", "budget": budget,
            "C_lower_if_X_le_B_ball": str(at_b),
            "displayed_C_lower": rounded,
            "C_equals_one_onset_log_witness": c1_witness,
            "C_equals_one_obstruction_ball": str(c1_lower),
            "C_equals_one_strict_gap_ball": str(c1_lower-1),
            "current_budget_endpoint_log_witness": gate_witness,
            "budget_obstruction_ball": str(gate_lower),
            "budget_constant_cap_ball": str(gate_cap),
            "budget_strict_gap_ball": str(gate_lower-gate_cap),
            "direct_sieve_main_coefficient_cap_ball": str(main_cap),
        }
        if main_witness is not None:
            m = arb(main_witness)
            lower, cap = obstruction(theta, m), main_cap*m*m
            assert theta*m >= arb(72).log() and lower > cap
            row.update({
                "direct_sieve_positivity_log_witness": main_witness,
                "direct_sieve_positivity_strict_gap_ball": str(lower-cap),
            })
        rows.append(row)
    exact_checks = finite_histogram_checks()
    out = {
        "status": "PASS: necessary distribution/onset constraints only",
        "precision_bits": ctx.prec, "B": "4000000000000000000",
        "formula": "F_theta(y)=y^3/(2theta)*(exp(-(1-theta)y)-exp(-y))-y^4*exp(-y/2)",
        "dependency": "prime-count lower bound pi(J)>J/log J for J>=71 (JS Lemma 3.1); manuscript integrality proof",
        "rows": rows,
        "finite_checks": {
            "exact_prime_modulus_histogram_cases": exact_checks,
            "integer_x_range": [3, 200],
            "scope": "prime-modulus subtotal at every integer H from 2 to x-1; does not prove the general analytic lemma",
        },
        "witness_scope": "Specific strict witnesses, not globally optimal lower limits; passing them is not sufficient for the distribution assumption.",
        "not_proved": ["any distribution upper bound", "a sufficient numerical pair C,X",
                       "a full Chen bridge", "a theoretical limit for other sieve or remainder models"],
        "input_sha256": {Path(__file__).name: hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
    }
    Path(__file__).with_suffix(".json").write_text(
        json.dumps(out, ensure_ascii=False, indent=2)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status":out["status"], "finite_checks":exact_checks,
                      "witness_rows":len(rows)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
