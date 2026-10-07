"""Check manuscript structure and agreement with the current numeric ledger.

This is an artifact-consistency audit. It does not prove the analytic lemmas,
authenticate external results, compile TeX, or verify a full Chen bridge.
"""
import ast
import hashlib
import json
import math
import re
from pathlib import Path
from flint import arb, ctx

ctx.prec = 192
ROOT = Path(__file__).resolve().parent


def escaped(text, position):
    count = 0
    while position and text[position-1] == "\\":
        position -= 1
        count += 1
    return count % 2 == 1


def without_comments(text):
    lines = []
    for line in text.splitlines():
        comment = next((i for i, ch in enumerate(line)
                        if ch == "%" and not escaped(line, i)), len(line))
        lines.append(line[:comment])
    return "\n".join(lines)


def audit():
    if not __debug__:
        raise RuntimeError("Run this audit without Python's -O option.")
    tex_path = ROOT / "conditional_chen_grh_bridge.tex"
    raw = tex_path.read_bytes()
    tex = raw.decode("utf-8")
    assert not [(i, b) for i, b in enumerate(raw)
                if b < 32 and b not in (9, 10)], "Unexpected source control bytes"
    clean = without_comments(tex)
    depth = 0
    for i, ch in enumerate(clean):
        if ch in "{}" and not escaped(clean, i):
            depth += 1 if ch == "{" else -1
            assert depth >= 0, "Unmatched closing brace"
    assert depth == 0, "Unmatched opening brace"
    stack = []
    for kind, env in re.findall(r"\\(begin|end)\{([^}]+)\}", clean):
        if kind == "begin":
            stack.append(env)
        else:
            assert stack and stack.pop() == env, f"Mismatched environment: {env}"
    assert not stack
    labels = re.findall(r"\\label\{([^}]+)\}", clean)
    assert len(labels) == len(set(labels)), "Duplicate labels"
    references = re.findall(r"\\(?:eqref|ref|pageref)\{([^}]+)\}", clean)
    assert set(references) <= set(labels), "Undefined references"
    bibkeys = re.findall(r"\\bibitem(?:\[[^\]]*\])?\{([^}]+)\}", clean)
    assert len(bibkeys) == len(set(bibkeys))
    citations = []
    for group in re.findall(r"\\cite(?:\[[^\]]*\])*\{([^}]+)\}", clean):
        citations.extend(key.strip() for key in group.split(","))
    assert set(citations) <= set(bibkeys), "Undefined citations"

    scripts = [
        "conditional_chen_composite_certificate.py",
        "conditional_chen_rankin_certificate.py",
        "conditional_chen_rosser_checks.py",
        "conditional_chen_refined_eh_search.py",
        "conditional_chen_uniform_product_certificate.py",
        "conditional_chen_uniform_product_checks.py",
        "conditional_chen_uniform_grh_probe.py",
        "conditional_chen_degreewise_rankin_certificate.py",
        "conditional_chen_degreewise_rankin_checks.py",
        "conditional_chen_degreewise_eh_search.py",
        "conditional_chen_local_rankin_probe.py",
        "conditional_chen_local_eh_probe.py",
        "conditional_chen_quadratic_rankin_certificate.py",
        "conditional_chen_quadratic_rankin_checks.py",
        "conditional_chen_quadratic_eh_search.py",
        "conditional_chen_quadratic_rankin_probe.py",
        "conditional_chen_prefix_band_certificate.py",
        "conditional_chen_prefix_band_checks.py",
        "conditional_chen_prefix_band_eh_search.py",
        "conditional_chen_prefix_band_probe.py",
        "conditional_chen_sharp_support_certificate.py",
        "conditional_chen_sharp_endpoint_certificate.py",
        "conditional_chen_distribution_tradeoff_certificate.py",
        "conditional_chen_weighted_distribution_certificate.py",
        "conditional_chen_dense_band_grh_certificate.py",
        "conditional_chen_sharp_grh_parameter_probe.py",
        "conditional_chen_prefix_block_certificate.py",
        "conditional_chen_prefix_block_checks.py",
        "conditional_chen_prefix_block_endpoint_certificate.py",
        "conditional_chen_prefix_block_refinement.py",
        "conditional_chen_sparse_modulus_certificate.py",
        "conditional_chen_prefix_block_endpoint_checks.py",
        "conditional_chen_prefix_block_update.py",
        "conditional_chen_binned_support_certificate.py",
        "conditional_chen_binned_support_checks.py",
        "conditional_chen_rough_interval_certificate.py",
        "conditional_chen_rough_interval_checks.py",
        "conditional_chen_rough_interval_update.py",
        "conditional_chen_accepted_support_certificate.py",
        "conditional_chen_accepted_support_checks.py",
        "conditional_chen_accepted_interval_certificate.py",
        "conditional_chen_accepted_interval_checks.py",
        "conditional_chen_accepted_interval_update.py",
        "conditional_chen_full_rosser_bridge_probe.py",
        "conditional_chen_weighted_rosser_support_certificate.py",
        "conditional_chen_full_rosser_w2_probe.py",
        "conditional_chen_full_rosser_bridge_checks.py",
        "conditional_chen_refresh_notes.py",
        "conditional_chen_source_audit.py",
    ]
    for name in scripts:
        ast.parse((ROOT / name).read_text(encoding="utf-8"), filename=name)
    base = json.loads((ROOT / "conditional_chen_composite_certificate.json")
                             .read_text(encoding="utf-8"))
    certificate = json.loads((ROOT / "conditional_chen_sharp_endpoint_certificate.json")
                             .read_text(encoding="utf-8"))
    standalone = json.loads((ROOT / "conditional_chen_prefix_band_certificate.json")
                            .read_text(encoding="utf-8"))
    checks = json.loads((ROOT / "conditional_chen_rosser_checks.json")
                       .read_text(encoding="utf-8"))
    search = json.loads((ROOT / "conditional_chen_prefix_band_eh_search.json")
                       .read_text(encoding="utf-8"))
    uniform = json.loads((ROOT / "conditional_chen_uniform_product_certificate.json")
                         .read_text(encoding="utf-8"))
    uniform_checks = json.loads((ROOT / "conditional_chen_uniform_product_checks.json")
                                .read_text(encoding="utf-8"))
    assert certificate["status"].startswith("PASS")
    assert certificate["finite_inputs_regenerated_in_this_run"] is False
    assert certificate["support_maxima_regenerated_in_this_run"] is True
    assert certificate["baseline_GRH_arithmetic_reproduction"].startswith("PASS")
    for name, digest in certificate["input_sha256"].items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == digest
    support = json.loads((ROOT / "conditional_chen_sharp_support_certificate.json")
                        .read_text(encoding="utf-8"))
    assert certificate["support_certificate"] == support
    assert support["finite_checks"]["coefficient_and_deleted_subset_cases"] == 52488
    assert support["finite_checks"]["accepted_coefficients_checked"] == 42226
    assert search["tested_candidates"] == 96
    assert search["not_an_effective_EH_only_numeric_record"] is True
    assert base["uniform_product_certificate"]["regenerated_in_this_run"] is True
    assert base["GRH"]["log_N0"] == 31650
    assert uniform["S"] == 100000000 and uniform["M_cutoff"] == 10000000
    assert uniform["complete_odd_prime_count"] == 5761454
    assert uniform["checked_prime_jumps"] == 5760226
    uniform_rows = {row["T"]: row for row in uniform["rows"]}
    for case_name in ["GRH", "EH_3_4", "EH_9_10", "EH_19_20", "EH_99_100"]:
        case = certificate[case_name]
        pre = case["presieve"]
        assert pre == base[case_name]["presieve"]
        assert pre["prime_generation"] == "complete integer sieve"
        assert pre["logarithm_base"] == "e"
        assert pre["band_optimizers_with_strict_signs"] > 0
        sr = case["sharp_support"]
        assert sr["complete_medium_prime_count"] == pre["medium_prime_count"]
        assert sr["kappa"] == ("11/169" if case_name == "GRH" else "5/49")
        assert abs(arb(case["effective_log_support_cost_ball"])
                   - arb(pre["log_support_cost_ball"])
                   - arb(sr["kappa"]).log()) < arb("1e-50")
        assert arb(pre["upper_relative_Rosser_gap_ball"]) < arb(pre["ordinary_upper_relative_gap_ball"])
        assert arb(pre["lower_relative_Rosser_gap_ball"]) < arb(pre["ordinary_lower_relative_gap_ball"])
        row = uniform_rows[pre["T"]]
        eps = arb(row["ordered_uniform_epsilon_ball"])
        assert abs(eps-arb(pre["uniform_epsilon_ball"])) < arb("1e-50")
        assert 0 < eps < arb(case["epsilon"])
        assert eps < arb(pre["coarse_uniform_epsilon_ball"])
        assert pre["uniform_product_input"] == "ordered finite certificate supplied by caller"
        bound = arb(row["ordered_logK_numerator"])/row["ordered_logK_denominator"]
        for key in ["finite_ordered_log_drawdown_ball",
                    "finite_to_far_log_ratio_ball", "far_to_far_log_ratio_ball"]:
            assert arb(row[key]) < bound
        table_fragment = f"{pre['T']}&{row['ordered_logK_numerator']}&"
        assert table_fragment in tex
    assert uniform_checks["status"].startswith("PASS")
    assert uniform_checks["integer_sequence_checks"] == 4096
    assert uniform_checks["direct_prime_product_endpoint_pairs"] == 1954
    assert uniform_checks["integer_and_prime_tail_checks"] == 9
    assert "4096" in tex and "1954" in tex
    assert r"\mathcal K_T<1+\epsilon" in tex
    assert r"\label{lem:ordered}" in tex
    assert r"\label{lem:degreewise}" in tex
    assert r"\label{lem:quadratic}" in tex
    assert r"\label{lem:prefixband}" in tex
    assert r"\label{lem:sharpsupport}" in tex
    degree_checks = json.loads((ROOT / "conditional_chen_degreewise_rankin_checks.json")
                               .read_text(encoding="utf-8"))
    assert degree_checks["status"].startswith("PASS")
    assert degree_checks["explicit_coefficient_configurations"] == 9
    assert degree_checks["first_failure_and_deleted_prime_checks"] == 2304
    quadratic_checks = json.loads((ROOT / "conditional_chen_quadratic_rankin_checks.json")
                                  .read_text(encoding="utf-8"))
    assert quadratic_checks["status"].startswith("PASS")
    assert quadratic_checks["explicit_coefficient_configurations"] == 9
    assert quadratic_checks["first_failure_and_deleted_prime_checks"] == 2304
    band_checks = json.loads((ROOT / "conditional_chen_prefix_band_checks.json")
                             .read_text(encoding="utf-8"))
    assert band_checks["status"].startswith("PASS")
    assert band_checks["atoms_are_primes"] is False
    assert band_checks["explicit_coefficient_configurations"] == 12
    assert band_checks["first_failure_and_deleted_subset_checks"] == 2688
    assert band_checks["first_failure_product_band_checks"] == 25
    assert band_checks["two_prime_first_failure_checks"] == 1
    grh = certificate["GRH"]
    assert grh["presieve"] == standalone
    assert grh["log_N0"] == 31000
    assert grh["alpha"] == "4215/1000000"
    assert grh["eta_plus"] == "1/2768" and grh["eta_minus"] == "1/2784"
    assert grh["presieve"]["u"] == "49/8"
    assert grh["presieve"]["sigma_candidates"] == [
        "0","1/10","7/50","9/50","11/50","7/25","9/25"]
    assert grh["presieve"]["separate_degree_cutoff"] == 12
    assert grh["presieve"]["band_sigma_candidates"] == [
        "2/25","3/25","4/25","1/5","6/25"]
    assert r"\alpha=4215/10^6" in tex
    assert r"N\ge e^{31000}" in tex
    assert r"\kappa=11/169" in tex and r"\kappa=5/49" in tex
    assert not re.search(r"\b(?:33200|36000|36800|38000|915|247|1036|289|1075|298|1054|293)\b", tex)
    assert arb(grh["support_slack_ball"]) > 0
    assert arb(grh["lower_composite_factor_ball"]) > 0
    assert arb(grh["margin_ball"]) > arb(".0006")
    assert arb(grh["F_s_ball"]) > 1
    assert 0 < arb(grh["h_integral_ball"]) < (arb(8)/3).log()

    advertised = re.search(r"Surplus & \$>([0-9.]+)\$", tex).group(1)
    assert arb(grh["margin_ball"]) > arb(advertised)
    for key in ["lower_ball", "half_upper_ball", "half_switch_ball", "AP_budget_ball"]:
        displayed = f"{float(arb(grh[key])):.12f}"
        assert displayed in tex
        assert abs(arb(grh[key])-arb(displayed)) < arb("5.1e-13")
    assert arb(grh["half_rectangle_error_ball"]) < arb(grh["rectangle_display_bound"])
    assert arb(grh["finite_error_ball"]) < arb(grh["finite_display_bound"])
    positive = certificate["GRH_positivity"]
    assert positive["presieve"] == grh["presieve"]
    assert positive["log_N0"] == 30824 and positive["alpha"] == "4238/1000000"
    assert positive["claimed_margin"] == ".00001"
    assert arb(positive["margin_ball"]) > arb(".0000337669646")
    assert arb(positive["support_slack_ball"]) > arb(".002084")
    assert arb(positive["half_rectangle_error_ball"]) < arb("1e-248")
    assert arb(positive["finite_error_ball"]) < arb("1e-1660")
    assert r"N\ge e^{30824}" in tex and r"\alpha=4238/10^6" in tex
    log_log_display = f"{float(arb(positive['log_log_N0_ball'])):.9f}"
    assert log_log_display in tex
    assert abs(arb(positive["log_log_N0_ball"])-arb(log_log_display)) < arb("5.1e-10")
    digit_quotient = arb(positive["log_N0"])/arb(10).log()
    digits = math.floor(float(digit_quotient))+1
    assert digits-1 < digit_quotient < digits
    assert "$"+str(digits)+"$" in tex
    for key, start in [("EH_3_4",774),("EH_9_10",210),("EH_19_20",165),("EH_99_100",141)]:
        weaker = certificate["distribution_positivity"][key]
        assert weaker["log_start"] == start and weaker["claimed_margin"] == ".00001"
        assert arb(weaker["surplus_ball"]) > arb(".00001")
        assert weaker["presieve"] == certificate[key]["presieve"]
        assert weaker["epsilon"] == certificate[key]["epsilon"]
        assert weaker["eta_plus"] == certificate[key]["eta_plus"]
        assert weaker["eta_minus"] == certificate[key]["eta_minus"]
    levels = {"EH_3_4":"3/4","EH_9_10":"9/10",
              "EH_19_20":"19/20","EH_99_100":"99/100"}
    for name, start in [("EH_3_4",786),("EH_9_10",212),
                        ("EH_19_20",166),("EH_99_100",142)]:
        case = certificate[name]
        assert case["log_start"] == start
        assert case["not_an_effective_EH_only_numeric_record"] is True
        assert "no RH or GRH" in case["assumptions"]
        assert arb(case["surplus_ball"]) > arb(case["claimed_margin"])
        pre = case["presieve"]
        selected = search["best"][name][0]
        for parameter in ["T", "exact_cutoff", "u"]:
            search_key = "w" if parameter == "exact_cutoff" else parameter
            assert str(pre[parameter]) == str(selected[search_key])
        assert selected["log_start"] == base[name]["log_start"]
        assert start < selected["log_start"]
        assert pre["sigma_candidates"] == selected["sigma_candidates"]
        assert pre["band_sigma_candidates"] == selected["band_sigma_candidates"]
        assert pre["separate_degree_cutoff"] == selected["degree_cutoff"]
        for parameter in ["epsilon", "eta_plus", "eta_minus"]:
            assert case[parameter] == f"1/{selected[parameter+'_inverse']}"
        assert (case["C1"], case["C2"]) == (selected["C1"], selected["C2"])
        row = "&".join([
            levels[name],
            str(pre["exact_cutoff"]), str(pre["T"]), str(pre["u"]),
            str(pre["separate_degree_cutoff"]), case["epsilon"].split("/")[1],
            case["eta_plus"].split("/")[1], case["eta_minus"].split("/")[1],
            f"({case['C1']},{case['C2']})",
        ])
        assert row in tex, "Manuscript parameter table disagrees with the certificate"
        for sign in ["plus", "minus"]:
            assert arb(pre[f"{'upper' if sign == 'plus' else 'lower'}_relative_Rosser_gap_ball"]) < arb(case[f"eta_{sign}"])
    assert checks["status"].startswith("PASS")
    assert checks["total_pointwise_pattern_checks"] == 602368
    assert checks["total_gap_identity_and_deleted_subset_checks"] == 3088
    assert all(case["nonzero_medium_gap"]
               for case in checks["fractional_power_universes"])
    assert {case["u"] for case in checks["fractional_power_universes"]} == {
        "5/1", "9/2", "35/8", "25/4", "17/4", "49/8", "39/8", "33/8",
    }
    assert "602368" in tex and "3088" in tex and "2304" in tex
    assert "2688" in tex and "52488" in tex

    documents = [
        "conditional_chen_results_20261003.md",
        "conditional_chen_audit_20261003.md",
        "conditional_chen_research_20261003.md",
    ]
    for name in documents:
        text = (ROOT / name).read_text(encoding="utf-8")
        if "research" in name:
            text = text.split("**追加更正")[0]  # The rest is marked historical.
        assert "exp(31000)" in text
        assert "exp(30824)" in text
        assert "786" in text and "212" in text
        assert "166" in text and "142" in text
        assert "774" in text and "210" in text
        assert "165" in text and "141" in text
    for name in documents[:2]:
        text = (ROOT / name).read_text(encoding="utf-8")
        headings = re.findall(r"^## (.+)$", text, flags=re.MULTILINE)
        assert len(headings) == len(set(headings)), "Duplicated current ledger sections"
        assert "4215/10^6" in text
        assert arb(grh["margin_ball"]) > arb("0.0006")
        assert "0.0107006692223" in text
        assert "0.0228594687806" in text
        assert arb(certificate["EH_3_4"]["surplus_ball"]) > arb("0.0107006692223")
        assert arb(certificate["EH_9_10"]["surplus_ball"]) > arb("0.0228594687806")

    tradeoff = json.loads((ROOT / "conditional_chen_distribution_tradeoff_certificate.json")
                         .read_text(encoding="utf-8"))
    assert tradeoff["status"].startswith("PASS")
    for name, digest in tradeoff["input_sha256"].items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == digest
    assert tradeoff["finite_checks"]["exact_prime_modulus_histogram_cases"] == 19701
    expected_tradeoff = [
        ("3/4", 20, "1.1697", 43, 22),
        ("9/10", 10, "602.35", 143, 57),
        ("19/20", 10, "4858.17", 336, 130),
        ("99/100", 10, "25860.94", 2246, 834),
    ]
    assert len(tradeoff["rows"]) == len(expected_tradeoff)
    for row, (theta, budget, lower, c1, gate) in zip(tradeoff["rows"], expected_tradeoff):
        assert (row["theta"], row["budget"], row["displayed_C_lower"],
                row["C_equals_one_onset_log_witness"],
                row["current_budget_endpoint_log_witness"]) == (theta, budget, lower, c1, gate)
        assert arb(row["C_lower_if_X_le_B_ball"]) > arb(lower)
        assert arb(row["C_equals_one_obstruction_ball"]) > 1
        assert arb(row["C_equals_one_strict_gap_ball"]) > 0
        assert arb(row["budget_obstruction_ball"]) > arb(gate)**2/budget
        assert arb(row["budget_strict_gap_ball"]) > 0
        assert lower in tex and str(c1) in tex and str(gate) in tex
    assert tradeoff["rows"][-1]["direct_sieve_positivity_log_witness"] == 492
    assert arb(tradeoff["rows"][-1]["direct_sieve_positivity_strict_gap_ball"]) > 0
    assert r"\label{prop:tradeoff}" in tex and r"\label{eq:integralitytradeoff}" in tex
    assert "19701" in tex and "492" in tex and "2246" in tex and "834" in tex
    for name in documents:
        note = (ROOT / name).read_text(encoding="utf-8").split("**追加更正")[0]
        assert "exp(834)" in note and "exp(2246)" in note and "492" in note

    weighted = json.loads((ROOT / "conditional_chen_weighted_distribution_certificate.json")
                          .read_text(encoding="utf-8"))
    assert weighted["status"].startswith("PASS") and weighted["precision_bits"] == 192
    assert weighted["C_and_X_are_unspecified_distribution_inputs"] is True
    assert weighted["finite_prime_product_and_gap_inputs_regenerated"] is False
    assert "in the sieve domain" in weighted["signed_input"]
    for name, digest in weighted["input_sha256"].items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == digest
    expected_weighted = [
        ("3/4", [722,723,740,1632]), ("9/10", [200,202,222,936]),
        ("19/20", [158,161,182,872]), ("99/100", [136,138,159,829]),
    ]
    assert len(weighted["endpoint_rows"]) == len(expected_weighted)
    for row, (theta, starts) in zip(weighted["endpoint_rows"], expected_weighted):
        assert row["theta"] == theta and row["reproduced_saved_main_and_surplus"]
        selected = [case for case in row["endpoint_cases"] if case["C_upper_budget"]]
        assert [case["C_upper_budget"] for case in selected] == [1,1000,10000,1000000]
        assert [case["log_start"] for case in selected] == starts
        assert theta+"&"+"&".join(str(x) for x in starts) in tex
        for case in selected:
            assert arb(case["margin_ball"]) > arb("0.00001")
            assert arb(case["preceding_integer_margin_ball"]) < arb("0.00001")
            assert 2 < arb(case["s_ball"]) < 3 and arb(case["factor_ball"]) > 0
        zero = row["endpoint_cases"][0]
        assert zero["C_upper_budget"] == 0
        assert arb(zero["preceding_integer_margin_ball"]) < 0
    small = weighted["exact_small_cases"]
    assert (small["even_N_count"], small["composite_cases"],
            small["pointwise_prime_divisor_patterns"], small["exact_signed_identity_checks"]) == (
                1156,10748,99158,21496)
    assert small["P2_transfer_checked_for_all_selected_N"] is True
    witness = weighted["support_class_checks"]["exact_prime_witness"]
    assert witness["coefficient"] == -1 and witness["p_cubed"] > witness["level_H"]
    assert witness["p_cubed"] == witness["p"]**3
    assert weighted["support_class_checks"]["two_factor_well_factorability_not_decided"] is True
    for label in ("lem:signed", "eq:signedinput", "eq:signedcriterion",
                  "prop:weighted", "prop:factorability"):
        assert "\\label{"+label+"}" in tex
    for count in ("10748","99158","21496","999983"):
        assert count in tex
    weighted_note = "conditional_chen_weighted_distribution_20261004.md"
    for name in documents+[weighted_note]:
        note = (ROOT/name).read_text(encoding="utf-8").split("**追加更正")[0]
        assert "max{X,exp(" in note and "829" in note
        for theta, starts in expected_weighted:
            assert "| "+theta+" | "+" | ".join(str(x) for x in starts)+" |" in note

    dense = json.loads((ROOT/"conditional_chen_dense_band_grh_certificate.json")
                      .read_text(encoding="utf-8"))
    assert dense["status"].startswith("PASS") and dense["precision_bits"] == 192
    assert dense["uniform_product_inputs_regenerated_in_this_run"] is False
    assert dense["baseline_30824_scalar_reproduction"].startswith("PASS")
    assert len(dense["direct_modified_budget_checks"]) == 2
    assert all(x.startswith("PASS") for x in dense["direct_modified_budget_checks"])
    for name, digest in dense["input_sha256"].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == digest, name
    finite = json.loads((ROOT/"conditional_chen_dense_band_grh_finite.json")
                       .read_text(encoding="utf-8"))
    for name, digest in finite["finite_generator_sha256"].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == digest, name
    dense_pre = dense["GRH_positivity"]["presieve"]
    assert {k:dense_pre[k] for k in finite["presieve"]} == finite["presieve"]
    assert dense_pre["prime_generation"] == "complete integer sieve"
    assert dense_pre["medium_prime_count"] == 348510
    assert dense_pre["sigma_candidates"] == standalone["sigma_candidates"]
    assert dense_pre["band_sigma_candidates"] == [
        "2/25","3/25","4/25","1/5","6/25","7/25","8/25","9/25"]
    assert dense_pre["separate_degree_cutoff"] == 12
    assert arb(dense_pre["upper_relative_Rosser_gap_ball"]) < arb(".000349096808")
    assert arb(dense_pre["lower_relative_Rosser_gap_ball"]) < arb(".000347781454")
    assert dense["GRH_dense_only_beta35"]["log_N0"] == 30779
    assert dense["beta_candidates"] == [
        "7/2","351/100","88/25","353/100","177/50","71/20","89/25","357/100"]
    assert len(dense["beta_endpoint_rows"]) == 8
    for row in dense["beta_endpoint_rows"]:
        assert arb(row["margin_ball"]) > arb(".00001")
        assert arb(row["preceding_margin_ball"]) < arb(".00001")
    for key,start,alpha,beta,claim in [
        ("GRH_positivity",30747,"4255/1000000","88/25",".00001"),
        ("GRH_at_31000",31000,"4218/1000000","351/100",".001")]:
        row = dense[key]
        assert (row["log_N0"],row["alpha"],row["beta"]) == (start,alpha,beta)
        assert row["epsilon"] == "1/41219"
        assert (row["eta_plus"],row["eta_minus"]) == ("1/2864","1/2875")
        assert row["presieve"] == dense_pre
        assert row["tight_square_free_AP_coefficient"] is True
        assert arb(row["margin_ball"]) > arb(claim)
        assert arb(row["support_slack_ball"]) > 0
        assert arb(row["half_rectangle_error_ball"]) < arb(row["rectangle_display_bound"])
        assert arb(row["finite_error_ball"]) < arb(row["finite_display_bound"])
        assert arb(row["AP_leading_coefficient_ball"]) < arb(".18")
        assert r"N\ge e^{"+str(start)+"}" in tex
    assert r"\label{sec:densegrh}" in tex and r"\label{eq:denseap}" in tex
    assert r"\beta$ & $88/25$ & $351/100$" in tex
    assert r"\alpha$ & $4255/10^6$ & $4218/10^6$" in tex
    assert "0.0000139212861" in tex and "0.0010866414727" in tex
    assert "$13354$" in tex and "10.333547707" in tex
    for name in documents:
        note = (ROOT/name).read_text(encoding="utf-8").split("**追加更正")[0]
        assert "exp(30747)" in note and "0.001" in note
        assert "88/25" in note and "351/100" in note
        assert "30779" in note and "1/41219" in note
        assert "1/2864" in note and "1/2875" in note

    prefix = json.loads((ROOT/'conditional_chen_prefix_block_refinement.json').read_text(encoding='utf-8'))
    sparse = json.loads((ROOT/'conditional_chen_sparse_modulus_certificate.json').read_text(encoding='utf-8'))
    endpoint_checks = json.loads((ROOT/'conditional_chen_prefix_block_endpoint_checks.json').read_text(encoding='utf-8'))
    prefix_checks = json.loads((ROOT/'conditional_chen_prefix_block_checks.json').read_text(encoding='utf-8'))
    initial = json.loads((ROOT/'conditional_chen_prefix_block_endpoint_certificate.json').read_text(encoding='utf-8'))
    for result in [prefix,sparse,endpoint_checks,prefix_checks,initial]:
        assert result['status'].startswith('PASS')
        for name,digest in result['input_sha256'].items():
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    assert (prefix_checks['deleted_prime_subsets'],prefix_checks['exact_first_failure_identities'],
            prefix_checks['pointwise_divisor_patterns'],prefix_checks['exact_block_coefficient_groups'])==(4608,9216,104976,132)
    current=sparse['GRH_positivity']; strong=sparse['GRH_at_31000']
    assert current['log_N0']==26820 and current['alpha']=='4267/1000000'
    assert strong['log_N0']==31000 and strong['alpha']=='3706/1000000'
    assert current['beta']==strong['beta']=='3'
    assert current['presieve']['u']=='87/16' and current['presieve']['log_bin_width']=='1/200'
    assert current['presieve']['maximum_actual_prefix_degree']==22
    assert current['presieve']['bin_limit']==16775 and current['presieve']['nonempty_prime_blocks']==1810
    assert (current['epsilon'],current['eta_plus'],current['eta_minus'])==('1/41219','1/2706','1/2709')
    assert current['sparse_modulus_count_budget']=='1/1653'
    assert arb(current['sparse_modulus_count']['selected']['count_fraction_ball'])<arb('1/1653')
    assert sparse['exact_small_Rankin_count_checks']==18
    assert len(sparse['complete_Rankin_products'])==3
    for row,claim,rect,finite_bound in [(current,'.00001','1e-212','1e-1446'),
                                       (strong,'.0179','1e-250','1e-1673')]:
        assert arb(row['margin_ball'])>arb(claim) and arb(row['support_slack_ball'])>0
        assert arb(row['half_rectangle_error_ball'])<arb(rect)
        assert arb(row['finite_error_ball'])<arb(finite_bound)
        assert arb(row['presieve']['upper_relative_Rosser_gap_ball'])<arb(row['eta_plus'])
        assert arb(row['presieve']['lower_relative_Rosser_gap_ball'])<arb(row['eta_minus'])
    assert [row['log_N0'] for row in endpoint_checks['GRH_checks']]==[26820,31000]
    assert all(row['quantities_independently_reproduced']==8 for row in endpoint_checks['GRH_checks'])
    assert len(endpoint_checks['distribution_checks'])==16
    expected_prefix=[('3/4',[627,629,648,1510]),('9/10',[171,174,198,885]),
                     ('19/20',[135,138,162,825]),('99/100',[115,118,143,791])]
    for theta,starts in expected_prefix:
        cells=[r for r in prefix['distribution_endpoints'] if r['theta']==theta]
        assert [r['log_start'] for r in cells]==starts
        assert [r['C_upper_budget'] for r in cells]==[1,1000,10000,1000000]
        assert all(arb(r['margin_ball'])>arb('.00001') and 2<arb(r['s_ball'])<3 for r in cells)
    for label in ['lem:prefixblocks','lem:sparsecount','eq:sparsefraction','eq:sparseap',
                  'sec:prefixgrh','prop:prefixdistribution','sec:prefixdistribution']:
        assert r'\label{'+label+'}' in tex
    assert '$26820$' in tex and '$0.0179$' in tex
    assert '$4267/1000000$' in tex and '$1853/500000$' in tex
    assert r'$s$ & $3.965864$ & $3.970352$' in tex
    assert '4608' in tex and '9216' in tex and '104976' in tex
    prefix_note='conditional_chen_prefix_block_20261004.md'
    for name in documents+[prefix_note]:
        note=(ROOT/name).read_text(encoding='utf-8').split('**追加更正')[0]
        assert 'exp(26820)' in note and '0.0179' in note and '1/1653' in note
        assert 'max{X,exp(A)}' in note and '31.052421204' in note
        for theta,starts in expected_prefix:
            assert '| '+theta+' | '+' | '.join(str(x) for x in starts)+' |' in note

    rough = json.loads((ROOT/'conditional_chen_rough_interval_certificate.json').read_text(encoding='utf-8'))
    rough_checks = json.loads((ROOT/'conditional_chen_rough_interval_checks.json').read_text(encoding='utf-8'))
    binned_checks = json.loads((ROOT/'conditional_chen_binned_support_checks.json').read_text(encoding='utf-8'))
    binned_comparison = json.loads((ROOT/'conditional_chen_binned_support_certificate.json').read_text(encoding='utf-8'))
    for name,h in binned_comparison['input_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h,name
    assert binned_comparison['sigma']=='4/5'
    assert arb(binned_comparison['ordinary_Rankin_to_bin_bound_ratio_ball'])>63
    for result in [rough,rough_checks,binned_checks]:
        assert result['status'].startswith('PASS')
        for name,h in result['input_sha256'].items():
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h,name
    previous=current;current=rough['GRH_positivity']
    assert (current['log_N0'],current['alpha'],current['beta'])==(24799,'4234/1000000','209/100')
    assert current['valid_log_interval']==[24799,26820]
    assert rough['analytic_tail_join_log_N']==previous['log_N0']==26820
    assert rough['GRH_at_31000']==strong
    assert current['presieve']==previous['presieve']
    assert current['binned_modulus_count_budget']=='1/118511'
    count=current['binned_support_count']
    for name,h in count['input_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h,name
    assert (count['sigma'],count['log_bin_width'],count['maximum_actual_degree'],
            count['label_limit'],count['nonempty_blocks'])==('0','1/100',22,8115,992)
    assert count['medium_prime_count']==348510
    assert arb(count['modulus_count_fraction_ball'])<arb(current['binned_modulus_count_budget'])
    assert arb(current['rough_density_ball'])<arb(current['rough_density_budget'])
    assert arb(current['margin_ball'])>arb('.00001') and arb(current['support_slack_ball'])>0
    assert arb(current['rough_floor_error_ratio_ball'])>0
    assert arb(current['rough_floor_error_ratio_ball'])<arb('1e-5302')
    assert arb(current['half_rectangle_error_ball'])<arb('1e-194')
    assert arb(current['finite_error_ball'])<arb('1e-1337')
    assert len(rough['candidate_rows'])==135
    assert rough_checks['analytic_join']['covered_log_interval']==current['valid_log_interval']
    assert rough_checks['analytic_join']['existing_half_line_start']==26820
    assert rough_checks['analytic_join']['inclusive_endpoints'] is True
    assert rough_checks['complete_integer_primes_independently_generated']==348513
    assert rough_checks['interval_checks'][0]['independent_endpoint_quantities']==12
    assert len(rough_checks['interval_checks'])==5
    assert all(arb(r['margin_ball'])>arb('.00001') for r in rough_checks['interval_checks'])
    assert binned_checks['binned_count_checks']['enumerated_subsets']==9216
    assert binned_checks['binned_count_checks']['exact_rational_histogram_coefficients_checked']==17172
    assert len(binned_checks['binned_count_checks']['configurations'])==48
    assert binned_checks['integer_roughness_checks']['exact_divisor_floor_terms_checked']==16384
    assert binned_checks['integer_roughness_checks']['integer_roughness_patterns_checked']==24316
    assert binned_checks['aggregate_prime_factor_checks']['modulus_patterns']==8608
    assert rough_checks['aggregate_prime_factor_modulus_patterns']==8608
    for label in ['sec:roughinterval','lem:binnedcount','lem:integerrough',
                  'eq:binnedcount','eq:integerrough','eq:roughsupport','eq:intervalplateau']:
        assert r'\label{'+label+'}' in tex
    front=tex.split('Let $B=4',1)[0]
    assert r'\alpha=2117/500000' in tex and r'\beta=209/100' in tex
    assert '17172' in tex and '24316' in tex and '118511' in tex
    rough_note='conditional_chen_rough_interval_20261004.md'
    for name in documents+[rough_note]:
        note=(ROOT/name).read_text(encoding='utf-8').split('**追加更正')[0]
        assert 'exp(24799)' in note and '0.0179' in note
        assert '24799≤L=log N≤26820' in note and '209/100' in note
        assert '1/118511' in note and '尚未拼接至有限验证' in note

    accepted=json.loads((ROOT/'conditional_chen_accepted_interval_certificate.json').read_text(encoding='utf-8'))
    accepted_checks=json.loads((ROOT/'conditional_chen_accepted_interval_checks.json').read_text(encoding='utf-8'))
    accepted_small=json.loads((ROOT/'conditional_chen_accepted_support_checks.json').read_text(encoding='utf-8'))
    standalone_accepted=json.loads((ROOT/'conditional_chen_accepted_support_certificate.json').read_text(encoding='utf-8'))
    for result in [accepted,accepted_checks,accepted_small,standalone_accepted]:
        if result is not standalone_accepted:assert result['status'].startswith('PASS')
        for name,h in result['input_sha256'].items():
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h,name
    preceding_rough=current;current=accepted['GRH_positivity'];strong=accepted['GRH_at_31000']
    assert (current['log_N0'],current['alpha'],current['beta'])==(24274,'4228/1000000','93/50')
    assert current['valid_log_interval']==[24274,24799]
    assert accepted['analytic_join_chain']==[24274,24799,26820,'infinity']
    assert accepted['analytic_tail_join_log_N']==preceding_rough['log_N0']
    assert current['presieve']==preceding_rough['presieve']
    assert current['binned_modulus_count_budget']=='1/992601'
    assert strong['log_N0']==31000 and strong['valid_log_interval']==[31000,None]
    assert strong['alpha']=='3737/1000000' and strong['beta']=='3'
    assert strong['claimed_margin']=='181/10000' and strong['fixed_beta_three_half_line']
    assert strong['binned_modulus_count_budget']=='1/1176659'
    assert len(accepted['candidate_rows'])==70 and len(accepted['stronger_candidates'])==9
    for row,degree,label_limit in [(current,22,8115),(strong,22,8211)]:
        count=row['binned_support_count']
        for name,h in count['input_sha256'].items():
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h,name
        assert count['sigma']=='0' and count['log_bin_width']=='1/100'
        assert count['medium_prime_count']==count['available_medium_prime_count']==348510
        assert count['maximum_actual_degree']==degree and count['label_limit']==label_limit
        assert count['nonempty_blocks']==992
        assert count['union_without_duplication'] and count['retains_all_preceding_checked_prefixes']
        assert count['complete_cutoffs_retained_after_deletion']
        assert len(count['state_count_balls'])==6
        integers=[arb(r['count_ball']).unique_fmpz() for r in count['state_count_balls']]
        assert None not in integers and sum(integers)==arb(count['count_upper_ball']).unique_fmpz()
        assert arb(count['modulus_count_fraction_ball'])<arb(row['binned_modulus_count_budget'])
        assert arb(row['rough_density_ball'])<arb(row['rough_density_budget'])
        assert arb(row['margin_ball'])>arb(row['claimed_margin'])
        assert arb(row['support_slack_ball'])>0 and arb(row['rough_floor_error_ratio_ball'])>0
        assert arb(row['half_rectangle_error_ball'])<arb(row['rectangle_display_bound'])
        assert arb(row['finite_error_ball'])<arb(row['finite_display_bound'])
    assert current['binned_support_count']==standalone_accepted
    assert accepted_small['enumerated_complete_subsets']==4608
    assert accepted_small['integer_histogram_coefficients_checked']==114876
    assert accepted_small['deleted_subset_patterns']==104976
    assert accepted_small['deleted_production_histograms_checked']==384
    assert accepted_checks['complete_integer_primes_independently_generated']==348513
    assert accepted_checks['aggregate_prime_factor_modulus_patterns']==8608
    assert accepted_checks['analytic_join']['covered_log_interval']==[24274,24799]
    assert accepted_checks['analytic_join']['prior_interval']==[24799,26820]
    assert accepted_checks['analytic_join']['existing_half_line_start']==26820
    assert accepted_checks['analytic_join']['inclusive_endpoints']
    assert len(accepted_checks['interval_checks'])==5
    assert len(accepted_checks['strong_half_line_checks'])==3
    for group in ['interval_checks','strong_half_line_checks']:
        assert accepted_checks[group][0]['independent_endpoint_quantities']==12
    for label in ['sec:acceptedinterval','lem:acceptedcount','eq:acceptedoperator',
                  'eq:acceptedblock','eq:acceptedfraction']:
        assert r'\label{'+label+'}' in tex
    front=tex.split('Let $B=4',1)[0]
    assert r'N\ge e^{24274}' in front and '$0.0181$' in front
    assert '24799' not in front and '26820' not in front
    assert '655257348020585286337311701158' in tex and '1176659' in tex
    accepted_note='conditional_chen_accepted_interval_20261004.md'
    for name in documents+[accepted_note]:
        note=(ROOT/name).read_text(encoding='utf-8').split('**追加更正')[0]
        assert 'exp(24274)' in note and '0.0181' in note
        assert '24274≤L=log N≤24799' in note and '93/50' in note
        assert '1/992601' in note and '六' in note and '尚未拼接至有限验证' in note

    local_bridge=json.loads((ROOT/'conditional_chen_full_rosser_bridge_checks.json').read_text(encoding='utf-8'))
    assert local_bridge['status'].startswith('PASS')
    assert local_bridge['complete_odd_prime_count']==148932
    assert local_bridge['local_scalar_reproductions'][-1]=={'w':2,'u':'61/20','C_budget':3489,'kappa':'3/25'}
    small_local=local_bridge['small_exact_checks']
    assert small_local['enumerated_complete_subsets']==9216
    assert small_local['weighted_histogram_coefficients_checked']==152316
    assert small_local['deleted_production_histograms_checked']==768
    assert small_local['pointwise_divisor_patterns_and_first_failure_identity_cases']==26244
    assert len(local_bridge['centered_pi_input_counterexamples'])==3
    for result_name in ['conditional_chen_full_rosser_bridge_probe.json',
                        'conditional_chen_weighted_rosser_support_certificate.json',
                        'conditional_chen_full_rosser_w2_probe.json',
                        'conditional_chen_full_rosser_bridge_checks.json']:
        result=json.loads((ROOT/result_name).read_text(encoding='utf-8'))
        for name,h in result['input_sha256'].items():
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==h,name
    for row in local_bridge['centered_pi_input_counterexamples']:
        assert row['status'].startswith('REFUTED')
        assert arb(row['necessary_A_epsilon_1_1000_lower_ball'])>arb(row['candidate_A_capacity_ball'])
        assert int(row['N'])>=4000000000000000000
    local_note='conditional_chen_full_rosser_bridge_20261004.md'
    for name in documents+[local_note]:
        note=(ROOT/name).read_text(encoding='utf-8').split('**追加更正')[0]
        assert '3489' in note and '0.255103' in note and '0.643081' in note
        assert '未证明' in note and 'exp(24274)' in note

    inputs = scripts + documents + [
        tex_path.name,
        local_note,
        'conditional_chen_full_rosser_bridge_probe.json',
        'conditional_chen_weighted_rosser_support_certificate.json',
        'conditional_chen_full_rosser_w2_probe.json',
        'conditional_chen_full_rosser_bridge_checks.json',
        weighted_note,
        "conditional_chen_composite_certificate.json",
        "conditional_chen_prefix_band_certificate.json",
        "conditional_chen_prefix_band_checks.json",
        "conditional_chen_prefix_band_eh_search.json",
        "conditional_chen_sharp_support_certificate.json",
        "conditional_chen_sharp_endpoint_certificate.json",
        "conditional_chen_distribution_tradeoff_certificate.json",
        "conditional_chen_weighted_distribution_certificate.json",
        "conditional_chen_dense_band_grh_certificate.json",
        "conditional_chen_dense_band_grh_finite.json",
        "conditional_chen_sharp_grh_parameter_probe.json",
        "conditional_chen_degreewise_rankin_certificate.json",
        "conditional_chen_rosser_checks.json",
        "conditional_chen_refined_eh_search.json",
        "conditional_chen_degreewise_eh_search.json",
        "conditional_chen_degreewise_rankin_checks.json",
        "conditional_chen_quadratic_rankin_certificate.json",
        "conditional_chen_quadratic_rankin_checks.json",
        "conditional_chen_quadratic_eh_search.json",
        "conditional_chen_uniform_product_certificate.json",
        "conditional_chen_uniform_product_checks.json",
        prefix_note,
        "conditional_chen_prefix_block_certificate.json",
        "conditional_chen_prefix_block_checks.json",
        "conditional_chen_prefix_block_grh_finite.json",
        "conditional_chen_prefix_block_endpoint_certificate.json",
        "conditional_chen_prefix_block_refinement.json",
        "conditional_chen_sparse_modulus_certificate.json",
        "conditional_chen_prefix_block_endpoint_checks.json",
        rough_note,
        "conditional_chen_binned_support_checks.json",
        "conditional_chen_binned_support_certificate.json",
        "conditional_chen_rough_interval_certificate.json",
        "conditional_chen_rough_interval_checks.json",
        accepted_note,
        "conditional_chen_accepted_support_certificate.json",
        "conditional_chen_accepted_support_checks.json",
        "conditional_chen_accepted_interval_certificate.json",
        "conditional_chen_accepted_interval_checks.json",
    ]
    return {
        "status": "PASS: source structure and artifact agreement only",
        "precision_bits": ctx.prec,
        "labels": len(labels),
        "reference_occurrences": len(references),
        "bibliography_keys": len(bibkeys),
        "GRH_log_N0": grh["log_N0"],
        "GRH_positivity_log_N0": current["log_N0"],
        "GRH_dense_previous_log_N0": dense["GRH_positivity"]["log_N0"],
        "GRH_previous_positivity_log_N0": positive["log_N0"],
        "GRH_dense_beta": [dense["GRH_positivity"]["beta"],dense["GRH_at_31000"]["beta"]],
        "GRH_stronger_claim": strong["claimed_margin"],
        "GRH_current_beta": current['beta'],
        "GRH_binned_medium_count_budget": current['binned_modulus_count_budget'],
        "GRH_rough_density_budget": current['rough_density_budget'],
        "GRH_analytic_tail_join": current['valid_log_interval'][1],
        "GRH_prior_sparse_start": previous['log_N0'],
        "GRH_prior_rough_start": preceding_rough['log_N0'],
        "GRH_analytic_join_chain": accepted['analytic_join_chain'],
        "GRH_support_count_states": 6,
        "local_full_rosser_signed_budget": 3489,
        "local_full_rosser_unproved_interval": [4000000000000000000,8000000000000000000],
        "refuted_local_centered_pi_inputs": 3,
        "prefix_distribution_log_starts": expected_prefix,
        "EH_log_starts": [786,212,166,142],
        "EH_positivity_log_starts": [774,210,165,141],
        "complete_product_prime_count": uniform["complete_odd_prime_count"],
        "finite_product_checks": uniform_checks,
        "distribution_tradeoff_witnesses": expected_tradeoff,
        "weighted_distribution_log_starts": expected_weighted,
        "signed_remainder_exact_small_cases": small,
        "sha256": {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
                   for name in inputs},
        "not_checked": [
            "general analytic lemmas and external theorem validity",
            "independent research review",
            "TeX compilation and PDF rendering",
            "connection to the finite verification range",
        ],
    }


if __name__ == "__main__":
    result = audit()
    Path(__file__).with_suffix(".json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
