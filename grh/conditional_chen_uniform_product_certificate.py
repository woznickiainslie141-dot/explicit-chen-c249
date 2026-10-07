"""Finite unconditional reciprocal-prime envelope, with an analytic tail.

The finite calculation checks every prime jump and every initial endpoint.
Between primes E(x)=sum_(p<x)1/p-log(log(x))-M decreases. The remaining
half-line uses the cited BJS/JS upper envelopes, not sampled points.
This certificate does not prove the cited far-range estimates.
"""
import json
import math
import time
from pathlib import Path
from flint import arb, ctx
from conditional_chen_rankin_certificate import odd_primes_to


def merge_ordered_statistics(left, right):
    """Enclose peak, trough and ordered drawdown on adjacent intervals."""
    p1, q1, a1 = left
    p2, q2, a2 = right
    return p1.max(p2), q1.min(q2), a1.max(a2).max(p1-q2)


def positive_decimal_upper(value, digits=12):
    """A float proposes a rational endpoint; Arb alone accepts it."""
    scale = 10**digits
    numerator = math.ceil(float(value)*scale)
    while not value < arb(numerator)/scale:
        numerator += 1
    return numerator, scale


def certify_uniform_products(
    cutoffs=(10000, 22000, 30000, 50000, 70000, 100000, 130000,
             200000, 500000, 1000000, 2000000, 3000000,
             5000000, 10000000, 20000000, 30000000),
    S=100000000, M_cutoff=10000000, progress=True,
):
    assert ctx.prec >= 192
    assert list(cutoffs) == sorted(set(cutoffs))
    assert 10000 <= cutoffs[0] <= cutoffs[-1] < S <= 10**12
    assert 2 <= M_cutoff <= S
    assert all(T % 2 == 0 for T in cutoffs) and S % 2 == 0
    start = time.monotonic()
    primes = odd_primes_to(S)
    if progress:
        print(f"Uniform envelope: generated {len(primes)} odd primes through {S}",
              flush=True)
    # M=gamma+sum_p(log(1-1/p)+1/p).
    # For n>M_cutoff, -log(1-1/n)-1/n <= 1/(2n(n-1)).
    # Summing over every integer bounds the omitted negative prime tail.
    M_finite = arb.const_euler()
    M_finite += (arb(1)/2).log()+arb(1)/2
    log_C_finite = arb(0)
    for p in primes:
        if p > M_cutoff:
            break
        inverse = arb(1)/p
        M_finite += (1-inverse).log()+inverse
        log_C_finite += (1-arb(1)/(p-1)**2).log()
    M_lower = M_finite-arb(1)/(2*M_cutoff)
    assert M_lower > arb(".2614971")
    # The product over every integer n>M_cutoff telescopes to
    # (M_cutoff-1)/M_cutoff. The prime-only tail is larger.
    log_D_infty_lower = arb(2).log()-arb.const_euler()+log_C_finite
    log_D_infty_lower += (1-arb(1)/M_cutoff).log()

    reciprocal = arb(1)/2
    log_product = arb(0)
    segment = -1
    maxima = []
    peaks, troughs, drawdowns = [], [], []
    checked_jumps = 0
    for i, p in enumerate(primes):
        while segment+1 < len(cutoffs) and p >= cutoffs[segment+1]:
            segment += 1
            T = cutoffs[segment]
            # T is even, hence not prime: strict and inclusive counts agree.
            endpoint = reciprocal-arb(T).log().log()-M_lower
            maxima.append(endpoint)
            logD = arb(T).log().log()+log_product
            peaks.append(logD)
            troughs.append(logD)
            drawdowns.append(arb(0))
        reciprocal += arb(1)/p
        log_factor = (1-arb(1)/(p-1)).log()
        if segment >= 0:
            postjump = reciprocal-arb(p).log().log()-M_lower
            maxima[segment] = maxima[segment].max(postjump)
            left_logD = arb(p).log().log()+log_product
            peaks[segment] = peaks[segment].max(left_logD)
            right_logD = left_logD+log_factor
            troughs[segment] = troughs[segment].min(right_logD)
            drawdowns[segment] = drawdowns[segment].max(
                peaks[segment]-right_logD)
            checked_jumps += 1
        log_product += log_factor
        if progress and (i+1) % 1000000 == 0:
            print(f"Uniform envelope: {i+1}/{len(primes)} primes, "
                  f"{time.monotonic()-start:.1f}s", flush=True)
    assert segment == len(cutoffs)-1
    final_logD = arb(S).log().log()+log_product
    peaks[-1] = peaks[-1].max(final_logD)
    troughs[-1] = troughs[-1].min(final_logD)
    drawdowns[-1] = drawdowns[-1].max(peaks[-1]-final_logD)
    # For v>=S, monotonicity of both published envelopes gives this bound.
    far = (2/(arb(S).sqrt()*arb(S).log())).max(
        arb(".000006836")/arb(10**12).log())
    running = far
    far_logK = far+arb(".000002964")/arb(S).log()
    far_logK += arb(1)/(S-1)+arb(".54")/S
    suffix_peak = suffix_trough = suffix_drawdown = None
    rows = []
    for j in reversed(range(len(cutoffs))):
        T, finite_segment = cutoffs[j], maxima[j]
        running = running.max(finite_segment)
        numerator, denominator = positive_decimal_upper(running)
        upper = arb(numerator)/denominator
        assert running < upper
        logK = upper+arb(".000002964")/arb(T).log()
        logK += arb(1)/(T-1)+arb(".54")/T
        epsilon = logK.exp()-1
        previous = (2/(arb(T).sqrt()*arb(T).log())).max(
            arb(".000006836")/arb(10**12).log())
        previous += arb(".000002964")/arb(T).log()
        previous += arb(1)/(T-1)+arb(".54")/T
        if suffix_peak is None:
            suffix_peak, suffix_trough, suffix_drawdown = (
                peaks[j], troughs[j], drawdowns[j])
        else:
            suffix_peak, suffix_trough, suffix_drawdown = merge_ordered_statistics(
                (peaks[j], troughs[j], drawdowns[j]),
                (suffix_peak, suffix_trough, suffix_drawdown))
        cross = suffix_peak-log_D_infty_lower+far
        combined_logK = suffix_drawdown.max(cross).max(far_logK)
        drawdown_numerator, drawdown_denominator = positive_decimal_upper(
            combined_logK)
        accepted_logK = arb(drawdown_numerator)/drawdown_denominator
        assert combined_logK < accepted_logK
        rows.append({
            "T": T,
            "reciprocal_upper_numerator": numerator,
            "reciprocal_upper_denominator": denominator,
            "reciprocal_upper_ball": str(upper),
            "all_v_ge_T_upper_enclosure_ball": str(running),
            "uniform_epsilon_ball": str(epsilon),
            "previous_uniform_epsilon_ball": str(previous.exp()-1),
            "finite_ordered_log_drawdown_ball": str(suffix_drawdown),
            "finite_to_far_log_ratio_ball": str(cross),
            "far_to_far_log_ratio_ball": str(far_logK),
            "ordered_logK_numerator": drawdown_numerator,
            "ordered_logK_denominator": drawdown_denominator,
            "ordered_uniform_epsilon_ball": str(accepted_logK.exp()-1),
        })
    rows.reverse()
    return {
        "status": "strict finite unconditional scalar certificate; BJS/JS far-range bounds remain cited dependencies",
        "precision_bits": ctx.prec,
        "S": S, "M_cutoff": M_cutoff,
        "M_finite_ball": str(M_finite), "M_lower_ball": str(M_lower),
        "log_D_infty_lower_ball": str(log_D_infty_lower),
        "complete_odd_prime_count": len(primes),
        "checked_prime_jumps": checked_jumps,
        "initial_endpoints_checked": list(cutoffs),
        "far_range_upper_ball": str(far),
        "rows": rows,
        "endpoint_convention": "E uses p<x; each post-jump supremum uses p<=x and covers the right limit; composite cutoff endpoints are checked separately",
        "continuum_coverage": "E decreases between prime jumps; no interior sampling",
        "product_formula": "K_T=exp(bbar_T+2.964e-6/log(T)+1/(T-1)+.54/T)",
        "ordered_product_formula": "log K_T <= max(finite ordered log D drawdown, max_(T<=u<=S) log D(u)-log D_infty_lower+b_S, original analytic log K_S)",
        "assumptions": "unconditional; no RH, GRH, or EH",
        "elapsed_seconds": round(time.monotonic()-start, 3),
    }


if __name__ == "__main__":
    ctx.prec = 192
    result = certify_uniform_products()
    Path(__file__).with_suffix(".json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ["M_finite_ball", "rows"]}, indent=2))
    print(json.dumps(result["rows"], indent=2))
