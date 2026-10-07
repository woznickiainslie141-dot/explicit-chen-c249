# Explicit Chen bounds: current research drafts

**最新整理：2026-10-07。** 本仓库公开分享 GPT 生成并协助内部核查的研究稿。首页以最新七分之八零点输入路线为主稿，区分 GRH 与历史无条件路线。原 c249 已移到 [历史目录](archive/README.md)，仓库地址沿用以保持旧链接有效。

For even N, D_{1,2}(N) counts primes p < N for which N-p >= 2 and Ω(N-p) <= 2, counting prime factors with multiplicity. Put

$$
U_N=2e^{-\gamma}\prod_{p>2}\left(1-\frac1{(p-1)^2}\right)
\prod_{\substack{p>2\\p\mid N}}\frac{p-1}{p-2}.
$$

## Current results and assumptions

| Route | Starting point for every even N | Draft's conclusion | Input |
|---|---|---|---|
| **Latest uniform 7/8 route** | N >= exp(51046), log log N ≈ 10.84048 | D_{1,2}(N) > 9×10^-7 U_N N/(log N)^2 | Uniform zero-free Re(s)>7/8 for **all Dirichlet L-functions**, all moduli and heights, with the principal pole allowed. No GRH estimate is used. |
| **GRH companion** | N >= exp(24274), log log N ≈ 10.09716 | D_{1,2}(N) > 10^-5 U_N N/(log N)^2 | GRH for Dirichlet L-functions, including RH for zeta. The coefficient improves to 0.0181 from exp(31000). |
| **Historical c=24.9** | N >= exp(exp(24.9)) | D_{1,2}(N)>0 | Older manuscript's unconditional Wu–BJS argument; preserved for provenance. |

The 7/8 assertion is exactly [Theorem 1.1 in the external OpenAI manuscript, fixed revision](https://github.com/openai/math/blob/adc7f1241b42e322a6451854ab7e4b4c146bf78a/preprints/The-Quasi-Riemann-Hypothesis-September-30-2026/build/paper.tex). **Adopting that external theorem removes the zero-free hypothesis of our reduction.** This repository does not reprove that theorem or verify its complete Lean dependency chain. The reduction under its explicit input remains separate from the external justification. A zeta-only QRH assumption is insufficient here.

These are the strongest stated results in the selected project manuscripts under their respective assumptions, not a claim of global optimality, an independently verified published record, or a completed bridge to the finite Goldbach verification range. exp(51046) is **not** exp(exp(51046)).

## Materials and reproduction

- **Main manuscript:** [QRH-route PDF](qrh/output/pdf/qrh_chen_51046.pdf), [LaTeX](qrh/qrh_chen_51046.tex), [final audit](qrh/qrh_final_audit.md), [reproduction guide](qrh/qrh_REPRODUCE.md).
- [GRH PDF](grh/output/pdf/conditional_chen_grh_bridge.pdf), [LaTeX](grh/conditional_chen_grh_bridge.tex), [certificate guide](grh/README.md).
- [Historical c249 package](archive/c249/README_REPRODUCE.md). Its historical private-review and authorship labels are superseded for the present publication by this README and [NOTICE](NOTICE.md).
- [Publication checks](PUBLICATION_CHECKS.json): integrity, current finite runs and PDF checks; these do not prove the analytic arguments.

With Python 3.12 and python-flint 0.9.0, run in the qrh directory:

~~~text
python -m pip install -r requirements.txt
python qrh_chen_reproduce.py
python qrh_chen_interface_checks.py
python qrh_chen_positivity_search.py
~~~

The first program regenerates the prime masses, ordered product and interval integrals without old JSON. The next two check interfaces and the declared finite search; that search does not establish global optimality. Analytic proofs and half-line propagation are in the manuscript. Compile each standalone TeX twice with pdflatex; third-party source snapshots are provenance, not TeX includes.

## AI use and review status

GPT generated the proof strategies and mathematical drafts and assisted with internal audits and certificate code. Cheng Huang publishes and maintains the materials as distributor, rather than claiming human authorship of the mathematical arguments. The QRH manuscript's personal author field is empty; the GRH manuscript uses the generic “Research manuscript”.

These are research drafts for discussion and independent review. Internal audits, interval computations and compilation do not constitute independent peer review or complete formal verification. Verification of the external 7/8 theorem is outside these certificates. No comprehensive novelty or priority determination is claimed.

## Version history

- **v2026.10.07-current-drafts:** current QRH and GRH manuscripts replace the old c249 homepage.
- **c249-before-20261007:** Git tag preserving the prior repository tree.

The historical repository name is retained for link continuity. Current statements are those linked above. See [NOTICE](NOTICE.md) for provenance and third-party rights.
