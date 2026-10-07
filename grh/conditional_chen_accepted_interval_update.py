"""Write the six-state support proof and its separately checked GRH endpoints."""
import json
from fractions import Fraction
from pathlib import Path
from flint import arb,ctx
from conditional_chen_prefix_block_update import bound,load_checked

ROOT=Path(__file__).resolve().parent
ctx.prec=192


def section(report,checks):
    weak,strong=report['GRH_positivity'],report['GRH_at_31000']
    start,anchor=weak['valid_log_interval']
    assert checks['analytic_join']['covered_log_interval']==[start,anchor]
    assert checks['strong_half_line']['claim']==strong['claimed_margin']
    table=[]
    for label,key,upper in [('$M_0$','lower_ball',False),('$M_1$','half_upper_ball',True),
        ('$M_2$','half_switch_ball',True),('$E_{\\rm AP}$','AP_budget_ball',True),
        ('Support slack','support_slack_ball',False),('Surplus','margin_ball',False)]:
        table.append(label+' & $'+('<' if upper else '>')+bound(weak[key],15,upper)+'$ & $'
                     +('<' if upper else '>')+bound(strong[key],15,upper)+'$'+r'\\')
    text=r'''% BEGIN CURRENT ACCEPTED INTERVAL
\subsection{The union of the accepted medium supports}\label{sec:acceptedinterval}
The preceding logarithmic count includes medium products that neither
Rosser sign accepts. The following positive count keeps every preceding
acceptance condition. It bounds the union of the two supports once,
without changing the sieve coefficients in Lemma \ref{lem:composite}.

\begin{lemma}[Six-state union count]\label{lem:acceptedcount}
Let $D_0=T^u$, $u>3$, and let $H_0=\kappa D_0>1$ bound both pure-medium
Rosser supports as in Lemma \ref{lem:sharpsupport}. Use the complete
medium-prime set, the degree bound $K$, labels $b(p)$ and cutoff $J$
of Lemma \ref{lem:binnedcount}, with $0<h<\log\min\mathcal M$.
Process the complete blocks in decreasing label order. In block
$\mathcal M_b$, put
\[
c_b=\left\lceil\frac{\log D_0-3\log\min\mathcal M_b}{h}\right\rceil,
\qquad m_b=|\mathcal M_b|.
\]
For a positive polynomial $F=\sum_B F_Bv^B$, define the positive operators
\[
\mathcal L_cF=\sum_{B<c}F_Bv^B,\quad
\mathcal H_cF=\sum_{B\ge c}F_Bv^B,\quad
\mathcal S_bF=v^bF\pmod{v^J}.
\]
Track the vector $X=(U_e,U_o,L_e,L_o,B_e,B_o)^\mathsf T$.
The first two coordinates mean only the upper sign is alive, the next
two mean only the lower sign is alive, and the final two mean both signs
are alive; subscripts record the parity of the number of selected primes.
Start with $B_e=1$ and every other coordinate zero. Within the block,
the append operator is
\begin{equation}\label{eq:acceptedoperator}
\mathcal A_bX=
\begin{pmatrix}
\mathcal S_b(U_o+\mathcal H_{c_b}B_o)\\
\mathcal S_b\mathcal L_{c_b}U_e\\
\mathcal S_b\mathcal L_{c_b}L_o\\
\mathcal S_b(L_e+\mathcal H_{c_b}B_e)\\
\mathcal S_b\mathcal L_{c_b}B_o\\
\mathcal S_b\mathcal L_{c_b}B_e
\end{pmatrix}.
\end{equation}
Apply the complete-block update
\begin{equation}\label{eq:acceptedblock}
X\ \longleftarrow\ \sum_{j=0}^{\min(K,m_b)}
                  \binom{m_b}{j}\mathcal A_b^jX.
\end{equation}
Then $A_{\rm acc}:=\sum_iX_i(1)$ bounds the number of distinct medium
products accepted by at least one Rosser sign. The same complete-set
bound holds after arbitrary prime deletion, retaining the complete
$K,J,c_b$ and replacing only $m_b$ by the number of available primes.
\end{lemma}
\begin{proof}
At new selected length $r$, the upper sign checks odd $r$ and the lower
sign checks even $r$. If $d_{\rm old}$ is the preceding product, a checked
append of prime $p$ passes only if $d_{\rm old}p^3<D_0$.
Its preceding label $B$ has $hB\le\log d_{\rm old}$, so actual passing
implies
\[
B<\frac{\log D_0-3\log p}{h}
 \le\frac{\log D_0-3\log\min\mathcal M_b}{h},\qquad B<c_b.
\]
Consequently, a true surviving sign is never killed by the relaxed
cutoff. Unchecked appends preserve that sign. For both signs alive,
exactly one sign is checked at each append: a high-label append leaves
the other sign alone. This gives \eqref{eq:acceptedoperator}.
Each label and flag state has at most one destination; the low and high
branches are disjoint. Induction shows that the relaxed alive set
contains the true alive set, and every retained relaxed subset occurs
in exactly one of the six exclusive states. Subsets with no surviving sign may
be dropped permanently.

All primes in a block have the same label and relaxed cutoff. Choosing
$j$ distinct primes therefore applies the same operator $j$ times, with
exact multiplicity $\binom{m_b}{j}$. Their actual decreasing order does
not affect this relaxed update. Every true final accepted product is
below $H_0$, hence contains at most $K$ primes and has label below $J$.
Every preceding label is no larger, so neither the block degree cutoff
nor intermediate truncation can remove that true contribution.
Excessive total degrees and possible passes retained by rounding are
additional positive terms. The resulting sum bounds the union, without
counting a subset twice when both signs survive.

After deletion, every surviving prime still belongs to its complete
block. The same complete minimum supplies a valid necessary passing
condition; the complete $K,J$ still cover every accepted product.
The operators have nonnegative coefficients and
$\binom{m_b'}{j}\le\binom{m_b}{j}$ for $m_b'\le m_b$.
Induction over the blocks therefore gives coefficientwise domination
in every state. No subtraction of uncertain intervals is required.
\end{proof}

Replace $A_\sigma$ by $A_{\rm acc}$ in \eqref{eq:roughsupport} and set
\begin{equation}\label{eq:acceptedfraction}
\nu_{\rm acc}=\frac{2^{\omega(Q_*)}A_{\rm acc}}{Q_*H_0}.
\end{equation}
The composite coefficient support is bounded by the union of the two
medium supports and a $T$-rough large factor, so the same proof gives
$\nu_{\rm acc}(a_T+b_T/D)H$ for both sieve signs.
For the whole $A_q$ sum, absorb $q$ into the large factor as in Section
\ref{sec:roughinterval}. The medium product is still in this same union;
the unique prime at least $z$ makes $(q,d)\mapsto qd$ injective over
the entire sum. Thus the new count also bounds its aggregate image.
The integer roughness boundary, prime-modulus centering, and bilinear
charges are unchanged.

For $T=5000000$, $w_0=5$, $u=87/16$, $\kappa=11/169$ and count width
$h=1/100$, the complete count has $K=22$, $J=8115$ and 992 blocks.
It gives
\[
A_{\rm acc}=655257348020585286337311701158,
\qquad \nu_{\rm acc}<1.007455\cdot10^{-6}<1/992601.
\]
The preceding unfiltered $A_0$ is more than $8.37$ times as large.
The gap budgets still use their independently certified width $1/200$.

\paragraph{Coverage of all larger integers.}
For the lower starting point use $L_0=__START__$, $L_*=__ANCHOR__$,
$\alpha=__WA__$ and $\beta=__WB__$; the medium family is $u=87/16$.
Equation \eqref{eq:intervalplateau}, with the new
$\nu_*:=1/992601$ and the same $a_*:=36411450583/10^{12}$,
supplies a constant leading AP expense throughout $[L_0,L_*]$.
The support condition $\alpha L-\beta\log L>\log(Q_*H_0)$,
$H\ge10^9$, $H\le\sqrt N$, $T<z$, $3<s<4$ and $k>1$ hold at
$L_0$ and persist. In particular $\alpha-\beta/L_0>0$ and
$1/2-\beta/L_0>0$. The finite charge $b_Te^{-tL_0}$ is positive
and below $10^{-5189}$; it has not been discarded.
The fixed-parameter monotonicity argument in Section \ref{sec:roughinterval}
proves the required surplus over this entire closed interval.
The preceding interval $[24799,26820]$ and the preceding half-line
$[26820,\infty)$ use the same $10^{-5}$ claim and normalization.
Together these inclusive ranges prove the first assertion of Theorem
\ref{thm:grh} for every $L\ge__START__$.

For the stronger claim, use $u=11/2$, the same $T,w_0,\kappa$,
$\alpha=__SA__$, $\beta=3$ and $L_0=31000$.
Its gap budgets are $\eta_+=1/3231$, $\eta_-=1/3242$, with
$\epsilon=1/41219$. The accepted count gives
$\nu_{\rm acc}<8.498637\cdot10^{-7}<\nu_*:=1/1176659$, and
$a_T<a_*:=36409264976/10^{12}$.
Now $L^{3-\beta}=1$, so the leading AP expense
\[
\frac{3\nu_*}{2U_*}(a_*+b_Te^{-tL})P(L)
\]
decreases on the whole half-line. Its value at $31000$ is a valid
constant upper bound for every $L\ge31000$.
All other losses decrease and the support slack increases.
This separate half-line argument proves the stronger coefficient
$__SC__$ of Theorem \ref{thm:grh}; it does not transfer a larger
claim through an earlier certificate with a smaller surplus.
The bilinear estimate continues to use the maximum modulus $H$.

\begin{center}
\begin{tabular}{lrr}\toprule
Term & $L_0=__START__$ & $L_0=31000$\\\midrule
$\alpha$ & $__WA__$ & $__SA__$\\
$\beta$ & $__WB__$ & $3$\\
$s$ & $__WS__$ & $__SS__$\\
__TABLE__
$E_{\rm rect}$ & $<10^{-188}$ & $<10^{-250}$\\
$E_{\rm fin}$ & $<10^{-1300}$ & $<10^{-1673}$\\\bottomrule
\end{tabular}
\end{center}
The separate-formula program reproduces 12 quantities for each new
endpoint, regenerates all 348513 integer primes through $T$, and
rechecks the earlier interval and tails. Exact small-set checks cover
4608 complete subsets, 114876 integer histogram coefficients,
104976 deleted-subset patterns and 384 deleted production histograms.
The earlier 8608 aggregate prime-factor patterns remain applicable.
The __CANDIDATES__ weak-endpoint candidates and nine stronger-count
candidates are a finite search, not a proof of global optimality.
These checks supplement the general proof above; they do not replace
the external analytic inputs or independent research review.
No distribution hypothesis beyond GRH has been introduced, and this
analytic interval join still does not reach finite verification.
% END CURRENT ACCEPTED INTERVAL
'''
    values={'START':start,'ANCHOR':anchor,'WA':str(Fraction(weak['alpha'])),
        'WB':weak['beta'],'SA':str(Fraction(strong['alpha'])),
        'SC':str(float(Fraction(strong['claimed_margin']))),
        'WS':f"{float(arb(weak['s_ball'])):.6f}",
        'SS':f"{float(arb(strong['s_ball'])):.6f}",
        'TABLE':'\n'.join(table),'CANDIDATES':len(report['candidate_rows'])}
    for k,v in values.items():text=text.replace('__'+k+'__',str(v))
    assert '__' not in text
    return text


def main():
    report=load_checked('conditional_chen_accepted_interval_certificate.json')
    checks=load_checked('conditional_chen_accepted_interval_checks.json')
    weak,strong=report['GRH_positivity'],report['GRH_at_31000']
    path=ROOT/'conditional_chen_grh_bridge.tex';tex=path.read_text(encoding='utf-8')
    front,rest=tex.split('Let $B=4',1)
    start=weak['log_N0'];coefficient=str(float(Fraction(strong['claimed_margin'])))
    front=front.replace('24799',str(start)).replace('0.0179',coefficient)
    front=front.replace('10.118558608',bound(arb(start).log(),9))
    digits=int((arb(start)/arb(10).log()).floor().unique_fmpz())+1
    front=front.replace('$10771$',f'${digits}$')
    front=front.replace('Positive prefix blocks retain all preceding checks. Logarithmic counts\n'
        'and an elementary rough-integer bound reduce the AP expense; a finite\n'
        'analytic interval joins an existing GRH tail.',
        'Positive prefix blocks retain preceding checks. A six-state count\n'
        'bounds the union of accepted supports; an integer roughness bound\n'
        'and a finite analytic interval reduce the GRH threshold.')
    tex=front+'Let $B=4'+rest
    tex=tex.replace('and the stronger coefficient $0.0179$ in Theorem \\ref{thm:grh}.',
        'and the preceding stronger coefficient $0.0179$.')
    tex=tex.replace('the lower starting point of Theorem \\ref{thm:grh} using GRH alone.\n'
        'Its stronger $0.0179$ coefficient at $L\\ge31000$ retains the preceding\n'
        'half-line certificate.',
        'the preceding starting point $e^{24799}$ using GRH alone.\n'
        'The preceding stronger $0.0179$ coefficient at $L\\ge31000$ retains\n'
        'its earlier half-line certificate.')
    proof=section(report,checks)
    if '% BEGIN CURRENT ACCEPTED INTERVAL' in tex:
        a,b=tex.split('% BEGIN CURRENT ACCEPTED INTERVAL',1)
        _,b=b.split('% END CURRENT ACCEPTED INTERVAL',1)
        tex=a+proof+b.lstrip('\n')
    else:
        tex=tex.replace('\\section{A parameterized EH-style reduction}',
            proof+'\n\\section{A parameterized EH-style reduction}',1)
    ledger=r'''% BEGIN ACCEPTED INTERVAL LEDGER
For the accepted-support union and current endpoints, run:
\begin{quote}\small
\path{conditional_chen_accepted_support_checks.py}\\
\path{conditional_chen_accepted_interval_certificate.py}\\
\path{conditional_chen_accepted_interval_checks.py}
\end{quote}
The counter uses six positive exclusive flag/parity states and exact
binomial block multiplicities. Its complete cutoffs remain fixed after
deletion. Counts are generated completely or reused only after generator
hash checks. The scalar checker uses its own endpoint formulas and
integer prime sieve, checks exact sums of the six state counts, and
rechecks the preceding interval and tails. The weak interval covers
$[24274,24799]$ and joins the earlier $[24799,26820]$ and GRH tail;
the stronger coefficient uses its own fixed-$\beta=3$ half-line.
Continuous coverage follows from Section \ref{sec:acceptedinterval},
not from finitely many sample points. Neither a finite-verification
splice nor independent research review is supplied by these programs.
% END ACCEPTED INTERVAL LEDGER
'''
    if '% BEGIN ACCEPTED INTERVAL LEDGER' in tex:
        a,b=tex.split('% BEGIN ACCEPTED INTERVAL LEDGER',1)
        _,b=b.split('% END ACCEPTED INTERVAL LEDGER',1)
        tex=a+ledger+b.lstrip('\n')
    else:
        tex=tex.replace('\\section*{Reproducibility and dependency ledger}',
            '\\section*{Reproducibility and dependency ledger}\n'+ledger,1)
    assert tex.count('\\label{lem:acceptedcount}')==1
    path.write_text(tex,encoding='utf-8',newline='\n')
    print(json.dumps({'GRH_log_N0':start,'stronger_coefficient':coefficient,'decimal_digits':digits}))


if __name__=='__main__':main()
