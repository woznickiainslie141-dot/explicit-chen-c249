"""Insert general support-count proofs and the checked finite-interval join."""
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from flint import arb,ctx
from conditional_chen_prefix_block_update import bound,load_checked

ROOT=Path(__file__).resolve().parent
ctx.prec=192


def section(report,checks):
    row=report['GRH_positivity'];count=row['binned_support_count']
    start,anchor=row['valid_log_interval']
    assert checks['analytic_join']['covered_log_interval']==[start,anchor]
    assert checks['interval_checks'][0]['independent_endpoint_quantities']==12
    table=[]
    for label,key,upper in [('$M_0$','lower_ball',False),('$M_1$','half_upper_ball',True),
        ('$M_2$','half_switch_ball',True),('$E_{\\rm AP}$','AP_budget_ball',True),
        ('Support slack','support_slack_ball',False),('Surplus','margin_ball',False)]:
        table.append(label+' & $'+('<' if upper else '>')+bound(row[key],15,upper)+'$'+r'\\')
    text=r'''% BEGIN CURRENT ROUGH INTERVAL
\subsection{Logarithmic support counts and rough large factors}\label{sec:roughinterval}
Two elementary count bounds improve the previous support estimate.
They introduce no distribution hypothesis. Throughout this subsection,
the actual sieve coefficients remain those of Lemma \ref{lem:composite}.

\begin{lemma}[Positive logarithmic-bin support count]\label{lem:binnedcount}
Let $\mathcal M$ be the complete medium-prime set, $H_0=\kappa T^u>1$,
and let $K$ be the largest number of its smallest primes whose product
is below $H_0$. Choose $0<h<\log\min\mathcal M$, put
$b(p)=\lfloor\log p/h\rfloor$ and $J=\lceil\log H_0/h\rceil$.
Partition $\mathcal M$ into blocks $\mathcal M_b$ of common label $b$.
For $\sigma\ge0$, define positive coefficients by
\[
\prod_{p\in\mathcal M_b}(1+p^{-\sigma}a)=\sum_j e_{b,j}a^j,
\qquad
P(v)=\prod_b\left(\sum_{0\le j\le K}e_{b,j}v^{jb}\right)\pmod{v^J}.
\]
Write $P_B=[v^B]P$. The number of square-free products of primes in
$\mathcal M$ that are below $H_0$ is at most
\begin{equation}\label{eq:binnedcount}
A_\sigma=\sum_{B<J}P_B
 \exp\!\left(\sigma\min\{\log H_0,h(B+K)\}\right).
\end{equation}
The same complete-set upper bound holds after arbitrary prime deletion.
For $\sigma=0$, $e_{b,j}=\binom{|\mathcal M_b|}{j}$ and $A_0=\sum_{B<J}P_B$.
\end{lemma}
\begin{proof}
A true product $d<H_0$ contains at most $K$ primes. Its label $B$ satisfies
$hB\le\log d<\log H_0$, hence $B<J$. Its block degrees are at most $K$,
so its positive contribution $d^{-\sigma}$ occurs in $P_B$.
Every intermediate label is at most its final label; truncation modulo
$v^J$ after each block cannot remove that contribution. Also
$\log d\le h(B+K)$, so its contribution to \eqref{eq:binnedcount} is at
least one. Products of excessive total degree that survive the block
truncations are harmless additional positive contributions.
Deleting primes decreases every elementary coefficient $e_{b,j}$,
while the complete-set degree and label bounds remain valid.
\end{proof}

\begin{lemma}[Counting rough integers]\label{lem:integerrough}
Let $V_T=\prod_{p\le T}(1-1/p)$, and let $\eta_+$ be the complete
normalized medium upper-gap budget supplied by Lemma \ref{lem:prefixblocks}
for $g(p)=1/(p-1)$. With the same $D_0$, $H_0$ and exact odd product $Q_*$,
the number $R_T(D)$ of positive integers $n<D$ having no prime factor
at most $T$ satisfies
\begin{equation}\label{eq:integerrough}
R_T(D)\le a_TD+b_T,\qquad
a_T=(1+\eta_+)V_T,\qquad b_T=2Q_*H_0.
\end{equation}
\end{lemma}
\begin{proof}
Use full inclusion--exclusion at all primes at most $w_0$, including $2$,
and the finite Rosser upper coefficients at $w_0<p\le T$.
For the integer sequence the density is $g_0(p)=1/p$.
The normalized positive recurrences in Lemma \ref{lem:prefixblocks}
also hold with $t_p=g(p)/(1-g(p))$: their passing and failing operators
depend only on the prime labels and $D_0$, and every coefficient of
$D(a)$ and $C(a)$ is a nonnegative polynomial in the $t_p$.
Thus decreasing $t_p$ decreases every propagated state and gap bound.
Here $t_{0,p}=1/(p-1)\le1/(p-2)$, so the normalized medium upper gap
at density $g_0$ is at most $\eta_+$. Multiplying by the exact small-prime
factor gives an upper main mass $\sum_d\lambda_d^+/d\le a_T$.

These integer upper weights have $|\lambda_d^+|\le1$ and support
$d<2Q_*H_0$ by Lemma \ref{lem:sharpsupport}; the factor $2$ is required
for the exact even prime. The number of positive multiples of $d$
strictly below $D$ is $\lceil D/d\rceil-1$, differing from $D/d$ by
an amount in $[-1,0)$. Summing the pointwise upper weights therefore gives
\[
R_T(D)\le\sum_d\lambda_d^+(\lceil D/d\rceil-1)
\le D\sum_d\frac{\lambda_d^+}{d}+\#\{d:\lambda_d^+\ne0\}
\le a_TD+b_T.
\]
This is an elementary integer-count bound, not an AP distribution input.
\end{proof}

Put $m=\omega(Q_*)$ for the complete exact odd-prime set and
\[
\nu_\sigma=\frac{2^m A_\sigma}{Q_*H_0}.
\]
Every large-prime factor of a composite coefficient has no prime factor
at most $T$. The coefficient support therefore has size at most
$2^mA_\sigma(a_TD+b_T)$. If $Q_*H_0D\le H$, its size is at most
\begin{equation}\label{eq:roughsupport}
\nu_\sigma(a_T+b_T/D)H.
\end{equation}
This counts all possible exact, medium and large factors; any vanishing
coefficient only reduces the support. Prime deletion preserves the bound.
For the aggregate $A_q$ remainder, write $d=e d_0d_1$ with large level
$D/q$ and absorb $q$ into $d_1$. The resulting large factor
$qd_1<D$ has no prime factor at most $T$. Because $d\mid P(z)$ and
$q\ge z$, $q$ is the unique prime factor of $qd$ that is at least $z$.
The map $(q,d)\mapsto qd$ is therefore injective over the whole sum,
and its complete image has at most $2^mA_\sigma R_T(D)$ elements.
Thus \eqref{eq:roughsupport} bounds the aggregate first AP remainder,
without summing a separate floor charge for each $q$.
The separate prime-modulus centering cost is retained.

For $T=5000000$, $w_0=5$, $u=87/16$ and $\kappa=11/169$,
use the preceding gap budgets and $\sigma=0$, $h=1/100$ for the count.
The complete count uses 348510 medium primes, $K=22$, $J=8115$,
and 992 nonempty blocks. Its positive polynomial calculation yields
\[
A_0<5.488142\cdot10^{30},\qquad
\nu_0<0.000008437984<\nu_*:=1/118511.
\]
The gap computation still uses $h=1/200$ as in Section \ref{sec:prefixgrh};
the two bin widths serve different bounds.
An independently regenerated complete integer prime set contains 348513
primes through $T$ and gives
\[
V_T<0.036397999733,\qquad
a_T<a_*:=36411450583/10^{12}.
\]

\paragraph{Finite analytic interval and the existing tail.}
Let $L_*=__ANCHOR__$, fix $\alpha=__ALPHA__$, $\beta=__BETA__$,
and put $D=e^{tL}$, $t=1/2-\alpha$, $H=e^{L/2}/L^\beta$.
Set $L_0=__START__$ and $c=\log(Q_*H_0)$.
On the closed interval $L_0\le L\le L_*$, the support condition
$\alpha L-\beta\log L>c$ persists since
$\alpha-\beta/L_0>0$. Also $H\ge10^9$, $H\le\sqrt N$,
$T<z$, $3<s<4$ and $k>1$ hold throughout.
Using \eqref{eq:roughsupport} in the pointwise GRH AP bound gives the
normalized leading expense
\[
\frac{3\nu_*}{2U_*}(a_*+b_Te^{-tL})P(L)L^{3-\beta}.
\]
Because $0<\beta<3$, this expression is not claimed to decrease on a
whole half-line. Instead, use the constant upper bound
\begin{equation}\label{eq:intervalplateau}
C_{\rm AP}=\frac{3\nu_*}{2U_*}
(a_*+b_Te^{-tL_0})P(L_0)L_*^{3-\beta}
\end{equation}
throughout this finite interval. This is valid because $P(L)$ and
$e^{-tL}$ decrease, and $L^{3-\beta}\le L_*^{3-\beta}$.
The positive floor charge is kept; at $L_0$,
$b_Te^{-tL_0}<10^{-5302}$.

Replace only the leading AP term by $C_{\rm AP}$.
Both other AP terms described after \eqref{eq:denseap}, every prime-modulus centering
charge, \eqref{eq:Erect}, and the finite loss are retained.
The bilinear estimate still uses the maximum modulus $H$.
Its proof in Section 6 needs $H\le\sqrt N$ and $H\ge10^9$, both verified
here; it does not use the smaller support count.
The endpoint certificate gives $s=__S__$ and
\begin{center}
\begin{tabular}{lr}\toprule
Term & Certified value\\\midrule
__TABLE__
$E_{\rm rect}$ & $<10^{-194}$\\
$E_{\rm fin}$ & $<10^{-1337}$\\\bottomrule
\end{tabular}
\end{center}
For these fixed parameters $M_0$ is constant, while $M_1$, $M_2$,
the other AP terms and all remaining expenses decrease for $L\ge L_0$.
Thus the strict surplus exceeding $10^{-5}$ at $L_0$ holds throughout
$[L_0,L_*]$. Section \ref{sec:prefixgrh} already proves the same
claim for every $L\ge L_*$, with common normalization $U_NN/L^2$.
These inclusive ranges cover the whole half-line $L\ge L_0$ and prove
the preceding starting point $e^{24799}$ using GRH alone.
The preceding stronger $0.0179$ coefficient at $L\ge31000$ retains
its earlier half-line certificate.

The separate-formula check reproduces 12 endpoint quantities and
regenerates the integer product with an independent prime sieve.
The histogram checks enumerate 9216 subsets and verify 17172 exact
rational coefficients at $\sigma=0,1$; the integer checks cover 16384
divisor-floor terms and 24316 roughness patterns.
An additional finite enumeration checks 8608 prime-factor modulus patterns,
including the injective aggregate map and its rough-factor support count.
These finite checks accompany the two general proofs above.
The 135 endpoint candidates are a finite parameter grid, not a global
optimality proof. This analytic-regime join does not connect to the
finite Goldbach verification endpoint $4\cdot10^{18}$.
% END CURRENT ROUGH INTERVAL
'''
    values={'START':start,'ANCHOR':anchor,'ALPHA':str(Fraction(row['alpha'])),
        'BETA':row['beta'],'S':f"{float(arb(row['s_ball'])):.6f}",'TABLE':'\n'.join(table)}
    for k,v in values.items():text=text.replace('__'+k+'__',str(v))
    assert '__' not in text
    return text


def main():
    report=load_checked('conditional_chen_rough_interval_certificate.json')
    checks=load_checked('conditional_chen_rough_interval_checks.json')
    proof=section(report,checks)
    path=ROOT/'conditional_chen_grh_bridge.tex';tex=path.read_text(encoding='utf-8')
    front,rest=tex.split('Let $B=4',1)
    start=report['GRH_positivity']['log_N0']
    for old in [26820,24799]:front=front.replace(str(old),str(start))
    front=front.replace('10.196903156',bound(arb(start).log(),9))
    digits=int((arb(start)/arb(10).log()).floor().unique_fmpz())+1
    front=front.replace('$11648$',f'${digits}$')
    if 'Logarithmic-bin support counts' not in front:
        front=front.replace('Rankin product bounds the number of actual moduli.',
            'Rankin product bounds the number of actual moduli.\n'
            'Logarithmic-bin support counts and an elementary rough-integer bound\n'
            'sharpen that count; a finite analytic interval is joined to an existing GRH tail.')
    front=front.replace('Positive prefix blocks retain all preceding checks, and a finite\n'
        'Rankin product bounds the number of actual moduli.\n'
        'Logarithmic-bin support counts and an elementary rough-integer bound\n'
        'sharpen that count; a finite analytic interval is joined to an existing GRH tail.',
        'Positive prefix blocks retain all preceding checks. Logarithmic counts\n'
        'and an elementary rough-integer bound reduce the AP expense; a finite\n'
        'analytic interval joins an existing GRH tail.')
    tex=front+'Let $B=4'+rest
    tex=tex.replace('The following calculation proves the current Theorem \\ref{thm:grh}.',
        'The following calculation proves the earlier starting point $e^{26820}$\n'
        'and the stronger coefficient $0.0179$ in Theorem \\ref{thm:grh}.')
    if '% BEGIN CURRENT ROUGH INTERVAL' in tex:
        a,b=tex.split('% BEGIN CURRENT ROUGH INTERVAL',1);_,b=b.split('% END CURRENT ROUGH INTERVAL',1)
        tex=a+proof+b.lstrip('\n')
    else:
        tex=tex.replace('\\section{A parameterized EH-style reduction}',proof+'\n\\section{A parameterized EH-style reduction}',1)
    ledger=r'''% BEGIN ROUGH INTERVAL LEDGER
For the logarithmic count and rough-integer interval, run:
\begin{quote}\small
\path{conditional_chen_binned_support_checks.py}\\
\path{conditional_chen_rough_interval_certificate.py}\\
\path{conditional_chen_rough_interval_checks.py}
\end{quote}
The support histograms are generated completely or reused from per-case
caches after checking generator hashes. The interval program computes the
complete integer Euler product freshly; the separate checker regenerates
its prime set independently. It rechecks the old tail, verifies the common
claim and inclusive join, and uses its own formulas for 12 new endpoint
quantities. Continuous coverage follows from the fixed-parameter argument
in Section \ref{sec:roughinterval}, not from the five sample points.
This work does not supply independent research review or a splice to the
finite verification range.
% END ROUGH INTERVAL LEDGER
'''
    if '% BEGIN ROUGH INTERVAL LEDGER' in tex:
        a,b=tex.split('% BEGIN ROUGH INTERVAL LEDGER',1);_,b=b.split('% END ROUGH INTERVAL LEDGER',1)
        tex=a+ledger+b.lstrip('\n')
    else:tex=tex.replace('\\section*{Reproducibility and dependency ledger}',
            '\\section*{Reproducibility and dependency ledger}\n'+ledger,1)
    assert tex.count('\\label{lem:binnedcount}')==1 and tex.count('\\label{lem:integerrough}')==1
    path.write_text(tex,encoding='utf-8',newline='\n')
    print(json.dumps({'GRH_log_N0':start,'analytic_tail_start':report['analytic_tail_join_log_N']}))


if __name__=='__main__':main()
