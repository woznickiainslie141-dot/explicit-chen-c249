"""Insert the current proved prefix-block and support-count sections in TeX.

This changes presentation only after hash-checked strict endpoint results
and the separate scalar checks are available. It does not assert a splice.
"""
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path
from flint import arb,ctx

ROOT=Path(__file__).resolve().parent
ctx.prec=192


def frac(value):
    x=Fraction(value)
    return str(x.numerator) if x.denominator==1 else f'{x.numerator}/{x.denominator}'


def bound(value,places=10,upper=False):
    x=arb(value); scale=10**places
    n=math.ceil(float(x)*scale) if upper else math.floor(float(x)*scale)
    while not (x<arb(n)/scale if upper else x>arb(n)/scale):n+=1 if upper else -1
    return f'{n//scale}.{n%scale:0{places}d}'


def load_checked(name):
    result=json.loads((ROOT/name).read_text(encoding='utf-8'))
    assert result['status'].startswith('PASS')
    for p,v in result['input_sha256'].items():
        assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==v,p
    return result


def sections():
    sparse=load_checked('conditional_chen_sparse_modulus_certificate.json')
    refined=load_checked('conditional_chen_prefix_block_refinement.json')
    checks=load_checked('conditional_chen_prefix_block_endpoint_checks.json')
    low,strong=sparse['GRH_positivity'],sparse['GRH_at_31000']
    assert checks['GRH_checks'][0]['log_N0']==low['log_N0']
    pre=low['presieve']; count=low['sparse_modulus_count']; chosen=count['selected']
    table=[]
    for key,formatter in [
        ('log_N0',lambda x:str(x)),('alpha',frac),('beta',frac),
        ('s_ball',lambda x:f'{float(arb(x)):.6f}'),('support_slack_ball',lambda x:'>'+bound(x)),
        ('AP_budget_ball',lambda x:'<'+bound(x,10,True)),
        ('margin_ball',lambda x:'>'+bound(x,12))]:
        label={'log_N0':'$L_0$','alpha':r'$\alpha$','beta':r'$\beta$',
               's_ball':'$s$','support_slack_ball':'Support slack',
               'AP_budget_ball':r'$E_{\rm AP}$','margin_ball':'Surplus'}[key]
        table.append(label+' & $'+formatter(low[key])+'$ & $'+formatter(strong[key])+'$'+r'\\')
    table='\n'.join(table)
    grh=r'''% BEGIN CURRENT PREFIX GRH
\subsection{All preceding checks and the actual modulus count}\label{sec:prefixgrh}
The following calculation proves the earlier starting point $e^{26820}$
and the preceding stronger coefficient $0.0179$.
Use Lemma \ref{lem:prefixblocks} with $T=__T__$, $w_0=5$,
$Q_*=15$, $u=__U__$, and $h=__H__$.
The complete integer sieve has __PRIMES__ medium primes;
the exact degree bound is $K=__K__$, the label limit is $J=__J__$,
and there are __BLOCKS__ nonempty blocks. The strict positive recurrences give
\[
\widehat R_+<__UP__<__EP__,\qquad
\widehat R_-<__LO__<__EM__.
\]
Set $\epsilon=__EPS__$, $\eta_+=__EP__$, $\eta_-=__EM__$,
$(C_1,C_2)=(106,108)$, and retain $\kappa=11/169$.
The actual support cost is
\[
c=\log(\kappa Q_*)+u\log T=__COST__\ldots.
\]
All uses of \eqref{eq:sieveparity} accept the new separate gap inputs;
its proof depends only on their domination of the two finite masses.

\begin{lemma}[Counting the actual modulus support]\label{lem:sparsecount}
For the composite coefficients in Lemma \ref{lem:composite}, put
$H_0=\kappa D_0$ and let $m$ be the number of complete exact primes.
For every $\sigma>0$, the number of their nonzero coefficients is at most
\[
2^m D H_0^\sigma
\prod_{w_0<p\le T}(1+p^{-\sigma}).
\]
This bound holds for every deletion of primes dividing $N$ and both signs.
If $Q_*H_0D\le H$, this number is at most $\nu H$, where
\begin{equation}\label{eq:sparsefraction}
\nu=\frac{2^m}{Q_*}H_0^{\sigma-1}
\prod_{w_0<p\le T}(1+p^{-\sigma}).
\end{equation}
\end{lemma}
\begin{proof}
Every contributing divisor has a unique exact, medium and large factorization.
There are at most $2^m$ exact divisors. Every medium factor is square-free,
uses only primes in $(w_0,T]$, and is below $H_0$ by Lemma
\ref{lem:sharpsupport}, for either medium sign. Rankin's elementary inequality
$1_{d_0<H_0}\le(H_0/d_0)^\sigma$ bounds their number by the Euler product
in the statement. A large factor is a positive integer below $D$,
so there are fewer than $D$ choices. Counting all such triples bounds
the coefficient support, even when some coefficients vanish.
Deleting primes decreases the positive Euler product and the available
exact factors, while the complete-set $H_0$ still bounds their support.
Dividing by $Q_*H_0D\le H$ proves the second assertion.
\end{proof}

For $\sigma=__SIGMA__$, a fresh 192-bit Euler-product sum over the complete
medium primes gives
\[
\log\prod_{5<p\le T}(1+p^{-\sigma})<__LOGPROD__,\qquad
\nu<__NU__<__NUBUDGET__.
\]
For the prime-factor sieves $A_q$, $q\ge z>T$ and the large level is $D/q$.
For their aggregate AP remainder, write $d=e d_0d_1$ and absorb $q$
into the large factor: $qd_1<D$. Since $d\mid P(z)$, the prime $q$
is the unique prime factor of $qd$ that is at least $z$. Thus
$(q,d)\mapsto qd$ is injective over the whole sum. All its images
belong to the same enlarged set of exact, medium and large triples,
so the same count bounds the aggregate first AP remainder.
The additional centering term
$r(q)/\varphi(d)$ retains its previous charge.

Set $H=\sqrt N/L^\beta$ with $\beta\ge3$. The pointwise AP bound and
centering argument used for \eqref{eq:R} give, for each actual support,
\[
R_{\rm support}\le\nu P(L)\frac{N}{L^{\beta-1}},\qquad
P(L)=q_+(L)+\frac1{16\pi}+\frac{2e^{-L/2}}{\log2}.
\]
Thus the lower sieve and half of the prime-factor sum charge a total of
$3R_{\rm support}/2$. Replace \eqref{eq:denseap} by
\begin{equation}\label{eq:sparseap}
\frac{3\nu P(L)}{2U_*L^{\beta-3}},
\end{equation}
using the fixed upper budget $\nu_*=__NUBUDGET__$ in this expense.
Both other AP terms, the prime-modulus centering cost, the full rectangle
charge, and every finite loss remain in the proof. The bilinear lemma
still uses the full maximum modulus $H$, rather than this count.

The strict endpoint certificate gives
\begin{center}
\begin{tabular}{lrr}\toprule
 & Lower starting point & Stronger count\\\midrule
__TABLE__
$E_{\rm rect}$ & $<10^{-212}$ & $<10^{-250}$\\
$E_{\rm fin}$ & $<10^{-1446}$ & $<10^{-1673}$\\\bottomrule
\end{tabular}
\end{center}
The respective surpluses exceed $10^{-5}$ and $__STRONG__$.
All analytic ranges used in the earlier calculation are checked again,
including $H\ge10^9$, $T<z$, $3<s<4$, $k>1$, and actual support
$\alpha L-\beta\log L>c$.
For each fixed parameter pair, support slack and $H$ increase thereafter.
With $\beta=3$, the leading AP expense decreases with $P(L)$;
every other expense decreases as before. The inequalities therefore
hold on the whole half-lines. This uses GRH alone.

The previous prefix-block budget using the full square-free modulus count
gave $N\ge e^{28279}$ and coefficient $0.0116$ at $e^{31000}$.
The original finite grid gave $e^{28334}$.
The new support count accounts for the further improvement.
The searched finite parameters and the four Rankin exponents do not
prove global optimality.
% END CURRENT PREFIX GRH
'''
    replacements={'T':pre['T'],'U':pre['u'],'H':pre['log_bin_width'],
        'PRIMES':pre['medium_prime_count'],'K':pre['maximum_actual_prefix_degree'],
        'J':pre['bin_limit'],'BLOCKS':pre['nonempty_prime_blocks'],
        'UP':bound(pre['upper_relative_Rosser_gap_ball'],12,True),
        'LO':bound(pre['lower_relative_Rosser_gap_ball'],12,True),
        'EP':low['eta_plus'],'EM':low['eta_minus'],'EPS':low['epsilon'],
        'COST':bound(low['effective_log_support_cost_ball'],9),
        'SIGMA':chosen['sigma'],'LOGPROD':bound(chosen['log_euler_product_ball'],10,True),
        'NU':bound(chosen['count_fraction_ball'],12,True),
        'NUBUDGET':low['sparse_modulus_count_budget'],'TABLE':table,
        'STRONG':str(float(Fraction(strong['claimed_margin'])))}
    for key,value in replacements.items():grh=grh.replace('__'+key+'__',str(value))
    families=[]; identifiers={}; cells=[]
    for row in refined['distribution_endpoints']:
        case=row['case']; p=case['presieve']
        key=(p['T'],p['exact_cutoff'],p['u'],p['log_bin_width'])
        if key not in identifiers:
            identifiers[key]=len(families)+1; families.append(case)
        cells.append({**row,'family':identifiers[key]})
    family_table=[]
    for i,case in enumerate(families,1):
        p=case['presieve']
        family_table.append(f"${i}$ & ${p['T']}$ & ${p['u']}$ & ${p['log_bin_width']}$ & "
            f"${case['epsilon']}$ & ${case['eta_plus']}$ & ${case['eta_minus']}$ & "
            f"$({case['C1']},{case['C2']})$"+r'\\')
    endpoint_table=[]
    for theta in ['3/4','9/10','19/20','99/100']:
        rows=[r for r in cells if r['theta']==theta]
        endpoint_table.append('$'+theta+'$ & '+ ' & '.join('$'+str(r['log_start'])
            +r'\ ('+str(r['family'])+')$' for r in rows)+r'\\')
    distribution=r'''% BEGIN CURRENT PREFIX DISTRIBUTION
\subsection{Prefix-block finite-constant conversion}\label{sec:prefixdistribution}
Recompute the finite medium masses using Lemma \ref{lem:prefixblocks}.
The following five families use $w_0=3$, $Q_*=3$, $\kappa=5/49$,
$D_0=T^u$, and the indicated gap and product budgets:
\begin{center}
\small
\begin{tabular}{rrrrrrrr}\toprule
Family & $T$ & $u$ & $h$ & $\epsilon$ & $\eta_+$ & $\eta_-$ & $(C_1,C_2)$\\\midrule
__FAMILIES__
\bottomrule\end{tabular}
\end{center}
For a chosen cell, use that family's composite lower weights at
$z=N^{1/3}$ and $D=N^\theta/(\kappa Q_*D_0)$.
If only the signed input \eqref{eq:signedinput} is assumed, it must hold
for this specified family on its sieve domain. The full absolute input
\eqref{eq:EHinput} implies it for every family with the same $C,X$.
No value of $C$ or $X$ is supplied by this construction.

\begin{proposition}[Prefix-block budget conversion]\label{prop:prefixdistribution}
Under that cell's signed input, or the full absolute input, every even
$N\ge\max\{X,e^A\}$ satisfies $\pi_2(N)>10^{-5}N/\log^2N$.
Here each entry is $A$ followed by its family number, and $C$ is bounded
by the column label:
\[
\begin{array}{c|rrrr}
\theta&C\le1&C\le10^3&C\le10^4&C\le10^6\\\hline
__ENDPOINTS__
\end{array}
\]
\end{proposition}
\begin{proof}
The complete-set bounds of Lemma \ref{lem:prefixblocks} supply the two
gaps in \eqref{eq:sieveparity}; all coefficient, support and signed
remainder identities continue to apply. Use exactly
\eqref{eq:signedcriterion} with each family's $c=\log(\kappa Q_*)+u\log T$.
The 192-bit certificates check every cell's strict surplus, $T<z$,
$L>3\log(10^{12})$, and $2<s=3\theta-3c/L<3$.
For fixed parameters, $l(s)$ increases and $u(s)$ decreases on $(2,3)$,
$s$ and $\rho_U(L)$ increase, and every deducted expense decreases.
Thus every cell extends to the claimed half-line.
The exact and separate-formula checks test finite arithmetic, not a
distribution upper bound or independent mathematical review.
\end{proof}
Ordinary EH still supplies only unspecified constants and starting points.
These entries improve the conversion of a supplied quantitative input,
and are not numeric EH-only records. In particular $X$ cannot be discarded.
The smallest support cost among these five families is __MINCOST__;
even $\theta<1$ allows only $L_B/3<14.278$ at $B$.
Thus these families still cannot reach the verified endpoint, even with
zero signed error. The complete absolute input also retains the integrality
constraints of Proposition \ref{prop:tradeoff}.
% END CURRENT PREFIX DISTRIBUTION
'''
    assert len(families)==5 and all(c['presieve']['exact_cutoff']==3 for c in families)
    distribution=distribution.replace('__FAMILIES__','\n'.join(family_table))
    distribution=distribution.replace('__ENDPOINTS__','\n'.join(endpoint_table))
    distribution=distribution.replace('__MINCOST__',bound(min(
        (arb(c['effective_log_support_cost_ball']) for c in families),key=float),9))
    return sparse,refined,grh,distribution


def main():
    sparse,refined,grh,distribution=sections()
    low,strong=sparse['GRH_positivity'],sparse['GRH_at_31000']
    path=ROOT/'conditional_chen_grh_bridge.tex'; tex=path.read_text(encoding='utf-8')
    current=low['log_N0']; coefficient=str(float(Fraction(strong['claimed_margin'])))
    if (ROOT/'conditional_chen_rough_interval_certificate.json').exists():
        current=load_checked('conditional_chen_rough_interval_certificate.json')['GRH_positivity']['log_N0']
    if (ROOT/'conditional_chen_accepted_interval_certificate.json').exists():
        accepted=load_checked('conditional_chen_accepted_interval_certificate.json')
        current=accepted['GRH_positivity']['log_N0']
        coefficient=str(float(Fraction(accepted['GRH_at_31000']['claimed_margin'])))
    # Only the statement/abstract preamble is replaced; earlier proof tables stay historical.
    front,rest=tex.split('Let $B=4',1)
    for old in [30747,26820]:front=front.replace(str(old),str(current))
    front=front.replace('surplus $0.001$','surplus $'+coefficient+'$')
    front=front.replace('replaced by $0.001$','replaced by $'+coefficient+'$')
    front=front.replace('10.333547707',bound(arb(current).log(),9))
    digits=int((arb(current)/arb(10).log()).floor().unique_fmpz())+1
    front=front.replace('$13354$',f'${digits}$')
    front=front.replace('The actual Rosser support is also smaller',
        'Positive prefix blocks retain all preceding checks, and a finite\n'
        'Rankin product bounds the number of actual moduli.\nThe actual Rosser support is also smaller') if 'Positive prefix blocks retain' not in front else front
    tex=front+'Let $B=4'+rest
    tex=tex.replace('For the two assertions of Theorem \\ref{thm:grh}, use the following',
        'That preceding threshold had $\\log(30747)=10.333547707\\ldots$\n'
        'and about $13354$ decimal digits.\n'
        'For the preceding assertions $N\\ge e^{30747}$ with coefficient $10^{-5}$\n'
        'and $N\\ge e^{31000}$ with coefficient $0.001$, use the following')
    tex=tex.replace('This proves both assertions for all $N$ in their respective half-lines.',
        'This proves those earlier assertions for all $N$ in their respective half-lines.')
    if '% BEGIN CURRENT PREFIX GRH' in tex:
        a,b=tex.split('% BEGIN CURRENT PREFIX GRH',1); _,b=b.split('% END CURRENT PREFIX GRH',1)
        tex=a+grh+b.lstrip('\n')
    else:tex=tex.replace('\\section{A parameterized EH-style reduction}',grh+'\n\\section{A parameterized EH-style reduction}',1)
    if '% BEGIN CURRENT PREFIX DISTRIBUTION' in tex:
        a,b=tex.split('% BEGIN CURRENT PREFIX DISTRIBUTION',1); _,b=b.split('% END CURRENT PREFIX DISTRIBUTION',1)
        tex=a+distribution+b.lstrip('\n')
    else:tex=tex.replace('\\begin{proposition}[A support obstruction to triple factorability]',
        distribution+'\n\\begin{proposition}[A support obstruction to triple factorability]',1)
    assert tex.count('\\label{prop:prefixdistribution}')==1
    tex=tex.replace('the current costs are $47.560368848',
                    'the earlier fixed-family costs are $47.560368848')
    tex=tex.replace('For the current GRH endpoints, then run',
                    'For the preceding dense-grid GRH endpoints, then run')
    ledger=r'''For the current prefix-block results, run the following files in order:
\begin{quote}\small
\path{conditional_chen_prefix_block_checks.py}\\
\path{conditional_chen_prefix_block_endpoint_certificate.py}\\
\path{conditional_chen_prefix_block_refinement.py}
\end{quote}
The refinement records
per-case caches tied to generator hashes. The initial grid includes a
previously generated $u=11/2$ finite sum, explicitly reused by hash.
The separate support-count program
\path{conditional_chen_sparse_modulus_certificate.py} freshly evaluates
the complete Euler products at four rational exponents for each of three
cutoffs; it also checks 18 exact small Rankin counts.
Then run \path{conditional_chen_prefix_block_endpoint_checks.py}.
It independently reproduces eight GRH quantities at each asserted endpoint
and the main mass and surplus in all 16 distribution cells, without calling
the production scalar evaluators. The exact prefix-block checks cover
4608 deleted prime subsets, 9216 first-failure identities, 104976 divisor
patterns, and 132 block-coefficient groups. These tests do not replace
Lemmas \ref{lem:prefixblocks} and \ref{lem:sparsecount}, the external inputs,
independent research review, or a verification-range splice.
'''
    if 'For the current prefix-block results, run' not in tex:
        tex=tex.replace('\\section*{Reproducibility and dependency ledger}',
                        '\\section*{Reproducibility and dependency ledger}\n'+ledger,1)
    assert 'For the current prefix-block results, run' in tex
    path.write_text(tex,encoding='utf-8',newline='\n')
    print(json.dumps({'updated_GRH_log_N0':current,'stronger_coefficient':coefficient,
                      'distribution_cells':len(refined['distribution_endpoints'])}))


if __name__=='__main__':main()
