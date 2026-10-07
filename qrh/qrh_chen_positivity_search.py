"""Exhaustive, explicitly bounded search for positivity, not a chosen coefficient.

Only the freshly reproduced u=11/2 construction is used.  This checks
endpoints of a fixed recipe, not the optimality of Chen's method.
"""
import hashlib
import json
from fractions import Fraction
from pathlib import Path

from flint import arb
import qrh_chen_certificate as q

ROOT=Path(__file__).resolve().parent


def run():
    reproduction=json.loads((ROOT/'qrh_chen_reproduction.json').read_text(encoding='utf-8'))
    for name,sha in reproduction['input_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==sha, name
    case=q.selected_case()
    pre=reproduction['finite_presieve']
    assert (pre['T'],pre['exact_cutoff'],pre['u'])==(5000000,5,'11/2')
    assert arb(pre['upper_relative_Rosser_gap_ball'])<arb(case['eta_plus'])
    assert arb(pre['lower_relative_Rosser_gap_ball'])<arb(case['eta_minus'])
    rows=[]
    outcomes=[]
    checked=0
    for n in [4,5,6,7,8]:
        delta=Fraction(1,10**n)
        before=None
        for L in range(50000,51201):
            row=q.endpoint(case,L,delta)
            checked+=1
            assert row is not None
            value=arb(row['margin'])
            if value>0:
                assert before is not None and arb(before['margin'])<0
                row['last_nonpositive_endpoint']=before
                rows.append(row)
                outcomes.append({'delta':str(delta),'first_positive_log_N':L,
                                 'preceding_integer_log_N':L-1,'comparisons':L-50000+1})
                print(f'delta={delta}: first positive integer endpoint in range is {L}',flush=True)
                break
            assert value<0, 'Unresolved sign; increase precision'
            before=row
        else:
            outcomes.append({'delta':str(delta),'first_positive_log_N':None,
                             'comparisons':1201,'last_checked_margin':before['margin']})
    assert rows
    best=min(rows,key=lambda row:row['log_N0'])
    # Check the L-dependent starting gates for the smaller fixed-parameter
    # half-line.  Universal contour and rectangle constants were separately
    # regenerated in the full reproduction.  None of these checks is a proof
    # of the derivative/continuum arguments recorded in the audit.
    L=arb(best['log_N0'])
    w=L-5*L.log()
    x3=L/8-arb(2).log()
    r=arb(3)/64
    prime_error=lambda t:10000000*(-r*t).exp()*(t+arb(50).log())**2
    lower_mass_gate=L**2*prime_error(L)+L**3*(-L).exp()/arb(2).log()
    upper_mass_gate=L**2*prime_error(L)+3/L
    pi_gate=L*(L/3)*prime_error(L/3)+27/L
    small_conductor_gate=arb('1.1e6')*w**20*(-r*w).exp()
    assert lower_mass_gate<1 and upper_mass_gate<1 and pi_gate<1
    assert small_conductor_gate<arb('3.2e-8')
    assert w>50000 and w.log()>arb('10.4')
    assert w/2-20*w.log()>0
    a=(L-5)/w
    b=a*(w-5)/(w-5*w.log())
    assert 1<a<arb('1.01') and b>1
    assert 2-6+arb('1.01')/w.log()<0
    assert 2+arb('6.5')*arb('1.01')-10<0
    assert 2+arb('1.5')*a-(L-5)/12<0
    assert 2+arb('7.5')*a-(L-5)/6<0
    assert 6+(L-5)/2-L-10*a-1/L.log()<0
    assert 17/x3<r and 16/x3<1 and 5-L/16<0
    fixed=q.endpoint(case,best['log_N0'],Fraction(best['delta']),Fraction(best['alpha']))
    assert abs(arb(fixed['margin'])-arb(best['margin']))<arb('1e-45')
    best['fixed_parameter_half_line_starting_gates']={
        'A_lower_gate':str(lower_mass_gate),'A_upper_gate':str(upper_mass_gate),
        'pi_upper_gate':str(pi_gate),'all_low_conductors_gate':str(small_conductor_gate),
        'w_ball':str(w),'AP_log_derivative_ratio_ball':str(a),
        'AP_nested_log_derivative_ratio_ball':str(b),
        'scope':'strict starting gates; written continuum arguments still required',
    }
    result={
        'status':'PASS: exhaustive finite positivity-recipe endpoint scan',
        'scope':{'integer_log_N_range':[50000,51200],'presieve_u':'11/2',
                 'delta_values':['1/'+str(10**n) for n in [4,5,6,7,8]],
                 'alpha_recipe':'ceiling_1e-9((cost+12.5 log L)/L)',
                 'comparisons_performed':checked,'target':'strictly positive, not >1e-5'},
        'selected':best,'first_positive_rows':rows,'scan_outcomes':outcomes,
        'not_checked':['general analytic half-line proof',
                       'other presieve constructions','L below 50000',
                       'other alpha choices or delta values','global optimality',
                       'external seven-eighths theorem proof'],
        'input_sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                         for name in [Path(__file__).name,'qrh_chen_certificate.py',
                                      'qrh_chen_reproduction.json']},
    }
    (ROOT/'qrh_chen_positivity_search.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(best,indent=2))
    return result


if __name__=='__main__':run()
