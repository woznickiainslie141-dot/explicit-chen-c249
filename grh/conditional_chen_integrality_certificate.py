"""Reject two proposed finite EH-style hypotheses by integer prime counts.

For P=pi(x), H<x, every prime q<=H with q-1 not dividing P-1
has max reduced-class error about P/(q-1) >=1/2. Exceptional q
are at most tau(P-1)<=2 sqrt(x). Thus the error sum is at least
pi(H)/2-sqrt(x). Use JS Lemma 3.1 at J=floor(H): pi(J)>J/log J.
This contradiction is unconditional, independent of GRH.
"""
from flint import arb,ctx
import json
from pathlib import Path
ctx.prec=192
B=arb(4)*arb(10)**18;L=B.log()
rows=[]
for theta,name,budget in [
    (arb(9)/10,'QEH_B level 9/10, constant 1/10, log power 2', B/(10*L**2)),
    (arb(3)/4,'old qEH level 3/4, constant 1, log power 4', B/L**4),
]:
    H=B**theta;J=H.floor()
    lower=J/(2*J.log())-B.sqrt()
    assert J>=71 and H<B and lower>budget
    rows.append({'candidate':name,'H_ball':str(H),'J_exact':str(J),
                 'integrality_lower_ball':str(lower),
                 'proposed_upper_ball':str(budget),'ratio_ball':str(lower/budget),
                 'verdict':'REFUTED at B; no nonvacuous bridge follows'})
out={'status':'unconditional arithmetic obstruction using the cited prime-count lower bound',
     'B':'4000000000000000000','precision_bits':ctx.prec,'candidates':rows}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps(out,indent=2))
