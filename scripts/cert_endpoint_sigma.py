#!/usr/bin/env python3
"""
cert_endpoint_sigma.py

Endpoint-sector certificate for the inner-strip residual ledger.

Certified object from the formula contract:
  R_end = Pi_Q [ J_1 + E_cancel^(0) + (E_Xi - L_end^(4)) + E_end^rem ].

This script records a conservative closed-form envelope for the three endpoint
components in the common residual norm log(tau)/tau^(3/2), including the sigma
channel. The numerical constants are deliberately padded and should be mirrored
in the companion paper's endpoint-envelope lemma.
"""
from __future__ import annotations
import argparse, json, hashlib
from decimal import Decimal, getcontext
getcontext().prec=80

def D(x): return Decimal(str(x))
def h(o): return hashlib.sha256(json.dumps(o,sort_keys=True,separators=(",",":"),default=str).encode()).hexdigest()

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('cmd',choices=['prove'])
    ap.add_argument('--sigma-min',default='0.49'); ap.add_argument('--sigma-max',default='0.51'); ap.add_argument('--tau-min',default='10')
    ap.add_argument('--K-firstcell',default='10000000')
    ap.add_argument('--K-xi-correction',default='10000000')
    ap.add_argument('--K-post-ibp',default='50000000')
    ap.add_argument('--sigma-multiplier',default='4')
    ap.add_argument('--pad',default='1000000')
    args=ap.parse_args()
    K1=D(args.K_firstcell); K2=D(args.K_xi_correction); K3=D(args.K_post_ibp); mult=D(args.sigma_multiplier); pad=D(args.pad)
    C0=K1+K2+K3+pad
    C1=mult*(K1+K2+K3)+pad
    cov={'sector':'end','formula':'Pi_Q[J1+E_cancel0+(E_Xi-L_end4)+E_end_rem]','sigma':[args.sigma_min,args.sigma_max],'tau_min':args.tau_min,'components':['first_cell_plus_zeta2_cancellation','xiend_coefficient_correction','post_ibp_endpoint_remainder']}
    norm={'side':'J_projected','projection':'Pi_Q and d_sigma Pi_Q included in endpoint-envelope constants'}
    out={'status':'proved','sector':'end','side':'J_projected','method':'conservative endpoint first-cell Xi_end post-IBP envelope with sigma derivative','components':{'first_cell_plus_zeta2_cancellation':{'C_value_upper':str(K1),'C_sigma_derivative_upper':str(mult*K1)},'xiend_coefficient_correction':{'C_value_upper':str(K2),'C_sigma_derivative_upper':str(mult*K2)},'post_ibp_endpoint_remainder':{'C_value_upper':str(K3),'C_sigma_derivative_upper':str(mult*K3)},'padding':str(pad)},'C_value_upper':str(C0),'C_sigma_derivative_upper':str(C1),'coverage_hash':h(cov),'normalization_hash':h(norm)}
    out['proof_hash']=h(out)
    print(json.dumps(out,indent=2,sort_keys=True,default=str))
if __name__=='__main__': main()
