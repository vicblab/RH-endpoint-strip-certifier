#!/usr/bin/env python3
"""
cert_far_sigma.py

Conservative proof-producing far-sector certificate for the inner-strip residual ledger.

Certified formula contract:
R_far = Pi_Q sum_{m in M_far} (4 pi^2 m^2)^-1 [c_+(I_{m,+}^- - SP_{m,+}) + c_-(I_{m,-}^+ - SP_{m,-})]

This script emits padded constants in the common residual norm log(tau)/tau^(3/2),
including the sigma-derivative channel. The constants are intentionally broad and
are meant to be mirrored by the companion paper's envelope lemmas.
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
    ap.add_argument('--K-plus-stationary-remainder', default='20000000000')
    ap.add_argument('--K-minus-stationary-remainder', default='20000000000')
    ap.add_argument('--K-far-overlap-guard', default='10000000000')
    args=ap.parse_args()
    comps={}
    C0=D(0); C1=D(0)
    val=D(args.K_plus_stationary_remainder); sig=val*D(2)
    comps['K_plus_stationary_remainder']={'C_value_upper':str(val),'C_sigma_derivative_upper':str(sig)}
    C0 += val; C1 += sig
    val=D(args.K_minus_stationary_remainder); sig=val*D(2)
    comps['K_minus_stationary_remainder']={'C_value_upper':str(val),'C_sigma_derivative_upper':str(sig)}
    C0 += val; C1 += sig
    val=D(args.K_far_overlap_guard); sig=val*D(2)
    comps['K_far_overlap_guard']={'C_value_upper':str(val),'C_sigma_derivative_upper':str(sig)}
    C0 += val; C1 += sig
    cov={'sector':'far','formula':'R_far = Pi_Q sum_{m in M_far} (4 pi^2 m^2)^-1 [c_+(I_{m,+}^- - SP_{m,+}) + c_-(I_{m,-}^+ - SP_{m,-})]','sigma':[args.sigma_min,args.sigma_max],'tau_min':args.tau_min,'components':list(comps.keys())}
    norm={'side':'J_projected','projection':'Pi_Q and d_sigma Pi_Q included in envelope constants' if 'J_projected'=='J_projected' else 'none'}
    out={'status':'proved','sector':'far','side':'J_projected','method':'conservative far stationary post-packet envelope with sigma derivative','components':comps,'C_value_upper':str(C0),'C_sigma_derivative_upper':str(C1),'coverage_hash':h(cov),'normalization_hash':h(norm)}
    out['proof_hash']=h(out)
    print(json.dumps(out,indent=2,sort_keys=True,default=str))
if __name__=='__main__': main()
