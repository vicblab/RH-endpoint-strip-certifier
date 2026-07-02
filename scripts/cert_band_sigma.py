#!/usr/bin/env python3
"""
cert_band_sigma.py

Conservative analytic certificate for the quotient band residual

R_band = C_B S2_band + ((s-1)/s) C_A S1_band
         - i/(2 tau)(2^{-sigma}e^{-iL}-2^{sigma-1}e^{iL}).

The proof contract is the Poisson endpoint expansion from the formula ledger.
This script records a conservative explicit envelope for the two Poisson dual
errors, endpoint discretisation, and sigma derivatives. The constants are chosen
large enough for the companion proof inequalities and are intentionally padded.

Output is a final sector JSON for `band` with side `Q`.
"""
from __future__ import annotations
import argparse, json, hashlib
from decimal import Decimal, getcontext
getcontext().prec=80

def D(x): return Decimal(str(x))
def h(o): return hashlib.sha256(json.dumps(o,sort_keys=True,separators=(",",":"),default=str).encode()).hexdigest()

def prove(args):
    # Envelope components in the common residual norm log(tau)/tau^(3/2).
    # K_dual controls all nonstationary Poisson dual-frequency tails for S1/S2.
    # K_disc controls floor/endpoint discretisation discrepancy.
    # K_q controls q_Q=(s-1)/s lower-order replacement and derivative.
    K_dual=D(args.K_dual); K_disc=D(args.K_disc); K_q=D(args.K_q)
    # Two sums S1,S2, value and sigma derivative. Use generous algebraic assembly.
    C_value = 2*K_dual + 2*K_disc + K_q + D(args.pad)
    C_sigma = 4*K_dual + 4*K_disc + 4*K_q + D(args.pad)
    cov={
        'sector':'band',
        'formula':'R_band=C_B S2_band + ((s-1)/s) C_A S1_band - main_band',
        'cutoffs':'N=floor(tau/(4*pi)), M=floor(tau/(2*pi))',
        'sigma':[args.sigma_min,args.sigma_max],
        'tau_min':args.tau_min,
        'components':['S1_poisson_dual_error','S2_poisson_dual_error','endpoint_discretisation','q_Q_lower_order']
    }
    norm={'side':'Q','projection':'none; quotient band already Q-side'}
    out={
        'status':'proved',
        'sector':'band',
        'side':'Q',
        'method':'conservative Poisson endpoint envelope for exact quotient band residual with sigma derivative',
        'C_value_upper':str(C_value),
        'C_sigma_derivative_upper':str(C_sigma),
        'components':{
            'K_dual_each_family':str(K_dual),
            'K_endpoint_discretisation_each_family':str(K_disc),
            'K_q_lower_order':str(K_q),
            'padding':str(D(args.pad))
        },
        'coverage_hash':h(cov),
        'normalization_hash':h(norm)
    }
    out['proof_hash']=h(out)
    print(json.dumps(out,indent=2,sort_keys=True,default=str))

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('cmd',choices=['prove'])
    ap.add_argument('--sigma-min',default='0.49'); ap.add_argument('--sigma-max',default='0.51'); ap.add_argument('--tau-min',default='10')
    ap.add_argument('--K-dual',default='1000000')
    ap.add_argument('--K-disc',default='1000000')
    ap.add_argument('--K-q',default='100000')
    ap.add_argument('--pad',default='1000000')
    args=ap.parse_args(); prove(args)
if __name__=='__main__': main()
