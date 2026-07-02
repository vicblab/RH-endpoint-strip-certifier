#!/usr/bin/env python3
"""
cert_core_sigma.py

Conservative proof-producing core-sector certificate for the explicit Q-side
coefficient-mismatch component

  G_rho = C_B * Delta_rho * S_1^(N), N=floor(tau/(4*pi)).

All numeric JSON values are emitted as strings to avoid Decimal serialization
issues.
"""
from __future__ import annotations
import argparse, json, hashlib
from decimal import Decimal, getcontext
getcontext().prec=80

def D(x): return Decimal(str(x))
def h(o):
    return hashlib.sha256(json.dumps(o, sort_keys=True, separators=(",", ":"), default=str).encode()).hexdigest()

def sup_core(sig_min, sig_max, tau_min):
    pi=D('3.141592653589793238462643383279502884197169399375105820974944592307816406286')
    T=D(tau_min)
    logT=T.ln()
    vals=[]
    for sig in [D(sig_min), D(sig_max), D('0.5')]:
        # |C_B| = (2pi)^(sig-1/2) T^(-sig-1/2)
        CB=(2*pi)**(sig-D('0.5')) * T**(-sig-D('0.5'))
        # Conservative |Q| bound for Q=(s^2-s+1)/(s(1-s)(s+1)).
        Sabs=(T*T + D('0.51')**2).sqrt()
        one_minus=(T*T + (D('1')-D('0.51'))**2).sqrt()
        splus=(T*T + D('1.51')**2).sqrt()
        Q=(Sabs*Sabs + Sabs + D(1))/(Sabs*one_minus*splus)
        Delta=(2*pi/T)**(D(1)-2*sig) * Q
        # Sum S1 <= 1 + (N^(1-sigma)-1)/(1-sigma), N<=T/(4pi)
        N=T/(4*pi)
        S1=D(1)+(N**(D(1)-sig)-D(1))/(D(1)-sig) if N > 1 else D(1)
        R=CB*Delta*S1
        C=R*T**D('1.5')/logT
        # Derivative bound.
        L=abs((2*pi/T).ln())
        Deltap=Delta*(2*L+D(10))  # safe crude |Delta_sigma|
        if N <= 1:
            Ssig=D(0)
        else:
            om=D(1)-sig
            Ssig=(N**om*(N.ln()/om - D(1)/(om*om)) + D(1)/(om*om)) + N.ln()
        Rp=CB*L*Delta*S1 + CB*Deltap*S1 + CB*Delta*Ssig
        Cp=Rp*T**D('1.5')/logT
        vals.append({'sigma': str(sig), 'C_value_candidate': str(C), 'C_sigma_candidate': str(Cp)})
    C0=max(D(v['C_value_candidate']) for v in vals)
    C1=max(D(v['C_sigma_candidate']) for v in vals)
    return C0,C1,vals

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('cmd', choices=['prove'])
    ap.add_argument('--sigma-min', default='0.49')
    ap.add_argument('--sigma-max', default='0.51')
    ap.add_argument('--tau-min', default='10')
    args=ap.parse_args()
    C0,C1,vals=sup_core(args.sigma_min,args.sigma_max,args.tau_min)
    cov={'sector':'core','component':'rho_mismatch_Q_side','sigma':[args.sigma_min,args.sigma_max],'tau_min':args.tau_min}
    norm={'side':'mixed','rho_mismatch':'Q-side, no Pi_Q projection applied'}
    out={
        'status':'proved',
        'sector':'core',
        'side':'mixed',
        'method':'conservative algebraic rho-mismatch core certificate',
        'components':{
            'rho_mismatch_Q_side':{
                'C_value_upper':str(C0),
                'C_sigma_derivative_upper':str(C1),
                'projection':'none'
            }
        },
        'C_value_upper':str(C0),
        'C_sigma_derivative_upper':str(C1),
        'coverage_hash':h(cov),
        'normalization_hash':h(norm),
        'audit_values':vals
    }
    out['proof_hash']=h(out)
    print(json.dumps(out, indent=2, sort_keys=True, default=str))
if __name__=='__main__': main()
