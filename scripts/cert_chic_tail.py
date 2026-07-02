#!/usr/bin/env python3
"""
cert_chic_tail.py

STRICT large-r compact-complement tail certificate.

It certifies the hook C_chic_tail PROVIDED that a separate, proved derivative-envelope JSON
is supplied. This script does not accept dummy raw constants as proof.

Tail model certified here
-------------------------
For the compact-complement/high-r tail, suppose the scaled integrand H(r,y) on
0 <= y <= 2Y satisfies

  H(0)=H(2Y)=0,
  |dH/dy| <= M1,
  |d^2H/dy^2| <= M2

uniformly over sigma in [0.01,0.49], w in [0,1/sqrt(10)], lambda in [0,1],
and r >= R.

For phase theta_r(y)=y^2/8+2*pi*r*y, theta'=y/4+2*pi*r, theta''=1/4,
we have for r>=R

  |int_0^{2Y} H(y) exp(i theta_r(y)) dy| <= A2 / r^2

with

  A2 = (2*M1 + 2Y*M2)/(2*pi)^2
       + (4Y*M1*theta2)/(2*pi)^3/R
       + (4Y^2*M1*theta2^2)/(2*pi)^4/R^2,
  theta2 = 1/4.

This gives the quotient-side tail hook

  C_chic_tail = 8*Kstar*A2*(1/R + 1/(R^2*sqrt(T0))) / log(T0)

using Kstar=1.262751, T0=10.

Input envelope JSON must contain:
  {
    "status": "proved",
    "M1_hprime": "...",
    "M2_hsecond": "...",
    "method": "...",
    "proof_hash": "..."
  }

Run:
  python cert_chic_tail.py spec
  python cert_chic_tail.py prove --envelope-json proved_tail_derivative_envelope.json
"""
import argparse, json, hashlib
from decimal import Decimal, getcontext
from pathlib import Path
getcontext().prec = 80

PI = Decimal('3.141592653589793238462643383279502884197169399375105820974944592307816406286')
T0 = Decimal('10')
KSTAR = Decimal('1.262751')


def D(x): return Decimal(str(x))


def compute_C(M1, M2, Y=Decimal(128), R=Decimal(32)):
    two_pi = 2*PI
    theta2 = Decimal(1)/4
    A2 = (2*M1 + 2*Y*M2)/(two_pi**2)
    A2 += (4*Y*M1*theta2)/(two_pi**3 * R)
    A2 += (4*(Y**2)*M1*(theta2**2))/(two_pi**4 * (R**2))
    C = 8*KSTAR*A2*(Decimal(1)/R + Decimal(1)/(R**2 * T0.sqrt()))/T0.ln()
    return A2, C


def load_envelope(path):
    data = json.loads(Path(path).read_text())
    if data.get('status') != 'proved':
        raise SystemExit(f"Envelope {path} is not proved; status={data.get('status')}")
    for key in ['M1_hprime','M2_hsecond','method']:
        if key not in data:
            raise SystemExit(f"Envelope {path} missing {key}")
    if not (data.get('proof_hash') or data.get('hash')):
        raise SystemExit(f"Envelope {path} missing proof_hash/hash")
    return D(data['M1_hprime']), D(data['M2_hsecond']), data


def prove(args):
    M1, M2, env = load_envelope(args.envelope_json)
    Y = D(args.Y); R = D(args.R)
    A2, C = compute_C(M1, M2, Y, R)
    out = {
        'status': 'proved',
        'C_chic_tail': str(C),
        'A2_rminus2_bound': str(A2),
        'inputs': {
            'Y': str(Y),
            'R': str(R),
            'T0': str(T0),
            'Kstar': str(KSTAR),
            'M1_hprime': str(M1),
            'M2_hsecond': str(M2),
            'envelope_hash': env.get('proof_hash') or env.get('hash')
        },
        'method': 'large-r compact-complement tail by two integrations by parts in y using proved derivative envelope',
    }
    out['proof_hash'] = hashlib.sha256(json.dumps(out, sort_keys=True).encode()).hexdigest()
    print(json.dumps(out, indent=2))


def spec(_args):
    print(json.dumps({
        'status': 'spec',
        'required_input': {
            'status': 'proved',
            'M1_hprime': 'uniform upper bound for |dH/dy|',
            'M2_hsecond': 'uniform upper bound for |d^2H/dy^2|',
            'method': 'Arb/analytic proof of derivative envelope',
            'proof_hash': '...'
        },
        'output_key': 'C_chic_tail',
        'formula': 'C=8*Kstar*A2*(1/R+1/(R^2*sqrt(T0)))/log(T0)',
        'warning': 'This proves only the large-r/y tail hook assuming the derivative envelope is proved.'
    }, indent=2))


def main():
    p=argparse.ArgumentParser()
    sub=p.add_subparsers(dest='cmd', required=True)
    sub.add_parser('spec')
    q=sub.add_parser('prove')
    q.add_argument('--envelope-json', required=True)
    q.add_argument('--Y', default='128')
    q.add_argument('--R', default='32')
    args=p.parse_args()
    if args.cmd=='spec': spec(args)
    else: prove(args)

if __name__ == '__main__': main()
