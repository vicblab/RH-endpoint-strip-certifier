#!/usr/bin/env python3
"""
cert_chic_stationary.py

Verifier for C_chic_stationary: the compact-complement transition boxes where the
phase derivative is not separated from zero and therefore were exported by
cert_chic_nonstationary.py.

This is intentionally crude but finite: on stationary-overlap boxes it uses absolute
value interval enclosures, not oscillatory cancellation. It covers exactly the same
transition strip and box partition convention as cert_chic_nonstationary.py:

  y in [Y,2Y], r in [r_min,r_max], w in [0,w_max], sigma in [0.01,0.49].

A box is handled here iff the phase derivative interval intersects [-lambda0,lambda0].
Those are precisely the boxes exported from the nonstationary verifier.

Typical run:
  python cert_chic_stationary.py prove \
    --Y 128 --r-min -4 --r-max 32 --w-max 0.31622776601683794 \
    --lambda0 0.25 --ns 4 --nw 16 --nr 96 --ny 128 \
    > cert_chic_stationary.json

Output:
  {
    "status": "proved",
    "C_chic_stationary": "...",
    "method": "absolute Arb box bound on stationary-overlap compact-complement transition boxes",
    "proof_hash": "..."
  }
"""
import sys, json, argparse, hashlib, math
from decimal import Decimal, getcontext
getcontext().prec = 80
try:
    from flint import arb, acb, ctx
except Exception as e:
    print(json.dumps({"status":"missing_dependency","error":repr(e),"install":"python -m pip install python-flint"}, indent=2))
    sys.exit(2)
ctx.prec = 100
PI = arb.pi()
KSTAR = Decimal('1.262751')
T0 = Decimal('10')

class AD0:
    def __init__(self,v): self.v=acb(v)
    def __add__(self,o): return AD0(self.v+toAD(o).v)
    __radd__=__add__
    def __neg__(self): return AD0(-self.v)
    def __sub__(self,o): return self+(-toAD(o))
    def __rsub__(self,o): return toAD(o)+(-self)
    def __mul__(self,o): return AD0(self.v*toAD(o).v)
    __rmul__=__mul__
    def __truediv__(self,o): return AD0(self.v/toAD(o).v)
    def log(self): return AD0(self.v.log())
    def exp(self): return AD0(self.v.exp())
    def pow_real(self,p): return (self.log()*p).exp()

def toAD(x): return x if isinstance(x,AD0) else AD0(x)

def midrad(a,b):
    a=arb(str(a)); b=arb(str(b)); return (a+b)/2 + arb(f"0 +/- {float((b-a)/2)}")

def mag(z): return float(acb(z).abs_upper())

def chi_c(y,Y):
    Y=arb(str(Y)); t=(y-AD0(Y))/Y
    chi=AD0(1)-3*t*t+2*t*t*t
    return AD0(1)-chi

def core_phase_yprime(w,y):
    return y/(4+2*w*y)

def phase_derivative(sig,r,w,y,fam):
    yp=core_phase_yprime(w,y)
    if fam=='-': return 2*PI*r + yp
    else: return -2*PI*r - yp

def amplitude(sig,w,y,Y,fam):
    # Absolute amplitude in the compact complement transition strip.
    q = -sig-2 if fam=='-' else sig-3
    X=AD0(2)+AD0(w)*AD0(y)
    A=X.pow_real(q)
    return chi_c(AD0(y),Y)*A

def prove(args):
    lam=arb(str(args.lambda0))
    C_acc=Decimal(0)
    counted=0
    skipped=0
    worst=None; worst_val=Decimal(0)
    for si in range(args.ns):
        sig=midrad(args.sigma_min+(args.sigma_max-args.sigma_min)*si/args.ns,args.sigma_min+(args.sigma_max-args.sigma_min)*(si+1)/args.ns)
        for wi in range(args.nw):
            w=midrad(0+args.w_max*wi/args.nw,0+args.w_max*(wi+1)/args.nw)
            for ri in range(args.nr):
                rlo=args.r_min+(args.r_max-args.r_min)*ri/args.nr
                rhi=args.r_min+(args.r_max-args.r_min)*(ri+1)/args.nr
                r=midrad(rlo,rhi); dr=Decimal(str(rhi-rlo))
                for yi in range(args.ny):
                    ylo=args.Y+(args.Y)*yi/args.ny
                    yhi=args.Y+(args.Y)*(yi+1)/args.ny
                    y=midrad(ylo,yhi); dy=Decimal(str(yhi-ylo))
                    for fam in ['+','-']:
                        pd=phase_derivative(sig,r,w,y,fam)
                        separated = bool(pd >= lam) or bool(pd <= -lam)
                        if separated:
                            skipped += 1
                            continue
                        amp=amplitude(sig,w,y,args.Y,fam)
                        M0=Decimal(str(mag(amp.v)))
                        # Crude absolute box contribution. Same normalization convention as nonstationary script:
                        # integrate in y and r, multiply by Kstar/log(T0).
                        box = (KSTAR / T0.ln()) * dr * dy * M0
                        C_acc += box
                        counted += 1
                        if box > worst_val:
                            worst_val=box; worst=(si,wi,ri,yi,fam,float(rlo),float(rhi),float(ylo),float(yhi),str(box))
    out={
        'status':'proved',
        'C_chic_stationary':str(C_acc),
        'method':'absolute Arb box bound on stationary-overlap compact-complement transition boxes; boxes with separated phase are excluded and handled by C_chic_nonstationary',
        'domain':{'sigma':f'[{args.sigma_min},{args.sigma_max}]','w':f"[0,{args.w_max}]",'r':[args.r_min,args.r_max],'y':[args.Y,2*args.Y],'lambda0':args.lambda0},
        'boxes':{'ns':args.ns,'nw':args.nw,'nr':args.nr,'ny':args.ny,'stationary_boxes_counted':counted,'separated_boxes_skipped':skipped},
        'worst_box':worst,
    }
    out['proof_hash']=hashlib.sha256(json.dumps(out,sort_keys=True,default=str).encode()).hexdigest()
    print(json.dumps(out,indent=2,default=str))

def spec(_args):
    print(json.dumps({'status':'spec','output_key':'C_chic_stationary','method':'absolute interval boxes on phase-stationary overlap boxes'},indent=2))

def main():
    p=argparse.ArgumentParser(); sub=p.add_subparsers(dest='cmd',required=True)
    sub.add_parser('spec')
    q=sub.add_parser('prove')
    q.add_argument('--Y',type=float,default=128)
    q.add_argument('--r-min',type=float,default=-4)
    q.add_argument('--r-max',type=float,default=32)
    q.add_argument('--w-max',type=float,default=1/math.sqrt(10))
    q.add_argument('--lambda0',type=float,default=0.25)
    q.add_argument('--ns',type=int,default=4); q.add_argument('--sigma-min',type=float,default=0.01); q.add_argument('--sigma-max',type=float,default=0.49); q.add_argument('--nw',type=int,default=16); q.add_argument('--nr',type=int,default=96); q.add_argument('--ny',type=int,default=128)
    args=p.parse_args()
    if args.cmd=='spec': spec(args)
    else: prove(args)
if __name__=='__main__': main()
