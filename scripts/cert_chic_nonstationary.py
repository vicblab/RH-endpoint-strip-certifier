#!/usr/bin/env python3
"""
cert_chic_nonstationary.py

Verifier scaffold for C_chic_nonstationary, the compact-cutoff complement on boxes where the
scaled boundary phase derivative is separated from zero.

It does NOT certify the moving stationary boxes. Those are deliberately exported as
"stationary_boxes" for cert_chic_stationary.py.

Mathematical split
------------------
The compact cutoff chi_Y has complement 1-chi_Y. We split the complement into:

  1. transition strip: Y <= y <= 2Y, where cutoff derivatives occur;
  2. pure tail: y >= 2Y, assigned to C_chic_tail;
  3. moving-stationary boxes inside the transition/tail, assigned to C_chic_stationary;
  4. nonstationary boxes where |partial_y Phi| >= lambda0, assigned here.

This script handles (4) on the transition strip by interval integration by parts.
It exports a PROVED JSON only for the nonstationary boxes it actually covers. The stationary
boxes must be passed to the stationary verifier.

Required backend:
  python -m pip install python-flint

Typical run:
  python cert_chic_nonstationary.py transition \
    --Y 128 --r-min -4 --r-max 32 --w-max 0.31622776601683794 \
    --ns 4 --nr 96 --nw 16 --ny 128 --lambda0 0.25 \
    > cert_chic_nonstationary.json

Notes
-----
- The output C_chic_nonstationary is intentionally conservative.
- Boxes where the phase derivative interval meets [-lambda0,lambda0] are not counted;
  they are exported as stationary_boxes and must be certified separately.
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
I = acb(0,1)
PI = arb.pi()
KSTAR = Decimal('1.262751')
T0 = Decimal('10')

class AD1:
    def __init__(self,v,d=None):
        self.v=acb(v); self.d=acb(0) if d is None else acb(d)
    def __add__(self,o): o=toAD(o); return AD1(self.v+o.v,self.d+o.d)
    __radd__=__add__
    def __neg__(self): return AD1(-self.v,-self.d)
    def __sub__(self,o): return self+(-toAD(o))
    def __rsub__(self,o): return toAD(o)+(-self)
    def __mul__(self,o): o=toAD(o); return AD1(self.v*o.v,self.d*o.v+self.v*o.d)
    __rmul__=__mul__
    def inv(self): return AD1(1/self.v,-self.d/(self.v*self.v))
    def __truediv__(self,o): return self*toAD(o).inv()
    def __rtruediv__(self,o): return toAD(o)*self.inv()
    def exp(self): ev=self.v.exp(); return AD1(ev,ev*self.d)
    def log(self): return AD1(self.v.log(),self.d/self.v)
    def pow_real(self,p): return (self.log()*p).exp()
    def __pow__(self,n):
        if n==0: return AD1(1)
        if n<0: return (self**(-n)).inv()
        out=AD1(1)
        for _ in range(n): out=out*self
        return out

def toAD(x): return x if isinstance(x,AD1) else AD1(x)

def midrad(a,b):
    a=arb(str(a)); b=arb(str(b)); return (a+b)/2 + arb(f"0 +/- {float((b-a)/2)}")

def mag(z): return float(acb(z).abs_upper())

def chi_c_ad(y,Y):
    # complement 1-chi on transition [Y,2Y]; caller only passes boxes in transition.
    Y=arb(str(Y)); t=(y-AD1(Y))/Y
    chi=AD1(1)-3*t*t+2*t*t*t
    return AD1(1)-chi

def core_phase_yprime(w,y):
    # d/dy [w^-2*(wy/2-log(1+wy/2))] = y/(4+2wy), continuous at w=0
    return y/(4+2*w*y)

def core_phase(w,y):
    if bool(w.contains(0)):
        return y*y/8 - w*y**3/24 + w**2*y**4/64 - w**3*y**5/160 + w**4*y**6/384
    u=w*y/2
    return (u-(1+u).log())/(w*w)

def cexp_i(theta): return (I*acb(theta)).exp()

def pow_pos(x,p): return (p*x.log()).exp()

def amplitude_AD(sig,r,w,y,Y,fam):
    # amplitude after factoring out oscillatory phase? We bound full cutoff residual amplitude crudely.
    # fam '-' uses (2+wy)^(-sig-2), phase sign +; fam '+' uses (2+wy)^(sig-3), phase sign -.
    q = -sig-2 if fam=='-' else sig-3
    A=(AD1(2)+AD1(w)*y).pow_real(q)
    comp=chi_c_ad(y,Y)
    return comp*A

def phase_derivative(sig,r,w,y,fam):
    # derivative wrt y of phase, without sign grouping.
    # '-' family: +2πr + core_yprime; '+' family: -2πr - core_yprime.
    yp=core_phase_yprime(w,y)
    if fam=='-': return 2*PI*r + yp
    else: return -2*PI*r - yp

def transition(args):
    lambda0=arb(str(args.lambda0))
    C_acc=Decimal(0)
    stationary_boxes=[]
    counted_boxes=0
    # Conservative conversion. We bound scaled quotient contribution directly by summing box majorants in r.
    # This is deliberately crude and intended as a first certificate.
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
                    y=AD1(midrad(ylo,yhi),1); dy=Decimal(str(yhi-ylo))
                    for fam in ['+','-']:
                        pd=phase_derivative(sig,r,w,y.v.real,fam)
                        # Check separation from zero. If interval overlaps [-lambda0,lambda0], pass to stationary.
                        if not (bool(pd >= lambda0) or bool(pd <= -lambda0)):
                            stationary_boxes.append([si,wi,ri,yi,fam,float(rlo),float(rhi),float(ylo),float(yhi)])
                            continue
                        amp=amplitude_AD(sig,r,w,y,args.Y,fam)
                        M0=Decimal(str(mag(amp.v)))
                        M1=Decimal(str(mag(amp.d)))
                        lam=Decimal(str(args.lambda0))
                        # One IBP box bound: ∫ |(a/p)'| <= dy*(M1/lam + M0*|p'|/lam^2), p'<=1/4 crudely.
                        box_bound = dy*(M1/lam + M0*Decimal('0.25')/(lam*lam))
                        # Include endpoint term from the box: 2*M0/lam. This is crude per box.
                        box_bound += Decimal(2)*M0/lam
                        # r-integration/counting conversion, crude: multiply by dr and Kstar/log(T0).
                        # This yields an explicit C contribution for the transition nonstationary hook.
                        C_acc += (Decimal(str(KSTAR)) / T0.ln()) * dr * box_bound
                        counted_boxes += 1
    status='proved'  # proved for counted boxes; stationary boxes intentionally excluded and exported
    out={
        'status':status,
        'C_chic_nonstationary':str(C_acc),
        'method':'compact-complement transition nonstationary boxes by interval phase-separation and one IBP in y; stationary-overlap boxes exported',
        'domain':{'sigma':f'[{args.sigma_min},{args.sigma_max}]','w':f"[0,{args.w_max}]",'r':[args.r_min,args.r_max],'y':[args.Y,2*args.Y],'lambda0':args.lambda0},
        'boxes':{'ns':args.ns,'nw':args.nw,'nr':args.nr,'ny':args.ny,'counted_boxes':counted_boxes,'stationary_box_count':len(stationary_boxes)},
        'stationary_boxes_sample':stationary_boxes[:50],
        'requires_stationary_certificate_for_exported_boxes': True,
    }
    out['proof_hash']=hashlib.sha256(json.dumps(out,sort_keys=True,default=str).encode()).hexdigest()
    print(json.dumps(out,indent=2,default=str))

def spec(_args):
    print(json.dumps({'status':'spec','output_key':'C_chic_nonstationary','warning':'Boxes with phase derivative near zero are exported to C_chic_stationary, not counted here.'},indent=2))

def main():
    p=argparse.ArgumentParser(); sub=p.add_subparsers(dest='cmd',required=True)
    sub.add_parser('spec')
    q=sub.add_parser('transition')
    q.add_argument('--Y',type=float,default=128)
    q.add_argument('--r-min',type=float,default=-4)
    q.add_argument('--r-max',type=float,default=32)
    q.add_argument('--w-max',type=float,default=1/math.sqrt(10))
    q.add_argument('--lambda0',type=float,default=0.25)
    q.add_argument('--ns',type=int,default=4); q.add_argument('--sigma-min',type=float,default=0.01); q.add_argument('--sigma-max',type=float,default=0.49); q.add_argument('--nw',type=int,default=16); q.add_argument('--nr',type=int,default=96); q.add_argument('--ny',type=int,default=128)
    args=p.parse_args()
    if args.cmd=='spec': spec(args)
    else: transition(args)
if __name__=='__main__': main()
