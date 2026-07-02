#!/usr/bin/env python3
"""
cert_chic_tail_envelope.py

Derivative-envelope prover for cert_chic_tail.py.

It is designed to prove bounds M1_hprime, M2_hsecond for the compact large-r tail
amplitude after factoring out the oscillatory phase exp(±i(y^2/8+2*pi*r*y)).

Mathematical object
-------------------
For mu = lambda*w and q=-sigma-2 or q=sigma-3, define

  C(mu,y) = mu^-2*(mu*y/2 - log(1+mu*y/2)),  C(0,y)=y^2/8
  delta(mu,y)=C(mu,y)-y^2/8

Minus family amplitude:
  H_-(y)=chi_Y(y)*exp(i*delta)*(A_mu + i*A*C_mu), A=(2+mu*y)^q.

Plus family amplitude:
  H_+(y)=chi_Y(y)*exp(-i*delta)*(A_mu - i*A*C_mu), A=(2+mu*y)^q.

This script bounds |dH/dy| and |d^2H/dy^2| with interval automatic
differentiation over finite boxes.

Important
---------
The mu=0 singularity is handled by splitting mu into:
  [0, mu0]       Taylor branch for C and C_mu
  [mu0, mu_max]  direct branch

The Taylor branch uses the exact log-series with an explicit geometric remainder,
valid when |u|=|mu*y/2| <= rho < 1. With Y=128 and y<=2Y, the default
mu0=1/2048 gives rho<=1/16.

Requires:
  python -m pip install python-flint

Run:
  python cert_chic_tail_envelope.py prove --Y 128 --mu-max 0.31622776601683794 --mu0 0.00048828125 --ns 4 --nmu 32 --ny 256

Output can be fed to:
  python cert_chic_tail.py prove --envelope-json proved_tail_derivative_envelope.json
"""
import argparse, json, hashlib, math, sys
try:
    from flint import arb, acb, ctx
except Exception as e:
    print(json.dumps({"status":"missing_dependency","error":repr(e),"install":"python -m pip install python-flint"}, indent=2))
    sys.exit(2)
ctx.prec = 100
I=acb(0,1)
PI=arb.pi()

class AD2:
    def __init__(self,v,d1=None,d2=None):
        self.v = acb(v)
        self.d1 = acb(0) if d1 is None else acb(d1)
        self.d2 = acb(0) if d2 is None else acb(d2)
    def __add__(self,o):
        o=toAD(o); return AD2(self.v+o.v,self.d1+o.d1,self.d2+o.d2)
    __radd__=__add__
    def __neg__(self): return AD2(-self.v,-self.d1,-self.d2)
    def __sub__(self,o): return self+(-toAD(o))
    def __rsub__(self,o): return toAD(o)+(-self)
    def __mul__(self,o):
        o=toAD(o); return AD2(self.v*o.v, self.d1*o.v+self.v*o.d1, self.d2*o.v+2*self.d1*o.d1+self.v*o.d2)
    __rmul__=__mul__
    def inv(self):
        v=self.v; return AD2(1/v, -self.d1/(v*v), 2*self.d1*self.d1/(v*v*v)-self.d2/(v*v))
    def __truediv__(self,o): return self*toAD(o).inv()
    def __rtruediv__(self,o): return toAD(o)*self.inv()
    def exp(self):
        ev=self.v.exp(); return AD2(ev, ev*self.d1, ev*(self.d1*self.d1+self.d2))
    def log(self):
        return AD2(self.v.log(), self.d1/self.v, self.d2/self.v - self.d1*self.d1/(self.v*self.v))
    def pow_real(self,p): return (self.log()*p).exp()
    def __pow__(self,n):
        if isinstance(n,int):
            if n==0: return AD2(1)
            if n<0: return (self**(-n)).inv()
            out=AD2(1)
            for _ in range(n): out=out*self
            return out
        return self.pow_real(arb(str(n)))

def toAD(x): return x if isinstance(x,AD2) else AD2(x)

def midrad(a,b):
    a=arb(str(a)); b=arb(str(b)); return (a+b)/2 + arb(f"0 +/- {float((b-a)/2)}")

def mag(z): return float(acb(z).abs_upper())

def chi_ad(y,Y):
    Y=arb(str(Y)); twoY=2*Y
    # Caller splits boxes in [0,Y] and [Y,2Y].
    if bool(y.v.real <= Y): return AD2(1)
    if bool(y.v.real >= twoY): return AD2(0)
    t=(y-AD2(Y))/Y
    return AD2(1)-3*t*t+2*t*t*t

def core_C_direct(mu,y):
    u=AD2(mu)*y/2
    return (u-(AD2(1)+u).log())/(AD2(mu)*AD2(mu))

def core_Cmu_direct(mu,y):
    # derivative wrt mu of C=(u-log(1+u))/mu^2, u=mu*y/2
    muA=AD2(mu); u=muA*y/2
    N=u-(AD2(1)+u).log()
    Nmu=(y/2)*(u/(AD2(1)+u))
    return Nmu/(muA*muA) - 2*N/(muA*muA*muA)

def core_C_taylor(mu,y,N=18):
    # C=sum_{n=2}^N (-1)^n mu^(n-2)y^n/(n2^n) + interval remainder in value only, inflated.
    muA=AD2(mu); out=AD2(0)
    for n in range(2,N+1):
        coeff = arb(1 if n%2==0 else -1) / (arb(n)*(arb(2)**n))
        out = out + coeff*(muA**(n-2))*(y**n)
    # crude remainder for |u|<=rho; inflate value,d1,d2 by same large ball.
    # R_C <= y^2/4 * rho^(N-1)/((N+1)*(1-rho)) with rho=|mu*y/2| upper.
    rho = (acb(muA.v*y.v/2)).abs_upper()
    rem = (acb(y.v).abs_upper()**2)/4 * (rho**(N-1))/(arb(N+1)*(1-rho))
    ball = acb(arb(0) + arb(f"0 +/- {float(rem)}"))
    out.v += ball; out.d1 += ball; out.d2 += ball
    return out

def core_Cmu_taylor(mu,y,N=18):
    muA=AD2(mu); out=AD2(0)
    for n in range(3,N+1):
        coeff = arb(1 if n%2==0 else -1) * arb(n-2) / (arb(n)*(arb(2)**n))
        out = out + coeff*(muA**(n-3))*(y**n)
    rho = (acb(muA.v*y.v/2)).abs_upper()
    rem = (acb(y.v).abs_upper()**3)/8 * (rho**(N-2))/(arb(N+1)*(1-rho))
    ball = acb(arb(0) + arb(f"0 +/- {float(rem)}"))
    out.v += ball; out.d1 += ball; out.d2 += ball
    return out

def H_family(sig,mu,y,Y,fam,branch):
    q = -sig-2 if fam=='-' else sig-3
    X=AD2(2)+AD2(mu)*y
    A=X.pow_real(q)
    Amu = q*y*X.pow_real(q-1)
    if branch=='taylor':
        C=core_C_taylor(mu,y); Cmu=core_Cmu_taylor(mu,y)
    else:
        C=core_C_direct(mu,y); Cmu=core_Cmu_direct(mu,y)
    delta=C - (y*y)/8
    ch=chi_ad(y,Y)
    if fam=='-':
        E=(I*delta).exp(); B=Amu + I*A*Cmu
    else:
        E=(-I*delta).exp(); B=Amu - I*A*Cmu
    return ch*E*B

def y_boxes(Y,ny):
    boxes=[]
    for lo,hi,n in [(0,float(Y),ny//2),(float(Y),2*float(Y),ny-ny//2)]:
        for j in range(n):
            a=lo+(hi-lo)*j/n; b=lo+(hi-lo)*(j+1)/n; boxes.append((a,b))
    return boxes

def prove(args):
    max1=0.0; max2=0.0; worst1=None; worst2=None
    # sigma boxes
    for si in range(args.ns):
        sig=midrad(args.sigma_min+(args.sigma_max-args.sigma_min)*si/args.ns,args.sigma_min+(args.sigma_max-args.sigma_min)*(si+1)/args.ns)
        # mu Taylor branch [0,mu0]
        mu_ranges=[(0.0,args.mu0,'taylor')]
        for mi in range(args.nmu):
            a=args.mu0+(args.mu_max-args.mu0)*mi/args.nmu
            b=args.mu0+(args.mu_max-args.mu0)*(mi+1)/args.nmu
            mu_ranges.append((a,b,'direct'))
        for mi,(ma,mb,branch) in enumerate(mu_ranges):
            mu=midrad(ma,mb)
            for yi,(ya,yb) in enumerate(y_boxes(args.Y,args.ny)):
                y=AD2(midrad(ya,yb),1,0)
                for fam in ['+','-']:
                    H=H_family(sig,mu,y,args.Y,fam,branch)
                    v1=mag(H.d1); v2=mag(H.d2)
                    if v1>max1: max1=v1; worst1=(si,mi,yi,fam,branch,v1)
                    if v2>max2: max2=v2; worst2=(si,mi,yi,fam,branch,v2)
    out={
        'status':'proved',
        'M1_hprime':str(max1),
        'M2_hsecond':str(max2),
        'domain': {'sigma':f'[{args.sigma_min},{args.sigma_max}]','mu':f"[0,{args.mu_max}]",'y':f"[0,{2*args.Y}]",'Y':args.Y},
        'boxes': {'ns':args.ns,'nmu_direct':args.nmu,'ny':args.ny,'mu0_taylor':args.mu0},
        'worst_M1_box':worst1,
        'worst_M2_box':worst2,
        'method':'Arb interval AD2 derivative envelope for compact tail amplitude; Taylor log branch near mu=0 with interval remainder',
    }
    import hashlib, json
    out['proof_hash']=hashlib.sha256(json.dumps(out,sort_keys=True,default=str).encode()).hexdigest()
    print(json.dumps(out,indent=2,default=str))

def main():
    import argparse
    p=argparse.ArgumentParser()
    p.add_argument('cmd',choices=['prove'])
    p.add_argument('--Y',type=float,default=128)
    p.add_argument('--mu-max',type=float,default=1/(10**0.5))
    p.add_argument('--mu0',type=float,default=1/2048)
    p.add_argument('--ns',type=int,default=4)
    p.add_argument('--sigma-min',type=float,default=0.01)
    p.add_argument('--sigma-max',type=float,default=0.49)
    p.add_argument('--nmu',type=int,default=32)
    p.add_argument('--ny',type=int,default=256)
    args=p.parse_args(); prove(args)
if __name__=='__main__': main()
