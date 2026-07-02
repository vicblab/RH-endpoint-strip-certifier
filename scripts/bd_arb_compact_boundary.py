#!/usr/bin/env python3
"""
bd_arb_compact_boundary.py

Compact-y boundary verifier for the outer-strip companion certificate.

This implements Option A: replace the non-compact boundary/Fresnel packet by a compact
packet in the scaled boundary variable y.

Cutoff convention:
  chi_Y(y) = 1                         for 0 <= y <= Y
           = 1 - 3t^2 + 2t^3           for Y <= y <= 2Y, t=(y-Y)/Y
           = 0                         for y >= 2Y

The verifier integrates only over [0,2Y]. There is therefore no boundary y-tail.
The complement (1-chi_Y) is NOT silently discarded: it must be assigned to the
far/nonstationary complement certificate. This script reports this explicitly.

Requires:
  python -m pip install python-flint

Examples:
  python bd_arb_compact_boundary.py selftest
  python bd_arb_compact_boundary.py grp-compact --tau-min 10 --tau-max 25 --Y 128 --ns 4 --nt 8 --ny 256
  python bd_arb_compact_boundary.py tail-compact --w-min 0.00390625 --w-max 0.31622776601683794 --R 32 --Y 128 --ns 4 --nw 8 --nr 48 --ny 256
  python bd_arb_compact_boundary.py central-taylor-compact --w0 0.00390625 --Y 128 --ns 4 --nw 4 --nr 48 --ny 256 --nl 8
  python bd_arb_compact_boundary.py tail-taylor-compact --w0 0.00390625 --R 32 --Y 128 --ns 4 --nw 4 --nr 48 --ny 256 --nl 8
"""
import sys, json, math, argparse, hashlib
try:
    from flint import arb, acb, ctx
except Exception as e:
    print(json.dumps({"status":"missing_dependency","error":repr(e),"install":"python -m pip install python-flint"}, indent=2))
    sys.exit(2)

ctx.prec = 100
I = acb(0, 1)
PI = arb.pi()
KSTAR10 = arb('1.262751')
CONV_TAIL = 8*KSTAR10*(arb(1)/50 + arb(1)/(125*arb(10).sqrt()))/arb(10).log()
BWIDTH = arb(4)


def midrad(a,b):
    a=arb(str(a)); b=arb(str(b)); rad=float((b-a)/2)
    return (a+b)/2 + arb(f"0 +/- {rad}")

def mag_acb(z):
    return float(z.abs_upper())

def cexp_i(theta):
    return (I*acb(theta)).exp()

def pow_pos(x,p):
    return (p*x.log()).exp()

def chi(y,Y):
    # Smooth cubic cutoff as an Arb expression. Branching is by interval location in caller.
    y0 = arb(str(Y)); y1 = 2*y0
    # Caller should pass boxes entirely in [0,Y], [Y,2Y], or outside.
    if bool(y <= y0):
        return arb(1)
    if bool(y >= y1):
        return arb(0)
    t=(y-y0)/y0
    return 1 - 3*t*t + 2*t*t*t

def theta_model(r,y):
    return y*y/8 + 2*PI*r*y

def phase_core(w,y):
    if bool(w.contains(0)):
        return y*y/8 - w*y**3/24 + w**2*y**4/64 - w**3*y**5/160 + w**4*y**6/384 - w**5*y**7/896
    u=w*y/2
    return (u-(1+u).log())/(w*w)

def phase_core_wprime(w,y):
    if bool(w.contains(0)):
        return -y**3/24 + w*y**4/32 - 3*w**2*y**5/160 + w**3*y**6/96 - 5*w**4*y**7/896
    u=w*y/2
    N=u-(1+u).log(); Np=(y/2)*(u/(1+u))
    return Np/(w*w) - 2*N/(w**3)

def D_minus_integrand(sig,r,w,y,Y):
    A=pow_pos(2+w*y, -sig-2); A0=pow_pos(arb(2), -sig-2)
    return acb(chi(y,Y))*(acb(A)*cexp_i(2*PI*r*y + phase_core(w,y)) - acb(A0)*cexp_i(theta_model(r,y)))

def D_plus_integrand(sig,r,w,y,Y):
    A=pow_pos(2+w*y, sig-3); A0=pow_pos(arb(2), sig-3)
    return acb(chi(y,Y))*(acb(A)*cexp_i(-2*PI*r*y - phase_core(w,y)) - acb(A0)*cexp_i(-theta_model(r,y)))

def fw_minus(sig,r,w,y,Y):
    q=-sig-2; X=2+w*y
    A=pow_pos(X,q); Aw=q*y*pow_pos(X,q-1); phiw=phase_core_wprime(w,y)
    return acb(chi(y,Y))*cexp_i(2*PI*r*y + phase_core(w,y))*(acb(Aw)+I*acb(A)*acb(phiw))

def fw_plus(sig,r,w,y,Y):
    q=sig-3; X=2+w*y
    A=pow_pos(X,q); Aw=q*y*pow_pos(X,q-1); phiw=-phase_core_wprime(w,y)
    return acb(chi(y,Y))*cexp_i(-2*PI*r*y - phase_core(w,y))*(acb(Aw)+I*acb(A)*acb(phiw))

def y_subintervals(Y,ny):
    # Split [0,2Y] with a break at Y so boxes don't straddle cutoff branch.
    out=[]
    for part_lo, part_hi, n in [(0.0, float(Y), ny//2), (float(Y), 2*float(Y), ny-ny//2)]:
        for j in range(n):
            a=part_lo+(part_hi-part_lo)*j/n; b=part_lo+(part_hi-part_lo)*(j+1)/n
            out.append((a,b))
    return out

def integrate_D_boxed(sig,r,w,fam,Y,ny):
    total=acb(0)
    for a,b in y_subintervals(Y,ny):
        ybox=midrad(a,b); dy=arb(str(b-a))
        val=D_plus_integrand(sig,r,w,ybox,Y) if fam=='+' else D_minus_integrand(sig,r,w,ybox,Y)
        total += val*acb(dy)
    return total

def integrate_H_boxed(sig,r,w,fam,Y,ny,nl):
    total=acb(0)
    for a,b in y_subintervals(Y,ny):
        ybox=midrad(a,b); dy=arb(str(b-a))
        for k in range(nl):
            la=k/nl; lb=(k+1)/nl
            lbox=midrad(la,lb); dl=arb(str(lb-la)); ww=lbox*w
            val=fw_plus(sig,r,ww,ybox,Y) if fam=='+' else fw_minus(sig,r,ww,ybox,Y)
            total += val*acb(dy*dl)
    return total

def tail_compact(args):
    maxv=0.0; worst=None
    for si in range(args.ns):
        sig=midrad(args.sigma_min+(args.sigma_max-args.sigma_min)*si/args.ns,args.sigma_min+(args.sigma_max-args.sigma_min)*(si+1)/args.ns)
        for wi in range(args.nw):
            w=midrad(args.w_min+(args.w_max-args.w_min)*wi/args.nw,args.w_min+(args.w_max-args.w_min)*(wi+1)/args.nw)
            for ri in range(args.nr):
                r=midrad(4+(args.R-4)*ri/args.nr,4+(args.R-4)*(ri+1)/args.nr)
                for fam in ['+','-']:
                    D=integrate_D_boxed(sig,r,w,fam,args.Y,args.ny)
                    scaled=D*acb(pow_pos(1+r, arb(3))/w)
                    val=mag_acb(scaled)
                    if val>maxv: maxv=val; worst=(si,wi,ri,fam,val)
    emit({"status":"compact_boundary_fragment","subsector":"bd_high_tail_compact_finite_w_r","w_range":[args.w_min,args.w_max],"r_range":[4,args.R],"Y":args.Y,"A_tail_sc_compact_fragment_upper":maxv,"worst_box":worst,"missing_for_full_certificate":["r>=R compact-complement tail reassignment"]})

def grp_compact(args):
    maxv=0.0; worst=None
    for si in range(args.ns):
        sig=midrad(args.sigma_min+(args.sigma_max-args.sigma_min)*si/args.ns,args.sigma_min+(args.sigma_max-args.sigma_min)*(si+1)/args.ns)
        for ti in range(args.nt):
            lo=args.tau_min+(args.tau_max-args.tau_min)*ti/args.nt; hi=args.tau_min+(args.tau_max-args.tau_min)*(ti+1)/args.nt
            tau=midrad(lo,hi); w=1/tau.sqrt(); alpha_lo=lo/(4*math.pi); alpha_hi=hi/(4*math.pi); root_hi=math.sqrt(hi)
            mmin=max(1,int(math.floor(alpha_lo-4*root_hi))-2); mmax=int(math.ceil(alpha_hi+4*root_hi))+2
            s=acb(sig,tau); total=acb(0)
            for m in range(mmin,mmax+1):
                ma=arb(m); r=(ma-tau/(4*PI))/tau.sqrt()
                Dp=integrate_D_boxed(sig,r,w,'+',args.Y,args.ny); Dm=integrate_D_boxed(sig,r,w,'-',args.Y,args.ny)
                Delta_p=acb(tau**(-arb('0.5')))*cexp_i(-4*PI*ma+tau*arb(2).log())*Dp
                Delta_m=acb(tau**(-arb('0.5')))*cexp_i( 4*PI*ma-tau*arb(2).log())*Dm
                total += (s*(s-2)*Delta_p + (1-s)*(s+1)*Delta_m)/(m*m)
            Rj=total/acb(4*PI*PI); norm=Rj.abs_upper()*tau.sqrt()/tau.log(); nf=float(norm)
            if nf>maxv: maxv=nf; worst=(si,ti,mmin,mmax,nf)
    emit({"status":"compact_boundary_fragment","subsector":"bd_grouped_central_compact_finite_tau","tau_range":[args.tau_min,args.tau_max],"Y":args.Y,"A_bd_grp_compact_fragment_upper":maxv,"worst_box":worst,"missing_for_full_certificate":["tau beyond tau_max or w->0 compact Taylor"]})

def tail_taylor_compact(args):
    maxv=0.0; worst=None
    for si in range(args.ns):
        sig=midrad(args.sigma_min+(args.sigma_max-args.sigma_min)*si/args.ns,args.sigma_min+(args.sigma_max-args.sigma_min)*(si+1)/args.ns)
        for wi in range(args.nw):
            w=midrad(0+args.w0*wi/args.nw,0+args.w0*(wi+1)/args.nw)
            for ri in range(args.nr):
                r=midrad(4+(args.R-4)*ri/args.nr,4+(args.R-4)*(ri+1)/args.nr)
                for fam in ['+','-']:
                    H=integrate_H_boxed(sig,r,w,fam,args.Y,args.ny,args.nl)
                    scaled=H*acb(pow_pos(1+r,arb(3)))
                    val=mag_acb(scaled)
                    if val>maxv: maxv=val; worst=(si,wi,ri,fam,val)
    emit({"status":"compact_boundary_fragment","subsector":"bd_high_tail_compact_taylor_w0","w_range":[0,args.w0],"r_range":[4,args.R],"Y":args.Y,"A_tail_sc_compact_taylor_upper":maxv,"A_bd_tail_compact_taylor_upper":float(CONV_TAIL)*maxv,"worst_box":worst,"missing_for_full_certificate":["r>=R compact-complement tail reassignment"]})

def central_taylor_compact(args):
    maxv=0.0; worst=None
    for si in range(args.ns):
        sig=midrad(args.sigma_min+(args.sigma_max-args.sigma_min)*si/args.ns,args.sigma_min+(args.sigma_max-args.sigma_min)*(si+1)/args.ns)
        for wi in range(args.nw):
            w=midrad(0+args.w0*wi/args.nw,0+args.w0*(wi+1)/args.nw)
            for ri in range(args.nr):
                r=midrad(-4+8*ri/args.nr,-4+8*(ri+1)/args.nr)
                for fam in ['+','-']:
                    H=integrate_H_boxed(sig,r,w,fam,args.Y,args.ny,args.nl)
                    val=mag_acb(H)
                    if val>maxv: maxv=val; worst=(si,wi,ri,fam,val)
    w0=arb(str(args.w0)); c0=1/(4*PI)-4*w0
    conv=(1/(4*PI*PI))*2*KSTAR10*arb(str(maxv))/(c0*c0)*(8+3*w0)/((1/(w0*w0)).log())
    emit({"status":"compact_boundary_fragment","subsector":"bd_grouped_central_compact_taylor_w0","w_range":[0,args.w0],"r_range":[-4,4],"Y":args.Y,"K_central_H_compact_sup_upper":maxv,"A_bd_grp_compact_taylor_upper":float(conv.abs_upper()),"worst_box":worst,"missing_for_full_certificate":[]})

def emit(obj):
    obj['cutoff']={"type":"cubic","chi":"1 on [0,Y], 1-3t^2+2t^3 on [Y,2Y], 0 on [2Y,infty]"}
    obj['warning']='Compact boundary removes boundary y-tail; complement must be covered by far/nonstationary reassignment certificate.'
    obj['hash']=hashlib.sha256(json.dumps(obj,sort_keys=True,default=str).encode()).hexdigest()
    print(json.dumps(obj,indent=2,default=str))

def selftest(): emit({"status":"ok","backend":"python-flint Arb/acb","CONV_TAIL":str(CONV_TAIL)})

def main():
    p=argparse.ArgumentParser(); sub=p.add_subparsers(dest='cmd',required=True)
    sub.add_parser('selftest')
    g=sub.add_parser('grp-compact'); g.add_argument('--tau-min',type=float,required=True); g.add_argument('--tau-max',type=float,required=True); g.add_argument('--Y',type=float,default=128); g.add_argument('--ns',type=int,default=4); g.add_argument('--sigma-min',type=float,default=0.01); g.add_argument('--sigma-max',type=float,default=0.49); g.add_argument('--nt',type=int,default=8); g.add_argument('--ny',type=int,default=256)
    t=sub.add_parser('tail-compact'); t.add_argument('--w-min',type=float,required=True); t.add_argument('--w-max',type=float,required=True); t.add_argument('--R',type=float,default=32); t.add_argument('--Y',type=float,default=128); t.add_argument('--ns',type=int,default=4); t.add_argument('--sigma-min',type=float,default=0.01); t.add_argument('--sigma-max',type=float,default=0.49); t.add_argument('--nw',type=int,default=8); t.add_argument('--nr',type=int,default=48); t.add_argument('--ny',type=int,default=256)
    ct=sub.add_parser('central-taylor-compact'); ct.add_argument('--w0',type=float,required=True); ct.add_argument('--Y',type=float,default=128); ct.add_argument('--ns',type=int,default=4); ct.add_argument('--sigma-min',type=float,default=0.01); ct.add_argument('--sigma-max',type=float,default=0.49); ct.add_argument('--nw',type=int,default=4); ct.add_argument('--nr',type=int,default=48); ct.add_argument('--ny',type=int,default=256); ct.add_argument('--nl',type=int,default=8)
    tt=sub.add_parser('tail-taylor-compact'); tt.add_argument('--w0',type=float,required=True); tt.add_argument('--R',type=float,default=32); tt.add_argument('--Y',type=float,default=128); tt.add_argument('--ns',type=int,default=4); tt.add_argument('--sigma-min',type=float,default=0.01); tt.add_argument('--sigma-max',type=float,default=0.49); tt.add_argument('--nw',type=int,default=4); tt.add_argument('--nr',type=int,default=48); tt.add_argument('--ny',type=int,default=256); tt.add_argument('--nl',type=int,default=8)
    args=p.parse_args()
    if args.cmd=='selftest': selftest()
    elif args.cmd=='grp-compact': grp_compact(args)
    elif args.cmd=='tail-compact': tail_compact(args)
    elif args.cmd=='central-taylor-compact': central_taylor_compact(args)
    elif args.cmd=='tail-taylor-compact': tail_taylor_compact(args)
if __name__=='__main__': main()
