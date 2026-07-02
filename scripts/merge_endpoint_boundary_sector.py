#!/usr/bin/env python3
"""
merge_endpoint_boundary_sector.py

Endpoint boundary-sector merger for the compact-boundary fragments and compact-
complement hooks produced by the endpoint rerun.

It is intentionally transparent:
  C_bd = A_grouped_compact + A_tail_compact + C_chic_tail
         + C_chic_nonstationary + C_chic_stationary

For finite compact high-tail fragments reporting A_tail_sc_compact_fragment_upper,
it applies the same high-tail conversion factor recorded in the companion paper.
For Taylor high-tail fragments reporting A_bd_tail_compact_taylor_upper, it uses
that converted value directly.

This script does not alter the mathematical proof; it only merges proved/fragment
outputs according to the existing compact-boundary ledger model.
"""
import argparse, json, glob, hashlib
from decimal import Decimal, getcontext
from pathlib import Path
getcontext().prec = 100

CONV_TAIL = Decimal('0.098843879459434312935638821476530184904991788627323')

def D(x): return Decimal(str(x))
def h(obj): return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":"), default=str).encode()).hexdigest()
def load(p): return json.loads(Path(p).read_text(encoding='utf-8'))

def expand(patterns):
    out=[]
    for pat in patterns:
        hits=glob.glob(pat)
        out.extend(hits if hits else [pat])
    return sorted(set(out))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--fragments', nargs='+', required=True)
    ap.add_argument('--tail-json', required=True)
    ap.add_argument('--nonstationary-json', required=True)
    ap.add_argument('--stationary-json', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--sigma-min', default='0')
    ap.add_argument('--sigma-max', default='0.01')
    ap.add_argument('--tau-min', default='10')
    ap.add_argument('--safe-eps', default='1e-12')
    ap.add_argument('--safe-pad', default='100')
    args=ap.parse_args()
    eps=D(args.safe_eps); pad=D(args.safe_pad)

    group_vals=[]; tail_vals=[]; deps={}
    for p in expand(args.fragments):
        obj=load(p)
        if obj.get('status') != 'compact_boundary_fragment':
            raise SystemExit(f'ERROR: {p}: expected compact_boundary_fragment, got {obj.get("status")}')
        dep={'file':p,'hash':obj.get('hash'),'subsector':obj.get('subsector')}
        if 'A_bd_grp_compact_fragment_upper' in obj:
            val=D(obj['A_bd_grp_compact_fragment_upper']); group_vals.append(val); dep['used_as']='grouped_compact'; dep['value']=str(val)
        if 'A_bd_grp_compact_taylor_upper' in obj:
            val=D(obj['A_bd_grp_compact_taylor_upper']); group_vals.append(val); dep['used_as']='grouped_compact_taylor'; dep['value']=str(val)
        if 'A_bd_tail_compact_taylor_upper' in obj:
            val=D(obj['A_bd_tail_compact_taylor_upper']); tail_vals.append(val); dep['used_as']='tail_compact_taylor_converted'; dep['value']=str(val)
        elif 'A_tail_sc_compact_fragment_upper' in obj:
            raw=D(obj['A_tail_sc_compact_fragment_upper']); val=CONV_TAIL*raw; tail_vals.append(val); dep['used_as']='tail_compact_finite_converted'; dep['raw_A_tail_sc']=str(raw); dep['conversion']=str(CONV_TAIL); dep['value']=str(val)
        deps[Path(p).name]=dep

    if not group_vals: raise SystemExit('ERROR: no grouped compact fragments found')
    if not tail_vals: raise SystemExit('ERROR: no compact tail fragments found')
    A_grp=max(group_vals); A_tail=max(tail_vals)

    hooks={}
    for name, path, keys in [
        ('chic_tail', args.tail_json, ['C_chic_tail','C_value_upper_safe','C_value_upper']),
        ('chic_nonstationary', args.nonstationary_json, ['C_chic_nonstationary','C_value_upper_safe','C_value_upper']),
        ('chic_stationary', args.stationary_json, ['C_chic_stationary','C_value_upper_safe','C_value_upper']),
    ]:
        obj=load(path)
        if obj.get('status') != 'proved':
            raise SystemExit(f'ERROR: {path}: hook status must be proved')
        val=None
        for k in keys:
            if k in obj:
                val=D(obj[k]); break
        if val is None: raise SystemExit(f'ERROR: {path}: no recognized constant key')
        hooks[name]={'file':path,'value':str(val),'proof_hash':obj.get('proof_hash')}

    C0 = A_grp + A_tail + sum(D(v['value']) for v in hooks.values())
    C0_safe = Decimal(0) if C0 == 0 else C0*(Decimal(1)+eps)+pad
    out={
        'status':'proved',
        'sector':'bd',
        'side':'J_projected',
        'method':'endpoint compact-boundary fragment max merge plus compact-complement hooks',
        'sigma_min':str(args.sigma_min),
        'sigma_max':str(args.sigma_max),
        'tau_min':str(args.tau_min),
        'C_value_upper':str(C0),
        'C_value_upper_safe':str(C0_safe),
        'C_sigma_derivative_upper':'0',
        'C_sigma_derivative_upper_safe':'0',
        'components':{
            'A_grouped_compact_upper':str(A_grp),
            'A_tail_compact_upper':str(A_tail),
            'compact_complement_hooks':hooks,
            'safe_rule':f'safe = raw*(1+{eps})+{pad}'
        },
        'dependencies':deps,
        'coverage_hash':h({'sector':'bd','sigma':[args.sigma_min,args.sigma_max],'tau_min':args.tau_min}),
        'normalization_hash':h({'sector':'bd','side':'J_projected','conversion_tail':str(CONV_TAIL)}),
    }
    out['proof_hash']=h(out)
    Path(args.output).write_text(json.dumps(out, indent=2, sort_keys=True, default=str), encoding='utf-8')
    print(json.dumps(out, indent=2, sort_keys=True, default=str))

if __name__=='__main__': main()
