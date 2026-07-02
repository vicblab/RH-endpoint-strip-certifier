#!/usr/bin/env python3
"""
merge_endpoint_certificate_safe.py
Strict one-channel endpoint merger for the endpoint-left quotient obstruction.

Purpose
-------
Merge seven proved endpoint sector JSON files for
    0 <= sigma <= 0.01, tau >= tau_min
and check the one-channel height inequality
    sqrt(T) > (C_endpoint / d_endpoint) * log(T).

This is NOT the inner two-channel merger. It sums only value residual constants.
Each sector input must have:
    status == "proved"
    sector in {band, bd, far, end, nonstat, floor, core}
    C_value_upper_safe
and must declare endpoint coverage compatible with sigma_min/sigma_max, either in
explicit fields or in a domain.sigma text field. The coverage check is deliberately
strict but allows old-style JSONs only when --allow-missing-coverage is supplied.
"""
from __future__ import annotations
import argparse, glob, json, hashlib, math, sys, re
from decimal import Decimal, getcontext
from pathlib import Path

getcontext().prec = 100
SECTORS = ['band', 'bd', 'far', 'end', 'nonstat', 'floor', 'core']
SIDES = {'Q', 'J_projected', 'zero', 'mixed', 'UNKNOWN'}

def D(x):
    return Decimal(str(x))

def hash_obj(obj) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()).hexdigest()

def load(p):
    with open(p, 'r', encoding='utf-8') as f:
        return json.load(f)

def req(o, k, p):
    if k not in o:
        raise SystemExit(f'ERROR: {p}: missing required field {k}')
    return o[k]

def expand(patterns):
    out=[]; missing=[]
    for pat in patterns:
        hits=glob.glob(pat)
        if hits:
            out += hits
        elif any(ch in pat for ch in '*?['):
            missing.append(pat)
        else:
            out.append(pat)
    if missing:
        raise SystemExit('ERROR: unmatched glob(s):\n  ' + '\n  '.join(missing))
    return sorted(set(out))

def parse_sigma_text(s):
    if not isinstance(s, str):
        return None
    nums = re.findall(r'[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?', s)
    if len(nums) >= 2:
        return D(nums[0]), D(nums[1])
    return None

def check_coverage(o, p, sigma_min, sigma_max, allow_missing=False):
    # Preferred explicit fields.
    sm = o.get('sigma_min')
    sx = o.get('sigma_max')
    if sm is not None and sx is not None:
        smD, sxD = D(sm), D(sx)
        if smD <= sigma_min and sxD >= sigma_max:
            return {'sigma_min': str(smD), 'sigma_max': str(sxD), 'coverage_check': 'explicit_fields_ok'}
        raise SystemExit(f'ERROR: {p}: explicit sigma coverage [{smD},{sxD}] does not cover requested [{sigma_min},{sigma_max}]')

    # Common domain object.
    dom = o.get('domain')
    if isinstance(dom, dict) and 'sigma' in dom:
        parsed = parse_sigma_text(dom['sigma'])
        if parsed:
            smD, sxD = parsed
            if smD <= sigma_min and sxD >= sigma_max:
                return {'sigma_min': str(smD), 'sigma_max': str(sxD), 'coverage_check': 'domain_sigma_ok'}
            raise SystemExit(f'ERROR: {p}: domain sigma [{smD},{sxD}] does not cover requested [{sigma_min},{sigma_max}]')

    if allow_missing:
        return {'sigma_min': 'MISSING', 'sigma_max': 'MISSING', 'coverage_check': 'MISSING_ALLOWED_FOR_DIAGNOSTIC_ONLY'}
    raise SystemExit(f'ERROR: {p}: missing endpoint sigma coverage fields; rerun certifier with sigma_min/sigma_max in JSON')

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--inputs', nargs='+', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--height', required=True, help='height T to test, e.g. 1e30')
    ap.add_argument('--d-endpoint', required=True, help='endpoint bracket lower bound')
    ap.add_argument('--sigma-min', default='0')
    ap.add_argument('--sigma-max', default='0.01')
    ap.add_argument('--tau-min', default='10')
    ap.add_argument('--normalization', default='N_Q')
    ap.add_argument('--allow-missing-coverage', action='store_true')
    args = ap.parse_args()

    sigma_min = D(args.sigma_min); sigma_max = D(args.sigma_max)
    T = D(args.height); d_endpoint = D(args.d_endpoint)
    paths = expand(args.inputs)
    if not paths:
        raise SystemExit('ERROR: no input files')

    by = {}; deps = {}
    for p in paths:
        if not Path(p).exists():
            raise SystemExit(f'ERROR: missing input {p}')
        o = load(p)
        if o.get('status') != 'proved':
            raise SystemExit(f"ERROR: {p}: status must be proved, got {o.get('status')!r}")
        sec = req(o, 'sector', p)
        if sec not in SECTORS:
            raise SystemExit(f'ERROR: {p}: invalid sector {sec}')
        if sec in by:
            raise SystemExit(f'ERROR: duplicate sector {sec}')
        side = o.get('side', 'UNKNOWN')
        if side not in SIDES:
            raise SystemExit(f'ERROR: {p}: invalid side {side}')
        c = D(req(o, 'C_value_upper_safe', p))
        if c < 0:
            raise SystemExit(f'ERROR: negative constant in {p}')
        cov = check_coverage(o, p, sigma_min, sigma_max, args.allow_missing_coverage)
        by[sec] = c
        deps[sec] = {
            'file': p,
            'side': side,
            'proof_hash': o.get('proof_hash', hash_obj(o)),
            'coverage_hash': o.get('coverage_hash', 'MISSING'),
            'normalization_hash': o.get('normalization_hash', 'MISSING'),
            'C_value_upper_safe': str(c),
            **cov,
        }

    missing = [s for s in SECTORS if s not in by]
    if missing:
        raise SystemExit('ERROR: missing sector certificates: ' + ', '.join(missing))

    C_endpoint = sum(by.values(), Decimal(0))
    lhs = math.sqrt(float(T))
    rhs = float(C_endpoint / d_endpoint) * math.log(float(T))
    ok = lhs > rhs

    out = {
        'status': 'endpoint_left_full_certificate' if ok else 'endpoint_left_height_failed',
        'sigma_min': str(sigma_min),
        'sigma_max': str(sigma_max),
        'tau_min': str(D(args.tau_min)),
        'normalization': args.normalization,
        'C_endpoint_residual_upper_safe': str(C_endpoint),
        'bracket_gap_lower': str(d_endpoint),
        'height_tested': str(T),
        'height_condition': 'sqrt(T) > (C_endpoint/d_endpoint)*log(T)',
        'height_condition_lhs': str(lhs),
        'height_condition_rhs': str(rhs),
        'height_condition_holds': ok,
        'sectors': deps,
    }
    out['proof_hash'] = hash_obj(out)
    Path(args.output).write_text(json.dumps(out, indent=2, sort_keys=True), encoding='utf-8')
    print(json.dumps(out, indent=2, sort_keys=True))
    if not ok:
        sys.exit(2)

if __name__ == '__main__':
    main()
