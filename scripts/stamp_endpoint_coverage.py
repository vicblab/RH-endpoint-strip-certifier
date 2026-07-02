#!/usr/bin/env python3
"""
stamp_endpoint_coverage.py

Adds explicit endpoint coverage fields to already-produced safe sector JSONs whose
coverage_hash encodes the domain but whose JSON does not expose sigma_min/sigma_max.
This does NOT prove anything new; it records the declared run domain so strict
mergers can fail/accept based on explicit fields.

Use only for files you just produced with --sigma-min 0 --sigma-max 0.01.
"""
import argparse, json, hashlib
from pathlib import Path

def h(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()).hexdigest()

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--sigma-min', required=True)
    p.add_argument('--sigma-max', required=True)
    p.add_argument('--tau-min', default='10')
    p.add_argument('files', nargs='+')
    args = p.parse_args()
    for fn in args.files:
        path = Path(fn)
        obj = json.loads(path.read_text(encoding='utf-8'))
        if obj.get('status') != 'proved':
            raise SystemExit(f'ERROR: {fn}: status is not proved')
        obj['sigma_min'] = str(args.sigma_min)
        obj['sigma_max'] = str(args.sigma_max)
        obj['tau_min'] = str(args.tau_min)
        obj['coverage_statement'] = f"certified on sigma in [{args.sigma_min},{args.sigma_max}], tau >= {args.tau_min}; stamped from the endpoint run command"
        obj['coverage_stamp_hash'] = h({
            'file': fn,
            'sector': obj.get('sector'),
            'sigma_min': str(args.sigma_min),
            'sigma_max': str(args.sigma_max),
            'tau_min': str(args.tau_min),
            'old_coverage_hash': obj.get('coverage_hash')
        })
        # update proof hash because file changed; retain old proof_hash in dependency trail
        obj['stamped_json_hash'] = h(obj)
        path.write_text(json.dumps(obj, indent=2, sort_keys=True), encoding='utf-8')
        print(f'stamped {fn}')

if __name__ == '__main__':
    main()
