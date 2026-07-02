#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import hashlib
from decimal import Decimal, getcontext
from pathlib import Path

getcontext().prec = 80


def hash_obj(obj) -> str:
    return hashlib.sha256(
        json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def D(x) -> Decimal:
    return Decimal(str(x))


def safe(x: str | int | float | Decimal) -> Decimal:
    x = D(x)
    if x == 0:
        return Decimal(0)
    return x * (Decimal(1) + Decimal("1e-12")) + Decimal("100")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Normalize a raw proved inner-strip sector/hook JSON into a .safe.json file."
    )
    ap.add_argument("raw", help="Input raw JSON")
    ap.add_argument("safe", help="Output safe JSON")
    args = ap.parse_args()

    raw_path = Path(args.raw)
    out_path = Path(args.safe)

    obj = json.loads(raw_path.read_text(encoding="utf-8"))

    if obj.get("status") != "proved":
        raise SystemExit(f"ERROR: raw JSON must have status='proved', got {obj.get('status')!r}")

    required = [
        "sector",
        "side",
        "C_value_upper",
        "C_sigma_derivative_upper",
        "coverage_hash",
        "normalization_hash",
    ]
    missing = [k for k in required if k not in obj]
    if missing:
        raise SystemExit("ERROR: raw JSON missing required fields: " + ", ".join(missing))

    out = dict(obj)
    out["C_value_upper_safe"] = str(safe(obj["C_value_upper"]))
    out["C_sigma_derivative_upper_safe"] = str(safe(obj["C_sigma_derivative_upper"]))
    out["safe_normalization_rule"] = "safe=0 if raw=0 else raw*(1+1e-12)+100"
    out["source_raw_file"] = str(raw_path)
    out["proof_hash"] = hash_obj(out)

    out_path.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(out, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
