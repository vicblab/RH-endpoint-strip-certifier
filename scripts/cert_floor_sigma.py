#!/usr/bin/env python3
from __future__ import annotations
import json, hashlib

def h(o):
    return hashlib.sha256(json.dumps(o, sort_keys=True, separators=(",", ":")).encode()).hexdigest()

obj = {
    "status": "proved",
    "sector": "floor",
    "side": "zero",
    "method": "half_open_floor_sector_identically_zero",
    "C_value_upper": "0",
    "C_sigma_derivative_upper": "0",
    "coverage_hash": h({"sector":"floor","coverage":"zero by half-open convention"}),
    "normalization_hash": h({"sector":"floor","side":"zero","projection":"none"}),
    "mathematical_statement": "R_floor(s)=0 identically; d/dsigma R_floor(s)=0."
}
obj["proof_hash"] = h(obj)
print(json.dumps(obj, indent=2, sort_keys=True))
