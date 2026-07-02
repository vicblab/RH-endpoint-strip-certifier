# RH Endpoint Strip Certifier

This repository reproduces my endpoint left certificate for the Mellin symmetric quotient obstruction project.

The purpose of this repository is precise:

```text
Certify the endpoint left residual bound on 0 <= sigma <= 0.01,
then verify the endpoint height inequality at T = 1e30.
```

The successful reference run produced:

```text
status = endpoint_left_full_certificate
height_condition_holds = True
proof_hash = 49a8a997898443ee15de999d3ace640854d8bddec2c20b408e0b2eda0cb3da20
```

This repository is a computational certificate repository. It is not by itself a proof of the Riemann Hypothesis. It supplies the endpoint left certificate that must be combined with the structural quotient obstruction theorem, the outer strip certificate, the inner strip certificate, and the reflection argument.

---

# 1. What this repository certifies

The endpoint left certifier works on the closed interval

```text
0 <= sigma <= 0.01
```

and base height

```text
tau >= 10
```

The final theorem application uses the open endpoint strip

```text
0 < Re(s) < 0.01
```

at high height

```text
absolute value of Im(s) >= 1e30
```

The closed interval

```text
0 <= sigma <= 0.01
```

is used because rigorous interval arithmetic certifies compact boxes.

The quotient normalization used here is

```text
N_Q(s) = B(s) + ((s - 1) / s) A(s)
```

The endpoint bracket lower bound used in the final height comparison is

```text
d_endpoint = 0.20405207717382234
```

The final height condition checked by the certificate is

```text
sqrt(T) > (C_endpoint / d_endpoint) log(T)
```

where `log` is the natural logarithm.

---

# 2. Reference successful output

The successful reference run produced the final file

```text
ENDPOINT_LEFT_CERTIFICATE_SAFE.json
```

with the decisive fields

```text
status = endpoint_left_full_certificate
sigma_min = 0
sigma_max = 0.01
tau_min = 10
normalization = N_Q

C_endpoint_residual_upper_safe =
92517448456.63735671709576210785011758253131058311532566421240408489986855178714577717052485812881457

bracket_gap_lower = 0.20405207717382234
height_tested = 1E+30
height_condition = sqrt(T) > (C_endpoint/d_endpoint) log(T)
height_condition_holds = True
height_condition_lhs = 1000000000000000.0
height_condition_rhs = 31319842553226.53

proof_hash =
49a8a997898443ee15de999d3ace640854d8bddec2c20b408e0b2eda0cb3da20
```

The key numerical comparison is

```text
1000000000000000.0 > 31319842553226.53
```

So the endpoint left certificate passes the height test at

```text
T = 1e30
```

---

# 3. What this repository does not claim

This repository does not by itself prove the Riemann Hypothesis.

It certifies only the endpoint left numerical interface. To obtain the intended high height off critical conclusion in the full critical strip, this endpoint result must be assembled with:

1. The structural quotient obstruction theorem.
2. The already certified outer strip result.
3. The already certified inner strip two channel result.
4. Reflection from the left endpoint to the right endpoint.
5. A finite height verification below the final analytic height threshold.

The endpoint result supplied here is for

```text
0 < Re(s) < 0.01, absolute value of Im(s) >= 1e30
```

and after reflection for

```text
0.99 < Re(s) < 1, absolute value of Im(s) >= 1e30
```

Critical line zeros are not excluded here, and the Riemann Hypothesis allows critical line zeros.

---

# 4. Repository structure

The repository should have this structure:

```text
RH-endpoint-strip-certifier/
├── README.md
├── requirements.txt
├── run_endpoint_certificate.sh
├── ENDPOINT_LEFT_CERTIFICATE_SAFE.json
├── SHA256SUMS.txt
├── scripts/
│   ├── bd_arb_compact_boundary.py
│   ├── cert_chic_tail_envelope.py
│   ├── cert_chic_tail.py
│   ├── cert_chic_nonstationary.py
│   ├── cert_chic_stationary.py
│   ├── cert_band_sigma.py
│   ├── cert_core_sigma.py
│   ├── cert_endpoint_sigma.py
│   ├── cert_far_sigma.py
│   ├── cert_floor_sigma.py
│   ├── cert_nonstat_sigma.py
│   ├── normalize_inner_certificate_json.py
│   ├── stamp_endpoint_coverage.py
│   ├── merge_endpoint_boundary_sector.py
│   └── merge_endpoint_certificate_safe.py
├── certs_endpoint/
│   ├── band_sector.safe.json
│   ├── bd_sector.safe.json
│   ├── far_sector.safe.json
│   ├── end_sector.safe.json
│   ├── nonstat_sector.safe.json
│   ├── floor_sector.safe.json
│   ├── core_sector.safe.json
│   └── raw certificate and hook JSON files
├── expected/
│   └── EXPECTED_FINAL_FIELDS.json
└── docs/
    └── ENDPOINT_PAPER_OUTLINE.md
```

The most important final file is

```text
ENDPOINT_LEFT_CERTIFICATE_SAFE.json
```

---

# 5. Install dependencies

Use a clean Linux Python virtual environment.

From the repository root, run:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

The `requirements.txt` file should contain

```text
python-flint
```

This package provides Arb and FLINT interval arithmetic through Python.

---

# 6. Quick check of the included reference certificate

Before recomputing anything, first check that the included reference certificate is valid JSON.

Run:

```bash
python -m json.tool ENDPOINT_LEFT_CERTIFICATE_SAFE.json >/dev/null && echo "final JSON OK"
```

Expected output:

```text
final JSON OK
```

Now inspect the decisive fields:

```bash
python - <<'PY'
import json
from pathlib import Path

obj = json.loads(Path("ENDPOINT_LEFT_CERTIFICATE_SAFE.json").read_text())

keys = [
    "status",
    "sigma_min",
    "sigma_max",
    "tau_min",
    "normalization",
    "C_endpoint_residual_upper_safe",
    "bracket_gap_lower",
    "height_tested",
    "height_condition",
    "height_condition_holds",
    "height_condition_lhs",
    "height_condition_rhs",
    "proof_hash"
]

for k in keys:
    print(k, "=", obj.get(k))
PY
```

Expected output should include:

```text
status = endpoint_left_full_certificate
sigma_min = 0
sigma_max = 0.01
tau_min = 10
normalization = N_Q
height_condition_holds = True
proof_hash = 49a8a997898443ee15de999d3ace640854d8bddec2c20b408e0b2eda0cb3da20
```

---

# 7. Full recomputation from scratch

To reproduce the certificate from scratch, run:

```bash
./run_endpoint_certificate.sh
```

This script rebuilds the `certs_endpoint/` calculation outputs, regenerates all sector certificates, merges them, and produces:

```text
ENDPOINT_LEFT_CERTIFICATE_SAFE.json
```

When the script succeeds, it should finish by printing something like:

```text
Endpoint certificate verified: 49a8a997898443ee15de999d3ace640854d8bddec2c20b408e0b2eda0cb3da20
```

If the script stops with an error, the certificate has not been reproduced.

---

# 8. What the full recomputation does

The reproduction script certifies seven residual sectors:

```text
band
bd
far
end
nonstat
floor
core
```

The non boundary sectors are generated by:

```text
cert_band_sigma.py
cert_core_sigma.py
cert_endpoint_sigma.py
cert_far_sigma.py
cert_floor_sigma.py
cert_nonstat_sigma.py
```

The boundary sector is generated from compact boundary fragments and compact complement hooks using:

```text
bd_arb_compact_boundary.py
cert_chic_tail_envelope.py
cert_chic_tail.py
cert_chic_nonstationary.py
cert_chic_stationary.py
merge_endpoint_boundary_sector.py
```

The final endpoint certificate is merged by:

```text
merge_endpoint_certificate_safe.py
```

The merge succeeds only if all required sectors have:

```text
status = proved
```

and the final endpoint height inequality passes.

---

# 9. Expected sector constants

In the successful reference run, the final safe sector constants were approximately:

```text
band     = 5.100100e6
bd       = 2.2441347956e10
far      = 5.0000000100e10
end      = 7.1000100000e7
nonstat  = 2.0000000100e10
floor    = 0
core     = 1.001507e2
```

The final merged endpoint residual constant was:

```text
C_endpoint_residual_upper_safe =
92517448456.63735671709576210785011758253131058311532566421240408489986855178714577717052485812881457
```

The dominant sectors are:

```text
far
bd
nonstat
```

This is useful for later constant refinement.

---

# 10. Check all sector files

After running the full reproduction script, verify all required sector files exist and are valid JSON:

```bash
python -m json.tool certs_endpoint/band_sector.safe.json >/dev/null && echo "band OK"
python -m json.tool certs_endpoint/bd_sector.safe.json >/dev/null && echo "bd OK"
python -m json.tool certs_endpoint/far_sector.safe.json >/dev/null && echo "far OK"
python -m json.tool certs_endpoint/end_sector.safe.json >/dev/null && echo "end OK"
python -m json.tool certs_endpoint/nonstat_sector.safe.json >/dev/null && echo "nonstat OK"
python -m json.tool certs_endpoint/floor_sector.safe.json >/dev/null && echo "floor OK"
python -m json.tool certs_endpoint/core_sector.safe.json >/dev/null && echo "core OK"
python -m json.tool ENDPOINT_LEFT_CERTIFICATE_SAFE.json >/dev/null && echo "endpoint final OK"
```

Expected output:

```text
band OK
bd OK
far OK
end OK
nonstat OK
floor OK
core OK
endpoint final OK
```

---

# 11. Check the final height comparison directly

Run:

```bash
python - <<'PY'
import json
from pathlib import Path

obj = json.loads(Path("ENDPOINT_LEFT_CERTIFICATE_SAFE.json").read_text())

lhs = float(obj["height_condition_lhs"])
rhs = float(obj["height_condition_rhs"])

print("lhs =", lhs)
print("rhs =", rhs)
print("lhs > rhs =", lhs > rhs)
PY
```

Expected output:

```text
lhs = 1000000000000000.0
rhs = 31319842553226.53
lhs > rhs = True
```

---

# 12. List final sector constants

Run:

```bash
python - <<'PY'
import json
from pathlib import Path

obj = json.loads(Path("ENDPOINT_LEFT_CERTIFICATE_SAFE.json").read_text())

print("Final endpoint status:", obj["status"])
print("Final endpoint C:", obj["C_endpoint_residual_upper_safe"])
print()

for name, sec in obj["sectors"].items():
    print(name)
    print("  C_value_upper_safe =", sec["C_value_upper_safe"])
    print("  side =", sec["side"])
    print("  sigma_min =", sec["sigma_min"])
    print("  sigma_max =", sec["sigma_max"])
    print("  proof_hash =", sec["proof_hash"])
    print()
PY
```

Expected result:

```text
Final endpoint status: endpoint_left_full_certificate
```

Each sector should report:

```text
sigma_min = 0
sigma_max = 0.01
```

---

# 13. Hash verification

The repository includes:

```text
SHA256SUMS.txt
```

To verify hashes, run:

```bash
sha256sum -c SHA256SUMS.txt
```

Expected result:

```text
Every checked file should report OK.
```

If you regenerate certificates, hashes may change. In that case, regenerate the hash file with:

```bash
find . -type f -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS.txt
```

Then commit the updated `SHA256SUMS.txt`.

---

# 14. How to make a release ZIP

From the parent directory of the repository, run:

```bash
zip -r RH-endpoint-strip-certifier.zip RH-endpoint-strip-certifier
```

The ZIP file can then be uploaded to GitHub or archived.

---

# 15. How this fits into the endpoint paper

The endpoint paper should cite:

```text
ENDPOINT_LEFT_CERTIFICATE_SAFE.json
status = endpoint_left_full_certificate
proof_hash = 49a8a997898443ee15de999d3ace640854d8bddec2c20b408e0b2eda0cb3da20
```

The paper should record:

```text
C_endpoint_residual_upper_safe =
92517448456.63735671709576210785011758253131058311532566421240408489986855178714577717052485812881457

d_endpoint =
0.20405207717382234
```

and the height comparison:

```text
1000000000000000.0 > 31319842553226.53
```

A clean endpoint theorem statement should be made at:

```text
absolute value of Im(s) >= 1e30
```

---

# 16. Publication caution

This repository is the computational certificate layer.

A rigorous endpoint paper must also provide the mathematical theorem layer:

1. Define the quotient numerator.
2. State the zero condition.
3. State the endpoint bracket.
4. Prove the endpoint bracket lower bound.
5. Define the residual decomposition.
6. Define the certificate contract.
7. Cite this repository and final JSON hash.
8. Prove the height contradiction.
9. Reflect the left endpoint result to the right endpoint.
10. Assemble with the already certified outer and inner strips.

