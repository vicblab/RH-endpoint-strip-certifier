# Endpoint rerun skeleton after sigma parameterization.
# Run from repository root after inspecting the patch/diff.

mkdir -p certs_endpoint

# 1) Compact grouped boundary, finite tau slabs
python scripts/bd_arb_compact_boundary.py grp-compact --sigma-min 0 --sigma-max 0.01 --tau-min 10 --tau-max 25 --Y 128 --ns 4 --nt 8 --ny 256 > certs_endpoint/compact_grp_10_25.raw.json
python scripts/bd_arb_compact_boundary.py grp-compact --sigma-min 0 --sigma-max 0.01 --tau-min 25 --tau-max 100 --Y 128 --ns 4 --nt 8 --ny 256 > certs_endpoint/compact_grp_25_100.raw.json
python scripts/bd_arb_compact_boundary.py grp-compact --sigma-min 0 --sigma-max 0.01 --tau-min 100 --tau-max 400 --Y 128 --ns 4 --nt 8 --ny 256 > certs_endpoint/compact_grp_100_400.raw.json
python scripts/bd_arb_compact_boundary.py grp-compact --sigma-min 0 --sigma-max 0.01 --tau-min 400 --tau-max 1600 --Y 128 --ns 4 --nt 8 --ny 256 > certs_endpoint/compact_grp_400_1600.raw.json
python scripts/bd_arb_compact_boundary.py grp-compact --sigma-min 0 --sigma-max 0.01 --tau-min 1600 --tau-max 6400 --Y 128 --ns 4 --nt 8 --ny 256 > certs_endpoint/compact_grp_1600_6400.raw.json
python scripts/bd_arb_compact_boundary.py grp-compact --sigma-min 0 --sigma-max 0.01 --tau-min 6400 --tau-max 65536 --Y 128 --ns 4 --nt 8 --ny 256 > certs_endpoint/compact_grp_6400_65536.raw.json

# 2) High-tau Taylor grouped compact-boundary slab
python scripts/bd_arb_compact_boundary.py central-taylor-compact --sigma-min 0 --sigma-max 0.01 --w0 0.00390625 --Y 128 --ns 4 --nw 4 --nr 48 --ny 256 --nl 8 > certs_endpoint/compact_grp_taylor_w0.raw.json

# 3) Compact high-tail slabs
python scripts/bd_arb_compact_boundary.py tail-compact --sigma-min 0 --sigma-max 0.01 --w-min 0.00390625 --w-max 0.00625 --R 32 --Y 128 --ns 4 --nw 8 --nr 48 --ny 256 > certs_endpoint/compact_tail_w_0p00390625_0p00625.raw.json
python scripts/bd_arb_compact_boundary.py tail-compact --sigma-min 0 --sigma-max 0.01 --w-min 0.00625 --w-max 0.0125 --R 32 --Y 128 --ns 4 --nw 8 --nr 48 --ny 256 > certs_endpoint/compact_tail_w_0p00625_0p0125.raw.json
python scripts/bd_arb_compact_boundary.py tail-compact --sigma-min 0 --sigma-max 0.01 --w-min 0.0125 --w-max 0.025 --R 32 --Y 128 --ns 4 --nw 8 --nr 48 --ny 256 > certs_endpoint/compact_tail_w_0p0125_0p025.raw.json
python scripts/bd_arb_compact_boundary.py tail-compact --sigma-min 0 --sigma-max 0.01 --w-min 0.025 --w-max 0.05 --R 32 --Y 128 --ns 4 --nw 8 --nr 48 --ny 256 > certs_endpoint/compact_tail_w_0p025_0p05.raw.json
python scripts/bd_arb_compact_boundary.py tail-compact --sigma-min 0 --sigma-max 0.01 --w-min 0.05 --w-max 0.31622776601683794 --R 32 --Y 128 --ns 4 --nw 8 --nr 48 --ny 256 > certs_endpoint/compact_tail_w_0p05_0p316.raw.json
python scripts/bd_arb_compact_boundary.py tail-taylor-compact --sigma-min 0 --sigma-max 0.01 --w0 0.00390625 --R 32 --Y 128 --ns 4 --nw 4 --nr 48 --ny 256 --nl 8 > certs_endpoint/compact_tail_taylor_w0.raw.json

# 4) Compact-complement tail envelope/hook
python scripts/cert_chic_tail_envelope.py prove --sigma-min 0 --sigma-max 0.01 --Y 128 --mu-max 0.31622776601683794 --mu0 0.00048828125 --ns 4 --nmu 32 --ny 256 > certs_endpoint/proved_tail_derivative_envelope.raw.json
python scripts/cert_chic_tail.py prove --envelope-json certs_endpoint/proved_tail_derivative_envelope.raw.json --Y 128 --R 32 > certs_endpoint/cert_chic_tail.raw.json

# 5) Compact-complement nonstationary/stationary hooks
python scripts/cert_chic_nonstationary.py transition --sigma-min 0 --sigma-max 0.01 --Y 128 --r-min -4 --r-max 32 --w-max 0.31622776601683794 --lambda0 0.25 --ns 4 --nw 16 --nr 96 --ny 128 > certs_endpoint/cert_chic_nonstationary.raw.json
python scripts/cert_chic_stationary.py prove --sigma-min 0 --sigma-max 0.01 --Y 128 --r-min -4 --r-max 32 --w-max 0.31622776601683794 --lambda0 0.25 --ns 4 --nw 16 --nr 96 --ny 128 > certs_endpoint/cert_chic_stationary.raw.json

# NOTE: You still need normalization/merge scripts that convert fragment outputs to proved sector JSONs.
# The dumped normalize_certificate_json.py is empty, so normalization must be supplied or restored.
