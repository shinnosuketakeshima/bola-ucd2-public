#!/usr/bin/env bash
# Verify every delivered FASTQ by content against the audit manifest, once, and
# record the result durably.
#
# This reads ~1.13 TiB and takes roughly 3-4 hours, so it deliberately does not
# live inside the preflight: the preflight expires in hours and would force a
# full re-hash after every long job. The seal it writes is bound to the bytes
# and to the audit manifest, not to a clock — it stays valid until either
# changes, and any change invalidates it automatically.
#
# Run this once, after the drive is mounted on this host. Every later step
# checks the seal cheaply, without re-reading the delivery.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV="$HOME/codex-envs/mapping"
STATUS="$ROOT/results/codex_reanalysis_v1/preflight_status.json"
MANIFEST="$ROOT/results/codex_reanalysis_v1/00_input_audit/input_manifest.tsv"
SEAL="$ROOT/results/codex_reanalysis_v1/00_input_audit/input_seal.json"

export PATH="$ENV/bin:$PATH"
cd "$ROOT"

if [ ! -f "$STATUS" ]; then
  echo "run run_remote/01_preflight.sh first" >&2
  exit 1
fi

# Take the delivery root from the preflight rather than guessing it, so the
# bytes hashed here are the bytes mapping will read. A fresh host can never
# get a "ready" preflight before this script has written the seal -- the seal
# is itself one of the checks -- so a "failed" verdict is accepted too, as
# long as it already located and confirmed a read-only delivery root.
DELIVER="$(python - "$STATUS" <<'PY'
import json, sys
status = json.load(open(sys.argv[1], encoding="utf-8"))
root = status.get("delivery_root")
if not root:
    sys.exit("preflight has not located a read-only delivery root yet; re-run run_remote/01_preflight.sh")
print(root)
PY
)"

echo "verifying delivered inputs under $DELIVER"
echo "this reads the whole delivery and takes about 3-4 hours; it runs once."
echo

python scripts/codex_reanalysis/input_seal.py \
  --deliver-root "$DELIVER" \
  --input-manifest "$MANIFEST" \
  --output "$SEAL"

echo
echo "seal written: $SEAL"
echo "next: bash run_remote/03_run_pilot.sh"
