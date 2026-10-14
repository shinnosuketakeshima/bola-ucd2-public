#!/usr/bin/env bash
# Verify this host can run the cohort. Nothing long starts until this passes.
# Optional argument: an explicit delivery root, if auto-detection cannot find it.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV="$HOME/codex-envs/mapping"
export PATH="$ENV/bin:$PATH"

if [ ! -x "$ENV/bin/python" ]; then
  echo "environment is missing; run run_remote/00_setup_env.sh first" >&2
  exit 1
fi

cd "$ROOT"
if [ $# -ge 1 ]; then
  python scripts/codex_reanalysis/remote_preflight.py --project-root "$ROOT" --deliver-root "$1"
else
  python scripts/codex_reanalysis/remote_preflight.py --project-root "$ROOT"
fi

echo
if [ -f "$ROOT/results/codex_reanalysis_v1/00_input_audit/input_seal.json" ]; then
  # The seal is content-bound, so it survives this preflight expiring. Once it
  # exists, the preflight only refreshes the machine state around it.
  echo "next: bash run_remote/03_run_pilot.sh   (or 05_run_cohort.sh once the pilot passed)"
else
  echo "next: bash run_remote/02_verify_inputs.sh"
fi
