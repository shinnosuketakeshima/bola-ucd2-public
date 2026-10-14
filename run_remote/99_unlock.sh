#!/usr/bin/env bash
# Clear the Snakemake directory lock left behind by a power loss or a kill.
# Safe to run when nothing is executing; refuses while the driver is alive.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV="$HOME/codex-envs/mapping"
export PATH="$ENV/bin:$PATH"

# Either driver can hold the lock, so both have to be dead before it is cleared.
for stage in pilot cohort; do
  pidfile="$HOME/codex_${stage}.pid"
  if [ -f "$pidfile" ] && kill -0 "$(cat "$pidfile")" 2>/dev/null; then
    echo "the $stage driver is still running (pid $(cat "$pidfile")); not unlocking" >&2
    exit 1
  fi
done

cd "$ROOT"
snakemake -s workflows/codex_reanalysis/Snakefile \
  --configfile configs/codex_reanalysis.example.yaml --unlock
rm -f "$HOME/codex_pilot.pid" "$HOME/codex_cohort.pid"

echo
echo "unlocked. Finished samples are kept; restart with:"
echo "  bash run_remote/03_run_pilot.sh    (before the gate)"
echo "  bash run_remote/05_run_cohort.sh   (after the gate passed)"
