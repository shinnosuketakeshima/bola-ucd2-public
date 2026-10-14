#!/usr/bin/env bash
# Map the two pilot animals, detached, before any cohort work is considered.
#
# The pilot is not a formality: it is the only place where the widened
# retention intervals, the -N ceiling and the assembly behaviour are tested on
# this host, under these settings, against this reference. The cohort refuses to
# start without a pilot that passed all of it.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV="$HOME/codex-envs/mapping"
STATUS="$ROOT/results/codex_reanalysis_v1/preflight_status.json"
SEAL="$ROOT/results/codex_reanalysis_v1/00_input_audit/input_seal.json"
MANIFEST="$ROOT/results/codex_reanalysis_v1/00_input_audit/input_manifest.tsv"
LOG="$ROOT/results/codex_reanalysis_v1/logs/pilot_driver.log"
PIDFILE="$HOME/codex_pilot.pid"
RUNNER="$HOME/codex_pilot_run.sh"

export PATH="$ENV/bin:$PATH"
cd "$ROOT"

# The host, the mount and the input seal are checked in one pass, and the
# settings used below are the ones that pass came from. No pilot gate is
# required here — this is the run that produces it.
VERIFIED="$(python scripts/codex_reanalysis/remote_gate_check.py \
  --preflight-status "$STATUS" \
  --input-seal "$SEAL" \
  --input-manifest "$MANIFEST")"

DELIVER="$(printf '%s' "$VERIFIED" | sed -n 1p)"
JOBS="${1:-2}"
CORES=$((JOBS * 8))

if [ -f "$PIDFILE" ] && kill -0 "$(cat "$PIDFILE")" 2>/dev/null; then
  echo "pilot already running with pid $(cat "$PIDFILE")"
  exit 0
fi

mkdir -p "$(dirname "$LOG")"

cat > "$RUNNER" <<RUNNER_EOF
#!/usr/bin/env bash
set -euo pipefail
export PATH="$ENV/bin:\$PATH"
export CODEX_DELIVER_ROOT="$DELIVER"
# Validated above, so provenance can record the delivered hash as verified on
# this host rather than carried over from the audit.
export CODEX_INPUT_SEAL="$SEAL"
cd "$ROOT"
snakemake -s workflows/codex_reanalysis/Snakefile \
  --configfile configs/codex_reanalysis.example.yaml \
  --cores $CORES --rerun-incomplete mapping_pilot
echo "PILOT_DRIVER_EXIT=\$?"
RUNNER_EOF
chmod +x "$RUNNER"

setsid nohup bash "$RUNNER" > "$LOG" 2>&1 < /dev/null &
echo $! > "$PIDFILE"
sleep 3

echo "pilot started"
echo "  pid            : $(cat "$PIDFILE")"
echo "  concurrent     : $JOBS samples (--cores $CORES)"
echo "  delivery root  : $DELIVER"
echo "  driver log     : $LOG"
echo
echo "progress: bash run_remote/status.sh"
echo "then    : bash run_remote/04_pilot_gate.sh"
