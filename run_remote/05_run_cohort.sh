#!/usr/bin/env bash
# Start cohort mapping, detached from this shell so an SSH disconnection cannot
# kill a multi-day job.
#
# Each sample is an independent job that publishes its outputs only once every
# artefact is complete and reconciled. A power loss therefore costs at most the
# samples in flight; finished samples stay finished and are skipped on restart.
#
# Optional argument: number of concurrent samples. Defaults to the count the
# preflight recommended for this host.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV="$HOME/codex-envs/mapping"
STATUS="$ROOT/results/codex_reanalysis_v1/preflight_status.json"
GATE="$ROOT/results/codex_reanalysis_v1/00_input_audit/pilot_status.json"
SEAL="$ROOT/results/codex_reanalysis_v1/00_input_audit/input_seal.json"
MANIFEST="$ROOT/results/codex_reanalysis_v1/00_input_audit/input_manifest.tsv"
LOG="$ROOT/results/codex_reanalysis_v1/logs/cohort_driver.log"
PIDFILE="$HOME/codex_cohort.pid"
RUNNER="$HOME/codex_cohort_run.sh"

export PATH="$ENV/bin:$PATH"
cd "$ROOT"

# The presence of a status file proves nothing: it may record a failure, may
# predate a remount, or may describe a machine state that has since changed.
# The verdict, its age, the current mount and the pilot gate are all checked in
# one pass, and the settings used below are the ones that pass came from.
VERIFIED="$(python scripts/codex_reanalysis/remote_gate_check.py \
  --preflight-status "$STATUS" \
  --input-seal "$SEAL" \
  --input-manifest "$MANIFEST" \
  --pilot-status "$GATE" \
  --require-pilot-gate)"

DELIVER="$(printf '%s' "$VERIFIED" | sed -n 1p)"
JOBS="${1:-$(printf '%s' "$VERIFIED" | sed -n 2p)}"
CORES=$((JOBS * 8))

if [ -f "$PIDFILE" ] && kill -0 "$(cat "$PIDFILE")" 2>/dev/null; then
  echo "cohort already running with pid $(cat "$PIDFILE")"
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
  --cores $CORES --rerun-incomplete mapping_all
echo "COHORT_DRIVER_EXIT=\$?"
RUNNER_EOF
chmod +x "$RUNNER"

setsid nohup bash "$RUNNER" > "$LOG" 2>&1 < /dev/null &
echo $! > "$PIDFILE"
sleep 3

echo "cohort started"
echo "  pid            : $(cat "$PIDFILE")"
echo "  concurrent     : $JOBS samples (--cores $CORES)"
echo "  delivery root  : $DELIVER"
echo "  driver log     : $LOG"
echo
echo "progress: bash run_remote/status.sh"
