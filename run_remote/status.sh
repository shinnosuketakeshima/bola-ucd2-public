#!/usr/bin/env bash
# Report progress without touching anything.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RESULTS="$ROOT/results/codex_reanalysis_v1"
MAPPING="$RESULTS/02_mapping"
PYTHON="$HOME/codex-envs/mapping/bin/python"

for stage in pilot cohort; do
  pidfile="$HOME/codex_${stage}.pid"
  if [ -f "$pidfile" ] && kill -0 "$(cat "$pidfile")" 2>/dev/null; then
    echo "$stage driver: RUNNING (pid $(cat "$pidfile"))"
  else
    echo "$stage driver: not running"
  fi
done

echo
printf '%-16s %-10s %12s %12s\n' sample status alignments 'MAPQ>=20'
for status in "$MAPPING"/*.mapping_status.json; do
  [ -e "$status" ] || continue
  "$PYTHON" - "$status" <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
total = d.get("alignment_count")
if total:
    q20 = d.get("mapq_20_59", 0) + d.get("mapq_60_plus", 0)
    print(f'{d.get("sample_id",""):16s} {d.get("status",""):10s} {total:12,d} {q20:12,d}')
else:
    print(f'{d.get("sample_id",""):16s} {d.get("status",""):10s} {"-":>12s} {"-":>12s}')
PY
done

echo
echo "mapped: $(grep -l '"status": "passed"' "$MAPPING"/*.mapping_status.json 2>/dev/null | wc -l) of 13"

# The gate verdict matters more than the count: a mapping that passed with a
# failed smoke still leaves the cohort closed.
GATE="$RESULTS/00_input_audit/pilot_status.json"
if [ -f "$GATE" ]; then
  echo
  "$PYTHON" - "$GATE" <<'PY'
import json, sys
d = json.load(open(sys.argv[1]))
print(f'pilot gate: {d.get("status")} (accepted={d.get("accepted")})')
for check in d.get("checks", []):
    if check.get("status") != "passed":
        print(f'  {check.get("name"):22s} {check.get("status")}')
if d.get("blocking_reason"):
    print(f'  blocking: {d["blocking_reason"]}')
PY
fi

echo
for log in pilot_driver cohort_driver; do
  [ -f "$RESULTS/logs/$log.log" ] || continue
  echo "--- $log ---"
  tail -3 "$RESULTS/logs/$log.log"
done
