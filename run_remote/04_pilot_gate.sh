#!/usr/bin/env bash
# Smoke-test the pilot mappings, then judge them.
#
# The smoke stage answers what mapping counts cannot: whether the -N ceiling was
# still truncating the secondary list, whether split reads keep every segment
# inside the widened retention intervals, and whether a local assembly of the
# contiguous MHC stretch reconstructs the region at the expected span, identity
# and gene content. The gate then folds those into the verdict the cohort reads.
#
# Nothing here is optional: a check with no evidence is reported not_run, which
# keeps the verdict out of 'passed' and keeps the cohort closed.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV="$HOME/codex-envs/mapping"
RESULTS="$ROOT/results/codex_reanalysis_v1"
MAPPING="$RESULTS/02_mapping"
REFERENCE="$RESULTS/01_reference"
SMOKE="$MAPPING/pilot_smoke"

ASSEMBLY="$HOME/codex-envs/assembly"
FLYE="$ASSEMBLY/bin/flye"

# Flye is invoked below by its full path, but it internally shells out to its
# own bundled `flye-minimap2` by name (`which("flye-minimap2")`), so its bin
# directory has to be on PATH too, not just findable via $FLYE.
export PATH="$ENV/bin:$ASSEMBLY/bin:$PATH"
cd "$ROOT"

if [ ! -x "$FLYE" ]; then
  echo "the smoke assembler is missing; run run_remote/00_setup_env.sh" >&2
  exit 1
fi

SECONDARY_MAX="$(python - <<'PY'
from scripts.codex_reanalysis.stream_target_mapping import DEFAULT_SECONDARY_MAX
print(DEFAULT_SECONDARY_MAX)
PY
)"

# The assembled stretch is the longest retained interval, which is the
# contiguous MHC span by construction. Deriving it keeps the smoke aligned with
# the BED the mapping actually used instead of a coordinate copied by hand.
REGION="$(python - "$REFERENCE/mapping_targets.bed" <<'PY'
import sys
best = None
for line in open(sys.argv[1], encoding="utf-8"):
    if not line.strip():
        continue
    contig, start, end = line.split("\t")[:3]
    span = int(end) - int(start)
    if best is None or span > best[0]:
        best = (span, contig, int(start), int(end))
if best is None:
    sys.exit("mapping BED contains no intervals")
_, contig, start, end = best
print(f"{contig}:{start + 1}-{end}")
PY
)"

SAMPLES="$(ls "$MAPPING"/*.mapping_status.json 2>/dev/null | xargs -r -n1 basename \
  | sed 's/\.mapping_status\.json$//')"
if [ -z "$SAMPLES" ]; then
  echo "no pilot mappings found; run run_remote/03_run_pilot.sh first" >&2
  exit 1
fi

mkdir -p "$SMOKE"
echo "smoke region  : $REGION"
echo "-N ceiling    : $SECONDARY_MAX"
echo "pilot samples : $(echo "$SAMPLES" | tr '\n' ' ')"
echo

for sample in $SAMPLES; do
  echo "=== $sample ==="
  # A non-zero exit means the smoke failed, not that the step is broken. Record
  # it and carry on, so the gate reports every animal rather than the first.
  python scripts/codex_reanalysis/pilot_smoke.py \
    --sample-id "$sample" \
    --bam "$MAPPING/$sample.targets.bam" \
    --mapping-bed "$REFERENCE/mapping_targets.bed" \
    --analytical-bed "$REFERENCE/targets.bed" \
    --reference "$REFERENCE/GCF_002263795.3_ARS-UCD2.0_genomic.fna" \
    --region "$REGION" \
    --secondary-max "$SECONDARY_MAX" \
    --output "$SMOKE/$sample.smoke_status.json" \
    --work-dir "$SMOKE/$sample" \
    --flye "$FLYE" \
    --threads 8 || echo "  ($sample smoke did not pass; recorded)"
  echo
done

echo "=== pilot gate ==="
set +e
python scripts/codex_reanalysis/pilot_gate.py \
  --root "$RESULTS" \
  --expected 2 \
  --report "$ROOT/docs/codex_reanalysis/PILOT_REPORT.md" \
  --config configs/codex_reanalysis.example.yaml \
  --samtools samtools
GATE_EXIT=$?
set -e

echo
if [ "$GATE_EXIT" -eq 0 ]; then
  echo "the pilot was accepted."
  echo "next: re-run bash run_remote/01_preflight.sh, then bash run_remote/05_run_cohort.sh"
else
  echo "the pilot was NOT accepted. Read:"
  echo "  $RESULTS/00_input_audit/pilot_status.json"
  echo "  $ROOT/docs/codex_reanalysis/PILOT_REPORT.md"
  echo "The cohort will refuse to start until every check passes."
fi
exit "$GATE_EXIT"
