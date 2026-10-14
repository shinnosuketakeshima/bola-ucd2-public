#!/usr/bin/env bash
# Create the pinned analysis environment on this host. Needs network access the
# first time. Installs entirely under $HOME, so no root is required.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV="$HOME/codex-envs/mapping"
export MAMBA_ROOT_PREFIX="$HOME/.codex-mamba"
LOCK="$ROOT/workflows/codex_reanalysis/envs/mapping.linux-64.lock"

mkdir -p "$MAMBA_ROOT_PREFIX" "$HOME/codex-envs" "$HOME/.local/bin"

if ! command -v micromamba >/dev/null 2>&1 && [ ! -x "$HOME/.local/bin/micromamba" ]; then
  echo "=== downloading micromamba ==="
  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
    | tar -xvj -C "$HOME/.local" bin/micromamba
fi
MM="$(command -v micromamba || echo "$HOME/.local/bin/micromamba")"

if [ -x "$ENV/bin/minimap2" ]; then
  echo "environment already present at $ENV"
else
  echo "=== creating environment from the explicit lock ==="
  # The lock pins exact builds, so this host gets the same packages the pilot
  # ran on rather than a fresh solve that could differ.
  "$MM" create -y -p "$ENV" -f "$LOCK"
fi

# Flye is only used by the pilot smoke assembly, and it is deliberately not in
# the mapping lock: adding it would change the environment the mapping results
# are produced in. It lives in its own prefix, pinned by version, and the
# version it resolves to is recorded in the smoke status.
ASSEMBLY="$HOME/codex-envs/assembly"
if [ -x "$ASSEMBLY/bin/flye" ]; then
  echo "assembly environment already present at $ASSEMBLY"
else
  echo "=== creating the smoke assembly environment ==="
  "$MM" create -y -p "$ASSEMBLY" -c conda-forge -c bioconda "flye=2.9.6"
fi

echo "=== versions ==="
"$ENV/bin/minimap2" --version
"$ENV/bin/samtools" --version | head -1
"$ENV/bin/python" --version
"$ENV/bin/pigz" --version
"$ENV/bin/git" --version
"$ENV/bin/snakemake" --version
"$ASSEMBLY/bin/flye" --version
echo
echo "environment ready: $ENV"
echo "smoke assembler : $ASSEMBLY"
echo "next: bash run_remote/01_preflight.sh"
