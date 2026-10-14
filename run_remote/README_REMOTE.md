# Running the cohort mapping on a remote host

The delivery drive is exFAT and irreplaceable, so it is mounted **read-only**
and nothing in this procedure writes to it. The project tree, all outputs and
all temporary files live on the remote host's local disk.

The steps are numbered and each one refuses to start until the previous one
passed. That ordering is enforced in code, not by convention: the cohort script
reads the pilot verdict and stops if it is anything other than a clean pass.

## What the host needs

- Linux, Python-capable user account, **no root required** for the analysis
- ~40 GB free on the local disk (project copy, outputs, sort temporaries and the
  pilot smoke assembly), or ~300 GB if the normalised inputs are copied across
  instead of linked — see step 1
- Network access once, to build the environment from the pinned lock file
- Root **only** to mount the drive, if the desktop does not auto-mount it

## Step 0 — mount the drive read-only (the one step that may need sudo)

If the drive auto-mounts under `/media/<user>/…` when plugged in while logged
into the desktop, nothing is needed. Otherwise:

```bash
lsblk -f                     # find the exFAT partition, e.g. /dev/sdb1
sudo mkdir -p /mnt/mybook
sudo mount -t exfat -o ro,uid=$(id -u),gid=$(id -g) /dev/sdb1 /mnt/mybook
```

`ro` is deliberate: a bug, a crash or a power loss then cannot damage the
delivered sequencing data.

## Step 1 — copy the project to local disk

From the mounted drive:

```bash
mkdir -p ~/bola_ucd2
rsync -a --info=progress2 \
  --exclude 'results/codex_reanalysis_v1/02_mapping' \
  --exclude 'results/codex_reanalysis_v1/00b_normalized' \
  /mnt/mybook/class1_analysis/ ~/bola_ucd2/

# The two normalised inputs are ~240 GB. Link them instead of copying: the
# register resolves its entries relative to the manifest after following
# symlinks, so a link to the read-only drive works and the local disk only
# needs ~40 GB.
ln -s /mnt/mybook/class1_analysis/results/codex_reanalysis_v1/00b_normalized \
      ~/bola_ucd2/results/codex_reanalysis_v1/00b_normalized
```

If the local disk has ~300 GB spare, dropping the second `--exclude` and
copying them instead is equally correct and slightly faster to read.

`.git` must come across — provenance records the commit SHA of every run and a
missing repository degrades that record to `unknown`.

## Step 2 — build the environment

```bash
bash ~/bola_ucd2/run_remote/00_setup_env.sh
```

Installs micromamba and the mapping environment under `$HOME` from
`workflows/codex_reanalysis/envs/mapping.linux-64.lock`, which pins exact builds
so this host gets the packages the pilot ran on rather than a fresh solve that
could differ.

It also creates a **separate** prefix holding Flye, used only by the pilot smoke
assembly. Flye is kept out of the mapping lock on purpose: adding it would change
the environment the mapping results are produced in. Its resolved version is
recorded in each smoke status.

## Step 3 — preflight

```bash
bash ~/bola_ucd2/run_remote/01_preflight.sh
# or, if the drive is somewhere unusual:
bash ~/bola_ucd2/run_remote/01_preflight.sh /mnt/mybook/Deliver_ADD
```

Finds the delivery root, confirms the mount is read-only, re-verifies the audit
and reference stage gates against their checksums, checks the tool versions
against the pinned configuration, and sizes the concurrency from available
memory and CPUs.

**The reference seal will usually report `resealed` here, and that is expected.**
The seal's cheap staleness guard includes modification time, and exFAT stores
timestamps as local time, so moving the drive between machines changes them
while the bytes stay identical. Re-sealing can only succeed if every byte still
matches the SHA-256 the reference stage recorded; changed content is a hard
stop, not a re-seal.

This verdict expires after 6 hours, because it describes machine state that can
change. Re-run it before each long job — it takes a few minutes.

## Step 4 — verify the delivered inputs (once, ~3-4 hours)

```bash
bash ~/bola_ucd2/run_remote/02_verify_inputs.sh
```

Hashes all 13 delivered FASTQ on this host and confirms them against the audit
manifest, writing `00_input_audit/input_seal.json`.

This is separate from the preflight on purpose. It reads ~1.13 TiB, so folding
it into a verdict that expires in 6 hours would force a full re-hash after every
long job. The seal is bound to the delivered bytes and to the SHA-256 of the
audit manifest, not to a clock: it stays valid until one of them changes, and
re-running the audit invalidates it automatically. Every later step validates it
cheaply, without re-reading the delivery.

The pilot and cohort runners export it as `CODEX_INPUT_SEAL`, so each mapping
records `input_hash_verification: input_seal` — the delivered hash was confirmed
on this host, not carried over from the audit. Without a seal that field says
`not_verified`, because mapping only re-checks the file's size.

## Step 5 — pilot (2 animals)

```bash
bash ~/bola_ucd2/run_remote/03_run_pilot.sh        # both animals at once
bash ~/bola_ucd2/run_remote/03_run_pilot.sh 1      # or one at a time
```

Detached with `setsid`, so closing the SSH session does not stop it. Watch with
`bash run_remote/status.sh`.

## Step 6 — smoke and gate

```bash
bash ~/bola_ucd2/run_remote/04_pilot_gate.sh
```

For each pilot animal this runs `scripts/codex_reanalysis/pilot_smoke.py`, which
answers what alignment counts cannot:

- **secondary saturation** — `minimap2 -N` truncates the secondary list at a
  fixed ceiling. A tail decaying into the ceiling means it was not binding; a
  spike sitting exactly on it means alignments were still being discarded and
  the multi-hit evidence is censored.
- **split-read retention** (`sv_smoke`) — the retention filter keeps an
  alignment when it overlaps a target interval, and the read's other segments
  only if they overlap one too. Reads whose `SA` segments land outside the
  retained intervals lose the evidence structural variants are called from.
- **assembly fidelity** — a Flye assembly of the longest retained interval (the
  contiguous MHC stretch, derived from the BED rather than hard-coded) is mapped
  back to that reference interval and judged on span, breadth, identity and how
  many analytical target genes it recovers. Collapse and chimerism show up here
  and nowhere else in mapping QC.

Every threshold is passed explicitly and recorded in the smoke status, so the
criterion a verdict was reached under is reviewable rather than implied.

The assembly writes into the project tree, which step 1 put on **local disk**.
That matters: an identical smoke run over a Windows-mounted drive failed with an
empty consensus and a vanished intermediate BAM. If the project was linked back
onto the delivery drive instead of copied, move it before running this step.

If the assembler dies, the smoke still writes its status with `assembly_smoke`
failed and the reason recorded, so the gate reports what blocked the pilot
rather than only that something did.

`pilot_gate.py` then folds these into the pilot verdict alongside the input,
reference, mapping, index, target-coverage, control-depth and cohort-contrast
checks, and records the SHA-256 of every BAM, index, reference artefact and the
configuration the verdict was reached on. A smoke result describing a different
BAM is a failure, not evidence — a re-map cannot inherit an acceptance it never
earned.

## Step 7 — cohort (remaining 11 animals)

```bash
bash ~/bola_ucd2/run_remote/01_preflight.sh        # refresh the 6-hour verdict
bash ~/bola_ucd2/run_remote/05_run_cohort.sh       # uses the recommended concurrency
bash ~/bola_ucd2/run_remote/05_run_cohort.sh 1     # or force one sample at a time
```

`05_run_cohort.sh` refuses to start unless, in a single validated pass: the
preflight is fresh and `ready`, the mount is still read-only *now*, the input
seal matches the current audit manifest, and the pilot verdict is a plain
`passed` with no excluded checks and every smoke check passed. The delivery root
and job count it uses are the ones that verification returned, so the values
acted on are the values that were checked.

## Watching

```bash
bash ~/bola_ucd2/run_remote/status.sh
```

Shows both drivers, per-sample mapping status, and any pilot check that is not
`passed`.

## After a power loss

Each sample publishes its outputs only once every artefact is complete and its
alignment count reconciles with `samtools flagstat`. An interrupted sample
leaves a `running` status and no BAM, never a plausible-looking partial result.
Finished samples are unaffected and are skipped on restart.

```bash
bash ~/bola_ucd2/run_remote/99_unlock.sh    # clear the Snakemake lock
bash ~/bola_ucd2/run_remote/01_preflight.sh
bash ~/bola_ucd2/run_remote/05_run_cohort.sh
```

The input seal survives a power loss: it describes bytes, not a session.

## Returning the results

Only `results/codex_reanalysis_v1/02_mapping/` is needed back. Copy it over the
network, or remount the drive read-write once at the end and write it there.
