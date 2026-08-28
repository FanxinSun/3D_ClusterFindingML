# First run — the detached sPHENIX TPC simulation pipeline (v6.1)

You have just cloned this repository onto a new machine. This page takes you
from nothing to (a) a working environment, (b) a two-minute demonstration that
the whole chain runs, and (c) a full production run and its post-production
battery. Read it top to bottom the first time; afterwards §4 is the only part
you will keep returning to.

**What this pipeline is.** It simulates the sPHENIX TPC end to end — proton-proton
collisions, Geant4 transport, charge drift, electronics, zero suppression,
clustering — and produces pixel and cluster datasets shaped exactly like the real
detector's output, for machine-learning cluster finding. It is *detached*: no
container, no CVMFS, no BNL account. The reference dataset it is calibrated
against is run 79507 (pp, 200 GeV).

---

## 0. Prerequisites in one command

```bash
bash setup_machine.sh            # Linux/WSL (x86_64 or arm64) and macOS (Apple Silicon)
bash setup_machine.sh --check    # report-only: tells you what is missing, changes nothing
```

It installs or builds ROOT 6.40.02, Geant4 11.4.2 (GDML + multithreaded) and
Pythia 8.317, ports the repository's hardcoded paths if your `$HOME` differs from
the reference machine, writes `env_v61.sh`, and builds the two in-repo binaries
(`standalone_tpc`, `gen_pp`). Budget an hour or two, mostly Geant4.

Then, **in every shell you work in**:

```bash
source env_v61.sh
```

### The four files git does not carry

They are large or private, so they are copied in by hand. `setup_machine.sh`
ends by checking for exactly these and printing OK/MISSING:

| file | size | what it is |
|---|---|---|
| `clusters_seeds_island_79507-0.root_ntuplizer.root` | 440 MB | **the real data** — the only authority every meter compares against |
| `CDB_offline/FIELDMAP_GAP/65/a9/…gap_rebuild_v2.root` | 299 MB | magnetic field map (26,144,996 grid points), used by Geant4 |
| `CDB_offline/TPC_DEADCHANNELMAP/ff/c3/…TPCDeadMap_79471.root` | small | dead-channel mask, used by the readout |
| `island_post/raw_lib_cmflash_w.root` | 3.8 MB | central-membrane laser-flash library, used by the composer |

Optionally also copy the v6.1 production artifacts (`island_post/*_v61.root`,
~2.5 GB) if you want to *analyse* v6.1 rather than *re-run* it.

> **Reproducibility note.** The chain runs identically on Apple Silicon, but a
> from-scratch generator+Geant4 run there is a *statistically equivalent, not
> identical* production: floating-point differences diverge the transport, so
> md5s will not match `island_post/production_manifest.txt`. Byte-level
> reproduction of a sealed artifact is only expected on x86_64 Linux with the
> reference stack. Copy the artifacts if you need *the* v6.1.

---

## 1. Two-minute demo — prove the chain works

Everything lands in `/tmp/pp_demo`; nothing in the repository is touched.

```bash
bash demo_pipeline.sh            # 150 collisions, 3 frames, ~90 seconds
```

It runs all six stages with the **real production configuration**, read out of
`pp_pipeline.sh` at runtime — only the event count is small. Expected output:

```
    collisions 150 | frames 3 | tune pT0Ref 1.85
    field "0:0|0|0.009|0|0|2.49|0.26" | twist twist_payload_v6.txt
[1/6] gen        →  150 events, fired frac 0.493
[2/6] g4         →  TpcSD: wrote 523745
[3/6] transport  →  rphi field IN-DIGI: SMOD 0.0090 … SPHI 2.4900 SCM 0.2600
                    TWIST in-digi from twist_payload_v6.txt (96 rows)
[4/6] readout    →  998377 pixels kept, 141744 dead-masked
[5/6] clustering →  95260 clusters [track 48715, looper 37487, noise 9058]
[6/6] prodclus   →  110157 clusters [track 53027, looper 48960, noise 8170]
```

Two lines are worth pausing on. The **transport** line is the position-distortion
model announcing itself — the field and the measured "twist" are applied to the
drifting charge *before* it is binned into pads, which is where a real distortion
acts. And 998,377 pixels over 3 frames is **333k pixels/frame**, which is the
production content number: the toy run is physically representative, not a
caricature.

### Look at what you made

```bash
ls -lh /tmp/pp_demo/*.root
root -l /tmp/pp_demo/demo_digi.root
```
then at the ROOT prompt:
```cpp
ntp_hit->Print("toponly")                                    // 10 branches
ntp_hit->Draw("adc","layer>=7&&adc>0")                       // ADC spectrum
ntp_hit->Draw("tbin:phi>>h(200,0,1,200,600,800)",
              "adc*(event==0&&layer==15)","colz")            // cluster shapes
```

`rm -rf /tmp/pp_demo` clears the 128 MB when you are done.

Two caveats so nobody over-reads the demo: with 150 collisions the composer
re-draws from a small pool, so the three frames share collisions (fine for a
demo, not for physics), and a single frame's content fluctuates far beyond the
±2–3 % realization spread that applies to a 250-frame production.

### If the demo fails

| symptom | cause and fix |
|---|---|
| `standalone_tpc missing` | `bash P5/build_standalone.sh` (needs Geant4 + xerces-c) |
| `cannot open …gap_rebuild_v2.root` / `TPCDeadMap…` | payload file not copied — see §0 |
| `GEO LOAD FAILED` | a ROOT macro was run from the wrong directory: stages 3–6 must run **from `island_post/`** (`tpc_geom_table.txt` is read relative to the working directory) |
| `unbound variable` right after Geant4 | `geant4.sh` is not `set -u`-clean; do not use `set -u` in that shell |

---

## 2. Full production run

Everything up to production lives in one driver, `pp_pipeline.sh`. The
post-production steps are deliberately outside it and run by hand.

### Step 0 — mint the new version names *before* launching

Edit the parameter block at the top of `pp_pipeline.sh`. At minimum set `VER`;
if you changed anything in transport or the field, also set a fresh `RAWPFX` so
the sealed libraries are not overwritten:

```bash
VER=v62                                  # every output file is named from this
RAWPFX=raw_lib_pp62                      # only if FIELD/TWIST/SIGMA0/KPRF changed
FIELD="0:0|0|0.009|0|0|2.49|0.26"        # v6.1 field model
TWIST="twist_payload_v6.txt"             # measured, byte-frozen — never rescale
COMP_SEED_BASE=2026095                   # composer seeds; see the σ note in Step 5
```

Do this **first**. The project learned it the hard way at v3.7: an inherited
script still wrote `_v36` names and nearly clobbered a sealed production.

### Step 1 — check the machine

```bash
source env_v61.sh
root-config --version                 # expect 6.40.02
ls $HOME/geant4/bin/geant4.sh         # Geant4 11.4.2
free -g                               # see Step 5 — the composer needs headroom
```

`gen_pp` rebuilds itself if `gen_pp.cc` is newer; `standalone_tpc` does not — the
g4 stage aborts with a message if it is missing.

### Step 2 — `gen`: the generator (~5–15 min, 10 chunks in parallel)

```bash
./pp_pipeline.sh gen
```

Pythia 8.317 pp minimum bias, 10 × 2,000 events → `P5/angantyr/pp_run_0..9.dat`
plus the MBD-trigger proxy consolidated into `island_post/pp_mbd.txt`.
**Verify:** it prints `fired frac` (expect ≈ 0.52) and hard-fails if any chunk
did not reach `GEN-DONE`.

### Step 3 — `g4`: Geant4 transport (~20 min, waves of 5)

```bash
./pp_pipeline.sh g4
```

`standalone_tpc` per chunk over the GDML geometry and the 26.1 M-point field map
→ `P5/PP_g4hit_0..9.root`. **Verify:** it checks every log for `TpcSD: wrote` and
that no output is empty. This is the expensive stage you almost never repeat — it
is species/geometry-level and blind to every detector knob below it.

### Step 4 — `gate`: transport into libraries + the anchor gate (~60 min, serial)

```bash
./pp_pipeline.sh gate
```

`tpc_transport` runs here **with the field and the twist**, producing
`${RAWPFX}_0..9.root`, plus a census readout used only to compute the gate.
**Stop and read the output — this is a decision point, not a formality:**

```
PP GATE: fired-mean kept px NNNNN vs anchor 10897 (+X.X%)
         uniform mean … | fired …/… | trigger bias xN.NNN
         PASS band +-15%: PASS / OUTSIDE
```

If it says OUTSIDE, the species or the tune is wrong, and continuing only makes a
bad production more expensive.

### Step 5 — `prod`: composer → readout → clustering (~60–80 min)

**Check memory first.** The composer's library load peaks near 27–30 GB; the v5.5
run died here when stale editor/session processes had pinned most of the swap.

```bash
free -g                       # want ~30 GB total available
./pp_pipeline.sh prod
```

Per batch it runs `frame_composer` → `tpc_readout` → `islandize91`, then merges
the five batches into `frames_pp_production_${VER}.root`,
`digi_frames_production_${VER}.root` (**the pixels**),
`island91_frames_production_${VER}.root`, and finally `islandize` →
`island_frames_${VER}.root`.

> **Discipline point.** Per-production content scatters by ±2–3 % from composer
> realization noise alone. Decide your seed-acceptance rule *before* you look at
> the number — the way the v5.5 re-roll was pre-declared — or you are seed-shopping.

### Step 6 — `prodclus`: truth-labeled production clusters (~minutes)

```bash
./pp_pipeline.sh prodclus
```

Writes `prodclus_${VER}.root`, the ported production clusterizer's output with
truth labels — the ML-facing cluster set.

### Running the whole thing

```bash
./pp_pipeline.sh all
# long runs:
nohup ./pp_pipeline.sh all > island_post/pp_logs/run_$(date +%F).log 2>&1 &
tail -f island_post/pp_logs/run_$(date +%F).log
```

The script uses `set -eo pipefail`: it halts at the first failure, and you resume
by naming the remaining stages (`./pp_pipeline.sh prod prodclus`). Per-stage logs
land in `island_post/pp_logs/` as `pp_gen_*.log`, `pp_g4_*.log`, `pp_gate_*.log`,
`pp_prod_*.log`.

---

## 3. Post-production, by hand from `island_post/`

```bash
cd island_post

# 3a. real-layout hit tree (the pixel-consumer/supervisor artifact), ~9 min
root -l -b -q 'hits69.C+("digi_frames_production_v62.root",
  "island91_frames_production_v62.root","hit69_frames_production_v62.root")'

# 3b. acceptance battery
root -l -b -q 'fieldmeter.C+("../clusters_seeds_island_79507-0.root_ntuplizer.root",
  "island91_frames_production_v62.root","v62")'
root -l -b -q 'cmcheck.C+("island91_frames_production_v62.root",1,"v62")'
root -l -b -q 'ms_real.C+' -e 'ms_real_split("../clusters_seeds_island_79507-0.root_ntuplizer.root",
  "island91_frames_production_v62.root","v62","v6.2")'
root -l -b -q '../sim_validation_plots/src/ms_nofinder.C+' -e 'nf_digipix("digi_frames_production_v62.root",60,
  "../clusters_seeds_island_79507-0.root_ntuplizer.root","v62")'

# 3c. figures, then seal
md5sum digi_frames_production_v62.root island91_frames_production_v62.root \
       prodclus_v62.root hit69_frames_production_v62.root    # → production_manifest.txt
```

Then write the `PIPELINE.md` entry and commit. **The seal is not bookkeeping** —
the manifest md5s are what let a later reviewer prove the file on disk is the file
the numbers came from.

A production is not considered delivered until its figures, residual table and
briefing have been regenerated on the new artifacts. Never defer replots because
a future retune "might" re-cut them.

---

## 4. The part that saves you the most time

Which stages a change actually forces:

| what you changed | re-run | cost |
|---|---|---|
| species, generator tune, seeds | `gen g4 gate prod prodclus` | ~2.5 h |
| field, twist, `SIGMA0`, `KPRF` (transport) | `gate prod prodclus` | ~2 h |
| composer: `RSPEC`, `ENV`, `iso`, seeds | `prod prodclus` | ~1.5 h |
| electronics: gain, ZS, tails, `σ_ped` | `prod prodclus` | ~1.5 h |
| export convention only | `islandize91` + `hits69` | minutes |

This staging is the whole point of the detached design: the Geant4 libraries are
frozen currency, so a readout question costs an hour, not a day. When you only
need to compare electronics settings, do not touch `g4` — and never re-run
`gen`/`g4` "to be safe", because new seeds there change the physics realization
and make your comparison non-comparable.

---

## 5. Where everything lives

| what | where |
|---|---|
| production driver | `pp_pipeline.sh` |
| demo | `demo_pipeline.sh` → `/tmp/pp_demo` |
| environment bootstrap | `setup_machine.sh` → `env_v61.sh` |
| simulation code (transport, composer, readout, clustering, exporters) | `island_post/*.C` |
| standalone Geant4 app + geometry + generator | `P5/` |
| figures, plotting macros, ledgers | `sim_validation_plots/` (**git submodule** — clone with `--recurse-submodules`) |
| the project ledger: every version, calibration and decision | `PIPELINE.md` |
| real-data findings | `REAL_DATA_PROBES.md` |
| independent review verdicts | `SIDE_REVIEW_VERDICTS.md` |
| cross-thread work requests | `island_post/*_request.md` |
| artifact checksums and provenance | `island_post/production_manifest.txt` |

Start with `PIPELINE.md` if you want the history of *why* the configuration is
what it is; it is the ledger the whole project is built on.
