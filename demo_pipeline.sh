#!/usr/bin/env bash
# =============================================================================
# demo_pipeline.sh — ~2-minute end-to-end demonstration of the detached pp
# pipeline at toy scale. Runs the SAME code and the SAME physics knobs as
# production (they are read out of pp_pipeline.sh at runtime), only with a
# handful of collisions and frames. EVERY output goes to /tmp/pp_demo — the
# repo, the sealed libraries and the production files are never touched.
#
#   usage:  bash demo_pipeline.sh [n_collisions] [n_frames]
#           defaults: 150 collisions, 3 frames
#
# Prerequisite: bash setup_machine.sh has been run and env sourced.
# See FIRST_RUN.md for the guided version of these steps.
# =============================================================================
set -eo pipefail

R="$(cd "$(dirname "$0")" && pwd)"
P5=$R/P5; GEN=$P5/angantyr; IP=$R/island_post; PIPE=$R/pp_pipeline.sh
OUT=/tmp/pp_demo; mkdir -p "$OUT"
NEV=${1:-150}      # collisions to generate + transport
NFR=${2:-3}        # frames to compose
SEED=20260816

# inherit the production configuration (so this demo is the real config, small)
ROOTDIR=$R
eval "$(grep -E '^(FIELD|TWIST|GEN_TUNE|SIGMA0|KPRF|GAIN|ISO_SCALE|TRIG_N|ENV|SPEC|FLASH|FM|DM)=' "$PIPE")"
RSPEC=$(grep -v '^#' "$IP/rspec99_v52.txt")
echo "=== demo config (inherited from pp_pipeline.sh) ==="
echo "    collisions $NEV | frames $NFR | tune pT0Ref $GEN_TUNE"
echo "    field \"$FIELD\" | twist $TWIST"
echo "    transport sigma $SIGMA0/$KPRF | gain $GAIN | iso $ISO_SCALE | trig $TRIG_N"
T0=$(date +%s); step(){ printf "    [%3ds] %s\n" $(( $(date +%s) - T0 )) "$1"; }

# ---------------- 1. generator ----------------
echo "=== [1/6] gen: Pythia8 pp MB, $NEV events ==="
cd "$GEN"
./gen_pp "$NEV" "$OUT/demo_run_0.dat" "$SEED" "$GEN_TUNE" \
    > "$OUT/log_gen.txt" 2> "$OUT/fired_demo_0.txt"
{ echo "# file event_index_in_file north south fired"
  awk -v f=demo_run_0.dat '$1=="FIRED"{print f, $2-1, $3, $4, $5}' "$OUT/fired_demo_0.txt"; } > "$OUT/demo_mbd.txt"
step "gen done: $(grep -c '^demo_run_0' "$OUT/demo_mbd.txt") events, fired frac $(awk '$1!~/^#/{n++;s+=$5} END{printf "%.3f", s/n}' "$OUT/demo_mbd.txt")"

# ---------------- 2. Geant4 ----------------
echo "=== [2/6] g4: standalone_tpc (local Geant4 + GDML + field map) ==="
cd "$P5"
source "$HOME/geant4/bin/geant4.sh"
./standalone_tpc sphenix_p5.gdml "$OUT/demo_run_0.dat" "$OUT/demo_g4hit_0.root" \
    "$NEV" "$FM" > "$OUT/log_g4.txt" 2>&1
step "g4 done: $(grep -o 'TpcSD: wrote [0-9]*' "$OUT/log_g4.txt" | tail -1)"

# ---------------- 3. transport (field + twist act here) ----------------
echo "=== [3/6] transport: ionization -> drift+diffusion -> FIELD+TWIST -> pads ==="
cd "$IP"
root -l -b -q -e "gROOT->ProcessLine(\".L tpc_digitize.C+\");
  tpc_transport(\"$OUT/demo_g4hit_0.root\",\"$OUT/demo_lib_0.root\",$NEV,$SIGMA0,$KPRF,\"$FIELD\",\"$TWIST\");" \
  > "$OUT/log_transport.txt" 2>&1
step "transport: $(grep -o 'rphi field IN-DIGI.*' "$OUT/log_transport.txt" | head -1)"
step "transport: $(grep -o 'TWIST in-digi.*' "$OUT/log_transport.txt" | head -1 | cut -c1-88)"

# ---------------- 4. composer + readout ----------------
echo "=== [4/6] composer -> readout: $NFR frames of pileup, then electronics ==="
root -l -b -q -e "gROOT->ProcessLine(\".L frame_composer.C+\");
  frame_composer(\"$OUT/demo_lib_0.root\",\"$OUT/demo_frames.root\",$NFR,600.,$SEED,0,
    \"$OUT/demo_g4hit_0.root\",\"$FLASH\",0.44,1.0,\"$SPEC\",1.5,2.0,1.1,0.75,\"$RSPEC\",$ISO_SCALE,
    \"$OUT/demo_mbd.txt|demo_run_0.dat\",$TRIG_N,0.0,\"$ENV\");
  gROOT->ProcessLine(\".L tpc_digitize.C+\");
  tpc_readout(\"$OUT/demo_frames.root\",\"$OUT/demo_digi.root\",$GAIN,20.0,1,1,4711,\"$DM\",
    11.0,0.39,0.55,0.81,0.016,7.0,36.0,70.0,11.0,940.0,2,0.29,10.0,1.24,1.06,-1.0,5.0,1.0);" \
  > "$OUT/log_prod.txt" 2>&1
step "readout done: $(grep -o 'tpc_readout:.*' "$OUT/log_prod.txt" | tail -1 | cut -c1-96)"

# ---------------- 5. clustering ----------------
echo "=== [5/6] islandize91: clusters in the real 91-branch layout (+ truth) ==="
root -l -b -q "islandize91.C+(\"$OUT/demo_digi.root\",\"$OUT/demo_island91.root\",1,\"\",
  \"$OUT/demo_frames.root\",\"\",\"\")" > "$OUT/log_island.txt" 2>&1
step "clustering: $(grep -o 'islandize91: .*clusters.*' "$OUT/log_island.txt" | tail -1 | cut -c1-96)"

# ---------------- 6. production clusterizer ----------------
echo "=== [6/6] prodclus: truth-labeled ML cluster set ==="
root -l -b -q -e "gROOT->ProcessLine(\".L prodclus.C+\");
  prodclus(\"$OUT/demo_digi.root\",\"$OUT/demo_prodclus.root\",1,5.0,3.0,5.0,0,10,20,1,\"\",
  \"$OUT/demo_frames.root\");" > "$OUT/log_prodclus.txt" 2>&1
step "prodclus: $(grep -o 'prodclus: .*clusters.*' "$OUT/log_prodclus.txt" | tail -1 | cut -c1-96)"

# ---------------- summary ----------------
echo "=== summary: what each stage produced (all in $OUT) ==="
root -l -b -q -e "
  auto n=[](const char* f,const char* t){ TFile* x=TFile::Open(f); if(!x||x->IsZombie()) return -1LL;
    TTree* y=(TTree*)x->Get(t); return y? y->GetEntries() : -1LL; };
  printf(\"  truth g4 steps (ntp_g4hit)  : %lld\n\", n(\"$OUT/demo_g4hit_0.root\",\"ntp_g4hit\"));
  printf(\"  transported pixels (raw_pix): %lld\n\", n(\"$OUT/demo_lib_0.root\",\"raw_pix\"));
  printf(\"  composed frame pixels       : %lld\n\", n(\"$OUT/demo_frames.root\",\"raw_pix\"));
  printf(\"  digitized pixels (ntp_hit)  : %lld   <- real-data-shaped output\n\", n(\"$OUT/demo_digi.root\",\"ntp_hit\"));
  printf(\"  island clusters             : %lld (+ ntp_truth %lld)\n\", n(\"$OUT/demo_island91.root\",\"ntp_cluster\"), n(\"$OUT/demo_island91.root\",\"ntp_truth\"));
  printf(\"  prodclus clusters (ntp_clus): %lld\n\", n(\"$OUT/demo_prodclus.root\",\"ntp_clus\"));
" 2>/dev/null | grep -E "g4 steps|transported|composed|digitized|island|prodclus" || true
ls -lh "$OUT"/*.root | awk '{printf "    %-28s %s\n", $NF, $5}'
step "TOTAL"
