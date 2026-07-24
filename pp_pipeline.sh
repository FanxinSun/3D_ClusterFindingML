#!/usr/bin/env bash
# =============================================================================
# pp_pipeline.sh — ALL-IN-ONE v5.0 pipeline driver: gen_pp.cc -> prodclus.C
# (species pivot 2026-07-22: real run 79507 is pp by scaler evidence;
#  the pAu v1-v4.0 line is closed — this is the species-correct chain.)
#
#   stage gen      Pythia8 pp MB (inelastic incl. diffraction) -> HepMC2 ASCII
#                  chunks + in-generator MBD proxy -> island_post/pp_mbd.txt
#   stage g4       standalone_tpc (local G4 11.4.2 + GDML + fieldmap) per chunk
#   stage gate     tpc_transport -> raw_lib_pp_* (kept AS the libs: NO THIN —
#                  v4.0's thinning was pAu species surgery, deliberately absent)
#                  + census readout -> PP GATE: MBD-fired mean kept-px vs the
#                  rate-free real anchor 8200. THE species validation number.
#   stage prod     frame_composer -> tpc_readout -> islandize91 per batch,
#                  hadd -> frames/digi/island91 _pp_production_${VER} + compact
#                  islandize
#   stage prodclus truth-labeled production clusters -> prodclus_${VER}.root
#
# usage:   ./pp_pipeline.sh all            # everything in order
#          ./pp_pipeline.sh gen g4         # just those stages
# Deliberately OUTSIDE this driver (post-production, run separately):
#   bestmatch_scan.C, funny_shapes.C, sim_to_csv.C, GenerateSimCSVs.C,
#   hits69.C, acceptance/figure macros.
#
# PP-RETUNE MARKERS: params tagged [PP-TBD] carry v4.0 values as starting
# points and are the retune campaign's knobs; everything untagged is
# detector-side (survives the species swap by construction).
# No `set -u`: geant4.sh is not nounset-clean (same reason as the old scripts).
# =============================================================================
set -eo pipefail

ROOTDIR=/home/rog/sPHENIX/3D_ClusterFindingML
P5=$ROOTDIR/P5
GEN=$P5/angantyr
IP=$ROOTDIR/island_post
LOGS=$IP/pp_logs; mkdir -p "$LOGS"

# ---------------- parameters ----------------
VER=v52
NCHUNK=10                 # generator/G4/transport chunks
GEN_PER=2000              # pp events per chunk (pp collision ~1/25 pAu content
                          #  -> 10x the v4.0 event count for comparable library charge)
GEN_SEED_BASE=20260722    # chunk i uses GEN_SEED_BASE+i
G4_PAR=5                  # parallel G4 jobs per wave (fieldmap ~300 MB resident each)
NB=5; PER=50              # production: NB batches x PER events = total sim events
COMP_SEED_BASE=2026093    # composer seed = ${COMP_SEED_BASE}${i}

SIGMA0=0.040; KPRF=0.0070 # transport PRF sigma(Ld) — detector-side (v4.0 values)
GAIN=1.005                # v5.1 LOCKED 2026-07-22 (scan pass B2: pixmean +0.8%)
ISO_SCALE=0.27            # [PP-TBD] iso singles scale — 0.27 was the Angantyr
                          #  content balance; re-derive at pp content
TRIG_N=1.0                # physical trigger (no fired-mean surgery)
RSPEC=$(grep -v '^#' "$IP/rspec99_v52.txt" 2>/dev/null)
[ -n "$RSPEC" ] || { echo "rspec99_v52.txt missing/empty in $IP"; exit 1; }
                          # PHYSICAL, zero fitted parameters (2026-07-22):
                          #  all 99 measured per-event rmbd values / eps 0.519
                          #  (proxy; replace file when team eps arrives).
                          #  Window accounting judged by acceptance, not pre-fit.
ENV="5.80:0.982,11.10:1.008,16.40:1.011,21.70:1.028,27.00:1.009,32.30:1.026,37.60:1.016,42.90:1.017,47.96:0.904"
                          # v5.2 nodes (ENV52, 2026-07-24): references at the
                          #  FIXED composer (correct fired flags — trigger now
                          #  genuinely MBD-fired) + v5.1 electronics + eps 0.588.
                          #  Trigger-subtracted ratio vs COMPLETE-62, unit-mean.
                          #  Notably flatter than ENV51: the correct (bigger)
                          #  trigger subtraction removes the early-window tilt.
                          #  Anchor recalib: 280.3 px/unit -> 10,897.
FLASH=raw_lib_cmflash_w.root                        # CM flash lib (species-free)
SPEC="0.008:1,0.009:1,0.011:1,0.012:1,0.013:1,0.014:1,0.018:1,0.021:1,0.027:1,0.037:1,2.2:0.25"
FM=$ROOTDIR/CDB_offline/FIELDMAP_GAP/65/a9/65a930ed6de9c0e049cd0f3ef226e6b4_sphenix3dbigmapxyz_gap_rebuild_v2.root
DM=$ROOTDIR/CDB_offline/TPC_DEADCHANNELMAP/ff/c3/ffc3f6498934c5a8ba31065292c6ebcc_TPCDeadMap_79471.root
ANCHOR=8200               # rate-free real fired-collision kept-px anchor

LIBS=""; EVALS=""; MBDDATS=""
for i in $(seq 0 $((NCHUNK-1))); do
  LIBS+="raw_lib_pp_$i.root,"; EVALS+="eval_pp_$i.root,"; MBDDATS+="pp_run_$i.dat,"
done
LIBS=${LIBS%,}; EVALS=${EVALS%,}; MBD="pp_mbd.txt|${MBDDATS%,}"

run_stage() { case " $STAGES " in *" $1 "*|*" all "*) return 0;; *) return 1;; esac; }
STAGES="$*"; [ -z "$STAGES" ] && { echo "usage: $0 [gen g4 gate prod prodclus | all]"; exit 1; }

# ---------------- [gen] ----------------
if run_stage gen; then
  echo "=== [gen] Pythia8 pp MB: $NCHUNK x $GEN_PER events ==="
  cd "$GEN"
  if [ ! -x gen_pp ] || [ gen_pp.cc -nt gen_pp ]; then
    g++ -O2 -std=c++17 gen_pp.cc -o gen_pp \
        -I install/include -L install/lib -lpythia8 -ldl \
        -Wl,-rpath,"$GEN/install/lib"
    echo "gen_pp compiled"
  fi
  for i in $(seq 0 $((NCHUNK-1))); do
    ./gen_pp $GEN_PER pp_run_$i.dat $((GEN_SEED_BASE+i)) \
        > "$LOGS/gen_pp_$i.log" 2> fired_pp_$i.txt &
  done
  wait
  # plain `wait` returns 0 even if a background job died — verify explicitly
  NDONE=$(grep -hl GEN-DONE "$LOGS"/gen_pp_*.log 2>/dev/null | wc -l)
  [ "$NDONE" -eq "$NCHUNK" ] || { echo "[gen] FAILED: only $NDONE/$NCHUNK chunks completed (see $LOGS/gen_pp_*.log)"; exit 1; }
  for i in $(seq 0 $((NCHUNK-1))); do
    [ -s pp_run_$i.dat ] || { echo "[gen] FAILED: pp_run_$i.dat empty/missing"; exit 1; }
  done
  grep -h GEN-DONE "$LOGS"/gen_pp_*.log
  # consolidate the MBD proxy in the composer's mapping format:
  #   <libfile> <event_index_0based> <north> <south> <fired>
  { echo "# file event_index_in_file north south fired"
    for i in $(seq 0 $((NCHUNK-1))); do
      awk -v f="pp_run_$i.dat" '$1=="FIRED"{print f, $2-1, $3, $4, $5}' fired_pp_$i.txt
    done; } > "$IP/pp_mbd.txt"
  echo "[gen] DONE: $(grep -c "^pp_run" "$IP/pp_mbd.txt") events in pp_mbd.txt," \
       "fired frac $(awk '$1!~/^#/{n++;s+=$5} END{printf "%.3f", s/n}' "$IP/pp_mbd.txt")"
fi

# ---------------- [g4] ----------------
if run_stage g4; then
  echo "=== [g4] standalone_tpc: $NCHUNK chunks, waves of $G4_PAR ==="
  [ -x "$P5/standalone_tpc" ] || { echo "standalone_tpc missing — run $P5/build_standalone.sh"; exit 1; }
  cd "$P5"
  source /home/rog/geant4/bin/geant4.sh
  for i in $(seq 0 $((NCHUNK-1))); do
    ./standalone_tpc sphenix_p5.gdml angantyr/pp_run_$i.dat PP_g4hit_$i.root \
        $GEN_PER "$FM" > "$LOGS/pp_g4_$i.log" 2>&1 &
    [ $(( (i+1) % G4_PAR )) -eq 0 ] && wait
  done
  wait
  for i in $(seq 0 $((NCHUNK-1))); do
    [ -s PP_g4hit_$i.root ] || { echo "[g4] FAILED: PP_g4hit_$i.root empty/missing (see $LOGS/pp_g4_$i.log)"; exit 1; }
    grep -q "TpcSD: wrote" "$LOGS/pp_g4_$i.log" || { echo "[g4] FAILED: chunk $i did not finish (see $LOGS/pp_g4_$i.log)"; exit 1; }
  done
  echo "[g4] DONE: $(ls -1 PP_g4hit_*.root | wc -l) g4hit files"
fi

# ---------------- [gate] ----------------
if run_stage gate; then
  echo "=== [gate] transport -> libs (NO THIN) + census -> anchor gate ==="
  cd "$IP"
  for i in $(seq 0 $((NCHUNK-1))); do
    root -l -b -q -e "gROOT->ProcessLine(\".L tpc_digitize.C+\");
      tpc_transport(\"$P5/PP_g4hit_$i.root\",\"raw_lib_pp_$i.root\",$GEN_PER,$SIGMA0,$KPRF);
      tpc_readout(\"raw_lib_pp_$i.root\",\"cens_pp_$i.root\",$GAIN,20.0,1,1,4711,\"$DM\",11.0,0.39,0.55,0.81,0.016,7.0,36.0,70.0,11.0,940.0,2,0.29,10.0,1.24,1.06,-1.0,5.0,0.0);" \
      > "$LOGS/pp_gate_$i.log" 2>&1
    echo "[gate] chunk $i transported+censused"
  done
  cat > "$LOGS/pp_census.C" <<'EOM'
void pp_census(const char *fn, const char *chunk, const char *out)
{
  TFile *f = TFile::Open(fn);
  TTree *t = (TTree *) f->Get("ntp_hit");
  float ev; t->SetBranchStatus("*", 0);
  t->SetBranchStatus("event", 1); t->SetBranchAddress("event", &ev);
  std::map<int, long> n;
  for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); n[(int) ev]++; }
  FILE *fo = fopen(out, "a");
  for (auto &kv : n) fprintf(fo, "%s %d %ld\n", chunk, kv.first, kv.second);
  fclose(fo);
}
EOM
  rm -f pp_census.txt
  for i in $(seq 0 $((NCHUNK-1))); do
    root -l -b -q "$LOGS/pp_census.C(\"cens_pp_$i.root\",\"$i\",\"pp_census.txt\")" > /dev/null 2>&1
  done
  python3 - "$ANCHOR" <<'EOP'
import sys
anchor = float(sys.argv[1])
fired = {}
for L in open("pp_mbd.txt"):
    if L.startswith("#"): continue
    f, e, n, s, fl = L.split()
    fired[(f.split("_")[-1].split(".")[0], int(e))] = int(fl)
cnt = {}
for L in open("pp_census.txt"):
    c, e, n = L.split()
    cnt[(c, int(e))] = int(n)
allv = [cnt.get(k, 0) for k in fired]
fv = [cnt.get(k, 0) for k, fl in fired.items() if fl]
mean = lambda v: sum(v) / max(1, len(v))
fm, um = mean(fv), mean(allv)
print(f"PP GATE: fired-mean kept px {fm:.0f} vs anchor {anchor:.0f} ({100*(fm/anchor-1):+.1f}%)")
print(f"         uniform mean {um:.0f} | fired {len(fv)}/{len(allv)} ({len(fv)/max(1,len(allv)):.3f}) | trigger bias x{fm/max(1e-9,um):.3f}")
print("         PASS band +-15%: " + ("PASS" if abs(fm/anchor-1) < 0.15 else "OUTSIDE — species/tune question BEFORE production"))
EOP
  rm -f cens_pp_*.root
fi

# ---------------- [prod] ----------------
if run_stage prod; then
  echo "=== [prod] composer -> readout -> islandize91, $NB x $PER events ==="
  cd "$IP"
  for i in $(seq 0 $((NCHUNK-1))); do
    ln -sf "$P5/PP_g4hit_$i.root" eval_pp_$i.root   # identity numbering: g4hit IS the eval
  done
  for i in $(seq 0 $((NB-1))); do
    root -l -b -q -e "
    gROOT->ProcessLine(\".L frame_composer.C+\");
    frame_composer(\"$LIBS\",\"fPP_$i.root\",$PER,600.,${COMP_SEED_BASE}$i,$((i*PER)),\"$EVALS\",\"$FLASH\",0.44,1.0,\"$SPEC\",1.5,2.0,1.1,0.75,\"$RSPEC\",$ISO_SCALE,\"$MBD\",$TRIG_N,0.0,\"$ENV\");
    gROOT->ProcessLine(\".L tpc_digitize.C+\");
    tpc_readout(\"fPP_$i.root\",\"dPP_$i.root\",$GAIN,20.0,1,1,4711,\"$DM\",11.0,0.39,0.55,0.81,0.016,7.0,36.0,70.0,11.0,940.0,2,0.29,10.0,1.24,1.06,-1.0,5.0,1.0);" \
      2>&1 | tee "$LOGS/pp_prod_$i.log" | grep -E "frame_composer: $PER|rate envelope|tpc_readout: fPP"
    root -l -b -q "islandize91.C+(\"dPP_$i.root\",\"i91pp_$i.root\",1,\"\",\"fPP_$i.root\")" 2>&1 | grep "islandize91: dPP"
  done
  hadd -f frames_pp_production_${VER}.root fPP_*.root  > /dev/null 2>&1
  hadd -f digi_frames_production_${VER}.root dPP_*.root > /dev/null 2>&1
  hadd -f island91_frames_production_${VER}.root i91pp_*.root > /dev/null 2>&1
  rm -f fPP_*.root dPP_*.root i91pp_*.root
  root -l -b -q -e "gROOT->ProcessLine(\".L islandize.C+\"); islandize(\"digi_frames_production_${VER}.root\",\"island_frames_${VER}.root\",1);" 2>&1 | grep islandize
  echo "[prod] DONE: frames/digi/island91 _pp_production_${VER} + island_frames_${VER}"
fi

# ---------------- [prodclus] ----------------
if run_stage prodclus; then
  echo "=== [prodclus] truth-labeled production clusters ==="
  cd "$IP"
  root -l -b -q -e "gROOT->ProcessLine(\".L prodclus.C+\");
    prodclus(\"digi_frames_production_${VER}.root\",\"prodclus_${VER}.root\",1,5.0,3.0,5.0,0,10,20,1,\"\",\"frames_pp_production_${VER}.root\");" \
    2>&1 | grep prodclus:
  echo "[prodclus] DONE: prodclus_${VER}.root"
fi

echo "=== PP PIPELINE (${VER}): requested stages complete ==="
