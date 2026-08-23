#!/usr/bin/env bash
# =============================================================================
# pp_pipeline_ideal.sh — IDEAL-FIELD CONTROL variant of pp_pipeline.sh
# (analysis session, 2026-08-20; user request: sim with NO distortion field,
#  steps g4hit->clusters rerun, everything else identical to V6).
# Byte-derived from pp_pipeline.sh with ONLY these deltas:
#   VER=v6ideal, FIELD="" (transport field OFF — the single physics change),
#   RAWPFX=raw_lib_ideal, and collision-free intermediate/log/census names
#   (fPPI_/dPPI_/i91ppI_/censI_/pp_census_ideal.txt) so nothing of the sealed
#   V6 production is touched. Export call unchanged (rowdr/field args empty
#   in V6 too -> identical). Run stages: gate prod (gen/g4 reuse V6 chunks).
# =============================================================================
set -eo pipefail

ROOTDIR=/home/rog/sPHENIX/3D_ClusterFindingML
P5=$ROOTDIR/P5
GEN=$P5/angantyr
IP=$ROOTDIR/island_post
LOGS=$IP/pp_logs; mkdir -p "$LOGS"

# ---------------- parameters ----------------
VER=v6ideal               # V6 (user naming, 2026-08-18): the laser-vetoed re-baseline; was provisionally v5.6
FIELD=""                  # IDEAL: no distortion field at transport
                          # v5.5 field-in-digitization (digi_field_request.md,
                          #  2026-08-13): the v5.4c empirical position-residual
                          #  field injected at transport (charge displaced in
                          #  r-phi BEFORE pad binning). Same measured params as
                          #  the v5.4c cluster overlay; the islandize91 export
                          #  then runs with the field OFF (rowdr only) so the
                          #  field enters exactly once. Empty string = off.
RAWPFX=raw_lib_ideal       # field-on library family — the sealed v5.3 libs
                          #  (raw_lib_pp_*) are left untouched.
NCHUNK=10                 # generator/G4/transport chunks
GEN_PER=2000              # pp events per chunk (pp collision ~1/25 pAu content
                          #  -> 10x the v4.0 event count for comparable library charge)
GEN_SEED_BASE=20260810    # chunk i uses GEN_SEED_BASE+i (v5.3 fresh seeds)
GEN_TUNE=1.85             # v5.3 soft-production tune: pT0Ref (Monash 2.28
                          #  default). Chosen 2026-07-25: 4-point census scan,
                          #  uniform closure -1.3% / fired anchor +2.9%.
G4_PAR=5                  # parallel G4 jobs per wave (fieldmap ~300 MB resident each)
NB=5; PER=50              # production: NB batches x PER events = total sim events
COMP_SEED_BASE=2026094    # composer seed = ${COMP_SEED_BASE}${i}
                          # v5.5 re-roll (2026-08-14, PRE-DECLARED): the
                          #  2026093x draw gave px/ev +5.8% vs target; the
                          #  measured composer realization sigma is ~3% per
                          #  50-frame block (~1.5-3% per production) — an
                          #  era-wide uncertainty never previously quoted.
                          #  Rule: accept this ONE re-roll if |px/ev-323799|
                          #  <= 2.5% (2 sigma); both draws ledgered.

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
ENV="5.80:0.968,11.10:1.064,16.40:0.965,21.70:1.055,27.00:1.055,32.30:1.036,37.60:0.913,42.90:1.012,47.96:0.931"
                          # ENV6 (2026-08-19): E2 residual correction of ENV53b against the
                          #  LASER-VETOED complete-61 real profile with the sim's own
                          #  window-length artefact removed (bin 10 compared over 874-960
                          #  only). ENV53b had been solved against profiles that held event
                          #  44's flash (real bin 3) and the sim's 44%-frequency injected
                          #  spike, and its last node against sim pixels at 961-970 that
                          #  the real DAQ never records. Ratios applied: 1.042 1.073 1.019
                          #  0.998 1.008 0.991 0.967 0.987 | node8 avg(0.976, 1.026).
                          #  Previous (ENV53b): 0.938,1.001,0.956,1.067,1.056,1.055,0.953,1.035,0.939
                          # v5.2 nodes (ENV52, 2026-07-24): references at the
                          #  FIXED composer (correct fired flags — trigger now
                          #  genuinely MBD-fired) + v5.1 electronics + eps 0.588.
                          #  Trigger-subtracted ratio vs COMPLETE-62, unit-mean.
                          #  Notably flatter than ENV51: the correct (bigger)
                          #  trigger subtraction removes the early-window tilt.
                          #  Anchor recalib: 280.3 px/unit -> 10,897.
FLASH=raw_lib_cmflash_w.root                        # CM flash lib (species-free)
FLASH_PROB=0.0            # v5.6 FINAL (2026-08-17, user): follow the official reco, which
                          #  vetoes GL1 laser events before TPC clustering -> NO flash
                          #  injection; the real reference drops event 44 (complete-61 /
                          #  all-99). Machinery parked: set >0 with SPEC to re-enable.
SPEC="2.2:1"              # (parked) giant-only spectrum, the yield fitted to event 44
FM=$ROOTDIR/CDB_offline/FIELDMAP_GAP/65/a9/65a930ed6de9c0e049cd0f3ef226e6b4_sphenix3dbigmapxyz_gap_rebuild_v2.root
DM=$ROOTDIR/CDB_offline/TPC_DEADCHANNELMAP/ff/c3/ffc3f6498934c5a8ba31065292c6ebcc_TPCDeadMap_79471.root
ANCHOR=10897              # rate-free fired anchor (complete-62 excess 38.87
                          #  x pp trig-only calib 280.3, v5.1 electronics)

LIBS=""; EVALS=""; MBDDATS=""
for i in $(seq 0 $((NCHUNK-1))); do
  LIBS+="${RAWPFX}_$i.root,"; EVALS+="eval_pp_$i.root,"; MBDDATS+="pp_run_$i.dat,"
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
    ./gen_pp $GEN_PER pp_run_$i.dat $((GEN_SEED_BASE+i)) $GEN_TUNE \
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
      tpc_transport(\"$P5/PP_g4hit_$i.root\",\"${RAWPFX}_$i.root\",$GEN_PER,$SIGMA0,$KPRF,\"$FIELD\");
      tpc_readout(\"${RAWPFX}_$i.root\",\"censI_pp_$i.root\",$GAIN,20.0,1,1,4711,\"$DM\",11.0,0.39,0.55,1.25,0.016,7.0,36.0,70.0,11.0,940.0,2,0.29,10.0,1.24,1.06,-1.0,5.0,0.0);" \
      > "$LOGS/pp_gateI_$i.log" 2>&1
    echo "[gate] chunk $i transported+censused"
  done
  cat > "$LOGS/pp_censusI.C" <<'EOM'
void pp_censusI(const char *fn, const char *chunk, const char *out)
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
  rm -f pp_census_ideal.txt
  for i in $(seq 0 $((NCHUNK-1))); do
    root -l -b -q "$LOGS/pp_censusI.C(\"censI_pp_$i.root\",\"$i\",\"pp_census_ideal.txt\")" > /dev/null 2>&1
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
for L in open("pp_census_ideal.txt"):
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
  rm -f censI_pp_*.root
fi

# ---------------- [prod] ----------------
if run_stage prod; then
  echo "=== [prod] composer -> readout -> islandize91, $NB x $PER events ==="
  cd "$IP"
  for i in $(seq 0 $((NCHUNK-1))); do
    ln -sf "$P5/PP_g4hit_$i.root" eval_pp_$i.root   # identity numbering: g4hit IS the eval
  done
  # LIB_SPLIT (2026-08-18): the composer holds every library pixel in memory
  # (~16 B x 1.36 G pixels ~ 22 GB for the 10-chunk pool) and dies at this
  # machine's ~27 GB ceiling (RAM 23 + swap 8; twice silently on 08-17/18).
  # Each 50-event batch now draws from HALF the pool (chunks 0-4 / 5-9,
  # alternating): 8.9k collisions per batch vs 3.25k draws — the pool is
  # still 2.7x oversampled, so per-batch reuse rises from ~0.18 to ~0.37
  # draws/collision; peak memory ~14 GB. The MBD map and truth tables follow
  # the same subset so (libsrc, event) bookkeeping stays aligned.
  LIB_SPLIT=${LIB_SPLIT:-2}
  for i in $(seq 0 $((NB-1))); do
    LIBS_B=""; EVALS_B=""; MBD_B=""
    for c in $(seq 0 $((NCHUNK-1))); do
      [ $(( c % LIB_SPLIT )) -eq $(( i % LIB_SPLIT )) ] || continue
      LIBS_B+="${RAWPFX}_$c.root,"; EVALS_B+="eval_pp_$c.root,"; MBD_B+="pp_run_$c.dat,"
    done
    LIBS_B=${LIBS_B%,}; EVALS_B=${EVALS_B%,}; MBD_B="pp_mbd.txt|${MBD_B%,}"
    root -l -b -q -e "
    gROOT->ProcessLine(\".L frame_composer.C+\");
    frame_composer(\"$LIBS_B\",\"fPPI_$i.root\",$PER,600.,${COMP_SEED_BASE}$i,$((i*PER)),\"$EVALS_B\",\"$FLASH\",$FLASH_PROB,1.0,\"$SPEC\",1.5,2.0,1.1,0.75,\"$RSPEC\",$ISO_SCALE,\"$MBD_B\",$TRIG_N,0.0,\"$ENV\");
    gROOT->ProcessLine(\".L tpc_digitize.C+\");
    tpc_readout(\"fPPI_$i.root\",\"dPPI_$i.root\",$GAIN,20.0,1,1,4711,\"$DM\",11.0,0.39,0.55,1.25,0.016,7.0,36.0,70.0,11.0,940.0,2,0.29,10.0,1.24,1.06,-1.0,5.0,1.0);" \
      2>&1 | tee "$LOGS/pp_prodI_$i.log" | grep -E "frame_composer: $PER|rate envelope|tpc_readout: fPPI"
    root -l -b -q "islandize91.C+(\"dPPI_$i.root\",\"i91ppI_$i.root\",1,\"\",\"fPPI_$i.root\",\"\",\"\")" 2>&1 | grep "islandize91: dPPI"
  done
  hadd -f frames_pp_production_${VER}.root fPPI_*.root  > /dev/null 2>&1
  hadd -f digi_frames_production_${VER}.root dPPI_*.root > /dev/null 2>&1
  hadd -f island91_frames_production_${VER}.root i91ppI_*.root > /dev/null 2>&1
  rm -f fPPI_*.root dPPI_*.root i91ppI_*.root
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
