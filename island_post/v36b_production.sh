#!/usr/bin/env bash
# P3 production v3.6 revB (2026-07-14): EMPIRICAL RATE ENVELOPE (leg d, final form).
# Supersedes the exp(-|t|/tau) envelope of the first v3.6 pass:
#   - The exponential form was misspecified: after the tau=2060 lock the
#     window-fraction gate (W1[60,240]/W2[270,600]/W3[650,950], flash excluded)
#     still sat at +5.2/-7.7/+5.8% - too peaked early, too flat late.
#   - Replaced by a PIECEWISE-LINEAR EMPIRICAL ENVELOPE (frame_composer arg 21,
#     "t_us:w,..." nodes, rejection-sampled): nodes derived from the
#     real/flat-sim arrivals ratio with the trigger bump (x1.119) divided out,
#     then iterated twice against the window gate at 100-frame flat-rate cal.
#   - Triggered collision (MBD-fired-weighted, t_c=0, trig_n=1.0) and
#     RSPEC x0.88 unchanged from the first v3.6 pass. rate_tau=0 (off).
# Readout: FULL revC parameter set (verbatim from v36_readout_fix.sh - the
# first v36 script had inherited the pre-revB line; stale-params incident #2).
# Env nodes passed as $1 (locked by v36_envcal iterations).
set -euo pipefail
cd /home/rog/sPHENIX/3D_ClusterFindingML/island_post
ENV="${1:?need env_spec nodes t_us:w,...}"
DM=/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/TPC_DEADCHANNELMAP/ff/c3/ffc3f6498934c5a8ba31065292c6ebcc_TPCDeadMap_79471.root
M=/home/rog/sPHENIX/3D_ClusterFindingML/macros-offline/detectors/sPHENIX
LIBS="raw_lib_pau.root,raw_lib_pauLa.root,raw_lib_pauLb.root,raw_lib_pauLc.root,raw_lib_pauLd.root"
FLASH="raw_lib_cmflash_w.root"
RSPEC="136:1,231:1,239:1,242:1,242:1,243:1,244:1,245:1,247:1,249:1"
MBD="mbd_weights.txt|pau200.dat,pau_chunk_a.dat,pau_chunk_b.dat,pau_chunk_c.dat,pau_chunk_d.dat"
SPEC="0.008:1,0.009:1,0.011:1,0.012:1,0.013:1,0.014:1,0.018:1,0.021:1,0.027:1,0.037:1,2.2:0.25"
EVALS="$M/pau_a_eval_g4svtx_eval.root,$M/pauL_a_eval_g4svtx_eval.root,$M/pauL_b_eval_g4svtx_eval.root,$M/pauL_c_eval_g4svtx_eval.root,$M/pauL_d_eval_g4svtx_eval.root"
NB=5; PER=50
for i in $(seq 0 $((NB-1))); do
  root -l -b -q -e "
  gROOT->ProcessLine(\".L frame_composer.C+\");
  frame_composer(\"$LIBS\",\"fP_$i.root\",$PER,275.,2026091$i,$((i*PER)),\"$EVALS\",\"$FLASH\",0.44,1.0,\"$SPEC\",1.5,2.0,1.1,0.75,\"$RSPEC\",1.0,\"$MBD\",1.0,0.0,\"$ENV\");
  gROOT->ProcessLine(\".L tpc_digitize.C+\");
  tpc_readout(\"fP_$i.root\",\"dP_$i.root\",0.93,20.0,1,1,4711,\"$DM\",11.0,0.39,0.55,1.0,0.021,7.0,36.0,70.0,11.0,940.0,2,0.29,10.0,1.24,1.06,-1.0,5.0,1.0);" 2>&1 | grep -E "frame_composer: $PER|rate envelope|tpc_readout: fP"
  root -l -b -q "islandize91.C+(\"dP_$i.root\",\"i91_$i.root\",1,\"\",\"fP_$i.root\")" 2>&1 | grep "islandize91: dP"
done
hadd -f frames_pau_production_v36.root fP_*.root > /dev/null 2>&1
hadd -f digi_frames_production_v36.root dP_*.root > /dev/null 2>&1
hadd -f island91_frames_production_v36.root i91_*.root > /dev/null 2>&1
rm -f fP_*.root dP_*.root i91_*.root
root -l -b -q -e 'gROOT->ProcessLine(".L islandize.C+"); islandize("digi_frames_production_v36.root","island_frames_v36.root",1);' 2>&1 | grep islandize
echo "=== V3.6revB PRODUCTION COMPLETE ==="
