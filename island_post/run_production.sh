#!/usr/bin/env bash
# P3 production v3.2 (B4.1: joint Pareto refit — gain 0.87, sigma_pad 0.55, tail 0.021/7,
# saturation re-clamp; cloud 0.12 unchanged)
# was v3.1 (adds B4 ion tail f=0.025 tau=7, rate re-trim) (P3.2 charge-spread CLOUD=0.12, gain 0.74, per-frame rate jitter
# ~275 kHz empirical spec, full CM flash) — see PIPELINE.md v3 calibration.
# Superseded header: v2.2 (flash realism: per-stripe+speckle dispersion, optical halo blur,
# region weights R1x1.25/R3x0.40, giant s=2.2, knee spec remapped): labeled p+Au frames @ 390 kHz (inelastic rate = rmbd/P_fire) from the
# full 320-collision library + CM LASER FLASH injection (collaboration PHG4TpcCentralMembrane
# model: delay 4346 ns, prob 0.44, empirical intensity spec, tbin jitter 1.5) -> cls=2 class.
set -euo pipefail
cd /home/rog/sPHENIX/3D_ClusterFindingML/island_post
DM=/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/TPC_DEADCHANNELMAP/ff/c3/ffc3f6498934c5a8ba31065292c6ebcc_TPCDeadMap_79471.root
M=/home/rog/sPHENIX/3D_ClusterFindingML/macros-offline/detectors/sPHENIX
LIBS="raw_lib_pau.root,raw_lib_pauLa.root,raw_lib_pauLb.root,raw_lib_pauLc.root,raw_lib_pauLd.root"
FLASH="raw_lib_cmflash_w.root"
RSPEC="154:1,263:1,272:1,275:1,275:1,276:1,277:1,278:1,281:1,283:1"
SPEC="0.008:1,0.009:1,0.011:1,0.012:1,0.013:1,0.014:1,0.018:1,0.021:1,0.027:1,0.037:1,2.2:0.25"
EVALS="$M/pau_a_eval_g4svtx_eval.root,$M/pauL_a_eval_g4svtx_eval.root,$M/pauL_b_eval_g4svtx_eval.root,$M/pauL_c_eval_g4svtx_eval.root,$M/pauL_d_eval_g4svtx_eval.root"
NB=5; PER=50   # 5 batches x 50 frames = 250 frames (pAu frames are light: ~1.3M raw px each)
for i in $(seq 0 $((NB-1))); do
  root -l -b -q -e "
  gROOT->ProcessLine(\".L frame_composer.C+\");
  frame_composer(\"$LIBS\",\"fP_$i.root\",$PER,275.,2026091$i,$((i*PER)),\"$EVALS\",\"$FLASH\",0.44,1.0,\"$SPEC\",1.5,2.0,1.1,0.75,\"$RSPEC\");
  gROOT->ProcessLine(\".L tpc_digitize.C+\");
  tpc_readout(\"fP_$i.root\",\"dP_$i.root\",0.87,20.0,1,1,4711,\"$DM\",11.0,0.39,0.55,0.0005,0.021,7.0);" 2>&1 | grep -E "frame_composer: $PER|tpc_readout: fP"
  root -l -b -q "islandize91.C+(\"dP_$i.root\",\"i91_$i.root\",1,\"\",\"fP_$i.root\")" 2>&1 | grep "islandize91: dP"
done
hadd -f frames_pau_production_v32.root fP_*.root > /dev/null 2>&1   # raw + sidecars (for reproducibility)
hadd -f digi_frames_production_v32.root dP_*.root > /dev/null 2>&1
hadd -f island91_frames_production_v32.root i91_*.root > /dev/null 2>&1
rm -f fP_*.root dP_*.root i91_*.root
echo "=== PRODUCTION COMPLETE ==="
