#!/usr/bin/env bash
# P3 production v3.5 (2026-07-12): composition-level realism pass on top of B4.3.
#   - CLOUD stays 0.12: the phi +4.9% "residual" was COMPOSITION (real's 10.5%
#     single-pixel islands dilute its mean phi; sim multi-px phi already matched
#     3.67 vs 3.68). CLOUD 0.106 was tried and REVERTED (fragmented merges, +35%
#     multi-px islands).
#   - NEW: isolated-hit background injection in frame_composer (iso_scale=1):
#     real singles 789/588/1143 per frame, adc ~ thr + Exp(4.5/6.5/7.0), 83%
#     diffuse + 17% on 100 hash-static hot pads; src 0xFE -> cls=2 noise labels.
#     Injected rates 974/460/1491 (pilot-fit effective-yield corrections).
#   - Flash region weights 0.36/0.87/0.43 (flash-only giant-yield inversion under
#     regional ZS; precision refit = dedicated laser_assess pass, declared).
#   - Readout = v3.4 winner with p_trig 0.26 -> 0.29 (singles dilute island-tail
#     denominators by ~x1.10).
# NOTE: fresh composition realization — iso draws shift the composer rng stream,
# so frames differ event-by-event from v3.2/v3.3/v3.4 (same seeds, same ensemble).
set -euo pipefail
cd /home/rog/sPHENIX/3D_ClusterFindingML/island_post
DM=/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/TPC_DEADCHANNELMAP/ff/c3/ffc3f6498934c5a8ba31065292c6ebcc_TPCDeadMap_79471.root
M=/home/rog/sPHENIX/3D_ClusterFindingML/macros-offline/detectors/sPHENIX
LIBS="raw_lib_pau.root,raw_lib_pauLa.root,raw_lib_pauLb.root,raw_lib_pauLc.root,raw_lib_pauLd.root"
FLASH="raw_lib_cmflash_w.root"
RSPEC="154:1,263:1,272:1,275:1,275:1,276:1,277:1,278:1,281:1,283:1"
SPEC="0.008:1,0.009:1,0.011:1,0.012:1,0.013:1,0.014:1,0.018:1,0.021:1,0.027:1,0.037:1,2.2:0.25"
EVALS="$M/pau_a_eval_g4svtx_eval.root,$M/pauL_a_eval_g4svtx_eval.root,$M/pauL_b_eval_g4svtx_eval.root,$M/pauL_c_eval_g4svtx_eval.root,$M/pauL_d_eval_g4svtx_eval.root"
NB=5; PER=50
for i in $(seq 0 $((NB-1))); do
  root -l -b -q -e "
  gROOT->ProcessLine(\".L frame_composer.C+\");
  frame_composer(\"$LIBS\",\"fP_$i.root\",$PER,275.,2026091$i,$((i*PER)),\"$EVALS\",\"$FLASH\",0.44,1.0,\"$SPEC\",1.5,2.0,1.1,0.75,\"$RSPEC\",1.0);
  gROOT->ProcessLine(\".L tpc_digitize.C+\");
  tpc_readout(\"fP_$i.root\",\"dP_$i.root\",0.93,20.0,1,1,4711,\"$DM\",11.0,0.39,0.55,0.00007,0.021,7.0,36.0,70.0,11.0,940.0,2,0.29,10.0,1.24,1.06);" 2>&1 | grep -E "frame_composer: $PER|tpc_readout: fP"
  root -l -b -q "islandize91.C+(\"dP_$i.root\",\"i91_$i.root\",1,\"\",\"fP_$i.root\")" 2>&1 | grep "islandize91: dP"
done
hadd -f frames_pau_production_v35.root fP_*.root > /dev/null 2>&1
hadd -f digi_frames_production_v35.root dP_*.root > /dev/null 2>&1
hadd -f island91_frames_production_v35.root i91_*.root > /dev/null 2>&1
rm -f fP_*.root dP_*.root i91_*.root
root -l -b -q -e 'gROOT->ProcessLine(".L islandize.C+"); islandize("digi_frames_production_v35.root","island_frames_v35.root",1);' 2>&1 | grep islandize
echo "=== V3.5 PRODUCTION COMPLETE ==="
