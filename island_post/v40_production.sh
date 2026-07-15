#!/usr/bin/env bash
# v4.0 RUN-THROUGH production (2026-07-16): SWAP+THIN library (Angantyr,
# alpha=3 thin), working point sigma(0.040,0.0070) g1.03 iso0.27, P2 refit,
# PHYSICAL rate (RSPEC ~623 kHz pileup) and PHYSICAL trigger (trig_n=1.0:
# thinned fired-mean == real 8.2k by construction). Envelope: v3.6 nodes as
# starting shape (one derivation iteration planned post-acceptance).
set -euo pipefail
cd /home/rog/sPHENIX/3D_ClusterFindingML/island_post
DM=/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/TPC_DEADCHANNELMAP/ff/c3/ffc3f6498934c5a8ba31065292c6ebcc_TPCDeadMap_79471.root
LIBS="raw_lib_v40_0.root,raw_lib_v40_1.root,raw_lib_v40_2.root,raw_lib_v40_3.root,raw_lib_v40_4.root,raw_lib_v40_5.root,raw_lib_v40_6.root,raw_lib_v40_7.root,raw_lib_v40_8.root,raw_lib_v40_9.root"
EVALS="eval_v40_0.root,eval_v40_1.root,eval_v40_2.root,eval_v40_3.root,eval_v40_4.root,eval_v40_5.root,eval_v40_6.root,eval_v40_7.root,eval_v40_8.root,eval_v40_9.root"
FLASH="raw_lib_cmflash_w.root"
RSPEC="369:1,628:1,651:1,662:1,662:1,663:1,666:1,667:1,674:1,680:1"
MBD="v40_mbd.txt|ang_run_0.dat,ang_run_1.dat,ang_run_2.dat,ang_run_3.dat,ang_run_4.dat,ang_run_5.dat,ang_run_6.dat,ang_run_7.dat,ang_run_8.dat,ang_run_9.dat"
SPEC="0.008:1,0.009:1,0.011:1,0.012:1,0.013:1,0.014:1,0.018:1,0.021:1,0.027:1,0.037:1,2.2:0.25"
ENV="$1"
NB=5; PER=50
for i in $(seq 0 $((NB-1))); do
  root -l -b -q -e "
  gROOT->ProcessLine(\".L frame_composer.C+\");
  frame_composer(\"$LIBS\",\"fP_$i.root\",$PER,600.,2026091$i,$((i*PER)),\"$EVALS\",\"$FLASH\",0.44,1.0,\"$SPEC\",1.5,2.0,1.1,0.75,\"$RSPEC\",0.27,\"$MBD\",1.0,0.0,\"$ENV\");
  gROOT->ProcessLine(\".L tpc_digitize.C+\");
  tpc_readout(\"fP_$i.root\",\"dP_$i.root\",1.03,20.0,1,1,4711,\"$DM\",11.0,0.39,0.55,0.70,0.021,7.0,36.0,70.0,11.0,940.0,2,0.29,10.0,1.24,1.06,-1.0,5.0,1.0);" 2>&1 | grep -E "frame_composer: $PER|rate envelope|tpc_readout: fP"
  root -l -b -q "islandize91.C+(\"dP_$i.root\",\"i91_$i.root\",1,\"\",\"fP_$i.root\")" 2>&1 | grep "islandize91: dP"
done
hadd -f frames_pau_production_v40.root fP_*.root > /dev/null 2>&1
hadd -f digi_frames_production_v40.root dP_*.root > /dev/null 2>&1
hadd -f island91_frames_production_v40.root i91_*.root > /dev/null 2>&1
rm -f fP_*.root dP_*.root i91_*.root
root -l -b -q -e 'gROOT->ProcessLine(".L islandize.C+"); islandize("digi_frames_production_v40.root","island_frames_v40.root",1);' 2>&1 | grep islandize
echo "=== V40 PRODUCTION COMPLETE ==="
