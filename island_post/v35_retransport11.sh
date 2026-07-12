#!/usr/bin/env bash
set -uo pipefail
# v3.5 retransport (2026-07-12): CLOUD 0.106 (recalibrated under B4.3 regional ZS;
# phi-size proxy 3.435 vs target 3.43) + flash region weights REFIT for thr-11 R1
# (excess-over-baseline real/v34: x0.24/x0.55/x0.67 on the old 1.25/1.00/0.40 ->
# 0.298/0.547/0.268, first pass; pilot-composed before full production).
cd /home/rog/sPHENIX/3D_ClusterFindingML/island_post
M=/home/rog/sPHENIX/3D_ClusterFindingML/macros-offline/detectors/sPHENIX
root -l -b -q -e 'gROOT->ProcessLine(".L tpc_digitize.C+");
tpc_transport("'$M'/pau_a_eval_g4svtx_eval.root","raw_lib11_pau.root",20,0.11);
tpc_transport("'$M'/pauL_a_eval_g4svtx_eval.root","raw_lib11_pauLa.root",75,0.11);
tpc_transport("'$M'/pauL_b_eval_g4svtx_eval.root","raw_lib11_pauLb.root",75,0.11);
tpc_transport("'$M'/pauL_c_eval_g4svtx_eval.root","raw_lib11_pauLc.root",75,0.11);
tpc_transport("'$M'/pauL_d_eval_g4svtx_eval.root","raw_lib11_pauLd.root",75,0.11);
tpc_transport("/home/rog/sPHENIX/3D_ClusterFindingML/laser_try/laser_cm_g4hits.root","raw_lib11_cmflash.root",1,0.11);' 2>&1 | grep transport:
echo "RETRANSPORT11 DONE"
