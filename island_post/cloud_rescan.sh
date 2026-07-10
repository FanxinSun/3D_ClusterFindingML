#!/usr/bin/env bash
set -uo pipefail
cd /home/rog/sPHENIX/3D_ClusterFindingML/island_post
M=/home/rog/sPHENIX/3D_ClusterFindingML/macros-offline/detectors/sPHENIX
for CL in 0.08 0.10; do
  sed -i "s|const double CLOUD = [0-9.]*;.*|const double CLOUD = $CL;|" tpc_digitize.C
  root -l -b -q -e "gROOT->ProcessLine(\".L tpc_digitize.C+\");
  tpc_transport(\"$M/pauL_a_eval_g4svtx_eval.root\",\"raw_scan_cl$CL.root\",75);" 2>&1 | grep transport:
done
sed -i "s|const double CLOUD = [0-9.]*;|const double CLOUD = 0.12;|" tpc_digitize.C
echo "CLOUD RESCAN TRANSPORTS DONE"
