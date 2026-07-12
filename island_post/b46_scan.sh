#!/bin/bash
# b46: near-threshold refinement — shaped trace (p2 scale=1) + R1 band lever
# (gR1 x sigma_pad_R1). Targets: sh/hi 0.233, shares 0.362/0.350/0.287,
# pixm 86.1/105.0/108.5, sub10 2.27e-4, run tails per B4.2.
set -e
cd "$(dirname "$0")"
DM=/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/TPC_DEADCHANNELMAP/ff/c3/ffc3f6498934c5a8ba31065292c6ebcc_TPCDeadMap_79471.root
OUT=b46_results.txt
for cfg in "1.24 0.40" "1.24 0.30" "1.30 0.40"; do
  set -- $cfg
  TAG="B46_g${1}_sp${2}"
  root -l -b -q -e "gROOT->ProcessLine(\".L tpc_digitize.C+\"); tpc_readout(\"raw_b45_scan.root\",\"b46_tmp.root\",0.93,20.0,1,1,4711,\"$DM\",11.0,0.39,0.55,1.0,0.021,7.0,36.0,70.0,11.0,940.0,2,0.29,10.0,$1,1.06,$2);" 2>&1 | grep "pixels kept"
  root -l -b -q -e "gROOT->ProcessLine(\".L b42_metrics.C+\"); b42_metrics(\"b46_tmp.root\",\"$TAG\");" 2>&1 | grep -E "B42METRIC|B43REGION" | tee -a $OUT
done
rm -f b46_tmp.root
echo "b46 scan complete"
