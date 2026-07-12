#!/bin/bash
# v3.5 step 1: CLOUD bracket under the corrected regional ZS (B4.3 winner readout).
# Single-collision proxy (pauL_a, 25 events) — the same proxy class that set
# CLOUD=0.12 originally (P3.2). Gate on the RELATIVE phi-size shift: frames give
# <phisize> 3.57 at CLOUD 0.12 vs real 3.40 -> want ~-5%.
set -e
cd "$(dirname "$0")"
DM=/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/TPC_DEADCHANNELMAP/ff/c3/ffc3f6498934c5a8ba31065292c6ebcc_TPCDeadMap_79471.root
M=/home/rog/sPHENIX/3D_ClusterFindingML/macros-offline/detectors/sPHENIX
OUT=b45_cloud_results.txt
: > $OUT
for CL in 0.106; do
  root -l -b -q -e "gROOT->ProcessLine(\".L tpc_digitize.C+\"); tpc_transport(\"$M/pauL_a_eval_g4svtx_eval.root\",\"b45_raw.root\",25,$CL);" 2>&1 | grep "tpc_transport:"
  root -l -b -q -e "gROOT->ProcessLine(\".L tpc_digitize.C+\"); tpc_readout(\"b45_raw.root\",\"b45_digi.root\",0.93,20.0,1,1,4711,\"$DM\",11.0,0.39,0.55,0.00007,0.021,7.0,36.0,70.0,11.0,940.0,2,0.26,10.0,1.24,1.06);" 2>&1 | grep "pixels kept"
  root -l -b -q -e "gROOT->ProcessLine(\".L islandize.C+\"); islandize(\"b45_digi.root\",\"b45_isl.root\",1);" > /dev/null 2>&1
  root -l -b -q -e '
  TFile*f=TFile::Open("b45_isl.root"); TTree*t=(TTree*)f->Get("island");
  TH1D hp("hp","",60,0.5,60.5); t->Draw("phisize>>hp","","goff");
  TH1D hz("hz","",200,0.5,200.5); t->Draw("zsize>>hz","","goff");
  double r1=t->GetEntries("layer<23");
  printf("B45CLOUD cloud='$CL' <phisize> %.3f | <zsize> %.2f | islands %.0f | R1 share %.3f\n",
    hp.GetMean(), hz.GetMean(), hz.Integral(), r1/hz.Integral());' 2>&1 | grep B45CLOUD | tee -a $OUT
done
rm -f b45_raw.root b45_digi.root b45_isl.root
echo "cloud bracket complete"
