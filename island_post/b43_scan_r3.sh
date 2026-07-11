#!/bin/bash
# b43 round 3: re-center z-tails (p_trig) and phi-size (sigma_pad) under regional ZS.
# Island targets: <phisize> 3.40, zsize P13/P20/P30 = 2.85/1.82/1.22e-2, isl/fr 23.9k.
set -e
cd "$(dirname "$0")"
DM=/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/TPC_DEADCHANNELMAP/ff/c3/ffc3f6498934c5a8ba31065292c6ebcc_TPCDeadMap_79471.root
OUT=b43_scan_results.txt
for cfg in "0.30 0.55" "0.26 0.55" "0.30 0.50" "0.26 0.50" "0.30 0.45" "0.26 0.45"; do
  set -- $cfg
  TAG="R6_pt${1}_sp${2}"
  root -l -b -q -e "gROOT->ProcessLine(\".L tpc_digitize.C+\"); tpc_readout(\"raw_b42_scan.root\",\"b43_tmp.root\",0.93,20.0,1,1,4711,\"$DM\",11.0,0.39,$2,0.00007,0.021,7.0,36.0,70.0,11.0,940.0,2,$1,10.0,1.24,1.06);" 2>&1 | grep "pixels kept"
  root -l -b -q -e "gROOT->ProcessLine(\".L b42_metrics.C+\"); b42_metrics(\"b43_tmp.root\",\"$TAG\");" 2>&1 | grep -E "B42METRIC|B43REGION" | tee -a $OUT
  root -l -b -q -e "gROOT->ProcessLine(\".L islandize.C+\"); islandize(\"b43_tmp.root\",\"b43_isl_tmp.root\",1);" > /dev/null 2>&1
  root -l -b -q -e '
  TFile*f=TFile::Open("b43_isl_tmp.root"); TTree*t=(TTree*)f->Get("island");
  TH1D h("h","",200,0.5,200.5); t->Draw("zsize>>h","","goff");
  TH1D hp("hp","",60,0.5,60.5); t->Draw("phisize>>hp","","goff");
  double tot=h.Integral();
  printf("B43ISLAND '$TAG' isl/fr %.0f | <zsize> %.2f | P13 %.2e P20 %.2e P30 %.2e | <phisize> %.2f\n",
    tot/40., h.GetMean(), h.Integral(h.FindBin(13),200)/tot, h.Integral(h.FindBin(20),200)/tot,
    h.Integral(h.FindBin(30),200)/tot, hp.GetMean());' 2>&1 | grep B43ISLAND | tee -a $OUT
done
rm -f b43_tmp.root b43_isl_tmp.root
echo "b43 round 3 complete"
