#!/usr/bin/env bash
# Leg (b) lever pilot: transport cloud sigma 0.12 (current) vs 0.04 (collaboration
# gem_cloud_sigma, PHG4TpcPadPlaneReadout). Identical pilot config otherwise:
# single lib (pau), 100 flat frames @275, iso on, flash/trig/env off, revC readout.
# Measures PHISPEC (phi2 target) + b42 anchors for the Pareto judgment.
set -euo pipefail
cd /home/rog/sPHENIX/3D_ClusterFindingML/island_post
M=/home/rog/sPHENIX/3D_ClusterFindingML/macros-offline/detectors/sPHENIX
DM=/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/TPC_DEADCHANNELMAP/ff/c3/ffc3f6498934c5a8ba31065292c6ebcc_TPCDeadMap_79471.root
for tag in ${TAGS:-c12 c04}; do
  LIB=raw_lib_pau_$tag.root; [ $tag = c12 ] && LIB=raw_lib_pau.root
  root -l -b -q -e "
  gROOT->ProcessLine(\".L frame_composer.C+\");
  frame_composer(\"$LIB\",\"fCL_$tag.root\",100,275.,31415,0,\"$M/pau_a_eval_g4svtx_eval.root\",\"\",0.0,1.0,\"\",1.5,2.0,1.1,0.75,\"\",1.0,\"\",0.0,0.0,\"\");
  gROOT->ProcessLine(\".L tpc_digitize.C+\");
  tpc_readout(\"fCL_$tag.root\",\"dCL_$tag.root\",0.93,20.0,1,1,4711,\"$DM\",11.0,0.39,0.55,1.0,0.021,7.0,36.0,70.0,11.0,940.0,2,0.29,10.0,1.24,1.06,-1.0,5.0,1.0);" 2>&1 | grep -E "frame_composer: 100|tpc_readout: fCL"
  root -l -b -q -e "gROOT->ProcessLine(\".L islandize.C+\"); islandize(\"dCL_$tag.root\",\"iCL_$tag.root\",1);" 2>&1 | grep islandize
  root -l -b -q -e "gROOT->ProcessLine(\".L b42_metrics.C+\"); b42_metrics(\"dCL_$tag.root\",\"$tag\");" 2>&1 | grep -E "B42METRIC|B43REGION"
done
root -l -b -q -e '
auto sp=[](const char*fn,const char*tag){
  TFile*f=TFile::Open(fn); TTree*t=(TTree*)f->Get("island");
  TH1D h("h","",8,0.5,8.5); t->Draw("phisize>>h","","goff"); h.Scale(1./h.Integral());
  double n=t->GetEntries();
  int e0=(int)t->GetMinimum("event"), e1=(int)t->GetMaximum("event");
  printf("PHISPEC %-4s",tag); for(int b=1;b<=8;++b) printf(" %.4f",h.GetBinContent(b));
  printf(" | isl/fr %.0f\n", n/(e1-e0+1)); };
for (const char*tg : {${PLIST:-"c12","c04"}}) sp(Form("iCL_%s.root",tg),tg); sp("island_real.root","real");' 2>&1 | grep PHISPEC
echo "=== CLOUD PILOT DONE ==="
