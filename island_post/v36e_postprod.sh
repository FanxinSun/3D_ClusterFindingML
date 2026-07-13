#!/usr/bin/env bash
# v3.6revE post-production batch: acceptance -> prodclus -> figures.
set -uo pipefail
cd /home/rog/sPHENIX/3D_ClusterFindingML/island_post

echo "=== [A] acceptance battery ==="
./v36b_accept.sh 2>&1 | grep -E "WINDOWS|BUMP|STEP|B42METRIC|B43REGION|ISLANDS"

echo "=== [B] production phi-spectrum + narrow-bright ==="
root -l -b -q -e '
auto sp=[](const char*fn,const char*tag){
  TFile*f=TFile::Open(fn); TTree*t=(TTree*)f->Get("island");
  TH1D h("h","",8,0.5,8.5); t->Draw("phisize>>h","","goff"); h.Scale(1./h.Integral());
  double nb=t->GetEntries("phisize==1&&adc>=80")/(double)t->GetEntries();
  printf("PHISPEC %-4s",tag); for(int b=1;b<=8;++b) printf(" %.4f",h.GetBinContent(b));
  printf(" | nbright %.5f\n", nb); };
sp("island_frames_v36.root","sim"); sp("island_real.root","real");' 2>&1 | grep PHISPEC

echo "=== [C] prodclus v36 (ML baseline) ==="
root -l -b -q -e 'gROOT->ProcessLine(".L prodclus.C+");
prodclus("digi_frames_production_v36.root","prodclus_v36.root",1,5.0,3.0,5.0,0,10,20,1,"","frames_pau_production_v36.root");' 2>&1 | grep prodclus:

echo "=== [D] figure suite replot ==="
# hits_profile.C holds TWO entry points: hits_profile() = exam5/AuAu-era record
# (same output filename!), hits_profile_pau() = living pAu figure. Call explicitly.
root -l -b -q -e 'gROOT->LoadMacro("hits_profile.C"); hits_profile_pau();' 2>&1 | grep -E "saved|Error in"
for m in cmp_pau p3_cmp residuals_evidence zs_shapes zsize_physics pau_day3 pau_final z_bridge cluster_gallery cmp_prodclus; do
  root -l -b -q "$m.C" 2>&1 | grep -E "saved|Error in" | head -3
done
root -l -b -q "asym_showcase.C" 2>&1 | grep -E "saved|Error in" | head -3
root -l -b -q "funny_shapes.C" 2>&1 | grep -E "saved|Error in" | head -2
cp -f /home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/funny_shapes.png \
      /home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/funny_shapes_v36.png

echo "=== [E] bestmatch scan (criterial, v36) ==="
root -l -b -q -e 'gROOT->ProcessLine(".L bestmatch_scan.C+"); bestmatch_scan();' 2>&1 | grep -E "BEST|best|score" | tail -8

echo "=== V3.6revE POSTPROD COMPLETE ==="
