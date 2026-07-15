#!/usr/bin/env bash
# v3.6revB acceptance battery: windows + trigger bump + pixel metrics + islands.
set -euo pipefail
cd /home/rog/sPHENIX/3D_ClusterFindingML/island_post
REAL=/home/rog/sPHENIX/3D_ClusterFindingML/clusters_seeds_island_79507-0.root_ntuplizer.root

echo "--- [1] arrivals: window fractions + trigger bump (production, RSPEC mix) ---"
root -l -b -q -e "
TFile*fr=TFile::Open(\"$REAL\"); TTree*tr=(TTree*)fr->Get(\"ntp_hit\");
TFile*fs=TFile::Open(\"digi_frames_production_v40b.root\"); TTree*ts=(TTree*)fs->Get(\"ntp_hit\");
TH1D hr(\"hr\",\"\",971,-0.5,970.5), hs(\"hs\",\"\",971,-0.5,970.5);
tr->Draw(\"tbin>>hr\",\"layer>=7&&layer<=54&&adc>0\",\"goff\"); ts->Draw(\"tbin>>hs\",\"\",\"goff\");
auto wf=[](TH1D&h,const char*t){ double tot=0,w1=0,w2=0,w3=0;
  for(int b=60;b<=950;++b){ if(b>=318&&b<=345) continue; double v=h.GetBinContent(b+1); tot+=v;
    if(b<=240) w1+=v; else if(b>=270&&b<=600) w2+=v; else if(b>=650) w3+=v; }
  printf(\"WINDOWS %-4s %.4f %.4f %.4f\n\", t, w1/tot, w2/tot, w3/tot); };
wf(hr,\"real\"); wf(hs,\"sim\");
auto bump=[](TH1D&h,const char*t){ double a=0,b=0;
  for(int x=60;x<=240;++x) a+=h.GetBinContent(x+1);
  for(int x=270;x<=450;++x) b+=h.GetBinContent(x+1);
  printf(\"BUMP %-4s pre/post %.4f\n\", t, (a/181.)/(b/181.)); };
bump(hr,\"real\"); bump(hs,\"sim\");
auto step=[](TH1D&h,const char*t){ double a=0,b=0;
  for(int x=230;x<=246;++x) a+=h.GetBinContent(x+1);
  for(int x=254;x<=270;++x) b+=h.GetBinContent(x+1);
  printf(\"STEP %-4s drift-edge %.4f\n\", t, a/b); };
step(hr,\"real\"); step(hs,\"sim\");" 2>&1 | grep -E "WINDOWS|BUMP|STEP"

echo "--- [2] pixel metrics (b42_metrics: shares, sub10, pixmean, runs) ---"
root -l -b -q -e "gROOT->ProcessLine(\".L b42_metrics.C+\"); b42_metrics(\"digi_frames_production_v40b.root\",\"v36revB\");" 2>&1 | grep -E "B42METRIC|B43REGION"
root -l -b -q -e "gROOT->ProcessLine(\".L b42_metrics.C+\"); b42_metrics(\"$REAL\",\"real\");" 2>&1 | grep -E "B42METRIC|B43REGION"

echo "--- [3] islands/frame + sizes ---"
root -l -b -q -e "
TFile*fr=TFile::Open(\"island_real.root\"); TTree*ir=(TTree*)fr->Get(\"island\");
TFile*fs=TFile::Open(\"island_frames_v40b.root\"); TTree*is=(TTree*)fs->Get(\"island\");
auto st=[](TTree*t,const char*tag){ double n=t->GetEntries();
  int evmin=(int)t->GetMinimum(\"event\"), evmax=(int)t->GetMaximum(\"event\");
  double nev=evmax-evmin+1;
  TH1D h1(\"h1\",\"\",200,0,200), h2(\"h2\",\"\",50,0,50);
  t->Draw(\"size>>h1\",\"\",\"goff\"); t->Draw(\"phisize>>h2\",\"\",\"goff\");
  printf(\"ISLANDS %-4s per-frame %.0f | size mean %.2f | phisize mean %.3f\n\",
    tag, n/nev, h1.GetMean(), h2.GetMean()); };
st(ir,\"real\"); st(is,\"sim\");" 2>&1 | grep ISLANDS
echo "=== ACCEPT BATTERY DONE ==="
