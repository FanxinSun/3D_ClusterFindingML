#!/usr/bin/env bash
# P3.2 charge-spread scan: CLOUD sigma in {0.06(current), 0.10, 0.14}
# one 75-collision lib re-transported per value -> 30 frames @390 -> B3 -> islands -> targets
set -uo pipefail
cd /home/rog/sPHENIX/3D_ClusterFindingML/island_post
DM=/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/TPC_DEADCHANNELMAP/ff/c3/ffc3f6498934c5a8ba31065292c6ebcc_TPCDeadMap_79471.root
EV=/home/rog/sPHENIX/3D_ClusterFindingML/macros-offline/detectors/sPHENIX/pauL_a_eval_g4svtx_eval.root
for CL in 0.06 0.10 0.14; do
  sed -i "s|const double CLOUD = [0-9.]*;|const double CLOUD = $CL;|" tpc_digitize.C
  root -l -b -q -e "
  gROOT->ProcessLine(\".L tpc_digitize.C+\");
  tpc_transport(\"$EV\",\"raw_scan_cl.root\",75);
  gROOT->ProcessLine(\".L frame_composer.C+\");
  frame_composer(\"raw_scan_cl.root\",\"f_scan.root\",30,390.,20260906,0,\"\");
  tpc_readout(\"f_scan.root\",\"d_scan.root\",0.70,20.0,1,1,4711,\"$DM\",11.0,0.39,0.0,0.0005);" > /dev/null 2>&1
  root -l -b -q 'islandize.C("d_scan.root","i_scan.root",1)' > /dev/null 2>&1
  root -l -b -q -e "
  TFile*f=TFile::Open(\"i_scan.root\"); TTree*t=(TTree*)f->Get(\"island\");
  TFile*fd=TFile::Open(\"d_scan.root\"); TTree*h=(TTree*)fd->Get(\"ntp_hit\");
  double px=h->GetEntries()/30., n=t->GetEntries()/30.;
  TH1D hx(\"hx\",\"\",121,-0.5,120.5); h->Draw(\"adc>>hx\",\"\",\"goff\");
  TH1D hm(\"hm\",\"\",120,0,1080); t->Draw(\"maxadc>>hm\",\"\",\"goff\");
  double m2060=hm.Integral(hm.FindBin(20),hm.FindBin(60))/hm.Integral();
  TH1D hs(\"hs\",\"\",300,0.5,300.5); t->Draw(\"size>>hs\",\"\",\"goff\");
  TH1D hp(\"hp\",\"\",30,0.5,30.5); t->Draw(\"phisize>>hp\",\"\",\"goff\");
  printf(\"CLOUD $CL: px/fr %.0fk pixmean %.1f | n/fr %.1fk | size %.2f phisize %.2f | maxadc frac20-60 %.3f\n\",
    px/1e3, hx.GetMean(), n/1e3, hs.GetMean(), hp.GetMean(), m2060);" 2>&1 | grep CLOUD
done
sed -i "s|const double CLOUD = [0-9.]*;|const double CLOUD = 0.06;|" tpc_digitize.C
echo "targets (real):   pixmean 99.15 | n/px-matched ratio 1 | size 10.7 phisize 3.41 | maxadc frac20-60 = $(root -l -b -q -e '
TFile*f=TFile::Open("island_real.root"); TTree*t=(TTree*)f->Get("island");
TH1D hm("hm","",120,0,1080); t->Draw("maxadc>>hm","","goff");
printf("%.3f", hm.Integral(hm.FindBin(20),hm.FindBin(60))/hm.Integral());' 2>/dev/null | tail -1)"
echo "SCAN DONE"
