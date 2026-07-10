#!/usr/bin/env bash
set -uo pipefail
cd /home/rog/sPHENIX/3D_ClusterFindingML/island_post
DM=/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/TPC_DEADCHANNELMAP/ff/c3/ffc3f6498934c5a8ba31065292c6ebcc_TPCDeadMap_79471.root
for CL in 0.10 0.12; do
  LIB="raw_scan_cl$CL.root"; [ "$CL" = "0.12" ] && LIB="raw_lib_pauLa.root"
  root -l -b -q -e "gROOT->ProcessLine(\".L frame_composer.C+\");
  frame_composer(\"$LIB\",\"f_j.root\",30,275.,20260911,0,\"\");" > /dev/null 2>&1
  for GS in "0.78 0.35" "0.84 0.50" "0.90 0.60"; do
    set -- $GS
    root -l -b -q -e "gROOT->ProcessLine(\".L tpc_digitize.C+\");
    tpc_readout(\"f_j.root\",\"d_j.root\",$1,20.0,1,1,4711,\"$DM\",11.0,0.39,$2,0.0005,0.025,7.0);" > /dev/null 2>&1
    root -l -b -q 'islandize.C("d_j.root","i_j.root",1)' > /dev/null 2>&1
    root -l -b -q -e "
    TFile*fd=TFile::Open(\"d_j.root\"); TTree*h=(TTree*)fd->Get(\"ntp_hit\");
    TFile*fi=TFile::Open(\"i_j.root\"); TTree*t=(TTree*)fi->Get(\"island\");
    TH1D hx(\"hx\",\"\",1101,-0.5,1100.5); h->Draw(\"adc>>hx\",\"\",\"goff\");
    TH1D hz(\"hz\",\"\",121,-0.5,120.5); h->Draw(\"adc>>hz\",\"\",\"goff\");
    double band=hz.Integral(hz.FindBin(11),hz.FindBin(19))/hz.Integral(hz.FindBin(20),hz.FindBin(40));
    TH1D hm(\"hm\",\"\",120,0,1080); t->Draw(\"maxadc>>hm\",\"\",\"goff\");
    double m2060=hm.Integral(hm.FindBin(20),hm.FindBin(60))/hm.Integral();
    TProfile pz(\"pz\",\"\",20,0,6000); t->Draw(\"zsize:adc>>pz\",\"\",\"goff\");
    TH1D hs(\"hs\",\"\",300,0.5,300.5); t->Draw(\"size>>hs\",\"\",\"goff\");
    TH1D hp(\"hp\",\"\",30,0.5,30.5); t->Draw(\"phisize>>hp\",\"\",\"goff\");
    double mid=hx.Integral(hx.FindBin(200),hx.FindBin(800))/hx.Integral();
    double sat=hx.Integral(hx.FindBin(930),hx.FindBin(960))/hx.Integral();
    printf(\"cl=$CL g=$1 sp=$2: pixm %.1f mid %.4f sat %.4f band %.3f | phi %.2f cliff %.3f | zg %.2f/%.2f | size %.2f px %.0fk\n\",
      hx.GetMean(), mid, sat, band, hp.GetMean(), m2060,
      pz.GetBinContent(pz.FindBin(3000)), pz.GetBinContent(pz.FindBin(5500)),
      hs.GetMean(), h->GetEntries()/30./1e3);" 2>&1 | grep "cl="
  done
done
echo "TARGETS: pixm 99.15 mid 0.0850 sat 0.0140 band 0.233 | phi 3.41 cliff 0.251 | zg 6.97/10.35 | size 10.7"
echo "SP GRID DONE"
