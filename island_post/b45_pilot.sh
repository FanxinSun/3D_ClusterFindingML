#!/bin/bash
# v3.5 pilot: 20 composed frames with CLOUD-0.106 libs + refit flash weights.
# Gates: flash excess-over-baseline per region (real: R1 0.00298 R2 0.00507 R3 0.00532),
# frames-level <phisize> (real 3.40), island R1 share (real 0.359).
set -e
cd "$(dirname "$0")"
DM=/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/TPC_DEADCHANNELMAP/ff/c3/ffc3f6498934c5a8ba31065292c6ebcc_TPCDeadMap_79471.root
M=/home/rog/sPHENIX/3D_ClusterFindingML/macros-offline/detectors/sPHENIX
LIBS="raw_lib_pau.root,raw_lib_pauLa.root,raw_lib_pauLb.root,raw_lib_pauLc.root,raw_lib_pauLd.root"
EVALS="$M/pau_a_eval_g4svtx_eval.root,$M/pauL_a_eval_g4svtx_eval.root,$M/pauL_b_eval_g4svtx_eval.root,$M/pauL_c_eval_g4svtx_eval.root,$M/pauL_d_eval_g4svtx_eval.root"
SPEC="0.008:1,0.009:1,0.011:1,0.012:1,0.013:1,0.014:1,0.018:1,0.021:1,0.027:1,0.037:1,2.2:0.25"
RSPEC="154:1,263:1,272:1,275:1,275:1,276:1,277:1,278:1,281:1,283:1"
root -l -b -q -e "
gROOT->ProcessLine(\".L frame_composer.C+\");
frame_composer(\"$LIBS\",\"fPILOT.root\",20,275.,20260910,0,\"$EVALS\",\"raw_lib_cmflash_w.root\",0.44,1.0,\"$SPEC\",1.5,2.0,1.1,0.75,\"$RSPEC\",1.0);
gROOT->ProcessLine(\".L tpc_digitize.C+\");
tpc_readout(\"fPILOT.root\",\"dPILOT.root\",0.93,20.0,1,1,4711,\"$DM\",11.0,0.39,0.55,1.0,0.021,7.0,36.0,70.0,11.0,940.0,2,0.26,10.0,1.24,1.06,-1.0,5.0,1.0);" 2>&1 | grep -E "frame_composer: 20|pixels kept"
[ dPILOT.root -nt fPILOT.root ] || { echo "GUARD: dPILOT stale — readout failed"; exit 1; }
root -l -b -q -e "gROOT->ProcessLine(\".L islandize.C+\"); islandize(\"dPILOT.root\",\"iPILOT.root\",1);" 2>&1 | tail -1
root -l -b -q -e '
// flash excess per region + island stats
TFile*fd=TFile::Open("dPILOT.root"); TTree*t=(TTree*)fd->Get("ntp_hit");
float lay,adc,tb; t->SetBranchStatus("*",0);
for(auto b:{"layer","adc","tbin"}) t->SetBranchStatus(b,1);
t->SetBranchAddress("layer",&lay); t->SetBranchAddress("adc",&adc); t->SetBranchAddress("tbin",&tb);
double win[3]={0,0,0}, side[3]={0,0,0}, tot[3]={0,0,0};
for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
  int rg=lay<23?0:(lay<39?1:2); tot[rg]+=1;
  if(tb>=322&&tb<=340) win[rg]+=1;
  else if((tb>=295&&tb<=313)||(tb>=350&&tb<=368)) side[rg]+=1; }
printf("B45PILOT flash excess: R1 %.5f (tgt .00298) R2 %.5f (tgt .00507) R3 %.5f (tgt .00532)\n",
  (win[0]-side[0]/2)/tot[0], (win[1]-side[1]/2)/tot[1], (win[2]-side[2]/2)/tot[2]);
TFile*fi=TFile::Open("iPILOT.root"); TTree*ti=(TTree*)fi->Get("island");
TH1D hp("hp","",60,0.5,60.5); ti->Draw("phisize>>hp","","goff");
TH1D hz("hz","",200,0.5,200.5); ti->Draw("zsize>>hz","","goff");
double n=hz.Integral(), r1=ti->GetEntries("layer<23");
printf("B45PILOT islands: /fr %.0f (tgt 23.9k) | <phisize> %.2f (tgt 3.40) | <zsize> %.2f (tgt 4.27) | R1 share %.3f (tgt 0.359)\n",
  n/20., hp.GetMean(), hz.GetMean(), r1/n);
double s1r1=ti->GetEntries("size==1&&layer<23"), s1r2=ti->GetEntries("size==1&&layer>=23&&layer<39"), s1r3=ti->GetEntries("size==1&&layer>=39");
printf("B45PILOT singles/fr: R1 %.0f (tgt 789) R2 %.0f (tgt 588) R3 %.0f (tgt 1143) | frac1px %.3f (tgt 0.105)\n",
  s1r1/20., s1r2/20., s1r3/20., (s1r1+s1r2+s1r3)/n);
double m2r1=ti->GetEntries("adc<80&&size==2&&layer<23")/20., m2r2=ti->GetEntries("adc<80&&size==2&&layer>=23&&layer<39")/20., m2r3=ti->GetEntries("adc<80&&size==2&&layer>=39")/20.;
double m34r1=ti->GetEntries("adc<80&&size>=3&&size<=4&&layer<23")/20., m34r2=ti->GetEntries("adc<80&&size>=3&&size<=4&&layer>=23&&layer<39")/20., m34r3=ti->GetEntries("adc<80&&size>=3&&size<=4&&layer>=39")/20.;
printf("B45PILOT lowadc 2px/fr: R1 %.0f (tgt 407) R2 %.0f (tgt 346) R3 %.0f (tgt 377) | 3-4px: R1 %.0f (tgt 237) R2 %.0f (tgt 37) R3 %.0f (tgt 38)\n", m2r1,m2r2,m2r3,m34r1,m34r2,m34r3);
TH1D hmx("hmx","",60,930.5,990.5); ti->Draw("maxadc>>hmx","maxadc>930&&maxadc<991","goff");
printf("B45PILOT maxadc bump: mean %.1f rms %.2f (real ~948.8 / ~5.8)\n", hmx.GetMean(), hmx.GetRMS());' 2>&1 | grep B45PILOT
echo "pilot complete"
