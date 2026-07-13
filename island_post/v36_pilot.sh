#!/bin/bash
# v3.6 pilot: frame time-structure (triggered collision + rate-decay envelope).
# Gates: arrivals bump excess (tgt x1.119 in [90,240] over plateau trend), post-250
# lnR slope (tgt ~0 after envelope), kept/fr + islands/fr guards.
set -e
cd "$(dirname "$0")"
DM=/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/TPC_DEADCHANNELMAP/ff/c3/ffc3f6498934c5a8ba31065292c6ebcc_TPCDeadMap_79471.root
M=/home/rog/sPHENIX/3D_ClusterFindingML/macros-offline/detectors/sPHENIX
LIBS="raw_lib_pau.root,raw_lib_pauLa.root,raw_lib_pauLb.root,raw_lib_pauLc.root,raw_lib_pauLd.root"
EVALS="$M/pau_a_eval_g4svtx_eval.root,$M/pauL_a_eval_g4svtx_eval.root,$M/pauL_b_eval_g4svtx_eval.root,$M/pauL_c_eval_g4svtx_eval.root,$M/pauL_d_eval_g4svtx_eval.root"
SPEC="0.008:1,0.009:1,0.011:1,0.012:1,0.013:1,0.014:1,0.018:1,0.021:1,0.027:1,0.037:1,2.2:0.25"
RSPEC="154:1,263:1,272:1,275:1,275:1,276:1,277:1,278:1,281:1,283:1"
MBD="mbd_weights.txt|pau200.dat,pau_chunk_a.dat,pau_chunk_b.dat,pau_chunk_c.dat,pau_chunk_d.dat"
TRIGN=${1:-2.0}; TAU=${2:-1754}
root -l -b -q -e "
gROOT->ProcessLine(\".L frame_composer.C+\");
frame_composer(\"$LIBS\",\"fPILOT.root\",20,275.,20260910,0,\"$EVALS\",\"raw_lib_cmflash_w.root\",0.44,1.0,\"$SPEC\",1.5,2.0,1.1,0.75,\"$RSPEC\",1.0,\"$MBD\",$TRIGN,$TAU);
gROOT->ProcessLine(\".L tpc_digitize.C+\");
tpc_readout(\"fPILOT.root\",\"dPILOT.root\",0.93,20.0,1,1,4711,\"$DM\",11.0,0.39,0.55,1.0,0.021,7.0,36.0,70.0,11.0,940.0,2,0.29,10.0,1.24,1.06,-1.0,5.0,1.0);" 2>&1 | grep -E "mbd map|frame_composer: 20|pixels kept"
root -l -b -q -e '
TFile*fr=TFile::Open("/home/rog/sPHENIX/3D_ClusterFindingML/clusters_seeds_island_79507-0.root_ntuplizer.root");
TTree*tr=(TTree*)fr->Get("ntp_hit");
TFile*fs=TFile::Open("dPILOT.root"); TTree*ts=(TTree*)fs->Get("ntp_hit");
TH1D hr("hr","",971,-0.5,970.5), hs("hs","",971,-0.5,970.5);
tr->Draw("tbin>>hr","layer>=7&&layer<=54&&adc>0","goff"); ts->Draw("tbin>>hs","","goff");
double kpf=hs.Integral()/20.; hr.Scale(1./hr.Integral()); hs.Scale(1./hs.Integral());
TGraph g; TF1 f2("f2","pol1",270,950);
for(int b=270;b<=950;++b){ if(b>=318&&b<=345) continue;
  double r=hr.GetBinContent(b+1), s=hs.GetBinContent(b+1);
  if(r>0&&s>0) g.SetPoint(g.GetN(), b, TMath::Log(r/s)); }
g.Fit(&f2,"RQ");
double bump=0; int nb=0;
for(int b=90;b<=240;++b){ double r=hr.GetBinContent(b+1), s=hs.GetBinContent(b+1);
  if(r>0&&s>0){ bump += TMath::Log(r/s)-f2.Eval(b); nb++; } }
printf("V36PILOT kept/fr %.0f | residual slope %.2e /tbin (tgt ~0) | residual bump x%.3f (tgt 1.00)\n",
  kpf, f2.GetParameter(1), TMath::Exp(bump/nb));' 2>&1 | grep V36PILOT
echo "pilot complete"
