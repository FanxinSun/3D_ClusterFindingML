#!/bin/bash
# tau transfer calibration: 100 flat-rate frames, same composition -> raw + readout slopes
set -e
cd "$(dirname "$0")"
DM=/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/TPC_DEADCHANNELMAP/ff/c3/ffc3f6498934c5a8ba31065292c6ebcc_TPCDeadMap_79471.root
M=/home/rog/sPHENIX/3D_ClusterFindingML/macros-offline/detectors/sPHENIX
LIBS="raw_lib_pau.root,raw_lib_pauLa.root,raw_lib_pauLb.root,raw_lib_pauLc.root,raw_lib_pauLd.root"
EVALS="$M/pau_a_eval_g4svtx_eval.root,$M/pauL_a_eval_g4svtx_eval.root,$M/pauL_b_eval_g4svtx_eval.root,$M/pauL_c_eval_g4svtx_eval.root,$M/pauL_d_eval_g4svtx_eval.root"
MBD="mbd_weights.txt|pau200.dat,pau_chunk_a.dat,pau_chunk_b.dat,pau_chunk_c.dat,pau_chunk_d.dat"
TAU=${1:-2330}
root -l -b -q -e "
gROOT->ProcessLine(\".L frame_composer.C+\");
frame_composer(\"$LIBS\",\"fTAU.root\",100,275.,20260914,0,\"$EVALS\",\"\",0.,1.0,\"\",0.,0.,0.,0.,\"\",1.0,\"$MBD\",1.0,$TAU);
gROOT->ProcessLine(\".L tpc_digitize.C+\");
tpc_readout(\"fTAU.root\",\"dTAU.root\",0.93,20.0,1,1,4711,\"$DM\",11.0,0.39,0.55,1.0,0.021,7.0,36.0,70.0,11.0,940.0,2,0.29,10.0,1.24,1.06,-1.0,5.0,1.0);" 2>&1 | grep -E "frame_composer: 100|pixels kept"
root -l -b -q -e '
auto slope=[](const char* fn, const char* tree, const char* tag){
  TFile*f=TFile::Open(fn); TTree*t=(TTree*)f->Get(tree);
  TH1D h("h","",971,-0.5,970.5); t->Draw("tbin>>h","","goff");
  h.Scale(1./h.Integral());
  TGraph g; TF1 fl("fl","pol1",270,950);
  for(int b=270;b<=950;++b){ double v=h.GetBinContent(b+1); if(v>0) g.SetPoint(g.GetN(),b,TMath::Log(v)); }
  g.Fit(&fl,"RQ");
  printf("TAUCAL %s slope %+.3e /tbin\n", tag, fl.GetParameter(1)); f->Close(); };
slope("fTAU.root","raw_pix","raw    ");
slope("dTAU.root","ntp_hit","readout");' 2>&1 | grep TAUCAL
echo "taucal complete"
