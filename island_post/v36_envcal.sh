#!/bin/bash
# tau transfer calibration: 100 flat-rate frames, same composition -> raw + readout slopes
set -e
cd "$(dirname "$0")"
DM=/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/TPC_DEADCHANNELMAP/ff/c3/ffc3f6498934c5a8ba31065292c6ebcc_TPCDeadMap_79471.root
M=/home/rog/sPHENIX/3D_ClusterFindingML/macros-offline/detectors/sPHENIX
LIBS="raw_lib_pau.root,raw_lib_pauLa.root,raw_lib_pauLb.root,raw_lib_pauLc.root,raw_lib_pauLd.root"
EVALS="$M/pau_a_eval_g4svtx_eval.root,$M/pauL_a_eval_g4svtx_eval.root,$M/pauL_b_eval_g4svtx_eval.root,$M/pauL_c_eval_g4svtx_eval.root,$M/pauL_d_eval_g4svtx_eval.root"
MBD="mbd_weights.txt|pau200.dat,pau_chunk_a.dat,pau_chunk_b.dat,pau_chunk_c.dat,pau_chunk_d.dat"
TAU=0
ENV="${1:-5.80:1.1537,11.10:1.1493,16.40:1.0862,21.70:0.9437,27.00:0.9853,32.30:0.9688,37.60:0.8798,42.90:0.7814,47.96:0.7445}"
root -l -b -q -e "
gROOT->ProcessLine(\".L frame_composer.C+\");
frame_composer(\"$LIBS\",\"fTAU.root\",100,275.,20260914,0,\"$EVALS\",\"\",0.,1.0,\"\",0.,0.,0.,0.,\"\",1.0,\"$MBD\",1.0,$TAU,\"$ENV\");
gROOT->ProcessLine(\".L tpc_digitize.C+\");
tpc_readout(\"fTAU.root\",\"dTAU.root\",0.93,20.0,1,1,4711,\"$DM\",11.0,0.39,0.55,1.0,0.021,7.0,36.0,70.0,11.0,940.0,2,0.29,10.0,1.24,1.06,-1.0,5.0,1.0);" 2>&1 | grep -E "frame_composer: 100|rate envelope|pixels kept"
root -l -b -q -e '
TFile*fr=TFile::Open("/home/rog/sPHENIX/3D_ClusterFindingML/clusters_seeds_island_79507-0.root_ntuplizer.root");
TTree*tr=(TTree*)fr->Get("ntp_hit");
TFile*fs=TFile::Open("dTAU.root"); TTree*ts=(TTree*)fs->Get("ntp_hit");
TH1D hr("hr","",971,-0.5,970.5), hs("hs","",971,-0.5,970.5);
tr->Draw("tbin>>hr","layer>=7&&layer<=54&&adc>0","goff"); ts->Draw("tbin>>hs","","goff");
auto wf=[](TH1D&h){ double tot=0,w1=0,w2=0,w3=0;
  for(int b=60;b<=950;++b){ if(b>=318&&b<=345) continue; double v=h.GetBinContent(b+1); tot+=v;
    if(b<=240) w1+=v; else if(b>=270&&b<=600) w2+=v; else if(b>=650) w3+=v; }
  printf("%.4f %.4f %.4f | ", w1/tot, w2/tot, w3/tot); };
printf("ENVCAL real: "); wf(hr); printf("sim: "); wf(hs); printf("\n");' 2>&1 | grep ENVCAL
echo "taucal complete"
