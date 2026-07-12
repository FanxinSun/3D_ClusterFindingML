#!/usr/bin/env bash
set -uo pipefail
# v3.5 retransport (2026-07-12): CLOUD 0.106 (recalibrated under B4.3 regional ZS;
# phi-size proxy 3.435 vs target 3.43) + flash region weights REFIT for thr-11 R1
# (excess-over-baseline real/v34: x0.24/x0.55/x0.67 on the old 1.25/1.00/0.40 ->
# 0.298/0.547/0.268, first pass; pilot-composed before full production).
cd /home/rog/sPHENIX/3D_ClusterFindingML/island_post
M=/home/rog/sPHENIX/3D_ClusterFindingML/macros-offline/detectors/sPHENIX
root -l -b -q -e 'gROOT->ProcessLine(".L tpc_digitize.C+");
tpc_transport("'$M'/pau_a_eval_g4svtx_eval.root","raw_lib_pau.root",20,0.106);
tpc_transport("'$M'/pauL_a_eval_g4svtx_eval.root","raw_lib_pauLa.root",75,0.106);
tpc_transport("'$M'/pauL_b_eval_g4svtx_eval.root","raw_lib_pauLb.root",75,0.106);
tpc_transport("'$M'/pauL_c_eval_g4svtx_eval.root","raw_lib_pauLc.root",75,0.106);
tpc_transport("'$M'/pauL_d_eval_g4svtx_eval.root","raw_lib_pauLd.root",75,0.106);
tpc_transport("/home/rog/sPHENIX/3D_ClusterFindingML/laser_try/laser_cm_g4hits.root","raw_lib_cmflash.root",1,0.106);' 2>&1 | grep transport:
# rebuild region-weighted flash lib
root -l -b -q -e '
TFile*fi=TFile::Open("raw_lib_cmflash.root"); TTree*t=(TTree*)fi->Get("raw_pix");
float ev,lay,sd,pad,tb,q,trk;
t->SetBranchAddress("event",&ev); t->SetBranchAddress("layer",&lay); t->SetBranchAddress("side",&sd);
t->SetBranchAddress("pad",&pad); t->SetBranchAddress("tbin",&tb); t->SetBranchAddress("q",&q); t->SetBranchAddress("trk",&trk);
TFile*fo=new TFile("raw_lib_cmflash_w.root","RECREATE");
TNtuple*o=new TNtuple("raw_pix","region-weighted CM flash","event:layer:side:pad:tbin:q:trk");
Long64_t N=t->GetEntries();
for(Long64_t i=0;i<N;++i){t->GetEntry(i);
  double w=(lay<=22)?0.298:(lay>=39?0.268:0.547);
  o->Fill(ev,lay,sd,pad,tb,q*w,trk);}
o->Write(); fo->Close(); printf("weighted flash lib rebuilt\n");' 2>&1 | grep rebuilt
echo "RETRANSPORT DONE"
