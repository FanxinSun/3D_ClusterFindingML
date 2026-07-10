#!/usr/bin/env bash
set -uo pipefail
cd /home/rog/sPHENIX/3D_ClusterFindingML/island_post
M=/home/rog/sPHENIX/3D_ClusterFindingML/macros-offline/detectors/sPHENIX
root -l -b -q -e 'gROOT->ProcessLine(".L tpc_digitize.C+");
tpc_transport("'$M'/pau_a_eval_g4svtx_eval.root","raw_lib_pau.root",20);
tpc_transport("'$M'/pauL_a_eval_g4svtx_eval.root","raw_lib_pauLa.root",75);
tpc_transport("'$M'/pauL_b_eval_g4svtx_eval.root","raw_lib_pauLb.root",75);
tpc_transport("'$M'/pauL_c_eval_g4svtx_eval.root","raw_lib_pauLc.root",75);
tpc_transport("'$M'/pauL_d_eval_g4svtx_eval.root","raw_lib_pauLd.root",75);
tpc_transport("/home/rog/sPHENIX/3D_ClusterFindingML/laser_try/laser_cm_g4hits.root","raw_lib_cmflash.root",1);' 2>&1 | grep transport:
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
  double w=(lay<=22)?1.25:(lay>=39?0.40:1.0);
  o->Fill(ev,lay,sd,pad,tb,q*w,trk);}
o->Write(); fo->Close(); printf("weighted flash lib rebuilt\n");' 2>&1 | grep rebuilt
echo "RETRANSPORT DONE"
