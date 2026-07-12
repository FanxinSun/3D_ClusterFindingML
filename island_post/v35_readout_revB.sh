#!/usr/bin/env bash
# v3.5 revision B (2026-07-12): readout-only rerun on frames_pau_production_v35
# fixing the near-threshold panel (user defect report):
#   - SHAPED sub-threshold trace: per-ADU keep-probability tables (per region),
#     fitted bin-wise as real-trace/raw and iterated once at the working point;
#     verified kept/target = 0.97-1.10 across all bins. p2 arg = table scale (1.0).
#   - band lever scanned (sigma_pad_R1 0.30-0.95 x gR1): sh/hi floor 0.273 vs
#     target 0.233 -> structurally locked (R1 raw spectral shape; external gain
#     map item). Kept uniform sigma_pad 0.55; sh/hi ~0.279 declared.
set -euo pipefail
cd /home/rog/sPHENIX/3D_ClusterFindingML/island_post
DM=/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/TPC_DEADCHANNELMAP/ff/c3/ffc3f6498934c5a8ba31065292c6ebcc_TPCDeadMap_79471.root
NB=5; PER=50
for i in $(seq 0 $((NB-1))); do
  LO=$((i*PER)); HI=$(((i+1)*PER))
  root -l -b -q -e "
  TFile*fi=TFile::Open(\"frames_pau_production_v35.root\");
  TNtuple*r=(TNtuple*)fi->Get(\"raw_pix\"); TNtuple*s=(TNtuple*)fi->Get(\"frame_truth\");
  TFile*fo=new TFile(\"fB_$i.root\",\"RECREATE\");
  TTree*rb=r->CopyTree(\"event>=$LO&&event<$HI\"); rb->SetName(\"raw_pix\"); rb->Write();
  TTree*sb=s->CopyTree(\"event>=$LO&&event<$HI\"); sb->SetName(\"frame_truth\"); sb->Write();
  printf(\"batch $i: %lld raw, %lld truth\n\", rb->GetEntries(), sb->GetEntries());
  fo->Close();" 2>&1 | grep "batch $i"
  root -l -b -q -e "
  gROOT->ProcessLine(\".L tpc_digitize.C+\");
  tpc_readout(\"fB_$i.root\",\"dP_$i.root\",0.93,20.0,1,1,4711,\"$DM\",11.0,0.39,0.55,1.0,0.021,7.0,36.0,70.0,11.0,940.0,2,0.29,10.0,1.24,1.06);" 2>&1 | grep "tpc_readout: fB"
  root -l -b -q "islandize91.C+(\"dP_$i.root\",\"i91_$i.root\",1,\"\",\"fB_$i.root\")" 2>&1 | grep "islandize91: dP"
  rm -f fB_$i.root
done
hadd -f digi_frames_production_v35.root dP_*.root > /dev/null 2>&1
hadd -f island91_frames_production_v35.root i91_*.root > /dev/null 2>&1
rm -f dP_*.root i91_*.root
root -l -b -q -e 'gROOT->ProcessLine(".L islandize.C+"); islandize("digi_frames_production_v35.root","island_frames_v35.root",1);' 2>&1 | grep islandize
echo "=== V3.5revB PRODUCTION COMPLETE ==="
