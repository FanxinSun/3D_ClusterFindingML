#!/usr/bin/env bash
# v3.6 readout fix (2026-07-14): the v36 production script inherited the
# PRE-revB readout call (p2 flat 7e-5 instead of table scale 1.0, no sigma_ped,
# no minis, p_trig 0.26) — sub10 collapsed to 1.4e-8. Readout-only rerun on the
# GOOD v3.6 frames with the full revC parameter set.
# was v3.5 revision C: near-threshold defect report, part 2.
# Adds on top of revB: per-pad pedestal spread sigma_ped=5 (real saturation BUMP
# mean/rms 948.8/5.8; sim was a delta at 949), readout-stage mini-cluster background
# (low-adc 2px {209,65,119}/fr + 3-4px {166,30,32}/fr, R2/R3 3px-only, pixel adu
# thr+1+Exp(7 / {5,2.5,2.5}), trk=-7 -> cls=2), both pilot-gated. Judged vs v3.4
# as the fallback baseline.
# revision B was: readout-only rerun on frames_pau_production_v35
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
  TFile*fi=TFile::Open(\"frames_pau_production_v36.root\");
  TNtuple*r=(TNtuple*)fi->Get(\"raw_pix\"); TNtuple*s=(TNtuple*)fi->Get(\"frame_truth\");
  TFile*fo=new TFile(\"fB_$i.root\",\"RECREATE\");
  TTree*rb=r->CopyTree(\"event>=$LO&&event<$HI\"); rb->SetName(\"raw_pix\"); rb->Write();
  TTree*sb=s->CopyTree(\"event>=$LO&&event<$HI\"); sb->SetName(\"frame_truth\"); sb->Write();
  printf(\"batch $i: %lld raw, %lld truth\n\", rb->GetEntries(), sb->GetEntries());
  fo->Close();" 2>&1 | grep "batch $i"
  root -l -b -q -e "
  gROOT->ProcessLine(\".L tpc_digitize.C+\");
  tpc_readout(\"fB_$i.root\",\"dP_$i.root\",0.93,20.0,1,1,4711,\"$DM\",11.0,0.39,0.55,1.0,0.021,7.0,36.0,70.0,11.0,940.0,2,0.29,10.0,1.24,1.06,-1.0,5.0,1.0);" 2>&1 | grep "tpc_readout: fB"
  root -l -b -q "islandize91.C+(\"dP_$i.root\",\"i91_$i.root\",1,\"\",\"fB_$i.root\")" 2>&1 | grep "islandize91: dP"
  rm -f fB_$i.root
done
hadd -f digi_frames_production_v36.root dP_*.root > /dev/null 2>&1
hadd -f island91_frames_production_v36.root i91_*.root > /dev/null 2>&1
rm -f dP_*.root i91_*.root
root -l -b -q -e 'gROOT->ProcessLine(".L islandize.C+"); islandize("digi_frames_production_v36.root","island_frames_v36.root",1);' 2>&1 | grep islandize
echo "=== V3.6 READOUT-FIX COMPLETE ==="
