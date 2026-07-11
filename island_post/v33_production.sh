#!/usr/bin/env bash
# P3 production v3.3 (B4.2: saturation-triggered baseline disturbance).
# READOUT-ONLY re-production: reuses the v3.2 composed raw frames verbatim
# (frames_pau_production_v32.root — same collisions, flash, rates, seeds), so v3.3
# differs from v3.2 ONLY in stage B. New B4.2 slow component on top of frozen B4.1:
#   fixed-duration LINEAR-ramp disturbance (ALICE slow-ion-component shape,
#   arXiv:2304.03881), fired by ADC saturation (>=940) with p_trig=0.33 (one
#   retriggerable disturbance per pad), amplitude 36 ADU (lognormal sigma 0.35,
#   mean-preserving), duration 70 tbins (gaussian sigma 0.30) — per-trigger
#   dispersion smooths the otherwise-bimodal island zsize tail — emit-only
#   (synthetic samples in empty tbins; no boost to real samples), emission
#   floor 11 (chains die through the retention band).
# Winner of the 8-round scan vs run 79507 per-pad run-length tail; island gate:
# P(zsize>=13/20/30) x0.91/x1.13/x1.07 of real (v3.2 was x0.35/x0.22/x0.17).
# Batched (50 frames) because tpc_readout sorts in memory.
set -euo pipefail
cd /home/rog/sPHENIX/3D_ClusterFindingML/island_post
DM=/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/TPC_DEADCHANNELMAP/ff/c3/ffc3f6498934c5a8ba31065292c6ebcc_TPCDeadMap_79471.root
NB=5; PER=50
for i in $(seq 0 $((NB-1))); do
  LO=$((i*PER)); HI=$(((i+1)*PER))
  root -l -b -q -e "
  TFile*fi=TFile::Open(\"frames_pau_production_v32.root\");
  TNtuple*r=(TNtuple*)fi->Get(\"raw_pix\"); TNtuple*s=(TNtuple*)fi->Get(\"frame_truth\");
  TFile*fo=new TFile(\"fB_$i.root\",\"RECREATE\");
  TTree*rb=r->CopyTree(\"event>=$LO&&event<$HI\"); rb->SetName(\"raw_pix\"); rb->Write();
  TTree*sb=s->CopyTree(\"event>=$LO&&event<$HI\"); sb->SetName(\"frame_truth\"); sb->Write();
  printf(\"batch $i: %lld raw, %lld truth\n\", rb->GetEntries(), sb->GetEntries());
  fo->Close();" 2>&1 | grep "batch $i"
  root -l -b -q -e "
  gROOT->ProcessLine(\".L tpc_digitize.C+\");
  tpc_readout(\"fB_$i.root\",\"dP_$i.root\",0.87,20.0,1,1,4711,\"$DM\",11.0,0.39,0.55,0.0005,0.021,7.0,36.0,70.0,11.0,940.0,2,0.33);" 2>&1 | grep "tpc_readout: fB"
  root -l -b -q "islandize91.C+(\"dP_$i.root\",\"i91_$i.root\",1,\"\",\"fB_$i.root\")" 2>&1 | grep "islandize91: dP"
  rm -f fB_$i.root
done
hadd -f digi_frames_production_v33.root dP_*.root > /dev/null 2>&1
hadd -f island91_frames_production_v33.root i91_*.root > /dev/null 2>&1
rm -f dP_*.root i91_*.root
root -l -b -q -e 'gROOT->ProcessLine(".L islandize.C+"); islandize("digi_frames_production_v33.root","island_frames_v33.root",1);' 2>&1 | grep islandize
echo "=== V3.3 PRODUCTION COMPLETE ==="
