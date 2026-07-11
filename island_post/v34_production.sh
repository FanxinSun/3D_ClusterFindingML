#!/usr/bin/env bash
# P3 production v3.4 (B4.3: REAL region-dependent deterministic ZS + region gains).
# READOUT-ONLY re-production off the frozen v3.2 composed frames (same collisions,
# flash, rates, seeds). Changes vs v3.3:
#   - ZS = CDB payload TpcADUThresholds10R1_20R23 (run 79507 per-region cliffs
#     confirmed: R1 keeps adc>=11, R2/R3 keep adc>=21); the B3 band/retention
#     machinery is OFF — the global "band 11-19 + step at 20" was purely
#     R1(thr10) stacked on R23(thr20).
#   - p2 = 7e-5 STANDALONE sub-threshold trace (real level).
#   - gaincal 0.93 (R2-anchored) x region gains gR1 1.24 / gR3 1.06
#     (phenomenological stand-in for the CAEN per-module corrections = B4.4).
#   - B4.2 saturation disturbance retained; emission floor follows the region
#     threshold; p_trig re-gated 0.33 -> 0.26 (regional ZS lengthens runs).
# Island gate (40 fr): zsize P13/P20/P30 x1.02/1.02/1.03, <zsize> -2.1%, region
# shares -1.7/-4.0/+7.0%, per-region pixmeans -5.1/-0.3/-3.1%, global pixmean
# -1.4%. Declared: <phisize> +4.9% (CLOUD recalib queued = stage A), sh/hi +18%
# (R1 spectral shape, B4.4 CAEN), R3-vs-R2 share imbalance ~+-4-7% (content).
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
  tpc_readout(\"fB_$i.root\",\"dP_$i.root\",0.93,20.0,1,1,4711,\"$DM\",11.0,0.39,0.55,0.00007,0.021,7.0,36.0,70.0,11.0,940.0,2,0.26,10.0,1.24,1.06);" 2>&1 | grep "tpc_readout: fB"
  root -l -b -q "islandize91.C+(\"dP_$i.root\",\"i91_$i.root\",1,\"\",\"fB_$i.root\")" 2>&1 | grep "islandize91: dP"
  rm -f fB_$i.root
done
hadd -f digi_frames_production_v34.root dP_*.root > /dev/null 2>&1
hadd -f island91_frames_production_v34.root i91_*.root > /dev/null 2>&1
rm -f dP_*.root i91_*.root
root -l -b -q -e 'gROOT->ProcessLine(".L islandize.C+"); islandize("digi_frames_production_v34.root","island_frames_v34.root",1);' 2>&1 | grep islandize
echo "=== V3.4 PRODUCTION COMPLETE ==="
