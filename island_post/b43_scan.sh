#!/bin/bash
# b43_scan.sh — B4.3 regional-ZS scan on 40 composed frames.
# Real payload TpcADUThresholds10R1_20R23 (run 79507 cliffs confirmed: R1>=11, R23>=21).
# Real targets: shares R1 0.362 R2 0.350 R3 0.287 | pixm R1 86.1 R2 105.0 R3 108.5
#               global pixm 99.2, kept/fr 257k; B4.2 run-tail targets unchanged.
set -e
cd "$(dirname "$0")"
DM=/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/TPC_DEADCHANNELMAP/ff/c3/ffc3f6498934c5a8ba31065292c6ebcc_TPCDeadMap_79471.root
RAW=raw_b42_scan.root
OUT=b43_scan_results.txt
if [ ! -f $RAW ]; then
  root -l -b -q -e '
  TFile*fi=TFile::Open("frames_pau_production_v32.root");
  TNtuple*r=(TNtuple*)fi->Get("raw_pix");
  TFile*fo=new TFile("raw_b42_scan.root","RECREATE");
  TTree*t=r->CopyTree("event<40"); t->Write("raw_pix");
  printf("subset: %lld rows\n", t->GetEntries()); fo->Close();' 2>&1 | grep subset
fi
run_cfg() {  # gain p2 ptrig gR1 gR3 tag
  root -l -b -q -e "gROOT->ProcessLine(\".L tpc_digitize.C+\"); tpc_readout(\"$RAW\",\"b43_tmp.root\",$1,20.0,1,1,4711,\"$DM\",11.0,0.39,0.55,$2,0.021,7.0,36.0,70.0,11.0,940.0,2,$3,10.0,$4,$5);" 2>&1 | grep "pixels kept"
  root -l -b -q -e "gROOT->ProcessLine(\".L b42_metrics.C+\"); b42_metrics(\"b43_tmp.root\",\"$6\");" 2>&1 | grep -E "B42METRIC|B43REGION" | tee -a $OUT
}
# round 2: region gains (R2-anchored g=0.93), p2 trimmed to real trace, p_trig re-gate
for cfg in "1.12 1.03" "1.18 1.03" "1.24 1.03" "1.12 1.06" "1.18 1.06" "1.24 1.06"; do
  set -- $cfg
  run_cfg 0.93 0.00007 0.36 $1 $2 "R5_gR1_${1}_gR3_${2}"
done
rm -f b43_tmp.root
echo "b43 scan complete -> $OUT"
