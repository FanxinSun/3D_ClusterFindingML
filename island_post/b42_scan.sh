#!/bin/bash
# b42_scan.sh — B4.2 slow-tail (tail2_frac, tail2_tau) grid scan on 40 composed frames.
# Baseline (0,-) first, then the grid. All other knobs frozen at v3.2 (B4.1).
# Targets (real, run 79507): run mean 3.066, P(>=10) 1.71e-2, P(>=20) 8.67e-3,
# P(>=30) 6.11e-3; guards: sub10 2.27e-4, sh/hi 0.233, pixmean/kept drift small.
set -e
cd "$(dirname "$0")"
DM=/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/TPC_DEADCHANNELMAP/ff/c3/ffc3f6498934c5a8ba31065292c6ebcc_TPCDeadMap_79471.root
RAW=raw_b42_scan.root
OUT=b42_scan_results.txt

run_cfg() {  # frac2 tau2 floor q0 emitonly ptrig tag
  root -l -b -q -e "gROOT->ProcessLine(\".L tpc_digitize.C+\"); tpc_readout(\"$RAW\",\"b42_tmp.root\",0.87,20.0,1,1,4711,\"$DM\",11.0,0.39,0.55,0.0005,0.021,7.0,$1,$2,$3,$4,$5,$6);" 2>&1 | grep "pixels kept"
  root -l -b -q -e "gROOT->ProcessLine(\".L b42_metrics.C+\"); b42_metrics(\"b42_tmp.root\",\"$7\");" 2>&1 | grep B42METRIC | tee -a $OUT
}

# round 8: winner-basin micro-tune — emission floor 11 (chains die through the
# retention band like real tails) and p_trig trim
for cfg in "36 70 0.40 11" "36 70 0.33 11" "40 70 0.35 11" "36 65 0.40 11"; do
  set -- $cfg
  run_cfg $1 $2 $4 940 2 $3 "R8_A${1}_T${2}_p${3}_fl${4}"
done
rm -f b42_tmp.root
echo "scan complete -> $OUT"
