#!/bin/bash
# world_scan.sh — factorized reshape scan on Monash (2026-09-02): (pT0Ref x StringPT:sigma) at default
# flavour, and (probStoUD x probQQtoQ) at a central (pT0Ref,sigma); 100k events/point, 3-way parallel.
cd "$(dirname "$0")"; N=100000; SEED=20260910; i=0; jobs=()
run() { i=$((i+1)); ./gen_world $N $((SEED+i)) "$1" set "$2" > "$1.log" 2>&1 & jobs+=($!); if [ ${#jobs[@]} -ge 3 ]; then wait -n; jobs=($(jobs -pr)); fi; }
for p in 2.28 2.45 2.60 2.80; do for s in 0.335 0.380 0.420; do
  run "scan_p${p}_s${s}_k0.217_q0.081" "MultipartonInteractions:pT0Ref=$p;StringPT:sigma=$s"; done; done
for k in 0.217 0.270 0.300; do for q in 0.081 0.095 0.110; do
  run "scan_p2.60_s0.380_k${k}_q${q}" "MultipartonInteractions:pT0Ref=2.60;StringPT:sigma=0.380;StringFlav:probStoUD=$k;StringFlav:probQQtoQ=$q"; done; done
wait; echo SCAN-DONE
