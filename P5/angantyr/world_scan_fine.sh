#!/bin/bash
# world_scan_fine.sh — refined reshape grid around the coarse-scan optimum (2026-09-02):
# pT0Ref {2.50,2.55,2.60} x StringPT:sigma {0.36,0.38,0.40}, flavour fixed at probStoUD 0.27 / probQQtoQ 0.11;
# 200k events/point, 3-way parallel.
cd "$(dirname "$0")"; N=200000; SEED=20260930; i=0; jobs=()
run() { i=$((i+1)); ./gen_world $N $((SEED+i)) "$1" set "$2" > "$1.log" 2>&1 & jobs+=($!); if [ ${#jobs[@]} -ge 3 ]; then wait -n; jobs=($(jobs -pr)); fi; }
for p in 2.50 2.55 2.60; do for s in 0.360 0.380 0.400; do
  run "fine_p${p}_s${s}_k0.270_q0.110" "MultipartonInteractions:pT0Ref=$p;StringPT:sigma=$s;StringFlav:probStoUD=0.27;StringFlav:probQQtoQ=0.11"; done; done
wait; echo FINE-DONE
