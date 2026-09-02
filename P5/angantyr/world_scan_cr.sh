#!/bin/bash
# world_scan_cr.sh — round 2 (last, pre-declared) on the winning ingredient: Monash + Detroit's colour
# reconnection (range 5.4, mode 2) + alphaS 0.18 flattens the pion slope; lower pT0Ref to restore world
# multiplicity. 200k ev/point.
cd "$(dirname "$0")"; N=200000; SEED=20260980; i=0; jobs=()
CR="ColourReconnection:range=5.4;ColourReconnection:mode=2;TimeShower:alphaSvalue=0.18"
run() { i=$((i+1)); ./gen_world $N $((SEED+i)) "$1" set "$2" > "$1.log" 2>&1 & jobs+=($!); if [ ${#jobs[@]} -ge 3 ]; then wait -n; jobs=($(jobs -pr)); fi; }
for p in 2.20 2.35; do for s in 0.300 0.360; do
  run "cr_p${p}_s${s}_k0.270_q0.110" "$CR;MultipartonInteractions:pT0Ref=$p;StringPT:sigma=$s;StringFlav:probStoUD=0.27;StringFlav:probQQtoQ=0.11"; done; done
run "cr_p2.20_s0.335_k0.300_q0.110" "$CR;MultipartonInteractions:pT0Ref=2.20;StringPT:sigma=0.335;StringFlav:probStoUD=0.30;StringFlav:probQQtoQ=0.11"
run "cr_p2.05_s0.335_k0.270_q0.110" "$CR;MultipartonInteractions:pT0Ref=2.05;StringPT:sigma=0.335;StringFlav:probStoUD=0.27;StringFlav:probQQtoQ=0.11"
wait; echo CR-DONE
