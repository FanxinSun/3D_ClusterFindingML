#!/bin/bash
# world_scan_det.sh — Detroit-based reshape grid (2026-09-02): official MDC2/Detroit cfg + overrides
# pT0Ref {1.15,1.25,1.35} x StringPT:sigma {0.335,0.30}, flavour probStoUD 0.27 / probQQtoQ 0.11; 200k ev/point.
cd "$(dirname "$0")"; N=200000; SEED=20260960; i=0; jobs=(); CFG=../../island_post/official_pp_mb_mdc2.cfg
run() { i=$((i+1)); ./gen_world $N $((SEED+i)) "$1" cfgset "$2" > "$1.log" 2>&1 & jobs+=($!); if [ ${#jobs[@]} -ge 3 ]; then wait -n; jobs=($(jobs -pr)); fi; }
for p in 1.15 1.25 1.35; do for s in 0.335 0.300; do
  run "det_p${p}_s${s}_k0.270_q0.110" "$CFG|MultipartonInteractions:pT0Ref=$p;StringPT:sigma=$s;StringFlav:probStoUD=0.27;StringFlav:probQQtoQ=0.11"; done; done
wait; echo DET-DONE
