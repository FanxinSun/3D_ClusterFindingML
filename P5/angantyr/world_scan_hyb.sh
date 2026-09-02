#!/bin/bash
# world_scan_hyb.sh — hybrid round (2026-09-02): isolate WHICH Detroit ingredient hardens the spectrum
# and starves soft multiplicity, to find a Monash/Detroit interpolation that fixes W1's soft tilt. 200k ev/pt.
cd "$(dirname "$0")"; N=200000; SEED=20260970; i=0; jobs=(); CFG=../../island_post/official_pp_mb_mdc2.cfg
FLAV="StringFlav:probStoUD=0.27;StringFlav:probQQtoQ=0.11"
run() { i=$((i+1)); ./gen_world $N $((SEED+i)) "$1" "$2" "$3" > "$1.log" 2>&1 & jobs+=($!); if [ ${#jobs[@]} -ge 3 ]; then wait -n; jobs=($(jobs -pr)); fi; }
run hyb_H1_det_pdf13_p1.40      cfgset "$CFG|PDF:pSet=13;MultipartonInteractions:pT0Ref=1.40;$FLAV"
run hyb_H2_det_pdf13_p1.25      cfgset "$CFG|PDF:pSet=13;MultipartonInteractions:pT0Ref=1.25;$FLAV"
run hyb_H3_mon_CR54_p2.55       set    "ColourReconnection:range=5.4;ColourReconnection:mode=2;TimeShower:alphaSvalue=0.18;MultipartonInteractions:pT0Ref=2.55;StringPT:sigma=0.36;$FLAV"
run hyb_H4_mon_bprof2_p2.55     set    "MultipartonInteractions:bProfile=2;MultipartonInteractions:coreRadius=0.56;MultipartonInteractions:coreFraction=0.78;MultipartonInteractions:pT0Ref=2.55;StringPT:sigma=0.36;$FLAV"
run hyb_H5_mon_bprof2_p2.30     set    "MultipartonInteractions:bProfile=2;MultipartonInteractions:coreRadius=0.56;MultipartonInteractions:coreFraction=0.78;MultipartonInteractions:pT0Ref=2.30;StringPT:sigma=0.36;$FLAV"
run hyb_H6_mon_pdf17_p2.55      set    "PDF:pSet=17;MultipartonInteractions:pT0Ref=2.55;StringPT:sigma=0.36;$FLAV"
wait; echo HYB-DONE
