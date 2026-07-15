#!/usr/bin/env bash
# v4.0 scale-up orchestrator: 1500 new events (total 2000) -> G4 x8 -> ready
# for census+thin+production. Stage markers: SCALE-A/B-DONE.
cd /home/rog/sPHENIX/3D_ClusterFindingML/P5/angantyr
for i in 4 5 6 7 8 9; do
  ./gen_pau 250 ang_run_$i.dat $((20260720+i)) > gen_$i.log 2> fired_$i.txt &
done
wait
echo "SCALE-A-DONE (generation)"
cd /home/rog/sPHENIX/3D_ClusterFindingML/P5
source /home/rog/geant4/bin/geant4.sh
FM=/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/FIELDMAP_GAP/65/a9/65a930ed6de9c0e049cd0f3ef226e6b4_sphenix3dbigmapxyz_gap_rebuild_v2.root
# 6 new files x 250 events, 6 parallel G4 jobs (existing 4x125 already done)
for i in 4 5 6 7 8 9; do
  ./standalone_tpc sphenix_p5.gdml angantyr/ang_run_$i.dat ANG_g4hit_$i.root 250 $FM > ang_g4_$i.log 2>&1 &
done
wait
echo "SCALE-B-DONE (G4)"
