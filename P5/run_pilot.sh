#!/usr/bin/env bash
# P5 pilot: 2 pau200 HepMC events through the MODERN G4 chain (alma9/ana.558,
# native). Geometry full-on (material budget under test), all cell/reco stages
# off, DST out with G4Hits. Measures the secondary-flood hypothesis.
set -o pipefail
cd /home/rog/sPHENIX/3D_ClusterFindingML/P5
source /cvmfs/sphenix.opensciencegrid.org/alma9.2-gcc-14.2.0/opt/sphenix/core/bin/sphenix_setup.sh -n
# EL9 compat: alma9 stack needs ICU 67 (host ships 74+); extracted from the
# libicu el9 rpm into el9libs (see PIPELINE.md P5 log)
export LD_LIBRARY_PATH=/home/rog/sPHENIX/3D_ClusterFindingML/P5/el9libs:$LD_LIBRARY_PATH
R=/home/rog/sPHENIX/3D_ClusterFindingML
export ROOT_INCLUDE_PATH=$R/modern_macros/detectors/sPHENIX:$R/modern_macros/common:$ROOT_INCLUDE_PATH
root.exe -q -b "Fun4All_P5_pilot.C(2, \"$R/hepmc_pau/pau200.dat\", \"P5_pilot.root\", \"\", 0, \".\")"
echo "=== P5 PILOT G4 EXIT $? ==="
ls -la *.root
