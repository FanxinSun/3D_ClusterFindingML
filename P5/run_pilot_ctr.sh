#!/usr/bin/env bash
# P5 pilot inside the EL9 container (native cling JIT needs EL9 system headers;
# host-run died on dictionary registration). alma9.sif + /cvmfs + home bound;
# el9libs provides libicu 67 that the minimal image lacks.
cd /home/rog/sPHENIX/3D_ClusterFindingML/P5
apptainer exec -B /cvmfs/sphenix.opensciencegrid.org -B /home/rog alma9_sb bash -c '
source /cvmfs/sphenix.opensciencegrid.org/alma9.2-gcc-14.2.0/opt/sphenix/core/bin/sphenix_setup.sh -n
export LD_LIBRARY_PATH=/home/rog/sPHENIX/3D_ClusterFindingML/P5/el9libs:$LD_LIBRARY_PATH
export PATH=/home/rog/sPHENIX/3D_ClusterFindingML/P5/el9bin:$PATH
R=/home/rog/sPHENIX/3D_ClusterFindingML
export ROOT_INCLUDE_PATH=$R/P5:$R/modern_macros/detectors/sPHENIX:$R/modern_macros/common:$ROOT_INCLUDE_PATH
cd $R/P5
root.exe -q -b "Fun4All_P5_pilot.C(1, \"$R/hepmc_pau/pau200.dat\", \"P5_pilot.root\", \"\", 0, \".\")"
echo "=== P5 PILOT G4 EXIT $? ==="
ls -la *.root 2>/dev/null'
