#!/bin/bash
# build.sh — rebuild the host-side laser-module runners (binaries are NOT in git).
#
# Both compile the collaboration module source VERBATIM against the stub headers
# in stubs/ (no Fun4All/Geant4 needed; only ROOT). Verified with the local
# coresoftware clone at the path below.
set -e
cd "$(dirname "$0")"
CS=/home/rog/sPHENIX/3D_ClusterFindingML/coresoftware/simulation/g4simulation/g4tpc

# CM diffuse-laser flash generator (PHG4TpcCentralMembrane).
# Usage: ./cm_gen <delay_ns> <n_events>   e.g. ./cm_gen 4346 300
g++ -O2 -std=c++17 -I stubs -I $CS $(root-config --cflags) \
    cm_gen.cc $CS/PHG4TpcCentralMembrane.cc $(root-config --libs) -o cm_gen

# DirectLaser generator (PHG4TpcDirectLaser) — trialed and EXCLUDED from the
# pipeline (endcap-peaked arrivals, no mid-drift spike); kept rebuildable for
# the record. -DSTUB_NO_UNITS because the module defines its own cm=1.0.
# Runtime needs CALIBRATIONROOT pointing at a dir with the fake pattern file
# (TTree "angles" with #theta/#phi branches) — see dl_gen.cc header comment.
# Usage: ./dl_gen <theta_deg> <phi_deg>   e.g. ./dl_gen 30 0
g++ -O2 -std=c++17 -DSTUB_NO_UNITS -I stubs -I $CS $(root-config --cflags) \
    dl_gen.cc $CS/PHG4TpcDirectLaser.cc $(root-config --libs) -o dl_gen

echo "built: cm_gen dl_gen"
