#!/usr/bin/env bash
# P5 event-display build: same host toolchain as build_standalone.sh, but linked
# WITH the GUI/vis components (build_standalone.sh does not need them).
set -e
cd "$(dirname "$0")"
source /home/rog/geant4/bin/geant4.sh
# The local 11.4.2 install advertises -DG4VIS_USE_OIX/-DG4VIS_USE_OI in
# `geant4-config --cflags` but never installed G4OpenInventorX.hh /
# G4OpenInventorXtExtended.hh, so G4VisExecutive.icc will not compile with them
# on. Undefine after the config flags — every other driver (OpenGL X/Xm/Qt,
# ToolsSG, RayTracer, VTK) is complete and stays registered.
g++ -O2 -std=c++$(root-config --cxxstandard) -o p5_display p5_display.cc \
  $(geant4-config --cflags --libs) \
  -UG4VIS_USE_OIX -UG4VIS_USE_OI \
  $(root-config --cflags --libs) \
  -lxerces-c
echo "build OK -> $(pwd)/p5_display"
