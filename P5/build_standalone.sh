#!/usr/bin/env bash
# P5 standalone app build: local Geant4 11.4.2 + local ROOT (host toolchain).
set -e
source /home/rog/geant4/bin/geant4.sh
g++ -O2 -o standalone_tpc standalone_tpc.cc \
  $(geant4-config --cflags --libs) \
  $(root-config --cflags --libs) \
  -lxerces-c
echo "build OK"
