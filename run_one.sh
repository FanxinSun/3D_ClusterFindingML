#!/usr/bin/env bash
# run_one.sh — non-interactive single sim job (bypasses the run_batch_sim.sh menu).
# Parallel-safe: does NOT rewrite localAlignmentParamsFile.txt (uses the one in place).
# Bakes in the validated exam5-era config: ana537 alignment (pre-transcoded file),
# gap-rebuild 3D field, merged static+module-edge distortion, SVTX_SCAN_ALL,
# ported TrkrNtuplizer, 37 us extended readout, TPC dead-channel masking.
#
# usage: run_one.sh <nevents> <output.root> [hepmc_file] [pileup_rate_hz]
#   no hepmc_file -> particle-gun mode (macro defaults: 60 pi-/evt, pT 0.1-20)
set -euo pipefail
NEV=$1; OUT=$2; HEPMC=${3:-}; PRATE=${4:-}

IMG=/cvmfs/sphenix.opensciencegrid.org/singularity/rhic_sl7_ext
SETUP=/cvmfs/sphenix.opensciencegrid.org/gcc-8.3/opt/sphenix/core/bin/sphenix_setup.sh
REPO=$HOME/sPHENIX/3D_ClusterFindingML
MACRODIR=$REPO/macros-offline/detectors/sPHENIX
COMMON=$REPO/macros-offline/common
PORT=$REPO/trkrntuplizer_port

export SINGULARITYENV_USE_CDB_OFFLINE=1
export SINGULARITYENV_MAG_FIELD=$REPO/CDB_offline/FIELDMAP_GAP/65/a9/65a930ed6de9c0e049cd0f3ef226e6b4_sphenix3dbigmapxyz_gap_rebuild_v2.root
export SINGULARITYENV_DISTORTION_MODE=merged
export SINGULARITYENV_SVTX_SCAN_ALL=1
export SINGULARITYENV_TRKRNTUP_SO=$PORT/libTrkrNtuplizerMin.so
export SINGULARITYENV_TRKRNTUP_H=$PORT/TrkrNtuplizerMin.h
export SINGULARITYENV_TRKRNTUP_OUT=${OUT%.root}_trkrntuple.root
# 14000 ns = ana.331 recording ceiling (padplane MaxT = 2x full drift = 27947 ns);
# larger values shift the clusterizer z-map via the shared geometry -> cluster corruption
export SINGULARITYENV_TPC_EXTENDED_READOUT_NS=${TPC_EXTENDED_READOUT_NS:-14000}
export SINGULARITYENV_TPC_DEADMAP=${TPC_DEADMAP:-$REPO/CDB_offline/TPC_DEADCHANNELMAP/ff/c3/ffc3f6498934c5a8ba31065292c6ebcc_TPCDeadMap_79471.root}
export SINGULARITYENV_TPC_DEADMASK_H=$PORT/TpcDeadMasker.h
if [ -n "$HEPMC" ]; then export SINGULARITYENV_HEPMC_FILE=$HEPMC; fi
if [ -n "$PRATE" ]; then
  export SINGULARITYENV_PILEUPRATE=$PRATE
  export SINGULARITYENV_PILEUP_FILE=${PILEUP_FILE:-$REPO/hepmc_auau/auau200_pileup.dat}
fi

exec singularity exec -B /cvmfs:/cvmfs "$IMG" bash -lc \
  "source '$SETUP' -n ana.331 >/dev/null 2>&1; export ROOT_INCLUDE_PATH='$COMMON':'$MACRODIR':\$ROOT_INCLUDE_PATH; cd '$MACRODIR' && root -b -q 'Fun4All_DistortionSim.C($NEV, \"\", \"$OUT\")'"
