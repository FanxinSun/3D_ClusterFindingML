#!/bin/bash
# v61_pilot.sh <SMOD_value> [tag] — width re-balance SOLVE pilot (chunk-0 CRN).
# GATED: run only after the probe thread's post-injection C(d) lands (or the
# user explicitly waives) — width_rebalance_request.md sequencing gate.
# Transport chunk 0 with the candidate smooth field (TWIST byte-frozen from
# twist_payload_v6.txt), gate-style readout (CRN seed 4711, mini off), then
# fc_pixels against the SMOD=0.0426 twisted baseline (v61pilot_base6t_digi).
# Baselines: v61pilot_base55_digi.root (no twist) / v61pilot_base6t_digi.root
# (twist, current SMOD) — produced 2026-08-23 by the harness A/B.
set -eo pipefail
cd "$(dirname "$0")"
SMOD=${1:?usage: v61_pilot.sh <SMOD_value> [tag]}
TAG=${2:-smod$SMOD}
DM=/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/TPC_DEADCHANNELMAP/ff/c3/ffc3f6498934c5a8ba31065292c6ebcc_TPCDeadMap_79471.root
P5=/home/rog/sPHENIX/3D_ClusterFindingML/P5
FIELD="0:0|0|${SMOD}|0|0|2.49|0.26"   # SPHI/SCM frozen at v5.4c values; SMOD = candidate
root -l -b -q -e "gROOT->ProcessLine(\".L tpc_digitize.C+\");
tpc_transport(\"$P5/PP_g4hit_0.root\",\"v61pilot_${TAG}_raw.root\",2000,0.040,0.0070,\"$FIELD\",\"twist_payload_v6.txt\");
tpc_readout(\"v61pilot_${TAG}_raw.root\",\"v61pilot_${TAG}_digi.root\",1.005,20.0,1,1,4711,\"$DM\",11.0,0.39,0.55,1.25,0.016,7.0,36.0,70.0,11.0,940.0,2,0.29,10.0,1.24,1.06,-1.0,5.0,0.0);" 2>&1 | grep -E "TWIST|rphi field|tpc_readout: v61" 
root -l -b -q -e "gROOT->ProcessLine(\".L ms_fieldcmp.C+\"); fc_pixels(\"v61pilot_base6t_digi.root\",\"v61pilot_${TAG}_digi.root\",2000,\"v61_${TAG}\");" 2>&1 | grep -E "ideal:|field:|FIELD SHARE"
rm -f "v61pilot_${TAG}_raw.root"
echo "PILOT ${TAG} DONE (candidate vs current-SMOD twisted baseline; negative share = width removed)"
