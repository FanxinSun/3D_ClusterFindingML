#!/bin/bash
# copy_cdb_smart.sh — Copy all 177 CDB domains from CVMFS, prioritizing
# tracking domains first, then the rest. Uses parallel rsync for speed.
#
# Total on CVMFS: 177 domains, ~205K files, ~21 GB
# Tracking-relevant domains are small & fast; calorimeter domains are huge.

SRC="/cvmfs/sphenix.opensciencegrid.org/calibrations/sphnxpro/cdb"
DST="/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline"

if [ ! -d "$SRC" ]; then
  echo "ERROR: CVMFS source not mounted at $SRC" >&2
  exit 1
fi

copy_domain() {
  local domain="$1"
  local domain_dir="$SRC/$domain/"
  if [ ! -d "$domain_dir" ]; then
    echo "  SKIP: $domain (not found)"
    return
  fi
  local nfiles=$(find "$domain_dir" -type f -not -name ".cvmfscatalog" 2>/dev/null | wc -l)
  echo "  Copying $domain ($nfiles files)..."
  mkdir -p "$DST/$domain"
  rsync -a --exclude='.cvmfscatalog' "$domain_dir" "$DST/$domain/" 2>/dev/null
  if [ $? -eq 0 ]; then
    echo "    ✓ $domain done"
  else
    echo "    ✗ $domain rsync failed, using cp..."
    find "$domain_dir" -type f -not -name ".cvmfscatalog" | while read f; do
      rel="${f#$domain_dir}"
      mkdir -p "$(dirname "$DST/$domain/$rel")"
      cp "$f" "$DST/$domain/$rel"
    done
    echo "    ✓ $domain done (cp)"
  fi
}

# === PHASE 1: Tracking-relevant domains (small, fast, high priority) ===
echo "=========================================="
echo "PHASE 1: Tracking & geometry domains"
echo "=========================================="

TRACKING_DOMAINS=(
  TRACKINGALIGNMENT
  Tracking_Geometry
  ACTSGEOMETRYCONFIG
  ACTSMATERIALMAP
  FIELDMAPTRACKING
  FIELDMAP_TRACKING
  MVTX_ALIGNMENT
  MVTX_DeadPixelMap
  MVTX_HotPixelMap
  InttFeeMap
  InttSurveyMap
  INTT_BCOMAP
  INTT_DACMAP
  INTT_Hotmap
  INTT_HotMap
  INTT_HotMap_default
  INTT_STREAMING_FEE_OFFSET
  TPC_DEADCHANNELMAP
  TPC_HOTCHANNELMAP
  TPC_FEE_CHANNEL_MAP
  TPC_DRIFT_VELOCITY
  TPC_DRIFT_VELOCITY_default
  TPC_STATIC_DISTORTION
  TPC_STATIC_CORRECTION_MODEL
  TPC_TIMEORDERED_DISTORTION
  TPC_ZS_THRESHOLDS
  TPC_TZERO_OFFSET
  TPC_TZERO_OFFSET_default
  TPC_LAMINATION_FIT_CORRECTION
  TPC_LAMINATION_FIT_CORRECTION_default
  TPC_Module_Edge
  TPC_DECODER_BAD_FEE
  TPC_DEDX_FITPARAM
  TPC_CAEN_CORRECTION_24974_1G
  TPC_CAEN_CORRECTION_24974_2G
  TPC_CAEN_CORRECTION_24974_5G
  TPOT_Pedestal
)

for d in "${TRACKING_DOMAINS[@]}"; do
  copy_domain "$d"
done

echo ""
echo "PHASE 1 COMPLETE. Tracking domains copied."
echo "Size so far: $(du -sh "$DST" | cut -f1)"
echo ""

# === PHASE 2: All remaining domains ===
echo "=========================================="
echo "PHASE 2: All remaining domains"
echo "=========================================="

# Build a set of already-copied domains
declare -A DONE
for d in "${TRACKING_DOMAINS[@]}"; do
  DONE[$d]=1
done

remaining=0
for domain_dir in "$SRC"/*/; do
  domain=$(basename "$domain_dir")
  if [ -z "${DONE[$domain]}" ]; then
    remaining=$((remaining + 1))
  fi
done
echo "$remaining remaining domains to copy"
echo ""

count=0
for domain_dir in "$SRC"/*/; do
  domain=$(basename "$domain_dir")
  if [ -z "${DONE[$domain]}" ]; then
    count=$((count + 1))
    echo "[$count/$remaining]"
    copy_domain "$domain"
  fi
done

echo ""
echo "=========================================="
echo "ALL DONE"
echo "=========================================="
echo "Total size: $(du -sh "$DST" | cut -f1)"
echo "Total files: $(find "$DST" -type f -not -name "*.sh" -not -name "*.yaml" -not -name ".cvmfscatalog" | wc -l)"
echo "Total domains: $(ls -d "$DST"/*/ 2>/dev/null | wc -l)"
