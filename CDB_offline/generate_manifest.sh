#!/bin/bash
# generate_manifest.sh — Create a YAML manifest of all CDB domains and their files.
# Runs against the CVMFS source to get accurate metadata.

SRC="/cvmfs/sphenix.opensciencegrid.org/calibrations/sphnxpro/cdb"
OUT="/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/cdb_manifest.yaml"

cat > "$OUT" <<'HEADER'
# sPHENIX Conditions Database (CDB) — Offline Payload Manifest
#
# Source: /cvmfs/sphenix.opensciencegrid.org/calibrations/sphnxpro/cdb/
# Generated from the publicly accessible CVMFS mirror of sPHENIX CDB payloads.
#
# IMPORTANT NOTES:
#   - These are the raw CDB payload files (ROOT histograms, text params, etc.)
#   - The CDB *lookup service* (run → payload mapping) requires BNL access,
#     but the payload files themselves are open on CVMFS.
#   - To use a TRACKINGALIGNMENT file directly, copy/symlink it to
#     ./localAlignmentParamsFile.txt in your working directory.
#   - For other domains, pass the file path directly to the relevant subsystem.
#
# Structure: domain_name/hash_prefix/hash_file
#   The hash directories (2-char hex) are content-addressed storage artifacts
#   from the CDB system. The actual filenames describe the payload contents.

HEADER

echo "domains:" >> "$OUT"

for domain_dir in "$SRC"/*/; do
  domain=$(basename "$domain_dir")
  nfiles=$(find "$domain_dir" -type f -not -name ".cvmfscatalog" 2>/dev/null | wc -l)
  total_size=$(du -sh "$domain_dir" 2>/dev/null | cut -f1)

  # Categorize by subsystem
  case "$domain" in
    CEMC*|cemc*|CALO*|calo*) category="Calorimeter_CEMC" ;;
    HCAL*|hcal*|ihcal*|ohcal*|IHCAL*|OHCAL*) category="Calorimeter_HCAL" ;;
    TPC_*) category="Tracking_TPC" ;;
    MVTX*) category="Tracking_MVTX" ;;
    INTT*|Intt*) category="Tracking_INTT" ;;
    TPOT*) category="Tracking_TPOT" ;;
    TRACKING*|Tracking*) category="Tracking_Global" ;;
    FIELD*) category="Tracking_FieldMap" ;;
    MBD*|mbd*|Centrality*|CENTRALITY*) category="EventCharacterization" ;;
    ZDC*) category="ForwardDetector" ;;
    sEPD*|SEPD*|sepd*) category="ForwardDetector" ;;
    GL1*|gl1*) category="DAQ_Trigger" ;;
    *) category="Other" ;;
  esac

  echo "  - name: \"$domain\"" >> "$OUT"
  echo "    category: \"$category\"" >> "$OUT"
  echo "    file_count: $nfiles" >> "$OUT"
  echo "    total_size: \"$total_size\"" >> "$OUT"
  echo "    files:" >> "$OUT"

  find "$domain_dir" -type f -not -name ".cvmfscatalog" -printf '%P %s %TY-%Tm-%Td\n' 2>/dev/null | sort | while read rel_path fsize fdate; do
    # Extract the human-readable filename (after the hash prefix)
    fname=$(basename "$rel_path")
    echo "      - path: \"$rel_path\"" >> "$OUT"
    echo "        filename: \"$fname\"" >> "$OUT"
    echo "        size_bytes: $fsize" >> "$OUT"
    echo "        date: \"$fdate\"" >> "$OUT"
  done

done

echo ""
echo "Manifest written to: $OUT"
echo "Total domains: $(grep -c 'name:' "$OUT")"
