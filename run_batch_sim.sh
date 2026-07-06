#!/usr/bin/env bash
#
# run_batch_sim.sh - sPHENIX tracker simulation in batch mode
#

set -euo pipefail

# Pin collation so the generated menus are byte-identical on every machine. Without it, the
# `sort` calls that order the alignment + field-map option lists below follow $LC_COLLATE:
# a C / C.UTF-8 box sorts by codepoint ('_' after letters) while an en_US.UTF-8 box uses
# dictionary order ('_' ignored) -> the SAME files get DIFFERENT option numbers per machine
# (e.g. "field choice 2" was gap-rebuild on the VM but the tracking map on WSL). C is chosen
# because it needs no installed locale and is identical everywhere.
export LC_ALL=C

IMG=/cvmfs/sphenix.opensciencegrid.org/singularity/rhic_sl7_ext
SETUP=/cvmfs/sphenix.opensciencegrid.org/gcc-8.3/opt/sphenix/core/bin/sphenix_setup.sh
MACRODIR="$HOME/sPHENIX/3D_ClusterFindingML/macros-offline/detectors/sPHENIX"
COMMON="$HOME/sPHENIX/3D_ClusterFindingML/macros-offline/common"
MACRO="${1:-Fun4All_DistortionSim.C}"
CDB_DIR="$HOME/sPHENIX/3D_ClusterFindingML/CDB_offline"

echo "========================================="
echo "      sPHENIX Batch Simulation Config    "
echo "========================================="
echo "Select tracking ALIGNMENT  (labels = TRUE content; duplicates removed;"
echo "auto-transcoded to ana.331 7-column so ANY build's file is crash-safe):"

# Build a content-fingerprinted, de-duplicated option list from TRACKINGALIGNMENT.
# Fingerprint = which DETECTORS carry non-zero alignment, decoded from each row's
# hitsetkey:  layer = (key>>16)&0xFF  ->  0-2 MVTX, 3-6 INTT, 7-54 TPC, 55-56 TPOT.
# A detector is "aligned" if any of its rows has a non-zero param (fields 2..7).
# The label IS the decoded truth (no separate decode needed); the ana331/494/537 in
# the filename is only an origin tag, shown as provenance. TPC alignment is called out
# because TPC ntp_hit shifts only when TPC rows are aligned.
ALIGN_DIR="./CDB_offline/TRACKINGALIGNMENT"
declare -A _seen_md5
entries=()
if [ -d "$ALIGN_DIR" ]; then
  while IFS= read -r -d $'\0' f; do
    md5=$(md5sum "$f" | cut -c1-32)
    if [ -n "${_seen_md5[$md5]:-}" ]; then continue; fi
    _seen_md5[$md5]=1
    # sig = '+'-joined list of detectors with any non-zero alignment, or IDEAL
    sig=$(awk '{k=$1+0; L=int(k/65536)%256; det=(L<=2?"MVTX":(L<=6?"INTT":(L<=54?"TPC":"TPOT")));
                z=1; for(i=2;i<=7;i++) if($i+0!=0) z=0; if(!z) a[det]=1}
               END{n=split("MVTX INTT TPC TPOT",r," "); s="";
                   for(i=1;i<=n;i++) if(r[i] in a) s=s (s==""?"":"+") r[i];
                   print (s==""?"IDEAL":s)}' "$f")
    case "$sig" in
      IDEAL) rank=0; cls="Ideal (no misalignment)" ;;
      *TPC*) rank=2; cls="Aligned: $sig  (TPC aligned -> TPC ntp_hit shifts)" ;;
      *)     rank=1; cls="Aligned: $sig  (TPC ideal -> TPC ntp_hit unchanged)" ;;
    esac
    prov=$(basename "$f" | sed -E 's/^[0-9a-f]{32}_//;s/\.txt$//')
    entries+=("${rank}|${f}|${cls}  --  ${prov}")
  done < <(find "$ALIGN_DIR" -type f -name "*.txt" -print0 | sort -z)
fi

apaths=()
idx=1
echo "  0) Truth-only (no alignment, no reco chain -- geometry/display only)"
if [ "${#entries[@]}" -gt 0 ]; then
  while IFS= read -r e; do
    apaths+=("$(printf '%s' "$e" | cut -d'|' -f2)")
    printf "  %d) %s\n" "$idx" "$(printf '%s' "$e" | cut -d'|' -f3-)"
    idx=$((idx+1))
  done < <(printf '%s\n' "${entries[@]}" | sort -t'|' -k1,1n -k3)
fi

read -p "Alignment choice [1]: " ACHOICE
ACHOICE=${ACHOICE:-1}
if [ "$ACHOICE" == "0" ]; then
  export USE_CDB_OFFLINE=0
  echo "--> truth-only mode (reco chain disabled; no ntp_hit)"
elif [ "$ACHOICE" -ge 1 ] 2>/dev/null && [ "$ACHOICE" -lt "$idx" ]; then
  export USE_CDB_OFFLINE=1
  selected_file="${apaths[$((ACHOICE-1))]}"
  # AlignmentTransformation reads ./localAlignmentParamsFile.txt from CWD (=$MACRODIR).
  # Transcode to 7 columns: newer dumps add 3 global-rotation columns ana.331 can't parse
  # (SIGSEGV); first 7 fields keep the real local alignment, keys are build-identical.
  echo "--> Transcoding $(basename "$selected_file") -> 7-column ana.331 alignment"
  awk 'NF>=7{print $1, $2, $3, $4, $5, $6, $7}' "$selected_file" > "$MACRODIR/localAlignmentParamsFile.txt"
  echo "    wrote ${MACRODIR}/localAlignmentParamsFile.txt ($(wc -l < "$MACRODIR/localAlignmentParamsFile.txt") rows)"
else
  echo "Invalid choice -> truth-only mode."
  export USE_CDB_OFFLINE=0
fi

echo ""
read -e -p "Number of Events [100]: " NEVENTS
NEVENTS=${NEVENTS:-100}

read -e -p "Particle Type (e.g. pi-, e-, gamma) [pi-]: " PARTICLE_TYPE
export PARTICLE_TYPE=${PARTICLE_TYPE:-pi-}

read -e -p "Number of particles per event (e.g. 10, 60, 100) [60]: " PARTICLE_COUNT
export PARTICLE_COUNT=${PARTICLE_COUNT:-60}

read -e -p "Particle pT min (GeV/c) [0.1]: " PARTICLE_PT_MIN
export PARTICLE_PT_MIN=${PARTICLE_PT_MIN:-0.1}

read -e -p "Particle pT max (GeV/c) [20.0]: " PARTICLE_PT_MAX
export PARTICLE_PT_MAX=${PARTICLE_PT_MAX:-20.0}

# ---------- Magnetic field (menu with visualized paths) ----------
# The ana.331 runtime is ROOT 6.24; it reads maps up to ~6.26 fine, but a far-newer file
# (e.g. 6.38) corrupts ROOT cleanup and SIGSEGVs. ensure_field_usable() auto re-encodes
# any too-new map to a native-6.24 copy (cached under field_maps_root624/) before use.
REPO="$HOME/sPHENIX/3D_ClusterFindingML"
REENC_DIR="$REPO/field_maps_root624"
SAFE_VER=62699   # <= ~6.26.x is safe for the 6.24 runtime; beyond needs re-encode
# ROOT file version = big-endian 4-byte int at offset 4 of the TFile header.
# Read it directly (no ROOT process -> instant, can't hang the menu).
field_version () { python3 -c "import sys;f=open(sys.argv[1],'rb');f.seek(4);print(int.from_bytes(f.read(4),'big'))" "$1" 2>/dev/null || echo 0; }
ensure_field_usable () {            # $1=map path -> echoes a path the 6.24 runtime can read
  local m="$1" v bn out
  v=$(field_version "$m"); v=${v:-0}
  if [ "$v" -le "$SAFE_VER" ]; then echo "$m"; return; fi
  mkdir -p "$REENC_DIR"
  bn=$(basename "$m" .root); out="$REENC_DIR/${bn}_root624.root"   # keeps xyz/cartesian in name for 3D detection
  if [ ! -f "$out" ]; then
    echo "    [bridge] map is ROOT $v (> 6.26); re-encoding to native 6.24 (one-time)..." >&2
    singularity exec -B /cvmfs:/cvmfs "$IMG" bash -lc \
      "source '$SETUP' -n ana.331 >/dev/null 2>&1; root -b -q '$REPO/reencode_field.C(\"$m\",\"$out\")'" >&2
  fi
  echo "$out"
}

echo ""
echo "Select MAGNETIC FIELD:"
echo "  0) Uniform 1.5 T  (constant solenoid, fast)"
fmaps=()
declare -A _seen_fhash
i=1
while IFS= read -r -d $'\0' m; do
  bn=$(basename "$m")
  h=$(printf '%s' "$bn" | cut -c1-32)              # CDB content hash prefix -> dedup identical maps
  if [ -n "${_seen_fhash[$h]:-}" ]; then continue; fi
  _seen_fhash[$h]=1
  case "$bn" in
    *steel*)            lbl="OHCAL steel return-field  (NOT the tracking field)";;
    *gap_rebuild*)      lbl="3D solenoid, gap-rebuild   (PRODUCTION default, recommended)";;
    *bigmapxyz.root)    lbl="3D solenoid, big map";;
    *trackingmapxyz*)   lbl="3D tracking map";;
    *cartesian*)        lbl="Measured 3D map  (PARTIAL: B=0 in unmeasured regions)";;
    *measured_smoothed*) lbl="3D solenoid, measured+smoothed  (yaw 2.2 mrad, 2022-12-02)";;
    *)                  lbl="3D field map";;
  esac
  v=$(field_version "$m"); v=${v:-0}
  flag=""; [ "$v" -gt "$SAFE_VER" ] && flag="  [ROOT $v -> auto re-encode to 6.24]"
  fmaps+=("$m")
  printf "  %d) %s%s\n        %s\n" "$i" "$lbl" "$flag" "$m"
  i=$((i+1))
done < <(find "$REPO"/CDB_offline/FIELDMAP* -name "*.root" -print0 2>/dev/null | sort -z)
echo "  c) Custom value (T) or .root path"
read -p "Field choice [0]: " FCHOICE
FCHOICE=${FCHOICE:-0}
if [ "$FCHOICE" == "0" ]; then
  MAG_FIELD="1.5"
elif [ "$FCHOICE" == "c" ]; then
  read -e -p "Enter field value (T) or .root path: " CUSTOM
  case "$CUSTOM" in
    *.root) MAG_FIELD=$(ensure_field_usable "$CUSTOM") ;;
    *)      MAG_FIELD="$CUSTOM" ;;
  esac
elif [ "$FCHOICE" -ge 1 ] 2>/dev/null && [ "$FCHOICE" -le "${#fmaps[@]}" ]; then
  MAG_FIELD=$(ensure_field_usable "${fmaps[$((FCHOICE-1))]}")
else
  echo "Invalid choice -> uniform 1.5 T"; MAG_FIELD="1.5"
fi
export MAG_FIELD
echo "--> MAG_FIELD = ${MAG_FIELD}"

# TPC distortion axis (applied during simulation in PHG4TpcDistortion). Read by the macro.
echo ""
echo "Select TPC DISTORTION (applied during simulation):"
echo "  1) static + module-edge (merged)        [recommended]"
echo "  2) static only"
echo "  3) static + module-edge + time-ordered"
echo "  4) none (ideal TPC)"
read -p "Distortion choice [1]: " DCHOICE
case "${DCHOICE:-1}" in
  2) export DISTORTION_MODE=static ;;
  3) export DISTORTION_MODE=merged_timeordered ;;
  4) export DISTORTION_MODE=none ;;
  *) export DISTORTION_MODE=merged ;;
esac
echo "--> DISTORTION_MODE = ${DISTORTION_MODE}"

read -p "Output file name [G4sPHENIX.root]: " OUTPUT_FILE
OUTPUT_FILE=${OUTPUT_FILE:-G4sPHENIX.root}

echo "========================================="
echo "Starting singularity container for batch mode..."

if [[ ! -e /cvmfs/sphenix.opensciencegrid.org ]]; then
  echo "ERROR: CVMFS not mounted. Try: cvmfs_config probe sphenix.opensciencegrid.org" >&2
  exit 1
fi
if [[ ! -f "$MACRODIR/$MACRO" ]]; then
  echo "ERROR: $MACRODIR/$MACRO not found." >&2
  exit 1
fi

BINDS=(-B /cvmfs:/cvmfs)

# Run ROOT in batch mode (-b -q).
# ROOT_INCLUDE_PATH must include macros-offline/common so the LOCAL G4_TPC.C (distortions
# ON, merged static+module-edge map) is used instead of the CVMFS default (distortions OFF).
exec singularity exec "${BINDS[@]}" "$IMG" \
  bash -lc "source '$SETUP' -n ana.331 >/dev/null 2>&1; export ROOT_INCLUDE_PATH='$COMMON':'$MACRODIR':\$ROOT_INCLUDE_PATH; export DISTORTION_MODE='${DISTORTION_MODE}'; cd '$MACRODIR' && root -b -q '${MACRO}(${NEVENTS}, \"\", \"${OUTPUT_FILE}\")'"
