#!/usr/bin/env bash
#
# p5_display.sh — interactive event display for the P5 pp chain (v5.5+).
#
# This is the visualisation counterpart of pp_pipeline.sh, and the replacement
# for sphenix_display.sh on this chain. The difference is not cosmetic:
#
#   sphenix_display.sh   CVMFS + singularity + a container Fun4All macro.
#                        Knows nothing about P5, standalone_tpc, or any batch
#                        product — it re-generates its own single-particle
#                        events inside the container. Still valid for the
#                        container reference chain; nothing here touches it.
#
#   p5_display.sh        host toolchain (local Geant4 11.4.2 + local ROOT), the
#                        SAME sphenix_p5.gdml the [g4] stage loads, and the
#                        batch FINAL OUTPUT drawn on top of it. No container.
#
# Two overlay modes, mutually exclusive because the pipeline makes them so:
#
#   frame     island_post/digi_frames_production_<VER>.root      pixels
#             island_post/island91_frames_production_<VER>.root  clusters+truth
#             truth cluster chains grouped by ntp_truth gtrackID = "tracks"
#             (this chain has no reconstructed tracker; every track drawn is a
#              truth trajectory, which is exactly what the ML labels are made of)
#
#   library   P5/PP_g4hit_<chunk>.root  raw Geant4 steps + per-track polylines
#             for ONE library collision. A production frame is a pileup draw
#             over many library collisions and frame_composer stores only the
#             remapped newid, not the source (file,event) — so a frame cannot be
#             traced back to its g4hits, and the display does not pretend it can.
#
# Geant4 stays live in the GUI: /run/beamOn works, on the pipeline's own HepMC
# chunk if you pick one, with the pipeline's field map.
#
# usage:  ./p5_display.sh                 interactive menu
#         ./p5_display.sh --frame 12      straight to a frame (no prompts)
#         ./p5_display.sh --g4event 7 --g4chunk 3
#         ./p5_display.sh --ver v55 --frame 3 --sector 5 --colour class
#         ./p5_display.sh --help          all p5_display options
#
# FULL REFERENCE (every flag + all 34 /p5/ commands): P5/DISPLAY.md
#
# Everything reachable from the session prompt is reachable from the shell too:
#         --geometry tracker --geomstyle surface 0.10 --view rphi
#         --cmd "/p5/show pix false" --cmd "/p5/trk/select 5803"
# --cmd takes ANY /p5/ or /vis/ command, is repeatable, and runs in the order
# given once the viewer exists — so a flag is only worth adding for the ones you
# type constantly. Handy for scripted stills:
#         ./p5_display.sh --batch --vis TSG_OFFSCREEN --frame 12 \
#             --view rphi --cmd "/vis/viewer/flush"
#
# Any option not listed below is passed through to the binary untouched.
set -eo pipefail

ROOTDIR=/home/rog/sPHENIX/3D_ClusterFindingML
P5=$ROOTDIR/P5
IP=$ROOTDIR/island_post
BIN=$P5/p5_display
SRC=$P5/p5_display.cc
G4ENV=/home/rog/geant4/bin/geant4.sh
DEFAULT_VER=v55
FM=$ROOTDIR/CDB_offline/FIELDMAP_GAP/65/a9/65a930ed6de9c0e049cd0f3ef226e6b4_sphenix3dbigmapxyz_gap_rebuild_v2.root

# ---------------- environment ----------------
[ -f "$G4ENV" ] || { echo "ERROR: local Geant4 not found at $G4ENV" >&2; exit 1; }
# geant4.sh is not nounset-clean — same reason pp_pipeline.sh avoids `set -u`
source "$G4ENV"
command -v root-config >/dev/null || { echo "ERROR: local ROOT not on PATH" >&2; exit 1; }

WANT_X=1
case " $* " in *" --batch "*) WANT_X=0;; esac
if [ "$WANT_X" = 1 ] && [ -z "${DISPLAY:-}" ] && [ -z "${WAYLAND_DISPLAY:-}" ]; then
  echo "ERROR: no DISPLAY. Under WSL2 this normally comes from WSLg; over ssh use -X." >&2
  echo "       For a windowless render instead:" >&2
  echo "         $0 --batch --vis TSG_OFFSCREEN --frame 0" >&2
  exit 1
fi

# ---------------- build if stale ----------------
if [ ! -x "$BIN" ] || [ "$SRC" -nt "$BIN" ]; then
  echo "[build] p5_display.cc is newer than the binary — rebuilding"
  "$P5/build_display.sh"
fi

# ---------------- non-interactive passthrough ----------------
if [ $# -gt 0 ]; then
  exec "$BIN" --repo "$ROOTDIR" "$@"
fi

# ---------------- interactive menu ----------------
echo "========================================================="
echo "   P5 event display — GDML geometry + batch final output "
echo "========================================================="

# versions actually present on disk, newest first
mapfile -t VERS < <(ls -1t "$IP"/digi_frames_production_*.root 2>/dev/null |
                    sed -E 's|.*digi_frames_production_(.*)\.root|\1|')
if [ ${#VERS[@]} -eq 0 ]; then
  echo "ERROR: no digi_frames_production_*.root in $IP — run ./pp_pipeline.sh prod first" >&2
  exit 1
fi
echo "Production versions found:"
for i in "${!VERS[@]}"; do
  v=${VERS[$i]}
  isl=$IP/island91_frames_production_$v.root
  printf "  %d) %-6s  digi %s%s\n" "$((i+1))" "$v" \
      "$(du -h "$IP/digi_frames_production_$v.root" | cut -f1)" \
      "$([ -f "$isl" ] && echo "  + island91 clusters" || echo "  (NO island91 clusters — pixels only)")"
done
read -e -p "Version [1 = ${VERS[0]}]: " VC
VC=${VC:-1}
VER=${VERS[$((VC-1))]:-$DEFAULT_VER}

echo ""
echo "Overlay mode:"
echo "  1) production frame  — digitised pixels + island91 clusters + truth chains"
echo "  2) library G4 event  — raw Geant4 steps + G4 tracks for one collision"
read -e -p "Mode [1]: " MODE
MODE=${MODE:-1}

ARGS=(--repo "$ROOTDIR" --ver "$VER")

if [ "$MODE" == "2" ]; then
  read -e -p "PP_g4hit chunk [0]: " CH
  CH=${CH:-0}
  # One pp minimum-bias collision is sparse by physics: the median library event
  # leaves ~3k G4 steps in the TPC gas and the spread is an order of magnitude.
  # If this chunk has been loaded before, its cached range index already knows
  # the per-event counts, so offer the busy ones instead of making you guess.
  GIDX=$P5/.p5disp/PP_g4hit_$CH.root.ntp_g4hit.idx
  if [ -f "$GIDX" ]; then
    python3 - "$GIDX" <<'EOP'
import sys, statistics
f = open(sys.argv[1]); f.readline(); nev = int(f.readline())
per = {}
for _ in range(nev):
    t = f.readline().split()
    per[int(t[0])] = sum(int(t[3 + 2 * k]) for k in range(int(t[1])))
v = sorted(per.values())
top = sorted(per.items(), key=lambda kv: -kv[1])[:8]
print(f"  this chunk: {len(v)} events with TPC hits, median {statistics.median(v):.0f} g4 hits")
print("  busiest:  " + "  ".join(f"{e}({n})" for e, n in top))
EOP
  else
    echo "  (per-event counts appear here once this chunk has been indexed;"
    echo "   in-session: /p5/g4list 10)"
  fi
  echo "  index 0..1999. A RANGE or LIST superimposes several collisions:"
  echo "    0-9      ten collisions at once      2,860,1354   these three"
  read -e -p "Library event index / range [0]: " EV
  ARGS+=(--g4chunk "$CH" --g4event "${EV:-0}")
  echo ""
  echo "Note: every library collision is generated at the origin (gen_pp runs with"
  echo "      vertex spread 0), so a range stacks them all on the same vertex. It is"
  echo "      N collisions superimposed, NOT a production frame — mode 1 is that,"
  echo "      with drift time, z spread and digitisation."
  read -e -p "G4 hit mark size in px [2.5]: " PS
  [ -n "$PS" ] && ARGS+=(--pixsize "$PS")
else
  # frame count comes free from the display's own range-index cache if it has
  # been built; asking ROOT for it would cost a full 83M-row branch scan
  IDX=$IP/.p5disp/digi_frames_production_$VER.root.ntp_hit.idx
  NEV=$([ -f "$IDX" ] && sed -n '2p' "$IDX")
  read -e -p "Frame index [0${NEV:+, 0..$((NEV-1))}]: " FR
  ARGS+=(--frame "${FR:-0}")
  read -e -p "TPC sector 0-11, blank = ALL 12: " SEC
  [ -n "$SEC" ] && ARGS+=(--sector "$SEC")
  echo "Colour mode:  1) level (hit-level pixels vs cluster-level vs truth chain)"
  echo "              2) class (ML truth label 0/1/2)      3) track (spotlight the leading chains)"
  read -e -p "Colour [1]: " CM
  case "${CM:-1}" in 2) ARGS+=(--colour class);; 3) ARGS+=(--colour track);; esac
  # A frame is ~200k pixels over the whole TPC. No colour scheme survives that
  # at full-detector zoom; the fix is to look at less of it. An ADC floor is the
  # single most effective knob (>=60 keeps roughly 1 pixel in 5).
  echo ""
  echo "Density: a frame is ~200k pixels / ~20k clusters. Blank keeps all of it;"
  echo "         an ADC floor thins the low-amplitude tail (60 keeps about 1 in 5)."
  read -e -p "Pixel ADC floor [none]: " AF
  [ -n "$AF" ] && ARGS+=(--adcmin "$AF")
fi

echo ""
echo "Geometry from the GDML:"
echo "  1) tpc      gas + field cage only (fastest, data unobstructed)      23 vol"
echo "  2) tracker  + TPC endcap, MVTX, INTT, TPOT, beam pipe              754 vol"
echo "  3) calo     + CEMC, HCal/cryostat, EPD                             794 vol"
echo "  4) all      + beam line: everything the export carries            1020 vol"
echo "     each system gets its own colour and a legend row (/p5/geom/tint false = grey)"
read -e -p "Geometry [1]: " GM
case "${GM:-1}" in 2) ARGS+=(--geometry tracker);; 3) ARGS+=(--geometry calo);; 4) ARGS+=(--geometry all);; esac
if [ "$MODE" != "2" ] && [ "${GM:-1}" != "1" ]; then
  echo "  note: same colours as mode 2. Over a full frame that competes with the"
  echo "        overlay — thin the data (ADC floor / sector) or /p5/geom/bright 0.5."
fi

echo ""
read -e -p "Wire up /run/beamOn with a HepMC chunk? (index 0-9, blank = gun only): " HM
if [ -n "$HM" ]; then
  H=$P5/angantyr/pp_run_$HM.dat
  [ -f "$H" ] || { echo "  $H not found — skipping"; H=""; }
  [ -n "$H" ] && ARGS+=(--hepmc "$H" --field "$FM")
fi

echo ""
echo "Viewer:  1) TSG  (ToolsSG + Qt, the Geant4 11 default)   2) OGLSQt  (classic OpenGL + Qt)"
echo "         3) OGLSX (classic OpenGL, plain X11, no Qt GUI panel)"
read -e -p "Viewer [1]: " VW
case "${VW:-1}" in
  2) ARGS+=(--vis OGLSQt);;
  3) ARGS+=(--vis OGLSX --session tcsh);;
  *) ARGS+=(--vis TSG);;
esac

echo ""
echo "launching: $BIN ${ARGS[*]}"
echo "  in the session:  /p5/print   /p5/frame <n>   /p5/view rphi|rz|3d|iso"
echo "                   /p5/sector <0-11>  (-1 = all 12)   /p5/zrange full"
echo "                   /p5/colour level|class|track   /p5/show pix|clus|trk <bool>"
echo "                   /p5/geometry tpc|tracker|calo|all   /p5/geom/show <name>"
echo "                   /p5/geom/style wireframe|surface [alpha]   /p5/geom/bright <f>"
echo "                   /p5/geom/tint <bool>"
echo "                   /p5/zoom <radius_cm>   help /p5/"
echo ""
exec "$BIN" "${ARGS[@]}"
