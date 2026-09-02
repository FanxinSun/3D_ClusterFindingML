#!/usr/bin/env bash
# =============================================================================
# setup_machine.sh — one-shot environment bootstrap for the v6.1 detached
# pipeline on a fresh machine. Linux/WSL (x86_64 or arm64) and macOS
# (Apple Silicon). Installs/builds the three external stacks, ports hardcoded
# paths if needed, writes the env file, and builds the in-repo binaries.
#
#   usage:   bash setup_machine.sh            # everything, skipping what exists
#            bash setup_machine.sh --check    # report-only, change nothing
#
# Reference versions (this is what v6.1 was produced with — newer minor
# versions work for running the chain, but md5-level reproduction of sealed
# artifacts is only expected on x86_64 Linux with the reference stack):
#   ROOT 6.40.02  ·  Geant4 11.4.2 (GDML, multithreaded)  ·  Pythia 8.317   [all C++20]
#
# Layout created (mirrors the reference machine so repo scripts run as-is):
#   $HOME/sw/root-<ver>/        ROOT install (or system/brew ROOT)
#   $HOME/geant4/               Geant4 install prefix (bin/geant4.sh)
#   <repo>/P5/angantyr/install  Pythia install
#   <repo>/env_v61.sh           source this in every working shell
#
# What this script does NOT fetch (private data — copy them yourself, see
# "payload check" at the end):
#   clusters_seeds_island_79507-0.root_ntuplizer.root       (repo root, 440M)
#   CDB_offline/FIELDMAP_GAP/65/a9/…gap_rebuild_v2.root     (299M)
#   CDB_offline/TPC_DEADCHANNELMAP/ff/c3/…TPCDeadMap_79471.root
#   island_post/raw_lib_cmflash_w.root                      (3.8M)
#   optional: island_post/*_v61.root production artifacts (see manifest)
# =============================================================================
set -eo pipefail

REPO="$(cd "$(dirname "$0")" && pwd)"
CHECK_ONLY=0
[ "$1" = "--check" ] && CHECK_ONLY=1

OS="$(uname -s)"; ARCH="$(uname -m)"
NCPU=4
if [ "$OS" = "Linux" ]; then NCPU=$(nproc); elif [ "$OS" = "Darwin" ]; then NCPU=$(sysctl -n hw.ncpu); fi
say()  { printf '\n=== %s\n' "$*"; }
have() { command -v "$1" >/dev/null 2>&1; }
run()  { if [ "$CHECK_ONLY" = 1 ]; then echo "[check] would run: $*"; else "$@"; fi; }

say "setup_machine: $OS/$ARCH, repo $REPO, ${NCPU} cores$( [ $CHECK_ONLY = 1 ] && echo ' (CHECK ONLY)')"

# ---------------------------------------------------------------- 1. OS packages
say "[1/8] OS packages"
if [ "$OS" = "Linux" ]; then
  PKGS="build-essential cmake git wget curl python3 python3-dev libxerces-c-dev \
        libx11-dev libxpm-dev libxft-dev libxext-dev libssl-dev libtbb-dev"
  if have apt-get; then
    run sudo apt-get update -y
    run sudo apt-get install -y $PKGS
  else
    echo "  non-apt Linux: install equivalents of: $PKGS"
  fi
elif [ "$OS" = "Darwin" ]; then
  have brew || { echo "  Homebrew required: https://brew.sh — install and re-run"; exit 1; }
  run brew install cmake xerces-c python coreutils gnu-sed bash wget
  echo "  note: macOS 'md5sum' is provided by coreutils (gmd5sum); env file aliases it."
else
  echo "  unsupported OS: $OS"; exit 1
fi

# ---------------------------------------------------------------- 2. ROOT
say "[2/8] ROOT (reference 6.40.02; any >=6.30 runs the chain)"
ROOT_OK=0
if have root-config; then
  echo "  found ROOT $(root-config --version) at $(command -v root-config) — keeping it"
  ROOT_OK=1
elif [ -x "$HOME/sw/root-6.40.02/bin/root-config" ]; then
  echo "  found $HOME/sw/root-6.40.02 — will use via env file"
  ROOT_OK=1
elif [ "$OS" = "Darwin" ]; then
  echo "  installing ROOT via Homebrew (current stable; fine for running the chain)"
  run brew install root && ROOT_OK=1
else
  echo "  no ROOT found. Two options (script tries neither automatically, both are long):"
  echo "    a) binary tarball for your distro from https://root.cern/install/all_releases/"
  echo "       (pick 6.40.02 if listed for your Ubuntu/gcc; untar to \$HOME/sw/root-6.40.02)"
  echo "    b) source build:"
  echo "       git clone --branch v6-40-02 --depth 1 https://github.com/root-project/root \$HOME/sw/src/root"
  echo "       cmake -S \$HOME/sw/src/root -B \$HOME/sw/build-root -DCMAKE_INSTALL_PREFIX=\$HOME/sw/root-6.40.02 \\"
  echo "             -DCMAKE_CXX_STANDARD=20 \\"
  echo "             -DCMAKE_BUILD_TYPE=RelWithDebInfo && cmake --build \$HOME/sw/build-root -j$NCPU --target install"
  echo "  re-run this script afterwards."
fi

# ---------------------------------------------------------------- 3. Geant4
say "[3/8] Geant4 11.4.2 (GDML ON, multithreaded — matches the v6.1 G4 stage)"
if [ -x "$HOME/geant4/bin/geant4.sh" ] || [ -f "$HOME/geant4/bin/geant4.sh" ]; then
  echo "  found $HOME/geant4 — keeping it"
else
  echo "  building Geant4 v11.4.2 into \$HOME/geant4 (includes ~2 GB data download)"
  if [ "$CHECK_ONLY" = 1 ]; then echo "[check] would clone+build geant4 v11.4.2"; else
    mkdir -p "$HOME/sw/src" "$HOME/sw/build-g4"
    [ -d "$HOME/sw/src/geant4" ] || git clone --branch v11.4.2 --depth 1 \
        https://github.com/Geant4/geant4.git "$HOME/sw/src/geant4"
    cmake -S "$HOME/sw/src/geant4" -B "$HOME/sw/build-g4" \
      -DCMAKE_INSTALL_PREFIX="$HOME/geant4" -DCMAKE_BUILD_TYPE=Release \
      -DGEANT4_USE_GDML=ON -DGEANT4_BUILD_MULTITHREADED=ON \
      -DGEANT4_BUILD_CXXSTD=20 \
      -DGEANT4_INSTALL_DATA=ON -DGEANT4_USE_SYSTEM_EXPAT=OFF
    cmake --build "$HOME/sw/build-g4" -j"$NCPU" --target install
  fi
fi

# ---------------------------------------------------------------- 4. Pythia
say "[4/8] Pythia 8.317 -> P5/angantyr/install"
if [ -f "$REPO/P5/angantyr/install/lib/libpythia8.a" ] || [ -f "$REPO/P5/angantyr/install/lib/libpythia8.so" ] || [ -f "$REPO/P5/angantyr/install/lib/libpythia8.dylib" ]; then
  echo "  found P5/angantyr/install — test-linking (arch ground truth)"
  if echo 'int main(){return 0;}' | g++ -x c++ - -L "$REPO/P5/angantyr/install/lib" -lpythia8 -ldl -o /tmp/.pythia_linktest 2>/dev/null; then
    rm -f /tmp/.pythia_linktest; echo "  link OK — keeping it"
  else
    echo "  link FAILED (wrong architecture or broken install) — rebuilding"
    run rm -rf "$REPO/P5/angantyr/install"
  fi
fi
if [ ! -d "$REPO/P5/angantyr/install" ]; then
  if [ "$CHECK_ONLY" = 1 ]; then echo "[check] would download+build pythia8317"; else
    cd "$REPO/P5/angantyr"
    [ -f pythia8317.tgz ] || wget -q https://pythia.org/download/pythia83/pythia8317.tgz \
      || { echo "  download failed — fetch pythia8317.tgz from https://pythia.org/releases/ manually into P5/angantyr/ and re-run"; exit 1; }
    tar xf pythia8317.tgz
    cd pythia8317 && ./configure --prefix="$REPO/P5/angantyr/install" \
      --cxx-common="-O2 -fPIC -pthread -std=c++20" && make -j"$NCPU" install
    cd "$REPO"
  fi
fi

# ---------------------------------------------------------------- 5. path port
say "[5/8] hardcoded-path port (repo scripts reference /home/rog/...)"
if [ "$HOME" = "/home/rog" ]; then
  echo "  HOME is /home/rog — nothing to port"
elif [ -f "$REPO/.paths_ported" ]; then
  echo "  already ported to $(cat "$REPO/.paths_ported") — skipping"
else
  N=$(grep -rl "/home/rog" "$REPO" \
        --include='*.sh' --include='*.C' --include='*.h' --include='*.py' --include='*.txt' --include='*.md' \
        2>/dev/null | grep -v "^$REPO/.git" | wc -l | tr -d ' ')
  echo "  porting /home/rog -> $HOME in $N tracked text files (originals: git is the backup)"
  if [ "$CHECK_ONLY" = 0 ]; then
    grep -rl "/home/rog" "$REPO" \
        --include='*.sh' --include='*.C' --include='*.h' --include='*.py' --include='*.txt' --include='*.md' \
        2>/dev/null | grep -v "^$REPO/.git" | while read -r f; do
      sed -i.pathbak "s|/home/rog|$HOME|g" "$f" && rm -f "$f.pathbak"
    done
    echo "$HOME" > "$REPO/.paths_ported"
    echo "  done. NOTE: this dirties the working tree relative to origin by design;"
    echo "  do not commit the ported paths back upstream."
  fi
fi

# ---------------------------------------------------------------- 6. env file
say "[6/8] env file -> $REPO/env_v61.sh"
if [ "$CHECK_ONLY" = 0 ]; then
  {
    echo "# source this in every shell that works on the pipeline"
    if have root-config && [ ! -x "$HOME/sw/root-6.40.02/bin/root-config" ]; then
      echo "# ROOT: system/brew install ($(root-config --version 2>/dev/null))"
    else
      echo "source \$HOME/sw/root-6.40.02/bin/thisroot.sh"
    fi
    echo "source \$HOME/geant4/bin/geant4.sh"
    echo "export PYTHIA8="$REPO"/P5/angantyr/install"
    if [ "$OS" = "Darwin" ]; then
      echo "command -v md5sum >/dev/null 2>&1 || alias md5sum='gmd5sum'"
      echo "command -v nproc  >/dev/null 2>&1 || alias nproc='sysctl -n hw.ncpu'"
    fi
  } > "$REPO/env_v61.sh"
  echo "  written."
fi

# ---------------------------------------------------------------- 7. in-repo builds
say "[7/8] in-repo binaries (standalone_tpc, gen_pp)"
if [ "$CHECK_ONLY" = 0 ] && [ -f "$HOME/geant4/bin/geant4.sh" ]; then
  ( cd "$REPO/P5" && bash build_standalone.sh ) || echo "  standalone_tpc build FAILED — check ROOT/Geant4/xerces above"
  ( cd "$REPO/P5/angantyr" && g++ -O2 -std=c++$(root-config --cxxstandard) gen_pp.cc -o gen_pp \
      -I install/include -L install/lib -lpythia8 -ldl \
      -Wl,-rpath,"$REPO/P5/angantyr/install/lib" ) \
    && echo "  gen_pp built" || echo "  gen_pp build FAILED"
else
  echo "  skipped (check-only or Geant4 missing)"
fi
echo "  note: all ROOT macros (.C) recompile themselves via ACLiC on first use —"
echo "  any *_C.so shipped from another machine is ignored/rebuilt automatically."

# ---------------------------------------------------------------- 8. payload check
say "[8/8] non-git payload check (copy these from the source machine)"
ok=1
for f in \
  "clusters_seeds_island_79507-0.root_ntuplizer.root" \
  "CDB_offline/FIELDMAP_GAP/65/a9/65a930ed6de9c0e049cd0f3ef226e6b4_sphenix3dbigmapxyz_gap_rebuild_v2.root" \
  "CDB_offline/TPC_DEADCHANNELMAP/ff/c3/ffc3f6498934c5a8ba31065292c6ebcc_TPCDeadMap_79471.root" \
  "island_post/raw_lib_cmflash_w.root"; do
  if [ -f "$REPO/$f" ]; then echo "  OK      $f"; else echo "  MISSING $f"; ok=0; fi
done
for f in digi island91; do
  p="$REPO/island_post/${f}_frames_production_v61.root"
  [ -f "$p" ] && echo "  OK      ${f}_frames_production_v61.root (analysis tier)" \
             || echo "  absent  ${f}_frames_production_v61.root (optional: ship, or re-run 'pp_pipeline.sh all')"
done
[ $ok = 1 ] || echo "  -> the chain cannot run until the MISSING files are copied in."

say "DONE. Next:  source $REPO/env_v61.sh   then verify with:"
echo "  root -b -q -e 'printf(\"ROOT %s\\n\", gROOT->GetVersion());'"
echo "  \$HOME/geant4/bin/geant4-config --version 2>/dev/null || true"
echo "  cd island_post && md5sum *_v61.root   # compare against production_manifest.txt"
echo "  (shipped-artifact md5s are only expected to match on x86_64 Linux; a fresh"
echo "   gen+g4 on Apple Silicon is a new, statistically-equivalent realization)"
