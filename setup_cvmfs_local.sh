#!/usr/bin/env bash
#
# setup_cvmfs_local.sh
#
# Replicates the proven sPHENIX CVMFS setup from the development VM
# (yaminocellist@192.168.255.129, Ubuntu 22.04, CVMFS 2.13.3) onto this
# local WSL2 box (Ubuntu 26.04 "resolute").
#
#   - CVMFS 2.13.3 + cvmfs-config-default 2.2-1   (CERN resolute-prod repo: same
#     version as the VM, natively built for 26.04)
#   - Apptainer from Ubuntu 26.04 universe        (provides the `singularity` cmd)
#   - autofs-mounted /cvmfs                        (systemd + autofs both present)
#   - openhtc.io Cloudflare CDN Stratum-1 pin      (copied verbatim from the VM)
#
# Run ONCE with root. Idempotent: safe to re-run.
#     sudo bash setup_cvmfs_local.sh
#
set -euo pipefail

if [[ ${EUID} -ne 0 ]]; then
  echo "Please run as root:  sudo bash $0" >&2
  exit 1
fi

echo "==> [1/6] Add CERN cvmfs-release apt repo (selects resolute-prod for 26.04)"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
curl -fsSL -o "$TMP/cvmfs-release.deb" \
  https://ecsft.cern.ch/dist/cvmfs/cvmfs-release/cvmfs-release-latest_all.deb
dpkg -i "$TMP/cvmfs-release.deb"

echo "==> [2/6] Install CVMFS 2.13.3 (CERN) + autofs + Apptainer (26.04 universe)"
apt-get update
DEBIAN_FRONTEND=noninteractive apt-get install -y \
  cvmfs cvmfs-config-default autofs apptainer

echo "==> [3/6] Write /etc/cvmfs/default.local"
# Mirrors the VM; adds sphenix to CVMFS_REPOSITORIES since this box exists to run it.
cat > /etc/cvmfs/default.local <<'CONF'
CVMFS_STRICT_MOUNT=no
CVMFS_REPOSITORIES=oasis.opensciencegrid.org,sphenix.opensciencegrid.org
CVMFS_HTTP_PROXY="DIRECT"
CVMFS_QUOTA_LIMIT=50000
CONF

echo "==> [4/6] Write /etc/cvmfs/domain.d/opensciencegrid.org.local (openhtc.io CDN pin)"
mkdir -p /etc/cvmfs/domain.d
cat > /etc/cvmfs/domain.d/opensciencegrid.org.local <<'CONF'
# Prefer the Cloudflare-backed openhtc.io CDN Stratum-1s (much faster than the
# direct :8000 servers, which gave ~150 KB/s). Mirrored from the dev VM.
CVMFS_SERVER_URL="http://s1fnal-cvmfs.openhtc.io/cvmfs/@fqrn@;http://s1bnl-cvmfs.openhtc.io/cvmfs/@fqrn@;http://s1ral-cvmfs.openhtc.io/cvmfs/@fqrn@;http://s1nikhef-cvmfs.openhtc.io/cvmfs/@fqrn@"
CONF

echo "==> [5/6] cvmfs_config setup + enable autofs"
cvmfs_config setup
systemctl enable --now autofs
# Nudge autofs to pick up the freshly written /cvmfs map.
systemctl reload autofs 2>/dev/null || systemctl restart autofs

echo "==> [6/6] Self-test"
cvmfs_config probe sphenix.opensciencegrid.org
echo "--- /cvmfs/sphenix.opensciencegrid.org (should list singularity/, gcc-8.3/, ...) ---"
ls /cvmfs/sphenix.opensciencegrid.org/ | head
echo
echo "SUCCESS: CVMFS + Apptainer ready. You can now run ./sphenix_display.sh"
