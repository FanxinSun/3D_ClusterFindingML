# sPHENIX CDB Offline Mirror

Offline copy of all **177 CDB payload domains** from the publicly accessible CVMFS mirror.

## Source

```
/cvmfs/sphenix.opensciencegrid.org/calibrations/sphnxpro/cdb/
```

**Total**: ~205,335 files, ~21 GB across 177 domains.

## Key Insight

The CDB **lookup service** (run → payload mapping) requires BNL access, but the
**payload files themselves are open** on CVMFS. This mirror allows offline use
without the CDB lookup — you just need to know which file to use.

---

## Quick Start: Using Tracking Alignment

```bash
# Copy the latest alignment file to your working directory:
cp CDB_offline/TRACKINGALIGNMENT/46/d0/46d02ad977567c26885ee9c9e9539c76_trackingalignment-ana537_2025p010-ana532.txt \
   ./localAlignmentParamsFile.txt

# AlignmentTransformation::createMap() will automatically pick it up
# before trying the CDB lookup.
```

---

## Domain Inventory by Subsystem

### Tracking — Global (36 files, ~12.5 MB)

| Domain | Files | Size | Description |
|--------|-------|------|-------------|
| TRACKINGALIGNMENT | 15 | 1.5M | Per-sensor Millepede alignment transforms (text) |
| Tracking_Geometry | 21 | 11M | TGeo geometry configs |

### Tracking — MVTX (8 files, ~4.8 MB)

| Domain | Files | Size | Description |
|--------|-------|------|-------------|
| MVTX_ALIGNMENT | 1 | 16K | MVTX-specific alignment |
| MVTX_DeadPixelMap | 3 | 56K | Dead pixel masks |
| MVTX_HotPixelMap | 4 | 4.7M | Hot pixel masks (ROOT) |

### Tracking — INTT (11,071 files, ~537 MB)

| Domain | Files | Size | Description |
|--------|-------|------|-------------|
| INTT_BCOMAP | 4,018 | 50M | BCO timing maps per run |
| INTT_HotMap | 7,040 | 486M | Hot channel masks per run |
| INTT_HotMap_default | 3 | 67K | Default hot channel mask |
| INTT_DACMAP | 4 | 220K | DAC calibration maps |
| InttFeeMap | 1 | 22K | FEE channel mapping |
| InttSurveyMap | 2 | 93K | Survey geometry |
| INTT_Hotmap | 1 | 36K | Legacy hot map |
| INTT_STREAMING_FEE_OFFSET | 2 | 34K | Streaming FEE timing offsets |

### Tracking — TPC (3,699 files, ~399 MB)

| Domain | Files | Size | Description |
|--------|-------|------|-------------|
| TPC_DEADCHANNELMAP | 642 | 63M | Dead channel masks per run |
| TPC_HOTCHANNELMAP | 642 | 27M | Hot channel masks per run |
| TPC_DRIFT_VELOCITY | 1,051 | 12M | Drift velocity per run |
| TPC_DECODER_BAD_FEE | 615 | 7.7M | Bad FEE decoder masks per run |
| TPC_LAMINATION_FIT_CORRECTION | 521 | 135M | Lamination distortion corrections |
| TPC_TZERO_OFFSET | 208 | 2.7M | T-zero timing offsets per run |
| TPC_STATIC_DISTORTION | 1 | 2.4M | Static distortion map (ROOT hist) |
| TPC_STATIC_CORRECTION_MODEL | 1 | 17M | Static distortion correction model |
| TPC_TIMEORDERED_DISTORTION | 1 | 120M | Time-ordered distortion map |
| TPC_FEE_CHANNEL_MAP | 6 | 677K | FEE-to-channel mapping |
| TPC_ZS_THRESHOLDS | 2 | 35K | Zero-suppression thresholds |
| TPC_Module_Edge | 2 | 165K | Module edge definitions |
| TPC_DEDX_FITPARAM | 1 | 55K | dE/dx fit parameters |
| TPC_CAEN_CORRECTION_24974_* | 3 | 1.5M | CAEN digitizer corrections |
| TPC_DRIFT_VELOCITY_default | 1 | 20K | Default drift velocity |
| TPC_LAMINATION_FIT_CORRECTION_default | 1 | 323K | Default lamination correction |
| TPC_TZERO_OFFSET_default | 1 | 20K | Default T-zero offset |
| PHGARFIELD_GAS | 1 | 9.9M | Magboltz gas table (Ar75/CF4-20/iso-5); **not read by the sim** — parametric, see note below |

### Tracking — TPOT/Micromegas (1 file, 80 KB)

| Domain | Files | Size | Description |
|--------|-------|------|-------------|
| TPOT_Pedestal | 1 | 80K | Pedestal values |

### Tracking — Field Maps (6 files, ~877 MB)

| Domain | Files | Size | Description |
|--------|-------|------|-------------|
| FIELDMAPBIG | 1 | 345M | Full magnetic field map |
| FIELDMAP_GAP | 1 | 279M | Gap field map |
| FIELDMAP_STEEL | 1 | 207M | Steel field map |
| FIELDMAPTRACKING | 1 | 19M | Tracking-region field map |
| FIELDMAP_TRACKING | 2 | 27M | Tracking field map (alt) |

### Acts Geometry (6 files, ~1.6 MB)

| Domain | Files | Size | Description |
|--------|-------|------|-------------|
| ACTSGEOMETRYCONFIG | 4 | 77K | TGeo → Acts JSON configs |
| ACTSMATERIALMAP | 2 | 1.5M | Acts material map JSON |

### Calorimeter — CEMC (~70K files, ~6 GB)

| Domain | Files | Size | Description |
|--------|-------|------|-------------|
| CEMC_BadTowerMap | 15,228 | 1.6G | Bad tower masks per run |
| CEMC_meanTime | 17,937 | 2.3G | Mean time calibration per run |
| CEMC_hotTowers_fracBadChi2 | 14,634 | 791M | Hot tower masks per run |
| CEMC_ZSCrossCalib | 6,064 | 722M | ZS cross-calibration per run |
| *(+ 30 smaller CEMC domains)* | | | |

### Calorimeter — HCAL (~75K files, ~6 GB)

| Domain | Files | Size | Description |
|--------|-------|------|-------------|
| HCALIN_hotTowers_fracBadChi2 | 16,414 | 236M | Inner HCAL hot towers per run |
| HCALIN_meanTime | 17,975 | 306M | Inner HCAL mean time per run |
| HCALOUT_hotTowers_fracBadChi2 | 16,414 | 450M | Outer HCAL hot towers per run |
| HCALOUT_meanTime | 17,975 | 556M | Outer HCAL mean time per run |
| *(+ many smaller HCAL domains)* | | | |

### Event Characterization — MBD/Centrality (~40K files, ~600 MB)

| Domain | Files | Size | Description |
|--------|-------|------|-------------|
| MBD_PED | 9,002 | 136M | MBD pedestal per run |
| MBD_PILEUP | 8,705 | 135M | MBD pileup correction |
| MBD_SAMPMAX | 9,030 | 100M | MBD sample max |
| MBD_T0CORR | 8,747 | 110M | MBD T0 correction |
| *(+ smaller MBD/Centrality domains)* | | | |

### Forward Detector — sEPD (~1.3K files, ~50 MB)

| Domain | Files | Size | Description |
|--------|-------|------|-------------|
| SEPD_EventPlaneCalib | 1,199 | 48M | sEPD event plane calibration |
| *(+ smaller sEPD domains)* | | | |

### Other (~140 files, ~11 GB)

| Domain | Files | Size | Description |
|--------|-------|------|-------------|
| EVENTPLANE_CORRECTION | 129 | 11G | Event plane correction (largest single domain!) |
| EVENTPLANE_CORRECTION_default | 2 | 120M | Default event plane correction |
| *(+ misc domains)* | | | |

---

## TRACKINGALIGNMENT Files — Detailed Timeline

| File | Date | Size | Notes |
|------|------|------|-------|
| `tracking_alignment_parameters_empty.txt` | 2023-04 | 37K | Zeros/baseline |
| `alignmentParamsFile_cdb_june28_2024.txt` | 2024-06 | 50K | First Run 2024 alignment |
| `AlignmentParamsFile_silicononly_july26_norotations.txt` | 2024-07 | 61K | Silicon-only, no rotations |
| `modified_version2_cdb_localAlignmentParamsFile.txt` | 2024-11 | 52K | Modified v2 |
| `mvtx_cosmic_interal_intt_realign_localAlignmentParamsFile.txt` | 2025-01 | 64K | MVTX cosmic + INTT realign |
| `silicon_tpot_tpc_offset_tpc_layers_localAlignmentParamsFile.txt` | 2025-02 | 98K | Silicon + TPOT + TPC |
| `Feb20_TPC_zshift_update_localAlignmentParamsFile.txt` | 2025-02 | 101K | TPC z-shift update |
| `Mar11_tpotupdate_localAlignmentParamsFile.txt` | 2025-03 | 101K | TPOT update |
| `shiftedMvtxGeomAlignmentParamsFile.txt` | 2025-03 | 103K | Shifted MVTX geometry |
| **`trackingalignment-ana537_2025p010-ana532.txt`** | **2025-03** | **120K** | **Latest: ana537 prod** |
| `tracking_alignment_ana494_2024p021_v001.txt` | 2025-07 | 113K | ana494 2024p021 v1 |
| `tracking_alignment_inttupdate_ana494_2024p021_v001.txt` | 2025-08 | 113K | + INTT update |
| `shiftedMvtxGeomAlignmentParamsFileInttUpdate.txt` | 2025-08 | 103K | Shifted MVTX + INTT |
| `trackingalignment-ana494_2024p022-ana500.txt` | 2025-09 | 114K | ana494 2024p022 |
| `trackingalignment-ana494_2024p022-ana516.txt` | 2025-10 | 113K | ana494 2024p022 update |

### Alignment File Format

Each line: `hitsetkey  alpha  beta  gamma  dx  dy  dz`

```
16 -0.000547 -0.002843 -0.001761  0.365910 -0.605773 -0.137230
48 -0.001986 -0.001279  0.000352  0.311148 -0.559624 -0.143669
80 -0.002220 -0.001339  0.001344  0.269082 -0.499655 -0.150248
```

- Columns 2-4: Millepede rotation angles (alpha, beta, gamma) in radians
- Columns 5-7: Translation offsets (dx, dy, dz) in cm
- Optional columns 8-10: Global rotation corrections (dgrx, dgry, dgrz)

---

## PHGARFIELD_GAS — Magboltz gas table is **not** consumed by the sim

*Research record, 2026-06-29. Scanned: local `coresoftware` @ `production_A_01-483-g21f9fed50`;
CVMFS `ana.220`–`ana.331` (gcc-8.3) and `ana.333`–`ana.479` + `new.*` (gcc-12.1.0).*

**Domain:** `PHGARFIELD_GAS/a8/51/a8516e920f926da938f7f138aed548e6_Ar75_CF20_iso5.gas`
(1 file, 9.9 MB ASCII Magboltz table — `CF4 20%, Ar 75%, iC4H10 5%, T=301.65 K, p=1.003 atm`).

**Verdict:** `run_batch_sim.sh` does **not** open this file, and no sPHENIX release does either.
sPHENIX consumes gas physics **parametrically**: the transport properties are extracted from the
gas table offline (once) and entered as parameters. This `.gas` payload is a **provenance/archive
artifact** — like the 1,051 `TPC_DRIFT_VELOCITY` files, it is present in the mirror but not read
by the simulation chain.

**Where the gas numbers actually live** (`macros-offline/common/G4_TPC.C`, matching the upstream
`PHG4TpcElectronDrift.cc` defaults for the `Ar75_CF20_iso5` mixture):

| Quantity | Value | Location |
|----------|-------|----------|
| drift velocity (sim) | 0.007550 cm/ns | `G4_TPC.C:69` |
| longitudinal diffusion | 0.014596 cm/√cm | `G4_TPC.C:118` |
| transverse diffusion | 0.005313 cm/√cm | `G4_TPC.C:119` |
| Ar / CF4 / isobutane fractions | 0.75 / 0.20 / 0.05 | `PHG4TpcElectronDrift.cc:943-946` |

The `*_frac` params feed only a dE/dx primary-ionization sum (`PHG4TpcElectronDrift.cc:267-301`),
never a gas-table read.

**Evidence (each independently conclusive):**

1. The only `.gas` reader in `coresoftware` is `offline/packages/PHGarfield/` (`LoadGasFile`/
   `MergeGasFile`). `libPHGarfield.so` is built in **zero** releases. It is a private dev tool:
   not registered in any macro, loads `PART_<i>.gas` from a hardcoded BNL path
   (`/direct/phenix+u/.../hemmick/.../gasfiles/`, not the CDB), and its `process_event()` is a stub.
2. Garfield++/`MediumMagboltz` is **absent** from both external stacks (gcc-8.3, gcc-12.1.0) — no
   release could load a gas table even if asked.
3. The string `PHGARFIELD_GAS` appears nowhere in `coresoftware`; nothing calls
   `getUrl("PHGARFIELD_GAS")`. `ana.479` key libs (`libg4tpc/libtpc/libffamodules/libcdbobjects`)
   and TPC headers contain no gas-table API.
4. `run_batch_sim.sh` references only `TRACKINGALIGNMENT`, `FIELDMAP*`, and the TPC distortion maps.

**Cross-release evolution (`ana.331` → `ana.479`):** method unchanged
(gas mixture → offline Magboltz, once → parametric constants → sim). The only new knob is
`use_PDG_gas_params()` (since **`ana.432`**, `PHG4TpcElectronDrift.cc:285-295`), which swaps the
per-gas dE/dx + primary-electron constants to PDG-2024 values feeding `electrons_per_gev`.
Diffusion and drift velocity remain macro-set parameters. This affects signal amplitude / hit
count only — **not** the distorted-hit geometry the distortion study targets — so it is cosmetic
here. Upgrading the image from `ana.331` would **not** make the `.gas` file get consumed.

---

## Directory Structure

```
CDB_offline/
├── README.md                    ← This file
├── cdb_manifest.yaml            ← Full YAML inventory (822K lines)
├── copy_cdb.sh                  ← Script to refresh from CVMFS
├── generate_manifest.sh         ← Script to regenerate manifest
│
├── TRACKINGALIGNMENT/           ← Per-sensor alignment transforms
│   ├── 46/d0/..._ana537_...txt  ← Latest production alignment
│   ├── ef/04/..._june28_...txt  ← First 2024 alignment
│   └── ...
├── TPC_STATIC_DISTORTION/       ← TPC static distortion map
│   └── 5a/3d/...root            ← 2.4 MB ROOT histogram
├── TPC_DEADCHANNELMAP/          ← Per-run TPC dead channel maps
│   ├── e6/70/...79394.root
│   ├── e6/44/...79509.root      ← Near your run 79507
│   └── ...
├── CEMC_BadTowerMap/            ← Per-run CEMC bad tower maps
│   └── (15,228 files)
└── ... (177 domains total)
```

---

## Files in This Directory

| File | Purpose |
|------|---------|
| `README.md` | This documentation |
| `cdb_manifest.yaml` | Complete YAML inventory of all domains and files |
| `copy_cdb.sh` | Shell script to copy/refresh all domains from CVMFS |
| `generate_manifest.sh` | Shell script to regenerate the YAML manifest |
