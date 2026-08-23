# Request: per-region azimuthal-offset ("twist") field term — real pixels carry a
# sawtooth r-phi profile that the v6 digi field does not

2026-08-21, from the "Real Data Probe" session to the pipeline session, relayed by the
user. Follow-up to `digi_field_request.md` (v5.5, field in digitization). This request
adds ONE missing component to that field model, measured on run 79507 pixels; it is the
first SIGNED field anchor (all v5.4/v5.5 anchors — circle-RMS width, d0 RMS, CM jump —
were width-level). Nothing in V6 is invalidated; V6 gates never looked at this axis.
Full probe records: REAL_DATA_PROBES.md (repo root), entries of 2026-08-21.

## What was measured (real vs sim digi v6, nf_digipix groups, trimmed global fits)

1. Split-half curvature (parameter level, blind to the cloud floor): fit rows 7-30 and
   31-54 of each track independently; Dsagitta = (k_in - k_out) x span^2/8.
   real median **+6.5 mm** (robust sigma 9.1; side 0 +6.0, side 1 +7.2; n = 4384)
   sim  median **+0.2 mm** (sigma 5.2; n = 21458; same estimator -> unbiased)
   CHARGE-INDEPENDENT: real bend<0 +6.05 mm (n 2024), bend>0 +7.12 mm (n 2360)
   -> a coherent r-phi displacement field, not a |curvature| or population effect.
   Producer: island_post/nonrms_probe.C (nonrms_probe_v6.txt).
2. The profile itself (island_post/twist_probe.C -> twist_probe_v6.txt, figure
   sim_validation_plots/twist_probe_v6.png): mean azimuthal displacement of kept pixels
   relative to the fitted circle, D(rphi) = r * wrap(phi_pixel - phi_fit(r)), per side
   and pad row, averaged over all tracks. REAL is a SAWTOOTH locked to the module
   regions (both sides, same sign):
     R1: +769/+922 um at row 7  -> -1268/-1408 um at row 22  (side 0 / side 1)
     R2: +449/+430 um at row 23 -> -568/-561 um at row 38
     R3: +982/+1022 um at row 39 -> -692/-1011 um at row 54
   Boundary jumps (real - sim): R1->R2 +1704 / +1857 um, R2->R3 +1311 / +1796 um.
   SIM digi v6 is flat (+-150 um noise) everywhere.
   Sector coherence: R1 (rows 7-14) per-sector means all positive, +182 +- 85 um
   (side 0), +265 +- 115 um (side 1) -> a common pattern with ~100 um sector scatter.
3. Same sign on BOTH sides: an E x B / drift distortion would flip sign with the drift
   direction; this does not -> geometric class (per-region azimuthal placement of the
   pad rows relative to the real tracks, as recorded in the real reco's pixel x,y).
   It is NOT the rowdr radial ramp: the V6 real-radius bake (tpc_geom_table.txt =
   GDML + rowdr, an R1-only +1.2 -> +0.4 mm radial ramp) is already common to BOTH
   sides' pixel geometry, so it cannot produce a real-vs-sim difference; the sawtooth
   is the missing AZIMUTHAL counterpart of that same inner-region structure.

## Interpretation (discussion)

A sawtooth in fit-orthogonal residuals is what per-region CONSTANT azimuthal offsets
look like after a smooth circle fit: the fit draws a smooth curve through three step
levels, so within each region the residual slopes down and jumps back at the boundary.
A 3-level step model built from the two jumps reproduces the bulk of it but not all
(see closure) — the regions also carry an extra within-region slope. The full measured
profile is therefore the spec; the step model is a 3-number approximation.

## Payload: island_post/twist_profile_v6.txt

Columns: layer side real_um sim_um delta_um. `delta` = real - sim (sim's own
field-model residual subtracted, so injecting `delta` lands sim ON real). Sign:
positive = pixel at LARGER azimuth (counter-clockwise in x-y) than the fitted circle.
Injection as an azimuthal displacement of the transported charge, per side and pad
row, at the digitization stage (the same place the v5.5 field enters, post-diffusion,
pre-binning): dphi(side,row) = delta(side,row) * 1e-4 / r [rad], r in cm.
3-number approximation (if a table is unwelcome): per-region constant r-phi offsets
0 / +1.78 / +3.34 mm (side 0: 0/+1.70/+3.02; side 1: 0/+1.86/+3.65), mean-subtracted.

## Closure on our side (twist_probe.C, sim pixels displaced in memory)

  sim + full delta table : Dsagitta med +6.27 mm (real +6.53; sides +6.18/+6.38 vs
                           +6.04/+7.15) | profile RMS-diff to real 126 um (sim alone 447)
  sim + 3-level steps    : Dsagitta +3.8 (const r-phi) / +4.4 mm (rigid rotation) |
                           RMS-diff 327 / 311 um
=> the table closes both observables to ~4%; the step model gets ~60-70% of the way.

## The ask

1. Add the per-(side, row) azimuthal displacement of twist_profile_v6.txt to the digi
   field model (alongside SMOD/SPHI/SCM), re-digitize the V6 production, re-export.
2. Acceptance this thread will run on the new digi (island_post/twist_probe.C +
   nonrms_probe.C, unchanged): split-half Dsagitta median per side within +-1 mm of
   real (+6.0 / +7.2); profile RMS-diff to real <= ~150 um; the pooled residual core
   (q68 1.88 mm) and the clipped GLOBAL (1.03) must NOT move (the term is
   fit-orthogonal: it should leave width-level meters alone — that is the check that
   nothing else was disturbed).
3. NOT requested: any change to response/ZS/cloud knobs, the road, or the sealed V6
   anchors; no retune of SMOD/SPHI/SCM (the autocorrelation-length mismatch — real
   granular vs sim smooth, REAL_DATA_PROBES 2026-08-21 — is a separate, second-order
   item; re-measure C(d) after this injection before touching it, the per-region steps
   may account for part of it).

## Decisions left to you (discussion)

- Geometry error vs physical distortion: for the ML dataset it does not matter — the
  real pixels AS RECORDED carry this pattern, and the ML trains on recorded pixels, so
  sim pixels should carry it too (rowdr precedent). If the pipeline prefers to keep
  "reco-geometry" effects at export level, note the rowdr lesson: pixel-level meters
  and any pixel-level training miss it unless it is in the digi.
- Per-sector scatter (~100 um) is not in the table; a per-(side, sector, row) table can
  be produced on request (twist_probe.C already accumulates it).
- The profile is measured on tracker-seed tracks (in-time) and is the fit-orthogonal
  part; injecting it as-is is what the closure test did.
