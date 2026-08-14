# Request: position-level realism for the detached transport (distortion + alignment)

**From:** the sim-output analysis / ML session (558f78d5)
**To:** the sim-pipeline session, relayed by the user
**Date:** 2026-08-05
**Status:** request with measured evidence; user decides priority and seals.

## Summary

Track-level analyses of real run 79507 against the v5.3 production have produced
the first observables that *see* the coherent position corrections dropped at
detachment (2026-07-08): a **+33% real excess in per-track circle-fit residuals
on stiff tracks** and a **~2.5 cm transverse impact-parameter offset of real
tracks from the coordinate origin**. The offset-vs-distortion diagnostic has
been RUN (see below): a frame/beam translation is ruled out at the 1 mm
level, the broadening is azimuthally isotropic with significant drift-time
structure — distortion-family. Request: apply the static distortion +
alignment content at the transport stage, using the numbers below as
acceptance meters. All pixel-level anchors should be invariant — that is the
regression guard, and it is also why this was invisible until now.

## Background (known, restated for completeness)

- The container-era chain applied `DISTORTION_MODE=merged` (static +
  module-edge maps) in-drift plus the ana537 `TRACKINGALIGNMENT`; the
  detached chain (`standalone_tpc` -> `tpc_transport`) applies neither.
  Declared ledger item since the 2026-07-20 side review, which also measured
  the related **R1 radius taper**: sim cluster radii low by 1.20 mm (L7) ->
  0.91 (L12) -> 0.40 (L22), R2/R3 agreeing — the signature of a coherent
  radial effect.
- This was parked deliberately: every pixel-level and statistical acceptance
  observable is blind to coherent position shifts. The criterion on record:
  decide by ML/analysis target, not fidelity ambition. The targets have now
  arrived.

## New measured evidence (2026-07-31 .. 08-05, this session)

Observables defined identically on real and sim. Real side: `ntp_clus_trk`
of `clusters_seeds_island_79507-0.root_ntuplizer.root` (the canonical real
source), tracks = (event, seedID) groups. Sim side: island91 v5.3
`ntp_cluster` + row-aligned `ntp_truth`, tracks = (event, gtrackID) groups
with `cls==0 && ntrks==1`. Full-crosser gates both sides (>=12 clusters,
lmin<=11, lmax>=50, whole-track circle fit rms <= 0.25 cm), pT windows by
**fitted curvature** (R_fit in [101,137] cm ~ pT 0.5; [357,596] cm ~ 1.5-2.5
GeV) — no truth in any selection.

1. **Per-track circle-fit residual RMS** (Kasa + Gauss-Newton circle fit):

   | window | real | sim v5.3 | data/MC |
   |---|---|---|---|
   | pT ~0.5 GeV (N=694/12,568) | 692 um | 733 um | **0.94** |
   | pT 1.5-2.5 GeV (N=55/899)  | 626 um | 469 um | **1.33** |

   The high-pT excess is ~415 um of real-only scatter in quadrature. It is
   **coherent**: after iterative 3 mm outlier cleaning (ms_realcheck) the
   ratio persists (607 vs 463 um -> 1.31), so it is not association junk.
   It is hidden at 0.5 GeV where crossing-angle smearing dominates both
   sides, and exposed on stiff tracks — exactly where a coherent effect
   must appear. Caveats: real high-pT N=55 (a ~3 sigma statement); the real
   tracker curates its clusters, biasing real widths LOW — the true gap is,
   if anything, larger.

2. **Transverse impact parameter** d0 = | |C| - R | from the whole-track fit
   (distance from the ntuple origin to the fitted circle): real median
   **2.51 cm** (2.69 cm in the audit's cleaned-accepted sample), distribution
   peaking at 2-3 cm, vs sim **0.26 cm** (origin vertices through ~700 um
   clusters). Real primaries do not point back to (0,0) in this frame.

3. Supporting: the previously declared R1 radius taper (item above) belongs
   to the same dropped-corrections family.

## Diagnostic RESULT (2026-08-05, `ms_d0diag`): translation ruled out —
## the effect is distortion-family

The offset-vs-distortion ambiguity is settled. Signed d0 (= |C| - R) vs the
circle-center azimuth phi_c on 3,756 real full-crossers:

- **No frame/beam translation:** cosine fit gives (x0, y0) =
  (+0.03, +0.08) cm — |offset| = 0.8 mm — and explains 0% of the d0
  variance. The GENERATOR branch (beam-spot offset) is CLOSED.
- **The broadening is azimuthally isotropic, per-track:** real d0s RMS
  2.99 cm vs sim 1.61 cm (|d0s| < 8 cm sample); medians 2.51 vs 0.26 cm.
- **Significant drift-time structure:** cosine-subtracted d0s medians per
  tbin bin sit within a few mm in the in-time band but jump to **+1.4 cm at
  tbin ~ 40 and +1.8 cm at tbin ~ 520** (per-bin median error ~0.2 cm);
  the sim control is flat everywhere. The largest biases are on
  OUT-OF-TIME tracks.
- Figure/ledger: `sim_validation_plots/ms_d0diag_v53f.png`,
  `ms_d0diag_v53f.txt` (entry `ms_d0diag` in `ms_real.C`).

Mechanism hypothesis (to be judged by the pipeline thread): the real
reconstruction applies z-dependent r-phi distortion corrections; for
out-of-time tracks those corrections are evaluated at the WRONG apparent z
and mis-apply, producing exactly a drift-time-dependent, azimuthally
unaligned d0 bias — on top of whatever residual static distortion /
alignment imperfection remains for in-time tracks (the high-pT 1.33
circle-RMS excess).

## Requested change (scoped)

1. Transport stage (tier 1): apply the **static distortion map**
   (+ module-edge; time-ordered optional later) as an in-drift distortion
   of electron arrival positions — the sim direction mirrors what the
   container's distortion module did — plus the **TRACKINGALIGNMENT**
   content where it moves TPC surfaces. Payloads are already mirrored in
   `CDB_offline`; the era-1 wiring (`Fun4All_DistortionSim.C`,
   `run_batch_sim.sh` menus) and the 2026-07-20 review's dropped-payload
   map show exactly what was consumed.
2. Optional tier 2 (for out-of-time realism, judged by the drift-time d0
   structure above): additionally emulate the real reco's z-dependent
   correction evaluated at APPARENT z, so out-of-time sim tracks acquire
   the same mis-correction the real ones carry.

## Acceptance meters (ready-made)

- High-pT circle-RMS data/MC: **1.33 -> 1.0** (meter: `ms_real.C`,
  `ms_real_split_v53f.txt`).
- d0: sim develops the isotropic broadening toward real (RMS 1.61 ->
  ~3 cm; medians 0.26 -> ~2.5 cm) and — if tier 2 is taken — the
  drift-time median structure (+1.4 / +1.8 cm bins) reproduces
  (meter: `ms_d0diag`).
- R1 radius taper 1.20 -> 0.40 mm closes (reviewer's per-layer radii,
  2026-07-20).
- **Invariance guard:** pixel-level anchors (pixmean, sub10, windows, step,
  bump, island counts/sizes) must NOT move beyond noise — coherent shifts
  are invisible there by construction; any movement flags an implementation
  error.

## Scope guard

Per the project's own criterion: this matters for consumers of absolute
positions and track-fit quality (track models, d0, any position-consuming
training target). Pure shape/statistics validation does not need it.
Priority vs other pipeline items is the user's call.

## Reproduction pointers

- Macros (self-contained, `island_post/`): `ms_real.C` (entries `ms_real`,
  `ms_realcheck`, `ms_real_split`, `ms_d0diag`), companion context
  `ms_split.C`, `truth_circle.C`.
- Ledgers/figures (v53f tags): `ms_real_v53f.txt`, `ms_realcheck_v53f.txt`,
  `ms_real_split_v53f.txt`, `ms_d0diag_v53f.txt`;
  `sim_validation_plots/ms_real_split_v53f.png` (supervisor-style summary),
  `ms_real_showcase_v53f.png`, `ms_realcheck_v53f.png`,
  `ms_d0diag_v53f.png` (the diagnostic verdict figure).
- Real-data-source rule respected throughout: the only real input is the
  ntuplizer root file.

---

# RESPONSE (sim-pipeline thread, 2026-08-05): measured verdict + partial fix

## Payload verdict: the mirrored distortion payloads CANNOT produce the
## measured residuals — mechanism reattributed

1. **Real ntuple cluster positions are UNCORRECTED and row-locked.**
   Measured: per-layer cluster-r p10-p90 spread 0.06-0.13 mm (discrete
   rows), so no distortion correction is applied to positions in this
   ntuple. The tier-2 "mis-applied correction" mechanism cannot operate on
   these positions.
2. **The static map's bulk is incompatible with run 79507.** The
   with_ModuleEdge map has mean radial deflection +9.7 -> +16.4 mm
   (r=31->71 cm, GROWING with r). A collection shift of that size would
   displace the layer-occupancy structure by ~2 rings (the L21/22 notch is
   aligned real-vs-sim to <1 ring) and contradicts the row-locked radii.
   Applying it in-drift would grossly break the invariance guard and every
   radial-share anchor. Its P (r-phi) component is sub-mm (rms 0.02-0.05
   mm) — far too small for the 415 um circle-RMS excess or cm-scale d0.
3. **Module-edge residual maps are um-scale** (P rms ~3-4 um, R/Z = 0 at
   all sampled radii) — negligible.
4. **The R1 radius taper is GEOMETRY-MODEL content, now measured exactly:**
   real cluster median radius minus our GDML nominal = +1.20 / +0.91 /
   +0.39 mm (L7/12/22), <0.05 mm for L23+ — the real reco's R1 row radii
   differ from our GDML; R2/R3 agree. Full 48-layer table:
   `island_post/real_row_radii_v54.txt` (measured from the run itself).
5. TRACKINGALIGNMENT TPC content (ana537 = ana494 for TPC, identical rows):
   alpha=beta=0; constant dz=-8.2 (global z convention — already absorbed
   in our fitted drift mapping); rigid global rotation (invisible to
   TPC-only meters); mm-scale per-surface dx,dy whose phi-mean radial
   projection is ~0 — does not produce the taper.

## Fix applied now (knob A): real row radii at cluster export

`islandize91.C` gained an optional per-layer radius-offset table
(`real_row_radii_v54.txt`); positions-only, pads/tbins/clustering
untouched -> pixel-level anchors invariant BY CONSTRUCTION (same digi
file). Output: `island91_frames_production_v54.root` (v5.4a). This closes
the taper meter exactly and feeds the circle-RMS meter with the corrected
R1 geometry. Meters running: ms_real_split + ms_d0diag at v54a tags.

## Remaining meters (circle-RMS 1.33, isotropic d0 broadening with
## drift-time structure): require a FITTED model, not a payload

What real clusters carry is phi/z-structured r-phi displacement content
that none of the mirrored payloads describes (bulk wrong, structure too
small). The pipeline-thread proposal: a phenomenological drift-structured
r-phi deflection field (parameterized vs z_true and r), fit to the
measured d0(tbin) profile and the stiff-track circle-RMS excess —
envelope-style methodology (fit family, CRN pilots, holdout gates,
invariance battery). This is a scan campaign (~days), and per this
request's own scope note the priority call is the user's.

---

# RESULTS (sim-pipeline thread, 2026-08-05): primary meters CLOSED (v5.4b)

Fitted two-knob field (probe-derived, one-shot solve, no iteration):
`islandize91` field arg `"0:0|0|0.0426|0|0|2.5"` on top of knob A =
SMOD 0.0426 cm per-(layer,sector,side) hashed rphi scatter (the cluster-
RMS mode) + SPHI 2.5 cm smooth k=2,3 azimuthal-harmonic translation field
(the d0 mode; k=1 excluded by the real null cosine). Official artifact:
`island91_frames_production_v54b.root` (positions-only overlay on the
SEALED v5.3b digi; cluster/label counts identical — invariance exact).

Acceptance (original ms_real entries, v54b tags):
- ms_real_split: stiff-window circle-RMS data/MC **1.33 -> 1.00**
  (real 626 um, sim 629); d0 median 2.51 vs 2.27 cm (was 0.27).
- ms_d0diag: sim RMS **3.01 vs real 2.99 cm**; sim cosine offset 0.09 cm
  vs real 0.08 — isotropy signature reproduced (no fake dipole).
- ms_realcheck (cleaned): stiff **1.00** (627 vs 626).
- R1 radius taper: closed exactly by knob A (+1.20/+0.91/+0.39 mm).

Declared residuals & findings:
1. 0.5-GeV window: 0.94 -> 0.83 (real 692, sim 838). The pT-blind
   geometric scatter EXPOSES a pre-existing sim low-pT over-smear (real
   base05 ~ 564 um vs sim 735 once field is accounted) — response-family
   item (transport step-3), not a field error.
2. The d0s-vs-tbin MEDIAN structure (+1.4 edge / +0.9 CM / +1.8 OOT) is
   UNREACHABLE by translation-class fields (orientation sign-mixing nulls
   medians; measured, P5). Tier-2 as scoped; likely requires an
   orientation-coupled mechanism (e.g. curvature-odd component) — open.
3. Methodology finding for the record: 3-parameter circle refits absorb
   per-track-smooth tangential fields to <1 um; only granular (layer/
   sector) and translation modes survive — payload-shaped smooth maps
   could never have produced these meters even if their bulk were right.

## ADDENDUM — v5.4b verification, truth_circle battery (2026-08-05, analysis session)

Independent rerun of the full truth_circle battery on
`island91_frames_production_v54b.root` (outputs `truth_circle_v54bx.*`):

- **Truth-hit side bit-identical** to v53f (16,626 fits, 20 um, R closure
  0.9990) — confirms v5.4b touched nothing upstream of the export.
- **Reco RMS growth is quadrature-consistent** with the SMOD knob:
  R-selected median 721 -> 816 um (721 (+) 426 = 838; R-window reshuffle
  explains the small gap). Consistent with ms_real numbers.
- **SURFACED A — multi-segment fraction doubled: 1.86% -> 3.69%**
  (cluster-weighted 2.82 -> 5.03%). NOT a label regression (label counts
  are exactly invariant): the per-SIDE hashed field is discontinuous at the
  CM crossing, so side-crossing tracks now split in the phi-chain
  diagnostic. Physically defensible (independent endplates) — but the REAL
  side-to-side discontinuity magnitude is unmeasured. Checkable on real
  data: run the same chaining on side-crossing ntp_clus_trk tracks (zelem
  available) and compare split rates. If real side-crossers do NOT split at
  the sim's rate, the per-side hash amplitude is too large at z ~ 0.
- **SURFACED B — curvature noise: R-window fraction 95.0% -> 87.4%**,
  median segment R_fit 117.7 -> 118.5 cm. The phi-varying SPHI harmonics
  are not a pure per-track translation, so they add few-percent-scale
  curvature (effective pT-resolution) noise. This mode was NOT in the
  acceptance battery. Whether real data carries the SAME curvature noise
  is measurable without truth: compare per-track half-arc curvature
  consistency (R_inner vs R_outer from the split fits) real vs sim —
  proposed as an additional meter before the field model is sealed.

Neither surfacing blocks v5.4b for label-level use (labels exact); both
are field-model shape questions for this thread to judge.

---

# ADDENDUM RESPONSE (pipeline thread, 2026-08-06): both shape questions
# answered by measurement; field revised -> v5.4c

Three-way cmcheck (real / v5.3 fit-floor / field-on):
1. CM discontinuity: REAL side-crossers show a GENUINE membrane jump —
   med|J| 0.480 +- ~0.13 cm (n=48; fit floor 0.304). v5.4b's independent
   per-side phases gave 0.985 (2.5x) and halved the crosser fraction.
   v5.4c: main harmonics side-SHARED + per-side SCM 0.26 cm ->
   med|J| 0.536 (inside real's error band), crossers 1.3% = real.
   Your split-rate diagnostic should recover accordingly — please re-run.
2. Curvature noise: real half-arc dk-width 0.402 (0.5 GeV; floor 0.241)
   vs sim 0.273-0.276 — the harmonics inject LESS curvature noise than
   real carries. No reduction warranted; your R-window drop flagged a
   real mode that is, if anything, undersized in sim.
Primary meters hold on v5.4c (stiff ratio 1.01, d0 RMS 3.05 vs 2.99,
med|d0| 2.27 vs 2.20); labels exact. FINAL ARTIFACT:
island91_frames_production_v54c.root (md5 5972cbe3...); v54b deleted.
cmcheck.C stays in island_post as the shared meter.

---

# RE-ACCEPTANCE (analysis session, 2026-08-06): v5.4c VERIFIED — CLOSED

Independent battery on island91_frames_production_v54c.root (md5 confirmed):
- split-rate diagnostic recovered as predicted: multi-segment 3.69% -> 2.26%
  (v5.3 baseline 1.86% + the real-sized membrane tail); R-window 87.4 -> 89.6%.
- cmcheck reproduced: real med|J| 0.480 (n=48), sim 0.536, crossers 1.3% both;
  half-arc dk05-width sim 0.276 vs real 0.402 (conservative, per Q2 verdict).
- primary meters: stiff circle-RMS 626/635 um (1.01); d0 RMS 2.99/3.05,
  cosine null both sides; truth-hit side bit-identical; labels bit-exact
  (4,524,011 / 3,207,672 / 300,462 unchanged v5.3 -> v5.4c).
v5.4c accepted as the final field artifact from this side. Open items carried:
med05 low-pT over-smear (response family) and the d0-median tbin profile
(tier-2). This document is closed.

---

# v5.5 DELIVERY (pipeline thread, 2026-08-14): field-in-digitization DONE
# (answers digi_field_request.md)

The v5.4c field now acts on transported charge BEFORE pad/tbin binning
(tpc_transport rphifield arg, hash-identical patterns; export runs rowdr-
only so the field enters exactly once). Artifacts (md5s in manifest):
  digi_frames_production_v55.root      (pixels WITH field)
  island91_frames_production_v55.root  (export: rowdr on, field off)
  prodclus_v55.root / hit69_frames_production_v55.root
Your acceptance meters, delivered side: stiff RMS ratio 0.98; d0 RMS
3.07 vs 2.99; CM jump 0.481 vs 0.480+-0.13 (dead-on), crossers 1.3%;
dk05 0.255 (conservative); med05 0.83 as you predicted (response-owned).
No SMOD re-solve was needed (quantization absorption ~ -2%, in-noise).
Labels: not bit-exact vs v53 by design; class balance drift <=0.4%
(island91 56.6/39.7/3.7; prodclus 53.7/44.6/1.6). Content note: px/ev
+2.9% vs target, quoted with the NEWLY MEASURED composer realization
sigma ~+-2-3% per production (era-wide correction — v5.3b's +0.2% was
partly a draw; ledger has the full investigation + the pre-declared
re-roll bookkeeping). Your nf_digipix differential predictions are ready
to test on digi_v55 — over to you.
