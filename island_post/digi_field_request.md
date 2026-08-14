# Request: field-in-digitization production (pixel-level field realism)

2026-08-13, from the ML/analysis session to the pipeline session, relayed by the user.
Follow-up to `distortion_remodel_request.md` (closed). That loop put the measured field
into CLUSTER positions (v5.4c, cluster-export overlay on the sealed v5.3 digi). This
request asks for the same field one stage earlier, so the PIXEL level carries it too.
Sections marked (discussion) wrap up the design conclusions agreed with the user on
2026-08-13 so the full rationale travels with the request.

## What the "field" is — agreed definition (discussion)

"v5.4 introduced an empirical position-residual field: the measured discrepancy
between reconstructed positions under the ideal GDML geometry and those observed in
run 79507 — anchored to the measured row radii, the stiff-track circle-fit width, the
d0 distribution, and the central-membrane jump. Its physical origin (mechanical
misalignment vs drift distortion) is left open by design; only the CDB static map is
excluded (real clusters uncorrected + row-locked; the map's cm-scale bulk matches
nothing observed)."

Anchor-by-component:
- rowdr (knob A): per-layer real reco row radius minus GDML nominal, measured
  directly (`real_row_radii_v54.txt`) — the one literally geometry-vs-geometry part.
- SMOD 0.0426 cm per-(layer,sector,side) hashed r-phi scatter: sized to the
  stiff-track circle-RMS meter (data/MC 1.33 -> 1.00).
- SPHI k=2,3 side-shared harmonics (amp 2.49): sized to the d0 RMS/median and its
  cosine isotropy; k=1 excluded BY a real-data null.
- SCM 0.26 cm per-side membrane term: sized to the measured CM jump 0.480 +- 0.13 cm
  (the only boundary-anchored component).

## Why (measured motivation)

`ms_nofinder.C::nf_digipix` (ledger `ms_nofinder_v54c.txt`, figure
`ms_nofinder_digipix_v54c.png`) ran the first matched pixel-level comparison:
real ntp_hit pixels (tracker-grouped) vs sim DIGI pixels (per-pixel gtrackID,
`digi_frames_production_v53.root` = the v5.4c input; 60 frames, 24,521 tracks):

- GLOBAL whole-track fit: real 2358 um vs sim 1497 um -> data/MC 1.57.
  In quadrature, real carries ~2.0 mm of long-range content vs sim ~0.67 mm.
  Part of the ~1.9 mm real-only excess is the distortion field, which sim pixels
  cannot carry by construction (the rest is road-association tails - declared).
- LOCAL 4-row short-sagitta fit: real 1205 vs sim 1337 um -> data/MC 0.90.
  Field-blind response meter: raw charge clouds agree to ~10%.

Wanted: field on BOTH sides of the global panel, and field-realistic sim pixels for
any future pixel-level training.

## Pipeline-stage map — what changes where (discussion)

NOW (v5.4c):
  charge from G4 -> [digitization: drift + diffusion + pad/tbin binning] -> pixels(v53)
                 -> clustering -> [export nudge: rowdr + FIELD on cluster xy] -> clusters(v54c)

REQUESTED:
  charge from G4 -> [digitization + FIELD displaces charge BEFORE binning] -> pixels(new)
                 -> clustering -> [export nudge: rowdr only, FIELD OFF] -> clusters(new)

Per stage: digitization gets the one real code change (displace each transported
electron's arrival point by the field before pad/tbin binning, so the warp expresses
itself as changed charge sharing between neighboring pads — as in reality); pixels are
regenerated, not edited; clustering code is untouched, just re-run; the export keeps
rowdr (a lookup convention about what radius a row is called — nothing physical to
inject) and switches the field components OFF so the field enters exactly once.
Untouched entirely: generator, G4 collision libraries, frame composer, clustering
algorithm.

## The ask

1. Apply the v5.4c field model as an r-phi displacement of the transported charge
   BEFORE pad/tbin quantization, in the digitization stage. Same model, same measured
   parameters as the v5.4c cluster overlay: islandize91 field arg
   `"0:0|0|0.0426|0|0|2.49|0.26"` semantics.
2. Re-digitize the 250-frame pp production from the existing G4 libraries, then
   re-export island91/prodclus/hit69 with the islandize91 field arg OFF.
3. NOT requested: rowdr (knob A) — stays at export/reader level as a lookup
   convention. (Noted for my own readers: pixel x,y should use the measured real row
   radii where mm-level correctness matters.)

## Why this stage and not elsewhere (discussion)

Why v5.4a/b/c never had a digi of their own: the overlay lives DOWNSTREAM of pixels
(it nudges cluster positions at export), so the pixel stage of v5.4c IS the sealed
v5.3 digi, byte-identical — a "digi_v54c" made by overlay is impossible in principle,
because pixels are integer (phibin, tbin) grid objects with no position between pads.
Real distortion acts on the drifting charge before it lands; faithful modeling must
therefore act inside digitization.

Why not bake the field into the GDML ("decode, add, recode — once and for all"):
1. No object to warp: pads/rows do not exist in the GDML — the sensitive volume is
   one featureless gas cylinder and the pad segmentation lives in readout software.
   Moving imaginary boundaries inside uniform gas changes no deposit.
2. Wrong object: neither candidate mechanism touches the particle — misaligned pads
   report shifted coordinates, deflected drift electrons land shifted; both live
   between trajectory and pixels. A GDML warp would distort TRUTH itself, corrupting
   every truth-anchored result (circle test, MS measurement, ML label semantics).
3. Cost: the GDML feeds the sealed G4 libraries and every content/shape tune since
   v3.6; changing it reopens all of that. Digitization injection reuses the G4
   libraries unchanged.
4. "Once and for all" is backwards: the field is run-79507-specific and provisional
   (space charge varies with conditions; a corrected official reprocessing would
   re-size it). It belongs in a versioned, swappable stage. And it is DEFINED as
   measured-real minus ideal: bake it into the ideal and the next measurement has no
   reference left.

## Expected acceptance consequences (plan for these, they are not failures)

- Labels will NOT be bit-exact vs v5.3/v5.4c: charge crossing pad boundaries changes
  cluster membership/fragmentation. Quantify label agreement + class-balance drift.
- Cluster-level meters re-verified against the v5.4c matched set (stiff circle-RMS
  ratio, d0 RMS/median + cosine null, cmcheck CM jump + crossers, split-rate). The
  effective cluster-level smear of an injected field may land slightly below the
  overlay's (quantization + centroiding absorb part of the displacement) -> a small
  SMOD re-solve against the stiff meter is plausible.
- The low-pT cluster RMS ratio (med05 family) will REMAIN ~0.8 — the response
  (charge-cloud width) is deliberately untouched by this request. Do not read that
  as the injection failing: v5.3's apparent 0.94 there was two errors canceling
  (missing field vs too-wide clouds; sqrt(733^2+426^2)=848 ~ the 842 now measured).
  The response fix is a separate, still-open axis (sim base ~735 vs real ~564 um at
  0.5 GeV; pixel-level local meter says clouds ~10% wide).
- ML datasets derived from cluster labels regenerate on the new production.

## Differential predictions (my acceptance meters, ready to run)

- nf_digipix GLOBAL: the real-only long-range excess shrinks by the field share;
  sim global RMS rises from 1497 um toward real; data/MC 1.57 moves toward 1.
- nf_digipix LOCAL: 0.90 stays put (local fits reject the smooth field; only the
  granular SMOD share, ~400 um in quadrature on mm clouds ~ +4%, may appear).
- Cluster-level battery: reproduces v5.4c values within errors after any re-solve.
A pass on all three = the field landed in the right stage, exactly once.

## Version scorecard for context (discussion; why the v5.4 direction stands)

| meter | v5.3 | v5.4c |
|---|---|---|
| stiff cluster-RMS data/MC | 1.31-1.33 | 0.99 |
| d0 (2.5 cm off-origin) | absent | matched (2.37 vs 2.51) |
| CM jump | absent | real-sized (0.536 vs 0.480) |
| half-arc curvature width | sim < real | still conservative (0.276 vs 0.402) |
| low-pT cluster RMS | 0.94 (cancellation) | 0.82 (honest, response-owned) |
| pixel-level local response | 0.90 | 0.90 (same pixels) |

The matched local fits graded the shared pixel stage, not the overlay; the response
flaw they surfaced predates v5.4 and is out of this request's scope.

## Cost/benefit framing for the decision

Local re-digitization of 250 frames + exports + acceptance battery (order hours of
compute, one verification campaign). Benefit is entirely at pixel level: matched
global comparisons and pixel-level ML realism. v5.4c remains valid and preferred for
cluster-level work regardless of this decision; nothing here reopens the closed
distortion loop. Version tag (v5.4d vs v5.5) is the pipeline session's call.
