# Side-Review Verdicts — sim-pipeline optimization project

Standing record of the bystanding-reviewer thread. Each entry states the claims
examined, the independent verification performed (all real-data checks use
`clusters_seeds_island_79507-0.root_ntuplizer.root` exclusively, per the
project's source rule), and the ruling. Newest first.

---

## 2026-08-05 — Adjudication: the circle-fit 1.33 discrepancy
**Parties.** The ML/analysis thread ("circle-fitting agent", request in
`island_post/distortion_remodel_request.md`) measured a ×1.33 real/sim excess
in stiff-track circle-fit residual RMS plus a 2.5 cm real d0 broadening, and
requested the mirrored distortion+alignment payloads be applied at transport.
The sim-pipeline thread ("optimizer") rejected the payload mechanism,
re-attributed the effects, and produced v5.4a (real row radii) + v5.4b
(fitted two-knob r-phi field).

### Independent verification (this review, 2026-08-05)

| claim | claimant | my measurement | verdict |
|---|---|---|---|
| Real `ntp_cluster` positions row-locked (uncorrected) | optimizer | per-layer r p10–p90 = 0.052/0.070/0.070/0.123 mm (L7/12/22/40) | **CONFIRMED** |
| Same for `ntp_clus_trk` — the tree the circle meters actually use (reviewer's own cross-examination) | — | identical medians and spreads to ntp_cluster | **CONFIRMED** — the tier-2 "mis-applied correction" mechanism cannot act on the meters' inputs |
| Static map bulk ⟨dR⟩ +9.7→+16.4 mm growing with r | optimizer | +8.4–9.8 mm (r=31) → +14.0–17.1 mm (r=71), per side | **CONFIRMED** — incompatible with row-locked radii + aligned L21/22 notch |
| Static map r-phi content sub-mm | optimizer | ⟨dP⟩ ≤ 0.12 mm, rms 0.017–0.092 mm | **CONFIRMED** — ~4× too small for the 415 µm excess, ~25× too small for cm-scale d0 |
| Real row radii table (knob A source) | optimizer | matches this reviewer's independent 2026-07-20 radii audit to the digit (L7 +1.20, L12 +0.91, L22 +0.40 mm) | **CONFIRMED** |
| v5.4b invariance guard | optimizer | v53 vs v54b: 8,032,145 clusters identical; label counts 4,524,011 / 3,207,672 / 300,462 identical | **CONFIRMED EXACT** (positions-only overlay) |
| Knob A closes the taper | optimizer | L7 mean r 31.3805 → 31.5002 vs real 31.5001; L30 shift −0.013 mm | **CONFIRMED** |
| v5.4b closes the stiff meter 1.33 → 1.00 | optimizer | re-ran `ms_real_split` on v54b: real 626 µm vs sim 629 → **1.00**; d0 median 2.51 vs 2.27 cm | **CONFIRMED** (digit-identical to the ledger) |
| Declared low-pT trade | optimizer | reproduced: 692 vs 838 µm → 0.83 (sim now 21% wide) | **CONFIRMED as declared** |
| The 1.33 excess itself (v5.3 baseline) | circle agent | baseline ledger + optimizer's independent `fieldmeter` reproduction (digit-exact per ledger); effect coherent (survives outlier cleaning, 1.31) | **ACCEPTED** |

### Ruling

1. **On the measurements: the circle-fitting agent is correct.** The 1.33
   stiff-track excess and the d0 broadening are real, coherent, and were
   measured with sound instruments (full-crosser gates, curvature-based pT
   windows, no truth in selection). Its invariance guard and scope guard were
   also correct and were honored.
2. **On the mechanism: the optimizer is correct.** The requested route is
   impossible on four independently verified counts: real positions are
   uncorrected and row-locked in *both* trees (tier-2 moot); the static map's
   radial bulk contradicts the run's own radii/occupancy; its r-phi content is
   an order of magnitude too small; the alignment payload's TPC content cannot
   produce the taper. The taper is geometry-model content (GDML vs real row
   radii), now closed exactly.
3. **On "not amendable": that premise mischaracterizes the record.** The
   optimizer rejected the *requested payload route*, not the fixability — it
   then amended the sim (knob A + fitted SMOD 0.0426 cm / SPHI 2.5 cm k=2,3
   field) and closed both primary meters (1.00; d0 RMS 3.008 vs 2.987 with
   the isotropy signature reproduced), which this review reproduced
   independently. What remains declared un-amendable *in the current model
   class* is only the d0-vs-tbin **median** structure (tier-2 residual,
   needs an orientation-coupled mechanism).

### Standing reviewer caveats

1. The v5.4b field is **descriptive, not derived**: two fitted amplitudes
   closing two primary meters, constrained by a genuine probe ladder (P1–P6
   mechanism decomposition) and passing untargeted checks (med|d0| 2.195 vs
   2.203; cosine isotropy 0.09 vs 0.08 cm). The physical identity of a
   2.5 cm k=2,3 azimuthal-harmonic displacement is unassigned — a question
   for the team (space charge? field cage? reco frame?), not for tuning.
2. The exposed **low-pT over-smear** (0.83; real base ~564 µm vs sim 735) is
   now the largest open track-level residual — response/transport family,
   pre-existing, previously masked by the missing field.
3. The real d0(tbin) median profile (+1.4 / +0.9 / +1.8 cm) remains open.
4. Statistics: the stiff real window is N=54 (~3σ statement); d0 median is
   90% closed (2.27 vs 2.51 cm), not exact.
5. The real tracker curates its clusters (real widths biased low), so the
   true pre-fix gap was, if anything, larger — unquantified.

Artifacts of this verification: `island_post/ms_real_split_v54brev.txt`,
`sim_validation_plots/ms_real_split_v54brev.png` (+showcase), review-tagged,
originals untouched.

---

## Appendix — prior verdicts of this review thread (summary)

**2026-07-20/21 — v4.0 supervisor-handoff audit.** All five manifest md5s,
hit69 branch-order identity (69/69), row closures (67,315,606 = digi−808;
event 43 = 271,654/25,523/15,141 = CSVs), z/phielem/seed conventions, R1
radii declaration, and all headline figures verified; bestmatch scan
reproduced exactly (f223 L12 score 0.14; watch 60/61). Findings: `hits69.C`
lacked the event-sort (fixed 07-22, verified on disk); sim ADC
float-continuous vs real integer ADU (traced to post-quantization ion-tail
addition at `tpc_digitize.C:640`; 88.7% of pixels fractional); CSV handoff
caveats (silicon rows, event-74-vs-median-43 occupancy); stale
funny_shapes_v40 title. Reviewer's own stage-map fieldmap misattribution
corrected by the optimizer (gap-rebuild map, 26,144,996 grid points —
verified in the G4 logs).

**2026-07-22 — v5.1 claims verification.** Species = pp: scaler evidence
verified (rmbd range 196.5–362.3 kHz exact; rawzdc/rawmbd = 0.0161; N/S
in-time symmetry 1.000); ledger imprecision noted (rzdc "11–28 kHz" is the
early/nonzero subset; full sample [0, 27.7] mean 7.0 kHz, declining to zero
mid-run — verdict unaffected). Truncation: 62 complete / 38 truncated
verified to the digit (endpoints [957,970] vs 94–838; complete-subset mean
323,799, CV 0.239). Mechanism measured: stops are detector-global,
activity-independent (corr −0.08), ~uniform — DAQ throttle-gate class
(~135 µs equivalent ON), team confirmation pending.

**2026-07-24 — ε_MBD analysis and storage audit.** ε_MBD classified across
the residual set (rate-scaled vs rate-shaped vs rate-free); internal meters
triangulate inconsistently (windows ε≲0.50, step ε≈0.6, CV ε≳0.52) —
supports the generator-multiplicity attribution; scan protocol specified
(fit on time-structure family only, rate-free anchors as holdouts).
Storage: none of coresoftware/macros-gcc12/modern_macros/macros-offline in
the v5.3 runtime; ~185 GB deletion menu delivered (optimizer executed a
~136 GB sweep on 07-30 per ledger).
