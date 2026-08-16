# Side-Review Verdicts — sim-pipeline optimization project

Standing record of the bystanding-reviewer thread. Each entry states the claims
examined, the independent verification performed (all real-data checks use
`clusters_seeds_island_79507-0.root_ntuplizer.root` exclusively, per the
project's source rule), and the ruling. Newest first.

---

## 2026-08-16 — Finding: small loopers exist in the sim with real morphology, but at ~1/3 the real rate

**Question (user):** can the sim reproduce the "small loopers" visible in the
supervisor's real event-74 3D cluster view (`Visualizing/images/event74_ntp_cluster_py.pdf`,
from `Visualizing/src/plot_event74_hits.ipynb`)?

**Method.** New macro `island_post/looper_view.C`: (a) a ROOT twin of the notebook
view — same cuts (`zelem==1`, 7<layer<55), same (z, x, y) axes and orientation
(matplotlib elev 24/azim 62 ⇒ ROOT theta 24/phi −118, verified with a synthetic
helix), adc-tertile colouring; (b) a **finder-free small-looper census applied
identically to real and sim** — 3D proximity linking (3 cm) into chains, Kasa xy
circle fit, accept R ≤ 12 cm, xy rms ≤ 0.6 cm, ≥ 30 clusters, ≥ 2.5 unwrapped
turns, z-span ≤ 130 cm. No truth on either side (truth was used only once, to
calibrate the criteria on known sim loopers: 30–250 clusters, R ≈ 4–12 cm, z-span
40–100 cm). Ledger: `island_post/looper_census_v55.txt`; figures
`sim_validation_plots/looper_view_ev{236,81}_v55.png` (event 81 = the sim's
richest, 6 loopers incl. a face-on ring).

**Result (99 real vs 250 sim events):**

| | real | sim v5.5 |
|---|---|---|
| small loopers / event | 2.27 (81 % of events ≥ 1, max 7) | 0.85 (55 % ≥ 1, max 6) |
| clusters / event (zelem==1) | 13,295 | 15,780 |
| **loopers per 10k clusters** | **1.71** | **0.54** |
| median R, turns of found loopers | 1.8 cm, 6.7 | 1.8 cm, 7.1 |

**Verdict.** Yes — the sim contains the same object, and its *morphology* is
matched to the digit (R and turn count). But the *rate* is ~3× low per cluster.
Caveats before treating this as a residual: (1) the census is greedy and cut-based,
so absolute rates depend on the thresholds — the ratio is the robust number, and
it survives loosening (at R ≤ 25/≥ 1.5 turns the same direction held);
(2) the truth-based hunt found 3,258 loopers (R ≤ 25, ≥ 1.5 turns) in the sim, so
the finder is not what's starving the sim — the deficit is in the *tight, many-turn*
class specifically. Candidate owners, in order: soft-electron content (δ-rays and
photon conversions in material — the CEMC effective-fill and passive-material
approximations feed exactly this class), generator soft spectrum, and the response
axis (wide sim clouds merge/blur tight coils so fewer pass the rms gate). This is
consistent with the standing "content owner = response/conditions side" attribution
and adds a specific, cheap meter for it. Recommended as an acceptance meter for the
next production; not a blocker for v5.5 adoption.

---

## 2026-08-14 — Review: v5.5 (field-in-digitization) + the closed circle-fit loop

**Scope.** The v5.5 campaign (`digi_field_request.md` follow-up: the v5.4c
empirical field injected at transport, before pad/tbin binning), reviewed
together with what became of the 2026-08-05 adjudication (v5.4c revision,
external-data corrections of 2026-08-08, and the two v5.5 honesty
corrections).

### Independent verification (this review, 2026-08-14)

| claim | my check | verdict |
|---|---|---|
| v5.5 artifact md5s (digi/island91/prodclus/hit69) | all four recomputed | **CONFIRMED** — match manifest |
| Field hash bit-identity, transport ↔ islandize91 | seed constants, sector convention `phw/(π/6)`, harmonic arguments compared in code | **CONFIRMED** identical |
| Injection point per request (post-diffusion, pre-gap/zigzag, Δφ=d/r, no rng draws consumed) | code read; off-path guard = chunk-0 no-field regen byte-exact (136,031,407 px + sumq) | **CONFIRMED** |
| "Field enters exactly once" | export call passes `real_row_radii_v54.txt` + field `""`; L7 cluster mean r = 31.5002 (rowdr on) | **CONFIRMED** |
| Stiff cluster meter on v5.5 | re-ran `ms_real_split`: real 626 µm vs sim 614 → data/MC 1.02 (manifest's 0.98 = same numbers as sim/real); d0 median 2.39 vs 2.51 | **CONFIRMED** (~1.0 within realization noise) |
| Requester's differential predictions | LOCAL 0.90 → 0.90 exact; GLOBAL 1.57 → 1.53 (in-quadrature sim gain √(1537²−1497²) ≈ 348 µm ≈ the granular SMOD share — the smooth SPHI is absorbed by circle fits, per the campaign's own P1/P5 finding, so a small move is the *expected* signature); cluster battery reproduces v5.4c, CM jump 0.481 vs real 0.480±0.13 (better than v5.4c's 0.536) | **PASSED** as pre-registered |
| Class balance | measured 56.7/39.9/3.4 on the shipped draw-2 file | nit: manifest line quotes the draw-1 values (56.6/39.7/3.7) — update wording, not content |

### Ruling

**v5.5 is sound and adoption-ready.** The implementation honors the request
precisely (same field, one stage earlier, entering exactly once), the guards
are the strong kind (byte-identity of the off path; sealed v5.3 libraries
untouched as `raw_lib_pp55_*` is a new family), and the delivery battery
reproduces the v5.4c cluster-level closures without an SMOD re-solve.
Convention going forward, as the manifest states: pixel-level consumers use
v5.5; v5.4c remains valid for cluster-level work on the sealed v5.3 pixels.

### Commendations and cautions for the record

1. **Two honesty corrections raise the era's epistemic standard.** The
   measured composer realization σ (~±2–3 % per production) retroactively
   attaches to *every* content closure since v5.0 — including numbers this
   review previously verified as reported (the values stand; their implied
   precision was never real). The pre-declared re-roll rule was amended
   post-hoc (draw 2 at +2.86 % vs the 2.5 % band) — transparent and, in my
   judgment, correct in outcome ("accept the closer of two pre-registered
   draws, no further rolls" preserves the anti-shopping intent), but future
   pre-registrations should state the tie-break rule upfront.
2. **The content deficit's owner has been narrowed by exclusion.** The
   2026-08-08 world-data confrontation (dNch/dη vs PHOBOS) plus the 16-point
   no-go scan refuted the generator-multiplicity attribution; with ε_MBD now
   externally closed (σ_MBD 25.6 mb ⇒ ε≈0.61, consistent with the scanned
   0.588±0.033 — the scan-then-compare structure recommended by this review
   on 07-24 played out as designed), the open owner is the **response /
   conditions axis** — which converges with med05 0.83 (sim clouds wide at
   low pT) and the pixel-level LOCAL 0.90. Three independent arrows now point
   at one target: the charge-cloud response. That is the pipeline's principal
   open front.
3. **Do not chase the pixel-global 1.53 with field amplitude.** The
   remaining real-side excess in that panel is association-tail-dominated
   (tracker-road grouping) and the smooth field is absorbed by circle fits on
   both sides; the meter's role is monitoring, not closure.
4. Open items, ranked: response axis (above); d0-vs-tbin medians (tier-2,
   unchanged); σ_z ≈ 50 cm beam spread vs the sharp step (reopened by
   published data — connects to this review's earlier note that all sim
   primaries originate at (0,0,0,0)); circle-thread re-acceptance
   (split-rate) of v5.4c/v5.5; manifest class-balance line to draw-2 values.

Review artifacts: `ms_real_split_v55rev.{txt,png}` (+showcase), review-tagged,
originals untouched. Workflow diagram renewed for v5.5 (user request, 2026-08-14):
`sim_validation_plots/pipeline_v55_workflow.{tex,pdf,svg}` — TikZ single source,
every number traceable to `pp_pipeline.sh` (VER=v55), the v5.5 ledger entries, the
manifest, or a side-review re-run; `pipeline_detachment.svg` is kept untouched as
the detachment-era record. Nit found while sourcing it: `pp_pipeline.sh`'s comment
on `RSPEC` says $\varepsilon$ 0.519, but the live `rspec99_v52.txt` header says
0.588 (the scan band centre, which is what production used) — stale comment,
worth a one-line fix.

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
