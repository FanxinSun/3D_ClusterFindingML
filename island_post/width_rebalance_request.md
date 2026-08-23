# Request: width-axis re-balance after the twist injection — SMOD/SPHI/SCM down-retune

2026-08-23, from the side-review thread to the pipeline session, relayed by the user.
Follow-up to `twist_field_request.md` (delivered 2026-08-23, reviewed same day:
SIDE_REVIEW_VERDICTS.md top entry — injection verified end-to-end, paired-library
check 19 um mean deviation over 96 cells). This request executes the down-retune
that the twist request itself pre-declared ("no retune ... re-measure C(d) after
this injection before touching it") and that the delivery parked as the WIDTH FLAG.

## User decision on record (2026-08-23)

Option chosen: **re-balance on top of the delivered V6-twist state** — NOT a rewind.
Rationale (review): the twist amplitude is measured, not fitted, so a "redo" could
only converge to the same end state at the cost of discarding a verified production
and the clean pp55/pp6t A/B pair; the intermediate state honestly separates measured
physics from fitted compensation (v5.3->v5.4 low-pT precedent).

## Why (measured motivation)

The v5.4b/c smooth terms were fitted to close the TOTAL real width — which, we now
know, already contained the twist. With the twist explicit, sim double-counts that
share. Delivered numbers (V6-twist vs pre-twist, real reference in parens):

| meter | pre-twist | post-twist | real |
|---|---|---|---|
| stiff circle RMS [um] | 618 (ratio 0.99) | 728 (ratio **1.16**) | ~626 |
| med05 [um] | 829 | 907 | 692 (response-owned axis) |
| CM jump [cm] | 0.484 | 0.529 | 0.480 +- 0.13 |
| half-arc dk (0.5 GeV) | 0.26 | 0.34 | 0.402 (still conservative) |

The excess is quadrature double-counting, not an implementation error: the twist's
stiff-meter share is sqrt(728^2 - 618^2) ~ 385 um, and the smooth terms still carry
the full pre-twist budget.

## Sequencing gate

The amplitude solve waits for the probe thread's **post-injection C(d)** (their
pending item) — the per-region steps may account for part of the granular-vs-smooth
correlation-length mismatch, which decides whether SMOD only shrinks or its
granularity model also changes. Prep (naming, seed rule, pilot harness) can proceed
in parallel. If C(d) is slow, the user may waive the gate explicitly — do not waive
it yourselves. The probe thread's formal acceptance of the twist delivery should be
collected in the same exchange.

## The ask

1. **Down-retune the smooth field terms on top of the frozen twist.** Free
   parameters: SMOD amplitude (primary); SCM trim (optional, round 1 only if the
   CM meter demands); SPHI frozen unless the d0-family holdouts move. The twist
   table `twist_payload_v6.txt` is **byte-frozen — never rescaled**: it is a
   measurement, not a knob.
2. Solve by the established pilot method (chunk-0 CRN pilots, Jacobian, one-shot
   solve — v5.4b P-ladder precedent), not by arithmetic. For orientation only:
   naive quadrature transfer ignoring the response floor gives SMOD' ~ 0.034;
   including the fieldcmp ideal floor (454 um) gives ~ 0.021 — the truth is in
   that bracket and the pilot decides. `ms_fieldcmp.C` (ideal-vs-field share
   decomposition) is the natural instrument.
3. One re-balance production + full delivery: exports (island91/prodclus/hits69),
   acceptance battery, figures + residual table + briefing regenerated on the new
   artifacts (standing user rule: a version is not delivered until they are),
   manifest md5s, PIPELINE.md entry.

## Acceptance (pre-registered)

Close:
- stiff circle RMS ratio back to 1.00 within the realization band;
- CM jump within 0.480 +- 0.13;
- med05 back to its pre-twist family (~829/0.83 — its remaining gap to real 692 is
  the response axis, explicitly OUT of scope here);
- d0 RMS / med|d0| hold their v5.4c-family values.

Frozen holdouts — must NOT move:
- split-half Dsagitta stays +6.3 mm within +-0.3 per side (the signed anchor: the
  guard that shrinking SMOD does not quietly eat the twist);
- twist profile RMS-diff to real stays <= ~150 um (`twist_probe.C`, unchanged);
- pixel pooled residual core q68 = 1.88 mm, bit-stable;
- spectral/content family (windows/bump/step/pixmean/sub10/shapes) within the
  declared realization sigma (+-2-3% content; no seed-shopping — pre-declare the
  seed rule before looking, one re-roll max, closer-of-two).

## Scope guards (unchanged from the twist request, restated)

NOT in scope: response/ZS/cloud knobs, the road, the generator and G4 libraries,
the geometry table, the twist payload, rowdr. If any acceptance meter appears to
require touching these, stop and report instead of tuning.

## Naming and provenance (hard requirement this time)

Mint **new version names before launching** (v6.1 / v6t1 — the tag is the pipeline
session's call) and a new lib family: no second in-place re-cut of `_v6` files.
The v3.7 lesson and the review's flag on the twist delivery both apply. Keep
pp55 and pp6t families on disk (the A/B pairs are working instruments now).

## Cost

Chunk-0 pilots + one production + delivery battery: ~2.5 h compute, one campaign.

## Decisions left to the pipeline session

- Version tag; whether SCM is trimmed in round 1 or left for the meter to decide;
  whether the pre-declared re-roll is spent if windows land at the band edge.
- Post-delivery: lift the interim flag "V6t clusters not for ML training" once the
  stiff meter is back at 1.00 (the ML thread should be told in the same handover).

---------------------------------------------------------------------------------
## PROBE-THREAD GATE RESPONSE (2026-08-23, Real Data Probe session)

The three gated items, measured independently on the delivered digi
(md5 ba837766), tags v6t; full record REAL_DATA_PROBES.md 2026-08-23.

1. TWIST FORMALLY ACCEPTED. Dsagitta sim +6.31 (side0 +6.24 / side1 +6.37) vs
   real +6.53 (+6.04/+7.15) — per-side deltas 0.20/0.78 mm, inside every gate;
   profile RMS-diff 127 um (was 447); boundary jumps residual 84-165 um (were
   ~1700-1850); charge independence preserved. Holdouts: pixel pooled q68
   1.88 mm bit-stable; clipped GLOBAL 1.04 -> 1.01 (toward real — the twist
   supplies previously-missing real LF: unexplained global surplus 418 -> 169
   um; this is the expected sign, not a disturbance).
2. C(d) RE-MEASURED: sim 0.57/0.50/0.39/0.22 at d=1/2/4/8 (was 0.77/0.67/
   0.55/0.25; real 0.30/0.18/0.08/-0.06). ~40% of the short-lag gap closed;
   sim still smooth-dominated (zero-crossing ~12-16 rows vs real ~6; sim
   C(24)=-0.32 = undiluted sawtooth, more evidence of the missing granular
   share).
3. RULING: **granularity/composition-change, NOT amplitude-only.** Post-twist,
   sim's granular share is BELOW real (cell coherence 145 vs 307 um; real
   surplus ~270 um; remaining LF surplus ~426 um correlated over ~2-3 rows)
   while sim's width excess (stiff +385 um in quadrature) sits in the SMOOTH
   share. An SMOD-amplitude-only shrink passes the pre-registered width gates
   but moves BOTH composition meters away from real (C(d)/C(0) up, cell
   coherence down). Options consistent with the data: (i) land the down-retune
   on SPHI (tension with the SPHI freeze / d0 anchors — pilot to reconcile), or
   (ii) re-model SMOD with ~2-3-row correlation (per-(side,sector,row-block)),
   letting it carry more width share at real's correlation shape. If v6.1
   proceeds amplitude-only, the composition axis stays open — record it; this
   thread will re-measure C(d)/cell coherence on v6.1 either way.
   Acceptance battery for v6.1 on this side: twist_probe.C + nonrms_probe.C +
   gtail_probe.C, unchanged, tag v61.

---------------------------------------------------------------------------------
## ADJUDICATED AMENDMENT (side-review thread, 2026-08-23) — scope updated per the
## gate response above; supersedes "SMOD amplitude (primary)" in The Ask

1. **Ruling accepted.** Amplitude-only is rejected as the primary path: it would
   pass the pre-registered width gates while moving BOTH composition meters away
   from real — a gate-passing cancellation of the kind this project has refused
   before (v5.3 low-pT 0.94 precedent). The composition axis is hereby promoted
   from probe-side diagnostic to acceptance gate.
2. **Chosen path: option (ii) — re-model SMOD with row-block correlation.**
   Grounds: SMOD as coded hashes independently per (layer, sector, side), so it is
   row-WHITE along a track — all C(0), nothing at d>=1 — while real's granular
   content is per-cell coherent (307 vs sim 145 um) and correlated over ~2-3 rows.
   A block-correlated SMOD manufactures exactly real's correlation shape and can
   therefore carry MORE width share while the total closes. Option (i) (land the
   cut on SPHI) is rejected as primary: SPHI is translation-mode, largely absorbed
   by circle fits (the v5.4b P1/P5 finding), so its stiff-RMS share is far too
   small to yield 385 um without gutting the d0 anchors — the tension the gate
   response itself flags — and it does nothing for the granular deficit. The
   pilot Jacobian may still measure an SPHI row as a cross-check.
3. **Implementation spec** (pipeline session refines): extend the field syntax
   with a per-(side, sector, row-block) hashed term (FSBLK; block length nb,
   pilot scans nb in {2, 3}), optionally retaining a reduced row-white SMOD
   residual. New hash-seed family distinct from all existing ones; deterministic,
   zero rng draws; **off-path byte-identity guard mandatory** (term-off regen of
   chunk 0 must match the sealed pp6t library exactly — the v5.5/twist
   discipline). In tpc_transport. Note for the record: islandize91's overlay
   field code will no longer be semantics-complete for the new component —
   mirror it or explicitly mark the overlay path deprecated (pipeline's call,
   ledger it either way).
4. **Solve**: the pilot ladder now solves the block amplitude (+ optional white
   residual) jointly against the width gates AND the composition meters. SCM
   trim optional as before; SPHI frozen; twist payload byte-frozen.
5. **Acceptance additions** (direction + minimum closure set here; the probe
   thread, as meter owner, fixes exact bands in one line BEFORE the solve pilot):
   C(1) from 0.57 toward real 0.30 by at least half the gap (<= ~0.44); cell
   coherence from 145 toward real 307 um by at least half (>= ~225 um). Measured
   by the probe battery (twist_probe + nonrms_probe + gtail_probe, tag v61).
   All original close-gates and frozen holdouts stand unchanged, including
   Dsagitta +6.3 +-0.3 per side and the pre-declared seed rule.
6. **Fallback clause**: if the pilot shows the block model cannot reach the width
   gates without violating a holdout, STOP AND REPORT — no silent scope changes.
   The amplitude-only fallback remains available only with the composition
   regression explicitly recorded and user sign-off obtained first.
7. **Cost delta**: ~20-line field-code extension + off-path guard + nb scan on
   the existing chunk-0 harness — still one campaign, ~3 h.

---------------------------------------------------------------------------------
## CLAUSE-6 RESOLUTION (side-review adjudication #2, 2026-08-23) — after the
## pipeline pilot's stop-and-report

1. **Falsification accepted; the amendment's composition gates are RETRACTED as
   v6.1 gates.** The pilot refuted the amendment's direction claim: replacing
   white (C=0) with block-correlated content (C=0.5-0.67) raises C(1) — the
   reachable band at fixed width budget is [0.63 all-white, 0.67 all-block]
   pilot-equivalent, and the gate (<=~0.49) is unreachable by construction. The
   error was shared by the probe ruling and this review's adjudication; the
   pilot's mixture cross-check (<0.01), dilution control (occupancy-degraded sim
   reads C(1)=0.70, up not down) and statistic-sensitivity check (pilot cell
   coherence floor-dominated at 122-139 um) are accepted as sound.
2. **Composition axis re-owned, not dropped.** The pilot locates the contrast in
   the RESPONSE FLOOR (C(1)~0.71 vs real 0.30) — frozen here by scope guard and
   already the pipeline's principal open front (med05 0.83, pixel LOCAL 0.90,
   content owner-by-exclusion). The axis moves to the response-family campaign,
   which now opens with four independent arrows at one target. The probe
   thread's forthcoming bands should be REDIRECTED there: floor-subtracted /
   decomposed composition meters (granular share; C_floor(d)) as that campaign's
   acceptance instruments — no bands on raw C(1) for v6.1.
3. **Path (ii) — amplitude-only fallback — RECOMMENDED for execution now**:
   confirm SMOD' (~0.009) with one pilot point, then the v6.1 production under
   every pre-declared rule (holdouts, seed rule, naming, exports, figures,
   manifest, briefing). Honesty rider: the raw C(1) on v6.1 will drift up by
   roughly the white-weight shift (order +0.05, bounded by the all-block edge);
   the probe battery measures and RECORDS it as a diagnostic — "best reachable
   configuration; axis owned by the response campaign" is the ledger language,
   not "no regression".
4. FSBLK code: keep behind its off-path guard as a shelved capability — it is
   the right instrument if the response campaign changes the floor's correlation
   structure and the granular re-mix becomes reachable.
5. **USER SIGN-OFF REQUIRED** (clause 6): execution of path (ii) awaits the
   user's explicit word to the pipeline session.

## METER-OWNER ACKNOWLEDGMENT (2026-08-23, Real Data Probe session)

Clause-6 resolution accepted: the composition ruling in the gate response above
is retracted as pilot-falsified on direction (shared error with adjudication #1;
the pilot's mixture, dilution and statistic-sensitivity controls are convincing).
Per the resolution: NO bands on raw C(1)/cell coherence for v6.1 — the unchanged
v61 battery (twist_probe + nonrms_probe + gtail_probe) runs on the v6.1 delivery
with the twist holdouts as GATES and raw C(1)/cell coherence RECORDED as
diagnostics (expected C(1) drift ~+0.05; ledger language per resolution §3).
Floor-subtracted/decomposed composition meters for the response campaign are
being built as island_post/cfloor_probe.C (C_sub(d) with the statistical row-mean
floor removed; white/short/smooth variance shares; thinned-sim sampling control;
floor-subtracted cell coherence) — baseline numbers on real + v6t to follow in
REAL_DATA_PROBES.md, bands proposed there for the response campaign, not here.
PROVENANCE FLAG for the v6.1 launch: island_post/*_v61.root files dated Aug 21
(digi/island91/prodclus) are a stale pre-twist draft already occupying the v61
filenames — remove or rename before the re-balance delivery mints v6.1, or the
md5 trail will collide.
