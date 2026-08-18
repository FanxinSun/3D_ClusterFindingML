# Request: laser-flash re-model — the run-79507 flash is ONE event, not 44 frames

2026-08-17, from the "Real Data Probe" session to the pipeline session, relayed by the
user. Everything below is measured on the sole real-data source
`clusters_seeds_island_79507-0.root_ntuplizer.root` (canon cut `layer>=7&&layer<=54&&adc>0`
unless stated) and is reproducible with `island_post/ev44_probe.C`
(`cd island_post && root -l -b -q 'ev44_probe.C+()'` → stdout tables +
`sim_validation_plots/ev44_probe.png`, 4 panels). Companion record: PIPELINE.md section
"REAL-DATA PROBE: the bright dot at (event 44, tbin 330) (2026-08-17)".

Trigger for the probe: the yellow dot at (event 44, tbin ~330) in the recording-extent
map (`trunc_check.C` panel 1). Question asked: laser, or a trigger / tbin=0 artifact?

---------------------------------------------------------------------------------------
## 1. Verdicts (each with the measurement behind it)

**V1 — The tbin-330 spike of the arrival curve is event 44 alone.**
At 1-tbin resolution the all-100 curve holds 64,398 hits at tbin 329 against ~26,800 in
the neighbouring tbins (the 2.4× spike quoted everywhere in PIPELINE.md). Without event
44: 26,784 vs 26,800 → 1.00×. Sum of the other 99 events over 327–331 = 133,660 vs
neighbours 134,040 → residual −0.28 %. The "REAL 62 complete" curve shows the same spike
because event 44 is one of the 62 complete events. So the SIM's injected spike
(2.07×) was tuned against a single event.

**V2 — Event 44 carries a genuine laser flash: 129,767 hits (adc>0) in tbin [322,340],
peak 37,614 at tbin 329**, on a pre-flash baseline of ~200 hits/tbin. Detector-wide:
side 0 / side 1 = 53.4k / 55.8k in 327–332, all 12 sectors on both sides (per-cell
variation ×2), all 48 layers; all 24 side×sector cells >3σ, the top-3 cells carry 17 % of
the excess (uniform = 12.5 %). Per-tbin core: 326:195 327:805 328:22,167 329:37,614
330:29,107 331:14,656 332:4,922 333:3,352 (background-subtracted centroid 330.4, rms
2.65 tbins).

**V3 — Two time markers exactly one full drift apart.** A prompt, detector-wide spike at
tbin 86–89 (+1,368 hits over neighbours; both sides, 10/12 sectors, inner-layer weighted)
= the laser fire time; the central-membrane flash at 329. Difference 242 tbins = the
reco's own tbin→z full drift (its z-mapping puts z=0 at tbin ~244). Apparent z of the
flash (physical convention, `z − 105.5·(zelem==0)`): ≈ −37 cm on side 1 and +34 cm on
side 0, i.e. ~35 cm "beyond" the membrane on each side, because the flash arrives ~85
tbins after a full drift → the laser fired ~85 tbins (~4.5 µs at 53 ns/tbin) after the
trigger. (The composer's CM-model delay of 4346 ns reproduces the arrival tbin; nothing
to change there.)

**V4 — Not a tbin=0 / trigger-start artifact.** No event has a tbin=0 spike (largest
tbin0 / mean(tbin1,2) ratio in the file: 1.30); tbin 0–6 sit flat at ~29.2k hits/tbin
(all events) with a smooth turn-on to 33.5k by tbin 11; no adc==0 hits below tbin 12;
event 44's tbin-0 content (450) is ordinary for its occupancy.

**V5 — It IS trigger-related in one sense: event 44 is a GL1 laser-triggered readout,
and the collaboration's reco recognised it.** GL1 = Global Level-1 trigger, sPHENIX's
central trigger processor (the `raw/live/scaled mbd/zdc` columns of `ntp_info` are its
scalers). Facts from coresoftware master (fetched 2026-08-17):
  - `offline/packages/tpc/LaserEventIdentifier.cc`: per-side histogram of TPC hits vs
    time sample, search range `SetRange(320, m_time_samples_max)`; laser event if
    `(peak/mean >= 7 && peak > 1000)` on either side; for `runnumber > 66153` it also
    sets `isGl1LaserEvent` from the GL1 packet (`getGTMAllBusyVector() & (1<<14)`) and
    `isGl1LaserPileupEvent` when the BCO gap to the previous laser event is short.
  - `offline/packages/tpc/DiffuseLaserEventSelector.cc`: for run > 66153 laser events
    are selected by `isGl1LaserEvent()` alone → the laser flash is a **triggered, GL1-
    tagged event class** in physics running (which is why it sits at a fixed sample).
  - `offline/packages/tpc/TpcClusterizer.cc`: `if (m_rejectEvent && laserInfo &&
    laserInfo->isLaserEvent()) return EVENT_OK;` — no TPC clustering for laser events;
    `macros/TrackingProduction/Fun4All_{Full,PRDF}Reconstruction.C` set
    `G4TPC::REJECT_LASER_EVENTS = true`.
  Applied to run 79507: event 44 passes the identifier by a mile (peak/mean = 63 with both
  sides combined; peak 37.6k) and NO other event passes (next-best peaks in [320,965] are
  1,271 and 1,036 hits with peak/mean 2.6 and 2.9). Consistently, `ntp_info` has for event
  44: `nhittpcall = 535,606` (largest of the run) but `ntrk = ntpcseed = nclustpc = 0`
  (nclusintt 1199, nclusmaps 2487 = normal → the event itself is fine, only TPC clustering
  was skipped).

**V6 — The 43 other events flagged in `laser_frames.txt` have no flash.** That list used
"hits[325–340] − hits[350–365] > 300", i.e. a 16-tbin window against a control 25 tbins
later — it flagged the local slope of each event's arrival curve, not a spike (single-event
arrival curves fluctuate by ±1000 hits/8 tbins from looper/track structure; Poisson would
be ±50). With a narrow statistic (hits[328–331] minus the adjacent 4+4 tbins) every other
event is < 400 and localised: top-3 of 24 side×sector cells carry 57–88 % of the excess
(event 18: one cell holds +340 of +281 net). The per-event tbin profile of the 43 "ordinary
flagged" events across 326–333 is flat (296, 292, 290, 285, 283 hits/tbin/event); their
layer and φ-fold "fingerprints" do not correlate with event 44 (r = −0.11 / +0.08) because
there is nothing there. Consequence: the "44 flashes", the "41-flash size distribution",
the "twin spikes" fine-tune and its "sampling noise of 44 flashes/run" argument all rest on
this list. **Reality: 1 flash in 100 events, giant.** (Poisson: one observation → a
per-event probability anywhere from ~0.2 % to ~3 % is compatible.)

**V7 — Collateral inside event 44 (all measured; none is physics of the light itself):**
  - Saturation: 3,707 pads have a flash sample ≥ 930 ADC (file maximum 963 = 10-bit
    ceiling minus pedestal; pile-up at 940–955); 6,077 (4.7 %) of the 129.8k flash hits
    are saturated. Window ADC quantiles 10/25/50/75/90 % = 21/29/57/198/685 (mean 187) vs
    20/25/38/83/182 outside the window.
  - Tail: the whole-event tail over 333–432 is ~62k hits above the event's own baseline
    (252/tbin), exponential with τ ≈ 26.5 tbins (fit 336–410). It is a per-channel
    recovery of the SATURATED pads: 0.205 hits/pad/tbin over 340–400 on saturated pads,
    0.020 on hard (500–930) pads, 0.0019 on the rest (×108); 3,707 × 0.205 × 60 ≈ 46k
    reproduces the 340–400 part. Not late light.
  - adc==0 hits: 15,817 in event 44 (43 % of the file's 36,573), all at tbins 328–330,
    concentrated in specific sectors (side 0: sectors 1,2,4,5,7 ~1–2k each, sectors 8,9
    zero; side 1: sectors 5,6 ~1.8k, sector 1 zero) and R2-heavy; 70 % sit on pads with no
    other hit within ±3 tbins → not saturation markers but a burst-only unpacker code
    (adc 1–10 hold ~500 hits each in the whole file, adc 0 holds 36.6k). In the other 24
    affected events they sit at random tbins (525, 595, 665, 790 …), so the old note
    "adc=0 concentrated at the flash" was event 44 too. Harmless: islandize drops adc==0.
  - Layer pattern of the flash: window/outside ratio per layer 0.14–0.61, rising with
    radius, with an odd/even alternation in R2/R3 (odd layers ≈ 2× even; e.g. L25/26/27/28
    = 4590/2214/5462/2215 window hits) — stripe geometry, already the CM model's domain.
  - The pre-/post-flash "deficits" of event 44 (250–320 at 0.58×, 400–480 at 0.67× of its
    late level) are ordinary pileup scatter: z = −0.8 / −0.9 among the 61 other complete
    events (window ratios scatter with rms 0.4–0.8, range 0.3–5).

**V8 — What stays valid.** `canon.h TPC_CUT_NOLASER` as a symmetric fiducial time cut
(it now simply drops event 44's flash + the first part of its tail on the real side and
the injected flash on the sim side); the truncation census (`trunc_check.C`, 62/38);
cluster-level real references (`ntp_cluster`, `prodclus_real_*`) never contained event
44's flash (the reco vetoed the event); hit-level real anchors (`ref_real.root`,
`island_real*.root`, `real_complete62`) DO contain event 44's 130k flash + ~62k tail +
15.8k adc==0 hits (0.5 % of the file's hits) — see §3.

---------------------------------------------------------------------------------------
## 2. Where the current flash model lives (what the re-model touches)

  - Composer call (pp_pipeline.sh line ~215; island_post/run_production.sh line ~23):
    `frame_composer(..., FLASH="raw_lib_cmflash_w.root", flash_prob 0.44, flash_scale 1.0,
    SPEC, jitter 1.5, pixdisp 2.0, stripedisp 1.1, blur 0.75, RSPEC, ...)` with
    `SPEC="0.008:1,0.009:1,0.011:1,0.012:1,0.013:1,0.014:1,0.018:1,0.021:1,0.027:1,0.037:1,2.2:0.25"`
    → 44 % of frames get a flash; of those 10/10.25 get a knee-level (0.008–0.037) flash and
    0.25/10.25 = 2.4 % the giant (2.2) → **~1.07 % of frames get the giant — that part
    already matches reality; the ~43 % of frames with knee-level flashes are the fiction.**
  - `island_post/laser_frames.txt` (flag list) — invalid as a flash list.
  - `island_post/canon.h` comment on TPC_CUT_NOLASER: "44/100 frames flagged in
    laser_frames.txt" → should read "1/100 events (event 44, GL1 laser-triggered)".
  - `island_post/laser_assess.C` / `sim_validation_plots/laser_assessment.png` (real vs CM
    model arrivals + φ-fold): the "REAL" curves there are event 44 in disguise; the
    per-flash comparison logic (spike ratio 2.07 vs 2.4, jitter fine-tune, "sampling noise
    of 44 flashes") is a one-event statistic → retire the spike-ratio gate.
  - Flash region weights in `raw_lib_cmflash_w.root` (0.36/0.87/0.43, "direct flash-only
    giant-yield inversion") were fitted against the giant, i.e. against event 44 → probably
    still valid; the pipeline agent can re-check against the numbers in V2/V7.
  - Truth labelling: flash pixels/islands → cls=2 (composer draw slot 0xFF / trk −1
    conventions) — mechanism unchanged; only its frequency is wrong (44 % of frames vs ~1 %).

---------------------------------------------------------------------------------------
## 3. The ask — two consistent options (decision open; recommendation given)

**Option A — model what the DAQ records (recommended).** Our ML dataset lives at the
hit/island level (islandize on `ntp_hit`), and the real hit-level data DO contain event
44's flash. So keep the CM-flash injection but at the real frequency and shape:
  1. `flash_prob` 0.44 → ~0.01 (defensible range 0.002–0.03), `SPEC` → giant only
     (the "2.2" bucket, or whatever scale reproduces the V2/V7 numbers): 129.8k adc>0
     hits in [322,340], ~4.7 % of them saturated (3.7k pads), ADC quantiles
     21/29/57/198/685, core per-tbin 805/22.2k/37.6k/29.1k/14.7k/4.9k at 327–332, side
     balance 53.4k/55.8k, sector spread ×2, layer alternation as in V7.
  2. Decide whether the readout stage should produce the saturated-channel recovery tail
     (0.205 hits/pad/tbin over 340–400 on saturated pads, τ ≈ 26 tbins, ~62k hits per
     flash event). The present ion tail (f 0.021, τ 7 tbins + saturation re-clamp) is far
     too short to make it; it affects ~1 frame in 100 — the pipeline agent's judgment.
     The prompt spike (+1.4k hits at tbin 86–89, inner-layer weighted) is optional.
  3. Do NOT model adc==0 (islandize drops it).
  4. Keep flash truth = cls=2 as now; expect the flash share of cls=2 to fall by ~×40.
     Optionally emulate the collaboration's tag on sim frames (per-side peak > 1000 with
     peak/mean ≥ 7 over samples ≥ 320) so a "laser event" flag travels with the frame.

**Option B — mimic the collaboration reco.** No flash injection at all (`flash_prob 0`)
and, on the real side, treat event 44 as a vetoed laser event (drop it from hit-level
anchors, or keep it flagged). The collaboration's production has no TPC clusters for
such events; cls=2 stays populated by the isolated-hit background (src 0xFE) + minis.

**Either way (housekeeping):** retire `laser_frames.txt` as a flash list; fix the canon.h
comment; drop the spike-ratio acceptance; re-frame or retire `laser_assess.C`; if hit-level
real anchors are rebuilt, decide event 44 in/out consistently with the option chosen (the
present anchors include it; the earlier "flash = 0.4 % of pixels, soft" note undercounts:
event 44's flash+tail+adc==0 ≈ 208k hits ≈ 0.75 % of the file, and it is not soft — window
mean ADC 187).

Nothing above was changed by the probe session: `canon.h`, `laser_frames.txt`, the
composer, the productions and all anchors are untouched. New files only:
`island_post/ev44_probe.C`, `sim_validation_plots/ev44_probe.png`, this request, and the
PIPELINE.md section.

---------------------------------------------------------------------------------------
## 4. How to verify a re-model (acceptance the probe would run)

  1. Per-frame narrow statistic on the sim (hits[328–331] − ½(hits[324–327]+hits[332–335])
     under the sim's tbin): every non-flash frame < 500; flash frames ~1 % of frames,
     each detector-wide (all 24 side×sector cells positive) with a core ≈ 100k hits.
  2. Ensemble arrival curve: the sim's spike height is then a Poisson-1 quantity — compare
     the flash-FRAME fraction and the per-flash shape (V2/V7 numbers), not the ensemble
     spike ratio.
  3. Real side unchanged: `ev44_probe.C` numbers above (129,767 / 37,614 / −0.28 % / 24/24).
