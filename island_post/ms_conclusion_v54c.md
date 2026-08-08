# Multiple scattering in the sPHENIX TPC — comprehensive conclusion (sim + data, v5.4c)

2026-08-07. Closes the circle-fit / multiple-scattering (MS) question at every level of
both simulation and real data, standing on the v5.4c production
(`island91_frames_production_v54c.root`, md5 5972cbe3…) and the canonical real file
(`clusters_seeds_island_79507-0.root_ntuplizer.root`, run 79507 pp, 100 events,
62 complete by hit-level endpoint; 61 by the cluster-level flag used in the symmetric scan).

## Verdict

**Multiple scattering in the TPC gas requires no explicit treatment, in simulation or in
data, at any level of the reconstruction chain.** The truth trajectory is a circle to
20 µm (5×10⁻⁴ of the 3.69 cm sagitta at 0.5 GeV); MS perturbs the tangent by
1.37 mrad at 0.5 GeV (Highland-order, 1/p-scaling verified); cluster resolution
produces 18–31 mrad of tangent noise on BOTH sim and real tracks — a factor ~20
larger, and pT-flat, which is the signature of resolution rather than scattering.
The real data itself confirms it: the real split-arc mismatch does not grow at low pT.
What DOES need modeling is not scattering but static distortion (d0 ≈ 2.4–2.5 cm
off-origin, 0.48 cm membrane jump) — which is what the v5.4 field overlay addresses.

## The ladder (one estimator family, five levels)

Circle fit = Kasa init + 6 Gauss–Newton iterations everywhere. Split-arc = independent
half fits with tangent mismatch Δψ at the border. Exhaustive finder = sectorized
conformal Hough + drift-coherence RANSAC + acceptance bar (≥12 clusters after 3 mm
cleaning, ≥13 layers, span ≥15, gap ≤6, RMS ≤ 0.20 cm, R ≥ 45 cm), identical
wherever it is used.

| level | sample | key number (real) | key number (sim v5.4c) | data/MC |
|---|---|---|---|---|
| truth trajectory | ntp_g4hit crossers, 0.5 GeV | — | circle RMS 20 µm; R_fit/R_exp 0.9990 | — |
| truth MS | split-arc on truth hits | — | σ(Δψ) 1.365±0.008 mrad @0.5 GeV; 0.367±0.005 stiff | — |
| tracker tracks | ntp_clus_trk vs truth-groups (r0=49) | σ(Δψ) 25.7±0.7 / 31.2±3.0 mrad; RMS 692 / 626 µm | 27.8 / 30.2 mrad; RMS 842 / 635 µm | σ: 0.92/1.03; RMS: 0.82/0.99 |
| all clusters, symmetric finder | full ntp_cluster both sides | 184 circles/ev; RMS 658 µm; in-time 0.34 | 385 circles/ev; RMS 815 µm; in-time 0.33 | rate 0.48; RMS 0.81 |
| pre-clustering | real: ALL pixels; sim: truth hits/collision | 293 circles/ev (median); RMS 805 µm; in-time 0.36 | 14.06 findable/collision; finder recall 0.96 | via bridge (below) |

## Level detail

### 1. Truth trajectory (sim) — the trajectory IS a circle
`truth_circle_v54cx.txt`: 16,626 full R1→R3 crossers at pT 0.45–0.55 GeV fit with
median RMS **20 µm** (p90 64 µm), median R_fit 118.62 cm against R_exp = 119.1 cm
(pT/0.3B at B = 1.4 T): ratio 0.9990. Deviation/sagitta = 5.4×10⁻⁴ over the gas.
Bit-identical to the v5.3 truth record (the v5.4 field is a positions-only overlay on
reconstructed clusters; truth hits untouched).

### 2. Truth-level MS (sim) — the only irreducible trajectory effect
`ms_split_v53f.txt` (truth side unchanged since), border r0 = 35 cm (adopted worst
point; `ms_r0scan_v53f.txt` is the justification):
- σ(Δψ) = **1.365 ± 0.008 mrad** at 0.5 GeV (fit noise 0.08); 0.367 ± 0.005 stiff.
- Cross-window ratio 3.72 vs Highland 3.63 (1/p scaling — it is scattering).
- Highland θ0 full crossing (order scale): 1.508 / 0.416 mrad; gas X0 = 112 m.
- Displacement equivalent ≈ 20–30 µm: 2.7% of cluster resolution.
- Mid-split (r0=49) values 1.073 / 0.307 mrad — quoted as the truth scale inside the
  r0=49 real figures.

### 3. Tracker-track level (real ntp_clus_trk vs sim truth-groups)
`ms_real_split_v54cx.txt` + `ms_real_split_v54cx.png` / `ms_real_showcase_v54cx.png`
(style-parallel to the supervisor-seen sim figure):
- σ(Δψ) real 25.70±0.70 mrad @0.5 GeV vs 31.21±3.00 stiff — **pT-flat**, i.e.
  resolution-driven; sim same shape (27.83 / 30.15). MS at the 1.07/0.31 mrad scale is
  invisible on both sides.
- Circle RMS: real 692 / 626 µm; sim 842 / 635 µm → data/MC **0.82** at 0.5 GeV
  (sim low-pT over-smear, response family, open) and **0.99** stiff (the v5.3-era +33%
  real excess is resolved by the v5.4 field).
- d0 median 2.51 (real) vs 2.37 cm (sim): off-origin scale matched by design (v5.4).
- Basis: 98.4% of ntp_clus_trk (event,seedID) groups are genuine circles
  (`ms_realcheck`: 86.1% clean + 12.3% rescued, 1.6% fake).

### 4. Symmetric cluster level (the grouping-bias-free comparison)
`ms_cluscmp_v54c.txt` + `.png` — the SAME exhaustive finder groups both sides' full
cluster trees (no tracker, no truth): real complete events 184 circles/ev vs sim
385/ev (**rate 0.48** at cluster-count ratio only 0.88, i.e. sim frames are ~2×
more track-organized per cluster — content/organization axis, not a resolution
statement); per-track RMS 658 vs 815 µm (**0.81**, same family as the 0.82 above);
in-time fraction **0.34 vs 0.33** — the time-structure of found circles matches.

### 5. Pre-clustering level
- **Sim, truth hits per collision** (`ms_g4scan_v54c.txt`): 14.06 findable
  track-class circles/collision under the same bar; finder recall **0.96**; raw found
  rate 34.3/collision decomposes into 39% first-finds + 42% second arcs of the same
  particle + 17% mixed-parent chains + 2% below-bar particles (dense-point
  double-counting — the caveat to carry to any raw pixel-level rate). 94% of found
  fits have RMS ≤ 0.05 cm.
- **Bridge**: 333 findable/frame (cluster-level truth calibration) ÷ 14.06/collision
  = 23.7 collisions/frame equivalent — the two sim levels close.
- **Real, ALL pixels** (`ms_pixscan_v54c.txt` + `.png`, 25 complete events): median
  **293 track-class circles/event straight from raw pixels** (1.7× the cluster-level
  rate — expect a large second-arc share per the truth-level decomposition); per-track
  RMS 805 µm (pad centers, no centroiding: cluster 658 µm + quantization); in-time
  0.36; R_fit median 87 cm. Tracks are findable pre-clustering when the (x,y) circle
  and the tbin-vs-layer line are demanded jointly; (x,y) alone is combinatorially
  saturated (purity 0.3% → 75% with drift coherence).
- Sim pixel trees carry no v5.4 field (the overlay is applied at island91 export), so
  the sim pre-clustering anchor is the truth-hit scan, connected to data through the
  cluster-level chain.

## Resolution ladder at 0.5 GeV (the whole story in µm)

MS displacement ~20–30 µm → truth-trajectory circle RMS 20 µm → real pixel-level
805 µm → real cluster-level 658–692 µm → sim cluster-level 815–842 µm. Sagitta to
measure: 21,400 µm over the pad rows. MS is 3% of the measurement noise and 0.1% of
the signal; **treat clusters as points on an ideal circle plus resolution; put the
modeling effort into distortions and response, not scattering.**

## Data/MC state after v5.4c (this campaign's scorecard)

| observable | data/MC | status |
|---|---|---|
| σ(Δψ) @0.5 / stiff | 0.92 / 1.03 | agree ≤8% |
| circle RMS stiff | 0.99 | fixed by v5.4 field |
| circle RMS @0.5 GeV | 0.82 | open — sim over-smears soft tracks (response family) |
| symmetric RMS (all clusters) | 0.81 | same family |
| in-time fraction of found circles | 1.03 | agree |
| d0 median | 1.06 | matched by design (v5.4) |
| track-class circles per event | 0.48 | open — organization/content axis (sim 2.1× at 0.88× occupancy); relay-worthy to the pipeline session |

## Receipts

- Finder calibration (sim truth, 50 frames): efficiency 77.7%, purity 75.2%
  (`missed_tracks_v53f.txt`, regenerated 2026-08-07 with number-identical output —
  doubles as the regression proof that the parameterized hunt() refactor is
  behavior-preserving).
- Real complete-event convention: hit/cluster-tbin p99.9 > 950 (62 hit-level /
  61 cluster-level of 100 events); all real rates quoted on complete events.
- v5.4c acceptance battery (2026-08-06): labels bit-exact v5.3→v5.4c; membrane jump
  real 0.480 vs sim 0.536 cm (inside errors); half-arc curvature width sim 0.276 vs
  real 0.402 (field conservative); multi-segment split-rate 2.26%.
- Figures QA'd per protocol: `ms_real_split_v54cx.png`, `ms_real_showcase_v54cx.png`,
  `ms_cluscmp_v54c.png`, `ms_pixscan_v54c.png`, `truth_circle_v54cx.png`.
