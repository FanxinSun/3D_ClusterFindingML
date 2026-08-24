# Real Data Probes — run 79507 (this thread's standing record)

Owner: the "Real Data Probe" session. Scope: probing the ONLY real-data source,
clusters_seeds_island_79507-0.root_ntuplizer.root (canon cut layer 7-54 && adc>0
unless stated). Convention (user directive 2026-08-19): probe records live HERE,
not in PIPELINE.md (that is the pipeline agent's ledger); anything that needs
pipeline ACTION goes out as an island_post/*_request.md handover doc, relayed by
the user. Producers are persistent macros in island_post/ with per-probe ledgers.

---------------------------------------------------------------------------------

## 2026-08-17 — the bright dot at (event 44, tbin ~330): the ONE laser flash

Producer: island_post/ev44_probe.C -> sim_validation_plots/ev44_probe.png.
Headline: the tbin-330 arrival spike is event 44 ALONE (129,767 TPC hits in
[322,340], peak 37,614 at 329; the other 99 events summed are flat to -0.28%).
Event 44 = a GL1 laser-triggered readout: prompt fire-time spike at tbin 86-89 +
CM flash at 329 (one full drift, 242 tbins, apart); LaserEventIdentifier tags it
(peak/mean 63 vs threshold 7) and TpcClusterizer skips it (ntp_info: 535,606 TPC
hits, 0 TPC clusters). NOT a tbin=0 artifact. laser_frames.txt's 43 other
"flashes" were arrival-curve slopes (16-tbin window vs a control 25 tbins later).
Collateral: 3,707 saturated pads (adc 940-963) whose recovery is the 340-400
tail (~62k hits, tau ~26.5 tbins); 15.8k adc==0 hits at 328-330 (burst-only
unpacker code, dropped by islandize).
ACTION (needed -> handed over): island_post/laser_flash_remodel_request.md
(verdicts V1-V8, options A/B). Consumed by the pipeline thread: PIPELINE.md
"REAL-DATA PROBE ... (2026-08-17)" + "v5.6 CAMPAIGN START" -> the V6 campaign
(FLASH_PROB 0.44 -> 0.01, giant-only spec, laser_frames.txt retired). That
PIPELINE.md section stays where it is as the consumed request's record; every
later probe record lives here instead.

---------------------------------------------------------------------------------

## 2026-08-19 — the nf_digipix GLOBAL long tail = road pickup (no pipeline action)

Question: why does the real GLOBAL whole-track pixel fit (ms_nofinder_digipix_v6,
med 2358 um) drag a flat tail to 6 mm while LOCAL is 1205 um and the same seeds at
cluster level fit to 780 um? Producer: island_post/gtail_probe.C ->
gtail_probe_v6.txt + sim_validation_plots/gtail_probe_v6.png (4 panels; nf_digipix
groups rebuilt exactly — road dxy<1.2 cm, |dtbin|<=6, ev44 veto — same MNF fitter).

ANSWER: the tail is the ASSOCIATION ROAD, not the tracks. One symmetric robust clip
(3x1.4826xMAD about the median residual, 3 iterations, floor 500 um) costs a median
3.1% of real pixels and 0.0% of sim pixels and gives:
  real med 2358 -> 1590 um | sim 1538 -> 1534 | data/MC GLOBAL 1.53 -> 1.04.
Sim v6 digi CARRIES the v5.5 field and loses nothing to the clip -> the clip does
not eat field content; what it removes on the real side is other tracks' pixels
inside the 1.2 cm x (+-6 tbin) road (no outlier cleaning by design of the RAW bar).
  - Tail (raw>3mm) = 29.2% of real fits vs 1.4% sim. Of the real tail: 57% trims to
    clean (med 14% of pixels dropped; max within-row peak-to-peak residual med
    13.3 mm = a second blob in the same pad rows — no distortion moves charge 13 mm
    within a row); the 43% survivors carry the same 12-15 mm two-blob rows
    (balanced blends the MAD clip cannot split; 20% of ALL real tracks lose >10% of
    pixels) plus a genuine low-pT non-circular population that sim shows too (sim
    tail R med 60-65 cm = loopers; >=6mm overflow 0.6% on BOTH sides).
  - Pileup fingerprint: tail fraction 23.0/29.8/34.9% across event-occupancy
    terciles raw; 5.9/6.5/8.0% after the clip.
  - Why only GLOBAL blows up: LOCAL is a median over 45 sliding windows (a
    contaminated row hits few windows), the cluster-level fit uses tracker-CLEANED
    clusters (780 um), the global pixel fit is one fit over everything.
    Quadrature closure: sim 1337(local) + LF 720 = 1518 ~ 1538 raw (closes); real
    1205 + 780 = 1435 << 2358 raw (does NOT close) but after the clip
    1205 + LF 914 = 1512 ~ 1590 (closes).
WHAT REMAINS after cleaning (the genuine real-only content, small): global surplus
sqrt(1590^2-1534^2) = 418 um; trimmed-LF surplus ~550 um at every span
(real 800/864/925 vs sim 589/657/762 um for spans 15-25/25-35/35-45 cm); cell
coherence (mean row residual per side,sector,layer, same estimator both sides):
real 307 um vs sim 122 um (~2.5x the visible v5.5 SMOD imprint). IMPLICATION for
the meters: GLOBAL pixel data/MC 1.53 is an association meter, not a field meter —
quote the clipped 1.04 (or a cleaned road) when the question is detector/field
content. nf_digipix/ms_nofinder.C left untouched; probe is standalone.

---------------------------------------------------------------------------------

## 2026-08-19 — the nf_digipix LOCAL bimodality = window-statistics mixture (no pipeline action)

Question: the real LOCAL 4-row sagitta RMS (ms_nofinder_digipix_v6 right panel) is
bimodal — shoulder ~500-800 um under the ~1100 um peak. Producer:
island_post/lsag_probe.C -> lsag_probe_v6.txt + sim_validation_plots/lsag_probe_v6.png
(windows rebuilt exactly as nf_digipix::doTrack; same fitter/gate/road/ev44-veto).

ANSWER: the shoulder is the SPARSE-WINDOW population, a statistics artifact of the
fit, not a detector feature.
  - Window RMS scales with pixel count n (both sides, same trend): real med by
    n-class 754 (n=5-6) / 1000 (7-9) / 1210 (10-14) / 1579 um (>=15); the two lowest
    classes carry 77% of the shoulder [300,900) and peak there themselves.
  - Real is bimodal and sim is not because of the n MIX: real windows med n 12
    (3.25 px/row), sim med n 24 (6.67 px/row) -> n<=9 windows are 35.3% of real vs
    1.8% of sim. 3-row windows are similar (34 vs 29%) — the sparsity is pixels per
    row (real clusters carry ~half the digi pixels of sim; same 2x as px/track
    124 vs 266), not missing rows.
  - CLOSURE: reweighting the sim windows to the real n-distribution reproduces the
    real shape — shoulder frac 23.7% vs real 21.9% (raw sim 4.9%), overlay traces
    the two-mode outline (fig P1). No hidden hardware needed.
  - Ruled out directly: side (med 1207/1203), sector (1132-1250 flat map), drift
    (1187-1212 flat in tbin), pad region (R1->R3 gradient 1019->1324 exists but is
    the SAME shape in sim 1278->1455, steps at R1/R2 and R2/R3 in both — fig P3),
    two-blob road pickup (rowPP>5mm windows: med 1948 um, 0.0% in the shoulder —
    they are the RIGHT tail + the 4000-overflow spike, 24.6% of real windows).
  - The near-zero bin (~1% of windows at <70 um) = n=5-6 overfit (5 points, 3 dof);
    reweighted sim shows it too.
METER IMPLICATION: the LOCAL data/MC 0.90 carries the n-mix inside it. Same-footing
versions: dof-corrected (rms/sqrt((n-3)/n)) med 1422 vs 1435 -> 0.99; matched clean
cell (rowPP<=5mm, n 10-14) 1145 vs 1069 -> 1.07. Honest response band ~0.99-1.07,
i.e. real clouds are NOT 10% narrower — the 0.90 was mostly sparse-window dof.
PIPELINE-FACING FLAG (user to relay if wanted): at the digi level sim makes ~2x the
above-threshold pixels per crossing of real (6.67 vs 3.25 px/row at matched window
population) — a content-level cloud/ZS difference the local meter was hiding.
ms_nofinder.C untouched; probe standalone.

---------------------------------------------------------------------------------

## 2026-08-20 — "a surprise of only halved": the LOCAL fit's floor is the charge
## cloud, and it sits within 6% of it (no pipeline action)

Question (supervisor, via user): the 4-row sagitta fit reduces the real median
only 2358 -> 1205 um — shouldn't it go MUCH lower? Producer:
island_post/lfloor_probe.C -> lfloor_probe_v6.txt +
sim_validation_plots/lfloor_probe_v6.png (same groups/fitter as nf_digipix;
global 3xMAD trim keeps road junk out of the floor measurements).

ANSWER: no — 1205 um is already the floor. Two facts:
  (1) "Only halved" is the linear-scale illusion. RMS combines in quadrature:
      the LOCAL fit removed sqrt(2358^2-1205^2) = 2026 um = 74% of the variance
      (road pickup + field + curvature). Sim for contrast removed 25%
      (1539 -> 1337): its global was already cloud-dominated.
  (2) What remains is NOT trajectory error. Per-window decomposition:
      rms med 1205 = HF (within-row pixel scatter) 1133 (+) LF (row-mean
      residuals) 219 um. The 4-row fit (3 dof on 3-4 rows) has already zeroed
      the sagitta content — a PERFECT trajectory model would still read ~1.1 mm,
      because the fit runs on raw pixels: pad-center positions repeated per tbin
      sample. Fit-free cloud width per row: real med 1105 um = 0.49 pad pitch
      (2 pads/row x 1.5 tbin samples/pad; ADC-weighted 997 um); the width
      distribution shows the discrete 2-pad / 3-pad comb in real AND sim
      (fig P3; the zero spike = single-pad rows). Shrinking the window further
      cannot cross this; it only loses dof (the lsag_probe shoulder).
HOW TO GO LOWER — change what a "point" is, not the window: ADC-centroid each
row and refit the whole track: same pixels -> real med 953 um, sim 809
(fig P4; tracker's own clusters: 780). Closure: 953 ~ centroid resolution (+)
the genuine long-range content (gtail trimmed LF 914 um) -> after centroiding,
what remains is resolution + field, i.e. the thing the supervisor wanted to
see. A LOCAL centroid sagitta cannot be quoted per 4-row window (3-4 points vs
3 dof leaves 0-1 residual dof), which is exactly why the pixel-level local
meter exists — and why its floor is the cloud width.
Cloud-width aside (feeds the same flag as the lsag probe): real 1105 um =
0.49 pitch (2 pads/row) vs sim 1301 um = 0.69 pitch (3 pads/row), but
ADC-WEIGHTED widths nearly match (997 vs 1056 um, 6%) — the sim-real cloud
difference is mostly marginal pads above threshold, not physical width.
ms_nofinder.C untouched; probe standalone.

---------------------------------------------------------------------------------

## 2026-08-20 — before/after the 4-row sagitta fit: metrics beyond one RMS ratio
## (no pipeline action)

Question (user): the before/after (2358 -> 1205) is assessed only as an RMS
ratio — what else? Legitimate critique sharpened by the probe: the two numbers
are UNPAIRED medians in different units (4,836 tracks vs 134,742 windows).
Producer: island_post/rscale_probe.C -> rscale_probe_v6.txt +
sim_validation_plots/rscale_probe_v6.png (same groups/fitter/gates as
nf_digipix; window gate at every L: >= max(3, L/2) rows, >= 5 px).

Four before/after metrics now measured:
  (1) RMS(L), the scale curve — median per-fit RMS vs window length L = 4..48
      rows, raw and per-window 3xMAD-clipped, real and sim:
        L:        4     8    16    32    48
        real raw  1205  1336  1575  2092  2330
        real clip 1165  1283  1392  1480  1587
        sim raw   1337  1430  1487  1523  1556
        sim clip  1321  1422  1480  1519  1552
      Reading: sim clip == sim raw at EVERY L (nothing to remove at any scale);
      real raw-clip gap grows monotonically 40 -> 743 um (road pickup
      accumulates with lever arm); real-clipped converges onto sim at large L
      (1587 vs 1552, +2%) from BELOW at small L (1165 vs 1321: real's narrower
      clouds). The single 4-vs-48 ratio is just the two endpoints of this curve.
  (2) PAIRED per track (same units at last): global med 2356 vs
      own-4-row-median med 1187; the 2D map shows the mechanism: a horizontal
      band at l ~ cloud width stretching to g = 6 mm = road-contaminated tracks.
  (3) Reduction-factor spectrum g/l per track: real med 1.82, q90 3.48
      (sim 1.15, q90 1.29); real tracks with no reduction: 0.8% (sim 5.7% —
      short tracks already at floor).
  (4) Removed content per track sqrt(g^2 - l^2): real med 1930 um, broad
      1-4 mm (sim 755 um, narrow) — the per-track distribution of exactly the
      quantity the one-number ratio compresses away.
These four separate the three contents by construction: cloud (small-L
plateau), road (raw-clip gap vs L), field/curvature (clip-sim gap vs L) —
answering "any other metrics" with measurements rather than a list.
ms_nofinder.C untouched; probe standalone.

---------------------------------------------------------------------------------

## 2026-08-20 — the same before/after metrics WITH the exhaustive finder:
## real and sim become identical (no pipeline action)

Question (user): were the before/after-sagitta probes with or without the
exhaustive finder? Answer: everything above ran WITHOUT it (ms_nofinder = the
no-finder branch: real grouped by the tracker road, sim by truth). Redone WITH:
MTK::hunt (missed_tracks.C: conformal Hough + drift-coherence RANSAC +
circularity bar) run on ALL clusters of BOTH sides (real ntp_cluster, 99
events; sim island91 v6, 50 frames — same finder config as ms_cluscmp), pixels
attached to member clusters by the same nf_digipix road, then the identical
rscale metrics. Producer: island_post/fscale_probe.C -> fscale_probe_v6.txt +
sim_validation_plots/fscale_probe_v6.png. Tracks: real 15,804 found -> 13,488
pass the pixel bar; sim 19,049 -> 16,921 (381/frame, matches ms_cluscmp 377).

RESULT: with symmetric grouping the before/after behavior of real and sim is
THE SAME within ~5-7% at every scale — and the no-finder asymmetry inverts:
   RMS(L) med [um]     L=4    8     16    32    48
   real raw           1357  1445  1613  1916  2063
   sim  raw           1466  1612  1813  2069  2200
   real clip          1338  1406  1461  1516  1545
   sim  clip          1437  1561  1631  1670  1705
   (no-finder was: real raw 1205->2330 vs sim 1337->1556, ratio 1.50 at L=48)
Paired: g/l med real 1.35 vs sim 1.34 (q90 2.36/2.18); removed content med
real 1341 vs sim 1408 um; the g/l and removed-content DISTRIBUTIONS overlap
(fig P3/P4).

READING: the dramatic real-vs-sim GLOBAL asymmetry of the no-finder branch
(data/MC 1.53; g/l 1.82 vs 1.15; removed 1930 vs 755 um) was the ASYMMETRIC
GROUPING, not the detector: grouped by the same algorithm with the same road,
sim develops the identical association pickup (sim raw 2200 -> clip 1705 at
L=48 - road junk it never had under truth grouping). Road pickup is a property
of the association method; truth grouping is what removed it on the sim side.
Residual with-finder differences are small and go the OTHER way: real sits
5-10% BELOW sim at all L (real narrower clouds; sim island91 cluster mix).
CAVEATS: the finder population is circularity-selected by construction (both
sides equally), 2.6x larger than the tracker-seed set (15.8k vs 6.1k) and
~2/3 out-of-time pileup tracks; absolute levels are not comparable 1:1 to the
no-finder numbers - the with/without comparison is about SYMMETRY, not scale.
ms_nofinder.C / missed_tracks.C untouched; probe standalone.

---------------------------------------------------------------------------------

## 2026-08-20 — road scan: the raw GLOBAL meter tracks the road, the clipped
## one does not (no pipeline action)

Question (user): "scan with the road see what will happen." Producer:
island_post/roadscan_probe.C -> roadscan_probe_v6.txt +
sim_validation_plots/roadscan_probe_v6.png. Method: one pass collects each
real pixel's candidate seed clusters within the widest road (1.8 cm, +-10 tb)
in bucket order; each (RXY, DT) config replays the exact nominal assignment
offline (first-in-road wins, hitID-dedup break semantics) — bit-exact
reproduction at the nominal (1.2, 6): Graw 2358 / Gclip 1590 / LOCAL 1205 /
124 px/track / 4836 tracks. Sim reference = truth-grouped digi (no road),
recomputed: GLOBAL 1539 / LOCAL 1337. Grid: RXY {0.3,0.6,0.9,1.2,1.8} cm x
DT {3,6,10} tbins.

RESULTS (DT=6 row; full 15-config table in the ledger):
   RXY [cm]        0.3    0.6    0.9    1.2*   1.8
   GLOBAL raw     1434   1691   1975   2358   3458   (data/MC 0.93 -> 2.25)
   GLOBAL clipped 1397   1558   1578   1590   1631   (data/MC 0.91 -> 1.06)
   pickup          325    658   1187   1741   3050
   LOCAL          1085   1168   1189   1205   1239   (data/MC 0.81 -> 0.93)
   px/track        108    118    121    124    130
Across the whole grid: data/MC GLOBAL raw spans 0.93-2.70; clipped spans
0.90-1.15. Track count is road-independent (4792-4853).

READING:
  1. The raw GLOBAL number is an association artifact almost entirely: the
     road knob alone dials data/MC from 0.93 to 2.70. The nominal 1.2 cm road
     buys only +15% pixels over 0.3 cm (clusters are mm-wide; there is nothing
     genuine left to collect beyond ~0.6 cm) while adding 1.4 mm of pickup in
     quadrature.
  2. The 3-sigma-clipped GLOBAL is road-INVARIANT (1397-1631 um over a 6x
     road-radius range): the clip recovers the same physics answer whatever
     the road — the road-robust way to quote the meter. Even at the tightest
     road ~325 um of pickup remains (overlapping clusters within 3 mm —
     unavoidable merges the clip still handles).
  3. The DT (time) window matters as much as dxy: at 1.2 cm, pickup goes
     1347 -> 1741 -> 2161 um for DT 3 -> 6 -> 10 (z-neighbours).
  4. LOCAL is mildly road-dependent too (1085-1239): a tight road trims
     cluster edges and biases the cloud width down (data/MC local 0.81 at
     0.3 cm) — so the local "response" ratio also carries the road inside it;
     the dof-corrected / matched-cell numbers (0.99-1.07, lfloor/lsag probes)
     remain the honest response statement.
RECOMMENDATION (meter, not pipeline): do not tune the road to beautify the
raw number — quote the clipped GLOBAL (road-invariant) or use symmetric
grouping (fscale_probe); if a raw pixel meter must exist, (0.6 cm, +-3 tb)
is the sweet spot (97% of the pixels, 1/3 of the pickup).
ms_nofinder.C untouched; probe standalone.

---------------------------------------------------------------------------------

## 2026-08-21 — non-RMS metric battery: two new findings (no pipeline action)

Question (user): the r/fscale comparisons are still all fit-RMS — really no
other metrics? Producer: island_post/nonrms_probe.C -> nonrms_probe_v6.txt +
sim_validation_plots/nonrms_probe_v6.png. Same nf_digipix groups (real =
nominal road, sim = truth digi), trimmed fits; four non-RMS axes:

  (1) POOLED signed residual shape (600,976 real / 6.57M sim pixels):
      |res| q68 real 1.88 mm = sim 1.88 mm EXACTLY — the response cores are
      identical as *shapes*, not just medians. Divergence only in tails:
      q95 6.38 vs 3.88, q99 10.9 vs 8.6 mm; frac |res|>3mm per track 11.9%
      vs 4.2%; per-track skewness +0.12 vs 0.00.
  (2) PULLS (res / track robust sigma, dimensionless): cores unit-width both
      sides; frac |pull|>3 per track: real 3.0% vs sim 0.0% — the scale-free
      association meter.
  (3) ROW-LAG AUTOCORRELATION of row-mean residuals — NEW FINDING: the
      correlation LENGTHS disagree. C(d)/C(0): real 0.30/0.18/0.08/-0.06 at
      d=1/2/4/8 (zero by ~6 rows) vs sim 0.77/0.67/0.55/0.25 (positive to
      ~15). Real's coherent content is SHORT-RANGE/granular (per-row/region
      offsets; a positive bump at lag 16 = the per-region row count, R1/R2/R3
      are 16 rows each); sim v6's baked field (SPHI harmonics et al.) is
      SMOOTH/long-range. Same magnitude class, wrong frequency content —
      invisible to any RMS. (Caveat: real row means are noisier, 3.25 vs 6.67
      px/row, which suppresses the real AMPLITUDE uniformly at all lags — but
      the decay length / zero crossing is dilution-invariant, so the
      short-vs-long-range contrast is robust.)
  (4) SPLIT-HALF CURVATURE (parameter-level, blind to the cloud floor) — NEW
      FINDING: real has a SYSTEMATIC signed offset. Dsagitta = (k_in - k_out)
      x span^2/8 over rows 7-30 vs 31-54: real median +6.5 mm, robust sigma
      9.1 mm (n=4384); sim median +0.2 mm, sigma 5.2 (n=21458; same estimator
      -> unbiased). ATTRIBUTED (2026-08-21, charge-resolved rerun): the offset
      is the SAME for both bending directions (bend<0 +6.05 mm n=2024, bend>0
      +7.12 mm n=2360; sim +0.3/+0.1) -> a coherent, charge-independent TWIST:
      an r-phi displacement field whose radial derivative differs inner vs
      outer half (space-charge / E x B class). NOT rowdr-the-radial-part: the
      V6 real-radius bake (tpc_geom_table.txt = GDML + rowdr, R1-only ramp
      +1.2 -> +0.4 mm) is common to BOTH sides' pixel geometry, so it cannot
      produce a real-vs-sim difference; the twist is rowdr's missing r-phi
      counterpart. The v5.4/v5.5 field model has no such term (SMOD random,
      SPHI phi-harmonics, SCM membrane step - all width-level anchors; this
      is the first SIGNED anchor). |d0|: med real 2.45 vs sim 2.56 cm (cores agree) but q90
      10.8 vs 6.7 cm (real far-d0 tail: secondaries/pileup/fakes).

SUMMARY: the non-RMS axes confirm the RMS-based story where they overlap
(cores identical, tails = association) and add two things RMS cannot see:
the field/alignment content's correlation length differs real-vs-sim
(granular vs smooth), and real carries a signed inner-vs-outer curvature
systematic (+6.5 mm sagitta median). ms_nofinder.C untouched; probe standalone.

---------------------------------------------------------------------------------

## 2026-08-21 — the twist profile: a per-region azimuthal SAWTOOTH, closed by
## injection -> handover island_post/twist_field_request.md (PIPELINE ACTION REQUESTED)

Producer: island_post/twist_probe.C -> twist_probe_v6.txt + twist_profile_v6.txt
(payload) + sim_validation_plots/twist_probe_v6.png. Measured the mean azimuthal
displacement D(rphi) = r * wrap(phi_pix - phi_fit(r)) of kept pixels vs pad row per
side (nf_digipix groups, trimmed fits). REAL: sawtooth locked to the module regions,
same sign both sides: R1 +769/+922 um (row 7) -> -1268/-1408 (row 22); R2 +449/+430
-> -568/-561; R3 +982/+1022 -> -692/-1011; jumps R1->R2 +1.7/+1.9 mm, R2->R3
+1.3/+1.8 mm; per-sector R1 means all positive (+182+-85 / +265+-115 um). SIM digi
v6: flat. Same sign on both sides -> geometric (not E x B); not rowdr (radial ramp
already common via the V6 bake). Interpretation: fit-orthogonal image of per-region
azimuthal offsets (+ within-region slope).
CLOSURE (sim pixels displaced in memory by the real-sim table): split-half Dsagitta
+0.21 -> +6.27 mm (real +6.53; per side +6.18/+6.38 vs +6.04/+7.15); profile
RMS-diff to real 447 -> 126 um. 3-level step models: +3.8/+4.4 mm, 327/311 um —
bulk but not all -> the table is the spec.
HANDOVER written: island_post/twist_field_request.md (payload, injection recipe
dphi = delta*1e-4/r per side/row at digi, acceptance = twist_probe + nonrms_probe
re-run, decisions left open). First SIGNED field anchor; V6 not invalidated.

---------------------------------------------------------------------------------

## 2026-08-23 — TWIST ACCEPTED on the probe meters (v6t) + post-injection C(d)
## + SMOD ruling: granularity/composition, not amplitude-only

Context: pipeline delivered V6-TWIST (digi_frames_production_v6.root re-cut in
place, md5 ba837766, manifest 2026-08-23; side review verified). Three gated
items measured here, independently, tags v6t (pre-twist _v6 ledgers/figures
kept as era records): twist_probe_v6t.txt, nonrms_probe_v6t.txt,
gtail_probe_v6t.txt (+ figures *_v6t.png).

(a) FORMAL ACCEPTANCE — PASS on all gates of twist_field_request.md:
    split-half Dsagitta sim +6.31 mm (side0 +6.24, side1 +6.37) vs real +6.53
    (+6.04/+7.15): per-side deltas 0.20 / 0.78 mm, inside +-1 mm (and inside
    the width-request's +6.3 +- 0.3 holdout); matches the delivery-run values
    independently. Profile RMS-diff to real 127 um (gate ~150; was 447);
    residual boundary jumps +84/+165 (s0), +90/+138 um (s1) — were ~1700-1850.
    Charge independence preserved (+5.72/+6.85). Holdouts: pooled q68 = 1.88 mm
    bit-stable BOTH sides; clipped GLOBAL moved 1.04 -> 1.01 — TOWARD real, as
    it must: the twist supplies previously-missing real LF content (real-only
    global surplus 418 -> 169 um; LF surplus 563 -> 426 um). Raw GLOBAL
    1.53 -> 1.49 (cosmetic). Note: the twist_probe "sim+injected" column on
    v6t is a double-injection diagnostic (+8.19 mm) — not a gate.
(b) POST-INJECTION C(d): sim 0.77/0.67/0.55/0.25 -> 0.57/0.50/0.39/0.22 at
    d = 1/2/4/8 (real 0.30/0.18/0.08/-0.06). Direction right, ~40% of the gap
    closed at short lags; sim zero-crossing still ~12-16 rows vs real ~6. New
    sim C(24) = -0.32 (real -0.02): the sawtooth's long-lag anticorrelation
    shows UNDILUTED in sim because sim lacks real's granular share — same
    diagnosis, not an injection artifact.
(c) RULING for the SMOD solve: **granularity/composition-change — amplitude-only
    is insufficient and moves the composition meters the wrong way.**
    Post-twist composition: C(1)/C(0) sim 0.57 vs real 0.30 (noise-corrected
    coherent signal fractions ~0.76 vs ~0.48); cell coherence real 307 vs sim
    145 um (was 122; real granular surplus ~ sqrt(307^2-145^2) ~ 270 um);
    remaining LF surplus sqrt(914^2-809^2) ~ 426 um, correlated over ~2-3 rows.
    Sim's EXCESS width (stiff 728 vs 626) sits in the SMOOTH share while sim's
    granular share is already BELOW real's. Shrinking SMOD (white, per-layer)
    closes the width gates but raises C(d)/C(0) and lowers cell coherence —
    both away from real. The data calls for either (i) landing the down-retune
    on the smooth terms (SPHI) — in tension with the request's SPHI freeze and
    its d0 anchors, their pilot to reconcile — or (ii) re-modeling SMOD with
    ~2-3-row correlation (per-(side,sector,row-block) granularity), which can
    carry a larger share of the width budget with real's correlation shape.
    If v6.1 proceeds amplitude-only regardless (pre-registered gates are
    width-axis and will pass), record explicitly that the composition axis
    (C(d), cell coherence) remains open and slightly worsens.
Gate response appended to island_post/width_rebalance_request.md.

---------------------------------------------------------------------------------

## 2026-08-23 — clause-6: composition ruling retracted (pilot-falsified);
## floor-subtracted meters built and baselined (response-campaign instruments)

The clause-6 resolution (width_rebalance_request.md) is accepted: the 2026-08-23
composition ruling above is RETRACTED — the pipeline pilot falsified its
direction (block-for-white RAISES C(1); reachable band 0.63-0.67 at fixed width;
occupancy-degradation control C(1)=0.70, up not down). Error shared with the
side review's adjudication. Per the resolution: no bands on raw C(1)/cell
coherence for v6.1 (recorded diagnostics only, expected C(1) drift ~+0.05);
the composition axis moves to the response campaign; meter-owner acknowledgment
appended to width_rebalance_request.md (with a provenance flag: stale Aug-21
*_v61.root draft files occupy the v6.1 filenames — rename/remove before the
re-balance delivery).

NEW INSTRUMENTS (island_post/cfloor_probe.C -> cfloor_probe_v6t.txt +
sim_validation_plots/cfloor_probe_v6t.png), baselined on real + v6t + a
thinned-sim sampling control (sim subsampled to real's px/row, keep p=0.49):
  - Statistical floor: real 60% of raw C(0) vs sim 21% -> after subtraction the
    lag-1 contrast VANISHES: C_sub(1) real 0.76 vs sim 0.72 (thinned 0.76) —
    the raw 0.30-vs-0.57 gap was the floor, exactly the pilot's thesis,
    confirmed on this harness by its own estimator.
  - What actually differs (the response campaign's targets):
      M1 non-stat signal amplitude sqrt(Csig(0)): real  795 vs sim 1414 um
         (the pixel-level image of the width double-count; drops at v6.1)
      M2 correlation-length ratio C_sub(4)/C_sub(1): real 0.26 vs sim 0.68
         (real decorrelates in 2-4 rows; zero-crossing 6 vs ~14)
      M3 floor-subtracted cell coherence: real 296 vs sim 156 um
      M4 C_sub(1): real 0.76 vs sim 0.72 — ALREADY MATCHED (holdout, not target)
    Stability under 2x thinning: +-0.04 on C_sub(1), +-0.01 on M2, +-2 um on M3.
  - Caveat baked into the figure: with the twist in both, real's long-lag
    plateau is NEGATIVE (-0.31, oscillatory sawtooth), so the 3-band
    WHITE/SHORT/SMOOTH shares overflow for real — judge real by the C_sub(d)
    curve and M1-M4; the shares remain useful for sim-vs-sim comparisons.
Bands on M1-M4 will be proposed when the response campaign's request doc
opens; v6.1 acceptance here stays the unchanged battery with twist holdouts
as gates. The v61 battery runs when the true v6.1 delivery lands (manifest
md5s — NOT the stale Aug-21 draft).

---------------------------------------------------------------------------------

## 2026-08-23 — v6.1 ACCEPTED on the probe battery (tag v61); diagnostics recorded

v6.1 (width re-balance, path (ii) amplitude-only, SMOD 0.0426 -> 0.009) shipped
and sealed by the pipeline; digi md5 verified against the manifest
(d0a5df64... — the stale Aug-21 draft was overwritten by the true delivery).
Battery: twist_probe + nonrms_probe + gtail_probe + cfloor_probe, tag v61,
ledgers island_post/*_v61.txt, figures sim_validation_plots/*_v61.png.

GATES (twist holdouts) — ALL PASS:
  - split-half Dsagitta sim +6.28 (side0 +6.23 / side1 +6.34) — inside
    +6.3 +- 0.3 per side; independently reproduces the pipeline delivery run;
    charge independence held (+5.73/+6.71); real reference unchanged
    (+6.53; +6.04/+7.15).
  - twist profile RMS-diff to real 140 um (gate <= ~150; matches their 140).
  - pixel pooled q68 = 1.88 mm bit-stable.
  - family checks: clipped GLOBAL data/MC 1.03 (v6t 1.01, pre-twist 1.04);
    raw 1.52; d0 med 2.52 / q90 6.87 (family).
DIAGNOSTICS (recorded per clause-6 adjudicated language — "best reachable
configuration; the composition axis is owned by the response campaign"):
  - raw C(1) = 0.61 on this harness (v6t 0.57; predicted drift ~+0.05;
    pipeline harness 0.662). raw cell coherence 141-154 um (v6t 145-158;
    real 307-310).
  - cfloor M-set: M1 sqrt(Csig0) sim 1385 um (v6t 1414; real 795) — the SMOD
    shrink removed only ~285 um in quadrature at pixel level; M2
    C_sub(4)/C_sub(1) = 0.68 (unchanged; real 0.26) — the correlation-length
    gap is intact, as an amplitude-only change must leave it; M4 C_sub(1)
    0.79 (v6t 0.72; real 0.76) and WHITE share 0.28 -> 0.21 — removing white
    SMOD raised the correlated fraction, the direction the clause-6 pilot
    predicted; all consistent, all diagnostics.
  - sim trimmed GLOBAL 1581 -> 1542, trimmed LF 809 -> 734 (SMOD width
    removed at pixel level, small — matches its small visible row-mean share).
Cosmetic: twist_probe P1 title hardcoded "(v6 digi)" -> now versioned;
v61 figure regenerated.

---------------------------------------------------------------------------------

## 2026-08-23 — DE-TWIST option for real-data fitting, measured (available, NOT adopted)

Question (user): any updates for real-data circle/sagitta fitting since v6.1?
v6.1 itself changes nothing on the real side (sim-only; all real fit numbers
bit-identical across v6/v6t/v61 ledgers). But the twist measurement created one
new OPTION, now quantified: correct real pixels by real's own measured profile
before fitting (subtract twist_probe's per-(side,row) mean D(rphi) from the
pixel azimuths). Producer: twist_probe.C "REAL DE-TWISTED" block (ledger tag
detw). Result: split-half Dsagitta +6.53 -> +1.85 mm (side0 +1.57 / side1
+2.09); residual mean profile RMS 431 -> 108 um. The remaining +1.85 mm =
per-sector scatter (~100 um) + higher-order structure the phi-averaged mean
profile cannot carry (a per-sector table and/or a second iteration would go
lower — measurable on request).
STATUS: available, NOT adopted. Adoption notes: (i) for sim-vs-real symmetric
meters it is unnecessary — both sides carry the twist since V6t; (ii) for
REAL-only physics fitting (ms_real / nf_tracks / exhaustive-finder curvature,
d0, sagitta studies) it removes a known +6.5 mm signed bias and is a one-line
correction; (iii) if adopted anywhere, adopt it as an explicit flag with the
profile file pinned (twist_profile_v6.txt), never silently — real-side numbers
change wherever it is on.

---------------------------------------------------------------------------------

## 2026-08-24 — de-twisted real in the nf_digipix meters: med RMS does NOT halve
## (GLOBAL -1.5%, LOCAL -0.1%) — and cannot, by the variance budget

Question (user): de-twist the real data, plot in the nf_digipix style, check
whether the med RMS drops by more than half. Producer:
island_post/detwist_digipix.C -> detwist_digipix_v61.txt +
sim_validation_plots/detwist_digipix_v61.png (nf_digipix-style two panels,
three curves: real / real DE-TWISTED / sim digi v6.1; groups, gates, fitter
verbatim nf_digipix; de-twist = subtract real's own measured per-(side,row)
profile from the pixel azimuths, measured in-macro exactly as twist_probe.C).

RESULT: GLOBAL med 2358 -> 2323 um (-1.5%), LOCAL 1205 -> 1204 (-0.1%);
data/MC 1.52 -> 1.50 and 0.90 -> 0.90. NOT more than half — and >50% is
impossible in principle: halving GLOBAL (2358 -> 1179) would land below the
1205 um charge-cloud floor, and the LOCAL panel IS that floor. The variance
budget says why: the twist's share of the whole-track RMS is
sqrt(2358^2 - 2323^2) ~ 405 um in quadrature (the previously-attributed twist
LF share), while the RMS medians are owned by road pickup (~1.7 mm in
quadrature, removable only by the 3-sigma clip) and the cloud width. Where
the de-twist IS large is the SIGNED meter: split-half Dsagitta +6.53 ->
+1.85 mm (previous entry) — a mm-scale bias can be a small RMS share; RMS
adds in quadrature, bias adds linearly in its own meter.
