# sPHENIX published results — local mirror (fetched 2026-08-08)

Public collaboration results only (papers from arXiv, conference notes from
sphenix.bnl.gov). No internal data. No HEPData records existed for these
at fetch time (hepdata.net: none found / bot-walled); no machine-readable
tables are published — figure PDFs + text only. Extracted text: *.txt.

## Files
- sPHENIX_2504.02240.pdf  sPH-BULK-2025-01  charged-hadron multiplicity, Au+Au 200 GeV (JHEP 08 (2025) 075)
- sPHENIX_2504.02242.pdf  sPH-BULK-2025-02  dE_T/deta, Au+Au 200 GeV (PRC)
- sPHENIX_2606.17184.pdf  sPH-JET-2026-01   dijet imbalance/acoplanarity, p+p 200 GeV, 41 pb-1  [2024 run = OUR period]
- sPHENIX_2607.03892.pdf  sPH-JET-2026-02   isolated prompt photon, p+p 200 GeV, 64.4 pb-1      [2024 run]
- sPH-CONF-JET-2025-03.pdf  inclusive jet cross-section, p+p 200 GeV (conf note)
- sPH-CONF-JET-2026-03.pdf  underlying event vs jet production, p+p 200 GeV (conf note)

## Numbers extracted for run-79507 pipeline cross-checks

1. **sigma_MBD (Vernier scan, 2024 RHIC pp period) = 25.6 mb ~ 61% of
   sigma_inel = 42 mb** (sPH-CONF-JET-2025-03; also the luminosity basis of
   sPH-JET-2026-02, L = 64.4 +5.9/-4.3 pb-1).
   -> OFFICIAL eps_MBD ~ 0.61. Our scan: eps = 0.588 +- 0.033. CONSISTENT
   (0.61 inside our 1-sigma band). The open "team eps_MBD" question is
   answered at collaboration level; rspec at 0.588 stands (0.61 would be a
   +3.7% content rescale = within our px closure resolution).
2. **pp MB trigger definition** (sPH-JET-2026-01): coincidence — at least
   one MBD tube above threshold on EACH side. Matches our gen proxy
   (mbdN>0 && mbdS>0) exactly. (The UE conf note uses an either-side
   definition for its MB sample — definition varies per analysis.)
3. **z-vertex width, 2024 pp: sigma_z ~ 50 cm Gaussian** (jet cross-section
   note: sim uses Gaussian sigma 50 cm "to match the z-vertex distribution
   in data"; photon paper: |z|<60 cm retains ~76% of events -> Gaussian
   sigma ~ 50 cm, consistent). SUPERSEDES the 16 cm nominal from the
   macros-repo beam params that our vz probe used. NOTE the standing
   tension: our CRN pilot showed vertex spread smears the drift-edge step
   (-4.8% per 16 cm) while the real step is sharp — at 50 cm the naive
   effect would be larger still; resolution must involve envelope refit
   absorption and/or z-window bookkeeping. OPEN modeling question,
   replaces the "ask team for sigma_z" item.
4. **sqrt(s) = 200 GeV** confirmed in every 2024-pp publication.
5. NOT available publicly: any event-level / hit-level data; no direct
   raw cross-check of run 79507 is possible from published material.
6. Follow-up material: sPH-CONF-JET-2026-03 (underlying-event densities vs
   Pythia-8 in pp 200) — candidate external test of the soft-content tune
   direction; sPH-CONF-JET-2026-02 (O+O jets) — species lever if OO data
   ever reaches us; Au+Au papers = generator-confrontation genre data
   (user's independent-paper idea), species-mismatched for this pipeline.

## World-data confrontation of the generator (gen_dnde, 2026-08-08)

dNch/deta(|eta|<0.5), INELASTIC pp 200, all charged finals (full decays),
6000 ev, vs PHOBOS inelastic 2.29 +- 0.08:
  ours (pT0Ref 1.85):   3.497   (+53%)
  Monash default 2.28:  2.723   (+19%)  [Monash known-high at RHIC]
  official MDC2 tune:   1.952   (-15%)
Definition caveat: strange-feed-down conventions differ (~10% scale); even
so +53% is far outside it. IMPLICATION: our +21% soft tune OVERSHOOTS
world multiplicity — pT0Ref=1.85 is an EFFECTIVE knob compensating
something else (spectrum shape: TPC charge is curler/track-length
weighted, so a softer spectrum at world-consistent multiplicity could
close the same TPC content; or residual response under-charge). v5.3b
remains valid as an effective detector-level calibration; its PHYSICS
interpretation is downgraded and the shape-vs-scale question is a
declared open item (binary: gen_dnde in P5/angantyr).

## Machine-readable tables recovered from arXiv LaTeX sources (2026-08-08)

No ancillary data files exist in any of the four tarballs (checked), and
no HEPData/CSV/TXT is published anywhere. But the papers' own result
tables were extracted from the .tex sources into CSV:
- sPH-BULK-2025-01_dNchdeta_AuAu200_table.csv — THE central result of
  2504.02240: dNch/deta||eta|<0.3 vs centrality (15 classes, 0-3% =
  723.4+-45.3) + <Npart> + per-participant-pair yields. This is the
  generator-confrontation target table (user's paper idea).
- sPH-BULK-2025-02_Npart_glauber_table.csv — Glauber <Npart> per
  centrality from 2504.02242 (its dE_T values are figure-only, no
  source table).
- The pp papers (2606.17184, 2607.03892) publish results as figures only;
  their key scalars (sigma_MBD 25.6 mb, L, vertex cuts) are in the
  extracted text above.
arxiv_src/ keeps the full LaTeX sources for provenance.

## Figure-data extraction from vector PDFs (2026-08-08) — VALIDATED

ROOT figure PDFs are vector graphics: marker draw-coordinates are in the
content streams, so extraction is exact digitization, not eyeballing.
- sPH-BULK-2025-01_dNchdeta_vs_eta_EXTRACTED.csv — the FULL dNch/deta(eta)
  distributions (121 points: 11 centralities x 11 eta bins, |eta|<1.05,
  both analysis methods per point) from
  arxiv_src/2504.02240/fig/Results/dNdEta_RHIC_woBRAHMS_altCent.pdf.
  Figure layout decoded: two pads (6+5 centralities), 2 method markers
  per point, log-y. One calibration constant (tick-label baseline offset)
  fixed on the paper's own mid-rapidity table; VALIDATION on the 6 table
  classes: ratios 0.996-1.005, scatter 0.33% -> sub-percent extraction.
  This dataset exists in machine-readable form NOWHERE else.
Same technique applies on demand to: the other Results figures (per-Npart
scaling, RHIC-combined), the pp papers' per-figure PDFs on the collab
site (dijet xJ/acoplanarity distributions, photon cross-section, MB
trigger-efficiency curve in the jet cross-section note).
Extractor: pdfminer.six (installed --user); parser inline in session logs.

## Full figure-extraction campaign (2026-08-08, "all available")

VALIDATED sub-percent method (see multiplicity benchmark above; Fig6
x-positions independently recover the Glauber Npart table to <=0.6).
- sPH-BULK-2025-01_dNchdeta_vs_eta_EXTRACTED.csv (121 pts, benchmark 0.33%)
- sPH-BULK-2025-02_dETdeta_perNpart_Fig6_EXTRACTED.csv — dET/deta/(0.5Npart)
  vs Npart, 8 centralities (0-5%: 3.792 -> dET/deta ~ 663 GeV). Fig-only data.
- sPH-BULK-2025-02_dETdeta_0-5pct_Fig5_EXTRACTED.csv — EMCal-only + HCal-only
  dET/deta(eta) at 0-5% (2 series x 25 eta bins).
- sPH-JET-2026-01_xJ_R04_range{0,1,2}_EXTRACTED.csv — unfolded dijet xJ
  distributions, R=0.4, three leading-jet pT ranges (19 data pts each;
  PYTHIA curves are line-drawn, reproducible, not extracted).
SKIPPED (documented, extractable on demand with the same method):
photon figure8/xtscaling (multi-experiment log compilation, no validation
anchor), dphi acoplanarity set, R=0.7 duplicates, remaining multiplicity
Results figs (content duplicates table + other-experiment overlays).
QA plots (tracklet/centrality checks) not result data - skipped.

## Retune attempt with these values (2026-08-08): NO-GO PROVEN
See PIPELINE.md "2-KNOB WORLD-CONSTRAINED RETUNE SCAN": at world dNch/deta
= 2.29, the TPC-charge proxy is 0.87-0.90 regardless of fragmentation
softening (required: 1.21). The content excess is NOT generator-ownable;
pT0Ref 1.85 stays as declared effective compensation. gen_scan2 kept.
