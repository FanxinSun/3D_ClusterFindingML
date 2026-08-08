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
