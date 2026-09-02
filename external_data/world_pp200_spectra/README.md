# World pp-200 pT spectra — manual download landing spot (2026-09)

Why: the spectrum-SHAPE axis of the generator scan (PIPELINE.md 2026-08-08
"WORLD-DATA CONFRONTATION": pT0Ref 1.85 overshoots world multiplicity +53%;
declared candidate owner of the charge-vs-multiplicity tension = a softer
spectrum at world multiplicity). These records are the missing input flagged
in the 2026-08-31 world-data note. Download by hand (HEPData blocks scripted
fetches from this box): open each record page -> "Download All" -> CSV, and
unpack into a folder named by the record id below.

1. PRIMARY — PHENIX, "Identified charged hadron production in p+p at
   sqrt(s)=200 and 62.4 GeV", PRC 83 (2011) 064903, arXiv:1102.0753.
   https://www.hepdata.net/record/ins886590   -> ins886590/
   pi+-, K+-, p, pbar invariant yields vs pT at |eta|<0.35 (pp 200:
   pi 0.3-3, K 0.4-2, p 0.5-4.5 GeV/c) + <pT> and dN/dy tables.
   This is the low-pT shape + species mix the TPC-charge proxy weights.

2. COMPANION — STAR, "Systematic measurements of identified particle spectra
   in pp, d+Au and Au+Au", PRC 79 (2009) 034909, arXiv:0808.2041.
   https://www.hepdata.net/record/ins793126   -> ins793126/
   Same species, independent detector (TPC dE/dx), pp 200 minbias +
   multiplicity classes; cross-checks PHENIX normalization and extends the
   shape lever arm.

3. OPTIONAL — PHENIX, "Mid-rapidity neutral pion production in pp at 200
   GeV", PRL 91 (2003) 241803, hep-ex/0304038.
   https://www.hepdata.net/record/ins617784   -> ins617784/
   pi0 cross section 1-14 GeV/c: anchors the spectrum tail; secondary for
   the TPC-charge question (curler-weighted charge is low-pT dominated).

Related (list only, fetch on need): STAR PLB 616 (2005) 8 pion/kaon/proton
pp+dAu 200 (nucl-ex/0309012); STAR high-pT identified pp 200
(arXiv:0901.0692); PHENIX 62.4 GeV inclusive charged hadrons
(arXiv:1202.4020) — energy-scaling cross-check only.

## Status (2026-09-02, after user download)

- ins886590 (PHENIX): DOWNLOADED + unpacked (HEPData-ins886590-v1.zip,
  yaml format). VERIFIED: Tables 1-4 = the pp-200 set — T1 pi+/pi-
  (0.3-3.0), T2 K+/K- (0.4-2.0), T3+T4 p/pbar (0.5-4.6 GeV/c), invariant
  cross sections E d3sig/dp3 at sqrt(s)=200, |eta|<0.35; T5-T8 = the same
  at 62.4 GeV. <pT>/dN/dy summary numbers are PRINTED TABLES in the paper
  (user's copy: /mnt/c/Users/ROG/Documents/Papers/adare2011.pdf), not in
  the HEPData record.
- ins793126 (STAR): DOWNLOADED + unpacked — but the record holds FIGURE
  data only (pp multiplicity distribution Fig03, tracking/vertex
  efficiencies Fig11-13, dAu/AuAu62/AuAu130 spectra Fig18-20, thermal-fit
  points): the pp pT SPECTRA of PRC 79 034909 are NOT in this HEPData
  record. For a STAR pp-spectra cross-check use STAR PLB 616 (2005) 8
  (nucl-ex/0309012) on HEPData, or digitize PRC 79's pp figures (user's
  copy: abelev2009.pdf). Kept anyway: the efficiency and multiplicity
  tables are useful context.
- ins617784 (PHENIX pi0): DOWNLOADED + unpacked (HEPData-ins617784-v1-csv.
  tar.gz). VERIFIED: Table 1 = pi0 invariant cross section E d3sig/dp3 vs pT,
  20 points, |eta|<0.35, sqrt(s)=200 (DOI 10.17182/hepdata.41956.v1/t1).
- ins628232 (STAR PLB 616 (2005) 8, nucl-ex/0309012): DOWNLOADED + unpacked
  (HEPData-ins628232-v1.zip; yaml + the collaboration's own .txt/.txt~ and
  PNG thumbnails, left as delivered). VERIFIED: fig2a/2b/2c = pi+-, K+-,
  p/pbar invariant yields 1/N_evt d2N/(2pi pT dpT dy) vs pT with an explicit
  "p+p NSD" series per species (pi 0.35-2.75, K 0.46-1.70, p 0.47-3.50
  GeV/c) beside d+Au minbias + 3 centrality classes; fig3 = identified
  R_dAu; fig4 = (p+pbar)/h for pp NSD and dAu minbias. This is the
  independent-detector (STAR TPC dE/dx) pp-spectra cross-check of PHENIX
  ins886590. NOTE: STAR pp yields are NSD-normalized (non-single-
  diffractive) while PHENIX is per BBC-triggered inelastic event —
  reconcile normalizations before comparing absolute yields; SHAPES
  compare directly. (Earlier guess ins628024 was wrong; HEPData's search
  box rejects slashes in arXiv ids — resolve via INSPIRE instead.)
- STAR high-pT identified pp (arXiv:0901.0692, INSPIRE 810388): NO HEPData
  record exists (404, checked 2026-09-02) — paper figures only; dropped.

Provenance rule (repo convention): keep the HEPData CSVs unmodified; any
derived/reformatted table goes beside them with a note naming the source
table + DOI. Public data only; independent of ~/Papers/n4-sphenix (which
holds none of these three).
