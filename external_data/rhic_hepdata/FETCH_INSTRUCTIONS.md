# HEPData .root bundles — manual fetch needed (Cloudflare blocks automation)

sPHENIX itself has NO public .root data (no HEPData uploads yet, no open-data
program, event-level files are collaboration-internal). Legacy-RHIC reference
records DO offer .root — but hepdata.net's bot protection blocks scripted
download from this machine. Fetching them is 3 browser clicks; drop the
.zip files in THIS directory and tell the sim-pipeline session — it will
unpack, organize, and wire the cross-checks.

1. PHOBOS multiplicity compilation (pp/dAu/CuCu/AuAu, incl. inelastic pp 200
   dNch/deta — the generator-confrontation reference):
   https://www.hepdata.net/record/ins876609    -> "Download All" -> ROOT
2. STAR identified spectra pp 200 (pi/K/p soft spectra shapes):
   https://www.hepdata.net/record/ins793126    -> ROOT
3. PHENIX pi+- cross sections pp 200:
   https://www.hepdata.net/record/ins1315330   -> ROOT

Key reference value used meanwhile (from literature, no download needed):
  inelastic pp 200 GeV: dNch/deta||eta|~0 = 2.29 +- 0.08 (PHOBOS)
