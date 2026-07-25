// attribution_pp.C — pp re-run of the pAu-era ancestry attribution (2026-07-25):
// TPC deposited energy classified by the PRIMARY ancestor's |eta| and pT
// (gancestor -> ntp_prim), plus the direct-vs-secondary split (gprimary flag
// on the DEPOSITING track). Purpose: split the -22% content deficit's
// ownership bound — primaries outside |eta|~1.1 can only contribute via
// MATERIAL interactions (splash-back / albedo, i.e. G4-modeled), while
// in-acceptance soft primaries are the generator's soft sector.
// pAu-era reference (10 HIJING pau200 events, energy shares):
//   |eta|: <0.5 19.3 | 0.5-1.1 26.3 | 1.1-2.0 26.4 | >2.0 28.0 (%)
//   pT:    <0.2 25.3 | 0.2-0.5 54.1 | 0.5-1.0 14.7 | >1.0 5.8  (%)
#include <TFile.h>
#include <TTree.h>
#include <cmath>
#include <cstdio>
#include <unordered_map>

void attribution_pp()
{
  double eeta[4] = {0, 0, 0, 0};   // |eta| <0.5, 0.5-1.1, 1.1-2.0, >2.0
  double ept[4] = {0, 0, 0, 0};    // pT <0.2, 0.2-0.5, 0.5-1.0, >1.0
  double edir = 0, esec = 0, eorph = 0, etot = 0;
  double esec_fwd = 0, esec_cen = 0;  // secondary-deposited, by primary |eta| vs 1.1
  long orph = 0, nhit = 0;
  for (int i = 0; i < 10; ++i)
  {
    TFile *f = TFile::Open(Form("/home/rog/sPHENIX/3D_ClusterFindingML/P5/PP_g4hit_%d.root", i));
    TTree *p = (TTree *) f->Get("ntp_prim");
    float pev, ptrk, ppt, peta;
    p->SetBranchStatus("*", 0);
    for (auto b : {"event", "ptrk", "ppt", "peta"}) p->SetBranchStatus(b, 1);
    p->SetBranchAddress("event", &pev);
    p->SetBranchAddress("ptrk", &ptrk);
    p->SetBranchAddress("ppt", &ppt);
    p->SetBranchAddress("peta", &peta);
    std::unordered_map<uint64_t, std::pair<float, float>> prim;  // (ev<<32|trk) -> (pt, eta)
    prim.reserve(p->GetEntries() * 2);
    for (Long64_t k = 0; k < p->GetEntries(); ++k)
    {
      p->GetEntry(k);
      prim[((uint64_t) (uint32_t) (int) pev << 32U) | (uint32_t) (int) ptrk] = {ppt, peta};
    }
    TTree *g = (TTree *) f->Get("ntp_g4hit");
    float gev, ged, ganc, gpri;
    g->SetBranchStatus("*", 0);
    for (auto b : {"event", "gedep", "gancestor", "gprimary"}) g->SetBranchStatus(b, 1);
    g->SetBranchAddress("event", &gev);
    g->SetBranchAddress("gedep", &ged);
    g->SetBranchAddress("gancestor", &ganc);
    g->SetBranchAddress("gprimary", &gpri);
    for (Long64_t k = 0; k < g->GetEntries(); ++k)
    {
      g->GetEntry(k);
      nhit++;
      etot += ged;
      auto it = prim.find(((uint64_t) (uint32_t) (int) gev << 32U) | (uint32_t) (int) ganc);
      if (it == prim.end()) { orph++; eorph += ged; continue; }
      double pt = it->second.first, aeta = std::fabs(it->second.second);
      int ie = aeta < 0.5 ? 0 : (aeta < 1.1 ? 1 : (aeta < 2.0 ? 2 : 3));
      int ip = pt < 0.2 ? 0 : (pt < 0.5 ? 1 : (pt < 1.0 ? 2 : 3));
      eeta[ie] += ged;
      ept[ip] += ged;
      if (gpri > 0.5) edir += ged;
      else { esec += ged; (aeta < 1.1 ? esec_cen : esec_fwd) += ged; }
    }
    f->Close();
    printf("ATTR file %d done\n", i);
  }
  double ea = etot - eorph;
  printf("ATTRPP hits %ld | orphan %.2f%% of energy\n", nhit, 100. * eorph / etot);
  printf("ATTRPP |eta|:  <0.5 %.1f%% | 0.5-1.1 %.1f%% | 1.1-2.0 %.1f%% | >2.0 %.1f%%   (pAu: 19.3/26.3/26.4/28.0)\n",
         100 * eeta[0] / ea, 100 * eeta[1] / ea, 100 * eeta[2] / ea, 100 * eeta[3] / ea);
  printf("ATTRPP pT:     <0.2 %.1f%% | 0.2-0.5 %.1f%% | 0.5-1.0 %.1f%% | >1.0 %.1f%%   (pAu: 25.3/54.1/14.7/5.8)\n",
         100 * ept[0] / ea, 100 * ept[1] / ea, 100 * ept[2] / ea, 100 * ept[3] / ea);
  printf("ATTRPP deposit: DIRECT-by-primary %.1f%% | via SECONDARIES %.1f%% "
         "(of which primary in-acceptance %.1f%%, primary forward %.1f%%)\n",
         100 * edir / ea, 100 * esec / ea, 100 * esec_cen / ea, 100 * esec_fwd / ea);
  printf("ATTRPP SPLASH-BACK BOUND (forward-primary share, G4-material-mediated): %.1f%%\n",
         100 * (eeta[2] + eeta[3]) / ea);
}
