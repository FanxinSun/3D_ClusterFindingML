// clustwist.C — CLUSTER-level twist profile (2026-08-21): mean fit-orthogonal r-phi
// displacement D = r * wrap(phi_cluster - phi_fit(r)) [um] per (layer, side), over
// full-crosser tracks (fieldmeter gates), real ntp_clus_trk seeds vs sim island91
// truth tracks. Purpose: decide whether the pixel-level sawtooth of
// twist_profile_v6.txt is also present in the real reco CLUSTERS (-> physical,
// inject everywhere) or only in the raw pixels (-> raw-geometry / alignment class:
// inject at digi, un-twist at cluster export). Sign: + = larger azimuth (CCW).
#include "fieldmeter.C"
#include <cstdio>
#include <map>
#include <vector>
#include <cmath>
namespace CT
{
struct Grp { std::vector<double> x, y; std::vector<int> lay, side; };
double sum[55][2], cnt[55][2];
// azimuth of the fitted circle at radius r (nearest intersection to phi_pix)
bool phiFit(const FM::Fit &F, double r, double phpix, double &phf)
{
  double d = std::hypot(F.a, F.b); if (d < 1e-9) return false;
  double cosang = (r * r + d * d - F.R * F.R) / (2 * r * d);
  if (cosang < -1 || cosang > 1) return false;
  double ang = std::acos(cosang), phc = std::atan2(F.b, F.a);
  double c1 = phc + ang, c2 = phc - ang;
  auto w = [](double v) { while (v > M_PI) v -= 2 * M_PI; while (v < -M_PI) v += 2 * M_PI; return v; };
  phf = (std::fabs(w(phpix - c1)) < std::fabs(w(phpix - c2))) ? c1 : c2;
  return true;
}
void accumulate(std::map<long long, Grp> &G, const char *tag, int &ntr)
{
  for (int L = 0; L < 55; ++L) sum[L][0] = sum[L][1] = cnt[L][0] = cnt[L][1] = 0;
  ntr = 0;
  for (auto &kv : G)
  {
    Grp &g = kv.second; size_t n = g.x.size(); if (n < 12) continue;
    int lmin = 99, lmax = -1; for (int l : g.lay) { lmin = std::min(lmin, l); lmax = std::max(lmax, l); }
    if (lmin > 11 || lmax < 50) continue;
    FM::Fit F = FM::fitCircle(g.x, g.y);
    if (!F.ok || F.rms > 0.25 || F.R < 35) continue;
    ntr++;
    for (size_t i = 0; i < n; ++i)
    {
      double r = std::hypot(g.x[i], g.y[i]), ph = std::atan2(g.y[i], g.x[i]), phf;
      if (!phiFit(F, r, ph, phf)) continue;
      double dphi = ph - phf; while (dphi > M_PI) dphi -= 2 * M_PI; while (dphi < -M_PI) dphi += 2 * M_PI;
      int s = g.side[i] ? 1 : 0, L = g.lay[i];
      sum[L][s] += r * dphi * 1e4; cnt[L][s] += 1;
    }
  }
  printf("CLUSTWIST %s: %d full-crosser tracks\n", tag, ntr);
}
}  // namespace CT
void clustwist(const char *realf = "/home/rog/sPHENIX/3D_ClusterFindingML/clusters_seeds_island_79507-0.root_ntuplizer.root",
               const char *simf = "island91_frames_production_v6.root")
{
  // ---- real: ntp_clus_trk grouped by (event, seedID), laser event vetoed
  std::map<long long, CT::Grp> G;
  {
    TFile *f = TFile::Open(realf); TTree *t = (TTree *) f->Get("ntp_clus_trk");
    float ev, sid, lay, x, y, ze; t->SetBranchStatus("*", 0);
    for (auto b : {"event", "seedID", "layer", "x", "y", "zelem"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("event", &ev); t->SetBranchAddress("seedID", &sid); t->SetBranchAddress("layer", &lay);
    t->SetBranchAddress("x", &x); t->SetBranchAddress("y", &y); t->SetBranchAddress("zelem", &ze);
    for (Long64_t i = 0; i < t->GetEntries(); ++i)
    {
      t->GetEntry(i); if ((int) ev == 44 || lay < 7 || lay > 54) continue;
      CT::Grp &g = G[((long long) ev << 32) | (long long) sid];
      g.x.push_back(x); g.y.push_back(y); g.lay.push_back((int) lay); g.side.push_back((int) ze);
    }
    f->Close();
  }
  int nreal; CT::accumulate(G, "REAL clusters", nreal);
  double rs[55][2], rc[55][2]; for (int L = 0; L < 55; ++L) for (int s = 0; s < 2; ++s) { rs[L][s] = CT::sum[L][s]; rc[L][s] = CT::cnt[L][s]; }
  // ---- sim: island91 ntp_cluster/ntp_truth, cls==0 && ntrks==1, grouped by (event, gtrackID)
  G.clear();
  {
    TFile *f = TFile::Open(simf); TTree *c = (TTree *) f->Get("ntp_cluster"); TTree *u = (TTree *) f->Get("ntp_truth");
    float ev, lay, x, y, ze, gtrk, cls, ntrks;
    c->SetBranchStatus("*", 0); for (auto b : {"event", "layer", "x", "y", "zelem"}) c->SetBranchStatus(b, 1);
    c->SetBranchAddress("event", &ev); c->SetBranchAddress("layer", &lay); c->SetBranchAddress("x", &x); c->SetBranchAddress("y", &y); c->SetBranchAddress("zelem", &ze);
    u->SetBranchStatus("*", 0); for (auto b : {"gtrackID", "cls", "ntrks"}) u->SetBranchStatus(b, 1);
    u->SetBranchAddress("gtrackID", &gtrk); u->SetBranchAddress("cls", &cls); u->SetBranchAddress("ntrks", &ntrks);
    for (Long64_t i = 0; i < c->GetEntries(); ++i)
    {
      c->GetEntry(i); u->GetEntry(i);
      if (cls != 0 || ntrks != 1 || gtrk <= 0 || lay < 7 || lay > 54) continue;
      CT::Grp &g = G[((long long) ev << 32) | (long long) gtrk];
      g.x.push_back(x); g.y.push_back(y); g.lay.push_back((int) lay); g.side.push_back((int) ze);
    }
    f->Close();
  }
  int nsim; CT::accumulate(G, "SIM island91 clusters", nsim);
  // ---- pixel-level real profile from the payload for side-by-side
  double px[55][2] = {{0}};
  { FILE *fp = fopen("twist_profile_v6.txt", "r"); char l[256];
    while (fp && fgets(l, 256, fp)) { int L, s; double re, si, de; if (l[0] == '#') continue; if (sscanf(l, "%d %d %lf %lf %lf", &L, &s, &re, &si, &de) == 5) px[L][s] = re; } if (fp) fclose(fp); }
  printf("CLUSTWIST table [um]: layer | side0: realPIX realCLU simCLU | side1: realPIX realCLU simCLU   (n real clu side0/side1)\n");
  for (int L = 7; L <= 54; ++L)
    printf("CLUSTWIST %2d | %+7.0f %+7.0f %+7.0f | %+7.0f %+7.0f %+7.0f   (%.0f/%.0f)\n", L,
           px[L][0], rc[L][0] > 0 ? rs[L][0] / rc[L][0] : 0, CT::cnt[L][0] > 0 ? CT::sum[L][0] / CT::cnt[L][0] : 0,
           px[L][1], rc[L][1] > 0 ? rs[L][1] / rc[L][1] : 0, CT::cnt[L][1] > 0 ? CT::sum[L][1] / CT::cnt[L][1] : 0, rc[L][0], rc[L][1]);
  // region summaries: first/last row and the boundary jumps, real clusters
  auto m = [&](int L, int s) { return rc[L][s] > 0 ? rs[L][s] / rc[L][s] : 0.; };
  printf("CLUSTWIST REALCLU summary side0: R1 %+.0f -> %+.0f | R2 %+.0f -> %+.0f | R3 %+.0f -> %+.0f | jumps %+.0f / %+.0f um\n",
         m(7,0), m(22,0), m(23,0), m(38,0), m(39,0), m(54,0), m(23,0)-m(22,0), m(39,0)-m(38,0));
  printf("CLUSTWIST REALCLU summary side1: R1 %+.0f -> %+.0f | R2 %+.0f -> %+.0f | R3 %+.0f -> %+.0f | jumps %+.0f / %+.0f um\n",
         m(7,1), m(22,1), m(23,1), m(38,1), m(39,1), m(54,1), m(23,1)-m(22,1), m(39,1)-m(38,1));
  printf("CLUSTWIST REALPIX summary side0: R1 %+.0f -> %+.0f | R2 %+.0f -> %+.0f | R3 %+.0f -> %+.0f | jumps %+.0f / %+.0f um (payload)\n",
         px[7][0], px[22][0], px[23][0], px[38][0], px[39][0], px[54][0], px[23][0]-px[22][0], px[39][0]-px[38][0]);
}
