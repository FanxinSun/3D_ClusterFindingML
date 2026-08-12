// ms_nofinder.C — the NO-exhaustive-finder branch of the MS / circle-fit
// comparison (user, 2026-08-10): each level fitted under its NATIVE grouping,
// sim vs real on the SAME canvas.
//   nf_hits     : sim ntp_g4hit (grouped by truth gtrackID, per collision)
//                 vs real ntp_hit pixels (grouped by the tracker: pixels
//                 road-matched to ntp_clus_trk seed clusters, same layer,
//                 dxy < 1.2 cm, |dtbin| <= 6). One canvas, log-x RMS.
//   nf_clusters : sim island91 ntp_cluster (grouped by row-aligned ntp_truth
//                 gtrackID) vs real ntp_clus_trk clusters (grouped by seedID).
//                 One canvas, linear RMS.
//   nf_tracks   : real ntp_clus_trk only (sim has no tracker output — skip):
//                 event display with fitted circles + RMS/stat panel.
// UNIFORM BAR everywhere (the only selection): >= 12 points, radial span
// >= 15 cm, 45 <= R_fit < 2e4 cm. RAW fit — no outlier cleaning, no pT
// windows. Fitter = Kasa + 6 Gauss-Newton, copied VERBATIM from ms_real.C
// (MSR::fitCircle) so every level shares one estimator.
// Real-pixel selection: CANON layer 7-54 && adc>0. The laser-flash window
// (tbin 322-340) is NOT excluded: the association road is tbin-limited
// around genuine seed clusters, so flash contamination is negligible, and
// the sim side (truth hits, no time axis) admits no symmetric cut.
// Outputs: ../sim_validation_plots/ms_nofinder_{hits,clusters,tracks_real}_<ver>.png
//          ms_nofinder_<ver>.txt (shared ledger, appended per entry)
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TEllipse.h>
#include <TPad.h>
#include <TROOT.h>
#include <cstdio>
#include <cstdarg>
#include <cmath>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>
#include <algorithm>

namespace MNF
{
struct Fit { double a = 0, b = 0, R = 0, rms = 0; int n = 0; bool ok = false; };

Fit fitCircle(const std::vector<double> &X, const std::vector<double> &Y)
{
  Fit F; F.n = (int) X.size();
  if (F.n < 5) return F;
  double Sx = 0, Sy = 0, Sxx = 0, Syy = 0, Sxy = 0, Sxz = 0, Syz = 0, Sz = 0;
  for (size_t i = 0; i < X.size(); ++i)
  {
    double x = X[i], y = Y[i], z = x * x + y * y;
    Sx += x; Sy += y; Sxx += x * x; Syy += y * y; Sxy += x * y;
    Sxz += x * z; Syz += y * z; Sz += z;
  }
  double n = F.n;
  double det = Sxx * (Syy * n - Sy * Sy) - Sxy * (Sxy * n - Sy * Sx) + Sx * (Sxy * Sy - Syy * Sx);
  if (std::fabs(det) < 1e-9) return F;
  double A = (Sxz * (Syy * n - Sy * Sy) - Sxy * (Syz * n - Sy * Sz) + Sx * (Syz * Sy - Syy * Sz)) / det;
  double B = (Sxx * (Syz * n - Sy * Sz) - Sxz * (Sxy * n - Sy * Sx) + Sx * (Sxy * Sz - Syz * Sx)) / det;
  double C = (Sxx * (Syy * Sz - Syz * Sy) - Sxy * (Sxy * Sz - Syz * Sx) + Sxz * (Sxy * Sy - Syy * Sx)) / det;
  F.a = A / 2; F.b = B / 2;
  double r2 = C + F.a * F.a + F.b * F.b;
  if (r2 <= 0) return F;
  F.R = std::sqrt(r2);
  for (int it = 0; it < 6; ++it)
  {
    double M[3][3] = {{0}}, v[3] = {0};
    for (size_t i = 0; i < X.size(); ++i)
    {
      double dx = X[i] - F.a, dy = Y[i] - F.b, rho = std::hypot(dx, dy);
      if (rho < 1e-9) continue;
      double res = rho - F.R, J[3] = {-dx / rho, -dy / rho, -1.};
      for (int p = 0; p < 3; ++p)
      {
        v[p] -= J[p] * res;
        for (int q = 0; q < 3; ++q) M[p][q] += J[p] * J[q];
      }
    }
    double d = M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
             - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
             + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
    if (std::fabs(d) < 1e-12) break;
    double d0 = (v[0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
               - M[0][1] * (v[1] * M[2][2] - M[1][2] * v[2])
               + M[0][2] * (v[1] * M[2][1] - M[1][1] * v[2])) / d;
    double d1 = (M[0][0] * (v[1] * M[2][2] - M[1][2] * v[2])
               - v[0] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
               + M[0][2] * (M[1][0] * v[2] - v[1] * M[2][0])) / d;
    double d2 = (M[0][0] * (M[1][1] * v[2] - v[1] * M[2][1])
               - M[0][1] * (M[1][0] * v[2] - v[1] * M[2][0])
               + v[0] * (M[1][0] * M[2][1] - M[1][1] * M[2][0])) / d;
    F.a += d0; F.b += d1; F.R += d2;
  }
  double s2 = 0;
  for (size_t i = 0; i < X.size(); ++i)
  {
    double res = std::hypot(X[i] - F.a, Y[i] - F.b) - F.R;
    s2 += res * res;
  }
  F.rms = std::sqrt(s2 / n);
  F.ok = true;
  return F;
}

double med(std::vector<double> v)
{
  if (v.empty()) return 0;
  std::sort(v.begin(), v.end());
  return v[v.size() / 2];
}

// the ONE selection of this file
struct Grp { std::vector<double> x, y, r; };
bool fitBar(const Grp &G, Fit &F)
{
  if ((int) G.x.size() < 12) return false;
  double rlo = 1e9, rhi = 0;
  for (double r : G.r) { rlo = std::min(rlo, r); rhi = std::max(rhi, r); }
  if (rhi - rlo < 15) return false;
  F = fitCircle(G.x, G.y);              // F is reference, so the fitted results stored;
  return F.ok && F.R >= 45 && F.R < 2e4;
}

FILE *openLedger(const char *ver)
{
  return fopen(Form("ms_nofinder_%s.txt", ver), "a");
}
}  // namespace MNF

// ---------------------------------------------------------------------------
// 1. HIT LEVEL, no finder: sim truth-grouped ntp_g4hit vs real tracker-grouped
//    ntp_hit pixels, same canvas (log-x RMS: the two are 40x apart).
void nf_hits(const char *g4pat = "../P5/PP_g4hit_%d.root", int nfiles = 10,
             const char *realf = "../clusters_seeds_island_79507-0.root_ntuplizer.root",
             const char *ver = "v54c")
{
  using namespace MNF;
  // --- sim: truth groups per collision -------------------------------------
  std::vector<double> srms, sR;
  long scoll = 0, sgrp = 0;
  for (int fi = 0; fi < nfiles; ++fi)
  {
    TFile *f = TFile::Open(Form(g4pat, fi));
    if (!f || f->IsZombie()) continue;
    TTree *t = (TTree *) f->Get("ntp_g4hit");
    float ev, gx, gy, tid;
    t->SetBranchStatus("*", 0);
    for (auto b : {"event", "gx", "gy", "gtrackID"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("event", &ev);
    t->SetBranchAddress("gx", &gx);
    t->SetBranchAddress("gy", &gy);
    t->SetBranchAddress("gtrackID", &tid);
    std::map<std::pair<int, int>, Grp> g;
    std::set<int> evs;
    for (Long64_t i = 0; i < t->GetEntries(); ++i)
    {
      t->GetEntry(i);
      double r = std::hypot(gx, gy);
      if (r < 20 || r > 78) continue;
      Grp &G = g[{(int) ev, (int) tid}];
      G.x.push_back(gx); G.y.push_back(gy); G.r.push_back(r);
      evs.insert((int) ev);
    }
    f->Close();
    scoll += (long) evs.size();
    for (auto &kv : g)
    {
      Fit F;
      if (!fitBar(kv.second, F)) continue;
      sgrp++;
      srms.push_back(F.rms * 1e4);          // um
      sR.push_back(F.R);
    }
  }
  // --- real: pixels road-matched to the tracker's seed clusters ------------
  struct SC { float x, y, tb; int seed; };
  std::vector<SC> sc;
  std::unordered_map<int, std::vector<int>> bucket;   // ev*100+layer -> sc idx
  std::map<std::pair<int, int>, int> seedIdx;
  {
    TFile *f = TFile::Open(realf);
    TTree *t = (TTree *) f->Get("ntp_clus_trk");
    float ev, sid, lay, x, y, tb;
    t->SetBranchStatus("*", 0);
    for (auto b : {"event", "seedID", "layer", "x", "y", "tbin"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("event", &ev);
    t->SetBranchAddress("seedID", &sid);
    t->SetBranchAddress("layer", &lay);
    t->SetBranchAddress("x", &x);
    t->SetBranchAddress("y", &y);
    t->SetBranchAddress("tbin", &tb);
    for (Long64_t i = 0; i < t->GetEntries(); ++i)
    {
      t->GetEntry(i);
      if (lay < 7 || lay > 54) continue;
      auto key = std::make_pair((int) ev, (int) sid);
      auto it = seedIdx.find(key);
      if (it == seedIdx.end()) it = seedIdx.insert({key, (int) seedIdx.size()}).first;
      bucket[(int) ev * 100 + (int) lay].push_back((int) sc.size());
      sc.push_back({x, y, tb, it->second});
    }
    f->Close();
  }
  int nseed = (int) seedIdx.size();
  std::vector<Grp> px(nseed);
  std::vector<std::set<int>> hid(nseed);
  {
    TFile *f = TFile::Open(realf);
    TTree *t = (TTree *) f->Get("ntp_hit");
    float ev, lay, x, y, tb, adc, hitID;
    t->SetBranchStatus("*", 0);
    for (auto b : {"event", "layer", "x", "y", "tbin", "adc", "hitID"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("event", &ev);
    t->SetBranchAddress("layer", &lay);
    t->SetBranchAddress("x", &x);
    t->SetBranchAddress("y", &y);
    t->SetBranchAddress("tbin", &tb);
    t->SetBranchAddress("adc", &adc);
    t->SetBranchAddress("hitID", &hitID);
    for (Long64_t i = 0; i < t->GetEntries(); ++i)
    {
      t->GetEntry(i);
      if (lay < 7 || lay > 54 || adc <= 0) continue;
      auto bit = bucket.find((int) ev * 100 + (int) lay);
      if (bit == bucket.end()) continue;
      for (int j : bit->second)
      {
        const SC &c = sc[j];
        if (std::fabs(tb - c.tb) > 6) continue;
        double dx = x - c.x, dy = y - c.y;
        if (dx * dx + dy * dy > 1.2 * 1.2) continue;
        if (!hid[c.seed].insert((int) hitID).second) break;   // already in this seed
        px[c.seed].x.push_back(x);
        px[c.seed].y.push_back(y);
        px[c.seed].r.push_back(std::hypot(x, y));
        break;                                                // one seed per pixel
      }
    }
    f->Close();
  }
  std::vector<double> rrms, rR;
  long rgrp = 0, rpx = 0;
  for (int s = 0; s < nseed; ++s)
  {
    Fit F;
    if (!fitBar(px[s], F)) continue;
    rgrp++;
    rpx += (long) px[s].x.size();
    rrms.push_back(F.rms * 1e4);
    rR.push_back(F.R);
  }
  // --- ledger + canvas ------------------------------------------------------
  FILE *fo = openLedger(ver);
  auto P = [&](const char *fmt, ...) {
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    va_start(ap, fmt); vfprintf(fo, fmt, ap); va_end(ap);
  };
  P("[nf_hits %s] HIT level, NO finder; bar >=12 pts, span >=15 cm, R_fit >=45 cm, raw fit\n", ver);
  P("  sim ntp_g4hit truth-grouped: %ld collisions -> %ld tracks | RMS med %.1f um | R med %.0f cm\n",
    scoll, sgrp, med(srms), med(sR));
  P("  real ntp_hit tracker-grouped (road dxy<1.2 cm, |dtbin|<=6): %d seeds -> %ld tracks "
    "(%.0f px/track mean) | RMS med %.0f um | R med %.0f cm\n",
    nseed, rgrp, rgrp ? (double) rpx / rgrp : 0, med(rrms), med(rR));
  P("  gap = the detector: diffusion + pad/tbin quantization (no centroiding at pixel level)\n");
  fclose(fo);
  gStyle->SetOptStat(0);
  double e[53];
  for (int i = 0; i <= 52; ++i) e[i] = 3. * std::pow(10., 3.523 * i / 52.);   // 3 um .. 10 mm
  TH1D *hs = new TH1D("nfh_s", ";per-track circle-fit RMS [#mum];tracks (unit area)", 52, e);
  TH1D *hr = new TH1D("nfh_r", ";per-track circle-fit RMS [#mum];tracks (unit area)", 52, e);
  for (double v : srms) hs->Fill(std::min(v, 9990.));
  for (double v : rrms) hr->Fill(std::min(v, 9990.));
  for (auto h : {hs, hr}) if (h->Integral() > 0) h->Scale(1. / h->Integral());
  TCanvas *cv = new TCanvas("cvnfh", "nf hits", 1100, 660);
  cv->SetLogx();
  hr->SetLineColor(kBlack); hr->SetLineWidth(2);
  hs->SetLineColor(kBlue + 1); hs->SetLineWidth(2); hs->SetLineStyle(2);
  hr->SetTitle("hit level, no finder: circle-fit quality under native grouping");
  hr->SetMaximum(1.35 * std::max(hr->GetMaximum(), hs->GetMaximum()));
  hr->Draw("hist");
  hs->Draw("hist same");
  TLegend *lg = new TLegend(0.14, 0.72, 0.60, 0.87);
  lg->SetBorderSize(0);
  lg->AddEntry(hr, Form("real ntp_hit, tracker-grouped: med %.0f #mum", med(rrms)), "l");
  lg->AddEntry(hs, Form("sim ntp_g4hit, truth-grouped: med %.1f #mum", med(srms)), "l");
  lg->Draw();
  TLatex tx; tx.SetNDC(); tx.SetTextSize(0.030);
  tx.DrawLatex(0.14, 0.66, Form("N tracks: real %ld, sim %ld", rgrp, sgrp));
  tx.DrawLatex(0.14, 0.61, "same fitter, same bar; no outlier removal, no p_{T} windows");
  tx.DrawLatex(0.14, 0.56, "the gap IS the detector: diffusion + pad/tbin quantization");
  cv->SaveAs(Form("../sim_validation_plots/ms_nofinder_hits_%s.png", ver));
  printf("wrote ms_nofinder_hits_%s.png\n", ver);
}

// ---------------------------------------------------------------------------
// 2. CLUSTER LEVEL, no finder: sim island91 ntp_cluster truth-grouped vs real
//    ntp_clus_trk tracker-grouped, same canvas.
void nf_clusters(const char *i91 = "island91_frames_production_v54c.root",
                 const char *realf = "../clusters_seeds_island_79507-0.root_ntuplizer.root",
                 const char *ver = "v54c")
{
  using namespace MNF;
  std::vector<double> rms[2], Rv[2];                  // [0]=real [1]=sim
  long ngrp[2] = {0, 0};
  {
    TFile *f = TFile::Open(realf);
    TTree *t = (TTree *) f->Get("ntp_clus_trk");
    float ev, sid, lay, x, y;
    t->SetBranchStatus("*", 0);
    for (auto b : {"event", "seedID", "layer", "x", "y"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("event", &ev);
    t->SetBranchAddress("seedID", &sid);
    t->SetBranchAddress("layer", &lay);
    t->SetBranchAddress("x", &x);
    t->SetBranchAddress("y", &y);
    std::map<std::pair<int, int>, Grp> g;
    for (Long64_t i = 0; i < t->GetEntries(); ++i)
    {
      t->GetEntry(i);
      if (lay < 7 || lay > 54) continue;
      Grp &G = g[{(int) ev, (int) sid}];
      G.x.push_back(x); G.y.push_back(y); G.r.push_back(std::hypot(x, y));
    }
    f->Close();
    for (auto &kv : g)
    {
      Fit F;
      if (!fitBar(kv.second, F)) continue;
      ngrp[0]++; rms[0].push_back(F.rms * 10); Rv[0].push_back(F.R);
    }
  }
  {
    TFile *f = TFile::Open(i91);
    TTree *c = (TTree *) f->Get("ntp_cluster");
    TTree *u = (TTree *) f->Get("ntp_truth");
    float ev, lay, x, y, tid;
    c->SetBranchStatus("*", 0);
    for (auto b : {"event", "layer", "x", "y"}) c->SetBranchStatus(b, 1);
    c->SetBranchAddress("event", &ev);
    c->SetBranchAddress("layer", &lay);
    c->SetBranchAddress("x", &x);
    c->SetBranchAddress("y", &y);
    u->SetBranchStatus("*", 0);
    u->SetBranchStatus("gtrackID", 1);
    u->SetBranchAddress("gtrackID", &tid);
    std::map<std::pair<int, int>, Grp> g;
    for (Long64_t i = 0; i < c->GetEntries(); ++i)
    {
      c->GetEntry(i); u->GetEntry(i);
      if (lay < 7 || lay > 54 || tid <= 0) continue;
      Grp &G = g[{(int) ev, (int) tid}];
      G.x.push_back(x); G.y.push_back(y); G.r.push_back(std::hypot(x, y));
    }
    f->Close();
    for (auto &kv : g)
    {
      Fit F;
      if (!fitBar(kv.second, F)) continue;
      ngrp[1]++; rms[1].push_back(F.rms * 10); Rv[1].push_back(F.R);
    }
  }
  FILE *fo = openLedger(ver);
  auto P = [&](const char *fmt, ...) {
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    va_start(ap, fmt); vfprintf(fo, fmt, ap); va_end(ap);
  };
  P("[nf_clusters %s] CLUSTER level, NO finder; same bar, raw fit\n", ver);
  P("  real ntp_clus_trk (event,seedID): %ld tracks | RMS med %.0f um | R med %.0f cm\n",
    ngrp[0], med(rms[0]) * 1000, med(Rv[0]));
  P("  sim island91 ntp_cluster (event,gtrackID>0): %ld tracks | RMS med %.0f um | R med %.0f cm\n",
    ngrp[1], med(rms[1]) * 1000, med(Rv[1]));
  P("  data/MC RMS %.2f\n", med(rms[1]) > 0 ? med(rms[0]) / med(rms[1]) : 0);
  fclose(fo);
  gStyle->SetOptStat(0);
  TH1D *h[2];
  const char *hn[2] = {"nfc_r", "nfc_s"};
  for (int s = 0; s < 2; ++s)
  {
    h[s] = new TH1D(hn[s], ";per-track circle-fit RMS [mm];tracks (unit area)", 50, 0, 2.5);
    for (double v : rms[s]) h[s]->Fill(std::min(v, 2.49));
    if (h[s]->Integral() > 0) h[s]->Scale(1. / h[s]->Integral());
  }
  TCanvas *cv = new TCanvas("cvnfc", "nf clusters", 1100, 660);
  h[0]->SetLineColor(kBlack); h[0]->SetLineWidth(2);
  h[1]->SetLineColor(kBlue + 1); h[1]->SetLineWidth(2); h[1]->SetLineStyle(2);
  h[0]->SetTitle("cluster level, no finder: circle-fit quality under native grouping");
  h[0]->SetMaximum(1.35 * std::max(h[0]->GetMaximum(), h[1]->GetMaximum()));
  h[0]->Draw("hist");
  h[1]->Draw("hist same");
  TLegend *lg = new TLegend(0.40, 0.72, 0.89, 0.87);
  lg->SetBorderSize(0);
  lg->AddEntry(h[0], Form("real ntp_clus_trk (tracker), med %.0f #mum", med(rms[0]) * 1000), "l");
  lg->AddEntry(h[1], Form("sim v5.4c ntp_cluster (truth), med %.0f #mum", med(rms[1]) * 1000), "l");
  lg->Draw();
  TLatex tx; tx.SetNDC(); tx.SetTextSize(0.030);
  tx.DrawLatex(0.40, 0.66, Form("N tracks: real %ld, sim %ld; data/MC %.2f",
                                ngrp[0], ngrp[1], med(rms[1]) > 0 ? med(rms[0]) / med(rms[1]) : 0));
  tx.DrawLatex(0.40, 0.61, "same fitter, same bar; no outlier removal, no p_{T} windows");
  cv->SaveAs(Form("../sim_validation_plots/ms_nofinder_clusters_%s.png", ver));
  printf("wrote ms_nofinder_clusters_%s.png\n", ver);
}

// ---------------------------------------------------------------------------
// 3. TRACK LEVEL, no finder: real-only (local sim reco yields no tracks).
//    ntp_clus_trk fitted; event display + statistics on one canvas.
void nf_tracks(const char *realf = "../clusters_seeds_island_79507-0.root_ntuplizer.root",
               const char *ver = "v54c", int showev = 7)
{
  using namespace MNF;
  std::map<std::pair<int, int>, Grp> g;
  {
    TFile *f = TFile::Open(realf);
    TTree *t = (TTree *) f->Get("ntp_clus_trk");
    float ev, sid, lay, x, y;
    t->SetBranchStatus("*", 0);
    for (auto b : {"event", "seedID", "layer", "x", "y"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("event", &ev);
    t->SetBranchAddress("seedID", &sid);
    t->SetBranchAddress("layer", &lay);
    t->SetBranchAddress("x", &x);
    t->SetBranchAddress("y", &y);
    for (Long64_t i = 0; i < t->GetEntries(); ++i)
    {
      t->GetEntry(i);
      if (lay < 7 || lay > 54) continue;
      Grp &G = g[{(int) ev, (int) sid}];
      G.x.push_back(x); G.y.push_back(y); G.r.push_back(std::hypot(x, y));
    }
    f->Close();
  }
  std::vector<double> rms, Rv, d0v;
  long ntot = 0, nfit = 0;
  for (auto &kv : g)
  {
    ntot++;
    Fit F;
    if (!fitBar(kv.second, F)) continue;
    nfit++;
    rms.push_back(F.rms * 10);
    Rv.push_back(F.R);
    d0v.push_back(std::fabs(std::hypot(F.a, F.b) - F.R));
  }
  FILE *fo = openLedger(ver);
  auto P = [&](const char *fmt, ...) {
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    va_start(ap, fmt); vfprintf(fo, fmt, ap); va_end(ap);
  };
  P("[nf_tracks %s] TRACK level, NO finder; real only (sim has no tracker output)\n", ver);
  P("  ntp_clus_trk seeds: %ld total, %ld pass bar (%.0f%%) | RMS med %.0f um | "
    "R med %.0f cm | |d0| med %.2f cm\n",
    ntot, nfit, ntot ? 100. * nfit / ntot : 0, med(rms) * 1000, med(Rv), med(d0v));
  fclose(fo);
  gStyle->SetOptStat(0);
  TCanvas *cv = new TCanvas("cvnft", "nf tracks", 1500, 700);
  cv->Divide(2, 1);
  cv->cd(1);
  {
    TH1 *fr = gPad->DrawFrame(-80, -80, 80, 80,
                              Form("real event %d: tracker tracks + fitted circles;x [cm];y [cm]", showev));
    fr->GetYaxis()->SetTitleOffset(1.25);
    int cols[5] = {kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 2, kOrange + 7};
    int k = 0;
    for (auto &kv : g)
    {
      if (kv.first.first != showev) continue;
      Fit F;
      bool pass = fitBar(kv.second, F);
      TGraph *gt = new TGraph();
      for (size_t i = 0; i < kv.second.x.size(); ++i)
        gt->SetPoint(gt->GetN(), kv.second.x[i], kv.second.y[i]);
      gt->SetMarkerStyle(pass ? 20 : 5);
      gt->SetMarkerSize(pass ? 0.45 : 0.4);
      gt->SetMarkerColor(pass ? cols[k % 5] : kGray + 1);
      gt->Draw("P same");
      if (pass && kv.second.x.size() >= 20)
      {
        // arc over the measured span only (full circles spirograph the pad)
        double sc = 0, ss = 0;
        for (size_t i = 0; i < kv.second.x.size(); ++i)
        {
          double th = std::atan2(kv.second.y[i] - F.b, kv.second.x[i] - F.a);
          sc += std::cos(th); ss += std::sin(th);
        }
        double th0 = std::atan2(ss, sc), lo = 0, hi = 0;
        for (size_t i = 0; i < kv.second.x.size(); ++i)
        {
          double d = std::atan2(kv.second.y[i] - F.b, kv.second.x[i] - F.a) - th0;
          while (d > M_PI) d -= 2 * M_PI;
          while (d < -M_PI) d += 2 * M_PI;
          lo = std::min(lo, d); hi = std::max(hi, d);
        }
        double r2d = 180. / M_PI;
        TEllipse *el = new TEllipse(F.a, F.b, F.R, F.R,
                                    (th0 + lo) * r2d - 4, (th0 + hi) * r2d + 4);
        el->SetFillStyle(0);
        el->SetNoEdges();
        el->SetLineColorAlpha(cols[k % 5], 0.55);
        el->SetLineWidth(1);
        el->Draw("same");
      }
      if (pass) k++;
    }
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.032);
    tx.DrawLatex(0.13, 0.86, Form("colored = pass bar (%d); gray = below bar", k));
    tx.DrawLatex(0.13, 0.81, "circles drawn for tracks with >= 20 clusters");
  }
  cv->cd(2);
  {
    TH1D *h = new TH1D("nft_r", ";per-track circle-fit RMS [mm];tracks", 50, 0, 2.5);
    for (double v : rms) h->Fill(std::min(v, 2.49));
    h->SetLineColor(kBlack); h->SetLineWidth(2);
    h->SetTitle("real tracker tracks: raw circle-fit RMS (all 100 events)");
    h->SetMaximum(1.35 * h->GetMaximum());
    h->Draw("hist");
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.032);
    tx.DrawLatex(0.48, 0.84, Form("seeds %ld, pass bar %ld (%.0f%%)", ntot, nfit, ntot ? 100. * nfit / ntot : 0));
    tx.DrawLatex(0.48, 0.78, Form("RMS median %.0f #mum", med(rms) * 1000));
    tx.DrawLatex(0.48, 0.72, Form("R_{fit} median %.0f cm", med(Rv)));
    tx.DrawLatex(0.48, 0.66, Form("|d0| median %.2f cm", med(d0v)));
    tx.SetTextSize(0.028);
    tx.DrawLatex(0.48, 0.58, "sim skipped: no tracker output in sim");
    tx.DrawLatex(0.48, 0.53, "(cluster-level mirror = nf_clusters)");
  }
  cv->SaveAs(Form("../sim_validation_plots/ms_nofinder_tracks_real_%s.png", ver));
  printf("wrote ms_nofinder_tracks_real_%s.png\n", ver);
}
