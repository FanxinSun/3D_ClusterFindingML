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
//   nf_ms_hits  : the ms_split.C split-arc MS (multiple scattering) tangent-
//                 mismatch test at HIT level (user, 2026-08-13): sim truth
//                 hits with the ms_split method VERBATIM (r0=35 reproduction
//                 + r0=49), real tracker-grouped pixels at r0=49 (rows start
//                 31.4 -> 35 impossible; ms_real precedent), one canvas.
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

// the ONE selection of this file. Values MEASURED as the working point
// (2026-08-20): ms_barscan_v6.txt = headline medians insensitive (n 8-20,
// span 10-25 cm, Rmin 40-60 move the med05-family ratio by <= 1%);
// ms_rocscan_v6.txt = in this given-grouping regime purity is ~0.98 at
// every variant and audit-fake seeds are rejected 100% everywhere, so
// tightening any cut only costs efficiency. R >= 45 cm (pT 0.189 GeV,
// margin above the 39.1 cm looper bound) is the one cut that also
// filters: R >= 40 admits looper leakage (-2.6 purity even here). NO gap
// or RMS gate here ON PURPOSE: gap kills stitched ghosts, which exist
// only under algorithmic grouping (that gate lives in
// missed_tracks.C::acceptTrack, where the ROC scan shows it IS the ghost
// filter); fits stay raw by the branch's design. Note the span here is
// RADIAL cm; acceptTrack's is layer units.
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

double wrapphi(double d)
{
  while (d > M_PI) d -= 2 * M_PI;
  while (d < -M_PI) d += 2 * M_PI;
  return d;
}
// tangent angle of circle (a,b,R) at its intersection with the cylinder
// r=r0, branch nearest (hx,hy), sign along (dx,dy) — verbatim ms_split.C.
bool tangentAtR(const Fit &F, double r0, double hx, double hy,
                double dx, double dy, double &psi)
{
  double d2 = F.a * F.a + F.b * F.b, d = std::sqrt(d2);
  if (d < 1e-6) return false;
  double alpha = (r0 * r0 - F.R * F.R + d2) / (2 * d);
  double h2 = r0 * r0 - alpha * alpha;
  if (h2 < 0) return false;
  double h = std::sqrt(h2);
  double bx = F.a * alpha / d, by = F.b * alpha / d;
  double px1 = bx - F.b / d * h, py1 = by + F.a / d * h;
  double px2 = bx + F.b / d * h, py2 = by - F.a / d * h;
  double px = px1, py = py1;
  if (std::hypot(px2 - hx, py2 - hy) < std::hypot(px1 - hx, py1 - hy)) { px = px2; py = py2; }
  double tx = -(py - F.b), ty = (px - F.a);
  if (tx * dx + ty * dy < 0) { tx = -tx; ty = -ty; }
  psi = std::atan2(ty, tx);
  return true;
}
double qrms(const std::vector<double> &v)      // sqrt(<v^2>): Highland ORDER scale
{
  if (v.empty()) return 0.;
  double s = 0;
  for (double q : v) s += q * q;
  return std::sqrt(s / v.size());
}

// real hit-level grouping shared by nf_hits and nf_ms_hits: pixels
// road-matched to the tracker's ntp_clus_trk seed clusters (same layer,
// dxy < 1.2 cm, |dtbin| <= 6, hitID-deduped, first seed wins).
int realPixGroups(const char *realf, std::vector<Grp> &px)
{
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
      if ((int) ev == 44) continue;   // V6 laser veto (canon.h)
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
  px.assign(nseed, Grp());
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
      if ((int) ev == 44) continue;   // V6 laser veto (canon.h)
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
  return nseed;
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
             const char *ver = "v6")
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
  std::vector<Grp> px;
  int nseed = realPixGroups(realf, px);
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
void nf_clusters(const char *i91 = "island91_frames_production_v6.root",
                 const char *realf = "../clusters_seeds_island_79507-0.root_ntuplizer.root",
                 const char *ver = "v6")
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
      if ((int) ev == 44) continue;   // V6 laser veto (canon.h)
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
  lg->AddEntry(h[1], Form("sim %s ntp_cluster (truth), med %.0f #mum", ver, med(rms[1]) * 1000), "l");
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
               const char *ver = "v6", int showev = 7)
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
      if ((int) ev == 44) continue;   // V6 laser veto (canon.h)
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

// ---------------------------------------------------------------------------
// 4. HIT-LEVEL MS CHECK, no finder: the ms_split.C split-arc tangent-mismatch
//    test at hit level, sim vs real on one canvas.
//    SIM ntp_g4hit (truth-grouped, ms_split method VERBATIM: pT windows from
//    truth momentum, full-crossing gates >=30 pts / rmin<=35 / rmax<=72 rule,
//    halves >=10 pts, |dpsi| < 8 mrad kink veto, Gluckstern-type fit noise,
//    per-track Highland ORDER scale) at TWO borders: r0=35 (the adopted
//    ms_split border — reproduction of the v53f record) and r0=49 (mid pad
//    rows, the only border the real side can share).
//    REAL ntp_hit pixels (tracker-grouped via realPixGroups) at r0=49;
//    windows by FITTED whole-track curvature (no truth exists); NO kink veto
//    (pixel-level fit noise is tens of mrad — the veto would clip the
//    distribution itself); sigma about the sample mean, 3sigma-clipped core.
void nf_ms_hits(const char *g4pat = "../P5/PP_g4hit_%d.root", int ng4 = 10,
                const char *realf = "../clusters_seeds_island_79507-0.root_ntuplizer.root",
                const char *ver = "v6")
{
  using namespace MNF;
  const double PTW[2][2] = {{0.45, 0.55}, {1.5, 2.5}};
  const double RCOEF = 100. / (0.299792458 * 1.4);          // R[cm] per GeV at 1.4 T
  const double RW[2][2] = {{101, 137}, {RCOEF * 1.5, RCOEF * 2.5}};
  const double X0LEN = 11200.;                              // gas X0 [cm]
  const double MPI = 0.13957;
  const double BRD[2] = {35., 49.};
  std::vector<double> sdp[2][2], sfn[2][2], hiw[2];         // sim: [border][window]
  std::vector<double> rdp[2], rfn[2];                       // real: [window]

  // ---------- sim: truth hits, ms_split method, two borders ----------
  struct Trk { std::vector<double> x, y, r; float pt = 0, p = 0; };
  for (int fi = 0; fi < ng4; ++fi)
  {
    TFile *f = TFile::Open(Form(g4pat, fi));
    if (!f || f->IsZombie()) { printf("missing %s\n", Form(g4pat, fi)); continue; }
    TTree *t = (TTree *) f->Get("ntp_g4hit");
    float ev, gx, gy, gpx, gpy, gpz, tid;
    t->SetBranchStatus("*", 0);
    for (auto b : {"event", "gx", "gy", "gpx", "gpy", "gpz", "gtrackID"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("event", &ev);
    t->SetBranchAddress("gx", &gx);
    t->SetBranchAddress("gy", &gy);
    t->SetBranchAddress("gpx", &gpx);
    t->SetBranchAddress("gpy", &gpy);
    t->SetBranchAddress("gpz", &gpz);
    t->SetBranchAddress("gtrackID", &tid);
    std::map<long, Trk> trks;
    for (Long64_t i = 0; i < t->GetEntries(); ++i)
    {
      t->GetEntry(i);
      if (tid <= 0) continue;
      double r = std::hypot(gx, gy);
      if (r < 20 || r > 78) continue;
      double pt = std::hypot(gpx, gpy);
      if (!((pt > PTW[0][0] && pt < PTW[0][1]) || (pt > PTW[1][0] && pt < PTW[1][1]))) continue;
      Trk &T = trks[(long) ev * 100000 + (long) tid];
      T.pt = pt; T.p = std::sqrt(pt * pt + gpz * gpz);
      T.x.push_back(gx); T.y.push_back(gy); T.r.push_back(r);
    }
    for (auto &kv : trks)
    {
      Trk &T = kv.second;
      int w = (T.pt < 1.0) ? 0 : 1;
      double rmin = 1e9, rmax = 0;
      for (double r : T.r) { rmin = std::min(rmin, r); rmax = std::max(rmax, r); }
      if ((int) T.x.size() < 30 || rmin > 35 || rmax < 72) continue;
      bool hlDone = false;
      for (int b = 0; b < 2; ++b)
      {
        const double R0 = BRD[b];
        std::vector<double> xi, yi, xo, yo;
        double hxi = 0, hyi = 0, hxo = 0, hyo = 0, dbi = 1e9, dbo = 1e9;
        for (size_t i = 0; i < T.x.size(); ++i)
        {
          if (T.r[i] < R0)
          {
            xi.push_back(T.x[i]); yi.push_back(T.y[i]);
            if (R0 - T.r[i] < dbi) { dbi = R0 - T.r[i]; hxi = T.x[i]; hyi = T.y[i]; }
          }
          else
          {
            xo.push_back(T.x[i]); yo.push_back(T.y[i]);
            if (T.r[i] - R0 < dbo) { dbo = T.r[i] - R0; hxo = T.x[i]; hyo = T.y[i]; }
          }
        }
        if ((int) xi.size() < 10 || (int) xo.size() < 10) continue;
        Fit Fi = fitCircle(xi, yi), Fo = fitCircle(xo, yo);
        if (!Fi.ok || !Fo.ok) continue;
        double dx = hxo - hxi, dy = hyo - hyi;
        double psi_i, psi_o;
        if (!tangentAtR(Fi, R0, hxi, hyi, dx, dy, psi_i)) continue;
        if (!tangentAtR(Fo, R0, hxo, hyo, dx, dy, psi_o)) continue;
        double dpsi = wrapphi(psi_o - psi_i) * 1e3;
        if (std::fabs(dpsi) > 8) continue;                  // ms_split kink veto
        sdp[b][w].push_back(dpsi);
        double Li = R0 - 20., Lo = 78. - R0;
        double s_i = Fi.rms / Li * std::sqrt(192. / (Fi.n + 4));
        double s_o = Fo.rms / Lo * std::sqrt(192. / (Fo.n + 4));
        sfn[b][w].push_back(std::sqrt(s_i * s_i + s_o * s_o) * 1e3);
        if (!hlDone)
        {
          double beta = T.p / std::sqrt(T.p * T.p + MPI * MPI);
          double path = 58.8 * (T.p / T.pt);
          double xX0 = path / X0LEN;
          hiw[w].push_back(13.6e-3 / (beta * T.p) * std::sqrt(xX0) * (1 + 0.038 * std::log(xX0)) * 1e3);
          hlDone = true;
        }
      }
    }
    f->Close();
    printf("g4 file %d: sim split fits %zu / %zu (r0=49)\n", fi, sdp[1][0].size(), sdp[1][1].size());
  }

  // ---------- real: tracker-grouped pixels, border 49 ----------
  {
    std::vector<Grp> px;
    realPixGroups(realf, px);
    const double R0 = 49.;
    for (auto &G : px)
    {
      if ((int) G.x.size() < 30) continue;
      double rmin = 1e9, rmax = 0;
      for (double r : G.r) { rmin = std::min(rmin, r); rmax = std::max(rmax, r); }
      if (rmin > 35 || rmax < 72) continue;
      Fit Fa = fitCircle(G.x, G.y);
      if (!Fa.ok) continue;
      int w = -1;
      for (int q = 0; q < 2; ++q)
        if (Fa.R > RW[q][0] && Fa.R < RW[q][1]) w = q;
      if (w < 0) continue;
      std::vector<double> xi, yi, xo, yo;
      double hxi = 0, hyi = 0, hxo = 0, hyo = 0, dbi = 1e9, dbo = 1e9;
      for (size_t i = 0; i < G.x.size(); ++i)
      {
        if (G.r[i] < R0)
        {
          xi.push_back(G.x[i]); yi.push_back(G.y[i]);
          if (R0 - G.r[i] < dbi) { dbi = R0 - G.r[i]; hxi = G.x[i]; hyi = G.y[i]; }
        }
        else
        {
          xo.push_back(G.x[i]); yo.push_back(G.y[i]);
          if (G.r[i] - R0 < dbo) { dbo = G.r[i] - R0; hxo = G.x[i]; hyo = G.y[i]; }
        }
      }
      if ((int) xi.size() < 10 || (int) xo.size() < 10) continue;
      Fit Fi = fitCircle(xi, yi), Fo = fitCircle(xo, yo);
      if (!Fi.ok || !Fo.ok) continue;
      double dx = hxo - hxi, dy = hyo - hyi;
      double psi_i, psi_o;
      if (!tangentAtR(Fi, R0, hxi, hyi, dx, dy, psi_i)) continue;
      if (!tangentAtR(Fo, R0, hxo, hyo, dx, dy, psi_o)) continue;
      rdp[w].push_back(wrapphi(psi_o - psi_i) * 1e3);       // NO veto at pixel level
      double Li = R0 - 31.4, Lo = 75.4 - R0;                // pad-row lever arms
      double s_i = Fi.rms / Li * std::sqrt(192. / (Fi.n + 4));
      double s_o = Fo.rms / Lo * std::sqrt(192. / (Fo.n + 4));
      rfn[w].push_back(std::sqrt(s_i * s_i + s_o * s_o) * 1e3);
    }
  }

  // ---------- stats: sigma about the mean, 3sigma-clipped core ----------
  auto stats = [](const std::vector<double> &v, double &sig, double &err,
                  double &core, double &tfrac) {
    sig = err = core = tfrac = 0;
    if (v.size() < 2) return;
    double mu = 0;
    for (double q : v) mu += q;
    mu /= v.size();
    double s2 = 0;
    for (double q : v) s2 += (q - mu) * (q - mu);
    sig = std::sqrt(s2 / v.size());
    err = sig / std::sqrt(2. * v.size());
    double sg = sig;
    long n = 0;
    for (int it = 0; it < 3; ++it)
    {
      double t2 = 0;
      n = 0;
      for (double q : v)
        if (std::fabs(q - mu) < 3 * sg) { t2 += (q - mu) * (q - mu); n++; }
      if (!n) break;
      sg = std::sqrt(t2 / n);
    }
    core = sg;
    tfrac = 1. - (double) n / v.size();
  };
  double ss[2][2], se[2][2], sc3[2][2], st[2][2];           // sim [border][w]
  double rs[2], re[2], rc[2], rt[2];                        // real [w]
  for (int b = 0; b < 2; ++b)
    for (int w = 0; w < 2; ++w) stats(sdp[b][w], ss[b][w], se[b][w], sc3[b][w], st[b][w]);
  for (int w = 0; w < 2; ++w) stats(rdp[w], rs[w], re[w], rc[w], rt[w]);

  FILE *fo = openLedger(ver);
  auto P = [&](const char *fmt, ...) {
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    va_start(ap, fmt); vfprintf(fo, fmt, ap); va_end(ap);
  };
  P("[nf_ms_hits %s] HIT-level split-arc MS check (ms_split method), %d g4 chunks\n", ver, ng4);
  for (int w = 0; w < 2; ++w)
    P("  sim r0=35 pT [%.2f,%.2f]: N=%zu sigma %.3f +- %.3f mrad (core %.3f, fit noise %.2f) "
      "| Highland order %.3f\n",
      PTW[w][0], PTW[w][1], sdp[0][w].size(), ss[0][w], se[0][w], sc3[0][w], med(sfn[0][w]), qrms(hiw[w]));
  P("  [r0=35 reproduction check vs ms_split_v53f record: 1.365 +- 0.008 / 0.367 +- 0.005]\n");
  for (int w = 0; w < 2; ++w)
    P("  sim r0=49 pT [%.2f,%.2f]: N=%zu sigma %.3f +- %.3f mrad (core %.3f, fit noise %.2f)\n",
      PTW[w][0], PTW[w][1], sdp[1][w].size(), ss[1][w], se[1][w], sc3[1][w], med(sfn[1][w]));
  for (int w = 0; w < 2; ++w)
    P("  real r0=49 window %d (%s): N=%zu sigma %.1f +- %.1f mrad (core %.1f, fit noise med %.1f)\n",
      w, w == 0 ? "R_fit 101-137 cm, pT~0.5" : "R_fit 357-595 cm, stiff",
      rdp[w].size(), rs[w], re[w], rc[w], med(rfn[w]));
  P("  1/p check: sim r0=49 ratio %.2f (Highland %.2f) | real ratio %.2f (pT-flat => resolution)\n",
    ss[1][1] > 0 ? ss[1][0] / ss[1][1] : 0,
    qrms(hiw[1]) > 0 ? qrms(hiw[0]) / qrms(hiw[1]) : 0,
    rs[1] > 0 ? rs[0] / rs[1] : 0);
  fclose(fo);

  // ---------- canvas: real (+-100 mrad) | sim (+-8 mrad) ----------
  gStyle->SetOptStat(0);
  TCanvas *cv = new TCanvas("cvnfm", "nf ms hits", 1500, 620);
  cv->Divide(2, 1);
  int wcol[2] = {kBlue + 1, kGreen + 2};
  cv->cd(1);
  {
    TH1D *h[2];
    for (int w = 0; w < 2; ++w)
    {
      h[w] = new TH1D(Form("nfm_r%d", w),
                      ";tangent mismatch at r = 49 cm  #Delta#psi [mrad];tracks (unit area)",
                      50, -100, 100);
      for (double q : rdp[w]) h[w]->Fill(std::max(-99.9, std::min(q, 99.9)));
      if (h[w]->Integral() > 0) h[w]->Scale(1. / h[w]->Integral());
      h[w]->SetLineColor(wcol[w]);
      h[w]->SetLineWidth(2);
    }
    h[0]->SetTitle("REAL ntp_hit pixels, tracker-grouped: split-arc mismatch");
    h[0]->SetMaximum(1.45 * std::max(h[0]->GetMaximum(), h[1]->GetMaximum()));
    h[0]->Draw("hist");
    h[1]->Draw("hist same");
    TLegend *lg = new TLegend(0.13, 0.70, 0.62, 0.87);
    lg->SetBorderSize(0);
    lg->AddEntry(h[0], Form("p_{T}^{fit} 0.45-0.55: #sigma = %.1f#pm%.1f (core %.1f)", rs[0], re[0], rc[0]), "l");
    lg->AddEntry(h[1], Form("p_{T}^{fit} 1.5-2.5:  #sigma = %.1f#pm%.1f (core %.1f)", rs[1], re[1], rc[1]), "l");
    lg->Draw();
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.032);
    tx.DrawLatex(0.13, 0.63, "pT-flat: resolution-driven, MS invisible");
    tx.DrawLatex(0.13, 0.57, "(raw pixels, no centroiding; edge bins collect tails)");
  }
  cv->cd(2);
  {
    TH1D *h[2];
    for (int w = 0; w < 2; ++w)
    {
      h[w] = new TH1D(Form("nfm_s%d", w),
                      ";tangent mismatch at r = 49 cm  #Delta#psi [mrad];tracks (unit area)",
                      81, -8.1, 8.1);
      for (double q : sdp[1][w]) h[w]->Fill(q);
      if (h[w]->Integral() > 0) h[w]->Scale(1. / h[w]->Integral());
      h[w]->SetLineColor(wcol[w]);
      h[w]->SetLineWidth(2);
    }
    h[0]->SetTitle("SIM ntp_g4hit truth hits: split-arc mismatch (ms_split method)");
    h[0]->SetMaximum(1.45 * std::max(h[0]->GetMaximum(), h[1]->GetMaximum()));
    h[0]->Draw("hist");
    h[1]->Draw("hist same");
    TLegend *lg = new TLegend(0.13, 0.70, 0.64, 0.87);
    lg->SetBorderSize(0);
    lg->AddEntry(h[0], Form("p_{T} 0.45-0.55: #sigma = %.3f#pm%.3f (core %.3f)", ss[1][0], se[1][0], sc3[1][0]), "l");
    lg->AddEntry(h[1], Form("p_{T} 1.5-2.5:  #sigma = %.3f#pm%.3f (core %.3f)", ss[1][1], se[1][1], sc3[1][1]), "l");
    lg->Draw();
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.032);
    tx.DrawLatex(0.13, 0.63, Form("1/p ratio %.2f (Highland %.2f): scattering",
                                  ss[1][1] > 0 ? ss[1][0] / ss[1][1] : 0,
                                  qrms(hiw[1]) > 0 ? qrms(hiw[0]) / qrms(hiw[1]) : 0));
    tx.DrawLatex(0.13, 0.57, Form("adopted border 35: #sigma = %.3f (record 1.365)", ss[0][0]));
    tx.DrawLatex(0.13, 0.51, "NOTE x-scale: 12x zoom vs left panel");
  }
  cv->SaveAs(Form("../sim_validation_plots/ms_nofinder_mshits_%s.png", ver));
  printf("wrote ms_nofinder_mshits_%s.png\n", ver);
}

// ---------------------------------------------------------------------------
// 5. HIT-LEVEL SHORT-SAGITTA FIT, no finder (user, 2026-08-13): replace the
//    global whole-track fit by LOCAL fits over 4 ADJACENT PAD ROWS, sliding
//    window start L = 7..51 (45 windows), sim vs real one canvas.
//    Algo corrections vs the proposal (verified): (a) sliding windows overlap
//    -> each point pushed into EVERY window containing its row (up to 4), a
//    single per-point key would leave one row per window; (b) ntp_g4hit has
//    no layer branch -> row = nearest pad-row radius from tpc_geom_table.txt
//    (valid for real pixels too: pixel r IS its row radius); (c) the global
//    bar cannot gate windows (span 1.7-3.3 cm << 15; curvature invisible:
//    sagitta ~100 um at 0.5 GeV under mm noise) -> sample = tracks passing
//    the GLOBAL nf_hits bar, then per-window >=3 distinct rows, n>=5.
//    Expected sim wall: G4 ~1 cm stepping -> 3-6 truth points/window, many
//    unfittable; truth local RMS = the um floor (locally exact arc).
void nf_sag_hits(const char *g4pat = "../P5/PP_g4hit_%d.root", int ng4 = 10,
                 const char *realf = "../clusters_seeds_island_79507-0.root_ntuplizer.root",
                 const char *ver = "v6")
{
  using namespace MNF;
  double rowR[55];
  {
    FILE *fp = fopen("tpc_geom_table.txt", "r");
    if (!fp) { printf("no tpc_geom_table.txt\n"); return; }
    char line[512];
    while (fgets(line, sizeof line, fp))
    {
      int L, nb; double r, sl, p0, p1;
      if (sscanf(line, "%d %d %lf %lf %lf %lf", &L, &nb, &r, &sl, &p0, &p1) == 6 && L >= 7 && L <= 54)
        rowR[L] = r;
    }
    fclose(fp);
  }
  auto nearRow = [&](double r) -> int {
    int best = -1; double bd = 1e9;
    for (int L = 7; L <= 54; ++L)
    {
      double d = std::fabs(r - rowR[L]);
      if (d < bd) { bd = d; best = L; }
    }
    return bd < 0.60 ? best : -1;               // beyond half-pitch of any row
  };
  // per accepted global track: subdivide into windows, fit each
  long nwin[2] = {0, 0}, nfitw[2] = {0, 0};     // [0]=real [1]=sim
  std::vector<double> wrms[2], wn[2];
  auto scanTrack = [&](const Grp &G, int side) {
    std::vector<std::vector<double>> wx(45), wy(45);
    std::vector<std::set<int>> wrow(45);
    for (size_t i = 0; i < G.x.size(); ++i)
    {
      int row = nearRow(G.r[i]);
      if (row < 0) continue;
      for (int w = std::max(7, row - 3); w <= std::min(51, row); ++w)
      {
        wx[w - 7].push_back(G.x[i]);
        wy[w - 7].push_back(G.y[i]);
        wrow[w - 7].insert(row);
      }
    }
    for (int w = 0; w < 45; ++w)
    {
      if (wx[w].empty()) continue;
      nwin[side]++;
      if ((int) wrow[w].size() < 3 || (int) wx[w].size() < 5) continue;
      Fit F = fitCircle(wx[w], wy[w]);
      if (!F.ok) continue;
      nfitw[side]++;
      wrms[side].push_back(F.rms * 1e4);        // um
      wn[side].push_back((double) wx[w].size());
    }
  };
  // --- sim: same sample as nf_hits (global bar), then windows --------------
  long sgrp = 0;
  for (int fi = 0; fi < ng4; ++fi)
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
    for (Long64_t i = 0; i < t->GetEntries(); ++i)
    {
      t->GetEntry(i);
      double r = std::hypot(gx, gy);
      if (r < 20 || r > 78) continue;
      Grp &G = g[{(int) ev, (int) tid}];
      G.x.push_back(gx); G.y.push_back(gy); G.r.push_back(r);
    }
    f->Close();
    for (auto &kv : g)
    {
      Fit F;
      if (!fitBar(kv.second, F)) continue;
      sgrp++;
      scanTrack(kv.second, 1);
    }
    printf("g4 file %d: sim windows fitted %ld\n", fi, nfitw[1]);
  }
  // --- real: same sample as nf_hits (global bar), then windows -------------
  long rgrp = 0;
  {
    std::vector<Grp> px;
    realPixGroups(realf, px);
    for (auto &G : px)
    {
      Fit F;
      if (!fitBar(G, F)) continue;
      rgrp++;
      scanTrack(G, 0);
    }
  }
  FILE *fo = openLedger(ver);
  auto P = [&](const char *fmt, ...) {
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    va_start(ap, fmt); vfprintf(fo, fmt, ap); va_end(ap);
  };
  P("[nf_sag_hits %s] HIT-level SHORT-SAGITTA local fits: 4 adjacent pad rows, 45 sliding windows\n", ver);
  P("  sample = tracks passing the nf_hits GLOBAL bar; window gate >=3 distinct rows, n>=5\n");
  P("  real: %ld tracks -> %ld/%ld windows fittable (%.0f%%) | local RMS med %.0f um | n/window med %.0f\n",
    rgrp, nfitw[0], nwin[0], nwin[0] ? 100. * nfitw[0] / nwin[0] : 0, med(wrms[0]), med(wn[0]));
  P("  sim : %ld tracks -> %ld/%ld windows fittable (%.0f%%) | local RMS med %.2f um | n/window med %.0f\n",
    sgrp, nfitw[1], nwin[1], nwin[1] ? 100. * nfitw[1] / nwin[1] : 0, med(wrms[1]), med(wn[1]));
  P("  (sim fittable fraction limited by G4 ~1 cm stepping: 3-6 truth pts per window;\n");
  P("   local RMS suppressed by sqrt((n-3)/n) dof factor: real ~0.87, sim ~0.63)\n");
  fclose(fo);
  gStyle->SetOptStat(0);
  double e[61];
  for (int i = 0; i <= 60; ++i) e[i] = 0.03 * std::pow(10., 5.523 * i / 60.);   // 0.03 um .. 10 mm
  TH1D *hs = new TH1D("nfg_s", ";per-window circle-fit RMS [#mum];windows (unit area)", 60, e);
  TH1D *hr = new TH1D("nfg_r", ";per-window circle-fit RMS [#mum];windows (unit area)", 60, e);
  for (double v : wrms[1]) hs->Fill(std::min(std::max(v, 0.031), 9990.));
  for (double v : wrms[0]) hr->Fill(std::min(std::max(v, 0.031), 9990.));
  for (auto h : {hs, hr}) if (h->Integral() > 0) h->Scale(1. / h->Integral());
  TCanvas *cv = new TCanvas("cvnfg", "nf sag hits", 1100, 660);
  cv->SetLogx();
  hr->SetLineColor(kBlack); hr->SetLineWidth(2);
  hs->SetLineColor(kBlue + 1); hs->SetLineWidth(2); hs->SetLineStyle(2);
  hr->SetTitle("hit level, no finder: SHORT-SAGITTA local fits (4 adjacent rows)");
  hr->SetMaximum(1.35 * std::max(hr->GetMaximum(), hs->GetMaximum()));
  hr->Draw("hist");
  hs->Draw("hist same");
  TLegend *lg = new TLegend(0.14, 0.72, 0.64, 0.87);
  lg->SetBorderSize(0);
  lg->AddEntry(hr, Form("real pixels, local: med %.0f #mum", med(wrms[0])), "l");
  lg->AddEntry(hs, Form("sim truth hits, local: med %.2f #mum", med(wrms[1])), "l");
  lg->Draw();
  TLatex tx; tx.SetNDC(); tx.SetTextSize(0.030);
  tx.DrawLatex(0.14, 0.66, Form("windows: real %ld, sim %ld (fittable %.0f%% / %.0f%%)",
                                nfitw[0], nfitw[1],
                                nwin[0] ? 100. * nfitw[0] / nwin[0] : 0,
                                nwin[1] ? 100. * nfitw[1] / nwin[1] : 0));
  tx.DrawLatex(0.14, 0.61, "same tracks as the global fit; window gate: >=3 rows, n>=5");
  tx.DrawLatex(0.14, 0.56, "local fit sees point scatter only; curvature invisible over 4 rows");
  cv->SaveAs(Form("../sim_validation_plots/ms_nofinder_saghits_%s.png", ver));
  printf("wrote ms_nofinder_saghits_%s.png\n", ver);
}

// ---------------------------------------------------------------------------
// 6. MATCHED PIXEL-LEVEL COMPARISON (user, 2026-08-13 "do it on v54c"):
//    sim DIGI pixels (digi_frames_production_v53.root — the sealed digi that
//    v5.4c is exported from; per-pixel gtrackID = truth grouping, x/y from
//    row radius x phi) vs real ntp_hit pixels (tracker-grouped). BOTH sides
//    now carry charge clouds -> the like-for-like the supervisor's
//    short-sagitta prediction assumes. One canvas: global fit | 4-row local.
//    DECLARED CAVEATS: (a) sim pixel positions are field-free (the v5.4
//    field is a cluster-export overlay; a pixel-level field needs
//    re-digitization) — local fits reject the smooth field anyway; real
//    granular per-layer alignment scatter (~400 um) survives locally, a ~5%
//    quadrature effect on the real side; (b) grouping asymmetry: sim = the
//    particle's own pixels (truth), real = road-matched pixels around seed
//    clusters (association tails included).
void nf_digipix(const char *digif = "digi_frames_production_v6.root", int nsim = 60,
                const char *realf = "../clusters_seeds_island_79507-0.root_ntuplizer.root",
                const char *ver = "v6", const char *dump = "")
{
  using namespace MNF;
  double rowR[55];
  {
    FILE *fp = fopen("tpc_geom_table.txt", "r");
    if (!fp) { printf("no tpc_geom_table.txt\n"); return; }
    char line[512];
    while (fgets(line, sizeof line, fp))
    {
      int L, nb; double r, sl, p0, p1;
      if (sscanf(line, "%d %d %lf %lf %lf %lf", &L, &nb, &r, &sl, &p0, &p1) == 6 && L >= 7 && L <= 54)
        rowR[L] = r;
    }
    fclose(fp);
  }
  auto nearRow = [&](double r) -> int {
    int best = -1; double bd = 1e9;
    for (int L = 7; L <= 54; ++L)
    {
      double d = std::fabs(r - rowR[L]);
      if (d < bd) { bd = d; best = L; }
    }
    return bd < 0.60 ? best : -1;
  };
  // one pass per side: global fit + 4-row sliding windows (same as nf_sag_hits)
  long ngrp[2] = {0, 0}, nwin[2] = {0, 0}, nfitw[2] = {0, 0}, npxs[2] = {0, 0};
  std::vector<double> grms[2], wrms[2];              // [0]=real [1]=sim, um
  auto doTrack = [&](const Grp &G, int side) {
    Fit F;
    if (!fitBar(G, F)) return;
    ngrp[side]++;
    npxs[side] += (long) G.x.size();
    grms[side].push_back(F.rms * 1e4);
    std::vector<std::vector<double>> wx(45), wy(45);
    std::vector<std::set<int>> wrow(45);
    for (size_t i = 0; i < G.x.size(); ++i)
    {
      int row = nearRow(G.r[i]);
      if (row < 0) continue;
      for (int w = std::max(7, row - 3); w <= std::min(51, row); ++w)
      {
        wx[w - 7].push_back(G.x[i]);
        wy[w - 7].push_back(G.y[i]);
        wrow[w - 7].insert(row);
      }
    }
    for (int w = 0; w < 45; ++w)
    {
      if (wx[w].empty()) continue;
      nwin[side]++;
      if ((int) wrow[w].size() < 3 || (int) wx[w].size() < 5) continue;
      Fit L = fitCircle(wx[w], wy[w]);
      if (!L.ok) continue;
      nfitw[side]++;
      wrms[side].push_back(L.rms * 1e4);
    }
  };
  // --- real: tracker-grouped pixels ----------------------------------------
  {
    std::vector<Grp> px;
    realPixGroups(realf, px);
    for (auto &G : px) doTrack(G, 0);
  }
  printf("real done: %ld tracks\n", ngrp[0]);
  // --- sim: digi pixels, per-pixel truth grouping --------------------------
  {
    TFile *f = TFile::Open(digif);
    if (!f || f->IsZombie()) { printf("no %s\n", digif); return; }
    TTree *t = (TTree *) f->Get("ntp_hit");
    float ev, lay, phi, adc, tid;
    t->SetBranchStatus("*", 0);
    for (auto b : {"event", "layer", "phi", "adc", "gtrackID"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("event", &ev);
    t->SetBranchAddress("layer", &lay);
    t->SetBranchAddress("phi", &phi);
    t->SetBranchAddress("adc", &adc);
    t->SetBranchAddress("gtrackID", &tid);
    std::map<std::pair<int, int>, Grp> g;
    for (Long64_t i = 0; i < t->GetEntries(); ++i)
    {
      t->GetEntry(i);
      if ((int) ev >= nsim) continue;
      if (lay < 7 || lay > 54 || adc <= 0 || tid <= 0) continue;
      double r = rowR[(int) lay];
      Grp &G = g[{(int) ev, (int) tid}];
      G.x.push_back(r * std::cos(phi));
      G.y.push_back(r * std::sin(phi));
      G.r.push_back(r);
    }
    f->Close();
    printf("sim digi: %zu truth groups in %d frames\n", g.size(), nsim);
    for (auto &kv : g) doTrack(kv.second, 1);
  }
  if (dump && dump[0])                       // per-fit RMS values for external (clean-label) figures
  {
    FILE *fd = fopen(dump, "w");
    if (fd)
    {
      for (int s = 0; s < 2; ++s)
      {
        for (double q : grms[s]) fprintf(fd, "G %d %.1f\n", s, q);
        for (double q : wrms[s]) fprintf(fd, "L %d %.1f\n", s, q);
      }
      fclose(fd);
      printf("dumped per-fit RMS values to %s\n", dump);
    }
  }
  FILE *fo = openLedger(ver);
  auto P = [&](const char *fmt, ...) {
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    va_start(ap, fmt); vfprintf(fo, fmt, ap); va_end(ap);
  };
  P("[nf_digipix %s] MATCHED pixel level: sim DIGI pixels (%s, truth-grouped) "
    "vs real pixels (tracker-grouped)\n", ver, digif);
  for (int s = 0; s < 2; ++s)
    P("  %s: %ld tracks (%.0f px/track) | GLOBAL RMS med %.0f um | LOCAL (4-row) med %.0f um "
      "(%ld/%ld windows, %.0f%%)\n",
      s == 0 ? "real" : "sim ", ngrp[s], ngrp[s] ? (double) npxs[s] / ngrp[s] : 0,
      med(grms[s]), med(wrms[s]), nfitw[s], nwin[s], nwin[s] ? 100. * nfitw[s] / nwin[s] : 0);
  P("  data/MC: GLOBAL %.2f | LOCAL %.2f   (supervisor test: matched local = response ratio)\n",
    med(grms[1]) > 0 ? med(grms[0]) / med(grms[1]) : 0,
    med(wrms[1]) > 0 ? med(wrms[0]) / med(wrms[1]) : 0);
  P("  global-vs-local gap: real %.0f -> %.0f | sim %.0f -> %.0f  "
    "(the removed part = each side's long-range share)\n",
    med(grms[0]), med(wrms[0]), med(grms[1]), med(wrms[1]));
  fclose(fo);
  gStyle->SetOptStat(0);
  TCanvas *cv = new TCanvas("cvnfd", "nf digipix", 1500, 620);
  cv->Divide(2, 1);
  const char *pt[2] = {"GLOBAL whole-track fit", "LOCAL 4-row short-sagitta fit"};
  for (int p = 0; p < 2; ++p)
  {
    cv->cd(p + 1);
    std::vector<double> *v[2] = {p == 0 ? &grms[0] : &wrms[0], p == 0 ? &grms[1] : &wrms[1]};
    double xhi = p == 0 ? 6000 : 4000;
    TH1D *h[2];
    for (int s = 0; s < 2; ++s)
    {
      h[s] = new TH1D(Form("nfd_%d_%d", p, s),
                      ";per-fit circle RMS [#mum];fits (unit area)", 60, 0, xhi);
      for (double q : *v[s]) h[s]->Fill(std::min(q, xhi - 1));
      if (h[s]->Integral() > 0) h[s]->Scale(1. / h[s]->Integral());
    }
    h[0]->SetLineColor(kBlack); h[0]->SetLineWidth(2);
    h[1]->SetLineColor(kBlue + 1); h[1]->SetLineWidth(2); h[1]->SetLineStyle(2);
    h[0]->SetTitle(Form("matched pixel level: %s", pt[p]));
    h[0]->SetMaximum(1.4 * std::max(h[0]->GetMaximum(), h[1]->GetMaximum()));
    h[0]->Draw("hist");
    h[1]->Draw("hist same");
    TLegend *lg = new TLegend(0.44, 0.72, 0.89, 0.87);
    lg->SetBorderSize(0);
    lg->AddEntry(h[0], Form("real pixels, med %.0f #mum", med(*v[0])), "l");
    lg->AddEntry(h[1], Form("sim digi pixels, med %.0f #mum", med(*v[1])), "l");
    lg->Draw();
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.030);
    tx.DrawLatex(0.44, 0.66, Form("data/MC %.2f",
                                  med(*v[1]) > 0 ? med(*v[0]) / med(*v[1]) : 0));
    if (p == 0)
    {
      tx.SetTextSize(0.024);
      tx.DrawLatex(0.44, 0.61, Form("sim digi: %s", digif));
      tx.SetTextSize(0.030);
    }
    else
    {
      tx.DrawLatex(0.44, 0.61, "local fits reject the smooth field on both sides:");
      tx.DrawLatex(0.44, 0.56, "this ratio is the RESPONSE, distortion-free");
    }
  }
  cv->SaveAs(Form("../sim_validation_plots/ms_nofinder_digipix_%s.png", ver));
  printf("wrote ms_nofinder_digipix_%s.png\n", ver);
}
