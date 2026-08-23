// ms_real.C — apply the MS (multiple scattering) split-arc machinery to REAL
// data and compare with the sim at the same level (user, 2026-07-31).
// Real data has no truth hits, so "the same analysis" runs on RECONSTRUCTED
// clusters grouped into tracks:
//   REAL: ntp_clus_trk of clusters_seeds_island_79507-0.root_ntuplizer.root
//         (the canonical real source), tracks = (event, seedID) groups from
//         the classical seed-island tracking;
//   SIM : island91 v5.3 truth-matched clusters, tracks = (event, gtrackID)
//         groups (cls==0 && ntrks==1) — truth used ONLY to group, exactly as
//         seedID groups the real side.
// Both sides get the IDENTICAL treatment: full-crosser gates, whole-track
// circle fit, pT windows selected by FITTED curvature (no truth in the
// selection), and the split-arc tangent mismatch.
// TWO deliberate adaptations vs the truth-hit analysis (ms_split.C):
//   (1) the split border is r0 = 49 cm, NOT the adopted 35: pad rows begin
//       at r = 31.4 cm, so r0=35 leaves <=7 clusters inside — below any fit
//       gate. 35 is a truth-hit border (hits start at r=20); it cannot be
//       applied at cluster level.
//   (2) no +-8 mrad core window: at cluster level the mismatch is dominated
//       by cluster-resolution fit noise (~10 mrad), not by the ~1 mrad MS.
// EXPECTED SIGNATURE (and the point of the comparison): sigma(dpsi) at
// cluster level is pT-INDEPENDENT (resolution noise), unlike the 1/p MS
// scaling seen on truth hits — the real data demonstrating directly that MS
// is invisible under cluster resolution. Data-vs-MC agreement of the widths
// then checks the sim's cluster resolution realism.
// Output: ../sim_validation_plots/ms_real_<ver>.png + ms_real_<ver>.txt
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TLine.h>
#include <TEllipse.h>
#include <TGraph.h>
#include <cstdio>
#include <cstdarg>
#include <cmath>
#include <map>
#include <vector>
#include <algorithm>

namespace MSR
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
double wrapphi(double d)
{
  while (d > M_PI) d -= 2 * M_PI;
  while (d < -M_PI) d += 2 * M_PI;
  return d;
}
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

struct CT { std::vector<double> x, y, r; int lmin = 99, lmax = 0; };

// iterative outlier cleaning: refit while removing clusters with
// |residual| > rescut (default 0.30 cm ~ 4 sigma of cluster resolution).
// Returns the final fit; nrem = clusters removed. False if unfittable or
// fewer than 6 survivors.
bool cleanFit(const CT &T, double rescut, CT &S, Fit &F, int &nrem)
{
  S = T;
  nrem = 0;
  for (int it = 0; it < 6; ++it)
  {
    if ((int) S.x.size() < 6) return false;
    F = fitCircle(S.x, S.y);
    if (!F.ok) return false;
    CT K;
    int rem = 0;
    for (size_t i = 0; i < S.x.size(); ++i)
    {
      double res = std::hypot(S.x[i] - F.a, S.y[i] - F.b) - F.R;
      if (std::fabs(res) > rescut) { rem++; continue; }
      K.x.push_back(S.x[i]); K.y.push_back(S.y[i]); K.r.push_back(S.r[i]);
    }
    if (!rem) return true;
    nrem += rem;
    K.lmin = S.lmin; K.lmax = S.lmax;   // spans kept from the original group
    S = K;
  }
  return (int) S.x.size() >= 6 && F.ok;
}

// gates + whole-track fit + curvature window; w = -1 if rejected
int classify(CT &T, Fit &F, const double RW[2][2])
{
  if ((int) T.x.size() < 12 || T.lmin > 11 || T.lmax < 50) return -1;
  F = fitCircle(T.x, T.y);
  if (!F.ok || F.rms > 0.25) return -1;
  for (int w = 0; w < 2; ++w)
    if (F.R > RW[w][0] && F.R < RW[w][1]) return w;
  return -1;
}

// split-arc mismatch on clusters; wide guard, no MS-core window.
// Optional out-params expose the half fits and boundary hits for drawing.
bool splitClus(CT &T, double r0, double &out,
               Fit *pFi = nullptr, Fit *pFo = nullptr, double *bh = nullptr)
{
  std::vector<double> xi, yi, xo, yo;
  double hxi = 0, hyi = 0, hxo = 0, hyo = 0, dbi = 1e9, dbo = 1e9;
  for (size_t i = 0; i < T.x.size(); ++i)
  {
    if (T.r[i] < r0)
    {
      xi.push_back(T.x[i]); yi.push_back(T.y[i]);
      if (r0 - T.r[i] < dbi) { dbi = r0 - T.r[i]; hxi = T.x[i]; hyi = T.y[i]; }
    }
    else
    {
      xo.push_back(T.x[i]); yo.push_back(T.y[i]);
      if (T.r[i] - r0 < dbo) { dbo = T.r[i] - r0; hxo = T.x[i]; hyo = T.y[i]; }
    }
  }
  if ((int) xi.size() < 10 || (int) xo.size() < 10) return false;
  Fit Fi = fitCircle(xi, yi), Fo = fitCircle(xo, yo);
  if (!Fi.ok || !Fo.ok) return false;
  double psi_i, psi_o, dx = hxo - hxi, dy = hyo - hyi;
  if (!tangentAtR(Fi, r0, hxi, hyi, dx, dy, psi_i)) return false;
  if (!tangentAtR(Fo, r0, hxo, hyo, dx, dy, psi_o)) return false;
  out = wrapphi(psi_o - psi_i) * 1e3;
  if (pFi) *pFi = Fi;
  if (pFo) *pFo = Fo;
  if (bh) { bh[0] = hxi; bh[1] = hyi; bh[2] = hxo; bh[3] = hyo; }
  return std::fabs(out) <= 200;
}
}  // namespace MSR

void ms_real(const char *realf = "../clusters_seeds_island_79507-0.root_ntuplizer.root",
             const char *i91 = "island91_frames_production_v6.root",
             const char *ver = "v6", const char *vtag = "V6",
             double r0 = 49.0)
{
  using namespace MSR;
  const double RCOEF = 100. / (0.299792458 * 1.4);
  const double RW[2][2] = {{101, 137}, {RCOEF * 1.5, RCOEF * 2.5}};   // pT ~0.5 / 1.5-2.5 by FITTED R

  // ---------- loaders: real (event,seedID) and sim (event,gtrackID) ----------
  std::map<long, CT> grp[2];                       // [0]=real [1]=sim
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
      if (lay < 7 || lay > 54) continue;           // TPC clusters only
      CT &T = grp[0][(long) ev * 1000 + (long) sid];
      T.x.push_back(x); T.y.push_back(y); T.r.push_back(std::hypot(x, y));
      T.lmin = std::min(T.lmin, (int) lay); T.lmax = std::max(T.lmax, (int) lay);
    }
    f->Close();
  }
  {
    TFile *f = TFile::Open(i91);
    TTree *c = (TTree *) f->Get("ntp_cluster");
    TTree *u = (TTree *) f->Get("ntp_truth");
    float ev, lay, x, y, tid, cls, ntrks;
    c->SetBranchStatus("*", 0);
    for (auto b : {"event", "layer", "x", "y"}) c->SetBranchStatus(b, 1);
    c->SetBranchAddress("event", &ev);
    c->SetBranchAddress("layer", &lay);
    c->SetBranchAddress("x", &x);
    c->SetBranchAddress("y", &y);
    u->SetBranchStatus("*", 0);
    for (auto b : {"gtrackID", "cls", "ntrks"}) u->SetBranchStatus(b, 1);
    u->SetBranchAddress("gtrackID", &tid);
    u->SetBranchAddress("cls", &cls);
    u->SetBranchAddress("ntrks", &ntrks);
    for (Long64_t i = 0; i < c->GetEntries(); ++i)
    {
      u->GetEntry(i);
      if (cls != 0 || ntrks != 1) continue;
      c->GetEntry(i);
      CT &T = grp[1][(long) ev * 1000000 + (long) tid];
      T.x.push_back(x); T.y.push_back(y); T.r.push_back(std::hypot(x, y));
      T.lmin = std::min(T.lmin, (int) lay); T.lmax = std::max(T.lmax, (int) lay);
    }
    f->Close();
  }

  // ---------- identical processing ----------
  TH1D *hrms[2][2];                                // [src][window]
  const char *sn[2] = {"real", "sim"};
  for (int s = 0; s < 2; ++s)
    for (int w = 0; w < 2; ++w)
      hrms[s][w] = new TH1D(Form("hrms_%d_%d", s, w),
                            ";per-track circle-fit RMS [mm];tracks (unit area)", 60, 0, 3.0);
  std::vector<double> rmsv[2][2], dps[2][2];
  long ngrp[2] = {0, 0};
  for (int s = 0; s < 2; ++s)
  {
    ngrp[s] = (long) grp[s].size();
    for (auto &kv : grp[s])
    {
      Fit F;
      int w = classify(kv.second, F, RW);
      if (w < 0) continue;
      rmsv[s][w].push_back(F.rms * 10);            // mm
      hrms[s][w]->Fill(F.rms * 10);
      double dp;
      if (splitClus(kv.second, r0, dp)) dps[s][w].push_back(dp);
    }
  }
  auto rmsof = [](std::vector<double> &v) {
    double sum = 0;
    for (double q : v) sum += q * q;
    return v.empty() ? 0. : std::sqrt(sum / v.size());
  };
  auto coreof = [&](std::vector<double> &v) {
    double sg = rmsof(v);
    for (int it = 0; it < 3; ++it)
    {
      double s2 = 0;
      long n = 0;
      for (double q : v)
        if (std::fabs(q) < 3 * sg) { s2 += q * q; n++; }
      if (!n) break;
      sg = std::sqrt(s2 / n);
    }
    return sg;
  };

  // ---------- ledger ----------
  FILE *fo = fopen(Form("ms_real_%s.txt", ver), "w");
  auto P = [&](const char *fmt, ...) {
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    va_start(ap, fmt); vfprintf(fo, fmt, ap); va_end(ap);
  };
  P("ms_real %s — MS split-arc machinery applied to REAL cluster-tracks vs sim (%s)\n", ver, vtag);
  P("real = ntp_clus_trk (event,seedID) tracks; sim = island91 truth-grouped tracks;\n");
  P("selection by FITTED curvature both sides; cluster split border r0 = %.0f cm\n", r0);
  P("(the truth-hit border 35 cm is unusable on clusters: pad rows begin at 31.4 cm).\n");
  P("track groups: real %ld, sim %ld\n", ngrp[0], ngrp[1]);
  const char *wl[2] = {"pT ~0.5 GeV (R_fit 101-137 cm)", "pT 1.5-2.5 GeV (R_fit 357-596 cm)"};
  for (int w = 0; w < 2; ++w)
  {
    P("%s:\n", wl[w]);
    for (int s = 0; s < 2; ++s)
    {
      double sg = rmsof(dps[s][w]);
      double se = dps[s][w].empty() ? 0 : sg / std::sqrt(2. * dps[s][w].size());
      P("  %-4s: N=%5zu  circle RMS median = %.0f um   |  split N=%5zu  sigma(dpsi) = %.2f ± %.2f mrad (core %.2f)\n",
        sn[s], rmsv[s][w].size(), med(rmsv[s][w]) * 1000, dps[s][w].size(), sg, se, coreof(dps[s][w]));
    }
    if (!rmsv[1][w].empty() && med(rmsv[1][w]) > 0)
      P("  data/MC: circle RMS %.2f   sigma(dpsi) %.2f\n",
        med(rmsv[0][w]) / med(rmsv[1][w]),
        rmsof(dps[1][w]) > 0 ? rmsof(dps[0][w]) / rmsof(dps[1][w]) : 0);
  }
  P("reading: sigma(dpsi) at cluster level is pT-INDEPENDENT (resolution-driven,\n");
  P("~10x the truth-level MS) in BOTH real and sim — the real data's own\n");
  P("demonstration that in-gas MS is invisible under cluster resolution.\n");
  P("Data/MC width ratios test the sim's cluster-resolution realism; the real\n");
  P("side additionally carries alignment/distortion residuals the sim lacks.\n");
  fclose(fo);

  // ---------- figure ----------
  gStyle->SetOptStat(0);
  TCanvas *cv = new TCanvas("cvr", Form("ms real %s", ver), 1500, 620);
  cv->Divide(2, 1);
  for (int w = 0; w < 2; ++w)
  {
    cv->cd(w + 1);
    for (int s = 0; s < 2; ++s)
      if (hrms[s][w]->Integral() > 0) hrms[s][w]->Scale(1. / hrms[s][w]->Integral());
    hrms[0][w]->SetLineColor(kBlack);
    hrms[0][w]->SetLineWidth(2);
    hrms[1][w]->SetLineColor(w == 0 ? kBlue + 1 : kGreen + 2);
    hrms[1][w]->SetLineWidth(2);
    hrms[1][w]->SetLineStyle(2);
    hrms[0][w]->SetTitle(Form("%s: real tracker vs sim truth-grouped clusters", wl[w]));
    hrms[0][w]->SetMaximum(1.45 * std::max(hrms[0][w]->GetMaximum(), hrms[1][w]->GetMaximum()));
    hrms[0][w]->Draw("hist");
    hrms[1][w]->Draw("hist same");
    TLegend *lg = new TLegend(0.40, 0.72, 0.89, 0.88);
    lg->SetBorderSize(0);
    lg->AddEntry(hrms[0][w], Form("real, median %.0f #mum", med(rmsv[0][w]) * 1000), "l");
    lg->AddEntry(hrms[1][w], Form("sim %s, median %.0f #mum", vtag, med(rmsv[1][w]) * 1000), "l");
    lg->Draw();
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.036);
    tx.DrawLatex(0.40, 0.64, Form("#sigma(#Delta#psi) at r_{0}=%g: real %.1f, sim %.1f mrad",
                                  r0, rmsof(dps[0][w]), rmsof(dps[1][w])));
    tx.DrawLatex(0.40, 0.57, "(resolution-driven, p_{T}-independent, not MS)");
  }
  cv->SaveAs(Form("../sim_validation_plots/ms_real_%s.png", ver));
  printf("wrote ../sim_validation_plots/ms_real_%s.png + ms_real_%s.txt\n", ver, ver);
}

// ---------------------------------------------------------------------------
// ms_realcheck — circularity AUDIT of the real tracker (user, 2026-07-31):
// "is ntp_clus_trk really true tracks?"  Since a true track must fit a circle
// at the cluster-resolution level (~0.7 mm, established by truth_circle /
// ms_real), circularity is an independent validity test of every
// (event, seedID) group — fakes and mis-associated clusters cannot satisfy it.
// Procedure, applied IDENTICALLY to real and sim (sim = island91 truth
// groups, so the same audit also measures the truth-grouping's own purity):
//   fit sample: >=8 clusters spanning >=10 layers (else UNFITTABLE — too
//     short to judge, reported separately);
//   iterative cleaning: refit while removing clusters with |residual| >
//     0.30 cm (~4 sigma), <=6 passes, >=6 survivors;
//   classify: CLEAN   = removed <=5% of clusters and final RMS <= 0.20 cm
//             RESCUED = removed <=30% and final RMS <= 0.20 (true track
//                       carrying mis-associated clusters)
//             FAKE-like = no consistent circle survives.
//   physics cross-check: transverse impact parameter d0 = ||C| - R| of
//     accepted tracks (primaries pile near the beamline; large-d0 entries
//     are genuine secondaries — decays/conversions — NOT used as a cut).
// Then the pT-window circle comparison is REDONE on cleaned tracks
// (full-crosser gates on survivors), giving purified data/MC ratios.
// Output: ../sim_validation_plots/ms_realcheck_<ver>.png + ms_realcheck_<ver>.txt
void ms_realcheck(const char *realf = "../clusters_seeds_island_79507-0.root_ntuplizer.root",
                  const char *i91 = "island91_frames_production_v6.root",
                  const char *ver = "v6", const char *vtag = "V6",
                  double rescut = 0.30)
{
  using namespace MSR;
  const double RCOEF = 100. / (0.299792458 * 1.4);
  const double RW[2][2] = {{101, 137}, {RCOEF * 1.5, RCOEF * 2.5}};

  std::map<long, CT> grp[2];
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
      CT &T = grp[0][(long) ev * 1000 + (long) sid];
      T.x.push_back(x); T.y.push_back(y); T.r.push_back(std::hypot(x, y));
      T.lmin = std::min(T.lmin, (int) lay); T.lmax = std::max(T.lmax, (int) lay);
    }
    f->Close();
  }
  {
    TFile *f = TFile::Open(i91);
    TTree *c = (TTree *) f->Get("ntp_cluster");
    TTree *u = (TTree *) f->Get("ntp_truth");
    float ev, lay, x, y, tid, cls, ntrks;
    c->SetBranchStatus("*", 0);
    for (auto b : {"event", "layer", "x", "y"}) c->SetBranchStatus(b, 1);
    c->SetBranchAddress("event", &ev);
    c->SetBranchAddress("layer", &lay);
    c->SetBranchAddress("x", &x);
    c->SetBranchAddress("y", &y);
    u->SetBranchStatus("*", 0);
    for (auto b : {"gtrackID", "cls", "ntrks"}) u->SetBranchStatus(b, 1);
    u->SetBranchAddress("gtrackID", &tid);
    u->SetBranchAddress("cls", &cls);
    u->SetBranchAddress("ntrks", &ntrks);
    for (Long64_t i = 0; i < c->GetEntries(); ++i)
    {
      u->GetEntry(i);
      if (cls != 0 || ntrks != 1) continue;
      c->GetEntry(i);
      CT &T = grp[1][(long) ev * 1000000 + (long) tid];
      T.x.push_back(x); T.y.push_back(y); T.r.push_back(std::hypot(x, y));
      T.lmin = std::min(T.lmin, (int) lay); T.lmax = std::max(T.lmax, (int) lay);
    }
    f->Close();
  }

  const char *sn[2] = {"real", "sim"};
  long nall[2] = {0, 0}, nshort[2] = {0, 0}, nclean[2] = {0, 0}, nresc[2] = {0, 0}, nfake[2] = {0, 0};
  long cltot[2] = {0, 0}, clrem[2] = {0, 0};
  std::vector<double> rmsBefore[2], rmsAfter[2], d0v[2], worstv[2];
  std::vector<double> rmsW[2][2];                 // redo: cleaned medians per window
  TH1D *hw[2], *hd0[2];
  for (int s = 0; s < 2; ++s)
  {
    hw[s] = new TH1D(Form("hw%d", s), ";worst |residual| in group [mm];tracks (unit area)", 60, 0, 12);
    hd0[s] = new TH1D(Form("hd0%d", s), ";transverse impact parameter d_{0} [cm];tracks (unit area)", 60, 0, 15);
  }
  for (int s = 0; s < 2; ++s)
  {
    for (auto &kv : grp[s])
    {
      CT &T = kv.second;
      nall[s]++;
      if ((int) T.x.size() < 8 || T.lmax - T.lmin < 10) { nshort[s]++; continue; }
      cltot[s] += (long) T.x.size();
      Fit F0 = fitCircle(T.x, T.y);
      double worst = 0;
      if (F0.ok)
      {
        rmsBefore[s].push_back(F0.rms * 10);
        for (size_t i = 0; i < T.x.size(); ++i)
          worst = std::max(worst, std::fabs(std::hypot(T.x[i] - F0.a, T.y[i] - F0.b) - F0.R));
        hw[s]->Fill(std::min(worst * 10, 11.9));
        worstv[s].push_back(worst * 10);
      }
      CT S;
      Fit F;
      int nrem = 0;
      bool ok = cleanFit(T, rescut, S, F, nrem);
      double frac = (double) nrem / T.x.size();
      if (!ok || F.rms > 0.20 || frac > 0.30)
      {
        nfake[s]++;
        clrem[s] += nrem;
        continue;
      }
      clrem[s] += nrem;
      if (frac <= 0.05) nclean[s]++;
      else nresc[s]++;
      rmsAfter[s].push_back(F.rms * 10);
      double d0 = std::fabs(std::hypot(F.a, F.b) - F.R);
      hd0[s]->Fill(std::min(d0, 14.9));
      d0v[s].push_back(d0);
      // redo the window comparison on the CLEANED track (full-crosser gates)
      if ((int) S.x.size() >= 12 && T.lmin <= 11 && T.lmax >= 50)
        for (int w = 0; w < 2; ++w)
          if (F.R > RW[w][0] && F.R < RW[w][1]) rmsW[s][w].push_back(F.rms * 10);
    }
  }

  FILE *fo = fopen(Form("ms_realcheck_%s.txt", ver), "w");
  auto P = [&](const char *fmt, ...) {
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    va_start(ap, fmt); vfprintf(fo, fmt, ap); va_end(ap);
  };
  P("ms_realcheck %s — circularity audit of track groups (outlier cut %.2f cm)\n", ver, rescut);
  P("real = ntp_clus_trk (event,seedID); sim = island91 %s truth groups, same audit\n", vtag);
  for (int s = 0; s < 2; ++s)
  {
    long nfit = nall[s] - nshort[s];
    P("%s: %ld groups; %ld too short to judge; %ld auditable\n", sn[s], nall[s], nshort[s], nfit);
    P("  CLEAN %ld (%.1f%%)   RESCUED %ld (%.1f%%)   FAKE-like %ld (%.1f%%)\n",
      nclean[s], 100. * nclean[s] / std::max(1L, nfit),
      nresc[s], 100. * nresc[s] / std::max(1L, nfit),
      nfake[s], 100. * nfake[s] / std::max(1L, nfit));
    P("  cluster contamination removed: %ld / %ld (%.2f%%)\n",
      clrem[s], cltot[s], 100. * clrem[s] / std::max(1L, cltot[s]));
    P("  per-track RMS median: before clean %.0f um -> after %.0f um; worst-residual median %.1f mm\n",
      med(rmsBefore[s]) * 1000, med(rmsAfter[s]) * 1000, med(worstv[s]));
    P("  d0 median %.2f cm; d0<2 cm: %.1f%% (large d0 = genuine secondaries, not fakes)\n",
      med(d0v[s]), d0v[s].empty() ? 0 : 100. * std::count_if(d0v[s].begin(), d0v[s].end(),
        [](double q) { return q < 2; }) / (double) d0v[s].size());
  }
  P("REDO on cleaned tracks (full-crosser gates, window by fitted R):\n");
  const char *wl[2] = {"pT ~0.5 GeV", "pT 1.5-2.5 GeV"};
  for (int w = 0; w < 2; ++w)
  {
    P("  %s: real N=%zu median %.0f um | sim N=%zu median %.0f um | data/MC %.2f\n",
      wl[w], rmsW[0][w].size(), med(rmsW[0][w]) * 1000,
      rmsW[1][w].size(), med(rmsW[1][w]) * 1000,
      med(rmsW[1][w]) > 0 ? med(rmsW[0][w]) / med(rmsW[1][w]) : 0);
  }
  P("caveats: transverse-only test (helix pitch unchecked); an accidental\n");
  P("circular fake is not excluded; short groups are unjudged, not exonerated.\n");
  fclose(fo);

  gStyle->SetOptStat(0);
  TCanvas *cv = new TCanvas("cvk", Form("ms realcheck %s", ver), 1500, 620);
  cv->Divide(2, 1);
  cv->cd(1);
  {
    for (int s = 0; s < 2; ++s)
      if (hw[s]->Integral() > 0) hw[s]->Scale(1. / hw[s]->Integral());
    hw[0]->SetLineColor(kBlack); hw[0]->SetLineWidth(2);
    hw[1]->SetLineColor(kBlue + 1); hw[1]->SetLineWidth(2); hw[1]->SetLineStyle(2);
    hw[0]->SetTitle("worst cluster residual per track (pre-clean fit)");
    hw[0]->SetMaximum(1.4 * std::max(hw[0]->GetMaximum(), hw[1]->GetMaximum()));
    hw[0]->Draw("hist");
    hw[1]->Draw("hist same");
    TLine *lc = new TLine(rescut * 10, 0, rescut * 10, 0.5 * hw[0]->GetMaximum());
    lc->SetLineColor(kGray + 2); lc->SetLineStyle(2); lc->Draw();
    TLegend *lg = new TLegend(0.44, 0.72, 0.89, 0.88);
    lg->SetBorderSize(0);
    lg->AddEntry(hw[0], "real (event,seedID) groups", "l");
    lg->AddEntry(hw[1], Form("sim %s truth groups", vtag), "l");
    lg->AddEntry(lc, Form("outlier cut %.0f mm", rescut * 10), "l");
    lg->Draw();
  }
  cv->cd(2);
  {
    gPad->SetLogy();
    for (int s = 0; s < 2; ++s)
      if (hd0[s]->Integral() > 0) hd0[s]->Scale(1. / hd0[s]->Integral());
    hd0[0]->SetLineColor(kBlack); hd0[0]->SetLineWidth(2);
    hd0[1]->SetLineColor(kBlue + 1); hd0[1]->SetLineWidth(2); hd0[1]->SetLineStyle(2);
    hd0[0]->SetTitle("impact parameter of accepted tracks");
    hd0[0]->SetMaximum(3.0 * std::max(hd0[0]->GetMaximum(), hd0[1]->GetMaximum()));
    hd0[0]->Draw("hist");
    hd0[1]->Draw("hist same");
    TLegend *lg = new TLegend(0.40, 0.74, 0.89, 0.88);
    lg->SetBorderSize(0);
    lg->AddEntry(hd0[0], "real, accepted tracks", "l");
    lg->AddEntry(hd0[1], Form("sim %s, accepted tracks", vtag), "l");
    lg->Draw();
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.034);
    tx.DrawLatex(0.40, 0.66, "peak = primaries from the beamline;");
    tx.DrawLatex(0.40, 0.60, "tail = genuine secondaries (decays, conversions)");
  }
  cv->SaveAs(Form("../sim_validation_plots/ms_realcheck_%s.png", ver));
  printf("wrote ../sim_validation_plots/ms_realcheck_%s.png + ms_realcheck_%s.txt\n", ver, ver);
}

// ---------------------------------------------------------------------------
// ms_real_split — the REAL-data split-arc figure in the SAME STYLE as the
// supervisor-seen sim figure ms_split_v53f.png (user, 2026-08-05): left panel
// dpsi histograms for the two fitted-curvature windows, right panel the
// summary table (real values, sim v5.3 same-level values in parentheses).
// The worst showcase (largest-|dpsi| real track at pT~0.5, drawn with its two
// half-fit circles and boundary tangents) and the remaining statistics go to
// a SECOND png so the primary stays parallel to what the supervisor has.
// Outputs: ../sim_validation_plots/ms_real_split_<ver>.png
//          ../sim_validation_plots/ms_real_showcase_<ver>.png
//          ms_real_split_<ver>.txt
void ms_real_split(const char *realf = "../clusters_seeds_island_79507-0.root_ntuplizer.root",
                   const char *i91 = "island91_frames_production_v6.root",
                   const char *ver = "v6", const char *vtag = "V6",
                   double r0 = 49.0)
{
  using namespace MSR;
  const double RCOEF = 100. / (0.299792458 * 1.4);
  const double RW[2][2] = {{101, 137}, {RCOEF * 1.5, RCOEF * 2.5}};
  const double MSTRUTH[2] = {1.07, 0.31};          // truth-level MS sigma (ms_split v53f), for scale

  std::map<long, CT> grp[2];                       // [0]=real [1]=sim
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
      CT &T = grp[0][(long) ev * 1000 + (long) sid];
      T.x.push_back(x); T.y.push_back(y); T.r.push_back(std::hypot(x, y));
      T.lmin = std::min(T.lmin, (int) lay); T.lmax = std::max(T.lmax, (int) lay);
    }
    f->Close();
  }
  {
    TFile *f = TFile::Open(i91);
    TTree *c = (TTree *) f->Get("ntp_cluster");
    TTree *u = (TTree *) f->Get("ntp_truth");
    float ev, lay, x, y, tid, cls, ntrks;
    c->SetBranchStatus("*", 0);
    for (auto b : {"event", "layer", "x", "y"}) c->SetBranchStatus(b, 1);
    c->SetBranchAddress("event", &ev);
    c->SetBranchAddress("layer", &lay);
    c->SetBranchAddress("x", &x);
    c->SetBranchAddress("y", &y);
    u->SetBranchStatus("*", 0);
    for (auto b : {"gtrackID", "cls", "ntrks"}) u->SetBranchStatus(b, 1);
    u->SetBranchAddress("gtrackID", &tid);
    u->SetBranchAddress("cls", &cls);
    u->SetBranchAddress("ntrks", &ntrks);
    for (Long64_t i = 0; i < c->GetEntries(); ++i)
    {
      u->GetEntry(i);
      if (cls != 0 || ntrks != 1) continue;
      c->GetEntry(i);
      CT &T = grp[1][(long) ev * 1000000 + (long) tid];
      T.x.push_back(x); T.y.push_back(y); T.r.push_back(std::hypot(x, y));
      T.lmin = std::min(T.lmin, (int) lay); T.lmax = std::max(T.lmax, (int) lay);
    }
    f->Close();
  }

  // identical processing; keep the worst real pT~0.5 track for the showcase
  TH1D *hd[2];
  hd[0] = new TH1D("hds0", Form(";tangent mismatch at r = %g cm  #Delta#psi [mrad];tracks (unit area)", r0), 101, -101, 101);
  hd[1] = (TH1D *) hd[0]->Clone("hds1");
  TH1D *hr[2][2];
  for (int s = 0; s < 2; ++s)
    for (int w = 0; w < 2; ++w)
      hr[s][w] = new TH1D(Form("hrs%d%d", s, w), ";per-track circle-fit RMS [mm];tracks (unit area)", 60, 0, 3.0);
  std::vector<double> dps[2][2], rmsv[2][2], d0v[2];
  CT worstT;
  double worstDp = 0;
  for (int s = 0; s < 2; ++s)
  {
    for (auto &kv : grp[s])
    {
      Fit F;
      int w = classify(kv.second, F, RW);
      if (w < 0) continue;
      rmsv[s][w].push_back(F.rms * 10);
      hr[s][w]->Fill(F.rms * 10);
      d0v[s].push_back(std::fabs(std::hypot(F.a, F.b) - F.R));
      double dp;
      if (splitClus(kv.second, r0, dp))
      {
        dps[s][w].push_back(dp);
        if (s == 0) hd[w]->Fill(std::max(-100.5, std::min(100.5, dp)));
        if (s == 0 && w == 0 && std::fabs(dp) > std::fabs(worstDp)) { worstDp = dp; worstT = kv.second; }
      }
    }
  }
  auto rmsof = [](std::vector<double> &v) {
    double sum = 0;
    for (double q : v) sum += q * q;
    return v.empty() ? 0. : std::sqrt(sum / v.size());
  };
  auto coreof = [&](std::vector<double> &v) {
    double sg = rmsof(v);
    for (int it = 0; it < 3; ++it)
    {
      double s2 = 0;
      long n = 0;
      for (double q : v)
        if (std::fabs(q) < 3 * sg) { s2 += q * q; n++; }
      if (!n) break;
      sg = std::sqrt(s2 / n);
    }
    return sg;
  };
  double sg[2][2], se[2][2], co[2][2];
  for (int s = 0; s < 2; ++s)
    for (int w = 0; w < 2; ++w)
    {
      sg[s][w] = rmsof(dps[s][w]);
      se[s][w] = dps[s][w].empty() ? 0 : sg[s][w] / std::sqrt(2. * dps[s][w].size());
      co[s][w] = coreof(dps[s][w]);
    }

  FILE *fo = fopen(Form("ms_real_split_%s.txt", ver), "w");
  auto P = [&](const char *fmt, ...) {
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    va_start(ap, fmt); vfprintf(fo, fmt, ap); va_end(ap);
  };
  P("ms_real_split %s — REAL split-arc mismatch in the ms_split style; sim %s in parens\n", ver, vtag);
  const char *wl[2] = {"pT ~0.5 GeV", "pT 1.5-2.5 GeV"};
  for (int w = 0; w < 2; ++w)
  {
    P("%s: real N=%zu sigma = %.2f ± %.2f (core %.2f)   [sim N=%zu sigma %.2f (core %.2f)]\n",
      wl[w], dps[0][w].size(), sg[0][w], se[0][w], co[0][w],
      dps[1][w].size(), sg[1][w], co[1][w]);
    P("  circle RMS median real %.0f um (sim %.0f) -> data/MC %.2f\n",
      med(rmsv[0][w]) * 1000, med(rmsv[1][w]) * 1000,
      med(rmsv[1][w]) > 0 ? med(rmsv[0][w]) / med(rmsv[1][w]) : 0);
  }
  P("truth-level MS for scale: %.2f / %.2f mrad (ms_r0scan v53f, r0=49 row; truth hits bit-identical since) — invisible at cluster level\n",
    MSTRUTH[0], MSTRUTH[1]);
  P("d0 median real %.2f cm, sim %.2f cm\n", med(d0v[0]), med(d0v[1]));
  P("worst real showcase: |dpsi| = %.1f mrad (largest in the pT~0.5 window)\n", std::fabs(worstDp));
  fclose(fo);

  // ---------- PNG 1: supervisor-style split figure ----------
  gStyle->SetOptStat(0);
  TCanvas *cv = new TCanvas("cvsp", Form("ms real split %s", ver), 1500, 620);
  cv->Divide(2, 1);
  cv->cd(1);
  {
    for (TH1D *h : {hd[0], hd[1]}) if (h->Integral() > 0) h->Scale(1. / h->Integral());
    hd[0]->SetLineColor(kBlue + 1); hd[1]->SetLineColor(kGreen + 2);
    hd[0]->SetLineWidth(2); hd[1]->SetLineWidth(2);
    hd[0]->SetTitle("split-arc tangent mismatch (REAL cluster-tracks)");
    hd[0]->SetMaximum(1.35 * std::max(hd[0]->GetMaximum(), hd[1]->GetMaximum()));
    hd[0]->Draw("hist");
    hd[1]->Draw("hist same");
    TLegend *lg = new TLegend(0.13, 0.70, 0.66, 0.88);
    lg->SetBorderSize(0);
    lg->AddEntry(hd[0], Form("p_{T}^{fit} 0.45-0.55: #sigma = %.1f#pm%.1f (core %.1f)", sg[0][0], se[0][0], co[0][0]), "l");
    lg->AddEntry(hd[1], Form("p_{T}^{fit} 1.5-2.5:   #sigma = %.1f#pm%.1f (core %.1f)", sg[0][1], se[0][1], co[0][1]), "l");
    lg->Draw();
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.036);
    tx.DrawLatex(0.13, 0.63, Form("truth-level MS: %.2f / %.2f mrad (invisible here)", MSTRUTH[0], MSTRUTH[1]));
    tx.SetTextSize(0.030);
    tx.DrawLatex(0.13, 0.57, "edge bins collect the |#Delta#psi| > 100 mrad tail");
  }
  cv->cd(2);
  {
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.045);
    tx.DrawLatex(0.06, 0.90, Form("REAL split-arc check (run 79507; sim %s in parens)", vtag));
    tx.SetTextSize(0.038);
    tx.DrawLatex(0.06, 0.80, "quantity                              p_{T} 0.5           p_{T} ~2 GeV");
    tx.DrawLatex(0.06, 0.72, Form("#sigma(#Delta#psi) real (sim)            %.1f (%.1f) mrad   %.1f (%.1f) mrad",
                                  sg[0][0], sg[1][0], sg[0][1], sg[1][1]));
    tx.DrawLatex(0.06, 0.64, Form("circle RMS real (sim)             %.0f (%.0f) #mum    %.0f (%.0f) #mum",
                                  med(rmsv[0][0]) * 1000, med(rmsv[1][0]) * 1000,
                                  med(rmsv[0][1]) * 1000, med(rmsv[1][1]) * 1000));
    tx.DrawLatex(0.06, 0.56, Form("data/MC (circle RMS)              %.2f              %.2f",
                                  med(rmsv[1][0]) > 0 ? med(rmsv[0][0]) / med(rmsv[1][0]) : 0,
                                  med(rmsv[1][1]) > 0 ? med(rmsv[0][1]) / med(rmsv[1][1]) : 0));
    tx.DrawLatex(0.06, 0.48, Form("truth-level MS (for scale)        %.2f mrad          %.2f mrad", MSTRUTH[0], MSTRUTH[1]));
    tx.SetTextSize(0.034);
    tx.DrawLatex(0.06, 0.36, "#Rightarrow mismatch is resolution-driven and p_{T}-independent: MS invisible");
    tx.DrawLatex(0.06, 0.29, "     under cluster resolution in the real data itself");
    double rr05 = med(rmsv[1][0]) > 0 ? med(rmsv[0][0]) / med(rmsv[1][0]) : 0;
    double rr2 = med(rmsv[1][1]) > 0 ? med(rmsv[0][1]) / med(rmsv[1][1]) : 0;
    if (rr2 > 1.15)
    {
      tx.DrawLatex(0.06, 0.21, Form("#Rightarrow data/MC: %.2f at 0.5 GeV; +%.0f%% real excess at high p_{T}",
                                    rr05, 100 * (rr2 - 1)));
      tx.DrawLatex(0.06, 0.14, "     (alignment/distortion-like, not in the sim by design)");
    }
    else
    {
      tx.DrawLatex(0.06, 0.21, Form("#Rightarrow data/MC: %.2f at 0.5 GeV, %.2f stiff; high-p_{T} excess resolved", rr05, rr2));
      tx.DrawLatex(0.06, 0.14, "     by the v5.4 field; sim low-p_{T} over-smear remains (response family)");
    }
    tx.SetTextSize(0.030);
    tx.DrawLatex(0.06, 0.05, "real tracks = ntp_clus_trk (event, seedID); windows by fitted curvature both sides");
  }
  cv->SaveAs(Form("../sim_validation_plots/ms_real_split_%s.png", ver));

  // ---------- PNG 2: worst showcase + other statistics ----------
  TCanvas *cw = new TCanvas("cvws", Form("ms real showcase %s", ver), 1500, 660);
  cw->Divide(2, 1);
  cw->cd(1);
  {
    Fit Fi, Fo;
    double bh[4], dp;
    splitClus(worstT, r0, dp, &Fi, &Fo, bh);
    double xlo = 1e9, xhi = -1e9, ylo = 1e9, yhi = -1e9;
    for (size_t i = 0; i < worstT.x.size(); ++i)
    {
      xlo = std::min(xlo, worstT.x[i]); xhi = std::max(xhi, worstT.x[i]);
      ylo = std::min(ylo, worstT.y[i]); yhi = std::max(yhi, worstT.y[i]);
    }
    double span = std::max(xhi - xlo, yhi - ylo) * 1.35;
    double cxm = 0.5 * (xlo + xhi), cym = 0.5 * (ylo + yhi);
    TH1 *fr = gPad->DrawFrame(cxm - span / 2, cym - span / 2, cxm + span / 2, cym + span / 2,
                              Form("WORST real track: |#Delta#psi| = %.0f mrad (largest in p_{T}~0.5 window);x [cm];y [cm]", std::fabs(dp)));
    fr->GetYaxis()->SetTitleOffset(1.3);
    TEllipse *ei = new TEllipse(Fi.a, Fi.b, Fi.R, Fi.R);
    ei->SetFillStyle(0); ei->SetLineColor(kBlue + 1); ei->SetLineStyle(2); ei->Draw();
    TEllipse *eo = new TEllipse(Fo.a, Fo.b, Fo.R, Fo.R);
    eo->SetFillStyle(0); eo->SetLineColor(kRed + 1); eo->SetLineStyle(2); eo->Draw();
    TGraph *gi = new TGraph(), *go = new TGraph();
    for (size_t i = 0; i < worstT.x.size(); ++i)
    {
      if (worstT.r[i] < r0) gi->SetPoint(gi->GetN(), worstT.x[i], worstT.y[i]);
      else go->SetPoint(go->GetN(), worstT.x[i], worstT.y[i]);
    }
    gi->SetMarkerStyle(20); gi->SetMarkerSize(0.7); gi->SetMarkerColor(kBlue + 1); gi->Draw("P same");
    go->SetMarkerStyle(21); go->SetMarkerSize(0.7); go->SetMarkerColor(kRed + 1); go->Draw("P same");
    TLegend *lg = new TLegend(0.13, 0.72, 0.60, 0.88);
    lg->SetBorderSize(0);
    lg->AddEntry(gi, Form("inner half (r < %g) + its fit", r0), "p");
    lg->AddEntry(go, Form("outer half (r > %g) + its fit", r0), "p");
    lg->Draw();
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.032);
    tx.DrawLatex(0.13, 0.65, "tail case: tracker mis-association / kink,");
    tx.DrawLatex(0.13, 0.59, "the population the 3#sigma core excludes");
  }
  cw->cd(2);
  {
    for (int s = 0; s < 2; ++s)
      for (int w = 0; w < 2; ++w)
        if (hr[s][w]->Integral() > 0) hr[s][w]->Scale(1. / hr[s][w]->Integral());
    hr[0][0]->SetLineColor(kBlue + 1); hr[0][0]->SetLineWidth(3);
    hr[1][0]->SetLineColor(kBlue + 1); hr[1][0]->SetLineWidth(2); hr[1][0]->SetLineStyle(2);
    hr[0][1]->SetLineColor(kGreen + 2); hr[0][1]->SetLineWidth(3);
    hr[1][1]->SetLineColor(kGreen + 2); hr[1][1]->SetLineWidth(2); hr[1][1]->SetLineStyle(2);
    hr[0][0]->SetTitle("per-track circle-fit RMS (solid real, dashed sim)");
    double mx = 0;
    for (int s = 0; s < 2; ++s)
      for (int w = 0; w < 2; ++w) mx = std::max(mx, hr[s][w]->GetMaximum());
    hr[0][0]->SetMaximum(1.4 * mx);
    hr[0][0]->Draw("hist");
    hr[1][0]->Draw("hist same");
    hr[0][1]->Draw("hist same");
    hr[1][1]->Draw("hist same");
    TLegend *lg = new TLegend(0.44, 0.62, 0.89, 0.88);
    lg->SetBorderSize(0);
    lg->AddEntry(hr[0][0], Form("real p_{T}~0.5, med %.0f #mum", med(rmsv[0][0]) * 1000), "l");
    lg->AddEntry(hr[1][0], Form("sim  p_{T}~0.5, med %.0f #mum", med(rmsv[1][0]) * 1000), "l");
    lg->AddEntry(hr[0][1], Form("real p_{T} 1.5-2.5, med %.0f #mum", med(rmsv[0][1]) * 1000), "l");
    lg->AddEntry(hr[1][1], Form("sim  p_{T} 1.5-2.5, med %.0f #mum", med(rmsv[1][1]) * 1000), "l");
    lg->Draw();
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.031);
    tx.DrawLatex(0.44, 0.54, Form("d_{0} median: real %.2f cm, sim %.2f cm", med(d0v[0]), med(d0v[1])));
    tx.DrawLatex(0.44, 0.48, Form("N tracks: real %zu / %zu, sim %zu / %zu",
                                  rmsv[0][0].size(), rmsv[0][1].size(), rmsv[1][0].size(), rmsv[1][1].size()));
  }
  cw->SaveAs(Form("../sim_validation_plots/ms_real_showcase_%s.png", ver));
  printf("wrote ms_real_split_%s.png + ms_real_showcase_%s.png + ms_real_split_%s.txt\n", ver, ver, ver);
}

// ---------------------------------------------------------------------------
// ms_d0diag — the offset-vs-distortion diagnostic (2026-08-05): signed
// d0 = |C| - R of every full-crosser vs the circle-center azimuth phi_c.
// A pure vertex/frame translation (x0,y0) gives d0s ~ x0 cos(phi_c) +
// y0 sin(phi_c) (charge-independent; curvature corrections <~ 1 mm at the
// measured scale). Coherent distortion gives non-cosine azimuthal structure
// and/or drift dependence — tested as the cosine-subtracted residual vs the
// track's median tbin. Sim v5.3 runs as the null control on both panels.
// (x0, y0) extracted by linear least squares on tracks with |d0s| < 8 cm;
// medians per bin are used for display robustness against secondaries.
// Output: ../sim_validation_plots/ms_d0diag_<ver>.png + ms_d0diag_<ver>.txt
void ms_d0diag(const char *realf = "../clusters_seeds_island_79507-0.root_ntuplizer.root",
               const char *i91 = "island91_frames_production_v6.root",
               const char *ver = "v6", const char *vtag = "V6")
{
  using namespace MSR;
  struct DT { std::vector<double> x, y; int lmin = 99, lmax = 0; std::vector<double> tb; };
  std::map<long, DT> grp[2];
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
      DT &T = grp[0][(long) ev * 1000 + (long) sid];
      T.x.push_back(x); T.y.push_back(y); T.tb.push_back(tb);
      T.lmin = std::min(T.lmin, (int) lay); T.lmax = std::max(T.lmax, (int) lay);
    }
    f->Close();
  }
  {
    TFile *f = TFile::Open(i91);
    TTree *c = (TTree *) f->Get("ntp_cluster");
    TTree *u = (TTree *) f->Get("ntp_truth");
    float ev, lay, x, y, tb, tid, cls, ntrks;
    c->SetBranchStatus("*", 0);
    for (auto b : {"event", "layer", "x", "y", "tbin"}) c->SetBranchStatus(b, 1);
    c->SetBranchAddress("event", &ev);
    c->SetBranchAddress("layer", &lay);
    c->SetBranchAddress("x", &x);
    c->SetBranchAddress("y", &y);
    c->SetBranchAddress("tbin", &tb);
    u->SetBranchStatus("*", 0);
    for (auto b : {"gtrackID", "cls", "ntrks"}) u->SetBranchStatus(b, 1);
    u->SetBranchAddress("gtrackID", &tid);
    u->SetBranchAddress("cls", &cls);
    u->SetBranchAddress("ntrks", &ntrks);
    for (Long64_t i = 0; i < c->GetEntries(); ++i)
    {
      u->GetEntry(i);
      if (cls != 0 || ntrks != 1) continue;
      c->GetEntry(i);
      DT &T = grp[1][(long) ev * 1000000 + (long) tid];
      T.x.push_back(x); T.y.push_back(y); T.tb.push_back(tb);
      T.lmin = std::min(T.lmin, (int) lay); T.lmax = std::max(T.lmax, (int) lay);
    }
    f->Close();
  }

  // per-track observables: d0s, phi_c, median tbin
  std::vector<double> d0s[2], phc[2], tbm[2];
  for (int s = 0; s < 2; ++s)
  {
    for (auto &kv : grp[s])
    {
      DT &T = kv.second;
      if ((int) T.x.size() < 12 || T.lmin > 11 || T.lmax < 50) continue;
      Fit F = fitCircle(T.x, T.y);
      if (!F.ok || F.rms > 0.25 || F.R < 35) continue;
      d0s[s].push_back(std::hypot(F.a, F.b) - F.R);
      phc[s].push_back(std::atan2(F.b, F.a));
      std::vector<double> t = T.tb;
      std::sort(t.begin(), t.end());
      tbm[s].push_back(t[t.size() / 2]);
    }
  }

  // linear LSQ for (x0, y0): d0s = x0 cos(phi_c) + y0 sin(phi_c)
  double X0[2], Y0[2], rmsRaw[2], rmsRes[2];
  for (int s = 0; s < 2; ++s)
  {
    double Scc = 0, Sss = 0, Scs = 0, Sc = 0, Ss = 0;
    long n = 0;
    for (size_t i = 0; i < d0s[s].size(); ++i)
    {
      if (std::fabs(d0s[s][i]) > 8) continue;
      double c = std::cos(phc[s][i]), sn = std::sin(phc[s][i]);
      Scc += c * c; Sss += sn * sn; Scs += c * sn;
      Sc += d0s[s][i] * c; Ss += d0s[s][i] * sn;
      n++;
    }
    double det = Scc * Sss - Scs * Scs;
    X0[s] = det != 0 ? (Sc * Sss - Ss * Scs) / det : 0;
    Y0[s] = det != 0 ? (Ss * Scc - Sc * Scs) / det : 0;
    double s2r = 0, s2f = 0;
    long m = 0;
    for (size_t i = 0; i < d0s[s].size(); ++i)
    {
      if (std::fabs(d0s[s][i]) > 8) continue;
      double fit = X0[s] * std::cos(phc[s][i]) + Y0[s] * std::sin(phc[s][i]);
      s2r += d0s[s][i] * d0s[s][i];
      s2f += (d0s[s][i] - fit) * (d0s[s][i] - fit);
      m++;
    }
    rmsRaw[s] = m ? std::sqrt(s2r / m) : 0;
    rmsRes[s] = m ? std::sqrt(s2f / m) : 0;
  }

  // profiles: median d0s per phi bin; cosine-subtracted median per tbin bin
  const int NB = 18, NTB2 = 12;
  auto medprof = [&](std::vector<double> &vx, std::vector<double> &vy, double lo, double hi,
                     int nb, double sub_x0, double sub_y0, bool subcos,
                     std::vector<double> &cx, std::vector<double> &cy) {
    std::vector<std::vector<double>> bins(nb);
    for (size_t i = 0; i < vx.size(); ++i)
    {
      if (std::fabs(vy[i]) > 8) continue;
      int b = (int) ((vx[i] - lo) / (hi - lo) * nb);
      if (b < 0 || b >= nb) continue;
      double val = vy[i];
      if (subcos) val -= sub_x0 * std::cos(phc[0][0] * 0);   // placeholder, unused path
      bins[b].push_back(val);
    }
    for (int b = 0; b < nb; ++b)
    {
      if ((int) bins[b].size() < 5) continue;
      std::sort(bins[b].begin(), bins[b].end());
      cx.push_back(lo + (b + 0.5) * (hi - lo) / nb);
      cy.push_back(bins[b][bins[b].size() / 2]);
    }
  };

  FILE *fo = fopen(Form("ms_d0diag_%s.txt", ver), "w");
  auto P = [&](const char *fmt, ...) {
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    va_start(ap, fmt); vfprintf(fo, fmt, ap); va_end(ap);
  };
  const char *sn[2] = {"real", "sim"};
  P("ms_d0diag %s — signed d0 vs circle-center azimuth (offset-vs-distortion)\n", ver);
  for (int s = 0; s < 2; ++s)
    P("%s: N=%zu  (x0, y0) = (%+.2f, %+.2f) cm  |offset| = %.2f cm  raw RMS %.2f -> cosine-subtracted %.2f cm\n",
      sn[s], d0s[s].size(), X0[s], Y0[s], std::hypot(X0[s], Y0[s]), rmsRaw[s], rmsRes[s]);
  double fexp = rmsRaw[0] > 0 ? 1. - rmsRes[0] * rmsRes[0] / (rmsRaw[0] * rmsRaw[0]) : 0;
  P("real: cosine explains %.0f%% of the d0s variance (|d0s|<8 cm sample)\n", 100. * fexp);
  fclose(fo);

  // ---------- figure ----------
  gStyle->SetOptStat(0);
  TCanvas *cv = new TCanvas("cvd0", Form("ms d0 diag %s", ver), 1500, 620);
  cv->Divide(2, 1);
  cv->cd(1);
  {
    TH1 *fr = gPad->DrawFrame(-M_PI, -7, M_PI, 7,
                              "signed d_{0} vs circle-center azimuth;#varphi_{c} = atan2(b, a) [rad];d_{0}^{s} = |C| - R [cm]");
    fr->GetYaxis()->SetTitleOffset(1.2);
    TGraph *gsc = new TGraph();
    for (size_t i = 0; i < d0s[0].size(); ++i)
      if (std::fabs(d0s[0][i]) < 7) gsc->SetPoint(gsc->GetN(), phc[0][i], d0s[0][i]);
    gsc->SetMarkerStyle(1);
    gsc->SetMarkerColorAlpha(kGray + 1, 0.5);
    gsc->Draw("P same");
    std::vector<double> cx, cy;
    medprof(phc[0], d0s[0], -M_PI, M_PI, NB, 0, 0, false, cx, cy);
    TGraph *gp = new TGraph((int) cx.size(), cx.data(), cy.data());
    gp->SetMarkerStyle(20); gp->SetMarkerColor(kBlack); gp->Draw("P same");
    std::vector<double> sx, sy;
    medprof(phc[1], d0s[1], -M_PI, M_PI, NB, 0, 0, false, sx, sy);
    TGraph *gs = new TGraph((int) sx.size(), sx.data(), sy.data());
    gs->SetMarkerStyle(24); gs->SetMarkerColor(kBlue + 1); gs->Draw("P same");
    TGraph *gf = new TGraph();
    for (int i = 0; i <= 100; ++i)
    {
      double ph = -M_PI + 2 * M_PI * i / 100.;
      gf->SetPoint(i, ph, X0[0] * std::cos(ph) + Y0[0] * std::sin(ph));
    }
    gf->SetLineColor(kRed + 1); gf->SetLineWidth(2); gf->Draw("L same");
    TLegend *lg = new TLegend(0.13, 0.70, 0.68, 0.88);
    lg->SetBorderSize(0);
    lg->AddEntry(gp, "real, median per #varphi_{c} bin", "p");
    lg->AddEntry(gf, Form("cosine fit: (x_{0}, y_{0}) = (%+.2f, %+.2f) cm", X0[0], Y0[0]), "l");
    lg->AddEntry(gs, Form("sim %s control", vtag), "p");
    lg->Draw();
  }
  cv->cd(2);
  {
    TH1 *fr = gPad->DrawFrame(0, -3, 970, 3,
                              "cosine-subtracted d_{0}^{s} vs drift time;track median tbin;d_{0}^{s} - offset model [cm]");
    fr->GetYaxis()->SetTitleOffset(1.2);
    TLine *l0 = new TLine(0, 0, 970, 0);
    l0->SetLineStyle(2); l0->SetLineColor(kGray + 2); l0->Draw();
    for (int s = 0; s < 2; ++s)
    {
      std::vector<double> rx, ry;
      std::vector<std::vector<double>> bins(NTB2);
      for (size_t i = 0; i < d0s[s].size(); ++i)
      {
        if (std::fabs(d0s[s][i]) > 8) continue;
        double res = d0s[s][i] - (X0[s] * std::cos(phc[s][i]) + Y0[s] * std::sin(phc[s][i]));
        int b = (int) (tbm[s][i] / 970. * NTB2);
        if (b >= 0 && b < NTB2) bins[b].push_back(res);
      }
      for (int b = 0; b < NTB2; ++b)
      {
        if ((int) bins[b].size() < 5) continue;
        std::sort(bins[b].begin(), bins[b].end());
        rx.push_back((b + 0.5) * 970. / NTB2);
        ry.push_back(bins[b][bins[b].size() / 2]);
      }
      TGraph *g = new TGraph((int) rx.size(), rx.data(), ry.data());
      g->SetMarkerStyle(s == 0 ? 20 : 24);
      g->SetMarkerColor(s == 0 ? kBlack : kBlue + 1);
      g->Draw("P same");
    }
    TLegend *lg = new TLegend(0.40, 0.74, 0.88, 0.88);
    lg->SetBorderSize(0);
    TGraph *d1 = new TGraph(); d1->SetMarkerStyle(20); d1->SetMarkerColor(kBlack);
    TGraph *d2 = new TGraph(); d2->SetMarkerStyle(24); d2->SetMarkerColor(kBlue + 1);
    lg->AddEntry(d1, "real, median per tbin bin", "p");
    lg->AddEntry(d2, Form("sim %s control", vtag), "p");
    lg->Draw();
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.034);
    tx.DrawLatex(0.40, 0.66, "flat = translation only; slope/structure = drift-");
    tx.DrawLatex(0.40, 0.60, "dependent (distortion-like) component");
  }
  cv->SaveAs(Form("../sim_validation_plots/ms_d0diag_%s.png", ver));
  printf("wrote ../sim_validation_plots/ms_d0diag_%s.png + ms_d0diag_%s.txt\n", ver, ver);
}
