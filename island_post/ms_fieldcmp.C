// ms_fieldcmp.C — SIM-ONLY ideal-field vs distortion-field comparison
// (user, 2026-08-20): the detached chain rerun from the G4 hits with NO
// distortion field (pp_pipeline_ideal.sh: FIELD="" at transport, everything
// else byte-identical to V6 — same G4 chunks, composer seeds, envelope,
// no-flash, tbin<=960, same export call) fitted against the V6 field-on
// production. No real data anywhere.
//   fc_clusters : island91 pair, truth-grouped (event,gtrackID) cluster
//                 tracks — whole-track circle RMS in the two fitted-pT
//                 windows, split-arc tangent mismatch at r0=49, signed d0.
//   fc_pixels   : digi pair, per-pixel-truth-grouped raw-pixel tracks —
//                 whole-track (global) fit and 4-adjacent-row local
//                 (short-sagitta) fit.
// Fitter/bars identical to the ms_nofinder battery (Kasa + 6 Gauss-Newton;
// global bar >=12 pts, span >=15 cm, 45<=R_fit<2e4; windows by fitted R).
// Truth hits themselves are the ideal trajectory by construction (20 um /
// 1.365 mrad) and are identical in both productions — quoted, not re-run.
// Outputs: ../sim_validation_plots/ms_fieldcmp_{clusters,pixels}_<ver>.png
//          ms_fieldcmp_<ver>.txt
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <cstdio>
#include <cstdarg>
#include <cmath>
#include <map>
#include <set>
#include <vector>
#include <algorithm>

namespace MFC
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

struct Grp { std::vector<double> x, y, r; };
bool fitBar(const Grp &G, Fit &F)
{
  if ((int) G.x.size() < 12) return false;
  double rlo = 1e9, rhi = 0;
  for (double r : G.r) { rlo = std::min(rlo, r); rhi = std::max(rhi, r); }
  if (rhi - rlo < 15) return false;
  F = fitCircle(G.x, G.y);
  return F.ok && F.R >= 45 && F.R < 2e4;
}
void stats(const std::vector<double> &v, double &sig, double &err, double &core)
{
  sig = err = core = 0;
  if (v.size() < 2) return;
  double mu = 0;
  for (double q : v) mu += q;
  mu /= v.size();
  double s2 = 0;
  for (double q : v) s2 += (q - mu) * (q - mu);
  sig = std::sqrt(s2 / v.size());
  err = sig / std::sqrt(2. * v.size());
  double sg = sig;
  for (int it = 0; it < 3; ++it)
  {
    double t2 = 0; long n = 0;
    for (double q : v)
      if (std::fabs(q - mu) < 3 * sg) { t2 += (q - mu) * (q - mu); n++; }
    if (!n) break;
    sg = std::sqrt(t2 / n);
  }
  core = sg;
}
bool loadRows(double rowR[55])
{
  FILE *fp = fopen("tpc_geom_table.txt", "r");
  if (!fp) { printf("no tpc_geom_table.txt\n"); return false; }
  char line[512];
  while (fgets(line, sizeof line, fp))
  {
    int L, nb; double r, sl, p0, p1;
    if (sscanf(line, "%d %d %lf %lf %lf %lf", &L, &nb, &r, &sl, &p0, &p1) == 6 && L >= 7 && L <= 54)
      rowR[L] = r;
  }
  fclose(fp);
  return true;
}
FILE *openLedger(const char *ver) { return fopen(Form("ms_fieldcmp_%s.txt", ver), "a"); }
const int CIDE = kBlue + 1, CDIS = kRed + 1;    // ideal solid blue, field-on dashed red
}  // namespace MFC

// ---------------------------------------------------------------------------
// CLUSTER LEVEL: island91 ideal vs field-on, truth-grouped tracks.
void fc_clusters(const char *ideal = "island91_frames_production_v6ideal.root",
                 const char *dist = "island91_frames_production_v6.root",
                 const char *ver = "v6")
{
  using namespace MFC;
  const double RCOEF = 100. / (0.299792458 * 1.4);
  const double RW[2][2] = {{101, 137}, {RCOEF * 1.5, RCOEF * 2.5}};
  const double R0 = 49.;
  const char *fn[2] = {ideal, dist};
  const char *sn[2] = {"ideal", "field"};
  std::vector<double> rmsw[2][2], dps[2][2], d0v[2];
  long ngrp[2] = {0, 0};
  for (int s = 0; s < 2; ++s)
  {
    TFile *f = TFile::Open(fn[s]);
    if (!f || f->IsZombie()) { printf("missing %s\n", fn[s]); return; }
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
      Grp &G = kv.second;
      Fit F;
      if (!fitBar(G, F)) continue;
      ngrp[s]++;
      int w = -1;
      for (int q = 0; q < 2; ++q)
        if (F.R > RW[q][0] && F.R < RW[q][1]) w = q;
      double d0s = std::hypot(F.a, F.b) - F.R;
      if (std::fabs(d0s) < 8) d0v[s].push_back(d0s);
      if (w < 0) continue;
      rmsw[s][w].push_back(F.rms * 1e4);
      // split-arc on full crossers
      double rmin = 1e9, rmax = 0;
      for (double r : G.r) { rmin = std::min(rmin, r); rmax = std::max(rmax, r); }
      if ((int) G.x.size() < 30 || rmin > 35 || rmax < 72) continue;
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
      double dx = hxo - hxi, dy = hyo - hyi, psi_i, psi_o;
      if (!tangentAtR(Fi, R0, hxi, hyi, dx, dy, psi_i)) continue;
      if (!tangentAtR(Fo, R0, hxo, hyo, dx, dy, psi_o)) continue;
      dps[s][w].push_back(wrapphi(psi_o - psi_i) * 1e3);
    }
    printf("fc_clusters %s: %ld tracks\n", sn[s], ngrp[s]);
  }
  double sg[2][2], se[2][2], sc[2][2], d0rms[2], d0e, d0c;
  for (int s = 0; s < 2; ++s)
  {
    for (int w = 0; w < 2; ++w) stats(dps[s][w], sg[s][w], se[s][w], sc[s][w]);
    stats(d0v[s], d0rms[s], d0e, d0c);
  }
  auto qsub = [](double a, double b) { return a > b ? std::sqrt(a * a - b * b) : 0.; };
  FILE *fo = openLedger(ver);
  auto P = [&](const char *fmt, ...) {
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    va_start(ap, fmt); vfprintf(fo, fmt, ap); va_end(ap);
  };
  P("[fc_clusters %s] SIM-ONLY ideal vs distortion field, island91 truth-grouped tracks\n", ver);
  P("  files: ideal=%s | field=%s (only delta: FIELD at transport)\n", ideal, dist);
  for (int s = 0; s < 2; ++s)
    P("  %s: %ld tracks | RMS med 0.5GeV %.0f um / stiff %.0f um | sigma(dpsi) %.1f+-%.1f / %.1f+-%.1f mrad "
      "| d0s: med|d0| %.2f cm, RMS %.2f cm (|d0s|<8)\n",
      sn[s], ngrp[s], med(rmsw[s][0]), med(rmsw[s][1]), sg[s][0], se[s][0], sg[s][1], se[s][1],
      [&]{ std::vector<double> a; for (double q : d0v[s]) a.push_back(std::fabs(q)); return med(a); }(),
      d0rms[s]);
  P("  FIELD SHARE (quadrature, dist (-) ideal): RMS 0.5GeV %.0f um | stiff %.0f um | dpsi 0.5GeV %.1f mrad\n",
    qsub(med(rmsw[1][0]), med(rmsw[0][0])), qsub(med(rmsw[1][1]), med(rmsw[0][1])),
    qsub(sg[1][0], sg[0][0]));
  P("  truth-hit level (identical in both by construction): 20 um circle RMS, 1.365 mrad MS\n");
  fclose(fo);
  gStyle->SetOptStat(0);
  TCanvas *cv = new TCanvas("cvfc", "fieldcmp clusters", 1500, 1150);
  cv->Divide(2, 2);
  const char *pt[2] = {"p_{T}^{fit} 0.45-0.55 GeV", "p_{T}^{fit} 1.5-2.5 GeV"};
  for (int w = 0; w < 2; ++w)
  {
    cv->cd(w + 1);
    TH1D *h[2];
    for (int s = 0; s < 2; ++s)
    {
      h[s] = new TH1D(Form("fcr%d%d", w, s), ";per-track circle-fit RMS [mm];tracks (unit area)", 50, 0, 2.0);
      for (double q : rmsw[s][w]) h[s]->Fill(std::min(q / 1000., 1.99));
      if (h[s]->Integral() > 0) h[s]->Scale(1. / h[s]->Integral());
    }
    h[0]->SetLineColor(CIDE); h[0]->SetLineWidth(2);
    h[1]->SetLineColor(CDIS); h[1]->SetLineWidth(2); h[1]->SetLineStyle(2);
    h[0]->SetTitle(Form("clusters, whole-track fit: %s", pt[w]));
    h[0]->SetMaximum(1.4 * std::max(h[0]->GetMaximum(), h[1]->GetMaximum()));
    h[0]->Draw("hist"); h[1]->Draw("hist same");
    TLegend *lg = new TLegend(0.44, 0.72, 0.89, 0.87);
    lg->SetBorderSize(0);
    lg->AddEntry(h[0], Form("ideal (no field), med %.0f #mum", med(rmsw[0][w])), "l");
    lg->AddEntry(h[1], Form("distortion field on, med %.0f #mum", med(rmsw[1][w])), "l");
    lg->Draw();
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.032);
    tx.DrawLatex(0.46, 0.65, Form("field share: %.0f #mum in quadrature",
                                  med(rmsw[1][w]) > med(rmsw[0][w])
                                  ? std::sqrt(std::pow(med(rmsw[1][w]), 2) - std::pow(med(rmsw[0][w]), 2)) : 0.));
  }
  cv->cd(3);
  {
    TH1D *h[2];
    for (int s = 0; s < 2; ++s)
    {
      h[s] = new TH1D(Form("fcd%d", s), ";tangent mismatch at r = 49 cm  #Delta#psi [mrad];tracks (unit area)", 50, -100, 100);
      for (double q : dps[s][0]) h[s]->Fill(std::max(-99.9, std::min(q, 99.9)));
      if (h[s]->Integral() > 0) h[s]->Scale(1. / h[s]->Integral());
    }
    h[0]->SetLineColor(CIDE); h[0]->SetLineWidth(2);
    h[1]->SetLineColor(CDIS); h[1]->SetLineWidth(2); h[1]->SetLineStyle(2);
    h[0]->SetTitle("split-arc mismatch, p_{T}~0.5 window");
    h[0]->SetMaximum(1.4 * std::max(h[0]->GetMaximum(), h[1]->GetMaximum()));
    h[0]->Draw("hist"); h[1]->Draw("hist same");
    TLegend *lg = new TLegend(0.40, 0.72, 0.89, 0.87);
    lg->SetBorderSize(0);
    lg->AddEntry(h[0], Form("ideal: #sigma = %.1f#pm%.1f mrad", sg[0][0], se[0][0]), "l");
    lg->AddEntry(h[1], Form("field on: #sigma = %.1f#pm%.1f mrad", sg[1][0], se[1][0]), "l");
    lg->Draw();
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.030);
    tx.DrawLatex(0.42, 0.65, Form("stiff window: %.1f vs %.1f mrad", sg[0][1], sg[1][1]));
    tx.DrawLatex(0.42, 0.59, "truth MS underneath: 1.07 mrad (identical both)");
  }
  cv->cd(4);
  {
    TH1D *h[2];
    for (int s = 0; s < 2; ++s)
    {
      h[s] = new TH1D(Form("fc0%d", s), ";signed d_{0} = |C| - R  [cm];tracks (unit area)", 60, -8, 8);
      for (double q : d0v[s]) h[s]->Fill(q);
      if (h[s]->Integral() > 0) h[s]->Scale(1. / h[s]->Integral());
    }
    h[0]->SetLineColor(CIDE); h[0]->SetLineWidth(2);
    h[1]->SetLineColor(CDIS); h[1]->SetLineWidth(2); h[1]->SetLineStyle(2);
    h[0]->SetTitle("impact parameter of the fitted circle");
    h[0]->SetMaximum(1.4 * std::max(h[0]->GetMaximum(), h[1]->GetMaximum()));
    h[0]->Draw("hist"); h[1]->Draw("hist same");
    TLegend *lg = new TLegend(0.55, 0.72, 0.89, 0.87);
    lg->SetBorderSize(0);
    lg->AddEntry(h[0], Form("ideal: RMS %.2f cm", d0rms[0]), "l");
    lg->AddEntry(h[1], Form("field on: RMS %.2f cm", d0rms[1]), "l");
    lg->Draw();
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.030);
    tx.DrawLatex(0.14, 0.80, "the d_{0} broadening IS the field:");
    tx.DrawLatex(0.14, 0.74, "trajectories point at the beamline;");
    tx.DrawLatex(0.14, 0.68, "their displaced images do not");
  }
  cv->SaveAs(Form("../sim_validation_plots/ms_fieldcmp_clusters_%s.png", ver));
  printf("wrote ms_fieldcmp_clusters_%s.png\n", ver);
}

// ---------------------------------------------------------------------------
// PIXEL LEVEL: digi ideal vs field-on, per-pixel-truth-grouped tracks;
// whole-track (global) vs 4-adjacent-row local (short-sagitta) fits.
void fc_pixels(const char *ideal = "digi_frames_production_v6ideal.root",
               const char *dist = "digi_frames_production_v6.root",
               int nsim = 60, const char *ver = "v6")
{
  using namespace MFC;
  double rowR[55];
  if (!loadRows(rowR)) return;
  auto nearRow = [&](double r) -> int {
    int best = -1; double bd = 1e9;
    for (int L = 7; L <= 54; ++L)
    {
      double d = std::fabs(r - rowR[L]);
      if (d < bd) { bd = d; best = L; }
    }
    return bd < 0.60 ? best : -1;
  };
  const char *fn[2] = {ideal, dist};
  const char *sn[2] = {"ideal", "field"};
  std::vector<double> grms[2], wrms[2];
  long ngrp[2] = {0, 0}, nfitw[2] = {0, 0}, nwin[2] = {0, 0};
  for (int s = 0; s < 2; ++s)
  {
    TFile *f = TFile::Open(fn[s]);
    if (!f || f->IsZombie()) { printf("missing %s\n", fn[s]); return; }
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
    for (auto &kv : g)
    {
      Grp &G = kv.second;
      Fit F;
      if (!fitBar(G, F)) continue;
      ngrp[s]++;
      grms[s].push_back(F.rms * 1e4);
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
        nwin[s]++;
        if ((int) wrow[w].size() < 3 || (int) wx[w].size() < 5) continue;
        Fit L = fitCircle(wx[w], wy[w]);
        if (!L.ok) continue;
        nfitw[s]++;
        wrms[s].push_back(L.rms * 1e4);
      }
    }
    printf("fc_pixels %s: %ld tracks\n", sn[s], ngrp[s]);
  }
  auto qsub = [](double a, double b) { return a > b ? std::sqrt(a * a - b * b) : 0.; };
  FILE *fo = openLedger(ver);
  auto P = [&](const char *fmt, ...) {
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    va_start(ap, fmt); vfprintf(fo, fmt, ap); va_end(ap);
  };
  P("[fc_pixels %s] SIM-ONLY ideal vs distortion field, digi raw-pixel tracks (%d frames)\n", ver, nsim);
  for (int s = 0; s < 2; ++s)
    P("  %s: %ld tracks | GLOBAL RMS med %.0f um | LOCAL 4-row med %.0f um (%ld/%ld windows)\n",
      sn[s], ngrp[s], med(grms[s]), med(wrms[s]), nfitw[s], nwin[s]);
  P("  FIELD SHARE: global %.0f um in quadrature | local shift %.0f um "
    "(short-sagitta fit rejects the smooth field by construction)\n",
    qsub(med(grms[1]), med(grms[0])), med(wrms[1]) - med(wrms[0]));
  fclose(fo);
  gStyle->SetOptStat(0);
  TCanvas *cv = new TCanvas("cvfp", "fieldcmp pixels", 1500, 620);
  cv->Divide(2, 1);
  const char *pt[2] = {"GLOBAL whole-track fit", "LOCAL 4-row short-sagitta fit"};
  for (int p = 0; p < 2; ++p)
  {
    cv->cd(p + 1);
    std::vector<double> *v[2] = {p == 0 ? &grms[0] : &wrms[0], p == 0 ? &grms[1] : &wrms[1]};
    double xhi = p == 0 ? 5000 : 4000;
    TH1D *h[2];
    for (int s = 0; s < 2; ++s)
    {
      h[s] = new TH1D(Form("fp%d%d", p, s), ";per-fit circle RMS [#mum];fits (unit area)", 60, 0, xhi);
      for (double q : *v[s]) h[s]->Fill(std::min(q, xhi - 1));
      if (h[s]->Integral() > 0) h[s]->Scale(1. / h[s]->Integral());
    }
    h[0]->SetLineColor(CIDE); h[0]->SetLineWidth(2);
    h[1]->SetLineColor(CDIS); h[1]->SetLineWidth(2); h[1]->SetLineStyle(2);
    h[0]->SetTitle(Form("raw pixels, truth-grouped: %s", pt[p]));
    h[0]->SetMaximum(1.4 * std::max(h[0]->GetMaximum(), h[1]->GetMaximum()));
    h[0]->Draw("hist"); h[1]->Draw("hist same");
    TLegend *lg = new TLegend(0.42, 0.72, 0.89, 0.87);
    lg->SetBorderSize(0);
    lg->AddEntry(h[0], Form("ideal (no field), med %.0f #mum", med(*v[0])), "l");
    lg->AddEntry(h[1], Form("distortion field on, med %.0f #mum", med(*v[1])), "l");
    lg->Draw();
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.030);
    if (p == 0)
      tx.DrawLatex(0.44, 0.65, Form("field share: %.0f #mum in quadrature",
                                    qsub(med(grms[1]), med(grms[0]))));
    else
    {
      tx.DrawLatex(0.44, 0.65, Form("shift: %.0f #mum", med(wrms[1]) - med(wrms[0])));
      tx.DrawLatex(0.44, 0.59, "local fit is field-blind by construction");
    }
  }
  cv->SaveAs(Form("../sim_validation_plots/ms_fieldcmp_pixels_%s.png", ver));
  printf("wrote ms_fieldcmp_pixels_%s.png\n", ver);
}
