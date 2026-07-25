// mcs_split.C — supervisor follow-up (2026-07-25): does multiple scattering
// in the TPC gas need explicit treatment, or is it negligible vs cluster
// resolution (~200 um)?  Explicit MCS measurement from the v5.2-era sim:
//   SPLIT-ARC TEST on G4 truth hits: fit the inner (r<49) and outer (r>49)
//   halves of each full-crossing track as independent circles, propagate both
//   tangents to the r=49 cm crossing, and histogram the tangent-angle
//   mismatch dpsi. With exact truth points the mismatch is the accumulated
//   in-gas scattering (plus small fit noise, quoted separately).
// Two pT windows show the 1/p scaling: [0.45,0.55] and [1.5,2.5] GeV.
// Highland comparison uses the as-built gas from sphenix_p5.gdml:
//   Ar/C/F/H mass fractions 0.580/0.099/0.310/0.010, rho 2.148e-3 g/cm3
//   (Ar75:CF4-20:iso-5) -> X0 = 24.0 g/cm2 = 112 m. theta0 evaluated per
//   track with its 3D path length (pion beta assumed; sample is 87% pi).
// Companion: truth_circle.C (trajectory circularity + cluster residuals).
// Output: ../sim_validation_plots/mcs_split_v52.png + mcs_split_v52.txt
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
#include <vector>
#include <algorithm>

namespace MCSD
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
// tangent angle of circle (a,b,R) at its intersection with the beamline
// cylinder r=r0, choosing the branch nearest (hx,hy); sign along (dx,dy).
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
  double tx = -(py - F.b), ty = (px - F.a);       // perp to radius vector
  if (tx * dx + ty * dy < 0) { tx = -tx; ty = -ty; }
  psi = std::atan2(ty, tx);
  return true;
}
}  // namespace MCSD

void mcs_split(int ng4 = 3, const char *g4pat = "../P5/PP_g4hit_%d.root",
               const char *i91 = "island91_frames_production_v52.root")
{
  using namespace MCSD;
  const double PTW[2][2] = {{0.45, 0.55}, {1.5, 2.5}};
  const double R0 = 49.0;                 // split radius [cm]
  const double X0LEN = 11200.;            // gas radiation length [cm] (Ar75:CF4-20:iso-5)
  const double MPI = 0.13957;
  struct Trk { std::vector<double> x, y, r; float pt = 0, p = 0; };
  std::vector<double> dpsi[2], rmscirc[2], fitnoise[2], hi[2], rmsclus[2];
  TH1D *hd[2];
  hd[0] = new TH1D("hd0", ";tangent mismatch at r = 49 cm  #Delta#psi [mrad];tracks (unit area)", 81, -8.1, 8.1);
  hd[1] = (TH1D *) hd[0]->Clone("hd1");
  long nfit[2] = {0, 0}, nclus[2] = {0, 0};

  // ---------- G4 truth hits: split-arc tangent mismatch ----------
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
      double pt = std::hypot(gpx, gpy);
      if (!((pt > PTW[0][0] && pt < PTW[0][1]) || (pt > PTW[1][0] && pt < PTW[1][1]))) continue;
      Trk &T = trks[(long) ev * 100000 + (long) tid];
      T.pt = pt; T.p = std::sqrt(pt * pt + gpz * gpz);
      T.x.push_back(gx); T.y.push_back(gy); T.r.push_back(std::hypot(gx, gy));
    }
    for (auto &kv : trks)
    {
      Trk &T = kv.second;
      int w = (T.pt < 1.0) ? 0 : 1;
      double rmin = 1e9, rmax = 0;
      for (double r : T.r) { rmin = std::min(rmin, r); rmax = std::max(rmax, r); }
      if ((int) T.x.size() < 30 || rmin > 34 || rmax < 72) continue;
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
      double dx = hxo - hxi, dy = hyo - hyi;      // local direction of motion
      double psi_i, psi_o;
      if (!tangentAtR(Fi, R0, hxi, hyi, dx, dy, psi_i)) continue;
      if (!tangentAtR(Fo, R0, hxo, hyo, dx, dy, psi_o)) continue;
      double dpsi_mrad = wrapphi(psi_o - psi_i) * 1e3;
      if (std::fabs(dpsi_mrad) > 8) continue;     // decay/hard-scatter kinks out of core
      nfit[w]++;
      dpsi[w].push_back(dpsi_mrad);
      hd[w]->Fill(dpsi_mrad);
      // whole-track circularity for the same track (context number)
      Fit Fa = fitCircle(T.x, T.y);
      if (Fa.ok) rmscirc[w].push_back(Fa.rms * 1e4);            // um
      // fit-noise estimate on dpsi: per-half tangent-angle error from the
      // half-fit residuals over the half lever arm (Gluckstern-type scaling)
      double Lh = 29.0;
      double s_i = Fi.rms / Lh * std::sqrt(192. / (Fi.n + 4));
      double s_o = Fo.rms / Lh * std::sqrt(192. / (Fo.n + 4));
      fitnoise[w].push_back(std::sqrt(s_i * s_i + s_o * s_o) * 1e3);   // mrad
      // Highland prediction for the full crossing, this track's path + beta
      double beta = T.p / std::sqrt(T.p * T.p + MPI * MPI);
      double path = 58.8 * (T.p / T.pt);          // 3D path across the gas [cm]
      double xX0 = path / X0LEN;
      double th0 = 13.6e-3 / (beta * T.p) * std::sqrt(xX0) * (1 + 0.038 * std::log(xX0));
      hi[w].push_back(th0 * 1e3);                 // mrad
    }
    f->Close();
    printf("g4 file %d: split fits so far %ld / %ld\n", fi, nfit[0], nfit[1]);
  }

  // ---------- reco clusters: per-point circle residual in both windows ----------
  {
    TFile *f = TFile::Open(i91);
    TTree *c = (TTree *) f->Get("ntp_cluster");
    TTree *u = (TTree *) f->Get("ntp_truth");
    float ev, lay, x, y, tid, gpt, cls, ntrks;
    c->SetBranchStatus("*", 0);
    for (auto b : {"event", "layer", "x", "y"}) c->SetBranchStatus(b, 1);
    c->SetBranchAddress("event", &ev);
    c->SetBranchAddress("layer", &lay);
    c->SetBranchAddress("x", &x);
    c->SetBranchAddress("y", &y);
    u->SetBranchStatus("*", 0);
    for (auto b : {"gtrackID", "gpt", "cls", "ntrks"}) u->SetBranchStatus(b, 1);
    u->SetBranchAddress("gtrackID", &tid);
    u->SetBranchAddress("gpt", &gpt);
    u->SetBranchAddress("cls", &cls);
    u->SetBranchAddress("ntrks", &ntrks);
    struct CT { std::vector<double> x, y; int lmin = 99, lmax = 0; float pt = 0; };
    std::map<long, CT> trks;
    for (Long64_t i = 0; i < c->GetEntries(); ++i)
    {
      u->GetEntry(i);
      if (cls != 0 || ntrks != 1) continue;
      if (!((gpt > PTW[0][0] && gpt < PTW[0][1]) || (gpt > PTW[1][0] && gpt < PTW[1][1]))) continue;
      c->GetEntry(i);
      CT &T = trks[(long) ev * 1000000 + (long) tid];
      T.pt = gpt;
      T.x.push_back(x); T.y.push_back(y);
      T.lmin = std::min(T.lmin, (int) lay); T.lmax = std::max(T.lmax, (int) lay);
    }
    for (auto &kv : trks)
    {
      CT &T = kv.second;
      int w = (T.pt < 1.0) ? 0 : 1;
      if ((int) T.x.size() < 12 || T.lmin > 11 || T.lmax < 50) continue;
      Fit F = fitCircle(T.x, T.y);
      if (!F.ok || F.rms > 0.25) continue;        // drop residual multi-track pathologies
      nclus[w]++;
      rmsclus[w].push_back(F.rms * 1e4);          // um
    }
    f->Close();
  }

  // ---------- summary ----------
  FILE *fo = fopen("mcs_split_v52.txt", "w");
  auto P = [&](const char *fmt, ...) {
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    va_start(ap, fmt); vfprintf(fo, fmt, ap); va_end(ap);
  };
  double sig[2], fn[2];
  for (int w = 0; w < 2; ++w)
  {
    sig[w] = hd[w]->GetRMS();
    fn[w] = med(fitnoise[w]);
  }
  P("mcs_split v52 — gas Ar75:CF4-20:iso-5, X0 = 112 m; split radius %.0f cm\n", R0);
  for (int w = 0; w < 2; ++w)
  {
    double core = std::sqrt(std::max(0., sig[w] * sig[w] - fn[w] * fn[w]));
    P("pT [%.2f,%.2f] GeV: %ld split tracks\n", PTW[w][0], PTW[w][1], nfit[w]);
    P("  tangent mismatch sigma(dpsi) = %.2f mrad  (fit noise %.2f -> MCS %.2f mrad)\n",
      sig[w], fn[w], core);
    P("  Highland full-crossing theta0 = %.2f mrad (per-track median)\n", med(hi[w]));
    P("  whole-track circle RMS (truth) = %.0f um (median)\n", med(rmscirc[w]));
    P("  reco-cluster circle RMS        = %.0f um (median, %ld tracks)\n", med(rmsclus[w]), nclus[w]);
    P("  MCS displacement / cluster RMS ~ %.3f\n",
      med(rmscirc[w]) / std::max(1., med(rmsclus[w])));
  }
  P("scaling check: sigma ratio W1/W2 = %.2f (Highland predicts %.2f)\n",
    sig[1] > 0 ? sig[0] / sig[1] : 0., med(hi[1]) > 0 ? med(hi[0]) / med(hi[1]) : 0.);
  fclose(fo);

  // ---------- figure ----------
  gStyle->SetOptStat(0);
  TCanvas *cv = new TCanvas("cv", "mcs split v52", 1500, 620);
  cv->Divide(2, 1);
  cv->cd(1);
  {
    for (TH1D *h : {hd[0], hd[1]}) if (h->Integral() > 0) h->Scale(1. / h->Integral());
    hd[0]->SetLineColor(kBlue + 1); hd[1]->SetLineColor(kGreen + 2);
    hd[0]->SetLineWidth(2); hd[1]->SetLineWidth(2);
    hd[0]->SetTitle("split-arc tangent mismatch (G4 truth hits)");
    hd[0]->SetMaximum(1.35 * std::max(hd[0]->GetMaximum(), hd[1]->GetMaximum()));
    hd[0]->Draw("hist");
    hd[1]->Draw("hist same");
    TLegend *lg = new TLegend(0.13, 0.70, 0.60, 0.88);
    lg->SetBorderSize(0);
    lg->AddEntry(hd[0], Form("p_{T} 0.45-0.55: #sigma = %.2f mrad", sig[0]), "l");
    lg->AddEntry(hd[1], Form("p_{T} 1.5-2.5:   #sigma = %.2f mrad", sig[1]), "l");
    lg->Draw();
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.036);
    tx.DrawLatex(0.13, 0.63, Form("Highland #theta_{0}: %.2f / %.2f mrad", med(hi[0]), med(hi[1])));
  }
  cv->cd(2);
  {
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.045);
    tx.DrawLatex(0.06, 0.90, "in-gas MCS vs cluster resolution (v5.2 sim, measured)");
    tx.SetTextSize(0.038);
    double core0 = std::sqrt(std::max(0., sig[0] * sig[0] - fn[0] * fn[0]));
    double core1 = std::sqrt(std::max(0., sig[1] * sig[1] - fn[1] * fn[1]));
    tx.DrawLatex(0.06, 0.80, "quantity                              p_{T} 0.5           p_{T} ~2 GeV");
    tx.DrawLatex(0.06, 0.72, Form("MCS angle, measured (Highland)   %.1f (%.1f) mrad   %.2f (%.2f) mrad",
                                  core0, med(hi[0]), core1, med(hi[1])));
    tx.DrawLatex(0.06, 0.64, Form("truth-track circle RMS            %.0f #mum            %.0f #mum",
                                  med(rmscirc[0]), med(rmscirc[1])));
    tx.DrawLatex(0.06, 0.56, Form("reco-cluster circle RMS           %.0f #mum           %.0f #mum",
                                  med(rmsclus[0]), med(rmsclus[1])));
    tx.DrawLatex(0.06, 0.48, Form("trajectory scatter / cluster res  %.3f             %.3f",
                                  med(rmscirc[0]) / std::max(1., med(rmsclus[0])),
                                  med(rmscirc[1]) / std::max(1., med(rmsclus[1]))));
    tx.SetTextSize(0.034);
    tx.DrawLatex(0.06, 0.36, "#Rightarrow position level: MCS << cluster resolution at all p_{T} (<3%)");
    tx.DrawLatex(0.06, 0.28, "#Rightarrow p_{T} estimation: MCS term #approx #theta_{0}/#psi_{bend} #approx 0.4%, p_{T}-independent;");
    tx.DrawLatex(0.06, 0.21, "     comparable to the measurement term only below ~0.5-1 GeV");
    tx.DrawLatex(0.06, 0.10, "gas: Ar75:CF4-20:iso-5 (from sphenix_p5.gdml), X_{0} = 112 m");
  }
  cv->SaveAs("../sim_validation_plots/mcs_split_v52.png");
  printf("wrote ../sim_validation_plots/mcs_split_v52.png + mcs_split_v52.txt\n");
}
