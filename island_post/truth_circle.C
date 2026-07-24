// truth_circle.C — supervisor question (2026-07-24): take the truth cluster
// coordinates of a ~500 MeV-pT track — how close is the trajectory to a circle?
// Three tiers of "truth", B=1.4 T -> R = pT/(0.3B) ~ 119 cm at 0.5 GeV:
//   T1 truth steps    : ntp_g4hit step midpoints (gx,gy) of primary tracks
//                       with vertex pT in [0.45,0.55] (P5/PP_g4hit_*.root)
//   T2 truth clusters : per-pad-row gedep-weighted centroids of the T1 steps
//                       (row = nearest tpc_geom_table radius, 0.7 cm gate —
//                        same assignment rule as tpc_transport)
//   T3 tagged reco    : island91 v5.1 production ntp_cluster (x,y) with
//                       row-aligned ntp_truth cls==0 && ntrks==1
// T3 CAVEAT (found by this study): frame_composer keys truth ids by
// (library FILE, raw per-collision trk) — frame_composer.C:70/101/442/484 —
// and ttab by first-in-file id (:279), so one ntp_truth id bundles ~several
// unrelated trajectories and carries a file-frozen representative gpt.
// T3 therefore splits each id group into single-trajectory ARCS by phi
// chaining and fits per arc; arcs with fitted R in [101,137] cm are the
// ~0.5 GeV subset (the gpt window itself selects scrambled labels).
// Fit: algebraic (Kasa) init + 6 Gauss-Newton geometric iterations.
// Quality: full R1->R3 crossers only (r<34 & r>72 cm; layers <=11 & >=50).
// Showcase panels use the MOST-SAMPLED qualifying track/arc (stated criterion).
// Output: ../sim_validation_plots/truth_circle_v51.png + truth_circle_v51.txt
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TEllipse.h>
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

namespace TCIRC
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
  for (int it = 0; it < 6; ++it)   // geometric refinement
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
double pct(std::vector<double> v, double p)
{
  if (v.empty()) return 0;
  std::sort(v.begin(), v.end());
  return v[std::min(v.size() - 1, (size_t)(p * v.size()))];
}
double wrapphi(double d)
{
  while (d > M_PI) d -= 2 * M_PI;
  while (d < -M_PI) d += 2 * M_PI;
  return d;
}
}  // namespace TCIRC

void truth_circle(double pt_lo = 0.45, double pt_hi = 0.55, int ng4 = 3,
                  const char *g4pat = "../P5/PP_g4hit_%d.root",
                  const char *i91 = "island91_frames_production_v51.root")
{
  using namespace TCIRC;
  const double BFIELD = 1.4;                                   // T (project currency)
  const double RCOEF = 100. / (0.299792458 * BFIELD);          // R[cm] = RCOEF * pT[GeV]
  const double RSEL_LO = 101, RSEL_HI = 137;                   // fitted-R window ~ pT 0.45-0.55
  double geoR[55]; int ngeo = 0;
  {
    FILE *ft = fopen("tpc_geom_table.txt", "r");
    if (!ft) { printf("no tpc_geom_table.txt (run from island_post/)\n"); return; }
    char line[512];
    while (fgets(line, sizeof line, ft))
    {
      int L, nb; double r, sl, p0, p1;
      if (line[0] == '#') continue;
      if (sscanf(line, "%d %d %lf %lf %lf %lf", &L, &nb, &r, &sl, &p0, &p1) == 6 && L >= 7 && L <= 54)
      { geoR[L] = r; ngeo++; }
    }
    fclose(ft);
    printf("geom table: %d rows\n", ngeo);
  }

  struct Trk { std::vector<double> x, y, r, w; float pt = 0; int flav = 0; };
  std::vector<double> rms1, rms2, rms3, rms3sel, rfit1, rfit3, rrat1;
  TH1D *hres1 = new TH1D("hres1", ";signed residual to fitted circle [mm];points (unit area)", 121, -3.025, 3.025);
  TH1D *hres2 = (TH1D *) hres1->Clone("hres2");
  TH1D *hres3 = (TH1D *) hres1->Clone("hres3");
  TH1D *hrms1 = new TH1D("hrms1", ";per-track residual RMS [mm];tracks (unit area)", 60, 0, 3.0);
  TH1D *hrms2 = (TH1D *) hrms1->Clone("hrms2");
  TH1D *hrms3 = (TH1D *) hrms1->Clone("hrms3");
  std::map<int, int> flavcnt;
  long nwin1 = 0, nfull1 = 0;
  Trk show1, show2;  // T1 showcase steps + its T2 centroids
  Fit showF1, showF2;
  std::vector<double> show3x, show3y, show3r;
  Fit showF3;
  double show1pt = 0;

  // ---------- T1 + T2: truth steps from single pp collisions ----------
  for (int fi = 0; fi < ng4; ++fi)
  {
    TFile *f = TFile::Open(Form(g4pat, fi));
    if (!f || f->IsZombie()) { printf("missing %s\n", Form(g4pat, fi)); continue; }
    TTree *t = (TTree *) f->Get("ntp_g4hit");
    float ev, gx, gy, gpx, gpy, gedep, tid, flav;
    t->SetBranchStatus("*", 0);
    for (auto b : {"event", "gx", "gy", "gpx", "gpy", "gedep", "gtrackID", "gflavor"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("event", &ev);
    t->SetBranchAddress("gx", &gx);
    t->SetBranchAddress("gy", &gy);
    t->SetBranchAddress("gpx", &gpx);
    t->SetBranchAddress("gpy", &gpy);
    t->SetBranchAddress("gedep", &gedep);
    t->SetBranchAddress("gtrackID", &tid);
    t->SetBranchAddress("gflavor", &flav);
    std::map<long, Trk> trks;
    for (Long64_t i = 0; i < t->GetEntries(); ++i)
    {
      t->GetEntry(i);
      if (tid <= 0) continue;                       // primaries only (secondaries carry -id)
      double pt = std::hypot(gpx, gpy);             // vertex momentum, constant per track
      if (pt < pt_lo || pt > pt_hi) continue;
      Trk &T = trks[(long) ev * 100000 + (long) tid];
      T.pt = pt; T.flav = (int) flav;
      T.x.push_back(gx); T.y.push_back(gy);
      T.r.push_back(std::hypot(gx, gy)); T.w.push_back(gedep);
    }
    for (auto &kv : trks)
    {
      Trk &T = kv.second;
      nwin1++;
      double rmin = 1e9, rmax = 0;
      for (double r : T.r) { rmin = std::min(rmin, r); rmax = std::max(rmax, r); }
      if ((int) T.x.size() < 25 || rmin > 34 || rmax < 72) continue;   // full crossers only
      Fit F = fitCircle(T.x, T.y);
      if (!F.ok) continue;
      nfull1++;
      flavcnt[std::abs(T.flav)]++;
      rms1.push_back(F.rms * 10);
      rfit1.push_back(F.R);
      rrat1.push_back(F.R / (RCOEF * T.pt));
      hrms1->Fill(F.rms * 10);
      for (size_t i = 0; i < T.x.size(); ++i)
        hres1->Fill((std::hypot(T.x[i] - F.a, T.y[i] - F.b) - F.R) * 10);
      // T2: gedep-weighted per-pad-row centroids (transport's 0.7 cm gate)
      double cx[55] = {0}, cy[55] = {0}, cw[55] = {0};
      for (size_t i = 0; i < T.x.size(); ++i)
      {
        int Lb = -1; double dbest = 0.7;
        for (int q = 7; q <= 54; ++q)
        {
          double d = std::fabs(T.r[i] - geoR[q]);
          if (d < dbest) { dbest = d; Lb = q; }
        }
        if (Lb < 0 || T.w[i] <= 0) continue;
        cx[Lb] += T.x[i] * T.w[i]; cy[Lb] += T.y[i] * T.w[i]; cw[Lb] += T.w[i];
      }
      Trk C;
      for (int q = 7; q <= 54; ++q)
        if (cw[q] > 0) { C.x.push_back(cx[q] / cw[q]); C.y.push_back(cy[q] / cw[q]); C.r.push_back(geoR[q]); }
      if ((int) C.x.size() >= 12)
      {
        Fit F2 = fitCircle(C.x, C.y);
        if (F2.ok)
        {
          rms2.push_back(F2.rms * 10);
          hrms2->Fill(F2.rms * 10);
          for (size_t i = 0; i < C.x.size(); ++i)
            hres2->Fill((std::hypot(C.x[i] - F2.a, C.y[i] - F2.b) - F2.R) * 10);
          // showcase: most-sampled CLEAN track (RMS<40 um, R within 2% of pT/(0.3B))
          if (F.rms < 0.004 && std::fabs(F.R / (RCOEF * T.pt) - 1) < 0.02 &&
              T.x.size() > show1.x.size()) { show1 = T; showF1 = F; show2 = C; showF2 = F2; show1pt = T.pt; }
        }
      }
    }
    f->Close();
    printf("T1 file %d done: %ld window primaries so far, %ld full crossers\n", fi, nwin1, nfull1);
  }

  // ---------- T3: truth-tagged reco clusters, split into trajectory arcs ----------
  long ngrp = 0, ngrp_multi = 0, narcs_tot = 0, nfit3 = 0, nfit3sel = 0;
  double wclus_tot = 0, wclus_multi = 0;
  {
    TFile *f = TFile::Open(i91);
    TTree *c = (TTree *) f->Get("ntp_cluster");
    TTree *u = (TTree *) f->Get("ntp_truth");
    if (c->GetEntries() != u->GetEntries()) { printf("ALIGN FAIL %lld vs %lld\n", c->GetEntries(), u->GetEntries()); return; }
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
    struct CL { float lay, x, y, phi; };
    std::map<long, std::vector<CL>> grp;   // (frame, conflated id) -> clusters
    for (Long64_t i = 0; i < c->GetEntries(); ++i)
    {
      u->GetEntry(i);
      if (cls != 0 || ntrks != 1) continue;
      if (gpt < pt_lo || gpt > pt_hi) continue;   // NB selects the class REPRESENTATIVE's pt
      c->GetEntry(i);
      grp[(long) ev * 1000000 + (long) tid].push_back({lay, x, y, (float) std::atan2(y, x)});
    }
    for (auto &kv : grp)
    {
      auto &v = kv.second;
      ngrp++;
      wclus_tot += v.size();
      std::sort(v.begin(), v.end(), [](const CL &a, const CL &b) { return a.lay < b.lay || (a.lay == b.lay && a.phi < b.phi); });
      // phi-chaining: one chain per trajectory arc
      struct Chain { std::vector<int> idx; double lL = 0, lP = 0, pL = 0, pP = 0; int np = 0; };
      std::vector<Chain> chains;
      for (int i = 0; i < (int) v.size(); ++i)
      {
        int best = -1; double bd = 0.06;
        for (int k = 0; k < (int) chains.size(); ++k)
        {
          Chain &ch = chains[k];
          double pred = ch.lP;
          if (ch.np >= 2 && ch.lL > ch.pL)
            pred += (ch.lP - ch.pP) / (ch.lL - ch.pL) * (v[i].lay - ch.lL);
          double d = std::fabs(wrapphi(v[i].phi - pred));
          if (d < bd) { bd = d; best = k; }
        }
        if (best < 0) { chains.push_back({}); best = (int) chains.size() - 1; }
        Chain &ch = chains[best];
        ch.idx.push_back(i);
        if (ch.np == 0) { ch.lL = v[i].lay; ch.lP = v[i].phi; }
        else { ch.pL = ch.lL; ch.pP = ch.lP; ch.lL = v[i].lay; ch.lP = v[i].phi; }
        ch.np++;
      }
      int nsub = 0;
      for (auto &ch : chains) if (ch.np >= 4) nsub++;
      narcs_tot += nsub;
      if (nsub >= 2) { ngrp_multi++; wclus_multi += v.size(); }
      for (auto &ch : chains)
      {
        if (ch.np < 12) continue;
        int lmin = 99, lmax = 0;
        std::vector<double> X, Y;
        for (int i : ch.idx)
        {
          lmin = std::min(lmin, (int) v[i].lay); lmax = std::max(lmax, (int) v[i].lay);
          X.push_back(v[i].x); Y.push_back(v[i].y);
        }
        if (lmin > 11 || lmax < 50) continue;
        Fit F = fitCircle(X, Y);
        if (!F.ok) continue;
        nfit3++;
        rms3.push_back(F.rms * 10);
        rfit3.push_back(F.R);
        if (F.R > RSEL_LO && F.R < RSEL_HI)         // the ~0.5 GeV subset by curvature
        {
          nfit3sel++;
          rms3sel.push_back(F.rms * 10);
          hrms3->Fill(F.rms * 10);
          for (size_t i = 0; i < X.size(); ++i)
            hres3->Fill((std::hypot(X[i] - F.a, Y[i] - F.b) - F.R) * 10);
          // showcase arc: most clusters among clean fits (RMS < 0.9 mm)
          if (F.rms < 0.09 && X.size() > show3x.size())
          {
            show3x = X; show3y = Y; showF3 = F;
            show3r.clear();
            for (size_t i = 0; i < X.size(); ++i) show3r.push_back(std::hypot(X[i], Y[i]));
          }
        }
      }
    }
    f->Close();
  }

  // ---------- summary ----------
  double Rexp = RCOEF * 0.5 * (pt_lo + pt_hi);
  // sagitta seen over a radial span [r_in, r_out] for a circle through the
  // vertex: point at distance r from origin sits at turning half-angle
  // alpha = asin(r/2R); chord = 2R sin(da), sagitta = R (1 - cos(da)).
  auto sagitta = [&](double r_in, double r_out) {
    double da = std::asin(r_out / (2 * Rexp)) - std::asin(r_in / (2 * Rexp));
    return Rexp * (1. - std::cos(da));
  };
  double sag = sagitta(20.0, 78.0);              // full gas volume (T1 steps)
  double sagrow = sagitta(geoR[7], geoR[54]);    // pad rows L7-L54 (clusters)
  FILE *fo = fopen("truth_circle_v51.txt", "w");
  auto P = [&](const char *fmt, ...) {
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    va_start(ap, fmt); vfprintf(fo, fmt, ap); va_end(ap);
  };
  P("truth_circle v51 — pT window [%.2f,%.2f] GeV, B=%.2f T, R_exp=%.1f cm\n", pt_lo, pt_hi, BFIELD, Rexp);
  P("sagitta: %.2f cm over the gas (r 20-78, T1)  |  %.2f cm over pad rows (r %.1f-%.1f, T2/T3)\n",
    sag, sagrow, geoR[7], geoR[54]);
  P("T1 truth steps   : %ld window primaries, %ld full R1->R3 crossers fitted\n", nwin1, nfull1);
  P("                   median RMS %.0f um  (p90 %.0f um), median R_fit %.2f cm, R_fit/R_exp(pT) %.4f\n",
    med(rms1) * 1000, pct(rms1, 0.90) * 1000, med(rfit1), med(rrat1));
  P("T2 truth clusters: median RMS %.0f um  (p90 %.0f um)\n", med(rms2) * 1000, pct(rms2, 0.90) * 1000);
  P("T3 tagged reco   : %ld id groups (gpt-window; ids are (file,trk)-conflated) -> %ld arcs\n", ngrp, narcs_tot);
  P("                   multi-arc groups %.1f%% (cluster-weighted %.1f%%) — composer truth-id conflation\n",
    100. * ngrp_multi / std::max(1L, ngrp), 100. * wclus_multi / std::max(1., wclus_tot));
  P("                   %ld full-span arcs fitted, median R_fit %.1f cm (~inclusive-pT image);\n",
    nfit3, med(rfit3));
  P("                   R-selected [%.0f,%.0f] cm (~0.5 GeV): %ld arcs, median RMS %.0f um (p90 %.0f um)\n",
    RSEL_LO, RSEL_HI, nfit3sel, med(rms3sel) * 1000, pct(rms3sel, 0.90) * 1000);
  P("deviation/sagitta: T1 %.2e (gas)  T2 %.2e  T3 %.2e (pad rows)\n",
    med(rms1) / 10 / sag, med(rms2) / 10 / sagrow, med(rms3sel) / 10 / sagrow);
  P("species of fitted T1 tracks: ");
  for (auto &kv : flavcnt) P("|pdg|=%d:%d  ", kv.first, kv.second);
  P("\n");
  fclose(fo);

  // ---------- figure ----------
  gStyle->SetOptStat(0);
  TCanvas *cv = new TCanvas("cv", "truth circle v51", 1600, 1200);
  cv->Divide(2, 2);

  cv->cd(1);
  {
    double xlo = 1e9, xhi = -1e9, ylo = 1e9, yhi = -1e9;
    for (size_t i = 0; i < show1.x.size(); ++i)
    {
      xlo = std::min(xlo, show1.x[i]); xhi = std::max(xhi, show1.x[i]);
      ylo = std::min(ylo, show1.y[i]); yhi = std::max(yhi, show1.y[i]);
    }
    double span = std::max(xhi - xlo, yhi - ylo) * 1.25;
    double cxm = 0.5 * (xlo + xhi), cym = 0.5 * (ylo + yhi);
    TH1 *fr = gPad->DrawFrame(cxm - span / 2, cym - span / 2, cxm + span / 2, cym + span / 2,
                              Form("showcase truth track (most-sampled, pT=%.3f GeV)  R_{fit}=%.1f cm;x [cm];y [cm]", show1pt, showF1.R));
    fr->GetYaxis()->SetTitleOffset(1.3);
    TEllipse *el = new TEllipse(showF1.a, showF1.b, showF1.R, showF1.R);
    el->SetFillStyle(0); el->SetLineColor(kGray + 2); el->SetLineStyle(2); el->Draw();
    TGraph *g1 = new TGraph(show1.x.size(), show1.x.data(), show1.y.data());
    g1->SetMarkerStyle(20); g1->SetMarkerSize(0.35); g1->SetMarkerColor(kBlue + 1); g1->Draw("P same");
    TGraph *g2 = new TGraph(show2.x.size(), show2.x.data(), show2.y.data());
    g2->SetMarkerStyle(24); g2->SetMarkerSize(1.0); g2->SetMarkerColor(kRed + 1); g2->Draw("P same");
    TLegend *lg = new TLegend(0.14, 0.72, 0.60, 0.88);
    lg->SetBorderSize(0);
    lg->AddEntry(g1, Form("T1 G4 steps (RMS %.0f #mum)", showF1.rms * 1e4), "p");
    lg->AddEntry(g2, Form("T2 truth clusters (RMS %.0f #mum)", showF2.rms * 1e4), "p");
    lg->AddEntry(el, "fitted circle", "l");
    lg->Draw();
  }

  cv->cd(2);
  {
    std::vector<double> rr1, dd1, rr2, dd2, rr3, dd3;
    for (size_t i = 0; i < show1.x.size(); ++i)
    { rr1.push_back(show1.r[i]); dd1.push_back((std::hypot(show1.x[i] - showF1.a, show1.y[i] - showF1.b) - showF1.R) * 10); }
    for (size_t i = 0; i < show2.x.size(); ++i)
    { rr2.push_back(show2.r[i]); dd2.push_back((std::hypot(show2.x[i] - showF2.a, show2.y[i] - showF2.b) - showF2.R) * 10); }
    for (size_t i = 0; i < show3x.size(); ++i)
    { rr3.push_back(show3r[i]); dd3.push_back((std::hypot(show3x[i] - showF3.a, show3y[i] - showF3.b) - showF3.R) * 10); }
    double dmax = 0.05;
    for (double d : dd1) dmax = std::max(dmax, std::fabs(d));
    for (double d : dd3) dmax = std::max(dmax, std::fabs(d));
    TH1 *fr = gPad->DrawFrame(18, -dmax * 1.3, 80, dmax * 1.3,
                              "residual vs radius (showcase tracks);cluster radius r [cm];residual to fitted circle [mm]");
    fr->GetYaxis()->SetTitleOffset(1.3);
    TGraph *q1 = new TGraph(rr1.size(), rr1.data(), dd1.data());
    q1->SetMarkerStyle(20); q1->SetMarkerSize(0.4); q1->SetMarkerColor(kBlue + 1); q1->Draw("P same");
    TGraph *q2 = new TGraph(rr2.size(), rr2.data(), dd2.data());
    q2->SetMarkerStyle(24); q2->SetMarkerSize(1.0); q2->SetMarkerColor(kRed + 1); q2->Draw("P same");
    TGraph *q3 = new TGraph(rr3.size(), rr3.data(), dd3.data());
    q3->SetMarkerStyle(21); q3->SetMarkerSize(0.7); q3->SetMarkerColor(kGreen + 2); q3->Draw("P same");
    TLegend *lg = new TLegend(0.14, 0.70, 0.74, 0.88);
    lg->SetBorderSize(0);
    lg->AddEntry(q1, "T1 G4 steps", "p");
    lg->AddEntry(q2, "T2 truth clusters", "p");
    lg->AddEntry(q3, Form("T3 reco-cluster arc (independent track, R_{fit}=%.0f cm)", showF3.R), "p");
    lg->Draw();
  }

  cv->cd(3);
  {
    gPad->SetLogy();
    for (TH1D *h : {hres1, hres2, hres3}) if (h->Integral() > 0) h->Scale(1. / h->Integral());
    hres3->SetLineColor(kGreen + 2); hres1->SetLineColor(kBlue + 1); hres2->SetLineColor(kRed + 1);
    hres3->SetLineWidth(2); hres1->SetLineWidth(2); hres2->SetLineWidth(2);
    hres3->SetTitle("per-point residuals, all fitted tracks");
    hres3->SetMaximum(1.5 * std::max({hres1->GetMaximum(), hres2->GetMaximum(), hres3->GetMaximum()}));
    hres3->Draw("hist");
    hres1->Draw("hist same");
    hres2->Draw("hist same");
    TLegend *lg = new TLegend(0.56, 0.68, 0.89, 0.88);
    lg->SetBorderSize(0);
    lg->AddEntry(hres1, Form("T1 steps, RMS %.0f #mum", hres1->GetRMS() * 1000), "l");
    lg->AddEntry(hres2, Form("T2 truth clus, RMS %.0f #mum", hres2->GetRMS() * 1000), "l");
    lg->AddEntry(hres3, Form("T3 reco arcs, RMS %.0f #mum", hres3->GetRMS() * 1000), "l");
    lg->Draw();
  }

  cv->cd(4);
  {
    gPad->SetLogy();
    for (TH1D *h : {hrms1, hrms2, hrms3}) if (h->Integral() > 0) h->Scale(1. / h->Integral());
    hrms3->SetLineColor(kGreen + 2); hrms1->SetLineColor(kBlue + 1); hrms2->SetLineColor(kRed + 1);
    hrms3->SetLineWidth(2); hrms1->SetLineWidth(2); hrms2->SetLineWidth(2);
    hrms3->SetTitle("per-track circle-fit RMS");
    hrms3->SetMaximum(3.0 * std::max({hrms1->GetMaximum(), hrms2->GetMaximum(), hrms3->GetMaximum()}));
    hrms3->SetMinimum(2e-4);
    hrms3->Draw("hist");
    hrms1->Draw("hist same");
    hrms2->Draw("hist same");
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.031);
    tx.DrawLatex(0.30, 0.84, Form("R_{exp}(0.5 GeV) = %.1f cm; sagitta %.2f (gas) / %.2f (rows) cm", Rexp, sag, sagrow));
    tx.DrawLatex(0.30, 0.79, Form("median RMS: T1 %.0f / T2 %.0f / T3 %.0f #mum",
                                  med(rms1) * 1000, med(rms2) * 1000, med(rms3sel) * 1000));
    tx.DrawLatex(0.30, 0.74, Form("T1 R_{fit}/R_{exp}(p_{T}) = %.4f;  T1 fits %ld, T3 arcs %ld", med(rrat1), nfull1, nfit3sel));
    tx.DrawLatex(0.30, 0.69, Form("T3 groups: %.0f%% multi-arc (truth-id conflation)", 100. * ngrp_multi / std::max(1L, ngrp)));
    TLegend *lg = new TLegend(0.30, 0.46, 0.70, 0.64);
    lg->SetBorderSize(0);
    lg->AddEntry(hrms1, "T1 G4 truth steps", "l");
    lg->AddEntry(hrms2, "T2 truth clusters (per pad row)", "l");
    lg->AddEntry(hrms3, "T3 reco-cluster arcs, R-selected", "l");
    lg->Draw();
  }

  cv->SaveAs("../sim_validation_plots/truth_circle_v51.png");
  printf("wrote ../sim_validation_plots/truth_circle_v51.png + truth_circle_v51.txt\n");
}
