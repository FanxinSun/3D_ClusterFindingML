// twist_probe.C — measure the coherent r-phi TWIST profile behind the real-only
// split-half curvature offset (nonrms_probe: Dsagitta med +6.5 mm, charge-
// independent), as a per-(side, pad-row) mean azimuthal displacement profile,
// and CLOSE it: inject the real-minus-sim profile into the sim pixels and
// check that sim reproduces the real split-half offset. The profile file
// twist_profile_<ver>.txt is the handover payload (twist_field_request.md).
// Definitions: per track (nf_digipix groups, trimmed global fit), for every
// kept pixel dphi = wrap(phi_pixel - phi_fit(r_pixel)) with phi_fit the
// azimuth of the fitted circle at that radius (branch nearest the pixel);
// D(rphi) = r*dphi [um]. Profile g(side,row) = mean over tracks and pixels.
// This is the FIT-ORTHOGONAL residual part of the field (what a circle fit
// cannot absorb) — injecting it as-is reproduces the observables.
// usage: root -l -b -q 'twist_probe.C+()'
// out: twist_probe_<ver>.txt, twist_profile_<ver>.txt,
//      ../sim_validation_plots/twist_probe_<ver>.png
#include "ms_nofinder.C"
#include <TROOT.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TLine.h>
#include <unordered_map>
#include <cstdint>

namespace TWP
{
struct GrpX { std::vector<double> x, y, r; std::vector<uint8_t> side, lay; };
double medv(std::vector<double> v)
{
  if (v.empty()) return 0;
  std::sort(v.begin(), v.end());
  return v[v.size() / 2];
}
bool trimFit(const std::vector<double> &X, const std::vector<double> &Y, MNF::Fit &F, std::vector<int> &kept)
{
  kept.resize(X.size());
  for (size_t k = 0; k < kept.size(); ++k) kept[k] = (int) k;
  for (int it = 0; it < 3; ++it)
  {
    std::vector<double> xs, ys;
    for (int i : kept) { xs.push_back(X[i]); ys.push_back(Y[i]); }
    F = MNF::fitCircle(xs, ys);
    if (!F.ok) return false;
    std::vector<double> res;
    for (int i : kept) res.push_back(std::hypot(X[i] - F.a, Y[i] - F.b) - F.R);
    std::vector<double> tmp = res;
    double md = medv(tmp);
    for (double &q : tmp) q = std::fabs(q - md);
    double thr = std::max(3 * 1.4826 * medv(tmp), 0.05);
    std::vector<int> nk;
    for (size_t k = 0; k < kept.size(); ++k)
      if (std::fabs(res[k] - md) <= thr) nk.push_back(kept[k]);
    if (nk.size() == kept.size()) break;
    if ((int) nk.size() < 12) break;
    kept.swap(nk);
  }
  std::vector<double> xs, ys;
  for (int i : kept) { xs.push_back(X[i]); ys.push_back(Y[i]); }
  F = MNF::fitCircle(xs, ys);
  return F.ok;
}
// azimuth of the fitted circle at radius r0, branch nearest (hx,hy)
bool phiAtR(const MNF::Fit &F, double r0, double hx, double hy, double &phi)
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
  phi = std::atan2(py, px);
  return true;
}
double signedK(const std::vector<double> &X, const std::vector<double> &Y,
               const std::vector<double> &Rr, const MNF::Fit &F)
{
  size_t i0 = 0, i1 = 0, i2 = 0;
  double rlo = 1e9, rhi = -1;
  for (size_t i = 0; i < Rr.size(); ++i)
  { if (Rr[i] < rlo) { rlo = Rr[i]; i0 = i; } if (Rr[i] > rhi) { rhi = Rr[i]; i2 = i; } }
  double rmid = 0.5 * (rlo + rhi), bd = 1e9;
  for (size_t i = 0; i < Rr.size(); ++i)
    if (std::fabs(Rr[i] - rmid) < bd) { bd = std::fabs(Rr[i] - rmid); i1 = i; }
  double cx = (X[i1] - X[i0]) * (Y[i2] - Y[i0]) - (Y[i1] - Y[i0]) * (X[i2] - X[i0]);
  return (cx >= 0 ? 1. : -1.) / F.R;
}
}  // namespace TWP

void twist_probe(const char *digif = "digi_frames_production_v6.root", int nsim = 60,
                 const char *realf = "../clusters_seeds_island_79507-0.root_ntuplizer.root",
                 const char *ver = "v6")
{
  using namespace TWP;
  gROOT->SetBatch(1);
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
    { double d = std::fabs(r - rowR[L]); if (d < bd) { bd = d; best = L; } }
    return bd < 0.60 ? best : -1;
  };

  // ---- groups (nominal road real / truth sim; + side + row) ----------------
  std::vector<GrpX> gr[2];
  {
    struct SC { float x, y, tb; int seed; };
    std::vector<SC> sc;
    std::unordered_map<int, std::vector<int>> bucket;
    std::map<std::pair<int, int>, int> seedIdx;
    TFile *f = TFile::Open(realf);
    TTree *t = (TTree *) f->Get("ntp_clus_trk");
    float ev, sid, lay, x, y, tb;
    t->SetBranchStatus("*", 0);
    for (auto b : {"event", "seedID", "layer", "x", "y", "tbin"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("event", &ev); t->SetBranchAddress("seedID", &sid);
    t->SetBranchAddress("layer", &lay); t->SetBranchAddress("x", &x);
    t->SetBranchAddress("y", &y); t->SetBranchAddress("tbin", &tb);
    for (Long64_t i = 0; i < t->GetEntries(); ++i)
    {
      t->GetEntry(i);
      if ((int) ev == 44) continue;
      if (lay < 7 || lay > 54) continue;
      auto key = std::make_pair((int) ev, (int) sid);
      auto it = seedIdx.find(key);
      if (it == seedIdx.end()) it = seedIdx.insert({key, (int) seedIdx.size()}).first;
      bucket[(int) ev * 100 + (int) lay].push_back((int) sc.size());
      sc.push_back({x, y, tb, it->second});
    }
    f->Close();
    gr[0].assign(seedIdx.size(), GrpX());
    std::vector<std::set<int>> hid(seedIdx.size());
    f = TFile::Open(realf);
    t = (TTree *) f->Get("ntp_hit");
    float adc, hitID, zel;
    t->SetBranchStatus("*", 0);
    for (auto b : {"event", "layer", "x", "y", "tbin", "adc", "hitID", "zelem"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("event", &ev); t->SetBranchAddress("layer", &lay);
    t->SetBranchAddress("x", &x); t->SetBranchAddress("y", &y);
    t->SetBranchAddress("tbin", &tb); t->SetBranchAddress("adc", &adc);
    t->SetBranchAddress("hitID", &hitID); t->SetBranchAddress("zelem", &zel);
    for (Long64_t i = 0; i < t->GetEntries(); ++i)
    {
      t->GetEntry(i);
      if ((int) ev == 44) continue;
      if (lay < 7 || lay > 54 || adc <= 0) continue;
      auto bit = bucket.find((int) ev * 100 + (int) lay);
      if (bit == bucket.end()) continue;
      for (int j : bit->second)
      {
        const SC &c = sc[j];
        if (std::fabs(tb - c.tb) > 6) continue;
        double dx = x - c.x, dy = y - c.y;
        if (dx * dx + dy * dy > 1.2 * 1.2) continue;
        if (!hid[c.seed].insert((int) hitID).second) break;
        GrpX &G = gr[0][c.seed];
        G.x.push_back(x); G.y.push_back(y); G.r.push_back(std::hypot(x, y));
        G.side.push_back((uint8_t) zel); G.lay.push_back((uint8_t) lay);
        break;
      }
    }
    f->Close();
    printf("real: %zu seeds grouped\n", gr[0].size());
  }
  {
    TFile *f = TFile::Open(digif);
    TTree *t = (TTree *) f->Get("ntp_hit");
    float ev, lay, phi, adc, tid, zel;
    t->SetBranchStatus("*", 0);
    for (auto b : {"event", "layer", "phi", "adc", "gtrackID", "zelem"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("event", &ev); t->SetBranchAddress("layer", &lay);
    t->SetBranchAddress("phi", &phi); t->SetBranchAddress("adc", &adc);
    t->SetBranchAddress("gtrackID", &tid); t->SetBranchAddress("zelem", &zel);
    std::map<std::pair<int, int>, int> gidx;
    for (Long64_t i = 0; i < t->GetEntries(); ++i)
    {
      t->GetEntry(i);
      if ((int) ev >= nsim) continue;
      if (lay < 7 || lay > 54 || adc <= 0 || tid <= 0) continue;
      auto key = std::make_pair((int) ev, (int) tid);
      auto it = gidx.find(key);
      if (it == gidx.end()) { it = gidx.insert({key, (int) gr[1].size()}).first; gr[1].push_back(GrpX()); }
      GrpX &G = gr[1][it->second];
      double r = rowR[(int) lay];
      G.x.push_back(r * std::cos(phi)); G.y.push_back(r * std::sin(phi)); G.r.push_back(r);
      G.side.push_back((uint8_t) zel); G.lay.push_back((uint8_t) lay);
    }
    f->Close();
    printf("sim: %zu truth groups\n", gr[1].size());
  }

  // ---- measurement engine: profile + split-half, on (optionally displaced) pixels ----
  // disp[side][L] in um of r-phi displacement added to sim pixels (closure test)
  auto measure = [&](int s, const double disp[2][55], double prof[2][55], double profSec[2][12][55],
                     std::vector<double> dsag[2], std::vector<double> &dsagAll, long &ntr) {
    double sum[2][55] = {{0}}, cnt[2][55] = {{0}};
    double sumS[2][12][55] = {{{0}}}, cntS[2][12][55] = {{{0}}};
    ntr = 0;
    for (auto &G : gr[s])
    {
      std::vector<double> X = G.x, Y = G.y;
      if (disp)
        for (size_t i = 0; i < X.size(); ++i)
        {
          int L = G.lay[i]; double r = G.r[i];
          double dphi = disp[G.side[i]][L] * 1e-4 / r;
          double ph = std::atan2(Y[i], X[i]) + dphi;
          X[i] = r * std::cos(ph); Y[i] = r * std::sin(ph);
        }
      MNF::Grp B; B.x = X; B.y = Y; B.r = G.r;
      MNF::Fit F0;
      if (!MNF::fitBar(B, F0)) continue;
      MNF::Fit FT; std::vector<int> kept;
      if (!trimFit(X, Y, FT, kept)) continue;
      ntr++;
      for (int i : kept)
      {
        int L = G.lay[i]; if (L < 7 || L > 54) continue;
        double phf;
        if (!phiAtR(FT, G.r[i], X[i], Y[i], phf)) continue;
        double dphi = MNF::wrapphi(std::atan2(Y[i], X[i]) - phf);
        double d = dphi * G.r[i] * 1e4;               // um
        if (std::fabs(d) > 8000) continue;            // guard: not a member of this track
        int sd = G.side[i];
        int sec = ((int) std::floor((std::atan2(Y[i], X[i]) + 2 * M_PI) / (M_PI / 6))) % 12;
        sum[sd][L] += d; cnt[sd][L]++;
        sumS[sd][sec][L] += d; cntS[sd][sec][L]++;
      }
      // split-half Dsagitta (same as nonrms_probe)
      std::vector<double> Xi, Yi, Ri, Xo, Yo, Ro; int nsd = 0;
      for (int i : kept)
      {
        int L = G.lay[i]; nsd += G.side[i];
        if (L <= 30) { Xi.push_back(X[i]); Yi.push_back(Y[i]); Ri.push_back(G.r[i]); }
        else { Xo.push_back(X[i]); Yo.push_back(Y[i]); Ro.push_back(G.r[i]); }
      }
      if ((int) Xi.size() >= 8 && (int) Xo.size() >= 8)
      {
        double rli = 1e9, rhi2 = 0;
        for (double q : Ri) rli = std::min(rli, q);
        for (double q : Ro) rhi2 = std::max(rhi2, q);
        MNF::Fit Fi = MNF::fitCircle(Xi, Yi), Fo2 = MNF::fitCircle(Xo, Yo);
        if (Fi.ok && Fo2.ok && Fi.R > 20 && Fo2.R > 20 && rhi2 - rli > 20)
        {
          double ki = signedK(Xi, Yi, Ri, Fi), ko = signedK(Xo, Yo, Ro, Fo2);
          double span = rhi2 - rli;
          double ds = (ki - ko) * span * span / 8 * 10;   // mm
          int sdT = (2 * nsd >= (int) kept.size()) ? 1 : 0;
          dsag[sdT].push_back(ds); dsagAll.push_back(ds);
        }
      }
    }
    for (int sd = 0; sd < 2; ++sd)
      for (int L = 7; L <= 54; ++L)
      {
        prof[sd][L] = cnt[sd][L] > 0 ? sum[sd][L] / cnt[sd][L] : 0;
        for (int sec = 0; sec < 12; ++sec)
          profSec[sd][sec][L] = cntS[sd][sec][L] > 0 ? sumS[sd][sec][L] / cntS[sd][sec][L] : 0;
      }
  };

  double profR[2][55] = {{0}}, profS[2][55] = {{0}}, profI[2][55] = {{0}};
  double secR[2][12][55] = {{{0}}}, secS[2][12][55] = {{{0}}}, secI[2][12][55] = {{{0}}};
  std::vector<double> dsR[2], dsS[2], dsI[2], dsRa, dsSa, dsIa;
  long ntR = 0, ntS = 0, ntI = 0;
  measure(0, nullptr, profR, secR, dsR, dsRa, ntR);
  measure(1, nullptr, profS, secS, dsS, dsSa, ntS);
  // closure: inject (real - sim) profile into sim
  double dprof[2][55] = {{0}};
  for (int sd = 0; sd < 2; ++sd) for (int L = 7; L <= 54; ++L) dprof[sd][L] = profR[sd][L] - profS[sd][L];
  measure(1, dprof, profI, secI, dsI, dsIa, ntI);
  // 3-level STEP models from the real profile's boundary jumps (rows 22->23, 38->39):
  double j1[2], j2[2];
  for (int sd = 0; sd < 2; ++sd) { j1[sd] = dprof[sd][23] - dprof[sd][22]; j2[sd] = dprof[sd][39] - dprof[sd][38]; }
  // B1: constant r-phi offset per region (levels 0 / j1 / j1+j2), minus the mean so the fit absorbs nothing extra
  double stepA[2][55] = {{0}}, stepB[2][55] = {{0}};
  for (int sd = 0; sd < 2; ++sd)
  {
    double lev[3] = {0, j1[sd], j1[sd] + j2[sd]};
    double mean = 0; for (int L = 7; L <= 54; ++L) mean += lev[L <= 22 ? 0 : L <= 38 ? 1 : 2]; mean /= 48;
    for (int L = 7; L <= 54; ++L) stepA[sd][L] = lev[L <= 22 ? 0 : L <= 38 ? 1 : 2] - mean;
    // B2: rigid ring rotation per region (constant dphi per region): dphi levels = jump / r_boundary
    double dph[3] = {0, j1[sd] / ((rowR[22] + rowR[23]) / 2), j1[sd] / ((rowR[22] + rowR[23]) / 2) + j2[sd] / ((rowR[38] + rowR[39]) / 2)};
    double mean2 = 0; for (int L = 7; L <= 54; ++L) mean2 += dph[L <= 22 ? 0 : L <= 38 ? 1 : 2] * rowR[L]; mean2 /= 48;
    for (int L = 7; L <= 54; ++L) stepB[sd][L] = dph[L <= 22 ? 0 : L <= 38 ? 1 : 2] * rowR[L] - mean2;
  }
  double profA[2][55] = {{0}}, profB[2][55] = {{0}}, secA[2][12][55] = {{{0}}}, secB[2][12][55] = {{{0}}};
  std::vector<double> dsA[2], dsB[2], dsAa, dsBa; long ntA = 0, ntB = 0;
  measure(1, stepA, profA, secA, dsA, dsAa, ntA);
  measure(1, stepB, profB, secB, dsB, dsBa, ntB);

  // ---- ledger + profile file ----------------------------------------------
  FILE *fo = fopen(Form("twist_probe_%s.txt", ver), "w");
  auto P = [&](const char *fmt, ...) {
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    va_start(ap, fmt); vfprintf(fo, fmt, ap); va_end(ap);
  };
  auto rsig = [&](std::vector<double> v) { double m = medv(v); for (double &q : v) q = std::fabs(q - m); return 1.4826 * medv(v); };
  P("[twist_probe %s] coherent r-phi twist profile: mean D(rphi) of kept pixels vs pad row, per side\n", ver);
  P("  tracks: real %ld | sim %ld | sim+injected %ld\n", ntR, ntS, ntI);
  P("  split-half Dsagitta med [mm] (robust sigma):\n");
  P("    real          all %+.2f (%.2f) | side0 %+.2f (n=%zu) | side1 %+.2f (n=%zu)\n", medv(dsRa), rsig(dsRa), medv(dsR[0]), dsR[0].size(), medv(dsR[1]), dsR[1].size());
  P("    sim           all %+.2f (%.2f) | side0 %+.2f (n=%zu) | side1 %+.2f (n=%zu)\n", medv(dsSa), rsig(dsSa), medv(dsS[0]), dsS[0].size(), medv(dsS[1]), dsS[1].size());
  P("    sim+injected  all %+.2f (%.2f) | side0 %+.2f (n=%zu) | side1 %+.2f (n=%zu)   <- CLOSURE vs real\n", medv(dsIa), rsig(dsIa), medv(dsI[0]), dsI[0].size(), medv(dsI[1]), dsI[1].size());
  P("  STEP MODELS (3 levels per side from the boundary jumps; j1 = R1->R2, j2 = R2->R3):\n");
  P("    jumps [um]: side0 j1 %+.0f j2 %+.0f | side1 j1 %+.0f j2 %+.0f\n", j1[0], j2[0], j1[1], j2[1]);
  auto rmsdiff = [&](double a[2][55]) { double s2 = 0; int n = 0; for (int sd = 0; sd < 2; ++sd) for (int L = 7; L <= 54; ++L) { double d = a[sd][L] - profR[sd][L]; s2 += d * d; n++; } return std::sqrt(s2 / n); };
  P("    B1 constant r-phi offset per region: sim+B1 Dsagitta med %+.2f mm | profile RMS diff to real %.0f um\n", medv(dsAa), rmsdiff(profA));
  P("    B2 rigid rotation per region (const dphi): sim+B2 Dsagitta med %+.2f mm | profile RMS diff to real %.0f um\n", medv(dsBa), rmsdiff(profB));
  P("    (full injected profile: Dsagitta %+.2f mm | profile RMS diff %.0f um ; sim alone: RMS diff %.0f um)\n", medv(dsIa), rmsdiff(profI), rmsdiff(profS));
  P("  profile [um]  (real | sim | real-sim = injected | sim+injected reproduces)\n");
  P("   row    side0: real   sim  delta  inj |  side1: real   sim  delta  inj\n");
  FILE *fp = fopen(Form("twist_profile_%s.txt", ver), "w");
  fprintf(fp, "# twist_profile_%s.txt — mean fit-orthogonal r-phi displacement [um] of kept pixels vs pad row,\n", ver);
  fprintf(fp, "# real (tracker-road groups, trimmed fits) minus sim digi v6 (truth groups). INJECT column 'delta'\n");
  fprintf(fp, "# as an azimuthal displacement D(rphi)=delta [um] at radius r: dphi = delta*1e-4/r, per side.\n");
  fprintf(fp, "# sign: positive = pixel at LARGER azimuth than the fitted circle (counter-clockwise in x-y).\n");
  fprintf(fp, "# layer side real_um sim_um delta_um\n");
  for (int L = 7; L <= 54; ++L)
  {
    P("   %2d   %7.0f %6.0f %6.0f %5.0f  |  %7.0f %6.0f %6.0f %5.0f\n", L,
      profR[0][L], profS[0][L], dprof[0][L], profI[0][L], profR[1][L], profS[1][L], dprof[1][L], profI[1][L]);
    for (int sd = 0; sd < 2; ++sd) fprintf(fp, "%d %d %.1f %.1f %.1f\n", L, sd, profR[sd][L], profS[sd][L], dprof[sd][L]);
  }
  fclose(fp);
  // sector coherence of the real profile: RMS across sectors of the R1 mean (rows 7-14) per side
  for (int sd = 0; sd < 2; ++sd)
  {
    std::vector<double> a;
    for (int sec = 0; sec < 12; ++sec)
    { double m = 0; int n = 0; for (int L = 7; L <= 14; ++L) if (secR[sd][sec][L] != 0) { m += secR[sd][sec][L]; n++; } if (n) a.push_back(m / n); }
    double mu = 0, s2 = 0; for (double q : a) mu += q; if (!a.empty()) mu /= a.size();
    for (double q : a) s2 += (q - mu) * (q - mu);
    P("  real side %d: R1 (rows 7-14) mean displacement per sector: mean %+.0f um, RMS across sectors %.0f um (values:", sd, mu, a.size() > 1 ? std::sqrt(s2 / (a.size() - 1)) : 0);
    for (double q : a) P(" %+.0f", q);
    P(")\n");
  }
  fclose(fo);
  printf("wrote twist_probe_%s.txt, twist_profile_%s.txt\n", ver, ver);

  // ---- figure ---------------------------------------------------------------
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.045);
  TCanvas *cv = new TCanvas("cvtw", "twist", 1700, 1250);
  cv->Divide(2, 2);
  TLatex tx; tx.SetNDC();
  auto profPanel = [&](int ipad, const char *title, double a0[2][55], double a1[2][55], const char *n0, const char *n1, double a2[2][55], const char *n2) {
    cv->cd(ipad);
    TH2D *fr = new TH2D(Form("twfr%d", ipad), Form("%s;pad row;#LT D(r#phi) #GT  [#mum]", title), 10, 6, 55, 10, -1600, 1600);
    fr->GetYaxis()->SetTitleOffset(1.25); fr->GetYaxis()->SetTitleSize(0.040);
    fr->Draw();
    TLine z; z.SetLineStyle(2); z.SetLineColor(kGray + 2); z.DrawLine(6, 0, 55, 0);
    TLine b1, b2; b1.SetLineStyle(3); b1.SetLineColor(kGray + 1); b1.DrawLine(22.5, -1600, 22.5, 1600); b1.DrawLine(38.5, -1600, 38.5, 1600);
    TLegend *L = new TLegend(0.40, 0.67, 0.89, 0.88); L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextSize(0.027);
    int col[2] = {kRed + 1, kBlue + 1};
    for (int sd = 0; sd < 2; ++sd)
    {
      TGraph *g0 = new TGraph, *g1 = new TGraph, *g2 = new TGraph;
      for (int Lr = 7; Lr <= 54; ++Lr) { g0->AddPoint(Lr, a0[sd][Lr]); g1->AddPoint(Lr, a1[sd][Lr]); if (a2) g2->AddPoint(Lr, a2[sd][Lr]); }
      g0->SetLineColor(col[sd]); g0->SetLineWidth(2); g0->SetMarkerStyle(20); g0->SetMarkerColor(col[sd]); g0->SetMarkerSize(0.7);
      g1->SetLineColor(col[sd]); g1->SetLineWidth(2); g1->SetLineStyle(2);
      g0->Draw("LP same"); g1->Draw("L same");
      L->AddEntry(g0, Form("%s, side %d", n0, sd), "lp");
      L->AddEntry(g1, Form("%s, side %d", n1, sd), "l");
      if (a2) { g2->SetLineColor(col[sd]); g2->SetLineWidth(3); g2->SetLineStyle(3); g2->Draw("L same"); L->AddEntry(g2, Form("%s, side %d", n2, sd), "l"); }
    }
    L->Draw();
  };
  profPanel(1, "twist profile: real vs sim (v6 digi)", profR, profS, "real", "sim", nullptr, "");
  profPanel(2, "closure: sim + injected (real#minussim) vs real", profR, profI, "real", "sim + injected", nullptr, "");
  // [3] Dsagitta distributions
  cv->cd(3);
  {
    gPad->SetLogy();
    TH1D *h0 = new TH1D("tw_ds0", "split-half #Delta sagitta: real / sim / sim+injected;#Delta sagitta [mm];tracks (unit area, log)", 100, -12.5, 12.5);
    TH1D *h1 = (TH1D *) h0->Clone("tw_ds1"); TH1D *h2 = (TH1D *) h0->Clone("tw_ds2");
    for (double q : dsRa) h0->Fill(std::max(-12.49, std::min(12.49, q)));
    for (double q : dsSa) h1->Fill(std::max(-12.49, std::min(12.49, q)));
    for (double q : dsIa) h2->Fill(std::max(-12.49, std::min(12.49, q)));
    for (auto h : {h0, h1, h2}) h->Scale(1. / h->Integral());
    h0->SetLineColor(kBlack); h0->SetLineWidth(2);
    h1->SetLineColor(kBlue + 1); h1->SetLineWidth(2); h1->SetLineStyle(2);
    h2->SetLineColor(kGreen + 2); h2->SetLineWidth(2);
    h0->SetMaximum(6 * std::max({h0->GetMaximum(), h1->GetMaximum(), h2->GetMaximum()})); h0->SetMinimum(1e-4);
    h0->Draw("hist"); h1->Draw("hist same"); h2->Draw("hist same");
    TLegend *L = new TLegend(0.14, 0.70, 0.60, 0.88); L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextSize(0.028);
    L->AddEntry(h0, Form("real, med %+.2f mm", medv(dsRa)), "l");
    L->AddEntry(h1, Form("sim, med %+.2f mm", medv(dsSa)), "l");
    L->AddEntry(h2, Form("sim + injected profile, med %+.2f mm", medv(dsIa)), "l");
    L->Draw();
  }
  // [4] sector coherence: R1 mean displacement per (side, sector)
  cv->cd(4);
  {
    gPad->SetRightMargin(0.14);
    TH2D *hm = new TH2D("tw_sec", "real R1 (rows 7-14) mean D(r#phi) by side #times sector [#mum];sector;side", 12, -0.5, 11.5, 2, -0.5, 1.5);
    for (int sd = 0; sd < 2; ++sd)
      for (int sec = 0; sec < 12; ++sec)
      { double m = 0; int n = 0; for (int L = 7; L <= 14; ++L) if (secR[sd][sec][L] != 0) { m += secR[sd][sec][L]; n++; } hm->SetBinContent(sec + 1, sd + 1, n ? m / n : 0); }
    gStyle->SetPaintTextFormat(".0f");
    hm->SetMarkerSize(1.4);
    hm->Draw("colz text");
  }
  cv->SaveAs(Form("../sim_validation_plots/twist_probe_%s.png", ver));
  printf("wrote ../sim_validation_plots/twist_probe_%s.png\n", ver);
}
