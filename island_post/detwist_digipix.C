// detwist_digipix.C — DE-TWIST the real pixels (subtract real's own measured
// per-(side,row) twist profile from the pixel azimuths), rerun the EXACT
// nf_digipix analysis, and plot in the nf_digipix STYLE (user, 2026-08-24):
// two panels, GLOBAL whole-track fit and LOCAL 4-row short-sagitta fit,
// per-fit circle RMS distributions (unit area), medians in the legend —
// now with three curves: real, real DE-TWISTED, sim digi (v6.1).
// Question answered in the ledger: does de-twisting reduce the real med RMS
// by more than half?
// Method: pass A measures the twist profile exactly as twist_probe.C (trimmed
// global fits, mean D(rphi) of kept pixels per side and pad row); pass B
// rotates every real pixel by -profile(side,row)/r and reruns the fits.
// Groups, gates, fitter = verbatim nf_digipix (road dxy<1.2 cm |dtbin|<=6,
// ev44 veto, bar >=12 px span>=15 45<=R<2e4; windows 4 rows, >=3 rows, >=5 px).
// usage: root -l -b -q 'detwist_digipix.C+()'
// out: detwist_digipix_<ver>.txt + ../sim_validation_plots/detwist_digipix_<ver>.png
#include "../sim_validation_plots/src/ms_nofinder.C"
#include <TROOT.h>
#include <cstdint>

namespace DTW
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
}  // namespace DTW

void detwist_digipix(const char *digif = "digi_frames_production_v61.root", int nsim = 60,
                     const char *realf = "../clusters_seeds_island_79507-0.root_ntuplizer.root",
                     const char *ver = "v61")
{
  using namespace DTW;
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

  // ---- groups ---------------------------------------------------------------
  std::vector<GrpX> gr[2];   // [0]=real [1]=sim
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
    float ev, lay, phi, adc, tid;
    t->SetBranchStatus("*", 0);
    for (auto b : {"event", "layer", "phi", "adc", "gtrackID"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("event", &ev); t->SetBranchAddress("layer", &lay);
    t->SetBranchAddress("phi", &phi); t->SetBranchAddress("adc", &adc);
    t->SetBranchAddress("gtrackID", &tid);
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
      G.side.push_back(0); G.lay.push_back((uint8_t) lay);
    }
    f->Close();
    printf("sim: %zu truth groups\n", gr[1].size());
  }

  // ---- pass A: measure the real twist profile (as twist_probe.C) ------------
  double prof[2][55] = {{0}};
  {
    double sum[2][55] = {{0}}, cnt[2][55] = {{0}};
    for (auto &G : gr[0])
    {
      MNF::Grp B; B.x = G.x; B.y = G.y; B.r = G.r;
      MNF::Fit F0;
      if (!MNF::fitBar(B, F0)) continue;
      MNF::Fit FT; std::vector<int> kept;
      if (!trimFit(G.x, G.y, FT, kept)) continue;
      for (int i : kept)
      {
        int L = G.lay[i]; if (L < 7 || L > 54) continue;
        double phf;
        if (!phiAtR(FT, G.r[i], G.x[i], G.y[i], phf)) continue;
        double d = MNF::wrapphi(std::atan2(G.y[i], G.x[i]) - phf) * G.r[i] * 1e4;
        if (std::fabs(d) > 8000) continue;
        sum[G.side[i]][L] += d; cnt[G.side[i]][L]++;
      }
    }
    for (int sd = 0; sd < 2; ++sd)
      for (int L = 7; L <= 54; ++L)
        prof[sd][L] = cnt[sd][L] > 0 ? sum[sd][L] / cnt[sd][L] : 0;
    printf("profile measured (e.g. side0 row7 %+.0f um, row22 %+.0f um)\n", prof[0][7], prof[0][22]);
  }
  // ---- de-twisted real groups ----------------------------------------------
  std::vector<GrpX> grD = gr[0];
  for (auto &G : grD)
    for (size_t i = 0; i < G.x.size(); ++i)
    {
      int L = G.lay[i]; if (L < 7 || L > 54) continue;
      double r = G.r[i];
      double ph = std::atan2(G.y[i], G.x[i]) - prof[G.side[i]][L] * 1e-4 / r;
      G.x[i] = r * std::cos(ph); G.y[i] = r * std::sin(ph);
    }

  // ---- nf_digipix analysis (verbatim gates) on the three samples ------------
  auto analyze = [&](std::vector<GrpX> &gv, std::vector<double> &grms, std::vector<double> &wrms, long &npx) {
    long ntr = 0; npx = 0;
    for (auto &G : gv)
    {
      MNF::Grp B; B.x = G.x; B.y = G.y; B.r = G.r;
      MNF::Fit F0;
      if (!MNF::fitBar(B, F0)) continue;
      ntr++; npx += (long) G.x.size();
      grms.push_back(F0.rms * 1e4);
      std::vector<std::vector<double>> wx(45), wy(45);
      std::vector<std::set<int>> wrow(45);
      for (size_t i = 0; i < G.x.size(); ++i)
      {
        int row = nearRow(G.r[i]);
        if (row < 0) continue;
        for (int w = std::max(7, row - 3); w <= std::min(51, row); ++w)
        { wx[w - 7].push_back(G.x[i]); wy[w - 7].push_back(G.y[i]); wrow[w - 7].insert(row); }
      }
      for (int w = 0; w < 45; ++w)
      {
        if ((int) wrow[w].size() < 3 || (int) wx[w].size() < 5) continue;
        MNF::Fit L = MNF::fitCircle(wx[w], wy[w]);
        if (L.ok) wrms.push_back(L.rms * 1e4);
      }
    }
    return ntr;
  };
  std::vector<double> gR, wR, gD, wD, gS, wS;
  long pxR, pxD, pxS;
  long nR = analyze(gr[0], gR, wR, pxR);
  long nD = analyze(grD, gD, wD, pxD);
  long nS = analyze(gr[1], gS, wS, pxS);

  // ---- ledger ---------------------------------------------------------------
  FILE *fo = fopen(Form("detwist_digipix_%s.txt", ver), "w");
  auto P = [&](const char *fmt, ...) {
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    va_start(ap, fmt); vfprintf(fo, fmt, ap); va_end(ap);
  };
  P("[detwist_digipix %s] nf_digipix meters on DE-TWISTED real (own measured profile subtracted)\n", ver);
  P("  real          : %ld tracks | GLOBAL med %.0f | LOCAL med %.0f um\n", nR, medv(gR), medv(wR));
  P("  real DE-TWIST : %ld tracks | GLOBAL med %.0f | LOCAL med %.0f um\n", nD, medv(gD), medv(wD));
  P("  sim digi %s   : %ld tracks | GLOBAL med %.0f | LOCAL med %.0f um\n", ver, nS, medv(gS), medv(wS));
  double redG = 100 * (1 - medv(gD) / medv(gR)), redL = 100 * (1 - medv(wD) / medv(wR));
  P("  REDUCTION: GLOBAL %.1f%%, LOCAL %.1f%%  ->  more than half? %s (would need >50%%)\n",
    redG, redL, (redG > 50 || redL > 50) ? "YES" : "NO");
  P("  data/MC: GLOBAL %.2f -> %.2f | LOCAL %.2f -> %.2f\n",
    medv(gR) / medv(gS), medv(gD) / medv(gS), medv(wR) / medv(wS), medv(wD) / medv(wS));
  fclose(fo);
  printf("wrote detwist_digipix_%s.txt\n", ver);

  // ---- figure: nf_digipix style, three curves per panel ---------------------
  gStyle->SetOptStat(0);
  TCanvas *cv = new TCanvas("cvdt", "detwist digipix", 1500, 620);
  cv->Divide(2, 1);
  const char *pt[2] = {"GLOBAL whole-track fit", "LOCAL 4-row short-sagitta fit"};
  for (int p = 0; p < 2; ++p)
  {
    cv->cd(p + 1);
    std::vector<double> *v[3] = {p == 0 ? &gR : &wR, p == 0 ? &gD : &wD, p == 0 ? &gS : &wS};
    double xhi = p == 0 ? 6000 : 4000;
    TH1D *h[3];
    for (int s = 0; s < 3; ++s)
    {
      h[s] = new TH1D(Form("dtw_%d_%d", p, s), ";per-fit circle RMS [#mum];fits (unit area)", 60, 0, xhi);
      for (double q : *v[s]) h[s]->Fill(std::min(q, xhi - 1));
      if (h[s]->Integral() > 0) h[s]->Scale(1. / h[s]->Integral());
    }
    h[0]->SetLineColor(kBlack); h[0]->SetLineWidth(2);
    h[1]->SetLineColor(kRed + 1); h[1]->SetLineWidth(2);
    h[2]->SetLineColor(kBlue + 1); h[2]->SetLineWidth(2); h[2]->SetLineStyle(2);
    h[0]->SetTitle(Form("matched pixel level: %s", pt[p]));
    h[0]->SetMaximum(1.4 * std::max({h[0]->GetMaximum(), h[1]->GetMaximum(), h[2]->GetMaximum()}));
    h[0]->Draw("hist");
    h[1]->Draw("hist same");
    h[2]->Draw("hist same");
    TLegend *lg = new TLegend(0.38, 0.66, 0.89, 0.88);
    lg->SetBorderSize(0);
    lg->AddEntry(h[0], Form("real pixels, med %.0f #mum", medv(*v[0])), "l");
    lg->AddEntry(h[1], Form("real DE-TWISTED, med %.0f #mum", medv(*v[1])), "l");
    lg->AddEntry(h[2], Form("sim digi pixels (%s), med %.0f #mum", ver, medv(*v[2])), "l");
    lg->Draw();
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.030);
    tx.DrawLatex(0.38, 0.60, Form("data/MC %.2f #rightarrow %.2f (de-twisted)",
                                  medv(*v[0]) / medv(*v[2]), medv(*v[1]) / medv(*v[2])));
    tx.SetTextSize(0.026);
    if (p == 0)
      tx.DrawLatex(0.38, 0.55, Form("de-twist reduction %.1f%% (half = 50%%)", 100 * (1 - medv(*v[1]) / medv(*v[0]))));
    else
      tx.DrawLatex(0.38, 0.55, Form("de-twist reduction %.1f%%; floor = cloud width", 100 * (1 - medv(*v[1]) / medv(*v[0]))));
  }
  cv->SaveAs(Form("../sim_validation_plots/detwist_digipix_%s.png", ver));
  printf("wrote ../sim_validation_plots/detwist_digipix_%s.png\n", ver);
}
