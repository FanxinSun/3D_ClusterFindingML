// cfloor_probe.C — FLOOR-SUBTRACTED / DECOMPOSED composition meters for the
// response campaign (clause-6 resolution 2026-08-23: no bands on raw C(1);
// the raw real-vs-sim C(1) contrast sits in the response floor, so the
// composition instruments must remove/decompose the floor first).
// Definitions (per side pooled over tracks; nf_digipix groups, trimmed fits):
//   m_L      = row-mean residual of kept pixels in pad row L (to the trimmed fit)
//   s2_L     = HF_track^2 / n_L = STATISTICAL floor of m_L (within-row scatter
//              propagated to the row mean; HF_track = within-row RMS pooled
//              over the track's >=2-px rows)
//   C(d)     = < m_i * m_{i+d} >  (raw autocovariance vs row lag, as nonrms)
//   Csig(0)  = C(0) - <s2>        (statistical floor removed; stat noise is
//              row-white so it lives ONLY in C(0))
//   C_sub(d) = C(d) / Csig(0), d>=1   ["C_floor(d)" of the resolution]
//   shares of the non-statistical signal variance:
//     SMOOTH = P = mean C_sub(d=10..14)   (long-range plateau incl. twist)
//     SHORT  = C_sub(1) - P               (granular, correlated ~1-3 rows)
//     WHITE  = 1 - C_sub(1)               (row-white NON-statistical content =
//                                          the response-floor systematic)
//   cellcoh_sub^2 = cellRMS^2 - <stat var of the cell mean>  (per side,sector,layer)
// CONTROL: "sim thinned" = sim pixels Bernoulli-subsampled (deterministic LCG)
// to real's px/row (keep p = 3.25/6.67) -> how much of the raw contrast my
// statistic attributes to sampling alone, on THIS harness.
// usage: root -l -b -q 'cfloor_probe.C+()'
// out: cfloor_probe_<ver>.txt + ../sim_validation_plots/cfloor_probe_<ver>.png
#include "../sim_validation_plots/src/ms_nofinder.C"
#include <TROOT.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TLine.h>
#include <unordered_map>
#include <cstdint>

namespace CFP
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
struct Acc
{
  double C[31] = {0}, Cn[31] = {0}, s2sum = 0, s2n = 0;
  std::map<int, std::pair<double, double>> cell;   // (side*12+sec)*100+L -> (sum m, sum s2), count in cn
  std::map<int, long> celln;
  long ntrk = 0;
};
}  // namespace CFP

void cfloor_probe(const char *digif = "digi_frames_production_v6.root", int nsim = 60,
                  const char *realf = "../clusters_seeds_island_79507-0.root_ntuplizer.root",
                  const char *ver = "v6t")
{
  using namespace CFP;
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

  // ---- groups: real (road) + sim (truth) ------------------------------------
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

  // ---- engine: accumulate C(d), stat floor, cells; optional thinning --------
  const double KEEPP = 3.25 / 6.67;   // thinned-sim keep probability (real/sim px per row)
  auto process = [&](int s, bool thin, Acc &A) {
    for (size_t gi = 0; gi < gr[s].size(); ++gi)
    {
      const GrpX &G0 = gr[s][gi];
      std::vector<double> X, Y, R2; std::vector<uint8_t> SD, LY;
      if (!thin) { X = G0.x; Y = G0.y; R2 = G0.r; SD = G0.side; LY = G0.lay; }
      else
      {
        unsigned lcg = 77777u + (unsigned) gi * 2654435761u;
        auto rnd = [&]() { lcg = lcg * 1664525u + 1013904223u; return (lcg >> 8) / 16777216.; };
        for (size_t i = 0; i < G0.x.size(); ++i)
          if (rnd() < KEEPP)
          { X.push_back(G0.x[i]); Y.push_back(G0.y[i]); R2.push_back(G0.r[i]); SD.push_back(G0.side[i]); LY.push_back(G0.lay[i]); }
      }
      if ((int) X.size() < 12) continue;
      double rlo = 1e9, rhi = 0;
      for (double q : R2) { rlo = std::min(rlo, q); rhi = std::max(rhi, q); }
      if (rhi - rlo < 15) continue;
      MNF::Fit F0 = MNF::fitCircle(X, Y);
      if (!F0.ok || F0.R < 45 || F0.R >= 2e4) continue;
      MNF::Fit FT; std::vector<int> kept;
      if (!trimFit(X, Y, FT, kept)) continue;
      A.ntrk++;
      // rows
      std::map<int, std::vector<int>> rows;
      for (int i : kept) { int L = LY[i]; if (L >= 7 && L <= 54) rows[L].push_back(i); }
      // track HF (within-row RMS pooled, rows >=2 px)
      double hf2 = 0; long nhf = 0;
      std::map<int, double> m;
      for (auto &kv : rows)
      {
        double mm = 0;
        for (int i : kv.second) mm += std::hypot(X[i] - FT.a, Y[i] - FT.b) - FT.R;
        mm /= kv.second.size();
        m[kv.first] = mm;
        if (kv.second.size() >= 2)
          for (int i : kv.second)
          { double d = (std::hypot(X[i] - FT.a, Y[i] - FT.b) - FT.R) - mm; hf2 += d * d; nhf++; }
      }
      double HF2 = nhf > 3 ? hf2 / nhf : 0;
      // C(d) + stat floor + cells
      for (auto &a : m)
        for (auto &b : m)
        {
          int d = b.first - a.first;
          if (d < 0 || d > 30) continue;
          A.C[d] += a.second * b.second * 1e8;   // um^2
          A.Cn[d]++;
        }
      for (auto &kv : rows)
      {
        double s2 = HF2 / kv.second.size() * 1e8;
        A.s2sum += s2; A.s2n++;
        // cell accumulation
        int i0 = kv.second[0];
        double mphi = std::atan2(Y[i0], X[i0]);
        int sec = ((int) std::floor((mphi + 2 * M_PI) / (M_PI / 6))) % 12;
        int sd = SD[i0];
        int key = (sd * 12 + sec) * 100 + kv.first;
        A.cell[key].first += m[kv.first] * 1e4;   // um
        A.cell[key].second += s2;
        A.celln[key]++;
      }
    }
  };
  Acc AR, AS, ATh;
  process(0, false, AR);
  process(1, false, AS);
  process(1, true, ATh);
  printf("tracks: real %ld sim %ld sim-thinned %ld\n", AR.ntrk, AS.ntrk, ATh.ntrk);

  // ---- derive metrics --------------------------------------------------------
  auto metrics = [&](Acc &A, double out_Csub[31], double &c0, double &s2m, double &W, double &SH, double &P,
                     double &cellraw, double &cellsub) {
    c0 = A.C[0] / A.Cn[0];
    s2m = A.s2sum / A.s2n;
    double csig0 = c0 - s2m;
    for (int d = 0; d <= 30; ++d) out_Csub[d] = (A.Cn[d] > 0 && csig0 > 0) ? (A.C[d] / A.Cn[d]) / csig0 : 0;
    P = 0; int np = 0;
    for (int d = 10; d <= 14; ++d) if (A.Cn[d] > 0) { P += out_Csub[d]; np++; }
    P = np ? P / np : 0;
    W = 1 - out_Csub[1];
    SH = out_Csub[1] - P;
    // cells (>=8 entries)
    double s2c = 0, sc2 = 0; long nc = 0; double statc = 0;
    for (auto &kv : A.cell)
    {
      long k = A.celln[kv.first];
      if (k < 8) continue;
      double mean = kv.second.first / k;
      s2c += mean; sc2 += mean * mean; nc++;
      statc += kv.second.second / (k * k);   // um^2 stat var of the cell mean
    }
    if (nc > 1)
    {
      double mu = s2c / nc;
      cellraw = std::sqrt(std::max(0., sc2 / nc - mu * mu));
      cellsub = std::sqrt(std::max(0., sc2 / nc - mu * mu - statc / nc));
    }
    else cellraw = cellsub = 0;
  };
  double CsR[31], CsS[31], CsT[31], c0[3], s2[3], W[3], SH[3], P[3], cr[3], cs[3];
  metrics(AR, CsR, c0[0], s2[0], W[0], SH[0], P[0], cr[0], cs[0]);
  metrics(AS, CsS, c0[1], s2[1], W[1], SH[1], P[1], cr[1], cs[1]);
  metrics(ATh, CsT, c0[2], s2[2], W[2], SH[2], P[2], cr[2], cs[2]);

  // ---- ledger ----------------------------------------------------------------
  FILE *fo = fopen(Form("cfloor_probe_%s.txt", ver), "w");
  auto Pr = [&](const char *fmt, ...) {
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    va_start(ap, fmt); vfprintf(fo, fmt, ap); va_end(ap);
  };
  Pr("[cfloor_probe %s] floor-subtracted composition meters (response-campaign instruments)\n", ver);
  Pr("  definitions in the header; thinned = sim subsampled to real px/row (keep p=%.2f)\n", KEEPP);
  const char *nm[3] = {"real       ", "sim        ", "sim thinned"};
  for (int k = 0; k < 3; ++k)
  {
    double *Cs = k == 0 ? CsR : k == 1 ? CsS : CsT;
    Pr("  %s: C(0) %.0f um^2 | stat floor %.0f um^2 (%.0f%%) | Csig(0) %.0f um^2 (= %0.f um rms)\n",
       nm[k], c0[k], s2[k], 100 * s2[k] / c0[k], c0[k] - s2[k], std::sqrt(std::max(0., c0[k] - s2[k])));
    Pr("      raw C(1)/C(0) %.2f -> FLOOR-SUB C_sub(1) %.2f | C_sub(2) %.2f | C_sub(4) %.2f | C_sub(8) %.2f | plateau(10-14) %.2f\n",
       (AR.C[1] / AR.Cn[1]) / c0[0] * 0 + (k == 0 ? (AR.C[1] / AR.Cn[1]) / c0[0] : k == 1 ? (AS.C[1] / AS.Cn[1]) / c0[1] : (ATh.C[1] / ATh.Cn[1]) / c0[2]),
       Cs[1], Cs[2], Cs[4], Cs[8], P[k]);
    Pr("      SHARES of non-stat signal: WHITE %.2f | SHORT(1-3 rows) %.2f | SMOOTH %.2f\n", W[k], SH[k], P[k]);
    Pr("      cell coherence: raw %.0f um -> floor-subtracted %.0f um\n", cr[k], cs[k]);
  }
  Pr("  reading guide: WHITE = response-floor systematic (row-white, non-statistical);\n");
  Pr("  SHORT = granular field/alignment; SMOOTH = long-range (SPHI+twist).\n");
  fclose(fo);
  printf("wrote cfloor_probe_%s.txt\n", ver);

  // ---- figure ----------------------------------------------------------------
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.05);
  TCanvas *cv = new TCanvas("cvcf", "cfloor", 1700, 640);
  cv->Divide(2, 1);
  cv->cd(1);
  {
    TH2D *fr = new TH2D("cf1", "floor-subtracted autocorrelation C_{sub}(d);row lag d;C_{sub}(d)", 10, 0, 26, 10, -0.45, 1.3);
    fr->GetYaxis()->SetTitleOffset(1.1);
    fr->Draw();
    TLine z; z.SetLineStyle(2); z.SetLineColor(kGray + 2); z.DrawLine(0, 0, 26, 0);
    TGraph *g0 = new TGraph, *g1 = new TGraph, *g2 = new TGraph;
    for (int d = 1; d <= 25; ++d)
    { g0->AddPoint(d, CsR[d]); g1->AddPoint(d, CsS[d]); g2->AddPoint(d, CsT[d]); }
    g0->SetLineColor(kBlack); g0->SetLineWidth(2); g0->SetMarkerStyle(20);
    g1->SetLineColor(kBlue + 1); g1->SetLineWidth(2); g1->SetLineStyle(2); g1->SetMarkerStyle(24); g1->SetMarkerColor(kBlue + 1);
    g2->SetLineColor(kOrange + 7); g2->SetLineWidth(2); g2->SetLineStyle(3); g2->SetMarkerStyle(26); g2->SetMarkerColor(kOrange + 7);
    g0->Draw("LP same"); g1->Draw("LP same"); g2->Draw("LP same");
    TLegend *L = new TLegend(0.44, 0.68, 0.89, 0.88); L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextSize(0.033);
    L->AddEntry(g0, Form("real, C_{sub}(1) = %.2f", CsR[1]), "lp");
    L->AddEntry(g1, Form("sim v6t, C_{sub}(1) = %.2f", CsS[1]), "lp");
    L->AddEntry(g2, Form("sim thinned to real px/row, %.2f", CsT[1]), "lp");
    L->Draw();
  }
  cv->cd(2);
  {
    TH1D *hb[3];
    const char *bn[3] = {"real", "sim v6t", "sim thinned"};
    int cc[3] = {kGray + 3, kBlue - 7, kOrange - 3};
    TH2D *fr = new TH2D("cf2", "non-statistical signal shares;;share of C_{sig}(0)", 3, 0, 3, 10, 0, 1.50);
    for (int k = 0; k < 3; ++k) fr->GetXaxis()->SetBinLabel(k + 1, bn[k]);
    fr->GetXaxis()->SetLabelSize(0.055);
    fr->Draw();
    // stacked bars: WHITE bottom, SHORT mid, SMOOTH top
    for (int k = 0; k < 3; ++k)
    {
      double x0 = k + 0.25, x1 = k + 0.75;
      double y0 = 0, vals[3] = {W[k], SH[k], P[k]};
      int vc[3] = {kAzure - 9, kOrange - 3, kRed - 7};
      for (int j = 0; j < 3; ++j)
      {
        TBox *b = new TBox(x0, y0, x1, y0 + std::max(0., vals[j]));
        b->SetFillColor(vc[j]); b->Draw();
        y0 += std::max(0., vals[j]);
      }
    }
    TLegend *L = new TLegend(0.14, 0.72, 0.60, 0.88); L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextSize(0.033);
    TBox *b1 = new TBox(); b1->SetFillColor(kAzure - 9);
    TBox *b2 = new TBox(); b2->SetFillColor(kOrange - 3);
    TBox *b3 = new TBox(); b3->SetFillColor(kRed - 7);
    L->AddEntry(b1, "WHITE (response-floor systematic)", "f");
    L->AddEntry(b2, "SHORT (granular, 1-3 rows)", "f");
    L->AddEntry(b3, "SMOOTH (plateau: SPHI + twist)", "f");
    L->Draw();
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.032);
    tx.DrawLatex(0.14, 0.60, Form("cell coherence raw #rightarrow sub [#mum]:"));
    tx.DrawLatex(0.14, 0.54, Form("real %.0f #rightarrow %.0f | sim %.0f #rightarrow %.0f | thin %.0f #rightarrow %.0f",
                                  cr[0], cs[0], cr[1], cs[1], cr[2], cs[2]));
    tx.SetTextSize(0.029); tx.SetTextColor(kRed + 1);
    tx.DrawLatex(0.14, 0.46, Form("caveat: real's plateau is NEGATIVE (%.2f, oscillatory", P[0]));
    tx.DrawLatex(0.14, 0.41, "twist) #Rightarrow 3-band shares overflow for real;");
    tx.DrawLatex(0.14, 0.36, "judge real by the C_{sub}(d) curve and M1-M4");
  }
  cv->SaveAs(Form("../sim_validation_plots/cfloor_probe_%s.png", ver));
  printf("wrote ../sim_validation_plots/cfloor_probe_%s.png\n", ver);
}
