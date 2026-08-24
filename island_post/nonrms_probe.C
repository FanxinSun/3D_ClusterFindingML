// nonrms_probe.C — metrics BEYOND fit RMS for the real-vs-sim pixel meter
// (user, 2026-08-20: "really no other metrics?"). Same groups as nf_digipix
// (real = tracker road at nominal 1.2 cm / +-6 tb, sim = truth-grouped digi),
// same fitter, global 3xMAD trim; then four non-RMS axes:
//   (1) POOLED signed residual distribution (every pixel, residual to the
//       trimmed fit): core width, tail fractions, skew — a shape, not a
//       per-fit scalar; association pickup appears as mm-scale shoulders.
//   (2) PULLS: residual / that track's own robust sigma (1.4826xMAD of kept
//       residuals) — dimensionless; core should be ~unit width both sides,
//       tails = association, independent of any scale difference.
//   (3) ROW-LAG AUTOCORRELATION C(d) of row-mean residuals along the track:
//       white noise -> ~0; coherent field/alignment -> positive, long-range;
//       the full frequency picture (LF/HF was its 2-bin version).
//   (4) FIT-PARAMETER physics: |d0| distribution; split-half curvature
//       consistency Delta-sagitta = (k_in - k_out) * span^2 / 8 (signed k via
//       endpoint cross product) — parameter level, blind to the cloud floor.
// usage: root -l -b -q 'nonrms_probe.C+()'
// out: nonrms_probe_<ver>.txt + ../sim_validation_plots/nonrms_probe_<ver>.png
#include "../sim_validation_plots/src/ms_nofinder.C"
#include <TROOT.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TLine.h>
#include <unordered_map>
#include <cstdint>

namespace NRP
{
struct GrpX { std::vector<double> x, y, r; };
double medv(std::vector<double> v)
{
  if (v.empty()) return 0;
  std::sort(v.begin(), v.end());
  return v[v.size() / 2];
}
bool trimFit(const GrpX &G, MNF::Fit &F, std::vector<int> &kept)
{
  kept.resize(G.x.size());
  for (size_t k = 0; k < kept.size(); ++k) kept[k] = (int) k;
  for (int it = 0; it < 3; ++it)
  {
    std::vector<double> X, Y;
    for (int i : kept) { X.push_back(G.x[i]); Y.push_back(G.y[i]); }
    F = MNF::fitCircle(X, Y);
    if (!F.ok) return false;
    std::vector<double> res;
    for (int i : kept) res.push_back(std::hypot(G.x[i] - F.a, G.y[i] - F.b) - F.R);
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
  std::vector<double> X, Y;
  for (int i : kept) { X.push_back(G.x[i]); Y.push_back(G.y[i]); }
  F = MNF::fitCircle(X, Y);
  return F.ok;
}
// signed curvature: sign from the cross product of (mid-first)x(last-first),
// points ordered by radius; k = sign / R
double signedK(const std::vector<double> &X, const std::vector<double> &Y,
               const std::vector<double> &Rr, const MNF::Fit &F)
{
  size_t i0 = 0, i1 = 0, i2 = 0;
  double rlo = 1e9, rhi = -1;
  for (size_t i = 0; i < Rr.size(); ++i)
  {
    if (Rr[i] < rlo) { rlo = Rr[i]; i0 = i; }
    if (Rr[i] > rhi) { rhi = Rr[i]; i2 = i; }
  }
  double rmid = 0.5 * (rlo + rhi), bd = 1e9;
  for (size_t i = 0; i < Rr.size(); ++i)
    if (std::fabs(Rr[i] - rmid) < bd) { bd = std::fabs(Rr[i] - rmid); i1 = i; }
  double cx = (X[i1] - X[i0]) * (Y[i2] - Y[i0]) - (Y[i1] - Y[i0]) * (X[i2] - X[i0]);
  return (cx >= 0 ? 1. : -1.) / F.R;
}
}  // namespace NRP

void nonrms_probe(const char *digif = "digi_frames_production_v6.root", int nsim = 60,
                  const char *realf = "../clusters_seeds_island_79507-0.root_ntuplizer.root",
                  const char *ver = "v6")
{
  using namespace NRP;
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

  // ---- groups (nominal road real / truth sim; verbatim as gtail_probe) -----
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
    float adc, hitID;
    t->SetBranchStatus("*", 0);
    for (auto b : {"event", "layer", "x", "y", "tbin", "adc", "hitID"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("event", &ev); t->SetBranchAddress("layer", &lay);
    t->SetBranchAddress("x", &x); t->SetBranchAddress("y", &y);
    t->SetBranchAddress("tbin", &tb); t->SetBranchAddress("adc", &adc);
    t->SetBranchAddress("hitID", &hitID);
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
    }
    f->Close();
    printf("sim: %zu truth groups\n", gr[1].size());
  }

  // ---- battery ---------------------------------------------------------------
  TH1D *hres[2], *hpull[2], *hd0[2], *hds[2];
  for (int s = 0; s < 2; ++s)
  {
    hres[s] = new TH1D(Form("nr_res%d", s), "pooled signed residual to the trimmed fit;residual [mm];pixels (unit area, log)", 120, -15, 15);
    hpull[s] = new TH1D(Form("nr_pull%d", s), "pooled pulls: residual / track robust #sigma;pull;pixels (unit area, log)", 100, -12.5, 12.5);
    hd0[s] = new TH1D(Form("nr_d0%d", s), "impact parameter of the trimmed fit;|d0| [cm];tracks (unit area)", 60, 0, 12);
    hds[s] = new TH1D(Form("nr_ds%d", s), "split-half curvature consistency;#Delta sagitta = (k_{in}#minus k_{out}) span^{2}/8 [mm];tracks (unit area)", 100, -12.5, 12.5);
  }
  double AC[2][31] = {{0}}, ACn[2][31] = {{0}};
  std::vector<double> tailr[2], skew[2], pull3[2], d0v[2], dsv[2], Rv[2];
  std::vector<double> dsq[2][2];   // [side][charge sign 0=neg 1=pos] split-half Dsagitta
  long npix[2] = {0, 0}, ntrk2[2] = {0, 0};
  for (int s = 0; s < 2; ++s)
    for (auto &G : gr[s])
    {
      MNF::Grp B; B.x = G.x; B.y = G.y; B.r = G.r;
      MNF::Fit F0;
      if (!MNF::fitBar(B, F0)) continue;
      MNF::Fit FT; std::vector<int> kept;
      if (!trimFit(G, FT, kept)) continue;
      ntrk2[s]++;
      Rv[s].push_back(FT.R);
      // (1)+(2) pooled residuals and pulls (ALL pixels, residual to trimmed fit)
      std::vector<double> resk;
      for (int i : kept) resk.push_back(std::hypot(G.x[i] - FT.a, G.y[i] - FT.b) - FT.R);
      std::vector<double> tmp = resk;
      double mdk = medv(tmp);
      for (double &q : tmp) q = std::fabs(q - mdk);
      double srob = std::max(1.4826 * medv(tmp), 1e-3);
      double s3 = 0, s2 = 0, ntl = 0;
      for (size_t i = 0; i < G.x.size(); ++i)
      {
        double res = std::hypot(G.x[i] - FT.a, G.y[i] - FT.b) - FT.R;
        hres[s]->Fill(std::max(-14.99, std::min(14.99, res * 10)));
        hpull[s]->Fill(std::max(-12.49, std::min(12.49, res / srob)));
        if (std::fabs(res / srob) > 3) ntl++;
        npix[s]++;
      }
      pull3[s].push_back(ntl / G.x.size());
      for (double q : resk) { s2 += q * q; s3 += q * q * q; }
      if (s2 > 0) skew[s].push_back(s3 / resk.size() / std::pow(s2 / resk.size(), 1.5));
      double tl = 0;
      for (size_t i = 0; i < G.x.size(); ++i)
      { double res = std::hypot(G.x[i] - FT.a, G.y[i] - FT.b) - FT.R; if (std::fabs(res) > 0.3) tl++; }
      tailr[s].push_back(tl / G.x.size());
      // (3) row-lag autocorrelation of row-mean residuals (kept)
      std::map<int, std::pair<double, int>> rowm;
      for (int i : kept)
      {
        int L = nearRow(G.r[i]);
        if (L < 0) continue;
        double res = std::hypot(G.x[i] - FT.a, G.y[i] - FT.b) - FT.R;
        rowm[L].first += res; rowm[L].second++;
      }
      std::map<int, double> m;
      for (auto &kv : rowm) m[kv.first] = kv.second.first / kv.second.second;
      for (auto &a : m)
        for (auto &b : m)
        {
          int d = b.first - a.first;
          if (d < 0 || d > 30) continue;
          AC[s][d] += a.second * b.second;
          ACn[s][d]++;
        }
      // (4) d0 + split-half curvature
      double d0 = std::fabs(std::hypot(FT.a, FT.b) - FT.R);
      hd0[s]->Fill(std::min(d0, 11.9));
      d0v[s].push_back(d0);
      std::vector<double> Xi, Yi, Ri, Xo, Yo, Ro;
      for (int i : kept)
      {
        int L = nearRow(G.r[i]);
        if (L < 0) continue;
        if (L <= 30) { Xi.push_back(G.x[i]); Yi.push_back(G.y[i]); Ri.push_back(G.r[i]); }
        else { Xo.push_back(G.x[i]); Yo.push_back(G.y[i]); Ro.push_back(G.r[i]); }
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
          hds[s]->Fill(std::max(-12.49, std::min(12.49, ds)));
          dsv[s].push_back(ds);
          double sfull = signedK(Xi, Yi, Ri, FT);   // full-fit bending sign (inner points, trimmed fit)
          dsq[s][sfull > 0 ? 1 : 0].push_back(ds);
        }
      }
    }

  // ---- ledger ----------------------------------------------------------------
  FILE *fo = fopen(Form("nonrms_probe_%s.txt", ver), "w");
  auto P = [&](const char *fmt, ...) {
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    va_start(ap, fmt); vfprintf(fo, fmt, ap); va_end(ap);
  };
  P("[nonrms_probe %s] non-RMS metric battery on the nf_digipix groups (real road / sim truth)\n", ver);
  for (int s = 0; s < 2; ++s)
  {
    // pooled residual quantiles
    double q68 = 0, q95 = 0, q99 = 0;
    {
      double qs[3] = {0.68, 0.95, 0.99}, xs[3];
      TH1D ha("ha", "", 3000, 0, 15);
      for (int b = 1; b <= 120; ++b)
      {
        double c = hres[s]->GetBinContent(b), r = std::fabs(hres[s]->GetBinCenter(b));
        ha.Fill(std::min(r, 14.9), c);
      }
      ha.GetQuantiles(3, xs, qs);
      q68 = xs[0]; q95 = xs[1]; q99 = xs[2];
    }
    double sg = 0;
    { std::vector<double> v = skew[s]; sg = medv(v); }
    P("  %s: %ld tracks, %ld pixels pooled\n", s ? "sim " : "real", ntrk2[s], npix[s]);
    P("    (1) |residual| quantiles q68/q95/q99: %.2f / %.2f / %.2f mm | frac |res|>3mm per track med %.1f%% | per-track skewness med %+.2f\n",
      q68, q95, q99, 100 * medv(tailr[s]), sg);
    P("    (2) pull tails: frac |pull|>3 per track med %.1f%%\n", 100 * medv(pull3[s]));
    P("    (3) row-lag autocorrelation C(d)/C(0): d=1 %.2f | 2 %.2f | 4 %.2f | 8 %.2f | 16 %.2f | 24 %.2f\n",
      AC[s][1] / ACn[s][1] / (AC[s][0] / ACn[s][0]), AC[s][2] / ACn[s][2] / (AC[s][0] / ACn[s][0]),
      AC[s][4] / ACn[s][4] / (AC[s][0] / ACn[s][0]), AC[s][8] / ACn[s][8] / (AC[s][0] / ACn[s][0]),
      AC[s][16] / ACn[s][16] / (AC[s][0] / ACn[s][0]), AC[s][24] / ACn[s][24] / (AC[s][0] / ACn[s][0]));
    std::vector<double> dv = dsv[s], tmp2 = dv;
    double mdd = medv(tmp2);
    for (double &q : tmp2) q = std::fabs(q - mdd);
    P("    (4) |d0| med %.2f cm (q90 %.2f) | R med %.0f cm | split-half Dsagitta: robust sigma %.2f mm, med %+.2f (n=%zu)\n",
      medv(d0v[s]), [&]{ std::vector<double> v = d0v[s]; std::sort(v.begin(), v.end()); return v.empty() ? 0. : v[(size_t)(0.9 * v.size())]; }(),
      medv(Rv[s]), 1.4826 * medv(tmp2), mdd, dsv[s].size());
    P("        charge-resolved Dsagitta med: bend<0 %+.2f mm (n=%zu) | bend>0 %+.2f mm (n=%zu)  [same sign both = coherent twist; opposite = |k| effect]\n",
      medv(dsq[s][0]), dsq[s][0].size(), medv(dsq[s][1]), dsq[s][1].size());
  }
  fclose(fo);
  printf("wrote nonrms_probe_%s.txt\n", ver);

  // ---- figure ----------------------------------------------------------------
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.045);
  TCanvas *cv = new TCanvas("cvnr", "nonrms", 1700, 1250);
  cv->Divide(2, 2);
  TLatex tx; tx.SetNDC();
  auto draw2 = [&](int ipad, TH1D *h0, TH1D *h1, bool logy, const char *a0, const char *a1) {
    cv->cd(ipad);
    if (logy) gPad->SetLogy();
    for (auto h : {h0, h1}) if (h->Integral() > 0) h->Scale(1. / h->Integral());
    h0->SetLineColor(kBlack); h0->SetLineWidth(2);
    h1->SetLineColor(kBlue + 1); h1->SetLineWidth(2); h1->SetLineStyle(2);
    h0->SetMaximum((logy ? 6 : 1.45) * std::max(h0->GetMaximum(), h1->GetMaximum()));
    if (logy) h0->SetMinimum(2e-7);
    h0->Draw("hist"); h1->Draw("hist same");
    TLegend *L = new TLegend(0.60, 0.74, 0.89, 0.87); L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextSize(0.030);
    L->AddEntry(h0, "real", "l"); L->AddEntry(h1, "sim", "l");
    L->Draw();
    tx.SetTextSize(0.028);
    if (a0[0]) tx.DrawLatex(0.14, 0.84, a0);
    if (a1[0]) tx.DrawLatex(0.14, 0.79, a1);
  };
  draw2(1, hres[0], hres[1], true, "same core; real mm-scale shoulders = road pickup", "");
  draw2(2, hpull[0], hpull[1], true, "scale-free: core ~unit width both sides;", "tails = association, not response");
  // autocorrelation panel
  cv->cd(3);
  {
    TH2D *fr = new TH2D("nr_fr3", "row-lag autocorrelation of row-mean residuals;row lag d;C(d) / C(0)", 10, 0, 26, 10, -0.25, 1.05);
    fr->Draw();
    TLine z; z.SetLineStyle(2); z.SetLineColor(kGray + 2); z.DrawLine(0, 0, 26, 0);
    TGraph *g0 = new TGraph, *g1 = new TGraph;
    for (int d = 0; d <= 25; ++d)
    {
      if (ACn[0][d] > 0) g0->AddPoint(d, AC[0][d] / ACn[0][d] / (AC[0][0] / ACn[0][0]));
      if (ACn[1][d] > 0) g1->AddPoint(d, AC[1][d] / ACn[1][d] / (AC[1][0] / ACn[1][0]));
    }
    g0->SetLineColor(kBlack); g0->SetLineWidth(2); g0->SetMarkerStyle(20);
    g1->SetLineColor(kBlue + 1); g1->SetLineWidth(2); g1->SetMarkerStyle(24); g1->SetMarkerColor(kBlue + 1); g1->SetLineStyle(2);
    g0->Draw("LP same"); g1->Draw("LP same");
    TLegend *L = new TLegend(0.60, 0.74, 0.89, 0.87); L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextSize(0.030);
    L->AddEntry(g0, "real", "lp"); L->AddEntry(g1, "sim", "lp");
    L->Draw();
    tx.SetTextSize(0.028);
    tx.DrawLatex(0.30, 0.60, "positive long-range C(d) = coherent bending");
    tx.DrawLatex(0.30, 0.55, "(field/alignment), invisible to any RMS");
  }
  draw2(4, hds[0], hds[1], true, "parameter-level: curvature agreement of the two", "track halves - blind to the cloud floor");
  cv->SaveAs(Form("../sim_validation_plots/nonrms_probe_%s.png", ver));
  printf("wrote ../sim_validation_plots/nonrms_probe_%s.png\n", ver);
}
