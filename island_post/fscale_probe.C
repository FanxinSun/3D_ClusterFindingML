// fscale_probe.C — the rscale_probe before/after-sagitta metrics REDONE WITH
// THE EXHAUSTIVE FINDER (user, 2026-08-20). The no-finder branch groups real
// pixels by the tracker road and sim by truth (asymmetric by design); here
// BOTH sides are grouped by the SAME algorithm: MTK::hunt (missed_tracks.C:
// sectorized conformal Hough + drift-coherence RANSAC + circularity bar
// >=12 clusters / >=13 layers / span>=15 / RMS<=0.20 cm / R>=45 / maxgap<=6),
// run on ALL clusters (real ntp_cluster, sim island91 v6), then pixels are
// attached to each found track's member clusters with the nf_digipix road
// (same event+layer, dxy < 1.2 cm, |dtbin| <= 6, first-track-wins).
// Metrics identical to rscale_probe: RMS(L) raw/3xMAD-clipped, paired
// per-track g/l, removed content sqrt(g^2-l^2). The no-finder real raw curve
// is overlaid from rscale_probe_<ver>.txt for direct comparison.
// usage: root -l -b -q 'fscale_probe.C+()'
// out: fscale_probe_<ver>.txt + ../sim_validation_plots/fscale_probe_<ver>.png
#include "../sim_validation_plots/src/missed_tracks.C"
#include <TROOT.h>
#include <TH2D.h>
#include <TLine.h>
#include <unordered_map>
#include <cstdint>

namespace FSP
{
double medv(std::vector<double> v)
{
  if (v.empty()) return 0;
  std::sort(v.begin(), v.end());
  return v[v.size() / 2];
}
struct GrpX { std::vector<double> x, y, r; };
double clippedRMS(const std::vector<double> &X, const std::vector<double> &Y)
{
  std::vector<int> kept(X.size());
  for (size_t k = 0; k < kept.size(); ++k) kept[k] = (int) k;
  MTK::Fit F;
  for (int it = 0; it < 3; ++it)
  {
    std::vector<double> xs, ys;
    for (int i : kept) { xs.push_back(X[i]); ys.push_back(Y[i]); }
    F = MTK::fitCircle(xs, ys);
    if (!F.ok) return -1;
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
    if ((int) nk.size() < 5) break;
    kept.swap(nk);
  }
  std::vector<double> xs, ys;
  for (int i : kept) { xs.push_back(X[i]); ys.push_back(Y[i]); }
  F = MTK::fitCircle(xs, ys);
  return F.ok ? F.rms : -1;
}
}  // namespace FSP

void fscale_probe(const char *i91 = "island91_frames_production_v6.root",
                  const char *digif = "digi_frames_production_v6.root", int nsim = 50,
                  const char *realf = "../clusters_seeds_island_79507-0.root_ntuplizer.root",
                  const char *ver = "v6")
{
  using namespace FSP;
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

  // ---- clusters + finder ----------------------------------------------------
  // member clusters per found track, keyed for the pixel road:
  // bucket[(ev, lay)] -> list of (clx, cly, cltb, trackIndex) in track order
  std::vector<GrpX> gpx[2];                        // pixel groups per found track
  long nfound[2] = {0, 0};
  std::vector<double> fmedtb[2];
  for (int s = 0; s < 2; ++s)
  {
    std::map<int, std::vector<MTK::Cl>> ev;
    if (s == 0)
    {
      TFile *f = TFile::Open(realf);
      TTree *c = (TTree *) f->Get("ntp_cluster");
      float evn, lay, x, y, tb;
      c->SetBranchStatus("*", 0);
      for (auto b : {"event", "layer", "x", "y", "tbin"}) c->SetBranchStatus(b, 1);
      c->SetBranchAddress("event", &evn); c->SetBranchAddress("layer", &lay);
      c->SetBranchAddress("x", &x); c->SetBranchAddress("y", &y); c->SetBranchAddress("tbin", &tb);
      for (Long64_t i = 0; i < c->GetEntries(); ++i)
      {
        c->GetEntry(i);
        if ((int) evn == 44) continue;             // V6 laser veto (canon.h)
        if (lay < 7 || lay > 54) continue;
        ev[(int) evn].push_back({x, y, lay, tb, -1});
      }
      f->Close();
    }
    else
    {
      TFile *f = TFile::Open(i91);
      TTree *c = (TTree *) f->Get("ntp_cluster");
      float evn, lay, x, y, tb;
      c->SetBranchStatus("*", 0);
      for (auto b : {"event", "layer", "x", "y", "tbin"}) c->SetBranchStatus(b, 1);
      c->SetBranchAddress("event", &evn); c->SetBranchAddress("layer", &lay);
      c->SetBranchAddress("x", &x); c->SetBranchAddress("y", &y); c->SetBranchAddress("tbin", &tb);
      for (Long64_t i = 0; i < c->GetEntries(); ++i)
      {
        c->GetEntry(i);
        if ((int) evn >= nsim) continue;
        if (lay < 7 || lay > 54) continue;
        ev[(int) evn].push_back({x, y, lay, tb, -1});
      }
      f->Close();
    }
    // run the finder, build road buckets
    struct RC { float x, y, tb; int trk; };
    std::unordered_map<long, std::vector<RC>> bucket;   // ev*100+lay
    long ntrk = 0;
    for (auto &kv : ev)
    {
      std::vector<MTK::Trk> t = MTK::hunt(kv.second);
      for (auto &T : t)
      {
        int gi = (int) gpx[s].size();
        gpx[s].push_back(GrpX());
        fmedtb[s].push_back(T.medtb);
        for (int i : T.idx)
        {
          const MTK::Cl &C = kv.second[i];
          bucket[(long) kv.first * 100 + (long) C.lay].push_back({C.x, C.y, C.tb, gi});
        }
        ntrk++;
      }
    }
    nfound[s] = ntrk;
    printf("%s: finder tracks %ld (events %zu)\n", s ? "sim" : "real", ntrk, ev.size());
    // ---- pixel road ---------------------------------------------------------
    if (s == 0)
    {
      TFile *f = TFile::Open(realf);
      TTree *t = (TTree *) f->Get("ntp_hit");
      float evn, lay, x, y, tb, adc;
      t->SetBranchStatus("*", 0);
      for (auto b : {"event", "layer", "x", "y", "tbin", "adc"}) t->SetBranchStatus(b, 1);
      t->SetBranchAddress("event", &evn); t->SetBranchAddress("layer", &lay);
      t->SetBranchAddress("x", &x); t->SetBranchAddress("y", &y);
      t->SetBranchAddress("tbin", &tb); t->SetBranchAddress("adc", &adc);
      for (Long64_t i = 0; i < t->GetEntries(); ++i)
      {
        t->GetEntry(i);
        if ((int) evn == 44) continue;
        if (lay < 7 || lay > 54 || adc <= 0) continue;
        auto bit = bucket.find((long) evn * 100 + (long) lay);
        if (bit == bucket.end()) continue;
        for (const RC &c : bit->second)
        {
          if (std::fabs(tb - c.tb) > 6) continue;
          double dx = x - c.x, dy = y - c.y;
          if (dx * dx + dy * dy > 1.2 * 1.2) continue;
          GrpX &G = gpx[0][c.trk];
          G.x.push_back(x); G.y.push_back(y); G.r.push_back(std::hypot(x, y));
          break;                                   // first track wins
        }
      }
      f->Close();
    }
    else
    {
      TFile *f = TFile::Open(digif);
      TTree *t = (TTree *) f->Get("ntp_hit");
      float evn, lay, phi, adc, tb;
      t->SetBranchStatus("*", 0);
      for (auto b : {"event", "layer", "phi", "adc", "tbin"}) t->SetBranchStatus(b, 1);
      t->SetBranchAddress("event", &evn); t->SetBranchAddress("layer", &lay);
      t->SetBranchAddress("phi", &phi); t->SetBranchAddress("adc", &adc);
      t->SetBranchAddress("tbin", &tb);
      for (Long64_t i = 0; i < t->GetEntries(); ++i)
      {
        t->GetEntry(i);
        if ((int) evn >= nsim) continue;
        if (lay < 7 || lay > 54 || adc <= 0) continue;
        auto bit = bucket.find((long) evn * 100 + (long) lay);
        if (bit == bucket.end()) continue;
        double r = rowR[(int) lay], x = r * std::cos(phi), y = r * std::sin(phi);
        for (const RC &c : bit->second)
        {
          if (std::fabs(tb - c.tb) > 6) continue;
          double dx = x - c.x, dy = y - c.y;
          if (dx * dx + dy * dy > 1.2 * 1.2) continue;
          GrpX &G = gpx[1][c.trk];
          G.x.push_back(x); G.y.push_back(y); G.r.push_back(r);
          break;
        }
      }
      f->Close();
    }
  }

  // ---- scale scan + paired (identical to rscale_probe) ----------------------
  const int NL = 8;
  int LS[NL] = {4, 6, 8, 12, 16, 24, 32, 48};
  std::vector<double> rmsL[2][NL], rmsLc[2][NL];
  std::vector<double> pg[2], pl[2], prat[2], prem[2];
  long nbar[2] = {0, 0};
  for (int s = 0; s < 2; ++s)
    for (auto &G : gpx[s])
    {
      if ((int) G.x.size() < 12) continue;
      double rlo = 1e9, rhi = 0;
      for (double q : G.r) { rlo = std::min(rlo, q); rhi = std::max(rhi, q); }
      if (rhi - rlo < 15) continue;
      MTK::Fit F0 = MTK::fitCircle(G.x, G.y);
      if (!F0.ok || F0.R < 45 || F0.R >= 2e4) continue;
      nbar[s]++;
      std::map<int, std::vector<int>> rows;
      for (size_t i = 0; i < G.x.size(); ++i)
      { int L = nearRow(G.r[i]); if (L >= 0) rows[L].push_back((int) i); }
      std::vector<double> loc4;
      for (int il = 0; il < NL; ++il)
      {
        int L = LS[il], step = (L == 4) ? 1 : std::max(1, L / 2);
        for (int w = 7; w + L - 1 <= 54; w += step)
        {
          std::vector<double> X, Y;
          int nr = 0;
          for (int rr = w; rr < w + L; ++rr)
          {
            auto it = rows.find(rr);
            if (it == rows.end()) continue;
            nr++;
            for (int i : it->second) { X.push_back(G.x[i]); Y.push_back(G.y[i]); }
          }
          if (nr < std::max(3, L / 2) || (int) X.size() < 5) continue;
          MTK::Fit F = MTK::fitCircle(X, Y);
          if (!F.ok) continue;
          rmsL[s][il].push_back(F.rms * 1e4);
          if (il == 0) loc4.push_back(F.rms * 1e4);
          double rc = clippedRMS(X, Y);
          if (rc >= 0) rmsLc[s][il].push_back(rc * 1e4);
        }
      }
      if (!loc4.empty())
      {
        double g = F0.rms * 1e4, l = medv(loc4);
        pg[s].push_back(g); pl[s].push_back(l);
        prat[s].push_back(g / l);
        prem[s].push_back(std::sqrt(std::max(0., g * g - l * l)));
      }
    }

  // ---- ledger ---------------------------------------------------------------
  FILE *fo = fopen(Form("fscale_probe_%s.txt", ver), "w");
  auto P = [&](const char *fmt, ...) {
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    va_start(ap, fmt); vfprintf(fo, fmt, ap); va_end(ap);
  };
  auto M = [&](std::vector<double> &v) { return medv(v); };
  P("[fscale_probe %s] before/after the 4-row sagitta fit WITH the exhaustive finder\n", ver);
  P("  grouping: MTK::hunt on ALL clusters both sides (real ntp_cluster, all 99 events; sim island91, %d frames),\n", nsim);
  P("  pixels attached by the nf_digipix road to member clusters; fit bar >=12 px, span>=15, 45<=R<2e4\n");
  P("  finder tracks: real %ld -> %ld pass pixel bar | sim %ld -> %ld\n", nfound[0], nbar[0], nfound[1], nbar[1]);
  P("  RMS(L): median per-fit RMS vs window length; clip = per-window 3xMAD\n");
  P("   L    real raw  real clip  sim raw  sim clip   (windows real/sim)\n");
  for (int il = 0; il < NL; ++il)
    P("  %2d    %6.0f    %6.0f    %6.0f   %6.0f     (%zu/%zu)\n", LS[il],
      M(rmsL[0][il]), M(rmsLc[0][il]), M(rmsL[1][il]), M(rmsLc[1][il]),
      rmsL[0][il].size(), rmsL[1][il].size());
  for (int s = 0; s < 2; ++s)
    P("  %s PAIRED (%zu tracks): global med %.0f | own-4-row med %.0f | g/l med %.2f (q90 %.2f) | removed med %.0f um\n",
      s ? "sim " : "real", pg[s].size(), M(pg[s]), M(pl[s]), M(prat[s]),
      [&]{ std::vector<double> v = prat[s]; std::sort(v.begin(), v.end()); return v.empty() ? 0. : v[(size_t)(0.9 * v.size())]; }(),
      M(prem[s]));
  P("  no-finder reference (rscale_probe_%s.txt): real raw 1205->2330, real clip 1165->1587, g/l 1.82, removed 1930\n", ver);
  fclose(fo);
  printf("wrote fscale_probe_%s.txt\n", ver);

  // ---- figure ----------------------------------------------------------------
  // no-finder real curves for overlay, parsed from rscale_probe ledger
  std::vector<std::pair<int, double>> nfraw, nfclip;
  {
    FILE *fr = fopen(Form("rscale_probe_%s.txt", ver), "r");
    if (fr)
    {
      char line[512];
      while (fgets(line, sizeof line, fr))
      {
        int L; double a, b, c, d;
        if (sscanf(line, " %d %lf %lf %lf %lf", &L, &a, &b, &c, &d) == 5 && L >= 4 && L <= 48)
        { nfraw.push_back({L, a}); nfclip.push_back({L, b}); }
      }
      fclose(fr);
    }
  }
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.045);
  TCanvas *cv = new TCanvas("cvfs", "fscale", 1700, 1250);
  cv->Divide(2, 2);
  TLatex tx; tx.SetNDC();
  // [1] RMS(L) with finder + no-finder overlay
  cv->cd(1);
  {
    TH2D *fr = new TH2D("fs1", "median RMS vs window length, WITH exhaustive finder;window length [pad rows];median per-fit RMS [#mum]", 10, 2, 52, 10, 0, 2800);
    fr->Draw();
    TGraph *g0 = new TGraph, *g1 = new TGraph, *g2 = new TGraph, *g3 = new TGraph, *g4 = new TGraph;
    for (int il = 0; il < NL; ++il)
    {
      g0->AddPoint(LS[il], M(rmsL[0][il])); g1->AddPoint(LS[il], M(rmsLc[0][il]));
      g2->AddPoint(LS[il], M(rmsL[1][il])); g3->AddPoint(LS[il], M(rmsLc[1][il]));
    }
    for (auto &q : nfraw) g4->AddPoint(q.first, q.second);
    for (auto g : {g0, g1, g2, g3}) g->SetLineWidth(2);
    g0->SetLineColor(kBlack); g0->SetMarkerStyle(20); g0->SetMarkerColor(kBlack);
    g1->SetLineColor(kRed + 1); g1->SetLineStyle(7); g1->SetMarkerStyle(24); g1->SetMarkerColor(kRed + 1);
    g2->SetLineColor(kBlue + 1); g2->SetMarkerStyle(21); g2->SetMarkerColor(kBlue + 1);
    g3->SetLineColor(kAzure + 6); g3->SetLineStyle(3); g3->SetMarkerStyle(25); g3->SetMarkerColor(kAzure + 6);
    g4->SetLineColor(kGray + 2); g4->SetLineWidth(3); g4->SetLineStyle(2);
    if (g4->GetN()) g4->Draw("L same");
    for (auto g : {g2, g3, g1, g0}) g->Draw("LP same");
    TLegend *L = new TLegend(0.14, 0.60, 0.66, 0.88); L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextSize(0.029);
    L->AddEntry(g0, "real raw (finder groups)", "lp");
    L->AddEntry(g1, "real, per-window 3#sigma clip", "lp");
    L->AddEntry(g2, "sim raw (finder groups)", "lp");
    L->AddEntry(g3, "sim, same clip", "lp");
    if (g4->GetN()) L->AddEntry(g4, "no-finder real raw (rscale_probe)", "l");
    L->Draw();
  }
  // [2] paired scatter (real)
  cv->cd(2);
  {
    TH2D *h = new TH2D("fs2", "REAL finder tracks, paired: before vs after;GLOBAL whole-track RMS [#mum];median of own 4-row windows [#mum]", 60, 0, 6000, 60, 0, 3000);
    for (size_t i = 0; i < pg[0].size(); ++i) h->Fill(std::min(pg[0][i], 5999.), std::min(pl[0][i], 2999.));
    h->Draw("colz");
    gPad->SetRightMargin(0.12);
    TLine d; d.SetLineStyle(2); d.SetLineColor(kGray + 2); d.DrawLine(0, 0, 3000, 3000);
    tx.SetTextSize(0.030);
    tx.DrawLatex(0.16, 0.84, "diagonal: no reduction");
  }
  // [3] reduction factor
  cv->cd(3);
  {
    TH1D *h0 = new TH1D("fs3_0", "reduction factor per track: global / own-local median;g / l;tracks (unit area)", 60, 0.5, 6.5);
    TH1D *h1 = (TH1D *) h0->Clone("fs3_1");
    for (double q : prat[0]) h0->Fill(std::min(q, 6.49));
    for (double q : prat[1]) h1->Fill(std::min(q, 6.49));
    for (auto h : {h0, h1}) if (h->Integral() > 0) h->Scale(1. / h->Integral());
    h0->SetLineColor(kBlack); h0->SetLineWidth(2);
    h1->SetLineColor(kBlue + 1); h1->SetLineWidth(2); h1->SetLineStyle(2);
    h0->SetMaximum(1.45 * std::max(h0->GetMaximum(), h1->GetMaximum()));
    h0->Draw("hist"); h1->Draw("hist same");
    TLegend *L = new TLegend(0.42, 0.70, 0.89, 0.87); L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextSize(0.030);
    L->AddEntry(h0, Form("real, med %.2f", M(prat[0])), "l");
    L->AddEntry(h1, Form("sim, med %.2f", M(prat[1])), "l");
    L->Draw();
    tx.SetTextSize(0.029);
    tx.DrawLatex(0.42, 0.64, "no-finder real med was 1.82");
  }
  // [4] removed content
  cv->cd(4);
  {
    TH1D *h0 = new TH1D("fs4_0", "removed content per track: #sqrt{g^{2}#minus l^{2}};[#mum];tracks (unit area)", 60, 0, 6000);
    TH1D *h1 = (TH1D *) h0->Clone("fs4_1");
    for (double q : prem[0]) h0->Fill(std::min(q, 5999.));
    for (double q : prem[1]) h1->Fill(std::min(q, 5999.));
    for (auto h : {h0, h1}) if (h->Integral() > 0) h->Scale(1. / h->Integral());
    h0->SetLineColor(kBlack); h0->SetLineWidth(2);
    h1->SetLineColor(kBlue + 1); h1->SetLineWidth(2); h1->SetLineStyle(2);
    h0->SetMaximum(1.45 * std::max(h0->GetMaximum(), h1->GetMaximum()));
    h0->Draw("hist"); h1->Draw("hist same");
    TLegend *L = new TLegend(0.42, 0.70, 0.89, 0.87); L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextSize(0.030);
    L->AddEntry(h0, Form("real, med %.0f #mum", M(prem[0])), "l");
    L->AddEntry(h1, Form("sim, med %.0f #mum", M(prem[1])), "l");
    L->Draw();
    tx.SetTextSize(0.029);
    tx.DrawLatex(0.42, 0.64, "no-finder real med was 1930 #mum");
  }
  cv->SaveAs(Form("../sim_validation_plots/fscale_probe_%s.png", ver));
  printf("wrote ../sim_validation_plots/fscale_probe_%s.png\n", ver);
}
