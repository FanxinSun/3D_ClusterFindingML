// rscale_probe.C — metrics for "before vs after the 4-row sagitta fit" beyond
// one RMS ratio (user concern 2026-08-20). The nf_digipix before/after compares
// two UNPAIRED medians in different units (4,836 tracks vs 134,742 windows).
// This probe measures the difference properly:
//   (1) RMS(L): median per-fit RMS vs window length L = 4..48 rows (raw and
//       per-window 3xMAD-clipped, real and sim) — the SHAPE of the growth
//       curve says what enters at which scale (cloud: flat; field/curvature:
//       smooth rise; road pickup: rise that the clip removes).
//   (2) PAIRED per-track: global RMS vs the median of that track's own 4-row
//       windows; reduction factor g/l per track; removed content
//       sqrt(g^2 - l^2) per track.
// Same groups / fitter / gates as nf_digipix (road, ev44 veto, fit bar);
// window gate at every L: >= max(3, L/2) distinct rows, >= 5 px.
// usage: root -l -b -q 'rscale_probe.C+()'
// out: rscale_probe_<ver>.txt + ../sim_validation_plots/rscale_probe_<ver>.png
#include "ms_nofinder.C"
#include <TGraph.h>
#include <TH2D.h>
#include <TLine.h>
#include <cstdint>

namespace RSP
{
struct GrpX
{
  std::vector<double> x, y, r;
  int ev = -1;
};
// per-window robust clip (same spec as gtail_probe::trimFit, min 5 px)
double clippedRMS(const std::vector<double> &X, const std::vector<double> &Y)
{
  std::vector<int> kept(X.size());
  for (size_t k = 0; k < kept.size(); ++k) kept[k] = (int) k;
  MNF::Fit F;
  for (int it = 0; it < 3; ++it)
  {
    std::vector<double> xs, ys;
    for (int i : kept) { xs.push_back(X[i]); ys.push_back(Y[i]); }
    F = MNF::fitCircle(xs, ys);
    if (!F.ok) return -1;
    std::vector<double> res;
    for (int i : kept) res.push_back(std::hypot(X[i] - F.a, Y[i] - F.b) - F.R);
    std::vector<double> tmp = res;
    double md = MNF::med(tmp);
    for (double &q : tmp) q = std::fabs(q - md);
    double thr = std::max(3 * 1.4826 * MNF::med(tmp), 0.05);
    std::vector<int> nk;
    for (size_t k = 0; k < kept.size(); ++k)
      if (std::fabs(res[k] - md) <= thr) nk.push_back(kept[k]);
    if (nk.size() == kept.size()) break;
    if ((int) nk.size() < 5) break;
    kept.swap(nk);
  }
  std::vector<double> xs, ys;
  for (int i : kept) { xs.push_back(X[i]); ys.push_back(Y[i]); }
  F = MNF::fitCircle(xs, ys);
  return F.ok ? F.rms : -1;
}
}  // namespace RSP

void rscale_probe(const char *digif = "digi_frames_production_v6.root", int nsim = 60,
                  const char *realf = "../clusters_seeds_island_79507-0.root_ntuplizer.root",
                  const char *ver = "v6")
{
  using namespace RSP;
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

  // ---- groups (verbatim nf_digipix construction) ---------------------------
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
    for (auto &kv : seedIdx) gr[0][kv.second].ev = kv.first.first;
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
    if (!f || f->IsZombie()) { printf("no %s\n", digif); return; }
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
      if (it == gidx.end()) { it = gidx.insert({key, (int) gr[1].size()}).first; gr[1].push_back(GrpX()); gr[1].back().ev = (int) ev; }
      GrpX &G = gr[1][it->second];
      double r = rowR[(int) lay];
      G.x.push_back(r * std::cos(phi)); G.y.push_back(r * std::sin(phi)); G.r.push_back(r);
    }
    f->Close();
    printf("sim: %zu truth groups\n", gr[1].size());
  }

  // ---- scale scan + paired per-track --------------------------------------
  const int NL = 8;
  int LS[NL] = {4, 6, 8, 12, 16, 24, 32, 48};
  std::vector<double> rmsL[2][NL], rmsLc[2][NL];       // per-window, um
  std::vector<double> pg[2], pl[2], prat[2], prem[2];  // paired per-track
  for (int s = 0; s < 2; ++s)
    for (auto &G : gr[s])
    {
      MNF::Grp B; B.x = G.x; B.y = G.y; B.r = G.r;
      MNF::Fit F0;
      if (!MNF::fitBar(B, F0)) continue;
      // rows -> pixel lists
      std::map<int, std::vector<int>> rows;
      for (size_t i = 0; i < G.x.size(); ++i)
      { int L = nearRow(G.r[i]); if (L >= 0) rows[L].push_back((int) i); }
      std::vector<double> loc4;                        // this track's 4-row window RMS
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
          MNF::Fit F = MNF::fitCircle(X, Y);
          if (!F.ok) continue;
          rmsL[s][il].push_back(F.rms * 1e4);
          if (il == 0) loc4.push_back(F.rms * 1e4);
          double rc = clippedRMS(X, Y);
          if (rc >= 0) rmsLc[s][il].push_back(rc * 1e4);
        }
      }
      if (!loc4.empty())
      {
        double g = F0.rms * 1e4, l = MNF::med(loc4);
        pg[s].push_back(g); pl[s].push_back(l);
        prat[s].push_back(g / l);
        prem[s].push_back(std::sqrt(std::max(0., g * g - l * l)));
      }
    }

  // ---- ledger ---------------------------------------------------------------
  FILE *fo = fopen(Form("rscale_probe_%s.txt", ver), "w");
  auto P = [&](const char *fmt, ...) {
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    va_start(ap, fmt); vfprintf(fo, fmt, ap); va_end(ap);
  };
  auto M = [&](std::vector<double> &v) { return MNF::med(v); };
  P("[rscale_probe %s] before/after the 4-row sagitta fit, beyond one RMS ratio\n", ver);
  P("  RMS(L): median per-fit RMS vs window length (rows); clip = per-window 3xMAD\n");
  P("   L    real raw  real clip  sim raw  sim clip   (windows real/sim)\n");
  for (int il = 0; il < NL; ++il)
    P("  %2d    %6.0f    %6.0f    %6.0f   %6.0f     (%zu/%zu)\n", LS[il],
      M(rmsL[0][il]), M(rmsLc[0][il]), M(rmsL[1][il]), M(rmsLc[1][il]),
      rmsL[0][il].size(), rmsL[1][il].size());
  for (int s = 0; s < 2; ++s)
  {
    long nge = 0;
    for (double q : prat[s]) if (q < 1.05) nge++;
    P("  %s PAIRED (%zu tracks): global med %.0f | own-4-row-median med %.0f | ratio g/l med %.2f "
      "(q90 %.2f) | removed sqrt(g^2-l^2) med %.0f um | tracks with no reduction (g<1.05*l): %.1f%%\n",
      s ? "sim " : "real", pg[s].size(), M(pg[s]), M(pl[s]), M(prat[s]),
      [&]{ std::vector<double> v = prat[s]; std::sort(v.begin(), v.end()); return v[(size_t)(0.9 * v.size())]; }(),
      M(prem[s]), 100. * nge / prat[s].size());
  }
  fclose(fo);
  printf("wrote rscale_probe_%s.txt\n", ver);

  // ---- figure ----------------------------------------------------------------
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.045);
  TCanvas *cv = new TCanvas("cvrs", "rscale", 1700, 1250);
  cv->Divide(2, 2);
  TLatex tx; tx.SetNDC();
  // [1] RMS(L)
  cv->cd(1);
  {
    TGraph *g0 = new TGraph, *g1 = new TGraph, *g2 = new TGraph, *g3 = new TGraph;
    for (int il = 0; il < NL; ++il)
    {
      g0->AddPoint(LS[il], M(rmsL[0][il])); g1->AddPoint(LS[il], M(rmsLc[0][il]));
      g2->AddPoint(LS[il], M(rmsL[1][il])); g3->AddPoint(LS[il], M(rmsLc[1][il]));
    }
    TH2D *fr = new TH2D("rs1", "median RMS vs fit window length;window length [pad rows];median per-fit RMS [#mum]", 10, 2, 52, 10, 0, 2800);
    fr->Draw();
    for (auto g : {g0, g1, g2, g3}) g->SetLineWidth(2);
    g0->SetLineColor(kBlack); g0->SetMarkerStyle(20); g0->SetMarkerColor(kBlack);
    g1->SetLineColor(kRed + 1); g1->SetLineStyle(7); g1->SetMarkerStyle(24); g1->SetMarkerColor(kRed + 1);
    g2->SetLineColor(kBlue + 1); g2->SetMarkerStyle(21); g2->SetMarkerColor(kBlue + 1);
    g3->SetLineColor(kAzure + 6); g3->SetLineStyle(3); g3->SetMarkerStyle(25); g3->SetMarkerColor(kAzure + 6);
    for (auto g : {g2, g3, g1, g0}) g->Draw("LP same");
    TLegend *L = new TLegend(0.14, 0.64, 0.62, 0.88); L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextSize(0.030);
    L->AddEntry(g0, "real raw", "lp");
    L->AddEntry(g1, "real, per-window 3#sigma clip", "lp");
    L->AddEntry(g2, "sim raw", "lp");
    L->AddEntry(g3, "sim, same clip", "lp");
    L->Draw();
    tx.SetTextSize(0.029);
    tx.DrawLatex(0.40, 0.30, "raw#minusclip gap = road pickup (grows with L);");
    tx.DrawLatex(0.40, 0.25, "clip#minussim gap = field/curvature share");
  }
  // [2] paired scatter (real)
  cv->cd(2);
  {
    TH2D *h = new TH2D("rs2", "REAL, paired per track: before vs after;GLOBAL whole-track RMS [#mum];median of own 4-row windows [#mum]", 60, 0, 6000, 60, 0, 3000);
    for (size_t i = 0; i < pg[0].size(); ++i) h->Fill(std::min(pg[0][i], 5999.), std::min(pl[0][i], 2999.));
    h->Draw("colz");
    gPad->SetRightMargin(0.12);
    TLine d; d.SetLineStyle(2); d.SetLineColor(kGray + 2); d.DrawLine(0, 0, 3000, 3000);
    tx.SetTextSize(0.030);
    tx.DrawLatex(0.16, 0.84, "diagonal: no reduction");
    tx.DrawLatex(0.16, 0.79, "horizontal band: local floor #approx cloud width");
  }
  // [3] reduction factor
  cv->cd(3);
  {
    TH1D *h0 = new TH1D("rs3_0", "reduction factor per track: global / own-local median;g / l;tracks (unit area)", 60, 0.5, 6.5);
    TH1D *h1 = (TH1D *) h0->Clone("rs3_1");
    for (double q : prat[0]) h0->Fill(std::min(q, 6.49));
    for (double q : prat[1]) h1->Fill(std::min(q, 6.49));
    for (auto h : {h0, h1}) h->Scale(1. / h->Integral());
    h0->SetLineColor(kBlack); h0->SetLineWidth(2);
    h1->SetLineColor(kBlue + 1); h1->SetLineWidth(2); h1->SetLineStyle(2);
    h0->SetMaximum(1.45 * std::max(h0->GetMaximum(), h1->GetMaximum()));
    h0->Draw("hist"); h1->Draw("hist same");
    TLegend *L = new TLegend(0.42, 0.70, 0.89, 0.87); L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextSize(0.030);
    L->AddEntry(h0, Form("real, med %.2f", M(prat[0])), "l");
    L->AddEntry(h1, Form("sim, med %.2f", M(prat[1])), "l");
    L->Draw();
    tx.SetTextSize(0.029);
    tx.DrawLatex(0.42, 0.64, "real tail (g/l > 2.5) = road-contaminated tracks;");
    tx.DrawLatex(0.42, 0.59, "sim narrow #approx pure dof + curvature");
  }
  // [4] removed content per track
  cv->cd(4);
  {
    TH1D *h0 = new TH1D("rs4_0", "removed content per track: #sqrt{g^{2}#minus l^{2}};[#mum];tracks (unit area)", 60, 0, 6000);
    TH1D *h1 = (TH1D *) h0->Clone("rs4_1");
    for (double q : prem[0]) h0->Fill(std::min(q, 5999.));
    for (double q : prem[1]) h1->Fill(std::min(q, 5999.));
    for (auto h : {h0, h1}) h->Scale(1. / h->Integral());
    h0->SetLineColor(kBlack); h0->SetLineWidth(2);
    h1->SetLineColor(kBlue + 1); h1->SetLineWidth(2); h1->SetLineStyle(2);
    h0->SetMaximum(1.45 * std::max(h0->GetMaximum(), h1->GetMaximum()));
    h0->Draw("hist"); h1->Draw("hist same");
    TLegend *L = new TLegend(0.42, 0.70, 0.89, 0.87); L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextSize(0.030);
    L->AddEntry(h0, Form("real, med %.0f #mum", M(prem[0])), "l");
    L->AddEntry(h1, Form("sim, med %.0f #mum", M(prem[1])), "l");
    L->Draw();
    tx.SetTextSize(0.029);
    tx.DrawLatex(0.42, 0.64, "= per-track long-range + road content");
    tx.DrawLatex(0.42, 0.59, "(the quantity the one-number ratio hides)");
  }
  cv->SaveAs(Form("../sim_validation_plots/rscale_probe_%s.png", ver));
  printf("wrote ../sim_validation_plots/rscale_probe_%s.png\n", ver);
}
