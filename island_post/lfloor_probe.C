// lfloor_probe.C — the "surprise of only halved": why does the LOCAL 4-row
// sagitta fit reduce the real GLOBAL pixel-fit median only from 2358 to
// 1205 um (nf_digipix v6)?
// Claim tested: RMS combines in quadrature, so LOCAL already removed
// sqrt(2358^2-1205^2) = 2027 um = 74% of the variance (road pickup + field +
// curvature); the 1205 um floor is the CHARGE-CLOUD WIDTH itself — the fit
// runs on raw pixels (pad-center positions repeated per tbin sample), so no
// trajectory fit, however local, can go below the cloud's r-phi spread. The
// way below the floor is CENTROIDING (= what clustering does): the same
// pixels, ADC-centroided per row and refit, land at the cluster-level scale
// (nf_tracks: 780 um).
// Measurements per track (real + sim digi, same groups/fitter as nf_digipix,
// global 3xMAD trim from gtail_probe to keep road junk out of the floor):
//   - window decomposition: rms vs HF (within-row) vs LF (row-mean) — the
//     4-row fit zeroes row-to-row content by construction (3 dof on 3-4 rows)
//   - fit-free cloud width: per-row r-phi sample RMS, in um and pad-pitch
//     units; distinct pads/row; tbin samples/pad (the comb duplication)
//   - ADC-weighted row centroids -> whole-track circle refit = RMS_centroid
// usage: root -l -b -q 'lfloor_probe.C+()'
// out: lfloor_probe_<ver>.txt + ../sim_validation_plots/lfloor_probe_<ver>.png
#include "../sim_validation_plots/src/ms_nofinder.C"
#include <THStack.h>
#include <TH2D.h>
#include <TLine.h>
#include <cstdint>

namespace LFP
{
struct GrpX
{
  std::vector<double> x, y, r;
  std::vector<float> adc;
  int ev = -1;
};
bool trimFit(const GrpX &G, MNF::Fit &F, std::vector<int> &kept,
             int iters = 3, double floorcm = 0.05)
{
  kept.resize(G.x.size());
  for (size_t k = 0; k < kept.size(); ++k) kept[k] = (int) k;
  for (int it = 0; it < iters; ++it)
  {
    std::vector<double> X, Y;
    for (int i : kept) { X.push_back(G.x[i]); Y.push_back(G.y[i]); }
    F = MNF::fitCircle(X, Y);
    if (!F.ok) return false;
    std::vector<double> res;
    for (int i : kept) res.push_back(std::hypot(G.x[i] - F.a, G.y[i] - F.b) - F.R);
    std::vector<double> tmp = res;
    double md = MNF::med(tmp);
    for (double &q : tmp) q = std::fabs(q - md);
    double thr = std::max(3 * 1.4826 * MNF::med(tmp), floorcm);
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
}  // namespace LFP

void lfloor_probe(const char *digif = "digi_frames_production_v6.root", int nsim = 60,
                  const char *realf = "../clusters_seeds_island_79507-0.root_ntuplizer.root",
                  const char *ver = "v6")
{
  using namespace LFP;
  gROOT->SetBatch(1);
  double rowR[55], pitch[55];
  {
    FILE *fp = fopen("tpc_geom_table.txt", "r");
    if (!fp) { printf("no tpc_geom_table.txt\n"); return; }
    char line[512];
    while (fgets(line, sizeof line, fp))
    {
      int L, nb; double r, sl, p0, p1;
      if (sscanf(line, "%d %d %lf %lf %lf %lf", &L, &nb, &r, &sl, &p0, &p1) == 6 && L >= 7 && L <= 54)
      { rowR[L] = r; pitch[L] = 2 * M_PI * r / nb; }
    }
    fclose(fp);
  }
  auto nearRow = [&](double r) -> int {
    int best = -1; double bd = 1e9;
    for (int L = 7; L <= 54; ++L)
    { double d = std::fabs(r - rowR[L]); if (d < bd) { bd = d; best = L; } }
    return bd < 0.60 ? best : -1;
  };

  // ---- groups (same construction as nf_digipix / gtail_probe, + adc) -------
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
        G.adc.push_back(adc);
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
      G.adc.push_back(adc);
    }
    f->Close();
    printf("sim: %zu truth groups\n", gr[1].size());
  }

  // ---- per-track measurements ----------------------------------------------
  // collections (um unless noted)
  std::vector<double> vGraw[2], vGtrim[2], vLoc[2], vHF[2], vLFw[2], vCent[2];
  std::vector<double> vWidth[2], vWpitch[2], vPads[2], vSpp[2], vWadc[2];
  long ntrk[2] = {0, 0};
  for (int s = 0; s < 2; ++s)
    for (auto &G : gr[s])
    {
      MNF::Grp B; B.x = G.x; B.y = G.y; B.r = G.r;
      MNF::Fit F0;
      if (!MNF::fitBar(B, F0)) continue;
      ntrk[s]++;
      vGraw[s].push_back(F0.rms * 1e4);
      MNF::Fit FT; std::vector<int> kept;
      if (!trimFit(G, FT, kept)) continue;
      vGtrim[s].push_back(FT.rms * 1e4);
      // ---- per-row structures on the trimmed set
      std::map<int, std::vector<int>> rows;
      for (int i : kept) { int L = nearRow(G.r[i]); if (L >= 0) rows[L].push_back(i); }
      std::vector<double> cx, cy;                      // ADC-weighted centroids
      for (auto &kv : rows)
      {
        int L = kv.first;
        auto &v = kv.second;
        double sw = 0, sx = 0, sy = 0;
        for (int i : v) { sw += G.adc[i]; sx += G.adc[i] * G.x[i]; sy += G.adc[i] * G.y[i]; }
        if (sw > 0) { cx.push_back(sx / sw); cy.push_back(sy / sw); }
        if ((int) v.size() < 2) continue;
        // fit-free r-phi cloud width of this row (around the unweighted mean angle)
        double mphi = 0;
        for (int i : v) mphi += std::atan2(G.y[i], G.x[i]);
        mphi /= v.size();
        double s1 = 0, s2 = 0, w1 = 0, w2 = 0;
        std::set<long> pads;
        for (int i : v)
        {
          double d = MNF::wrapphi(std::atan2(G.y[i], G.x[i]) - mphi) * rowR[L];
          s1 += d; s2 += d * d;
          w1 += G.adc[i] * d; w2 += G.adc[i] * d * d;
          pads.insert(std::lround(std::atan2(G.y[i], G.x[i]) / (pitch[L] / rowR[L])));
        }
        double m = s1 / v.size(), wid = std::sqrt(std::max(0., s2 / v.size() - m * m));
        double sw2 = 0; for (int i : v) sw2 += G.adc[i];
        double wm = w1 / sw2, wadc = std::sqrt(std::max(0., w2 / sw2 - wm * wm));
        vWidth[s].push_back(wid * 1e4);
        vWadc[s].push_back(wadc * 1e4);
        vWpitch[s].push_back(wid / pitch[L]);
        vPads[s].push_back((double) pads.size());
        vSpp[s].push_back((double) v.size() / pads.size());
      }
      // ---- centroid-level whole-track refit
      if ((int) cx.size() >= 12)
      {
        double rlo = 1e9, rhi = 0;
        for (size_t k = 0; k < cx.size(); ++k)
        { double rr = std::hypot(cx[k], cy[k]); rlo = std::min(rlo, rr); rhi = std::max(rhi, rr); }
        if (rhi - rlo >= 15)
        {
          MNF::Fit FC = MNF::fitCircle(cx, cy);
          if (FC.ok && FC.R >= 45 && FC.R < 2e4) vCent[s].push_back(FC.rms * 1e4);
        }
      }
      // ---- local windows (verbatim nf_digipix) + HF/LF decomposition
      std::vector<std::vector<int>> widx(45);
      std::vector<std::set<int>> wrow(45);
      for (size_t i = 0; i < G.x.size(); ++i)
      {
        int row = nearRow(G.r[i]);
        if (row < 0) continue;
        for (int w = std::max(7, row - 3); w <= std::min(51, row); ++w)
        { widx[w - 7].push_back((int) i); wrow[w - 7].insert(row); }
      }
      for (int w = 0; w < 45; ++w)
      {
        if ((int) wrow[w].size() < 3 || (int) widx[w].size() < 5) continue;
        std::vector<double> X, Y;
        for (int i : widx[w]) { X.push_back(G.x[i]); Y.push_back(G.y[i]); }
        MNF::Fit L = MNF::fitCircle(X, Y);
        if (!L.ok) continue;
        vLoc[s].push_back(L.rms * 1e4);
        std::map<int, std::vector<double>> rr;
        for (int i : widx[w])
        { int Lr = nearRow(G.r[i]); if (Lr >= 0) rr[Lr].push_back(std::hypot(G.x[i] - L.a, G.y[i] - L.b) - L.R); }
        double lf2 = 0, hf2 = 0; long nT = 0;
        for (auto &kv : rr)
        {
          double m = 0; for (double q : kv.second) m += q; m /= kv.second.size();
          lf2 += m * m * kv.second.size();
          for (double q : kv.second) hf2 += (q - m) * (q - m);
          nT += (long) kv.second.size();
        }
        if (nT) { vHF[s].push_back(std::sqrt(hf2 / nT) * 1e4); vLFw[s].push_back(std::sqrt(lf2 / nT) * 1e4); }
      }
    }

  // ---- ledger ----------------------------------------------------------------
  FILE *fo = fopen(Form("lfloor_probe_%s.txt", ver), "w");
  auto P = [&](const char *fmt, ...) {
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    va_start(ap, fmt); vfprintf(fo, fmt, ap); va_end(ap);
  };
  auto M = [&](std::vector<double> &v) { return MNF::med(v); };
  P("[lfloor_probe %s] the 'only halved' question: what the LOCAL fit removed and what its floor is\n", ver);
  for (int s = 0; s < 2; ++s)
  {
    P("  %s (%ld tracks):\n", s ? "sim " : "real", ntrk[s]);
    P("    LADDER: pixel GLOBAL raw %4.0f -> 3sig-clipped %4.0f -> pixel LOCAL 4-row %4.0f -> row-centroid GLOBAL refit %4.0f um (%d tracks)\n",
      M(vGraw[s]), M(vGtrim[s]), M(vLoc[s]), M(vCent[s]), (int) vCent[s].size());
    P("    variance removed by LOCAL vs raw GLOBAL: sqrt(%0.f^2-%0.f^2)=%4.0f um = %.0f%% of the variance\n",
      M(vGraw[s]), M(vLoc[s]), std::sqrt(std::max(0., M(vGraw[s]) * M(vGraw[s]) - M(vLoc[s]) * M(vLoc[s]))),
      100. * (1. - (M(vLoc[s]) * M(vLoc[s])) / (M(vGraw[s]) * M(vGraw[s]))));
    P("    window decomposition: rms med %4.0f = HF (within-row) %4.0f (+) LF (row means) %4.0f um -> the 4-row fit zeroes trajectory content\n",
      M(vLoc[s]), M(vHF[s]), M(vLFw[s]));
    P("    fit-free cloud width per row: sample RMS med %4.0f um (= %.2f pad pitch) | ADC-weighted %4.0f um | pads/row med %.0f | tbin samples/pad med %.1f\n",
      M(vWidth[s]), M(vWpitch[s]), M(vWadc[s]), M(vPads[s]), M(vSpp[s]));
  }
  P("  reference: nf_tracks cluster-level whole-track fit (tracker clusters) = 780 um -> centroid refit closes onto it\n");
  fclose(fo);
  printf("wrote lfloor_probe_%s.txt\n", ver);

  // ---- figure -----------------------------------------------------------------
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.045);
  TCanvas *cv = new TCanvas("cvlf", "lfloor", 1700, 1250);
  cv->Divide(2, 2);
  TLatex tx; tx.SetNDC();
  // [1] the ladder (real)
  cv->cd(1);
  {
    std::vector<double> *vv[4] = {&vGraw[0], &vGtrim[0], &vLoc[0], &vCent[0]};
    const char *nm[4] = {"pixel GLOBAL raw", "pixel GLOBAL 3#sigma-clipped", "pixel LOCAL 4-row", "row-centroid GLOBAL refit"};
    int cc[4] = {kBlack, kGray + 2, kBlue + 1, kGreen + 2};
    int ls[4] = {1, 7, 1, 1};
    TH1D *h[4];
    for (int k = 0; k < 4; ++k)
    {
      h[k] = new TH1D(Form("lf1_%d", k), "REAL: the estimator ladder;per-fit circle RMS [#mum];fits (unit area)", 75, 0, 6000);
      for (double q : *vv[k]) h[k]->Fill(std::min(q, 5999.));
      h[k]->Scale(1. / h[k]->Integral());
      h[k]->SetLineColor(cc[k]); h[k]->SetLineWidth(2); h[k]->SetLineStyle(ls[k]);
    }
    double mx = 0; for (int k = 0; k < 4; ++k) mx = std::max(mx, h[k]->GetMaximum());
    h[0]->SetMaximum(1.45 * mx);
    h[0]->Draw("hist"); for (int k = 1; k < 4; ++k) h[k]->Draw("hist same");
    TLegend *L = new TLegend(0.37, 0.60, 0.89, 0.87); L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextSize(0.029);
    for (int k = 0; k < 4; ++k) L->AddEntry(h[k], Form("%s, med %.0f #mum", nm[k], M(*vv[k])), "l");
    L->Draw();
    tx.SetTextSize(0.030);
    tx.DrawLatex(0.37, 0.54, "no trajectory fit beats the cloud width;");
    tx.DrawLatex(0.37, 0.49, "centroiding (= clustering) does");
  }
  // [2] window decomposition (real)
  cv->cd(2);
  {
    TH1D *h0 = new TH1D("lf2_0", "REAL windows: rms vs its two components;[#mum];windows (unit area)", 60, 0, 3000);
    TH1D *h1 = (TH1D *) h0->Clone("lf2_1"); TH1D *h2 = (TH1D *) h0->Clone("lf2_2");
    for (double q : vLoc[0]) h0->Fill(std::min(q, 2999.));
    for (double q : vHF[0]) h1->Fill(std::min(q, 2999.));
    for (double q : vLFw[0]) h2->Fill(std::min(q, 2999.));
    for (auto h : {h0, h1, h2}) h->Scale(1. / h->Integral());
    h0->SetLineColor(kBlack); h0->SetLineWidth(3);
    h1->SetLineColor(kOrange + 7); h1->SetLineWidth(2); h1->SetLineStyle(7);
    h2->SetLineColor(kAzure + 2); h2->SetLineWidth(2);
    h0->SetMaximum(1.45 * std::max({h0->GetMaximum(), h1->GetMaximum(), h2->GetMaximum()}));
    h0->Draw("hist"); h1->Draw("hist same"); h2->Draw("hist same");
    TLegend *L = new TLegend(0.40, 0.64, 0.89, 0.87); L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextSize(0.029);
    L->AddEntry(h0, Form("window RMS, med %.0f #mum", M(vLoc[0])), "l");
    L->AddEntry(h1, Form("HF: within-row scatter, med %.0f #mum", M(vHF[0])), "l");
    L->AddEntry(h2, Form("LF: row means after fit, med %.0f #mum", M(vLFw[0])), "l");
    L->Draw();
    tx.SetTextSize(0.030);
    tx.DrawLatex(0.40, 0.58, "RMS #approx HF alone: the sagitta content");
    tx.DrawLatex(0.40, 0.53, "is already zeroed - the floor is the cloud");
  }
  // [3] fit-free cloud width real vs sim
  cv->cd(3);
  {
    TH1D *h0 = new TH1D("lf3_0", "fit-free r#phi cloud width per pad row;within-row sample RMS [#mum];rows (unit area)", 60, 0, 4000);
    TH1D *h1 = (TH1D *) h0->Clone("lf3_1");
    for (double q : vWidth[0]) h0->Fill(std::min(q, 3999.));
    for (double q : vWidth[1]) h1->Fill(std::min(q, 3999.));
    for (auto h : {h0, h1}) h->Scale(1. / h->Integral());
    h0->SetLineColor(kBlack); h0->SetLineWidth(2);
    h1->SetLineColor(kBlue + 1); h1->SetLineWidth(2); h1->SetLineStyle(2);
    h0->SetMaximum(1.45 * std::max(h0->GetMaximum(), h1->GetMaximum()));
    h0->Draw("hist"); h1->Draw("hist same");
    TLegend *L = new TLegend(0.40, 0.68, 0.89, 0.87); L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextSize(0.029);
    L->AddEntry(h0, Form("real, med %.0f #mum = %.2f pad pitch", M(vWidth[0]), M(vWpitch[0])), "l");
    L->AddEntry(h1, Form("sim digi, med %.0f #mum = %.2f pad pitch", M(vWidth[1]), M(vWpitch[1])), "l");
    L->Draw();
    tx.SetTextSize(0.029);
    tx.DrawLatex(0.40, 0.62, Form("real: %.0f pads/row #times %.1f tbin samples/pad", M(vPads[0]), M(vSpp[0])));
    tx.DrawLatex(0.40, 0.57, Form("ADC-weighted width: real %.0f / sim %.0f #mum", M(vWadc[0]), M(vWadc[1])));
    tx.DrawLatex(0.40, 0.52, "pixels = pad centers repeated per tbin:");
    tx.DrawLatex(0.40, 0.47, "the fit sees the pad comb, not a point");
  }
  // [4] centroid refit real vs sim + the 780 reference
  cv->cd(4);
  {
    TH1D *h0 = new TH1D("lf4_0", "row-centroid whole-track refit (same pixels, centroided);per-fit circle RMS [#mum];tracks (unit area)", 60, 0, 3000);
    TH1D *h1 = (TH1D *) h0->Clone("lf4_1");
    for (double q : vCent[0]) h0->Fill(std::min(q, 2999.));
    for (double q : vCent[1]) h1->Fill(std::min(q, 2999.));
    for (auto h : {h0, h1}) h->Scale(1. / h->Integral());
    h0->SetLineColor(kGreen + 2); h0->SetLineWidth(2);
    h1->SetLineColor(kBlue + 1); h1->SetLineWidth(2); h1->SetLineStyle(2);
    h0->SetMaximum(1.45 * std::max(h0->GetMaximum(), h1->GetMaximum()));
    h0->Draw("hist"); h1->Draw("hist same");
    TLine ln; ln.SetLineStyle(2); ln.SetLineColor(kGray + 2);
    ln.DrawLine(780, 0, 780, h0->GetMaximum() / 1.45);
    TLegend *L = new TLegend(0.40, 0.68, 0.89, 0.87); L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextSize(0.029);
    L->AddEntry(h0, Form("real, med %.0f #mum", M(vCent[0])), "l");
    L->AddEntry(h1, Form("sim digi, med %.0f #mum", M(vCent[1])), "l");
    L->Draw();
    tx.SetTextSize(0.029);
    tx.DrawLatex(0.40, 0.62, "dashed: nf_tracks cluster-level fit, 780 #mum");
    tx.DrawLatex(0.40, 0.57, "centroiding removes the cloud width;");
    tx.DrawLatex(0.40, 0.52, "what remains is resolution + field");
  }
  cv->SaveAs(Form("../sim_validation_plots/lfloor_probe_%s.png", ver));
  printf("wrote ../sim_validation_plots/lfloor_probe_%s.png\n", ver);
}
