// lsag_probe.C — what makes the REAL LOCAL 4-row short-sagitta RMS distribution
// (nf_digipix right panel) BIMODAL: shoulder ~500-800 um under the main peak
// ~1100 um. Windows rebuilt exactly as nf_digipix::doTrack (4 adjacent pad
// rows, 45 sliding windows, gate >=3 distinct rows & >=5 px, same MNF fitter,
// same road / ev44 veto on the real side). Splits tested per window:
//   n-class (5-6 / 7-9 / 10-14 / 15+ px) -> dof factor sqrt((n-3)/n)
//   nrows (3 vs 4) | pad region (R1 w<=19 / MIX 20-22 / R2 23-35 / MIX 36-38 /
//   R3 >=39) | side | sector | mean tbin (drift) | two-blob rows (rowPP>5mm)
// plus a dof-corrected overlay rms/sqrt((n-3)/n).
// usage: root -l -b -q 'lsag_probe.C+()'   out: lsag_probe_<ver>.txt +
// ../sim_validation_plots/lsag_probe_<ver>.png
#include "ms_nofinder.C"
#include <THStack.h>
#include <TH2D.h>
#include <TLine.h>
#include <TProfile.h>
#include <cstdint>

namespace LSP
{
struct GrpX
{
  std::vector<double> x, y, r;
  std::vector<float> tb;
  std::vector<uint8_t> side;
  int ev = -1;
};
struct Win
{
  float rms, rmsc, tb;   // um, dof-corrected um, mean tbin
  short n, nrows, w, sec, side, reg;  // reg 0=R1 1=MIX 2=R2 3=MIX 4=R3
  float rowPP;           // max within-row peak-to-peak residual, um
};
}  // namespace LSP

void lsag_probe(const char *digif = "digi_frames_production_v6.root", int nsim = 60,
                const char *realf = "../clusters_seeds_island_79507-0.root_ntuplizer.root",
                const char *ver = "v6")
{
  using namespace LSP;
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

  // ---- groups (verbatim road / truth grouping, as gtail_probe) -------------
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
      // printf("key values: %d, %d;   ", key.first, key.second);
      auto it = seedIdx.find(key);
      // printf("Value before insert: %d, %d, %d;   ", it->first.first, it->first.second, it->second);
      if (it == seedIdx.end()) it = seedIdx.insert({key, (int) seedIdx.size()}).first;
      // printf("Value after insert: %d, %d, %d   ", it->first.first, it->first.second, it->second);
      bucket[(int) ev * 100 + (int) lay].push_back((int) sc.size());
      sc.push_back({x, y, tb, it->second});
      // printf("Value in sc: %d\n", it->second);
      // printf("**********************************************\n");
    }
    f->Close();
    gr[0].assign(seedIdx.size(), GrpX());
    for (auto &kv : seedIdx) gr[0][kv.second].ev = kv.first.first;
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
        G.tb.push_back(tb); G.side.push_back((uint8_t) zel);
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
    float ev, lay, phi, adc, tid, tb, zel;
    t->SetBranchStatus("*", 0);
    for (auto b : {"event", "layer", "phi", "adc", "gtrackID", "tbin", "zelem"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("event", &ev); t->SetBranchAddress("layer", &lay);
    t->SetBranchAddress("phi", &phi); t->SetBranchAddress("adc", &adc);
    t->SetBranchAddress("gtrackID", &tid); t->SetBranchAddress("tbin", &tb);
    t->SetBranchAddress("zelem", &zel);
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
      G.tb.push_back(tb); G.side.push_back((uint8_t) zel);
    }
    f->Close();
    printf("sim: %zu truth groups\n", gr[1].size());
  }

  // ---- windows (verbatim nf_digipix::doTrack construction) -----------------
  std::vector<Win> wn[2];
  for (int s = 0; s < 2; ++s)
    for (auto &G : gr[s])
    {
      MNF::Grp B; B.x = G.x; B.y = G.y; B.r = G.r;
      MNF::Fit F0;
      if (!MNF::fitBar(B, F0)) continue;                 // same track population
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
        Win W;
        W.n = (short) widx[w].size(); W.nrows = (short) wrow[w].size(); W.w = (short) (w + 7);
        W.rms = (float) (L.rms * 1e4);
        W.rmsc = (float) (W.rms / std::sqrt((W.n - 3.) / W.n));
        double mtb = 0, mphi = 0; int nsd = 0;
        std::map<int, std::vector<double>> rows;
        for (int i : widx[w])
        {
          mtb += G.tb[i]; nsd += G.side[i];
          mphi += std::atan2(G.y[i], G.x[i]);
          rows[nearRow(G.r[i])].push_back(std::hypot(G.x[i] - L.a, G.y[i] - L.b) - L.R);
        }
        W.tb = (float) (mtb / W.n);
        W.sec = (short) (((int) std::floor((mphi / W.n + 2 * M_PI) / (M_PI / 6))) % 12);
        W.side = (short) (2 * nsd >= W.n ? 1 : 0);
        int wl = w + 7, wh = wl + 3;
        W.reg = wh <= 22 ? 0 : (wl <= 22 ? 1 : (wh <= 38 ? 2 : (wl <= 38 ? 3 : 4)));
        W.rowPP = 0;
        for (auto &kv : rows)
          if (kv.second.size() >= 2)
          {
            double lo = 1e9, hi = -1e9;
            for (double q : kv.second) { lo = std::min(lo, q); hi = std::max(hi, q); }
            W.rowPP = std::max(W.rowPP, (float) ((hi - lo) * 1e4));
          }
        wn[s].push_back(W);
      }
    }
  printf("windows: real %zu sim %zu\n", wn[0].size(), wn[1].size());

  // ---- ledger ---------------------------------------------------------------
  FILE *fo = fopen(Form("lsag_probe_%s.txt", ver), "w");
  auto P = [&](const char *fmt, ...) {
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    va_start(ap, fmt); vfprintf(fo, fmt, ap); va_end(ap);
  };
  auto medof = [&](std::vector<double> &v) { return MNF::med(v); };
  P("[lsag_probe %s] anatomy of the LOCAL 4-row sagitta RMS bimodality (nf_digipix windows)\n", ver);
  const char *sn2[2] = {"real", "sim "};
  // overall + shoulder/peak occupancy
  for (int s = 0; s < 2; ++s)
  {
    std::vector<double> a;
    long nsh = 0, npk = 0;
    for (auto &W : wn[s]) { a.push_back(W.rms); if (W.rms > 300 && W.rms < 900) nsh++; if (W.rms >= 900 && W.rms < 1600) npk++; }
    P("  %s: %zu windows | med %.0f um | in shoulder [300,900) %.1f%% | in peak [900,1600) %.1f%%\n",
      sn2[s], wn[s].size(), medof(a), 100. * nsh / wn[s].size(), 100. * npk / wn[s].size());
  }
  // split tables
  auto table = [&](const char *name, std::function<int(const Win &)> cls, int ncls, std::function<const char *(int)> lbl) {
    for (int s = 0; s < 2; ++s)
    {
      P("  -- %s (%s):\n", name, sn2[s]);
      for (int c = 0; c < ncls; ++c)
      {
        std::vector<double> v; long nsh = 0, ntot = 0;
        for (auto &W : wn[s]) if (cls(W) == c) { v.push_back(W.rms); ntot++; if (W.rms > 300 && W.rms < 900) nsh++; }
        if (!ntot) continue;
        P("     %-14s %7ld (%4.1f%%) med %4.0f um | shoulder frac %4.1f%%\n",
          lbl(c), ntot, 100. * ntot / wn[s].size(), medof(v), 100. * nsh / ntot);
      }
    }
  };
  table("n-class", [](const Win &W) { return W.n <= 6 ? 0 : W.n <= 9 ? 1 : W.n <= 14 ? 2 : 3; }, 4,
        [](int c) { return (const char *[]){"n=5-6", "n=7-9", "n=10-14", "n>=15"}[c]; });
  table("nrows", [](const Win &W) { return W.nrows - 3; }, 2,
        [](int c) { return (const char *[]){"3 rows", "4 rows"}[c]; });
  table("pad region", [](const Win &W) { return (int) W.reg; }, 5,
        [](int c) { return (const char *[]){"R1", "R1/R2 mix", "R2", "R2/R3 mix", "R3"}[c]; });
  table("side", [](const Win &W) { return (int) W.side; }, 2,
        [](int c) { return (const char *[]){"side 0", "side 1"}[c]; });
  table("drift (tbin)", [](const Win &W) { return W.tb < 120 ? 0 : W.tb < 235 ? 1 : W.tb < 350 ? 2 : 3; }, 4,
        [](int c) { return (const char *[]){"tb<120 long", "tb 120-235", "tb 235-350", "tb>350 mixed"}[c]; });
  table("two-blob", [](const Win &W) { return W.rowPP > 5000 ? 1 : 0; }, 2,
        [](int c) { return (const char *[]){"rowPP<=5mm", "rowPP>5mm"}[c]; });
  // sparsity source: px per row, n medians
  for (int s = 0; s < 2; ++s)
  {
    std::vector<double> ppr, nn;
    for (auto &W : wn[s]) { ppr.push_back((double) W.n / W.nrows); nn.push_back(W.n); }
    P("  %s window stats: med n %.0f | med px-per-row %.2f\n", sn2[s], medof(nn), medof(ppr));
  }
  // matched clean-response cell: rowPP<=5mm & n 10-14
  for (int s = 0; s < 2; ++s)
  {
    std::vector<double> v;
    for (auto &W : wn[s]) if (W.rowPP <= 5000 && W.n >= 10 && W.n <= 14) v.push_back(W.rms);
    P("  %s matched clean cell (rowPP<=5mm, n 10-14): %zu windows, med %.0f um\n", sn2[s], v.size(), medof(v));
  }
  // closure: reweight sim windows to the real n-mix
  double rwMed = 0, rwSh = 0;
  {
    std::map<int, double> cR, cS;
    for (auto &W : wn[0]) cR[std::min((int) W.n, 40)]++;
    for (auto &W : wn[1]) cS[std::min((int) W.n, 40)]++;
    double totR = wn[0].size(), totS = wn[1].size(), sw = 0, swsh = 0;
    TH1D hq("hq_rw", "", 400, 0, 4000);
    for (auto &W : wn[1])
    {
      int nb = std::min((int) W.n, 40);
      if (!cS.count(nb) || cS[nb] == 0 || !cR.count(nb)) continue;
      double w = (cR[nb] / totR) / (cS[nb] / totS);
      hq.Fill(std::min((double) W.rms, 3999.), w);
      sw += w; if (W.rms > 300 && W.rms < 900) swsh += w;
    }
    double q = 0.5, x = 0; hq.GetQuantiles(1, &x, &q);
    rwMed = x; rwSh = 100 * swsh / sw;
    P("  CLOSURE: sim reweighted to the real n-mix: med %.0f um | shoulder frac %.1f%% (real 21.9%%)\n", rwMed, rwSh);
  }
  // sector medians (real only, hunting hardware)
  {
    P("  -- sector x side median rms (real):\n");
    for (int sd = 0; sd < 2; ++sd)
    {
      P("     side %d:", sd);
      for (int sec = 0; sec < 12; ++sec)
      {
        std::vector<double> v;
        for (auto &W : wn[0]) if (W.side == sd && W.sec == sec) v.push_back(W.rms);
        P(" %4.0f", medof(v));
      }
      P("\n");
    }
  }
  fclose(fo);
  printf("wrote lsag_probe_%s.txt\n", ver);

  // ---- figure ----------------------------------------------------------------
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.045);
  TCanvas *cv = new TCanvas("cvls", "lsag", 1700, 1250);
  cv->Divide(2, 2);
  // [1] real stacked by n-class
  cv->cd(1);
  {
    THStack *st = new THStack("ls_st1", "real LOCAL RMS by window pixel count;per-fit circle RMS [#mum];windows");
    int cc[4] = {kAzure - 9, kAzure + 2, kOrange - 3, kRed - 7};
    const char *ln[4] = {"n = 5-6 px", "n = 7-9", "n = 10-14", "n #geq 15"};
    TH1D *hh[4];
    for (int c = 0; c < 4; ++c)
    {
      hh[c] = new TH1D(Form("ls1_%d", c), "", 60, 0, 4000);
      hh[c]->SetFillColor(cc[c]); hh[c]->SetLineColor(cc[c]);
    }
    for (auto &W : wn[0]) hh[W.n <= 6 ? 0 : W.n <= 9 ? 1 : W.n <= 14 ? 2 : 3]->Fill(std::min((double) W.rms, 3999.));
    for (int c = 0; c < 4; ++c) st->Add(hh[c]);
    st->Draw("hist");
    // closure overlay: sim reweighted to the real n-mix, scaled to the real window count
    TH1D *hrw = new TH1D("ls1_rw", "", 60, 0, 4000);
    {
      std::map<int, double> cR, cS;
      for (auto &W : wn[0]) cR[std::min((int) W.n, 40)]++;
      for (auto &W : wn[1]) cS[std::min((int) W.n, 40)]++;
      double totR = wn[0].size(), totS = wn[1].size();
      for (auto &W : wn[1])
      {
        int nb = std::min((int) W.n, 40);
        if (!cS.count(nb) || cS[nb] == 0 || !cR.count(nb)) continue;
        hrw->Fill(std::min((double) W.rms, 3999.), (cR[nb] / totR) / (cS[nb] / totS));
      }
      hrw->Scale(wn[0].size() / hrw->Integral());
      hrw->SetLineColor(kGray + 3); hrw->SetLineWidth(3); hrw->SetLineStyle(2);
      hrw->Draw("hist same");
    }
    TLegend *L = new TLegend(0.52, 0.60, 0.89, 0.87); L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextSize(0.030);
    for (int c = 0; c < 4; ++c)
    {
      std::vector<double> v; for (auto &W : wn[0]) if ((W.n <= 6 ? 0 : W.n <= 9 ? 1 : W.n <= 14 ? 2 : 3) == c) v.push_back(W.rms);
      L->AddEntry(hh[c], Form("%s, med %.0f #mum", ln[c], MNF::med(v)), "f");
    }
    L->AddEntry((TObject *) gROOT->FindObject("ls1_rw"), "sim reweighted to the real n-mix", "l");
    L->Draw();
  }
  // [2] dof-corrected overlay
  cv->cd(2);
  {
    TH1D *h0 = new TH1D("ls2_0", "dof correction test: rms / #sqrt{(n-3)/n};per-fit circle RMS [#mum];windows (unit area)", 60, 0, 4000);
    TH1D *h1 = new TH1D("ls2_1", "", 60, 0, 4000);
    TH1D *h2 = new TH1D("ls2_2", "", 60, 0, 4000);
    for (auto &W : wn[0]) { h0->Fill(std::min((double) W.rms, 3999.)); h1->Fill(std::min((double) W.rmsc, 3999.)); }
    for (auto &W : wn[1]) h2->Fill(std::min((double) W.rmsc, 3999.));
    for (auto h : {h0, h1, h2}) h->Scale(1. / h->Integral());
    h0->SetLineColor(kBlack); h0->SetLineWidth(2);
    h1->SetLineColor(kRed + 1); h1->SetLineWidth(2); h1->SetLineStyle(7);
    h2->SetLineColor(kBlue + 1); h2->SetLineWidth(2); h2->SetLineStyle(3);
    h0->SetMaximum(1.4 * std::max({h0->GetMaximum(), h1->GetMaximum(), h2->GetMaximum()}));
    h0->Draw("hist"); h1->Draw("hist same"); h2->Draw("hist same");
    std::vector<double> v0, v1, v2;
    for (auto &W : wn[0]) { v0.push_back(W.rms); v1.push_back(W.rmsc); }
    for (auto &W : wn[1]) v2.push_back(W.rmsc);
    TLegend *L = new TLegend(0.42, 0.66, 0.89, 0.87); L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextSize(0.030);
    L->AddEntry(h0, Form("real raw, med %.0f #mum", MNF::med(v0)), "l");
    L->AddEntry(h1, Form("real dof-corrected, med %.0f #mum", MNF::med(v1)), "l");
    L->AddEntry(h2, Form("sim dof-corrected, med %.0f #mum", MNF::med(v2)), "l");
    L->Draw();
  }
  // [3] median rms vs window start row, real & sim
  cv->cd(3);
  {
    TH1D *p0 = new TH1D("ls3_0", "median window RMS vs start row;window start pad row;median RMS [#mum]", 45, 6.5, 51.5);
    TH1D *p1 = (TH1D *) p0->Clone("ls3_1");
    for (int w = 7; w <= 51; ++w)
      for (int s = 0; s < 2; ++s)
      {
        std::vector<double> v;
        for (auto &W : wn[s]) if (W.w == w) v.push_back(W.rms);
        (s ? p1 : p0)->SetBinContent(w - 6, MNF::med(v));
      }
    p0->SetLineColor(kBlack); p0->SetLineWidth(2);
    p1->SetLineColor(kBlue + 1); p1->SetLineWidth(2); p1->SetLineStyle(2);
    p0->SetMinimum(0); p0->SetMaximum(2200);
    p0->Draw("hist"); p1->Draw("hist same");
    TLine bl; bl.SetLineStyle(2); bl.SetLineColor(kGray + 2);
    bl.DrawLine(22.5 - 3, 0, 22.5 - 3, 2200); bl.DrawLine(38.5 - 3, 0, 38.5 - 3, 2200);
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.030);
    tx.DrawLatex(0.30, 0.20, "dashed: last all-R1 / all-R2 window starts");
    TLegend *L = new TLegend(0.55, 0.74, 0.89, 0.87); L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextSize(0.030);
    L->AddEntry(p0, "real", "l"); L->AddEntry(p1, "sim digi", "l");
    L->Draw();
  }
  // [4] side x sector median map (real)
  cv->cd(4);
  {
    gPad->SetRightMargin(0.14);
    TH2D *hm = new TH2D("ls4", "real median window RMS by side #times sector;sector;side", 12, -0.5, 11.5, 2, -0.5, 1.5);
    for (int sd = 0; sd < 2; ++sd)
      for (int sec = 0; sec < 12; ++sec)
      {
        std::vector<double> v;
        for (auto &W : wn[0]) if (W.side == sd && W.sec == sec) v.push_back(W.rms);
        hm->SetBinContent(sec + 1, sd + 1, MNF::med(v));
      }
    hm->GetZaxis()->SetTitle("median RMS [#mum]");
    gStyle->SetPaintTextFormat(".0f");
    hm->SetMarkerSize(1.5);
    hm->Draw("colz text");
  }
  cv->SaveAs(Form("../sim_validation_plots/lsag_probe_%s.png", ver));
  printf("wrote ../sim_validation_plots/lsag_probe_%s.png\n", ver);
}
