// roadscan_probe.C — scan the pixel-association ROAD (nominal dxy < 1.2 cm,
// |dtbin| <= 6; ms_nofinder.C::realPixGroups) and watch every meter move
// (user, 2026-08-20: "scan with the road see what will happen").
// Method: one pass over ntp_hit collecting, per pixel, ALL candidate seed
// clusters within the WIDEST road (dxy < 1.8 cm, |dtbin| <= 10) in the
// original bucket order; each (RXY, DT) config then replays the exact
// nominal assignment offline (first in-road cluster wins; hitID-dedup with
// break semantics) — bit-exact with realPixGroups at (1.2, 6).
// Per config: GLOBAL raw / 3xMAD-clipped medians, LOCAL 4-row median,
// pickup = sqrt(Graw^2 - Gclip^2), px/track, clip cost; data/MC against the
// truth-grouped sim digi reference (no road on sim), recomputed here.
// usage: root -l -b -q 'roadscan_probe.C+()'
// out: roadscan_probe_<ver>.txt + ../sim_validation_plots/roadscan_probe_<ver>.png
#include "../sim_validation_plots/src/ms_nofinder.C"
#include <TROOT.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TLine.h>
#include <unordered_map>
#include <cstdint>

namespace RDS
{
struct Cand { int trk; float dxy, dtb; };
struct Pix { float x, y; int hitID; std::vector<Cand> cand; };
double medv(std::vector<double> v)
{
  if (v.empty()) return 0;
  std::sort(v.begin(), v.end());
  return v[v.size() / 2];
}
double trimmedRMS(const std::vector<double> &X, const std::vector<double> &Y, double &fdrop)
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
  fdrop = 1. - (double) kept.size() / X.size();
  std::vector<double> xs, ys;
  for (int i : kept) { xs.push_back(X[i]); ys.push_back(Y[i]); }
  F = MNF::fitCircle(xs, ys);
  return F.ok ? F.rms : -1;
}
}  // namespace RDS

void roadscan_probe(const char *digif = "digi_frames_production_v6.root", int nsim = 60,
                    const char *realf = "../clusters_seeds_island_79507-0.root_ntuplizer.root",
                    const char *ver = "v6")
{
  using namespace RDS;
  gROOT->SetBatch(1);
  const double RXYMAX = 1.8, DTMAX = 10;
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

  // ---- pass 1: seed clusters ------------------------------------------------
  struct SC { float x, y, tb; int seed; };
  std::vector<SC> sc;
  std::unordered_map<int, std::vector<int>> bucket;
  int nseed = 0;
  {
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
    nseed = (int) seedIdx.size();
    printf("real: %d seeds\n", nseed);
  }
  // ---- pass 2: pixels with candidate lists (widest road) --------------------
  std::vector<Pix> px;
  {
    TFile *f = TFile::Open(realf);
    TTree *t = (TTree *) f->Get("ntp_hit");
    float ev, lay, x, y, tb, adc, hitID;
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
      Pix P; P.x = x; P.y = y; P.hitID = (int) hitID;
      for (int j : bit->second)
      {
        const SC &c = sc[j];
        float dtb = std::fabs(tb - c.tb);
        if (dtb > DTMAX) continue;
        double dx = x - c.x, dy = y - c.y, dxy = std::sqrt(dx * dx + dy * dy);
        if (dxy > RXYMAX) continue;
        P.cand.push_back({c.seed, (float) dxy, dtb});
      }
      if (!P.cand.empty()) px.push_back(std::move(P));
    }
    f->Close();
    printf("real: %zu pixels with candidates (widest road)\n", px.size());
  }
  // ---- sim truth-grouped reference (no road) --------------------------------
  double simG = 0, simL = 0;
  {
    TFile *f = TFile::Open(digif);
    TTree *t = (TTree *) f->Get("ntp_hit");
    float ev, lay, phi, adc, tid;
    t->SetBranchStatus("*", 0);
    for (auto b : {"event", "layer", "phi", "adc", "gtrackID"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("event", &ev); t->SetBranchAddress("layer", &lay);
    t->SetBranchAddress("phi", &phi); t->SetBranchAddress("adc", &adc);
    t->SetBranchAddress("gtrackID", &tid);
    std::map<std::pair<int, int>, MNF::Grp> g;
    for (Long64_t i = 0; i < t->GetEntries(); ++i)
    {
      t->GetEntry(i);
      if ((int) ev >= nsim) continue;
      if (lay < 7 || lay > 54 || adc <= 0 || tid <= 0) continue;
      double r = rowR[(int) lay];
      MNF::Grp &G = g[{(int) ev, (int) tid}];
      G.x.push_back(r * std::cos(phi)); G.y.push_back(r * std::sin(phi)); G.r.push_back(r);
    }
    f->Close();
    std::vector<double> vg, vl;
    for (auto &kv : g)
    {
      MNF::Fit F0;
      if (!MNF::fitBar(kv.second, F0)) continue;
      vg.push_back(F0.rms * 1e4);
      auto &G = kv.second;
      std::map<int, std::vector<int>> rows;
      for (size_t i = 0; i < G.x.size(); ++i)
      { int L = nearRow(G.r[i]); if (L >= 0) rows[L].push_back((int) i); }
      for (int w = 7; w <= 51; ++w)
      {
        std::vector<double> X, Y; int nr = 0;
        for (int rr = w; rr < w + 4; ++rr)
        { auto it = rows.find(rr); if (it == rows.end()) continue; nr++;
          for (int i : it->second) { X.push_back(G.x[i]); Y.push_back(G.y[i]); } }
        if (nr < 3 || (int) X.size() < 5) continue;
        MNF::Fit F = MNF::fitCircle(X, Y);
        if (F.ok) vl.push_back(F.rms * 1e4);
      }
    }
    simG = medv(vg); simL = medv(vl);
    printf("sim reference (truth-grouped): GLOBAL med %.0f | LOCAL med %.0f\n", simG, simL);
  }

  // ---- the scan --------------------------------------------------------------
  const int NR = 5, NDT = 3;
  double RXY[NR] = {0.3, 0.6, 0.9, 1.2, 1.8};
  double DTW[NDT] = {3, 6, 10};
  FILE *fo = fopen(Form("roadscan_probe_%s.txt", ver), "w");
  auto P = [&](const char *fmt, ...) {
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    va_start(ap, fmt); vfprintf(fo, fmt, ap); va_end(ap);
  };
  P("[roadscan_probe %s] road scan of the real pixel association (nominal 1.2 cm, +-6 tbins)\n", ver);
  P("  sim reference (truth-grouped digi, no road): GLOBAL %.0f | LOCAL %.0f um\n", simG, simL);
  P("  RXY  DT | tracks px/trk | Graw  Gclip  pickup | LOCAL | dM/C: Graw Gclip LOCAL | clip%%\n");
  double gr_[NDT][NR], gc_[NDT][NR], lo_[NDT][NR], pk_[NDT][NR], pt_[NDT][NR];
  for (int idt = 0; idt < NDT; ++idt)
    for (int ir = 0; ir < NR; ++ir)
    {
      double rxy = RXY[ir], dtw = DTW[idt];
      // assignment (exact nominal semantics)
      std::vector<MNF::Grp> grp(nseed);
      std::vector<std::set<int>> hid(nseed);
      for (auto &Pp : px)
        for (const Cand &c : Pp.cand)
        {
          if (c.dxy > rxy || c.dtb > dtw) continue;
          if (!hid[c.trk].insert(Pp.hitID).second) break;
          grp[c.trk].x.push_back(Pp.x);
          grp[c.trk].y.push_back(Pp.y);
          grp[c.trk].r.push_back(std::hypot(Pp.x, Pp.y));
          break;
        }
      // metrics
      std::vector<double> vg, vgc, vl, vfd;
      long npx2 = 0, ntr = 0;
      for (auto &G : grp)
      {
        MNF::Fit F0;
        if (!MNF::fitBar(G, F0)) continue;
        ntr++; npx2 += (long) G.x.size();
        vg.push_back(F0.rms * 1e4);
        double fd = 0, rt = trimmedRMS(G.x, G.y, fd);
        if (rt >= 0) { vgc.push_back(rt * 1e4); vfd.push_back(fd); }
        std::map<int, std::vector<int>> rows;
        for (size_t i = 0; i < G.x.size(); ++i)
        { int L = nearRow(G.r[i]); if (L >= 0) rows[L].push_back((int) i); }
        for (int w = 7; w <= 51; ++w)
        {
          std::vector<double> X, Y; int nr = 0;
          for (int rr = w; rr < w + 4; ++rr)
          { auto it = rows.find(rr); if (it == rows.end()) continue; nr++;
            for (int i : it->second) { X.push_back(G.x[i]); Y.push_back(G.y[i]); } }
          if (nr < 3 || (int) X.size() < 5) continue;
          MNF::Fit F = MNF::fitCircle(X, Y);
          if (F.ok) vl.push_back(F.rms * 1e4);
        }
      }
      double g = medv(vg), gc = medv(vgc), l = medv(vl);
      double pick = std::sqrt(std::max(0., g * g - gc * gc));
      gr_[idt][ir] = g; gc_[idt][ir] = gc; lo_[idt][ir] = l; pk_[idt][ir] = pick;
      pt_[idt][ir] = ntr ? (double) npx2 / ntr : 0;
      P("  %.1f %3.0f | %5ld  %5.0f | %5.0f %5.0f  %5.0f | %5.0f | %10.2f %5.2f %5.2f | %4.1f\n",
        rxy, dtw, ntr, pt_[idt][ir], g, gc, pick, l,
        simG > 0 ? g / simG : 0, simG > 0 ? gc / simG : 0, simL > 0 ? l / simL : 0,
        100 * medv(vfd));
    }
  fclose(fo);
  printf("wrote roadscan_probe_%s.txt\n", ver);

  // ---- figure ----------------------------------------------------------------
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.045);
  TCanvas *cv = new TCanvas("cvrd", "roadscan", 1700, 1250);
  cv->Divide(2, 2);
  TLatex tx; tx.SetNDC();
  int cols[NDT] = {kAzure + 2, kBlack, kOrange + 7};
  int mks[NDT] = {24, 20, 26};
  const char *dl[NDT] = {"|dtbin| #leq 3", "|dtbin| #leq 6 (nominal)", "|dtbin| #leq 10"};
  auto panel = [&](int ipad, const char *title, double arr[NDT][NR], double ymax, double refline, const char *refname) {
    cv->cd(ipad);
    TH2D *fr = new TH2D(Form("rdfr%d", ipad), Form("%s;road radius dxy [cm];[#mum]", title), 10, 0.15, 1.95, 10, 0, ymax);
    fr->Draw();
    TLegend *L = new TLegend(0.14, 0.68, 0.60, 0.88); L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextSize(0.030);
    for (int idt = 0; idt < NDT; ++idt)
    {
      TGraph *g = new TGraph;
      for (int ir = 0; ir < NR; ++ir) g->AddPoint(RXY[ir], arr[idt][ir]);
      g->SetLineColor(cols[idt]); g->SetMarkerColor(cols[idt]);
      g->SetLineWidth(2); g->SetMarkerStyle(mks[idt]);
      g->Draw("LP same");
      L->AddEntry(g, dl[idt], "lp");
    }
    if (refline > 0)
    {
      TLine *ln = new TLine(0.15, refline, 1.95, refline);
      ln->SetLineStyle(2); ln->SetLineColor(kGray + 2); ln->Draw();
      L->AddEntry(ln, refname, "l");
    }
    TLine *nom = new TLine(1.2, 0, 1.2, ymax);
    nom->SetLineStyle(3); nom->SetLineColor(kGray + 1); nom->Draw();
    L->Draw();
  };
  panel(1, "GLOBAL median vs road: raw fans out, clipped flat", gr_, 4400, simG, "sim (truth-grouped)");
  cv->cd(1);
  {
    TLegend *L2 = new TLegend(0.47, 0.16, 0.89, 0.30); L2->SetBorderSize(0); L2->SetFillStyle(0); L2->SetTextSize(0.030);
    for (int idt = 0; idt < NDT; ++idt)
    {
      TGraph *g = new TGraph;
      for (int ir = 0; ir < NR; ++ir) g->AddPoint(RXY[ir], gc_[idt][ir]);
      g->SetLineColor(cols[idt]); g->SetMarkerColor(cols[idt]);
      g->SetLineWidth(2); g->SetLineStyle(7); g->SetMarkerStyle(mks[idt]); g->SetMarkerSize(0.8);
      g->Draw("LP same");
      if (idt == 1) L2->AddEntry(g, "same, 3#sigma-clipped (dashed): road-invariant", "lp");
    }
    L2->Draw();
  }
  panel(2, "pickup #sqrt{G_{raw}^{2} #minus G_{clip}^{2}} vs road", pk_, 2600, -1, "");
  panel(3, "LOCAL 4-row median vs road", lo_, 1800, simL, "sim (truth-grouped)");
  // p4: px/track
  cv->cd(4);
  {
    TH2D *fr = new TH2D("rdfr4", "pixels per track vs road (association efficiency);road radius dxy [cm];px / track", 10, 0.15, 1.95, 10, 0, 200);
    fr->Draw();
    TLegend *L = new TLegend(0.55, 0.16, 0.89, 0.36); L->SetBorderSize(0); L->SetFillStyle(0); L->SetTextSize(0.030);
    for (int idt = 0; idt < NDT; ++idt)
    {
      TGraph *g = new TGraph;
      for (int ir = 0; ir < NR; ++ir) g->AddPoint(RXY[ir], pt_[idt][ir]);
      g->SetLineColor(cols[idt]); g->SetMarkerColor(cols[idt]);
      g->SetLineWidth(2); g->SetMarkerStyle(mks[idt]);
      g->Draw("LP same");
      L->AddEntry(g, dl[idt], "lp");
    }
    TLine *nom = new TLine(1.2, 0, 1.2, 200); nom->SetLineStyle(3); nom->SetLineColor(kGray + 1); nom->Draw();
    L->Draw();
  }
  cv->SaveAs(Form("../sim_validation_plots/roadscan_probe_%s.png", ver));
  printf("wrote ../sim_validation_plots/roadscan_probe_%s.png\n", ver);
}
