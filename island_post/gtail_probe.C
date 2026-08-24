// gtail_probe.C — why does the REAL GLOBAL whole-track pixel fit (nf_digipix)
// have a long tail (med 2358 um, plateau to 6 mm) while LOCAL is 1205 um and
// the SAME seeds fitted at cluster level give 780 um?
// Method: rebuild the exact nf_digipix real groups (road: same layer,
// dxy<1.2 cm, |dtbin|<=6, hitID-dedup, first-seed-wins, ev44 vetoed) and the
// sim digi truth groups, fit with the SAME estimator (MNF::fitCircle), then
// per track:
//   - robust re-fit: 3 iterations dropping |res-med| > 3*1.4826*MAD (floor
//     500 um) -> rms_trim, fdrop  [contamination collapses, physics survives]
//   - exact pixel-weighted decomposition of rms^2 into LF (row-mean residuals
//     = long-range shape) + HF (within-row scatter = charge-cloud width)
//   - per-row peak-to-peak residual (two-blob rows = a second track in road)
//   - fitted R, span, d0, side/sector/layer cell means (alignment coherence)
// Classes (thresholds printed): CONTAM rms_trim<1.3*med_trim ; BENT LF>HF ;
// WIDE rest. usage: root -l -b -q 'gtail_probe.C+(...)'
// out: ../sim_validation_plots/gtail_probe_<ver>.png + gtail_probe_<ver>.txt
#include "../sim_validation_plots/src/ms_nofinder.C"
#include <THStack.h>
#include <TH2D.h>
#include <TLine.h>
#include <TGraph.h>
#include <cstdint>

namespace GTP
{
struct GrpX
{
  std::vector<double> x, y, r;
  std::vector<float> tb;
  std::vector<uint8_t> lay, side;
  int ev = -1;
};
struct Diag
{
  int ev, npx, nrows, cls, gidx;  // cls 1=trims-clean 2=LF-survivor 3=HF-survivor
  double rms0, R0, span, rmsT, RT, d0T, fdrop, LF, HF, rowPP, keptfrac;
};
double medv(std::vector<double> v) { return MNF::med(v); }

// robust trim: iterate MAD-clip about the median residual, refit
bool trimFit(const GrpX &G, const std::vector<int> &idx0, MNF::Fit &F,
             std::vector<int> &kept, int iters = 3, double floorcm = 0.05)
{
  kept = idx0;
  for (int it = 0; it < iters; ++it)
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
    double thr = std::max(3 * 1.4826 * medv(tmp), floorcm);
    std::vector<int> nk;
    for (size_t k = 0; k < kept.size(); ++k)
      if (std::fabs(res[k] - md) <= thr) nk.push_back(kept[k]);
    if (nk.size() == kept.size()) break;
    if ((int) nk.size() < 12) break;   // don't trim below the bar
    kept.swap(nk);
  }
  std::vector<double> X, Y;
  for (int i : kept) { X.push_back(G.x[i]); Y.push_back(G.y[i]); }
  F = MNF::fitCircle(X, Y);
  return F.ok;
}
}  // namespace GTP

void gtail_probe(const char *digif = "digi_frames_production_v6.root", int nsim = 60,
                 const char *realf = "../clusters_seeds_island_79507-0.root_ntuplizer.root",
                 const char *ver = "v6")
{
  using namespace GTP;
  gROOT->SetBatch(1);
  // ---- row radii (same table as nf_digipix) ------------------------------
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

  // ---- real groups: EXACT nf_digipix road, extended payload ---------------
  std::vector<GrpX> gr[2];                        // [0]=real [1]=sim
  std::map<int, long> evocc;                      // real TPC px per event
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
      evocc[(int) ev]++;
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
        G.tb.push_back(tb); G.lay.push_back((uint8_t) lay); G.side.push_back((uint8_t) zel);
        break;
      }
    }
    f->Close();
    printf("real: %zu seeds grouped\n", gr[0].size());
  }
  // ---- sim digi groups (nf_digipix selection + payload) --------------------
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
      G.tb.push_back(tb); G.lay.push_back((uint8_t) lay); G.side.push_back((uint8_t) zel);
    }
    f->Close();
    printf("sim: %zu truth groups in %d frames\n", gr[1].size(), nsim);
  }

  // ---- per-track diagnostics ----------------------------------------------
  std::vector<Diag> dg[2];
  // cell coherence: mean row residual per (side,sector,layer), tracks rmsT<4000um
  std::map<int, std::pair<double, long>> cell[2];
  for (int s = 0; s < 2; ++s)
    for (auto &G : gr[s])
    {
      MNF::Grp B; B.x = G.x; B.y = G.y; B.r = G.r;
      MNF::Fit F0;
      if (!MNF::fitBar(B, F0)) continue;
      Diag D; D.ev = G.ev; D.npx = (int) G.x.size();
      D.gidx = (int) (&G - &gr[s][0]);
      D.rms0 = F0.rms * 1e4; D.R0 = F0.R;
      double rlo = 1e9, rhi = 0;
      for (double q : G.r) { rlo = std::min(rlo, q); rhi = std::max(rhi, q); }
      D.span = rhi - rlo;
      std::vector<int> all(G.x.size()); for (size_t k = 0; k < all.size(); ++k) all[k] = (int) k;
      MNF::Fit FT; std::vector<int> kept;
      if (!trimFit(G, all, FT, kept)) continue;
      D.rmsT = FT.rms * 1e4; D.RT = FT.R; D.fdrop = 1. - (double) kept.size() / G.x.size();
      D.keptfrac = (double) kept.size() / G.x.size();
      D.d0T = std::fabs(std::hypot(FT.a, FT.b) - FT.R);
      // decomposition + row structure on the TRIMMED set (LF/HF) and RAW set (rowPP)
      std::map<int, std::vector<double>> rowsT, rowsR;
      std::map<int, std::pair<double, int>> rowmeta;  // sum tb , n  (unused now; side/sector below)
      for (int i : kept)
      {
        int L = nearRow(G.r[i]); if (L < 0) continue;
        rowsT[L].push_back(std::hypot(G.x[i] - FT.a, G.y[i] - FT.b) - FT.R);
      }
      for (size_t i = 0; i < G.x.size(); ++i)
      {
        int L = nearRow(G.r[i]); if (L < 0) continue;
        rowsR[L].push_back(std::hypot(G.x[i] - FT.a, G.y[i] - FT.b) - FT.R);
      }
      double lf2 = 0, hf2 = 0; long nT = 0;
      for (auto &kv : rowsT)
      {
        double m = 0; for (double q : kv.second) m += q; m /= kv.second.size();
        lf2 += m * m * kv.second.size();
        for (double q : kv.second) hf2 += (q - m) * (q - m);
        nT += (long) kv.second.size();
      }
      D.LF = nT ? std::sqrt(lf2 / nT) * 1e4 : 0;
      D.HF = nT ? std::sqrt(hf2 / nT) * 1e4 : 0;
      D.nrows = (int) rowsR.size();
      D.rowPP = 0;
      for (auto &kv : rowsR)
        if (kv.second.size() >= 2)
        {
          double lo = 1e9, hi = -1e9;
          for (double q : kv.second) { lo = std::min(lo, q); hi = std::max(hi, q); }
          D.rowPP = std::max(D.rowPP, (hi - lo) * 1e4);
        }
      D.cls = -1;
      dg[s].push_back(D);
      // cell accumulation (trimmed rows, per side/sector of the row's pixels)
      if (D.rmsT < 4000)
      {
        std::map<int, std::vector<std::pair<int,double>>> dummy;
        std::map<int, std::vector<int>> rowpix;
        for (int i : kept) { int L = nearRow(G.r[i]); if (L >= 0) rowpix[L].push_back(i); }
        for (auto &kv : rowpix)
        {
          double m = 0, mphi = 0; int nside = 0;
          for (int i : kv.second)
          {
            m += std::hypot(G.x[i] - FT.a, G.y[i] - FT.b) - FT.R;
            mphi += std::atan2(G.y[i], G.x[i]);
            nside += G.side[i];
          }
          m /= kv.second.size(); mphi /= kv.second.size();
          int sec = ((int) std::floor((mphi + 2 * M_PI) / (M_PI / 6))) % 12;
          int sd = (2 * nside >= (int) kv.second.size()) ? 1 : 0;
          auto &C = cell[s][(sd * 12 + sec) * 100 + kv.first];
          C.first += m; C.second++;
        }
      }
    }

  // ---- classification ------------------------------------------------------
  std::vector<double> allT[2];
  for (int s = 0; s < 2; ++s) for (auto &D : dg[s]) allT[s].push_back(D.rmsT);
  double medT[2] = {medv(allT[0]), medv(allT[1])};
  for (int s = 0; s < 2; ++s)
    for (auto &D : dg[s])
    {
      if (D.rmsT < 1.3 * medT[s]) D.cls = 1;        // collapses to clean -> outliers were the story
      else if (D.LF > D.HF) D.cls = 2;              // long-range shape dominates
      else D.cls = 3;                                // wide scatter everywhere
    }

  // ---- ledger ---------------------------------------------------------------
  FILE *fo = fopen(Form("gtail_probe_%s.txt", ver), "w");
  auto P = [&](const char *fmt, ...) {
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    va_start(ap, fmt); vfprintf(fo, fmt, ap); va_end(ap);
  };
  auto Q = [&](std::vector<double> v, double q) -> double {
    if (v.empty()) return 0; std::sort(v.begin(), v.end());
    return v[std::min(v.size() - 1, (size_t)(q * v.size()))];
  };
  P("[gtail_probe %s] anatomy of the GLOBAL whole-track pixel-fit tail (nf_digipix population)\n", ver);
  for (int s = 0; s < 2; ++s)
  {
    std::vector<double> r0, rt, fd, lf, hf;
    for (auto &D : dg[s]) { r0.push_back(D.rms0); rt.push_back(D.rmsT); fd.push_back(D.fdrop); lf.push_back(D.LF); hf.push_back(D.HF); }
    P("  %s: %zu tracks | RAW rms med %.0f (75%% %.0f, 90%% %.0f, 95%% %.0f) | >3mm %.1f%% | >=6mm %.1f%%\n",
      s ? "sim " : "real", dg[s].size(), Q(r0, .50), Q(r0, .75), Q(r0, .90), Q(r0, .95),
      100. * std::count_if(r0.begin(), r0.end(), [](double q) { return q > 3000; }) / r0.size(),
      100. * std::count_if(r0.begin(), r0.end(), [](double q) { return q >= 5999; }) / r0.size());
    P("        TRIMMED rms med %.0f | med fdrop %.1f%% | med LF %.0f HF %.0f um\n",
      Q(rt, .50), 100 * Q(fd, .50), Q(lf, .50), Q(hf, .50));
  }
  std::vector<double> raw0, raw1;
  for (auto &D : dg[0]) raw0.push_back(D.rms0);
  for (auto &D : dg[1]) raw1.push_back(D.rms0);
  double medR0 = medv(raw0), medR1 = medv(raw1);
  P("  data/MC GLOBAL: raw %.2f -> trimmed %.2f   (trim = 3x MAD clip, 3 iters, floor 500 um)\n",
    medR0 / medR1, medT[0] / medT[1]);
  {
    long nfd = std::count_if(dg[0].begin(), dg[0].end(), [](const Diag &D) { return D.fdrop > 0.10; });
    P("  real tracks with >10%% of pixels clipped: %.1f%% (sim: %.1f%%)\n",
      100. * nfd / dg[0].size(),
      100. * std::count_if(dg[1].begin(), dg[1].end(), [](const Diag &D) { return D.fdrop > 0.10; }) / (double) dg[1].size());
  }
  // genuine long-range share: trimmed LF by radial span
  for (int s2 = 0; s2 < 2; ++s2)
  {
    std::vector<double> lfb[3];
    for (auto &D : dg[s2])
    { int b = D.span < 25 ? 0 : D.span < 35 ? 1 : 2; lfb[b].push_back(D.LF); }
    P("  %s trimmed LF med by span [15,25)/[25,35)/[35,45] cm: %.0f / %.0f / %.0f um\n",
      s2 ? "sim " : "real", medv(lfb[0]), medv(lfb[1]), medv(lfb[2]));
  }
  // tail composition
  for (int s = 0; s < 2; ++s)
  {
    long n[4] = {0, 0, 0, 0}; std::vector<double> fdc[4], rtc[4], Rc[4], ppc[4];
    long ntail = 0;
    for (auto &D : dg[s])
      if (D.rms0 > 3000)
      { ntail++; n[D.cls]++; fdc[D.cls].push_back(D.fdrop); rtc[D.cls].push_back(D.rmsT); Rc[D.cls].push_back(D.R0); ppc[D.cls].push_back(D.rowPP); }
    P("  %s TAIL (raw rms>3mm): %ld tracks (%.1f%%)\n", s ? "sim " : "real", ntail, 100. * ntail / std::max((size_t) 1, dg[s].size()));
    const char *cn[4] = {"", "CONTAM (trims to clean)", "BENT   (LF>HF survives)", "WIDE   (HF survives)"};
    for (int c = 1; c <= 3; ++c)
      P("    %-24s %5ld (%4.1f%% of tail) | med fdrop %4.1f%% rmsT %4.0f um R %.0f cm rowPP %4.0f um\n",
        cn[c], n[c], ntail ? 100. * n[c] / ntail : 0, 100 * medv(fdc[c]), medv(rtc[c]), medv(Rc[c]), medv(ppc[c]));
  }
  // cell coherence
  double cellRMS[2] = {0, 0};
  for (int s = 0; s < 2; ++s)
  {
    std::vector<double> cm;
    for (auto &kv : cell[s]) if (kv.second.second >= 8) cm.push_back(kv.second.first / kv.second.second * 1e4);
    double m = 0, s2 = 0; for (double q : cm) m += q; if (!cm.empty()) m /= cm.size();
    for (double q : cm) s2 += (q - m) * (q - m);
    cellRMS[s] = cm.size() > 1 ? std::sqrt(s2 / (cm.size() - 1)) : 0;
    P("  %s cell coherence (side,sector,layer row-mean residual, >=8 tracks): %zu cells, RMS %.0f um\n",
      s ? "sim " : "real", cm.size(), cellRMS[s]);
  }
  // figure-facing computed stats
  double medFd[2], tailCleanPct = 0, rowPPlo = 0, rowPPhi = 0, tTail[3] = {0, 0, 0};
  for (int s = 0; s < 2; ++s)
  { std::vector<double> fd; for (auto &D : dg[s]) fd.push_back(D.fdrop); medFd[s] = 100 * Q(fd, .5); }
  {
    long nt = 0, nc = 0; std::vector<double> pp[4];
    for (auto &D : dg[0]) if (D.rms0 > 3000) { nt++; if (D.cls == 1) nc++; pp[D.cls].push_back(D.rowPP); }
    tailCleanPct = nt ? 100. * nc / nt : 0;
    rowPPlo = 1e9; rowPPhi = 0;
    for (int c = 1; c <= 3; ++c) { double m2 = medv(pp[c]); rowPPlo = std::min(rowPPlo, m2); rowPPhi = std::max(rowPPhi, m2); }
  }
  // occupancy correlation (real)
  {
    std::vector<std::pair<long, double>> eo; std::vector<std::pair<long, double>> eoT;
    for (auto &D : dg[0]) { eo.push_back({evocc[D.ev], D.rms0}); eoT.push_back({evocc[D.ev], D.rmsT}); }
    std::sort(eo.begin(), eo.end()); std::sort(eoT.begin(), eoT.end());
    long n3 = eo.size() / 3;
    std::vector<double> vt[3];
    for (int tz = 0; tz < 3; ++tz)
    { long a = tz * n3, b = (tz == 2) ? eoT.size() : a + n3; for (long i = a; i < b; ++i) vt[tz].push_back(eoT[i].second); }
    for (int tzone = 0; tzone < 3; ++tzone)
    {
      std::vector<double> v; long a = tzone * n3, b = (tzone == 2) ? eo.size() : a + n3;
      for (long i = a; i < b; ++i) v.push_back(eo[i].second);
      P("  real occupancy tercile %d (ev px %ldk-%ldk): med raw rms %.0f um | raw>3mm %.1f%% | trimmed>3mm %.1f%%\n",
        tzone + 1, eo[a].first / 1000, eo[b - 1].first / 1000, Q(v, .5),
        100. * std::count_if(v.begin(), v.end(), [](double q) { return q > 3000; }) / v.size(),
        100. * std::count_if(vt[tzone].begin(), vt[tzone].end(), [](double q) { return q > 3000; }) / std::max((size_t) 1, vt[tzone].size()));
    }
  }
  // worst list
  {
    std::vector<const Diag *> w;
    for (auto &D : dg[0]) w.push_back(&D);
    std::sort(w.begin(), w.end(), [](const Diag *a, const Diag *b) { return a->rms0 > b->rms0; });
    P("  real worst 12: ev/npx/R0/span | rms0 -> rmsT (fdrop%%) LF/HF rowPP cls\n");
    for (int i = 0; i < 12 && i < (int) w.size(); ++i)
      P("    ev%3d n%4d R%6.0f sp%3.0f | %5.0f -> %4.0f (%4.1f%%) %4.0f/%4.0f %5.0f %s\n",
        w[i]->ev, w[i]->npx, w[i]->R0, w[i]->span, w[i]->rms0, w[i]->rmsT, 100 * w[i]->fdrop,
        w[i]->LF, w[i]->HF, w[i]->rowPP, w[i]->cls == 1 ? "CONTAM" : w[i]->cls == 2 ? "BENT" : "WIDE");
  }
  // exemplar: largest raw->trim collapse among trims-to-clean tail tracks
  const Diag *EX = nullptr;
  for (auto &D : dg[0])
    if (D.rms0 > 3000 && D.cls == 1 && (!EX || D.rms0 - D.rmsT > EX->rms0 - EX->rmsT)) EX = &D;
  if (EX) P("  exemplar (panel 3): ev %d, %d px, R0 %.0f cm, rms %0.f -> %.0f um, fdrop %.1f%%\n",
            EX->ev, EX->npx, EX->R0, EX->rms0, EX->rmsT, 100 * EX->fdrop);
  fclose(fo);
  printf("wrote gtail_probe_%s.txt\n", ver);

  // ---- figure ---------------------------------------------------------------
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.045);
  TCanvas *cv = new TCanvas("cvgt", "gtail", 1700, 1250);
  cv->Divide(2, 2);
  // [1] raw vs trimmed distributions
  cv->cd(1);
  TH1D *h[4];
  const char *hn[4] = {"real raw (road, no cleaning)", "real after 3#sigma clip", "sim digi raw", "sim after same clip"};
  int hc[4] = {kBlack, kRed + 1, kBlue + 1, kAzure + 6};
  for (int k = 0; k < 4; ++k)
  {
    h[k] = new TH1D(Form("gt_h%d", k), "GLOBAL whole-track pixel fit: raw vs 3#sigma-clipped;per-fit circle RMS [#mum];fits (unit area)", 60, 0, 6000);
    int s2 = k / 2;
    for (auto &D : dg[s2]) h[k]->Fill(std::min(k % 2 ? D.rmsT : D.rms0, 5999.));
    h[k]->Scale(1. / h[k]->Integral());
    h[k]->SetLineColor(hc[k]); h[k]->SetLineWidth(2);
    if (k == 1) h[k]->SetLineStyle(7);
    if (k == 3) h[k]->SetLineStyle(3);
  }
  h[0]->SetMaximum(1.45 * std::max(h[2]->GetMaximum(), h[3]->GetMaximum()));
  h[0]->Draw("hist"); for (int k = 1; k < 4; ++k) h[k]->Draw("hist same");
  TLegend *L1 = new TLegend(0.40, 0.62, 0.89, 0.87);
  L1->SetBorderSize(0); L1->SetFillStyle(0); L1->SetTextSize(0.030);
  double medsv[4] = {medR0, medT[0], medR1, medT[1]};
  for (int k = 0; k < 4; ++k) L1->AddEntry(h[k], Form("%s, med %.0f #mum", hn[k], medsv[k]), "l");
  L1->Draw();
  TLatex tx; tx.SetNDC(); tx.SetTextSize(0.034);
  tx.DrawLatex(0.40, 0.55, Form("data/MC %.2f #rightarrow %.2f", medR0 / medR1, medT[0] / medT[1]));
  tx.SetTextSize(0.027);
  tx.DrawLatex(0.40, 0.50, Form("clip cost: real %.1f%%, sim %.1f%% of pixels (median)", medFd[0], medFd[1]));
  // [2] real raw stacked by class
  cv->cd(2);
  THStack *st = new THStack("gt_st", "real raw RMS by what survives the clip;per-fit circle RMS [#mum];tracks");
  TH1D *hs[4];
  const char *sn[3] = {"trims to clean (association outliers or already clean)",
                       "long-range shape survives (LF>HF: blends, kinks, bent)",
                       "wide scatter survives (HF: blends, loopers)"};
  int sc2[3] = {kAzure - 9, kOrange - 3, kRed - 7};
  for (int c = 1; c <= 3; ++c)
  {
    hs[c] = new TH1D(Form("gt_s%d", c), "", 60, 0, 6000);
    for (auto &D : dg[0]) if (D.cls == c) hs[c]->Fill(std::min(D.rms0, 5999.));
    hs[c]->SetFillColor(sc2[c - 1]); hs[c]->SetLineColor(sc2[c - 1]);
    st->Add(hs[c]);
  }
  st->Draw("hist");
  TLegend *L2 = new TLegend(0.33, 0.68, 0.89, 0.87);
  L2->SetBorderSize(0); L2->SetFillStyle(0); L2->SetTextSize(0.027);
  for (int c = 1; c <= 3; ++c) L2->AddEntry(hs[c], sn[c - 1], "f");
  L2->Draw();
  tx.SetTextSize(0.030);
  tx.DrawLatex(0.33, 0.62, Form("tail (>3 mm): %.0f%% trims to clean; survivors carry", tailCleanPct));
  tx.DrawLatex(0.33, 0.57, Form("%.0f-%.0f mm two-blob rows #rightarrow unresolved second tracks", rowPPlo / 1000, rowPPhi / 1000));
  // [3] exemplar residual vs r
  cv->cd(3);
  if (EX)
  {
    const GrpX &G = gr[0][EX->gidx];
    std::vector<int> all((size_t) G.x.size()); for (size_t k = 0; k < all.size(); ++k) all[k] = (int) k;
    MNF::Fit FT; std::vector<int> kept;
    trimFit(G, all, FT, kept);
    std::set<int> K(kept.begin(), kept.end());
    TGraph *gk = new TGraph, *gd = new TGraph;
    for (size_t i = 0; i < G.x.size(); ++i)
    {
      double res = (std::hypot(G.x[i] - FT.a, G.y[i] - FT.b) - FT.R) * 10;  // mm
      (K.count((int) i) ? gk : gd)->AddPoint(G.r[i], res);
    }
    TH2D *fr = new TH2D("gt_fr", Form("exemplar CONTAM track (ev %d): the road picked up a second track;pad-row radius r [cm];residual to clipped fit [mm]", EX->ev),
                        10, 28, 80, 10, -12, 25);
    fr->Draw();
    TLine z0; z0.SetLineStyle(2); z0.SetLineColor(kGray + 2); z0.DrawLine(28, 0, 80, 0);
    gk->SetMarkerStyle(20); gk->SetMarkerSize(0.6); gk->SetMarkerColor(kBlack); gk->Draw("P same");
    gd->SetMarkerStyle(24); gd->SetMarkerSize(0.8); gd->SetMarkerColor(kRed + 1); gd->Draw("P same");
    TLegend *L3 = new TLegend(0.14, 0.72, 0.66, 0.87);
    L3->SetBorderSize(0); L3->SetFillStyle(0); L3->SetTextSize(0.030);
    L3->AddEntry(gk, Form("kept pixels (%zu): the track, RMS %.0f #mum", kept.size(), EX->rmsT), "p");
    L3->AddEntry(gd, Form("clipped pixels (%d): inside the 1.2 cm road", (int) (G.x.size() - kept.size())), "p");
    L3->Draw();
    tx.SetTextSize(0.030);
    tx.DrawLatex(0.14, 0.66, Form("raw fit through everything: RMS %.0f #mum", EX->rms0));
  }
  // [4] tail fraction vs occupancy
  cv->cd(4);
  {
    std::vector<std::pair<long, double>> eo, eoT;
    for (auto &D : dg[0]) { eo.push_back({evocc[D.ev], D.rms0}); eoT.push_back({evocc[D.ev], D.rmsT}); }
    std::sort(eo.begin(), eo.end()); std::sort(eoT.begin(), eoT.end());
    long n3 = eo.size() / 3;
    TH1D *hr = new TH1D("gt_or", "tail fraction vs event occupancy (real);event-occupancy tercile;fraction of fits with RMS > 3 mm", 3, 0.5, 3.5);
    TH1D *ht = (TH1D *) hr->Clone("gt_ot");
    for (int tz = 0; tz < 3; ++tz)
    {
      long a = tz * n3, b = (tz == 2) ? eo.size() : a + n3;
      double fr2 = 0, ft = 0;
      for (long i = a; i < b; ++i) { if (eo[i].second > 3000) fr2++; if (eoT[i].second > 3000) ft++; }
      hr->SetBinContent(tz + 1, fr2 / (b - a)); ht->SetBinContent(tz + 1, ft / (b - a)); tTail[tz] = 100. * ft / (b - a);
      hr->GetXaxis()->SetBinLabel(tz + 1, Form("%s (%ld-%ldk px)", tz == 0 ? "quiet" : tz == 1 ? "mid" : "busy", eo[a].first / 1000, eo[b - 1].first / 1000));
    }
    hr->SetFillColor(kGray + 1); hr->SetLineColor(kGray + 3);
    ht->SetFillColor(kRed - 7); ht->SetLineColor(kRed + 1);
    hr->SetMaximum(0.45); hr->SetMinimum(0); hr->SetBarWidth(0.42); hr->SetBarOffset(0.06); hr->GetXaxis()->SetLabelSize(0.045);
    ht->SetBarWidth(0.42); ht->SetBarOffset(0.52);
    hr->Draw("bar"); ht->Draw("bar same");
    TLegend *L4 = new TLegend(0.14, 0.74, 0.72, 0.87);
    L4->SetBorderSize(0); L4->SetFillStyle(0); L4->SetTextSize(0.030);
    L4->AddEntry(hr, "raw fits: tail grows with pileup", "f");
    L4->AddEntry(ht, Form("after 3#sigma clip: %.0f-%.0f%% remnant", tTail[0], tTail[2]), "f");
    L4->Draw();
    tx.SetTextSize(0.028);
    tx.DrawLatex(0.14, 0.66, Form("cell coherence (side,sector,layer): real %.0f #mum vs sim %.0f #mum", cellRMS[0], cellRMS[1]));
    tx.DrawLatex(0.14, 0.61, "= the genuine (small) alignment/field share, not the tail");
  }
  cv->SaveAs(Form("../sim_validation_plots/gtail_probe_%s.png", ver));
  printf("wrote ../sim_validation_plots/gtail_probe_%s.png\n", ver);
}
