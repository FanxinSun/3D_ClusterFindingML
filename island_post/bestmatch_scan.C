// bestmatch_scan.C — PERSISTENT best-match window search (was an inline throwaway on
// 2026-07-11 when it found v3.2 f125/L12; provenance rule: producers live as macros).
//
// Finds the sim window (frame x R1 layer, notebook window class phi[0,1) tbin[600,800])
// whose island-level features best match a REAL target window (default e74/L15, the
// funny_shapes notebook view). Features: n islands, <size>, <maxadc>, protrusion
// fraction (size <= 0.75*phisize*zsize), frac(size>=10). Score = sum of squared pulls
// with fixed normalizations (30% relative on n/<size>/<maxadc>; 0.12/0.10 absolute on
// the two fractions). Count prefilter via one branch-address pass (TTree::Draw per
// window crashed ROOT at 4000 windows — hence the accumulator).
//
// Render the winner with: funny_shapes(real, digi, 74, <L>, <f>, 15, 0,1,600,800,
//   ".../funny_shapes_vXX_bestmatch.png", "vXX BEST-MATCH")
#include <TFile.h>
#include <TH1D.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <vector>

void bestmatch_scan(const char *simisl91 = "island91_frames_production_v40b.root",
                    const char *realisl91 = "island91_real.root",
                    int realevent = 74, int reallayer = 15,
                    const char *simdigi = "digi_frames_production_v40b.root",
                    const char *realdigi = "/home/rog/sPHENIX/3D_ClusterFindingML/clusters_seeds_island_79507-0.root_ntuplizer.root")
{
  auto feats = [&](TTree *t, int ev, int lay, double *out) -> bool {
    const char *W = Form("event==%d&&layer==%d&&phi>=0&&phi<1&&tbin>=600&&tbin<=800", ev, lay);
    double n = t->GetEntries(W);
    if (n < 10) return false;
    TH1D h1("h1", "", 300, 0, 300);
    t->Draw("size>>h1", W, "goff");
    TH1D h2("h2", "", 200, 0, 1080);
    t->Draw("maxadc>>h2", W, "goff");
    double prot = t->GetEntries(Form("%s&&phisize>=2&&zsize>=2&&size<=0.75*phisize*zsize", W));
    double big = t->GetEntries(Form("%s&&size>=10", W));
    out[0] = n;
    out[1] = h1.GetMean();
    out[2] = h2.GetMean();
    out[3] = prot / n;
    out[4] = big / n;
    return true;
  };
  TFile *fr = TFile::Open(realisl91);
  TTree *r = (TTree *) fr->Get("ntp_cluster");
  double R[5];
  if (!feats(r, realevent, reallayer, R))
  {
    printf("bestmatch_scan: real target window too empty\n");
    return;
  }
  printf("REAL e%d L%d window: n %.0f | <size> %.2f | <maxadc> %.0f | prot %.2f | frac>=10 %.2f\n",
         realevent, reallayer, R[0], R[1], R[2], R[3], R[4]);
  TFile *fs = TFile::Open(simisl91);
  TTree *s = (TTree *) fs->Get("ntp_cluster");
  // count prefilter: one pass, no per-window Draw
  std::map<int, int> cnt;
  {
    float ev, lay, phi, tb;
    s->SetBranchStatus("*", 0);
    for (const char *b : {"event", "layer", "phi", "tbin"}) s->SetBranchStatus(b, 1);
    s->SetBranchAddress("event", &ev);
    s->SetBranchAddress("layer", &lay);
    s->SetBranchAddress("phi", &phi);
    s->SetBranchAddress("tbin", &tb);
    Long64_t N = s->GetEntries();
    for (Long64_t i = 0; i < N; ++i)
    {
      s->GetEntry(i);
      if (lay < 7 || lay > 22 || phi < 0 || phi >= 1 || tb < 600 || tb > 800) continue;
      cnt[((int) ev) * 100 + (int) lay]++;
    }
    s->SetBranchStatus("*", 1);
  }
  std::vector<std::pair<double, int>> best;
  for (auto &kv : cnt)
  {
    if (kv.second < R[0] * 0.6 || kv.second > R[0] * 1.6) continue;
    double F[5];
    if (!feats(s, kv.first / 100, kv.first % 100, F)) continue;
    double sc = 0;
    sc += std::pow((F[0] - R[0]) / (0.3 * R[0]), 2);
    sc += std::pow((F[1] - R[1]) / (0.3 * R[1]), 2);
    sc += std::pow((F[2] - R[2]) / (0.3 * R[2]), 2);
    sc += std::pow((F[3] - R[3]) / 0.12, 2);
    sc += std::pow((F[4] - R[4]) / 0.10, 2);
    best.push_back({sc, kv.first});
  }
  std::sort(best.begin(), best.end());
  // ---- CRITERIAL feature 6 (user 2026-07-13): count of tu-shaped ~7px clusters ----
  // (asym_showcase pixel test: top tbin row exactly ONE non-corner pixel ON a >=3-wide
  // base; cluster size 6-8). Windows re-islandized from their own pixels (4-connected
  // components) so the test is branch-semantics independent. Applied to the
  // statistical top-30, score += ((n_sim - n_real)/1)^2, then re-ranked.
  auto tucount = [](std::vector<std::array<int, 3>> &px) {  // (pad, tbin, adc)
    if (px.empty()) return 0;
    int n = (int) px.size();
    std::vector<int> comp(n, -1);
    int nc = 0;
    for (int i = 0; i < n; ++i)
    {
      if (comp[i] >= 0) continue;
      std::vector<int> stack = {i};
      comp[i] = nc;
      while (!stack.empty())
      {
        int j = stack.back();
        stack.pop_back();
        for (int k = 0; k < n; ++k)
        {
          if (comp[k] >= 0) continue;
          if (std::abs(px[j][0] - px[k][0]) + std::abs(px[j][1] - px[k][1]) == 1)
          {
            comp[k] = nc;
            stack.push_back(k);
          }
        }
      }
      nc++;
    }
    int ntu = 0;
    for (int c = 0; c < nc; ++c)
    {
      std::vector<std::array<int, 3>> cp;
      for (int i = 0; i < n; ++i)
        if (comp[i] == c) cp.push_back(px[i]);
      if (cp.size() < 6 || cp.size() > 8) continue;
      // 4-orientation tu test (user 2026-07-14: exact, 90-left/right, 180
      // rotations all in scope). Same upright test applied in 4 transforms.
      auto passes = [](const std::vector<std::array<int, 2>> &q) {
        int tmax = -1 << 30;
        for (auto &p : q) tmax = std::max(tmax, p[1]);
        int ntop = 0, xtop = -1 << 30, plo = 1 << 30, phi2 = -1 << 30, nbelow = 0;
        bool onbase = false;
        for (auto &p : q)
        {
          plo = std::min(plo, p[0]);
          phi2 = std::max(phi2, p[0]);
          if (p[1] == tmax) { ntop++; xtop = p[0]; }
          if (p[1] == tmax - 1) nbelow++;
        }
        for (auto &p : q)
          if (p[1] == tmax - 1 && p[0] == xtop) onbase = true;
        return ntop == 1 && nbelow >= 3 && onbase && xtop != plo && xtop != phi2;
      };
      bool istu = false;
      for (int o = 0; o < 4 && !istu; ++o)
      {
        std::vector<std::array<int, 2>> q;
        q.reserve(cp.size());
        for (auto &p : cp)
        {
          if (o == 0) q.push_back({p[0], p[1]});         // up
          else if (o == 1) q.push_back({p[0], -p[1]});   // 180
          else if (o == 2) q.push_back({p[1], p[0]});    // 90 right (+pad)
          else q.push_back({p[1], -p[0]});               // 90 left (-pad)
        }
        istu = passes(q);
      }
      if (!istu) continue;
      // definition v3 (user 2026-07-14): the tank must be DISTINCT - no foreign
      // window pixel within Chebyshev distance 1 of any component pixel
      // (diagonal attachments split by 4-connectivity were counting as tu).
      bool isolated = true;
      for (int i = 0; i < n && isolated; ++i)
      {
        if (comp[i] == c) continue;
        for (auto &p : cp)
          if (std::abs(px[i][0] - p[0]) <= 1 && std::abs(px[i][1] - p[1]) <= 1)
          {
            isolated = false;
            break;
          }
      }
      if (!isolated) continue;
      // definition v4 (user 2026-07-14): BRIGHTNESS UNIFORMITY - a hot core
      // pixel (max/median > 3) reads as two overlapped clusters to the eye,
      // not one clean tank.
      {
        std::vector<int> adcs;
        for (auto &p : cp) adcs.push_back(p[2]);
        std::sort(adcs.begin(), adcs.end());
        double med = adcs[adcs.size() / 2];
        if (med <= 0 || adcs.back() / med > 3.0) continue;
      }
      ntu++;
    }
    return ntu;
  };
  auto winpix = [&](const char *fn, bool isReal, std::map<int, std::vector<std::array<int, 3>>> &out,
                    const std::map<int, char> &want) {
    TFile *fd = TFile::Open(fn);
    TTree *t = (TTree *) fd->Get("ntp_hit");
    float ev, lay, pb, tb, adc, ph;
    t->SetBranchStatus("*", 0);
    for (auto bn : {"event", "layer", "phibin", "adc", "phi"}) t->SetBranchStatus(bn, 1);
    t->SetBranchStatus(isReal ? "tbin" : "zbin", 1);
    t->SetBranchAddress("event", &ev);
    t->SetBranchAddress("layer", &lay);
    t->SetBranchAddress("phibin", &pb);
    t->SetBranchAddress(isReal ? "tbin" : "zbin", &tb);
    t->SetBranchAddress("adc", &adc);
    t->SetBranchAddress("phi", &ph);
    for (Long64_t i = 0; i < t->GetEntries(); ++i)
    {
      t->GetEntry(i);
      if (adc <= 0 || ph < 0 || ph >= 1 || tb < 600 || tb > 800) continue;
      int key = ((int) ev) * 100 + (int) lay;
      if (!want.count(key)) continue;
      out[key].push_back({(int) pb, (int) tb, (int) adc});
    }
    fd->Close();
  };
  std::map<int, char> wantR{{realevent * 100 + reallayer, 1}};
  std::map<int, std::vector<std::array<int, 3>>> pxR;
  winpix(realdigi, true, pxR, wantR);
  int tuR = tucount(pxR[realevent * 100 + reallayer]);
  printf("REAL tu-count (6-8px, window): %d   [CRITERIAL feature 6]\n", tuR);
  std::map<int, char> wantS;
  const int WATCH = 12 * 100 + 15;  // f12 L15: user-judged visually better (2026-07-14)
  wantS[WATCH] = 1;
  for (size_t i = 0; i < 60 && i < best.size(); ++i) wantS[best[i].second] = 1;
  std::map<int, std::vector<std::array<int, 3>>> pxS;
  winpix(simdigi, false, pxS, wantS);
  std::vector<std::pair<double, std::pair<int, int>>> rescored;  // (score, (id, ntu))
  bool watchInPool = false;
  for (size_t i = 0; i < 60 && i < best.size(); ++i)
  {
    if (best[i].second == WATCH) watchInPool = true;
    int ntu = tucount(pxS[best[i].second]);
    double sc = best[i].first + std::pow((double) (ntu - tuR) / 1.0, 2);
    rescored.push_back({sc, {best[i].second, ntu}});
  }
  // user-watch window: force-score even if outside the statistical pool
  if (!watchInPool)
  {
    for (size_t i = 0; i < best.size(); ++i)
      if (best[i].second == WATCH)
      {
        int ntu = tucount(pxS[WATCH]);
        printf("WATCH f12 L15: stat rank %zu, stat score %.2f, tu %d (real %d)\n",
               i + 1, best[i].first, ntu, tuR);
        rescored.push_back({best[i].first + std::pow((double) (ntu - tuR), 2), {WATCH, ntu}});
        break;
      }
  }
  std::sort(rescored.begin(), rescored.end());
  for (size_t i = 0; i < rescored.size(); ++i)
    if (rescored[i].second.first == WATCH)
      printf("WATCH f12 L15 final rank: %zu/%zu (score %.2f, tu %d)\n",
             i + 1, rescored.size(), rescored[i].first, rescored[i].second.second);
  best.clear();
  for (auto &r : rescored) best.push_back({r.first, r.second.first});
  std::map<int, int> tumap;
  for (auto &r : rescored) tumap[r.second.first] = r.second.second;
  // ---- layout score: NON-CRITERIAL (user 2026-07-12) — visual-similarity aid only.
  // The five statistical features above remain the pipeline criteria; this extra
  // term quantifies spatial LAYOUT of the window (what the eye keys on): mean
  // nearest-neighbour centroid distance (window-normalised phi/1.0, tbin/200) and
  // clumpiness (std/mean of island counts over a 4x4 grid). Ad-hoc 15% pulls.
  auto layout = [&](TTree *t, int ev, int lay, double &nnd, double &clump) {
    const char *W = Form("event==%d&&layer==%d&&phi>=0&&phi<1&&tbin>=600&&tbin<=800", ev, lay);
    Long64_t n = t->Draw("phi:tbin", W, "goff");
    if (n < 3) { nnd = -1; clump = -1; return; }
    double *px = t->GetV1(), *pt = t->GetV2();
    double sum = 0;
    double grid[16] = {0};
    for (Long64_t i = 0; i < n; ++i)
    {
      double dmin = 1e9;
      for (Long64_t j = 0; j < n; ++j)
      {
        if (i == j) continue;
        double dx = px[i] - px[j], dt = (pt[i] - pt[j]) / 200.;
        double d = sqrt(dx * dx + dt * dt);
        if (d < dmin) dmin = d;
      }
      sum += dmin;
      int gx = std::min(3, (int) (px[i] * 4)), gt = std::min(3, (int) ((pt[i] - 600) / 50));
      grid[gx * 4 + gt] += 1;
    }
    nnd = sum / n;
    double m = n / 16., v = 0;
    for (double g : grid) v += (g - m) * (g - m);
    clump = sqrt(v / 16.) / m;
  };
  double Rn, Rc;
  layout(r, realevent, reallayer, Rn, Rc);
  printf("REAL layout: nnd %.3f clump %.2f   [NON-CRITERIAL aid]\n", Rn, Rc);
  printf("top matches (score | frame layer | n <size> <mx> prot fr10 | layout):\n");
  double bestlay = 1e9;
  int bestlayid = -1;
  for (size_t i = 0; i < 5 && i < best.size(); ++i)
  {
    double F[5];
    feats(s, best[i].second / 100, best[i].second % 100, F);
    double Sn, Sc;
    layout(s, best[i].second / 100, best[i].second % 100, Sn, Sc);
    double lsc = std::pow((Sn - Rn) / (0.15 * Rn), 2) + std::pow((Sc - Rc) / (0.15 * Rc), 2);
    if (lsc < bestlay) { bestlay = lsc; bestlayid = best[i].second; }
    printf("  %.2f | f%d L%d | %.0f %.2f %.0f %.2f %.2f | tu %d (real %d) | layout %.2f (nnd %.3f clump %.2f)\n",
           best[i].first, best[i].second / 100, best[i].second % 100, F[0], F[1], F[2], F[3], F[4],
           tumap[best[i].second], tuR, lsc, Sn, Sc);
  }
  if (bestlayid >= 0)
  {
    printf("layout-best of top-5 (NON-CRITERIAL): f%d L%d (layout %.2f)\n",
           bestlayid / 100, bestlayid % 100, bestlay);
  }
}
