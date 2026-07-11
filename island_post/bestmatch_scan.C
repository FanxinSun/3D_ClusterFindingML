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

void bestmatch_scan(const char *simisl91 = "island91_frames_production_v33.root",
                    const char *realisl91 = "island91_real.root",
                    int realevent = 74, int reallayer = 15)
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
  printf("top matches (score | frame layer | n <size> <mx> prot fr10):\n");
  for (size_t i = 0; i < 5 && i < best.size(); ++i)
  {
    double F[5];
    feats(s, best[i].second / 100, best[i].second % 100, F);
    printf("  %.2f | f%d L%d | %.0f %.2f %.0f %.2f %.2f\n",
           best[i].first, best[i].second / 100, best[i].second % 100, F[0], F[1], F[2], F[3], F[4]);
  }
}
