// tc_diag.C — one-off diagnostic for truth_circle T3 pathology (scratch)
#include <TFile.h>
#include <TTree.h>
#include <cstdio>
#include <cmath>
#include <map>
#include <vector>
#include <algorithm>

void tc_diag(double pt_lo = 0.45, double pt_hi = 0.55)
{
  TFile *f = TFile::Open("island91_frames_production_v51.root");
  TTree *c = (TTree *) f->Get("ntp_cluster");
  TTree *u = (TTree *) f->Get("ntp_truth");
  float ev, lay, x, y, z, phi, tb, zel, tid, gpt, cls, ntrks;
  c->SetBranchStatus("*", 0);
  for (auto b : {"event", "layer", "x", "y", "z", "phi", "tbin", "zelem"}) c->SetBranchStatus(b, 1);
  c->SetBranchAddress("event", &ev);
  c->SetBranchAddress("layer", &lay);
  c->SetBranchAddress("x", &x);
  c->SetBranchAddress("y", &y);
  c->SetBranchAddress("z", &z);
  c->SetBranchAddress("phi", &phi);
  c->SetBranchAddress("tbin", &tb);
  c->SetBranchAddress("zelem", &zel);
  u->SetBranchStatus("*", 0);
  for (auto b : {"gtrackID", "gpt", "cls", "ntrks"}) u->SetBranchStatus(b, 1);
  u->SetBranchAddress("gtrackID", &tid);
  u->SetBranchAddress("gpt", &gpt);
  u->SetBranchAddress("cls", &cls);
  u->SetBranchAddress("ntrks", &ntrks);
  struct Row { float lay, x, y, z, phi, tb, zel; };
  std::map<long, std::vector<Row>> trks;
  std::map<long, float> tpt;
  for (Long64_t i = 0; i < c->GetEntries(); ++i)
  {
    u->GetEntry(i);
    if (cls != 0 || ntrks != 1) continue;
    if (gpt < pt_lo || gpt > pt_hi) continue;
    c->GetEntry(i);
    long k = (long) ev * 1000000 + (long) tid;
    trks[k].push_back({lay, x, y, z, phi, tb, zel});
    tpt[k] = gpt;
  }
  // simple circle rms via 3-point curvature proxy is overkill; reuse full fit
  auto fitrms = [](std::vector<Row> &v, double &Rout) {
    double Sx = 0, Sy = 0, Sxx = 0, Syy = 0, Sxy = 0, Sxz = 0, Syz = 0, Sz = 0;
    for (auto &r : v)
    { double zz = r.x * r.x + r.y * r.y;
      Sx += r.x; Sy += r.y; Sxx += r.x * r.x; Syy += r.y * r.y; Sxy += r.x * r.y;
      Sxz += r.x * zz; Syz += r.y * zz; Sz += zz; }
    double n = v.size();
    double det = Sxx * (Syy * n - Sy * Sy) - Sxy * (Sxy * n - Sy * Sx) + Sx * (Sxy * Sy - Syy * Sx);
    if (std::fabs(det) < 1e-9) { Rout = 0; return 1e9; }
    double A = (Sxz * (Syy * n - Sy * Sy) - Sxy * (Syz * n - Sy * Sz) + Sx * (Syz * Sy - Syy * Sz)) / det;
    double B = (Sxx * (Syz * n - Sy * Sz) - Sxz * (Sxy * n - Sy * Sx) + Sx * (Sxy * Sz - Syz * Sx)) / det;
    double C = (Sxx * (Syy * Sz - Syz * Sy) - Sxy * (Sxy * Sz - Syz * Sx) + Sxz * (Sxy * Sy - Syy * Sx)) / det;
    double a = A / 2, b = B / 2, R = std::sqrt(std::max(1e-12, C + a * a + b * b));
    double s2 = 0;
    for (auto &r : v) { double d = std::hypot(r.x - a, r.y - b) - R; s2 += d * d; }
    Rout = R;
    return std::sqrt(s2 / n);
  };
  struct Q { long k; double rms, R; int n; };
  std::vector<Q> qs;
  for (auto &kv : trks)
  {
    auto &v = kv.second;
    int lmin = 99, lmax = 0;
    for (auto &r : v) { lmin = std::min(lmin, (int) r.lay); lmax = std::max(lmax, (int) r.lay); }
    if ((int) v.size() < 12 || lmin > 11 || lmax < 50) continue;
    double R; double rms = fitrms(v, R);
    qs.push_back({kv.first, rms, R, (int) v.size()});
  }
  std::sort(qs.begin(), qs.end(), [](const Q &a, const Q &b) { return a.rms > b.rms; });
  printf("fitted %zu tracks; worst 3 and median 3:\n", qs.size());
  auto dump = [&](Q &q) {
    auto v = trks[q.k];
    std::sort(v.begin(), v.end(), [](const Row &a, const Row &b) { return a.lay < b.lay; });
    printf("--- ev %ld tid %ld pt %.3f n %d rms(cm) %.3f Rfit %.1f\n",
           q.k / 1000000, q.k % 1000000, tpt[q.k], q.n, q.rms, q.R);
    for (auto &r : v)
      printf("  L%2.0f zel %1.0f phi %+7.4f x %+8.3f y %+8.3f z %+8.2f tb %6.1f\n",
             r.lay, r.zel, r.phi, r.x, r.y, r.z, r.tb);
  };
  for (int i = 0; i < 3 && i < (int) qs.size(); ++i) dump(qs[i]);
  printf("=== median region ===\n");
  for (size_t i = qs.size() / 2; i < qs.size() / 2 + 3 && i < qs.size(); ++i) dump(qs[i]);
}
