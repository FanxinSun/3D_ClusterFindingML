// fieldmeter.C — compact meter for the v5.4 r-phi field campaign (2026-08-05).
// One pass over real (ntp_clus_trk) + sim (island91 ntp_cluster/ntp_truth,
// cls==0 && ntrks==1), ms_real-identical track gates and Kasa+GN circle fit:
//   FM RMS  : per-track circle-RMS medians in the two R_fit windows
//   FM D0   : d0s RMS (|d0s|<8) and median |d0s|
//   FM PROF : 12-bin median d0s vs track-median tbin over [0,960]
// This is the ITERATION instrument; acceptance runs use ms_real.C entries.
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <vector>

namespace FM
{
struct Fit
{
  double a = 0, b = 0, R = 0, rms = 9e9;
  bool ok = false;
};
Fit fitCircle(const std::vector<double> &x, const std::vector<double> &y)
{
  Fit F;
  size_t n = x.size();
  if (n < 3) return F;
  double sx = 0, sy = 0;
  for (size_t i = 0; i < n; ++i) { sx += x[i]; sy += y[i]; }
  double mx = sx / n, my = sy / n;
  double Suu = 0, Suv = 0, Svv = 0, Suuu = 0, Svvv = 0, Suvv = 0, Svuu = 0;
  for (size_t i = 0; i < n; ++i)
  {
    double u = x[i] - mx, v = y[i] - my;
    Suu += u * u; Suv += u * v; Svv += v * v;
    Suuu += u * u * u; Svvv += v * v * v; Suvv += u * v * v; Svuu += v * u * u;
  }
  double det = Suu * Svv - Suv * Suv;
  if (std::fabs(det) < 1e-12) return F;
  double uc = 0.5 * (Svv * (Suuu + Suvv) - Suv * (Svvv + Svuu)) / det;
  double vc = 0.5 * (Suu * (Svvv + Svuu) - Suv * (Suuu + Suvv)) / det;
  double a = uc + mx, b = vc + my;
  double R = 0;
  for (size_t i = 0; i < n; ++i) R += std::hypot(x[i] - a, y[i] - b);
  R /= n;
  for (int it = 0; it < 3; ++it)  // Gauss-Newton refinement
  {
    double JTJ[3][3] = {{0}}, JTr[3] = {0};
    for (size_t i = 0; i < n; ++i)
    {
      double dx = x[i] - a, dy = y[i] - b, di = std::hypot(dx, dy);
      if (di < 1e-9) continue;
      double res = di - R, J0 = -dx / di, J1 = -dy / di, J2 = -1;
      double J[3] = {J0, J1, J2};
      for (int p = 0; p < 3; ++p)
      {
        JTr[p] += J[p] * res;
        for (int q = 0; q < 3; ++q) JTJ[p][q] += J[p] * J[q];
      }
    }
    double A0[3][4];
    for (int p = 0; p < 3; ++p)
    {
      for (int q = 0; q < 3; ++q) A0[p][q] = JTJ[p][q];
      A0[p][3] = -JTr[p];
    }
    for (int c = 0; c < 3; ++c)
    {
      int pv = c;
      for (int rr = c + 1; rr < 3; ++rr)
        if (std::fabs(A0[rr][c]) > std::fabs(A0[pv][c])) pv = rr;
      if (std::fabs(A0[pv][c]) < 1e-12) { pv = -1; break; }
      for (int q = 0; q < 4; ++q) std::swap(A0[c][q], A0[pv][q]);
      for (int rr = 0; rr < 3; ++rr)
      {
        if (rr == c) continue;
        double f = A0[rr][c] / A0[c][c];
        for (int q = 0; q < 4; ++q) A0[rr][q] -= f * A0[c][q];
      }
    }
    if (std::fabs(A0[0][0]) < 1e-12 || std::fabs(A0[1][1]) < 1e-12 || std::fabs(A0[2][2]) < 1e-12) break;
    a += A0[0][3] / A0[0][0];
    b += A0[1][3] / A0[1][1];
    R += A0[2][3] / A0[2][2];
  }
  double s2 = 0;
  for (size_t i = 0; i < n; ++i)
  {
    double d = std::hypot(x[i] - a, y[i] - b) - R;
    s2 += d * d;
  }
  F.a = a; F.b = b; F.R = R;
  F.rms = std::sqrt(s2 / n);
  F.ok = true;
  return F;
}
double med(std::vector<double> v)
{
  if (v.empty()) return 0;
  std::sort(v.begin(), v.end());
  return v[v.size() / 2];
}
}  // namespace FM
using namespace FM;

void fieldmeter(const char *realf = "../clusters_seeds_island_79507-0.root_ntuplizer.root",
                const char *i91 = "island91_frames_production_v54.root",
                const char *tag = "v54a")
{
  struct DT
  {
    std::vector<double> x, y, tb;
    int lmin = 99, lmax = 0;
  };
  std::map<long long, DT> grp[2];

  {  // real
    TFile *f = TFile::Open(realf);
    TTree *t = (TTree *) f->Get("ntp_clus_trk");
    float ev, sid, lay, x, y, tb;
    t->SetBranchStatus("*", 0);
    for (auto b : {"event", "seedID", "layer", "x", "y", "tbin"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("event", &ev);
    t->SetBranchAddress("seedID", &sid);
    t->SetBranchAddress("layer", &lay);
    t->SetBranchAddress("x", &x);
    t->SetBranchAddress("y", &y);
    t->SetBranchAddress("tbin", &tb);
    for (Long64_t i = 0; i < t->GetEntries(); ++i)
    {
      t->GetEntry(i);
      int L = (int) lay;
      if (L < 7 || L > 54) continue;
      DT &T = grp[0][((long long) ev << 24) | (long long) sid];
      T.x.push_back(x); T.y.push_back(y); T.tb.push_back(tb);
      T.lmin = std::min(T.lmin, L); T.lmax = std::max(T.lmax, L);
    }
    f->Close();
  }
  {  // sim
    TFile *f = TFile::Open(i91);
    TTree *c = (TTree *) f->Get("ntp_cluster");
    TTree *u = (TTree *) f->Get("ntp_truth");
    float ev, lay, x, y, tb;
    float gtrk, cls, ntrks;
    c->SetBranchStatus("*", 0);
    for (auto b : {"event", "layer", "x", "y", "tbin"}) c->SetBranchStatus(b, 1);
    c->SetBranchAddress("event", &ev);
    c->SetBranchAddress("layer", &lay);
    c->SetBranchAddress("x", &x);
    c->SetBranchAddress("y", &y);
    c->SetBranchAddress("tbin", &tb);
    u->SetBranchStatus("*", 0);
    for (auto b : {"gtrackID", "cls", "ntrks"}) u->SetBranchStatus(b, 1);
    u->SetBranchAddress("gtrackID", &gtrk);
    u->SetBranchAddress("cls", &cls);
    u->SetBranchAddress("ntrks", &ntrks);
    for (Long64_t i = 0; i < c->GetEntries(); ++i)
    {
      c->GetEntry(i);
      u->GetEntry(i);
      int L = (int) lay;
      if (L < 7 || L > 54) continue;
      if ((int) cls != 0 || (int) ntrks != 1) continue;
      DT &T = grp[1][((long long) ev << 24) | (long long) (int) gtrk];
      T.x.push_back(x); T.y.push_back(y); T.tb.push_back(tb);
      T.lmin = std::min(T.lmin, L); T.lmax = std::max(T.lmax, L);
    }
    f->Close();
  }

  const int NTB2 = 12;
  for (int s = 0; s < 2; ++s)
  {
    std::vector<double> rms05, rmsHi, d0v, tbmv;
    std::vector<std::vector<double>> prof(NTB2);
    for (auto &kv : grp[s])
    {
      DT &T = kv.second;
      if ((int) T.x.size() < 12 || T.lmin > 11 || T.lmax < 50) continue;
      Fit F = fitCircle(T.x, T.y);
      if (!F.ok || F.rms > 0.25 || F.R < 35) continue;
      if (F.R >= 101 && F.R <= 137) rms05.push_back(F.rms * 1e4);
      if (F.R >= 357 && F.R <= 596) rmsHi.push_back(F.rms * 1e4);
      double d0 = std::hypot(F.a, F.b) - F.R;
      std::vector<double> t = T.tb;
      std::sort(t.begin(), t.end());
      double tm = t[t.size() / 2];
      if (std::fabs(d0) < 8)
      {
        d0v.push_back(d0);
        int b = (int) (tm / 960. * NTB2);
        if (b >= 0 && b < NTB2) prof[b].push_back(d0);
      }
    }
    double s2 = 0;
    for (double d : d0v) s2 += d * d;
    std::vector<double> ad;
    for (double d : d0v) ad.push_back(std::fabs(d));
    printf("FM%d RMS %s: n05 %zu med05 %.1f | nHi %zu medHi %.1f\n", s, tag,
           rms05.size(), med(rms05), rmsHi.size(), med(rmsHi));
    printf("FM%d D0  %s: n %zu rms %.3f med|d0| %.3f\n", s, tag, d0v.size(),
           d0v.empty() ? 0 : std::sqrt(s2 / d0v.size()), med(ad));
    printf("FM%d PROF %s:", s, tag);
    for (int b = 0; b < NTB2; ++b)
      printf(" %s", prof[b].size() < 5 ? "nan" : Form("%+.2f", med(prof[b])));
    printf("   (bin centers 40,120,...,920; n:");
    for (int b = 0; b < NTB2; ++b) printf(" %zu", prof[b].size());
    printf(")\n");
  }
}
