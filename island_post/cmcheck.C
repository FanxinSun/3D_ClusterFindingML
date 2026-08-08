// cmcheck.C — field-model shape checks (2026-08-05, answering the circle-fit
// thread's two v5.4b questions):
//   CMJUMP: side-to-side rphi discontinuity at the central membrane —
//     circle fit on side-A clusters of a side-crossing track, coherent
//     signed rphi offset of its side-B clusters; per-track J, ensemble
//     median|J| and robust width.
//   HALFARC: curvature-noise meter — R_in (L<=30) vs R_out (L>30) circle
//     fits per full-crosser; dk = ln(R_in/R_out); robust width of dk.
// Runs identically on real (ntp_clus_trk, side = zelem) and sim island91
// (ntp_cluster+ntp_truth cls==0&&ntrks==1, side = zelem).
#include <TFile.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <vector>

namespace CMC
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
  if (n < 4) return F;
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
  F.a = uc + mx; F.b = vc + my;
  double R = 0;
  for (size_t i = 0; i < n; ++i) R += std::hypot(x[i] - F.a, y[i] - F.b);
  F.R = R / n;
  double s2 = 0;
  for (size_t i = 0; i < n; ++i)
  {
    double d = std::hypot(x[i] - F.a, y[i] - F.b) - F.R;
    s2 += d * d;
  }
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
double q68(std::vector<double> v)  // robust width: 68th pct of |x - median|
{
  if (v.size() < 5) return 0;
  double m = med(v);
  std::vector<double> a;
  for (double x : v) a.push_back(std::fabs(x - m));
  std::sort(a.begin(), a.end());
  return a[(size_t) (0.68 * (a.size() - 1))];
}
}  // namespace CMC
using namespace CMC;

void cmcheck(const char *file, int isSim, const char *tag)
{
  struct DT
  {
    std::vector<double> x, y, L, sd;
    int lmin = 99, lmax = 0;
  };
  std::map<long long, DT> grp;
  TFile *f = TFile::Open(file);
  TTree *t = (TTree *) f->Get(isSim ? "ntp_cluster" : "ntp_clus_trk");
  TTree *u = isSim ? (TTree *) f->Get("ntp_truth") : nullptr;
  float ev, sid = 0, lay, x, y, ze;
  float gtrk = 0, cls = 0, ntrks = 0;
  t->SetBranchStatus("*", 0);
  for (auto b : {"event", "layer", "x", "y", "zelem"}) t->SetBranchStatus(b, 1);
  t->SetBranchAddress("event", &ev);
  t->SetBranchAddress("layer", &lay);
  t->SetBranchAddress("x", &x);
  t->SetBranchAddress("y", &y);
  t->SetBranchAddress("zelem", &ze);
  if (!isSim)
  {
    t->SetBranchStatus("seedID", 1);
    t->SetBranchAddress("seedID", &sid);
  }
  else
  {
    u->SetBranchStatus("*", 0);
    for (auto b : {"gtrackID", "cls", "ntrks"}) u->SetBranchStatus(b, 1);
    u->SetBranchAddress("gtrackID", &gtrk);
    u->SetBranchAddress("cls", &cls);
    u->SetBranchAddress("ntrks", &ntrks);
  }
  for (Long64_t i = 0; i < t->GetEntries(); ++i)
  {
    t->GetEntry(i);
    if (isSim)
    {
      u->GetEntry(i);
      if ((int) cls != 0 || (int) ntrks != 1) continue;
    }
    int L = (int) lay;
    if (L < 7 || L > 54) continue;
    DT &T = grp[((long long) ev << 24) | (long long) (isSim ? (int) gtrk : (int) sid)];
    T.x.push_back(x); T.y.push_back(y); T.L.push_back(L); T.sd.push_back(ze);
    T.lmin = std::min(T.lmin, L); T.lmax = std::max(T.lmax, L);
  }
  f->Close();

  std::vector<double> J, DK, DK05;
  long ncross = 0, nfull = 0;
  for (auto &kv : grp)
  {
    DT &T = kv.second;
    size_t n = T.x.size();
    if ((int) n < 12 || T.lmin > 11 || T.lmax < 50) continue;
    Fit Fall = fitCircle(T.x, T.y);
    if (!Fall.ok || Fall.rms > 0.4 || Fall.R < 35) continue;
    nfull++;
    // HALFARC: inner vs outer half fits
    std::vector<double> xi, yi, xo, yo;
    for (size_t i = 0; i < n; ++i)
      (T.L[i] <= 30 ? xi : xo).push_back(T.x[i]),
      (T.L[i] <= 30 ? yi : yo).push_back(T.y[i]);
    if (xi.size() >= 6 && xo.size() >= 6)
    {
      Fit Fi = fitCircle(xi, yi), Fo = fitCircle(xo, yo);
      if (Fi.ok && Fo.ok && Fi.R > 20 && Fo.R > 20 && Fi.R < 5000 && Fo.R < 5000)
      {
        double dk = std::log(Fi.R / Fo.R);
        DK.push_back(dk);
        if (Fall.R >= 101 && Fall.R <= 137) DK05.push_back(dk);
      }
    }
    // CMJUMP: side-crossers only
    std::vector<double> xa, ya, xb, yb;
    for (size_t i = 0; i < n; ++i)
      (T.sd[i] < 0.5 ? xa : xb).push_back(T.x[i]),
      (T.sd[i] < 0.5 ? ya : yb).push_back(T.y[i]);
    if (xa.size() >= 6 && xb.size() >= 3)
    {
      ncross++;
      const std::vector<double> &fx = xa.size() >= xb.size() ? xa : xb;
      const std::vector<double> &fy = xa.size() >= xb.size() ? ya : yb;
      const std::vector<double> &px = xa.size() >= xb.size() ? xb : xa;
      const std::vector<double> &py = xa.size() >= xb.size() ? yb : ya;
      if (px.size() < 3) continue;
      Fit FA = fitCircle(fx, fy);
      if (!FA.ok || FA.rms > 0.4) continue;
      double s = 0;
      for (size_t i = 0; i < px.size(); ++i)
        s += std::hypot(px[i] - FA.a, py[i] - FA.b) - FA.R;
      // signed mean radial residual of the minority side vs majority-side
      // circle: for near-radial tracks the field jump appears here (rphi
      // displacement ~ radial residual of the fitted arc at these radii)
      J.push_back(s / px.size());
    }
  }
  std::vector<double> aJ;
  for (double j : J) aJ.push_back(std::fabs(j));
  printf("CMC %s: full-crossers %ld | side-crossers %ld (%.1f%%)\n", tag, nfull, ncross,
         nfull ? 100. * ncross / nfull : 0);
  printf("CMC %s CMJUMP: n %zu | med|J| %.3f cm | width(q68) %.3f cm\n", tag, J.size(),
         med(aJ), q68(J));
  printf("CMC %s HALFARC: n %zu dk-width %.4f | n05 %zu dk05-width %.4f\n", tag,
         DK.size(), q68(DK), DK05.size(), q68(DK05));
}
