// ms_barscan.C — statistics and sensitivity of the MNF::fitBar cuts
// (>=12 points, radial span >=15 cm, 45<=R_fit<2e4 cm), asked 2026-08-20:
// per-cut distributions, cut-flow and N-1 marginals, and the stability of
// the headline medians under varying each cut, on BOTH sides (real
// ntp_clus_trk seeds, V6 laser veto; sim island91 V6 truth groups).
#include <TFile.h>
#include <TTree.h>
#include <cstdio>
#include <cstdarg>
#include <cmath>
#include <map>
#include <vector>
#include <algorithm>

namespace MBS
{
struct Fit { double a = 0, b = 0, R = 0, rms = 0; int n = 0; bool ok = false; };
Fit fitCircle(const std::vector<double> &X, const std::vector<double> &Y)
{
  Fit F; F.n = (int) X.size();
  if (F.n < 5) return F;
  double Sx = 0, Sy = 0, Sxx = 0, Syy = 0, Sxy = 0, Sxz = 0, Syz = 0, Sz = 0;
  for (size_t i = 0; i < X.size(); ++i)
  {
    double x = X[i], y = Y[i], z = x * x + y * y;
    Sx += x; Sy += y; Sxx += x * x; Syy += y * y; Sxy += x * y;
    Sxz += x * z; Syz += y * z; Sz += z;
  }
  double n = F.n;
  double det = Sxx * (Syy * n - Sy * Sy) - Sxy * (Sxy * n - Sy * Sx) + Sx * (Sxy * Sy - Syy * Sx);
  if (std::fabs(det) < 1e-9) return F;
  double A = (Sxz * (Syy * n - Sy * Sy) - Sxy * (Syz * n - Sy * Sz) + Sx * (Syz * Sy - Syy * Sz)) / det;
  double B = (Sxx * (Syz * n - Sy * Sz) - Sxz * (Sxy * n - Sy * Sx) + Sx * (Sxy * Sz - Syz * Sx)) / det;
  double C = (Sxx * (Syy * Sz - Syz * Sy) - Sxy * (Sxy * Sz - Syz * Sx) + Sxz * (Sxy * Sy - Syy * Sx)) / det;
  F.a = A / 2; F.b = B / 2;
  double r2 = C + F.a * F.a + F.b * F.b;
  if (r2 <= 0) return F;
  F.R = std::sqrt(r2);
  for (int it = 0; it < 6; ++it)
  {
    double M[3][3] = {{0}}, v[3] = {0};
    for (size_t i = 0; i < X.size(); ++i)
    {
      double dx = X[i] - F.a, dy = Y[i] - F.b, rho = std::hypot(dx, dy);
      if (rho < 1e-9) continue;
      double res = rho - F.R, J[3] = {-dx / rho, -dy / rho, -1.};
      for (int p = 0; p < 3; ++p)
      {
        v[p] -= J[p] * res;
        for (int q = 0; q < 3; ++q) M[p][q] += J[p] * J[q];
      }
    }
    double d = M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
             - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
             + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
    if (std::fabs(d) < 1e-12) break;
    double d0 = (v[0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
               - M[0][1] * (v[1] * M[2][2] - M[1][2] * v[2])
               + M[0][2] * (v[1] * M[2][1] - M[1][1] * v[2])) / d;
    double d1 = (M[0][0] * (v[1] * M[2][2] - M[1][2] * v[2])
               - v[0] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
               + M[0][2] * (M[1][0] * v[2] - v[1] * M[2][0])) / d;
    double d2 = (M[0][0] * (M[1][1] * v[2] - v[1] * M[2][1])
               - M[0][1] * (M[1][0] * v[2] - v[1] * M[2][0])
               + v[0] * (M[1][0] * M[2][1] - M[1][1] * M[2][0])) / d;
    F.a += d0; F.b += d1; F.R += d2;
  }
  double s2 = 0;
  for (size_t i = 0; i < X.size(); ++i)
  {
    double res = std::hypot(X[i] - F.a, Y[i] - F.b) - F.R;
    s2 += res * res;
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
double pct(std::vector<double> v, double p)
{
  if (v.empty()) return 0;
  std::sort(v.begin(), v.end());
  return v[(size_t) (p * (v.size() - 1))];
}
struct Rec { int n; double span, R, rms; bool okfit; };
}  // namespace MBS

void ms_barscan(const char *i91 = "island91_frames_production_v6.root",
                const char *realf = "../clusters_seeds_island_79507-0.root_ntuplizer.root",
                const char *ver = "v6")
{
  using namespace MBS;
  std::vector<Rec> rec[2];                            // [0]=real [1]=sim
  {
    TFile *f = TFile::Open(realf);
    if (!f || f->IsZombie()) { printf("missing %s\n", realf); return; }
    TTree *t = (TTree *) f->Get("ntp_clus_trk");
    float ev, sid, lay, x, y;
    t->SetBranchStatus("*", 0);
    for (auto b : {"event", "seedID", "layer", "x", "y"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("event", &ev);
    t->SetBranchAddress("seedID", &sid);
    t->SetBranchAddress("layer", &lay);
    t->SetBranchAddress("x", &x);
    t->SetBranchAddress("y", &y);
    std::map<std::pair<int, int>, std::vector<std::array<double, 2>>> g;
    for (Long64_t i = 0; i < t->GetEntries(); ++i)
    {
      t->GetEntry(i);
      if ((int) ev == 44) continue;                   // V6 laser veto (no-op here, convention)
      if (lay < 7 || lay > 54) continue;
      g[{(int) ev, (int) sid}].push_back({(double) x, (double) y});
    }
    f->Close();
    for (auto &kv : g)
    {
      std::vector<double> X, Y;
      double rlo = 1e9, rhi = 0;
      for (auto &p : kv.second)
      {
        X.push_back(p[0]); Y.push_back(p[1]);
        double r = std::hypot(p[0], p[1]);
        rlo = std::min(rlo, r); rhi = std::max(rhi, r);
      }
      Fit F = fitCircle(X, Y);
      rec[0].push_back({(int) X.size(), rhi - rlo, F.R, F.rms * 1e4, F.ok});
    }
  }
  {
    TFile *f = TFile::Open(i91);
    if (!f || f->IsZombie()) { printf("missing %s\n", i91); return; }
    TTree *c = (TTree *) f->Get("ntp_cluster");
    TTree *u = (TTree *) f->Get("ntp_truth");
    float ev, lay, x, y, tid;
    c->SetBranchStatus("*", 0);
    for (auto b : {"event", "layer", "x", "y"}) c->SetBranchStatus(b, 1);
    c->SetBranchAddress("event", &ev);
    c->SetBranchAddress("layer", &lay);
    c->SetBranchAddress("x", &x);
    c->SetBranchAddress("y", &y);
    u->SetBranchStatus("*", 0);
    u->SetBranchStatus("gtrackID", 1);
    u->SetBranchAddress("gtrackID", &tid);
    std::map<std::pair<int, int>, std::vector<std::array<double, 2>>> g;
    for (Long64_t i = 0; i < c->GetEntries(); ++i)
    {
      c->GetEntry(i); u->GetEntry(i);
      if (lay < 7 || lay > 54 || tid <= 0) continue;
      g[{(int) ev, (int) tid}].push_back({(double) x, (double) y});
    }
    f->Close();
    for (auto &kv : g)
    {
      std::vector<double> X, Y;
      double rlo = 1e9, rhi = 0;
      for (auto &p : kv.second)
      {
        X.push_back(p[0]); Y.push_back(p[1]);
        double r = std::hypot(p[0], p[1]);
        rlo = std::min(rlo, r); rhi = std::max(rhi, r);
      }
      Fit F = fitCircle(X, Y);
      rec[1].push_back({(int) X.size(), rhi - rlo, F.R, F.rms * 1e4, F.ok});
    }
  }
  FILE *fo = fopen(Form("ms_barscan_%s.txt", ver), "w");
  auto P = [&](const char *fmt, ...) {
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    va_start(ap, fmt); vfprintf(fo, fmt, ap); va_end(ap);
  };
  const char *sn[2] = {"real seeds ", "sim truthgr"};
  P("ms_barscan %s — statistics + sensitivity of the fitBar cuts (n>=12, span>=15, 45<=R<2e4)\n", ver);
  P("R=45 cm at 1.4 T <-> pT = 0.189 GeV (looper/track class boundary 0.164 GeV <-> R = 39.1 cm)\n");
  for (int s = 0; s < 2; ++s)
  {
    std::vector<double> nn, sp;
    for (auto &r : rec[s]) { nn.push_back(r.n); sp.push_back(r.span); }
    P("%s: %zu raw groups | n p10/med %.0f/%.0f | span p10/med %.1f/%.1f cm\n",
      sn[s], rec[s].size(), pct(nn, 0.10), med(nn), pct(sp, 0.10), med(sp));
    long f_n = 0, f_sp = 0, f_fit = 0, f_R = 0, pass = 0;
    long m_n = 0, m_sp = 0, m_R = 0, base_n = 0, base_sp = 0, base_R = 0;
    for (auto &r : rec[s])
    {
      bool cn = r.n >= 12, cs = r.span >= 15, cf = r.okfit, cR = r.okfit && r.R >= 45 && r.R < 2e4;
      if (!cn) f_n++;
      else if (!cs) f_sp++;
      else if (!cf) f_fit++;
      else if (!cR) f_R++;
      else pass++;
      if (cs && cR) { base_n++; if (!cn) m_n++; }
      if (cn && cR) { base_sp++; if (!cs) m_sp++; }
      if (cn && cs && cf) { base_R++; if (!cR) m_R++; }
    }
    P("  cut-flow: fail n %.1f%% -> fail span %.1f%% -> unfittable %.2f%% -> fail R-window %.1f%% -> PASS %.1f%%\n",
      100. * f_n / rec[s].size(), 100. * f_sp / rec[s].size(), 100. * f_fit / rec[s].size(),
      100. * f_R / rec[s].size(), 100. * pass / rec[s].size());
    P("  N-1 marginals: n<12 %.1f%% | span<15 %.1f%% | R outside %.1f%%\n",
      base_n ? 100. * m_n / base_n : 0, base_sp ? 100. * m_sp / base_sp : 0,
      base_R ? 100. * m_R / base_R : 0);
  }
  auto scan = [&](int nmin, double spmin, double rmin) {
    double m[2][2];                                   // [side][all/w05]
    for (int s = 0; s < 2; ++s)
    {
      std::vector<double> a, w;
      for (auto &r : rec[s])
      {
        if (r.n < nmin || r.span < spmin || !r.okfit || r.R < rmin || r.R >= 2e4) continue;
        a.push_back(r.rms);
        if (r.R > 101 && r.R < 137) w.push_back(r.rms);
      }
      m[s][0] = med(a); m[s][1] = med(w);
    }
    P("  n>=%2d span>=%2.0f R>=%2.0f : RMSmed all %4.0f / %4.0f um (ratio %.3f) | 0.5GeV-win %4.0f / %4.0f (ratio %.3f)\n",
      nmin, spmin, rmin, m[0][0], m[1][0], m[1][0] > 0 ? m[0][0] / m[1][0] : 0,
      m[0][1], m[1][1], m[1][1] > 0 ? m[0][1] / m[1][1] : 0);
  };
  P("sensitivity (real / sim medians, real:sim ratio) — NOMINAL first:\n");
  scan(12, 15, 45);
  for (int n : {8, 16, 20}) scan(n, 15, 45);
  for (double s : {10., 20., 25.}) scan(12, s, 45);
  for (double r : {40., 50., 60.}) scan(12, 15, r);
  {
    long nb = 0, tot = 0;
    for (int s = 0; s < 2; ++s)
    {
      nb = tot = 0;
      for (auto &r : rec[s])
      {
        if (r.n < 12 || r.span < 15 || !r.okfit || r.R < 45 || r.R >= 2e4) continue;
        tot++;
        if (r.R < 60) nb++;
      }
      P("%s: accepted with R in [45,60) (near-boundary): %.1f%%\n", sn[s], tot ? 100. * nb / tot : 0);
    }
  }
  fclose(fo);
  printf("wrote ms_barscan_%s.txt\n", ver);
}
