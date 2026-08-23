// missed_tracks.C — exhaustive circle search over ALL clusters (user,
// 2026-07-31): find track-class circles the classical tracker missed.
// Motivation: ~95% of real clusters are on no ntp_clus_trk track; are they
// all loopers/noise, or do they hide real tracks?
// Method (identical for real and sim):
//   1. conformal map (u,v) = (x,y)/r^2: circles through ~the beamline become
//      ~straight lines (d0 up to a few cm tolerated by wide roads);
//   2. Hough accumulator over line (theta, rho) proposes candidates
//      EXHAUSTIVELY; accepted tracks are subtracted from the accumulator and
//      the search continues until no candidate is left;
//   3. every candidate must pass the ESTABLISHED circularity bar (same as
//      ms_realcheck): iterative 3 mm outlier cleaning, >=10 survivors,
//      layer span >= 15, RMS <= 0.20 cm, R_fit >= 35 cm (track-class
//      curvature; loopers deliberately out of scope).
// Runs:
//   REAL-UNUSED: clusters NOT on any (event,seedID) track  -> missed tracks;
//   REAL-FULL  : all clusters -> re-find rate of the tracker's own tracks
//                (finder-efficiency cross-check on known tracks);
//   SIM        : island91 v5.3 clusters -> absolute finder efficiency and
//                purity vs truth, and the truth count of findable tracks.
// Per found track the median cluster tbin classifies IN-TIME (band 60-360)
// vs OUT-OF-TIME — the expected home of tracker-missed but genuine tracks.
// Output: ../sim_validation_plots/missed_tracks_<ver>.png + missed_tracks_<ver>.txt
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TEllipse.h>
#include <TGraph.h>
#include <cstdio>
#include <cstdarg>
#include <cmath>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <numeric>

namespace MTK
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

struct Cl { float x, y, lay, tb; int id; };       // id: sim gtrackID (-1 real)
struct Trk { std::vector<int> idx; Fit F; double medtb = 0; };

// the established acceptance bar (same circularity criteria as ms_realcheck)
struct CandRec
{
  int n, nlay, span, maxgap;
  double R, rms;
  bool acc;
  std::vector<int> idx;
};
std::vector<CandRec> *gCand = nullptr;   // armed by mt_rocscan around hunt()

bool acceptTrack(const std::vector<Cl> &C, std::vector<int> &idx, Fit &F,
                 double rescut = 0.30)
{
  for (int it = 0; it < 6; ++it)
  {
    if ((int) idx.size() < 12) return false;
    std::vector<double> X, Y;
    for (int i : idx) { X.push_back(C[i].x); Y.push_back(C[i].y); }
    F = fitCircle(X, Y);
    if (!F.ok) return false;
    std::vector<int> keep;
    for (int i : idx)
    {
      double res = std::hypot(C[i].x - F.a, C[i].y - F.b) - F.R;
      if (std::fabs(res) <= rescut) keep.push_back(i);
    }
    if (keep.size() == idx.size()) break;
    idx = keep;
  }
  if ((int) idx.size() < 12) return false;
  int lmin = 99, lmax = 0;
  std::set<int> lays;
  for (int i : idx)
  {
    lmin = std::min(lmin, (int) C[i].lay);
    lmax = std::max(lmax, (int) C[i].lay);
    lays.insert((int) C[i].lay);
  }
  // contiguity: real tracks have no big radial holes; stitched ghosts do
  std::vector<int> lv(lays.begin(), lays.end());
  int maxgap = 0;
  for (size_t i = 1; i < lv.size(); ++i)
    maxgap = std::max(maxgap, lv[i] - lv[i - 1]);
  bool pass = F.ok && F.rms <= 0.20 && lmax - lmin >= 15 && (int) lays.size() >= 13 &&
              F.R >= 45 && F.R < 2e4 && maxgap <= 6;
  // ROC-scan collector (2026-08-20): when armed, record every candidate that
  // survives the FIXED cleaning stage, with its gate scalars and the nominal
  // verdict. Claiming behavior is untouched (return value unchanged), so the
  // nominal candidate stream stays bit-exact; the offline cut scan re-applies
  // gates to these records (exact for tighter-than-nominal claiming, second-
  // order approximate for looser).
  if (gCand)
    gCand->push_back({(int) idx.size(), (int) lays.size(), lmax - lmin, maxgap,
                      F.R, F.rms, pass, idx});
  return pass;
}

// exhaustive circle search on one event's cluster list.
// A single global Hough saturates (29k clusters x 360 theta stamps / 144k
// bins -> mean ~72 per bin, drowning 10-cluster peaks), so the search is
// SECTORIZED: overlapping 30-degree azimuthal wedges, each with a local
// (theta, rho) Hough in conformal space (background ~6/bin, peaks of >=10
// stand out). Roads are collected event-wide; validation = acceptTrack.
// Two passes over all sectors; accepted clusters are flagged used.
std::vector<Trk> hunt(const std::vector<Cl> &C, double cohband = 8, double cohslope = 8)
{
  const int NSEC = 36, NT = 80, NRHO = 400;
  const double RHOMAX = 0.04, ROAD = 2.0e-4, HALFW = 15. * M_PI / 180.;
  const double ROADPHI = 40. * M_PI / 180.;   // roads are azimuthally LOCAL
  int n = (int) C.size();
  std::vector<double> u(n), v(n), phi(n);
  std::vector<char> used(n, 0);
  for (int i = 0; i < n; ++i)
  {
    double r2 = (double) C[i].x * C[i].x + (double) C[i].y * C[i].y;
    u[i] = C[i].x / r2;
    v[i] = C[i].y / r2;
    phi[i] = std::atan2(C[i].y, C[i].x);
  }
  auto wrap = [](double d) {
    while (d > M_PI) d -= 2 * M_PI;
    while (d < -M_PI) d += 2 * M_PI;
    return d;
  };
  std::vector<Trk> out;
  std::vector<int> acc(NT * NRHO);
  std::vector<int> sec;
  for (int pass = 0; pass < 2; ++pass)
  {
    for (int s = 0; s < NSEC; ++s)
    {
      double phic = -M_PI + (s + 0.5) * 2 * M_PI / NSEC;
      double t0 = phic + M_PI / 2 - 40. * M_PI / 180.;   // local theta range: sector normal +-40 deg
      sec.clear();
      for (int i = 0; i < n; ++i)
        if (!used[i] && std::fabs(wrap(phi[i] - phic)) < HALFW) sec.push_back(i);
      if ((int) sec.size() < 10) continue;
      std::fill(acc.begin(), acc.end(), 0);
      for (int i : sec)
        for (int t = 0; t < NT; ++t)
        {
          double th = t0 + t * (M_PI / 180.);
          int b = (int) ((u[i] * std::cos(th) + v[i] * std::sin(th) + RHOMAX) / (2 * RHOMAX) * NRHO);
          if (b >= 0 && b < NRHO) acc[t * NRHO + b]++;
        }
      for (int tries = 0; tries < 60; ++tries)
      {
        int best = -1, bc = 8;                     // require >= 9 in a bin
        for (int b = 0; b < NT * NRHO; ++b)
          if (acc[b] > bc) { bc = acc[b]; best = b; }
        if (best < 0) break;
        int tb = best / NRHO, rb = best % NRHO;
        double th = t0 + tb * (M_PI / 180.);
        double rho = -RHOMAX + (rb + 0.5) * (2 * RHOMAX) / NRHO;
        double cth = std::cos(th), sth = std::sin(th);
        // road: conformal-line proximity AND azimuthal locality
        std::vector<int> road;
        for (int i = 0; i < n; ++i)
        {
          if (used[i] || std::fabs(wrap(phi[i] - phic)) > ROADPHI) continue;
          if (std::fabs(u[i] * cth + v[i] * sth - rho) > ROAD) continue;
          road.push_back(i);
        }
        // DRIFT COHERENCE: a real track is also a straight line in
        // tbin-vs-layer (z drifts smoothly along the crossing); random road
        // clusters spread over ~970 tbins and cannot hold a +-12-bin linear
        // consensus. RANSAC over cluster pairs; without this the x-y
        // projection is combinatorially saturated at this pileup (measured:
        // 0.3% finder purity in sim without it).
        std::vector<int> cand;
        if ((int) road.size() >= 12)
        {
          unsigned lcg = 12345u + (unsigned) (tb * 131 + rb * 7 + tries);
          auto rnd = [&]() { lcg = lcg * 1664525u + 1013904223u; return lcg; };
          size_t bestn = 0;
          double bA = 0, bB = 0;
          int npair = std::min((size_t) 40, road.size() * (road.size() - 1) / 2);
          for (int q = 0; q < npair; ++q)
          {
            int i1 = road[rnd() % road.size()], i2 = road[rnd() % road.size()];
            if (i1 == i2 || std::fabs(C[i1].lay - C[i2].lay) < 6) continue;
            double A = (C[i2].tb - C[i1].tb) / (C[i2].lay - C[i1].lay);
            if (std::fabs(A) > cohslope) continue;          // coherence slope sanity
            double B = C[i1].tb - A * C[i1].lay;
            size_t nin = 0;
            for (int i : road)
              if (std::fabs(C[i].tb - (A * C[i].lay + B)) < cohband) nin++;
            if (nin > bestn) { bestn = nin; bA = A; bB = B; }
          }
          if ((int) bestn >= 13)
          {
            // keep inliers, then best-1 cluster per layer (by conformal distance)
            std::map<int, std::pair<double, int>> bylay;
            for (int i : road)
            {
              if (std::fabs(C[i].tb - (bA * C[i].lay + bB)) >= cohband) continue;
              double d = std::fabs(u[i] * cth + v[i] * sth - rho);
              int L = (int) C[i].lay;
              auto it = bylay.find(L);
              if (it == bylay.end() || d < it->second.first) bylay[L] = {d, i};
            }
            for (auto &q : bylay) cand.push_back(q.second.second);
          }
        }
        bool ok = false;
        Fit F;
        if ((int) cand.size() >= 12)
        {
          std::vector<int> idx = cand;
          if (acceptTrack(C, idx, F))
          {
            Trk T;
            T.idx = idx;
            T.F = F;
            std::vector<double> tbs;
            for (int i : idx) { used[i] = 1; tbs.push_back(C[i].tb); }
            std::sort(tbs.begin(), tbs.end());
            T.medtb = tbs[tbs.size() / 2];
            out.push_back(T);
            ok = true;
          }
        }
        if (ok)
        {
          // remove the accepted clusters from this sector's accumulator
          for (int i : out.back().idx)
          {
            bool insec = std::fabs(wrap(phi[i] - phic)) < HALFW;
            if (!insec) continue;
            for (int t = 0; t < NT; ++t)
            {
              double th2 = t0 + t * (M_PI / 180.);
              int b = (int) ((u[i] * std::cos(th2) + v[i] * std::sin(th2) + RHOMAX) / (2 * RHOMAX) * NRHO);
              if (b >= 0 && b < NRHO) acc[t * NRHO + b]--;
            }
          }
        }
        else
        {
          for (int dt = -1; dt <= 1; ++dt)
            for (int dr = -1; dr <= 1; ++dr)
            {
              int t2 = tb + dt, r2 = rb + dr;
              if (t2 >= 0 && t2 < NT && r2 >= 0 && r2 < NRHO) acc[t2 * NRHO + r2] = 0;
            }
        }
      }
    }
  }
  return out;
}
}  // namespace MTK

void missed_tracks(const char *realf = "../clusters_seeds_island_79507-0.root_ntuplizer.root",
                   const char *i91 = "island91_frames_production_v53.root",
                   const char *ver = "v53f", const char *vtag = "v5.3",
                   int nsimev = 50)
{
  using namespace MTK;

  // ---------- REAL: clusters + used-flag from ntp_clus_trk ----------
  std::map<int, std::vector<Cl>> rev;             // event -> clusters
  std::map<int, std::vector<char>> ronrk;         // event -> on-track flag
  {
    std::set<long long> trkkey;
    TFile *f = TFile::Open(realf);
    TTree *t = (TTree *) f->Get("ntp_clus_trk");
    float ev, lay, x, y;
    t->SetBranchStatus("*", 0);
    for (auto b : {"event", "layer", "x", "y"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("event", &ev);
    t->SetBranchAddress("layer", &lay);
    t->SetBranchAddress("x", &x);
    t->SetBranchAddress("y", &y);
    long ntrkrows = 0;
    for (Long64_t i = 0; i < t->GetEntries(); ++i)
    {
      t->GetEntry(i);
      if ((int) ev == 44) continue;   // V6 laser veto (canon.h)
      if (lay < 7 || lay > 54) continue;
      trkkey.insert(((long long) ev << 40) ^ ((long long) llround(x * 1e3) << 20) ^ (long long) llround(y * 1e3));
      ntrkrows++;
    }
    TTree *c = (TTree *) f->Get("ntp_cluster");
    float tb;
    c->SetBranchStatus("*", 0);
    for (auto b : {"event", "layer", "x", "y", "tbin"}) c->SetBranchStatus(b, 1);
    c->SetBranchAddress("event", &ev);
    c->SetBranchAddress("layer", &lay);
    c->SetBranchAddress("x", &x);
    c->SetBranchAddress("y", &y);
    c->SetBranchAddress("tbin", &tb);
    long nmatch = 0, nc = 0;
    for (Long64_t i = 0; i < c->GetEntries(); ++i)
    {
      c->GetEntry(i);
      if ((int) ev == 44) continue;   // V6 laser veto (canon.h)
      if (lay < 7 || lay > 54) continue;
      bool on = trkkey.count(((long long) ev << 40) ^ ((long long) llround(x * 1e3) << 20) ^ (long long) llround(y * 1e3)) > 0;
      rev[(int) ev].push_back({x, y, lay, tb, -1});
      ronrk[(int) ev].push_back(on ? 1 : 0);
      if (on) nmatch++;
      nc++;
    }
    f->Close();
    printf("real: %ld TPC clusters, %ld flagged on-track (clus_trk rows %ld -> match check)\n",
           nc, nmatch, ntrkrows);
  }

  // ---------- run finder: REAL-UNUSED and REAL-FULL ----------
  long nMiss = 0, nMissOOT = 0, nFull = 0, nKnown = 0, nKnownRefound = 0;
  std::vector<double> missR, missN, missTb;
  int showev = -1;
  std::vector<Trk> showTrks;
  std::vector<Cl> showCl;
  std::vector<char> showOn;
  for (auto &kv : rev)
  {
    std::vector<Cl> &C = kv.second;
    std::vector<char> &on = ronrk[kv.first];
    // unused set
    std::vector<Cl> Cu;
    for (size_t i = 0; i < C.size(); ++i)
      if (!on[i]) Cu.push_back(C[i]);
    std::vector<Trk> miss = hunt(Cu);
    for (auto &T : miss)
    {
      nMiss++;
      missR.push_back(T.F.R);
      missN.push_back((double) T.idx.size());
      missTb.push_back(T.medtb);
      if (T.medtb < 60 || T.medtb > 360) nMissOOT++;
    }
    if (showev < 0 && miss.size() >= 8) { showev = kv.first; showTrks = miss; showCl = C; showOn = on; }
    // full-set run: re-find rate of the tracker's own tracks
    std::vector<Trk> full = hunt(C);
    nFull += (long) full.size();
    // known tracks of this event = connected clusters flagged on-track,
    // grouped is unnecessary: count a known track re-found if >=50% of some
    // found track's clusters are on-track flagged (proxy at event level)
    for (auto &T : full)
    {
      int non = 0;
      for (int i : T.idx) if (on[i]) non++;
      if (non >= (int) T.idx.size() / 2) nKnownRefound++;
    }
  }
  // count known full-crosser tracker tracks (from ms_realcheck acceptance: ~ use groups >= 10 on-track clusters per seed — proxy: clus_trk groups already audited; use audited number 6038 accepted)
  nKnown = 6038;                                   // 5286 clean + 752 rescued (ms_realcheck)

  // ---------- SIM: finder efficiency + purity vs truth ----------
  long sTruth = 0, sTruthMatched = 0, sFound = 0, sPure = 0;
  {
    TFile *f = TFile::Open(i91);
    TTree *c = (TTree *) f->Get("ntp_cluster");
    TTree *u = (TTree *) f->Get("ntp_truth");
    float ev, lay, x, y, tb, tid, cls, ntrks;
    c->SetBranchStatus("*", 0);
    for (auto b : {"event", "layer", "x", "y", "tbin"}) c->SetBranchStatus(b, 1);
    c->SetBranchAddress("event", &ev);
    c->SetBranchAddress("layer", &lay);
    c->SetBranchAddress("x", &x);
    c->SetBranchAddress("y", &y);
    c->SetBranchAddress("tbin", &tb);
    u->SetBranchStatus("*", 0);
    for (auto b : {"gtrackID", "cls", "ntrks"}) u->SetBranchStatus(b, 1);
    u->SetBranchAddress("gtrackID", &tid);
    u->SetBranchAddress("cls", &cls);
    u->SetBranchAddress("ntrks", &ntrks);
    std::map<int, std::vector<Cl>> sev;
    for (Long64_t i = 0; i < c->GetEntries(); ++i)
    {
      u->GetEntry(i);
      c->GetEntry(i);
      if ((int) ev >= nsimev) continue;
      int id = (cls == 0 && ntrks == 1) ? (int) tid : -1;
      sev[(int) ev].push_back({x, y, lay, tb, id});
    }
    f->Close();
    for (auto &kv : sev)
    {
      std::vector<Cl> &C = kv.second;
      // truth full-crossers (the findable denominator), same bar as finder
      std::map<int, std::vector<int>> tg;
      for (size_t i = 0; i < C.size(); ++i)
        if (C[i].id >= 0) tg[C[i].id].push_back((int) i);
      std::map<int, int> truthN;
      for (auto &g : tg)
      {
        if ((int) g.second.size() < 10) continue;
        std::vector<int> idx = g.second;
        Fit F;
        if (acceptTrack(C, idx, F)) { sTruth++; truthN[g.first] = (int) g.second.size(); }
      }
      std::vector<Trk> found = hunt(C);
      std::set<int> matched;
      for (auto &T : found)
      {
        sFound++;
        std::map<int, int> vote;
        for (int i : T.idx) if (C[i].id >= 0) vote[C[i].id]++;
        int mid = -1, mv = 0;
        for (auto &q : vote) if (q.second > mv) { mv = q.second; mid = q.first; }
        if (mid >= 0 && mv >= (int) (0.7 * T.idx.size())) sPure++;
        if (mid >= 0 && truthN.count(mid) && mv >= truthN[mid] / 2) matched.insert(mid);
      }
      sTruthMatched += (long) matched.size();
    }
  }

  // ---------- ledger ----------
  FILE *fo = fopen(Form("missed_tracks_%s.txt", ver), "w");
  auto P = [&](const char *fmt, ...) {
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    va_start(ap, fmt); vfprintf(fo, fmt, ap); va_end(ap);
  };
  auto med = [](std::vector<double> v) {
    if (v.empty()) return 0.;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
  };
  P("missed_tracks %s — exhaustive Hough+circle search for tracks the classical\n", ver);
  P("tracker missed. Bar: >=13 distinct layers after 3 mm cleaning, span >=15,\n");
  P("no layer gap >6, RMS <= 0.20 cm, R_fit >= 45 cm; roads azimuthally local\n");
  P("(+-40 deg), drift-coherent (tbin-vs-layer RANSAC +-8), best-1 per layer.\n");
  P("FINDER CALIBRATION (sim %s, %d events, truth-known):\n", vtag, nsimev);
  P("  truth findable full-crossers: %ld (%.1f/event)\n", sTruth, (double) sTruth / nsimev);
  P("  finder efficiency on them: %.1f%%   purity of found tracks: %.1f%% (%ld found)\n",
    sTruth ? 100. * sTruthMatched / sTruth : 0, sFound ? 100. * sPure / sFound : 0, sFound);
  P("REAL, full cluster set: finder returns %ld tracks/100 ev; %.1f%% of them are\n",
    nFull, nFull ? 100. * nKnownRefound / nFull : 0);
  P("  majority-on-track (re-finds of the tracker's %ld audited tracks).\n", nKnown);
  P("REAL, UNUSED clusters only (the actual question):\n");
  P("  MISSED track-class circles found: %ld (%.2f/event)\n", nMiss, nMiss / 100.);
  P("  of these OUT-OF-TIME (median tbin outside [60,360]): %ld (%.1f%%)\n",
    nMissOOT, nMiss ? 100. * nMissOOT / nMiss : 0);
  P("  medians: clusters/track %.0f, R_fit %.0f cm (pT ~%.2f GeV), tbin %.0f\n",
    med(missN), med(missR), med(missR) * 0.299792458 * 1.4 / 100, med(missTb));
  double effs = sTruth ? (double) sTruthMatched / sTruth : 0;
  double purs = sFound ? (double) sPure / sFound : 0;
  P("SIM-CALIBRATED estimate: true missed tracks/event ~ found x purity / eff\n");
  P("  = %.1f x %.2f / %.2f = %.0f per event (in-time share %.0f%%)\n",
    nMiss / 100., purs, effs, effs > 0 ? nMiss / 100. * purs / effs : 0,
    nMiss ? 100. * (nMiss - nMissOOT) / nMiss : 0);
  P("reading: the tracker saves ~60 tracks/event; the in-time missed subset is\n");
  P("recoverable inefficiency, the out-of-time subset is invisible to an\n");
  P("in-time tracker by design but is real charge that cluster-level ML must\n");
  P("classify. All rates carry the finder calibration above.\n");
  fclose(fo);

  // ---------- figure: showcase event + summary ----------
  gStyle->SetOptStat(0);
  TCanvas *cv = new TCanvas("cvm", Form("missed tracks %s", ver), 1500, 700);
  cv->Divide(2, 1);
  cv->cd(1);
  {
    TH1 *fr = gPad->DrawFrame(-80, -80, 80, 80,
                              Form("real event %d: unused clusters + missed tracks;x [cm];y [cm]", showev));
    fr->GetYaxis()->SetTitleOffset(1.25);
    TGraph *gu = new TGraph();
    for (size_t i = 0; i < showCl.size(); ++i)
      if (!showOn[i]) gu->SetPoint(gu->GetN(), showCl[i].x, showCl[i].y);
    gu->SetMarkerStyle(1);
    gu->SetMarkerColorAlpha(kGray + 1, 0.5);
    gu->Draw("P same");
    int cols[3] = {kBlue + 1, kRed + 1, kGreen + 2};
    // T.idx indexes the UNUSED list — rebuild it to draw member clusters
    std::vector<Cl> Cu;
    for (size_t i = 0; i < showCl.size(); ++i)
      if (!showOn[i]) Cu.push_back(showCl[i]);
    for (size_t k = 0; k < showTrks.size(); ++k)
    {
      Trk &T = showTrks[k];
      TGraph *gt = new TGraph();
      for (int i : T.idx) gt->SetPoint(gt->GetN(), Cu[i].x, Cu[i].y);
      gt->SetMarkerStyle(20);
      gt->SetMarkerSize(0.45);
      gt->SetMarkerColor(cols[k % 3]);
      gt->Draw("P same");
      TEllipse *el = new TEllipse(T.F.a, T.F.b, T.F.R, T.F.R);
      el->SetFillStyle(0);
      el->SetLineColorAlpha(cols[k % 3], 0.6);
      el->SetLineStyle(3);
      el->Draw();
    }
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.034);
    tx.DrawLatex(0.13, 0.86, Form("%zu missed tracks in this event", showTrks.size()));
    tx.DrawLatex(0.13, 0.80, "(first event with >= 8 finds)");
  }
  cv->cd(2);
  {
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.040);
    tx.DrawLatex(0.06, 0.90, "exhaustive circle search: summary");
    tx.SetTextSize(0.034);
    tx.DrawLatex(0.06, 0.80, Form("finder calibration on sim %s (truth-known):", vtag));
    tx.DrawLatex(0.10, 0.74, Form("efficiency %.0f%%, purity %.0f%% on %.1f findable/event",
                                  sTruth ? 100. * sTruthMatched / sTruth : 0,
                                  sFound ? 100. * sPure / sFound : 0, (double) sTruth / nsimev));
    tx.DrawLatex(0.06, 0.64, Form("real, full set: %.1f found/event; %.0f%% re-find the tracker's tracks",
                                  nFull / 100., nFull ? 100. * nKnownRefound / nFull : 0));
    tx.DrawLatex(0.06, 0.54, Form("real, UNUSED clusters: %.2f missed tracks/event", nMiss / 100.));
    tx.DrawLatex(0.10, 0.48, Form("%.0f%% out-of-time (median tbin outside [60,360])",
                                  nMiss ? 100. * nMissOOT / nMiss : 0));
    tx.DrawLatex(0.10, 0.42, Form("median R_{fit} %.0f cm  (p_{T} ~ %.2f GeV), %.0f clusters/track",
                                  med(missR), med(missR) * 0.299792458 * 1.4 / 100, med(missN)));
    tx.SetTextSize(0.030);
    tx.DrawLatex(0.06, 0.28, "in-time finds = recoverable tracker inefficiency;");
    tx.DrawLatex(0.06, 0.22, "out-of-time finds = genuine tracks an in-time tracker skips by design,");
    tx.DrawLatex(0.06, 0.16, "but real charge that cluster-level ML must classify.");
  }
  cv->SaveAs(Form("../sim_validation_plots/missed_tracks_%s.png", ver));
  printf("wrote ../sim_validation_plots/missed_tracks_%s.png + missed_tracks_%s.txt\n", ver, ver);
}

// ---------------------------------------------------------------------------
// mt_cluscmp — ITEM 2 (2026-08-06): algorithm-symmetric CLUSTER-level
// comparison. The SAME finder groups BOTH sides' full ntp_cluster trees
// (real: all TPC clusters; sim v5.4c: all island91 clusters, NO truth
// filter) — removing the tracker-vs-truth grouping asymmetry of ms_real.
// Real rates quoted on COMPLETE events (cluster-tbin p99.9 > 950) per the
// dual-reference convention; sim frames are always complete.
void mt_cluscmp(const char *realf = "../clusters_seeds_island_79507-0.root_ntuplizer.root",
                const char *i91 = "island91_frames_production_v6.root",
                const char *ver = "v6", int nsimev = 50)
{
  using namespace MTK;
  std::map<int, std::vector<Cl>> ev[2];
  std::map<int, bool> complete;
  {
    TFile *f = TFile::Open(realf);
    TTree *c = (TTree *) f->Get("ntp_cluster");
    float evn, lay, x, y, tb;
    c->SetBranchStatus("*", 0);
    for (auto b : {"event", "layer", "x", "y", "tbin"}) c->SetBranchStatus(b, 1);
    c->SetBranchAddress("event", &evn);
    c->SetBranchAddress("layer", &lay);
    c->SetBranchAddress("x", &x);
    c->SetBranchAddress("y", &y);
    c->SetBranchAddress("tbin", &tb);
    std::map<int, std::vector<int>> tbc;           // event -> 5-tbin counts
    for (Long64_t i = 0; i < c->GetEntries(); ++i)
    {
      c->GetEntry(i);
      if ((int) evn == 44) continue;   // V6 laser veto (canon.h)
      if (lay < 7 || lay > 54) continue;
      ev[0][(int) evn].push_back({x, y, lay, tb, -1});
      auto &h = tbc[(int) evn];
      if (h.empty()) h.assign(200, 0);
      int b = (int) (tb / 5);
      if (b >= 0 && b < 200) h[b]++;
    }
    f->Close();
    for (auto &kv : tbc)
    {
      long tot = 0;
      for (int q : kv.second) tot += q;
      long acc2 = 0;
      int endp = 0;
      for (int b = 0; b < 200; ++b)
      {
        acc2 += kv.second[b];
        if (acc2 >= (long) (0.999 * tot)) { endp = b * 5; break; }
      }
      complete[kv.first] = endp > 950;
    }
  }
  {
    TFile *f = TFile::Open(i91);
    TTree *c = (TTree *) f->Get("ntp_cluster");
    float evn, lay, x, y, tb;
    c->SetBranchStatus("*", 0);
    for (auto b : {"event", "layer", "x", "y", "tbin"}) c->SetBranchStatus(b, 1);
    c->SetBranchAddress("event", &evn);
    c->SetBranchAddress("layer", &lay);
    c->SetBranchAddress("x", &x);
    c->SetBranchAddress("y", &y);
    c->SetBranchAddress("tbin", &tb);
    for (Long64_t i = 0; i < c->GetEntries(); ++i)
    {
      c->GetEntry(i);
      if ((int) evn >= nsimev) continue;
      ev[1][(int) evn].push_back({x, y, lay, tb, -1});
    }
    f->Close();
  }
  const char *sn[2] = {"real", "sim"};
  std::vector<double> nfound[2], rmsv[2], rv[2];
  long nin[2] = {0, 0}, ntot[2] = {0, 0}, nevC = 0;
  TH1D *hn[2], *hr[2], *hR[2], *hT[2];
  for (int s = 0; s < 2; ++s)
  {
    hn[s] = new TH1D(Form("hcn%d", s), ";track-class circles found / event;events (unit area)", 30, 0, 600);
    hr[s] = new TH1D(Form("hcr%d", s), ";per-track circle-fit RMS [mm];tracks (unit area)", 50, 0, 2.5);
    hR[s] = new TH1D(Form("hcR%d", s), ";fitted radius [cm]  (p_{T} = 0.0042 R);circles / event / bin", 30, 45, 345);
    hT[s] = new TH1D(Form("hcT%d", s), ";median drift-time bin of the circle;circles / event / bin", 40, 0, 1000);
  }
  // per-track dump (v5.5 figure work): side, event, n, R, rms, medtb
  FILE *fdump = fopen(Form("ms_cluscmp_%s_tracks.txt", ver), "w");
  fprintf(fdump, "# side(0=real,1=sim) event ncl R_cm rms_cm medtb\n");
  for (int s = 0; s < 2; ++s)
  {
    for (auto &kv : ev[s])
    {
      if (s == 0 && !complete[kv.first]) continue;   // real: complete events only
      if (s == 0) nevC++;
      std::vector<Trk> t = hunt(kv.second);
      nfound[s].push_back((double) t.size());
      hn[s]->Fill(std::min((double) t.size(), 599.));
      for (auto &T : t)
      {
        ntot[s]++;
        if (T.medtb >= 60 && T.medtb <= 360) nin[s]++;
        rmsv[s].push_back(T.F.rms * 10);
        rv[s].push_back(T.F.R);
        hr[s]->Fill(std::min(T.F.rms * 10, 2.49));
        hR[s]->Fill(std::min(T.F.R, 344.9));
        hT[s]->Fill(T.medtb);
        fprintf(fdump, "%d %d %d %.2f %.4f %.1f\n", s, kv.first, (int) T.idx.size(), T.F.R, T.F.rms, T.medtb);
      }
    }
  }
  fclose(fdump);
  auto med = [](std::vector<double> v) {
    if (v.empty()) return 0.;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
  };
  FILE *fo = fopen(Form("ms_cluscmp_%s.txt", ver), "w");
  auto P = [&](const char *fmt, ...) {
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    va_start(ap, fmt); vfprintf(fo, fmt, ap); va_end(ap);
  };
  P("ms_cluscmp %s — SYMMETRIC cluster-level exhaustive circle comparison\n", ver);
  P("same finder + acceptance both sides; real = complete events (%ld), sim = %d frames\n",
    nevC, nsimev);
  for (int s = 0; s < 2; ++s)
    P("%s: found/event median %.0f (mean %.1f) | in-time frac %.2f | "
      "per-track RMS med %.0f um | R_fit med %.0f cm\n",
      sn[s], med(nfound[s]),
      nfound[s].empty() ? 0 : std::accumulate(nfound[s].begin(), nfound[s].end(), 0.) / nfound[s].size(),
      ntot[s] ? (double) nin[s] / ntot[s] : 0, med(rmsv[s]) * 1000, med(rv[s]));
  P("data/MC: rate %.2f | RMS %.2f\n",
    med(nfound[1]) > 0 ? med(nfound[0]) / med(nfound[1]) : 0,
    med(rmsv[1]) > 0 ? med(rmsv[0]) / med(rmsv[1]) : 0);
  fclose(fo);
  gStyle->SetOptStat(0);
  TCanvas *cv = new TCanvas("cvcc", Form("cluscmp %s", ver), 1500, 620);
  cv->Divide(2, 1);
  for (int p = 0; p < 2; ++p)
  {
    cv->cd(p + 1);
    TH1D **h = p == 0 ? hn : hr;
    for (int s = 0; s < 2; ++s)
      if (h[s]->Integral() > 0) h[s]->Scale(1. / h[s]->Integral());
    h[0]->SetLineColor(kBlack); h[0]->SetLineWidth(2);
    h[1]->SetLineColor(kBlue + 1); h[1]->SetLineWidth(2); h[1]->SetLineStyle(2);
    h[0]->SetTitle(p == 0 ? "exhaustive finder: tracks per event (symmetric grouping)"
                          : "per-track circle-fit RMS (symmetric grouping)");
    h[0]->SetMaximum(1.4 * std::max(h[0]->GetMaximum(), h[1]->GetMaximum()));
    h[0]->Draw("hist");
    h[1]->Draw("hist same");
    TLegend *lg = new TLegend(0.50, 0.72, 0.89, 0.88);
    lg->SetBorderSize(0);
    lg->AddEntry(h[0], p == 0 ? Form("real complete, med %.0f/ev", med(nfound[0]))
                              : Form("real, med %.0f #mum", med(rmsv[0]) * 1000), "l");
    lg->AddEntry(h[1], p == 0 ? Form("sim %s, med %.0f/ev", ver, med(nfound[1]))
                              : Form("sim %s, med %.0f #mum", ver, med(rmsv[1]) * 1000), "l");
    lg->Draw();
  }
  {
    // per-event normalized radius + drift-time spectra of found circles (report figure input)
    double nR = nevC > 0 ? nevC : 1, nS = nsimev > 0 ? nsimev : 1;
    hR[0]->Scale(1. / nR); hR[1]->Scale(1. / nS); hT[0]->Scale(1. / nR); hT[1]->Scale(1. / nS);
    TFile fh(Form("ms_cluscmp_%s_hists.root", ver), "RECREATE");
    for (int s = 0; s < 2; ++s) { hn[s]->Write(); hr[s]->Write(); hR[s]->Write(); hT[s]->Write(); }
    fh.Close();
  }
  cv->SaveAs(Form("../sim_validation_plots/ms_cluscmp_%s.png", ver));
  printf("wrote ms_cluscmp_%s outputs\n", ver);
}

// ---------------------------------------------------------------------------
// mt_g4scan — ITEM 3a (2026-08-06): exhaustive circle search on sim TRUTH
// hits (ntp_g4hit), per single pp collision. Pseudo-layer = int(r) [cm]
// (truth steps carry no pad row); drift-coherence coordinate = gz
// (band 5 cm, slope 6 cm per pseudo-layer). Reports the from-first-
// principles findable track-class rate per collision, against the direct
// truth-group count under the same acceptance.
void mt_g4scan(const char *g4pat = "../P5/PP_g4hit_%d.root", int nfiles = 1,
               const char *ver = "v6")
{
  using namespace MTK;
  long ncoll = 0, nfound = 0, ntruth = 0, nrec = 0, ndupe = 0, nsub = 0, nghost = 0, nrms05 = 0;
  std::vector<double> rv, nv, pv;
  for (int fi = 0; fi < nfiles; ++fi)
  {
    TFile *f = TFile::Open(Form(g4pat, fi));
    if (!f || f->IsZombie()) continue;
    TTree *t = (TTree *) f->Get("ntp_g4hit");
    float ev, gx, gy, gz, tid;
    t->SetBranchStatus("*", 0);
    for (auto b : {"event", "gx", "gy", "gz", "gtrackID"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("event", &ev);
    t->SetBranchAddress("gx", &gx);
    t->SetBranchAddress("gy", &gy);
    t->SetBranchAddress("gz", &gz);
    t->SetBranchAddress("gtrackID", &tid);
    std::map<int, std::vector<Cl>> evs;
    for (Long64_t i = 0; i < t->GetEntries(); ++i)
    {
      t->GetEntry(i);
      double r = std::hypot(gx, gy);
      if (r < 20 || r > 78) continue;
      evs[(int) ev].push_back({gx, gy, (float) (int) r, gz, (int) tid});
    }
    f->Close();
    for (auto &kv : evs)
    {
      std::vector<Cl> &C = kv.second;
      ncoll++;
      // direct truth count under the same acceptance
      std::map<int, std::vector<int>> tg;
      for (size_t i = 0; i < C.size(); ++i) tg[C[i].id].push_back((int) i);
      std::set<int> acc;
      for (auto &g : tg)
      {
        if ((int) g.second.size() < 12) continue;
        std::vector<int> idx = g.second;
        Fit F;
        if (acceptTrack(C, idx, F)) { ntruth++; acc.insert(g.first); }
      }
      std::vector<Trk> found = hunt(C, 5, 6);
      nfound += (long) found.size();
      nv.push_back((double) found.size());
      std::set<int> claimed;
      for (auto &T : found)
      {
        rv.push_back(T.F.R);
        if (T.F.rms <= 0.05) nrms05++;
        std::map<int, int> cnt;
        for (int i : T.idx) cnt[C[i].id]++;
        int mtid = -1, mc = 0;
        for (auto &q : cnt)
          if (q.second > mc) { mc = q.second; mtid = q.first; }
        double pur = T.idx.empty() ? 0 : (double) mc / T.idx.size();
        pv.push_back(pur);
        if (pur < 0.6) { nghost++; continue; }            // mixed-parent chain
        if (!acc.count(mtid)) { nsub++; continue; }       // real particle below the acceptance bar
        if (claimed.count(mtid)) { ndupe++; continue; }   // second arc of an already-found particle
        claimed.insert(mtid);
        nrec++;
      }
    }
  }
  auto med = [](std::vector<double> v) {
    if (v.empty()) return 0.;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
  };
  FILE *fo = fopen(Form("ms_g4scan_%s.txt", ver), "w");
  auto P = [&](const char *fmt, ...) {
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    va_start(ap, fmt); vfprintf(fo, fmt, ap); va_end(ap);
  };
  P("ms_g4scan %s — exhaustive circle search on sim TRUTH hits (per collision)\n", ver);
  P("collisions %ld | truth findable %.3f/coll | finder found %.3f/coll\n",
    ncoll, (double) ntruth / std::max(1L, ncoll), (double) nfound / std::max(1L, ncoll));
  double fn = std::max(1L, nfound);
  P("composition of found: recovered-findable %.3f/coll (recall %.2f) | extra arc of same particle %.2f | "
    "real particle below bar %.2f | mixed-parent %.2f\n",
    (double) nrec / std::max(1L, ncoll), ntruth ? (double) nrec / ntruth : 0,
    ndupe / fn, nsub / fn, nghost / fn);
  P("single-parent purity median %.2f | fit RMS<=0.05 cm fraction %.2f\n", med(pv), nrms05 / fn);
  P("found R_fit median %.0f cm; found/collision median %.0f\n", med(rv), med(nv));
  fclose(fo);
  printf("wrote ms_g4scan_%s.txt\n", ver);
}

// ---------------------------------------------------------------------------
// mt_pixscan — ITEM 3b (2026-08-06): exhaustive circle search on ALL real
// ntp_hit PIXELS (pre-clustering level; canonical cut layer 7-54 && adc>0),
// on the first nev COMPLETE events (pixel-tbin p99.9 > 950). Same finder,
// same acceptance, tbin coherence. Pixel positions are pad centers
// (no charge weighting), so the per-point RMS is expected slightly above
// the cluster level.
void mt_pixscan(const char *realf = "../clusters_seeds_island_79507-0.root_ntuplizer.root",
                int nev = 25, const char *ver = "v6")
{
  using namespace MTK;
  TFile *f = TFile::Open(realf);
  TTree *t = (TTree *) f->Get("ntp_hit");
  float ev, lay, x, y, tb, adc;
  t->SetBranchStatus("*", 0);
  for (auto b : {"event", "layer", "x", "y", "tbin", "adc"}) t->SetBranchStatus(b, 1);
  t->SetBranchAddress("event", &ev);
  t->SetBranchAddress("layer", &lay);
  t->SetBranchAddress("x", &x);
  t->SetBranchAddress("y", &y);
  t->SetBranchAddress("tbin", &tb);
  t->SetBranchAddress("adc", &adc);
  // pass 1: endpoints
  std::map<int, std::vector<long>> tbc;
  for (Long64_t i = 0; i < t->GetEntries(); ++i)
  {
    t->GetEntry(i);
    if ((int) ev == 44) continue;   // V6 laser veto (canon.h)
    if (lay < 7 || lay > 54 || adc <= 0) continue;
    auto &h = tbc[(int) ev];
    if (h.empty()) h.assign(200, 0);
    int b = (int) (tb / 5);
    if (b >= 0 && b < 200) h[b]++;
  }
  std::set<int> chosen;
  for (auto &kv : tbc)
  {
    if ((int) chosen.size() >= nev) break;
    long tot = 0;
    for (long q : kv.second) tot += q;
    long acc2 = 0;
    int endp = 0;
    for (int b = 0; b < 200; ++b)
    {
      acc2 += kv.second[b];
      if (acc2 >= (long) (0.999 * tot)) { endp = b * 5; break; }
    }
    if (endp > 950) chosen.insert(kv.first);
  }
  // pass 2: load chosen
  std::map<int, std::vector<Cl>> evs;
  for (Long64_t i = 0; i < t->GetEntries(); ++i)
  {
    t->GetEntry(i);
    if ((int) ev == 44) continue;   // V6 laser veto (canon.h)
    if (lay < 7 || lay > 54 || adc <= 0) continue;
    if (!chosen.count((int) ev)) continue;
    evs[(int) ev].push_back({x, y, lay, tb, -1});
  }
  f->Close();
  std::vector<double> nfound, rmsv, rv;
  long nin = 0, ntot = 0;
  int showev = -1;
  std::vector<Trk> showT;
  std::vector<Cl> showC;
  for (auto &kv : evs)
  {
    std::vector<Trk> tr = hunt(kv.second);
    nfound.push_back((double) tr.size());
    printf("pixscan ev %d: %zu px -> %zu tracks\n", kv.first, kv.second.size(), tr.size());
    for (auto &T : tr)
    {
      ntot++;
      if (T.medtb >= 60 && T.medtb <= 360) nin++;
      rmsv.push_back(T.F.rms * 10);
      rv.push_back(T.F.R);
    }
    if (showev < 0 && tr.size() >= 40) { showev = kv.first; showT = tr; showC = kv.second; }
  }
  auto med = [](std::vector<double> v) {
    if (v.empty()) return 0.;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
  };
  FILE *fo = fopen(Form("ms_pixscan_%s.txt", ver), "w");
  auto P = [&](const char *fmt, ...) {
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    va_start(ap, fmt); vfprintf(fo, fmt, ap); va_end(ap);
  };
  P("ms_pixscan %s — exhaustive circle search on ALL real ntp_hit pixels\n", ver);
  P("complete events used: %zu | found/event median %.0f (mean %.1f)\n",
    evs.size(), med(nfound),
    nfound.empty() ? 0 : std::accumulate(nfound.begin(), nfound.end(), 0.) / nfound.size());
  P("in-time fraction %.2f | per-track RMS med %.0f um | R_fit med %.0f cm\n",
    ntot ? (double) nin / ntot : 0, med(rmsv) * 1000, med(rv));
  fclose(fo);
  gStyle->SetOptStat(0);
  TCanvas *cv = new TCanvas("cvpx", Form("pixscan %s", ver), 1500, 700);
  cv->Divide(2, 1);
  cv->cd(1);
  {
    TH1 *fr = gPad->DrawFrame(-80, -80, 80, 80,
                              Form("real event %d: ALL pixels + found tracks;x [cm];y [cm]", showev));
    fr->GetYaxis()->SetTitleOffset(1.25);
    TGraph *gu = new TGraph();
    for (size_t i = 0; i < showC.size(); i += 7)     // thin the pixel cloud for drawing
      gu->SetPoint(gu->GetN(), showC[i].x, showC[i].y);
    gu->SetMarkerStyle(1);
    gu->SetMarkerColorAlpha(kGray + 1, 0.35);
    gu->Draw("P same");
    int cols[3] = {kBlue + 1, kRed + 1, kGreen + 2};
    for (size_t k = 0; k < showT.size(); ++k)
    {
      TGraph *gt = new TGraph();
      for (int i : showT[k].idx) gt->SetPoint(gt->GetN(), showC[i].x, showC[i].y);
      gt->SetMarkerStyle(20);
      gt->SetMarkerSize(0.4);
      gt->SetMarkerColor(cols[k % 3]);
      gt->Draw("P same");
    }
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.034);
    tx.DrawLatex(0.13, 0.86, Form("%zu tracks found from raw pixels", showT.size()));
    tx.DrawLatex(0.13, 0.80, "(first complete event with >= 40 finds)");
    tx.DrawLatex(0.13, 0.74, "(pixel cloud drawn 1:7)");
  }
  cv->cd(2);
  {
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.038);
    tx.DrawLatex(0.06, 0.88, "pixel-level exhaustive search (real, complete events)");
    tx.SetTextSize(0.034);
    tx.DrawLatex(0.06, 0.78, Form("events %zu; found/event median %.0f", evs.size(), med(nfound)));
    tx.DrawLatex(0.06, 0.70, Form("in-time fraction %.2f", ntot ? (double) nin / ntot : 0));
    tx.DrawLatex(0.06, 0.62, Form("per-track RMS median %.0f #mum (pad centers, no centroiding)", med(rmsv) * 1000));
    tx.DrawLatex(0.06, 0.54, Form("R_{fit} median %.0f cm", med(rv)));
    tx.SetTextSize(0.030);
    tx.DrawLatex(0.06, 0.40, "same finder + acceptance as the cluster level;");
    tx.DrawLatex(0.06, 0.33, "tracks are findable straight from raw pixels when");
    tx.DrawLatex(0.06, 0.26, "the (x, y) circle and tbin-line are demanded JOINTLY");
  }
  cv->SaveAs(Form("../sim_validation_plots/ms_pixscan_%s.png", ver));
  printf("wrote ms_pixscan_%s outputs\n", ver);
}

// ---------------------------------------------------------------------------
// mt_rocscan — efficiency/purity working curves for the acceptance-bar cuts
// (user method, 2026-08-20): ONE exhaustive finder pass per dataset with the
// candidate collector armed (claiming stays bit-exact nominal), then an
// OFFLINE scan re-applies gate variants to the cached candidates.
//   Regime FINDER (the fake-circle filter the user asked about): sim = hunt()
//   over ALL island91-V6 clusters of 50 frames, candidates labeled by truth
//   majority; eff = findable truth tracks recovered (findable = nominal
//   acceptTrack on truth groups, FIXED across the scan; dupes not double-
//   counted); purity = accepted candidates with majority purity >= 0.6 on a
//   findable track (dupes count as real) / accepted. Real = hunt() over ALL
//   clusters of the complete-61 events; reference = audited-genuine seeds
//   (cleanFit-style, >=12 clusters); refind eff + composition proxies.
//   Regime GIVEN-GROUPING (the no-finder bar's home, my comparison): the
//   same collector run on truth groups (sim; genuine = gpt >= 0.19) and on
//   tracker seeds (real; genuine = cleanFit ok), same offline scan.
// Scan restricted to the post-cleaning gates; the 0.30 cm cleaning itself is
// fixed (rms variants above 0.30 would be no-ops by construction). Looser-
// than-nominal points carry a second-order claiming approximation (declared).
void mt_rocscan(const char *realf = "../clusters_seeds_island_79507-0.root_ntuplizer.root",
                const char *i91 = "island91_frames_production_v6.root",
                const char *ver = "v6", int nsimev = 50)
{
  using namespace MTK;
  struct SR { int n, nlay, span, maxgap; double R, rms, medtb, pur; bool nom; int lab; };
  // lab: FINDER sim: 1 = majority (pur>=0.6) on findable truth, 0 = else (ghost/sub-bar)
  //      FINDER real: 1 = majority overlap on reference seed, 0 = else
  //      GIVEN sim: 1 = genuine track-class (gpt>=0.19), 0 = else
  //      GIVEN real: 1 = cleanFit-genuine seed, 0 = fake
  std::vector<SR> fs, fr, gs, gr;                     // finder-sim/real, given-sim/real
  std::vector<int> fsTid;                             // finder-sim: majority tid key (dedup)
  long findableS = 0, refSeeds = 0;
  std::vector<int> fsFindKey;                         // per-candidate findable-truth key or -1
  std::vector<int> frSeedKey;                         // finder-real: matched reference seed key or -1

  // ---------------- SIM ----------------
  {
    TFile *f = TFile::Open(i91);
    if (!f || f->IsZombie()) { printf("missing %s\n", i91); return; }
    TTree *c = (TTree *) f->Get("ntp_cluster");
    TTree *u = (TTree *) f->Get("ntp_truth");
    float ev, lay, x, y, tb, tid, gpt;
    c->SetBranchStatus("*", 0);
    for (auto b : {"event", "layer", "x", "y", "tbin"}) c->SetBranchStatus(b, 1);
    c->SetBranchAddress("event", &ev);
    c->SetBranchAddress("layer", &lay);
    c->SetBranchAddress("x", &x);
    c->SetBranchAddress("y", &y);
    c->SetBranchAddress("tbin", &tb);
    u->SetBranchStatus("*", 0);
    for (auto b : {"gtrackID", "gpt"}) u->SetBranchStatus(b, 1);
    u->SetBranchAddress("gtrackID", &tid);
    u->SetBranchAddress("gpt", &gpt);
    std::map<int, std::vector<Cl>> evc;
    std::map<int, std::vector<int>> evtid;
    std::map<int, std::vector<float>> evgpt;
    for (Long64_t i = 0; i < c->GetEntries(); ++i)
    {
      c->GetEntry(i); u->GetEntry(i);
      if ((int) ev >= nsimev) continue;
      if (lay < 7 || lay > 54) continue;
      evc[(int) ev].push_back({x, y, lay, tb, -1});
      evtid[(int) ev].push_back((int) tid);
      evgpt[(int) ev].push_back(gpt);
    }
    f->Close();
    for (auto &kv : evc)
    {
      int e = kv.first;
      std::vector<Cl> &C = kv.second;
      std::vector<int> &T = evtid[e];
      std::vector<float> &G = evgpt[e];
      // truth groups
      std::map<int, std::vector<int>> tg;
      for (size_t i = 0; i < C.size(); ++i)
        if (T[i] > 0) tg[T[i]].push_back((int) i);
      // findable set (FIXED nominal definition) + GIVEN-regime records
      std::set<int> findable;
      gCand = nullptr;
      for (auto &g : tg)
      {
        if ((int) g.second.size() < 12) continue;
        std::vector<int> idx = g.second;
        Fit F;
        if (acceptTrack(C, idx, F)) { findable.insert(g.first); findableS++; }
      }
      std::vector<CandRec> cand;
      gCand = &cand;
      for (auto &g : tg)
      {
        if ((int) g.second.size() < 12) continue;
        std::vector<int> idx = g.second;
        Fit F;
        acceptTrack(C, idx, F);
      }
      gCand = nullptr;
      for (auto &r : cand)
      {
        double gp = 0;
        std::vector<double> tbv;
        for (int i : r.idx) tbv.push_back(C[i].tb);
        if (!r.idx.empty()) gp = G[r.idx[0]];
        std::sort(tbv.begin(), tbv.end());
        gs.push_back({r.n, r.nlay, r.span, r.maxgap, r.R, r.rms,
                      tbv.empty() ? 0 : tbv[tbv.size() / 2], 1.0, r.acc,
                      gp >= 0.19 ? 1 : 0});
      }
      // FINDER regime
      cand.clear();
      gCand = &cand;
      hunt(C);
      gCand = nullptr;
      for (auto &r : cand)
      {
        std::map<int, int> cnt;
        std::vector<double> tbv;
        for (int i : r.idx)
        {
          if (T[i] > 0) cnt[T[i]]++;
          tbv.push_back(C[i].tb);
        }
        int mtid = -1, mc = 0;
        for (auto &q : cnt)
          if (q.second > mc) { mc = q.second; mtid = q.first; }
        double pur = r.idx.empty() ? 0 : (double) mc / r.idx.size();
        bool onfind = pur >= 0.6 && findable.count(mtid);
        std::sort(tbv.begin(), tbv.end());
        fs.push_back({r.n, r.nlay, r.span, r.maxgap, r.R, r.rms,
                      tbv.empty() ? 0 : tbv[tbv.size() / 2], pur, r.acc,
                      onfind ? 1 : 0});
        fsFindKey.push_back(onfind ? e * 1000000 + mtid : -1);
      }
      printf("rocscan sim frame %d: %zu finder cands\n", e, cand.size());
    }
  }
  // ---------------- REAL ----------------
  {
    TFile *f = TFile::Open(realf);
    if (!f || f->IsZombie()) { printf("missing %s\n", realf); return; }
    // seeds + audit
    struct SCl { float x, y; int lay; };
    std::map<std::pair<int, int>, std::vector<SCl>> seeds;
    {
      TTree *t = (TTree *) f->Get("ntp_clus_trk");
      float ev, sid, lay, x, y;
      t->SetBranchStatus("*", 0);
      for (auto b : {"event", "seedID", "layer", "x", "y"}) t->SetBranchStatus(b, 1);
      t->SetBranchAddress("event", &ev);
      t->SetBranchAddress("seedID", &sid);
      t->SetBranchAddress("layer", &lay);
      t->SetBranchAddress("x", &x);
      t->SetBranchAddress("y", &y);
      for (Long64_t i = 0; i < t->GetEntries(); ++i)
      {
        t->GetEntry(i);
        if ((int) ev == 44) continue;                 // V6 laser veto
        if (lay < 7 || lay > 54) continue;
        seeds[{(int) ev, (int) sid}].push_back({x, y, (int) lay});
      }
    }
    // genuine reference seeds (cleanFit-style: iterative 0.30 cm, ok && >=6)
    std::map<std::pair<int, int>, int> refkey;        // seed -> reference index
    {
      int k = 0;
      for (auto &kv : seeds)
      {
        if ((int) kv.second.size() < 12) continue;
        std::vector<double> X, Y;
        for (auto &p : kv.second) { X.push_back(p.x); Y.push_back(p.y); }
        bool ok = true;
        Fit F;
        for (int it = 0; it < 6; ++it)
        {
          if ((int) X.size() < 6) { ok = false; break; }
          std::vector<double> Xc = X, Yc = Y;
          F = fitCircle(Xc, Yc);
          if (!F.ok) { ok = false; break; }
          std::vector<double> KX, KY;
          for (size_t i = 0; i < X.size(); ++i)
          {
            double res = std::hypot(X[i] - F.a, Y[i] - F.b) - F.R;
            if (std::fabs(res) <= 0.30) { KX.push_back(X[i]); KY.push_back(Y[i]); }
          }
          if (KX.size() == X.size()) break;
          X = KX; Y = KY;
        }
        ok = ok && (int) X.size() >= 6;
        // GIVEN-regime real record via collector on the seed's clusters
        std::vector<Cl> SC;
        for (auto &p : kv.second) SC.push_back({p.x, p.y, (float) p.lay, 0, -1});
        std::vector<int> idx(SC.size());
        for (size_t i = 0; i < SC.size(); ++i) idx[i] = (int) i;
        std::vector<CandRec> cand;
        gCand = &cand;
        Fit FF;
        acceptTrack(SC, idx, FF);
        gCand = nullptr;
        for (auto &r : cand)
          gr.push_back({r.n, r.nlay, r.span, r.maxgap, r.R, r.rms, 0, 1.0, r.acc, ok ? 1 : 0});
        if (ok) { refkey[kv.first] = k++; refSeeds++; }
      }
    }
    // complete-61 events + finder pass
    TTree *c = (TTree *) f->Get("ntp_cluster");
    float ev, lay, x, y, tb;
    c->SetBranchStatus("*", 0);
    for (auto b : {"event", "layer", "x", "y", "tbin"}) c->SetBranchStatus(b, 1);
    c->SetBranchAddress("event", &ev);
    c->SetBranchAddress("layer", &lay);
    c->SetBranchAddress("x", &x);
    c->SetBranchAddress("y", &y);
    c->SetBranchAddress("tbin", &tb);
    std::map<int, std::vector<Cl>> evc;
    std::map<int, std::vector<int>> tbc;
    for (Long64_t i = 0; i < c->GetEntries(); ++i)
    {
      c->GetEntry(i);
      if ((int) ev == 44) continue;                   // V6 laser veto
      if (lay < 7 || lay > 54) continue;
      evc[(int) ev].push_back({x, y, lay, tb, -1});
      auto &h = tbc[(int) ev];
      if (h.empty()) h.assign(200, 0);
      int b = (int) (tb / 5);
      if (b >= 0 && b < 200) h[b]++;
    }
    f->Close();
    for (auto &kv : evc)
    {
      long tot = 0, acc2 = 0;
      int endp = 0;
      for (int q : tbc[kv.first]) tot += q;
      for (int b = 0; b < 200; ++b)
      {
        acc2 += tbc[kv.first][b];
        if (acc2 >= (long) (0.999 * tot)) { endp = b * 5; break; }
      }
      if (endp <= 950) continue;                      // complete events only
      std::vector<Cl> &C = kv.second;
      // per-(layer) position lookup of this event's seed clusters
      std::map<int, std::vector<std::array<double, 3>>> lut;  // lay -> (x, y, seedkey)
      for (auto &sv : seeds)
      {
        if (sv.first.first != kv.first) continue;
        auto it = refkey.find(sv.first);
        if (it == refkey.end()) continue;
        for (auto &p : sv.second)
          lut[p.lay].push_back({(double) p.x, (double) p.y, (double) it->second});
      }
      std::vector<CandRec> cand;
      gCand = &cand;
      hunt(C);
      gCand = nullptr;
      for (auto &r : cand)
      {
        std::map<int, int> cnt;
        std::vector<double> tbv;
        for (int i : r.idx)
        {
          tbv.push_back(C[i].tb);
          auto lit = lut.find((int) C[i].lay);
          if (lit == lut.end()) continue;
          for (auto &p : lit->second)
            if (std::hypot(C[i].x - p[0], C[i].y - p[1]) < 0.1) { cnt[(int) p[2]]++; break; }
        }
        int msk = -1, mc = 0;
        for (auto &q : cnt)
          if (q.second > mc) { mc = q.second; msk = q.first; }
        double ov = r.idx.empty() ? 0 : (double) mc / r.idx.size();
        std::sort(tbv.begin(), tbv.end());
        fr.push_back({r.n, r.nlay, r.span, r.maxgap, r.R, r.rms,
                      tbv.empty() ? 0 : tbv[tbv.size() / 2], ov, r.acc,
                      ov >= 0.5 ? 1 : 0});
        frSeedKey.push_back(ov >= 0.5 ? msk : -1);
      }
      printf("rocscan real ev %d: %zu finder cands\n", kv.first, cand.size());
    }
  }
  // ---------------- offline scan ----------------
  struct Var { const char *fam; int n, nlay, span, gap; double rms, rmin; };
  std::vector<Var> vars = {
    {"NOMINAL", 12, 13, 15, 6, 0.20, 45},
    {"rms", 12, 13, 15, 6, 0.10, 45}, {"rms", 12, 13, 15, 6, 0.15, 45}, {"rms", 12, 13, 15, 6, 0.30, 45},
    {"gap", 12, 13, 15, 3, 0.20, 45}, {"gap", 12, 13, 15, 10, 0.20, 45}, {"gap", 12, 13, 15, 99, 0.20, 45},
    {"nlay", 12, 9, 15, 6, 0.20, 45}, {"nlay", 12, 17, 15, 6, 0.20, 45},
    {"span", 12, 13, 10, 6, 0.20, 45}, {"span", 12, 13, 20, 6, 0.20, 45},
    {"n", 16, 13, 15, 6, 0.20, 45}, {"n", 20, 13, 15, 6, 0.20, 45},
    {"Rmin", 12, 13, 15, 6, 0.20, 40}, {"Rmin", 12, 13, 15, 6, 0.20, 50}, {"Rmin", 12, 13, 15, 6, 0.20, 60}};
  auto pass = [](const SR &r, const Var &v) {
    return r.n >= v.n && r.nlay >= v.nlay && r.span >= v.span && r.maxgap <= v.gap &&
           r.rms <= v.rms && r.R >= v.rmin && r.R < 2e4;
  };
  FILE *fo = fopen(Form("ms_rocscan_%s.txt", ver), "w");
  auto P = [&](const char *fmt, ...) {
    va_list ap; va_start(ap, fmt); vprintf(fmt, ap); va_end(ap);
    va_start(ap, fmt); vfprintf(fo, fmt, ap); va_end(ap);
  };
  P("ms_rocscan %s — eff/purity working curves of the acceptance cuts (user method)\n", ver);
  P("finder regime: sim %zu cands / findable %ld | real %zu cands / reference seeds %ld (complete-61)\n",
    fs.size(), findableS, fr.size(), refSeeds);
  P("given-grouping regime: sim %zu truth groups (genuine = gpt>=0.19) | real %zu seeds (genuine = cleanFit ok)\n",
    gs.size(), gr.size());
  P("columns: FINDER sim eff | sim purity || real refind-eff | acc/ev | in-time || GIVEN sim eff | pur || real eff | pur\n");
  std::vector<std::array<double, 4>> rocF, rocG;      // for the figure: (effS, purS) x 2 regimes
  for (auto &v : vars)
  {
    // finder sim: eff via unique findable keys; purity over accepted
    std::set<int> rec;
    long accS = 0, realS = 0;
    for (size_t i = 0; i < fs.size(); ++i)
    {
      if (!pass(fs[i], v)) continue;
      accS++;
      if (fs[i].lab) { realS++; if (fsFindKey[i] >= 0) rec.insert(fsFindKey[i]); }
    }
    double effS = findableS ? (double) rec.size() / findableS : 0;
    double purS = accS ? (double) realS / accS : 0;
    // finder real
    std::set<int> refound;
    long accR = 0, intime = 0;
    for (size_t i = 0; i < fr.size(); ++i)
    {
      if (!pass(fr[i], v)) continue;
      accR++;
      if (fr[i].medtb >= 60 && fr[i].medtb <= 360) intime++;
      if (frSeedKey[i] >= 0) refound.insert(frSeedKey[i]);
    }
    double effR = refSeeds ? (double) refound.size() / refSeeds : 0;
    // given regimes
    auto gep = [&](std::vector<SR> &g, double &e, double &p) {
      long gen = 0, acc = 0, accgen = 0;
      for (auto &r : g)
      {
        if (r.lab) gen++;
        if (!pass(r, v)) continue;
        acc++;
        if (r.lab) accgen++;
      }
      e = gen ? (double) accgen / gen : 0;
      p = acc ? (double) accgen / acc : 0;
    };
    double eGS, pGS, eGR, pGR;
    gep(gs, eGS, pGS); gep(gr, eGR, pGR);
    P("%-8s n%2d L%2d s%2d g%2d rms%.2f R%2.0f : %.3f %.3f || %.3f %5.1f %.2f || %.3f %.3f || %.3f %.3f\n",
      v.fam, v.n, v.nlay, v.span, v.gap, v.rms, v.rmin,
      effS, purS, effR, (double) accR / 61., accR ? (double) intime / accR : 0,
      eGS, pGS, eGR, pGR);
    rocF.push_back({effS, purS, eGS, pGS});
  }
  P("sealed-calibration cross-check at NOMINAL: eff 0.777 / purity 0.752 (2026-07-31, same frames class)\n");
  P("NOTE: looser-than-nominal points (rms 0.30, gap 10/99, nlay 9, span 10, Rmin 40) carry the\n");
  P("claiming approximation; tighter points are exact subsets of the sealed candidate stream.\n");
  fclose(fo);
  printf("wrote ms_rocscan_%s.txt\n", ver);
}
