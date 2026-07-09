// scan_day2.C — day-2 electronics calibration scan (stage-B replica, in-memory).
// Loads raw_pix once, sweeps (gaincal, thr, retention) against frozen real anchors
// from ref_real.root:  per-pixel ADC mean / mean(>30) / near-threshold shape,
// run-length mean & <=3 fraction. Prints ranked table; no files written.

#include <TFile.h>
#include <TH1D.h>
#include <TNtuple.h>
#include <TRandom3.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

namespace SCN
{
const double MV_PER_E = 7.68e-3;
const double ADU_PER_MV = 1024. / 2200.;
const double NOISE_ADU = 670. * (7.68e-3 / 2.4) * (1024. / 2200.);
const double PEDESTAL = 74.4;

struct Row
{
  uint64_t col;
  int tb;
  float q;
};
std::vector<Row> V;

void loadRaw(const char *rawin)
{
  TFile *fi = TFile::Open(rawin);
  TNtuple *r = (TNtuple *) fi->Get("raw_pix");
  float ev, L, sd, pad, tb, q, trk;
  r->SetBranchAddress("event", &ev);
  r->SetBranchAddress("layer", &L);
  r->SetBranchAddress("side", &sd);
  r->SetBranchAddress("pad", &pad);
  r->SetBranchAddress("tbin", &tb);
  r->SetBranchAddress("q", &q);
  r->SetBranchAddress("trk", &trk);
  V.reserve(r->GetEntries());
  for (Long64_t i = 0; i < r->GetEntries(); ++i)
  {
    r->GetEntry(i);
    uint64_t col = ((uint64_t) (uint32_t) ev << 40U) | ((uint64_t) (uint32_t) L << 32U) |
                   ((uint64_t) (((int) sd) ? 1 : 0) << 24U) | (uint64_t) (uint32_t) pad;
    V.push_back({col, (int) tb, q});
  }
  std::sort(V.begin(), V.end(), [](const Row &a, const Row &b) {
    return a.col == b.col ? a.tb < b.tb : a.col < b.col;
  });
  printf("scan: loaded %zu raw pixels\n", V.size());
}

struct Res
{
  double pixmean, mean30, runmean, runle3, nkept, chi2z, band, sub10, flat;
};

Res eval(double gaincal, double thr, int rpre, int rpost, TH1D *r_adcz, double thr2 = -1, TH1D **hout = nullptr, double sigma_pad = 0., double p_keep = 1.0, double sigma_thr = 0., double p2 = 0.)
{
  TRandom3 rng(4711);
  TH1D hz("hz", "", 121, -0.5, 120.5);
  double sum = 0, n = 0, sum30 = 0, n30 = 0;
  double runsum = 0, runn = 0, runle3 = 0;
  std::vector<double> adu;
  std::vector<char> keep;
  size_t i = 0;
  while (i < V.size())
  {
    size_t j = i;
    while (j < V.size() && V[j].col == V[i].col)
    {
      ++j;
    }
    adu.assign(j - i, 0.);
    keep.assign(j - i, 0);
    double gpad = 1.0;
    if (sigma_pad > 0)
    {
      TRandom3 gr((UInt_t) (V[i].col & 0xFFFFFFFFFFULL));  // static per (layer,side,pad)
      gpad = std::exp(gr.Gaus(0., sigma_pad) - 0.5 * sigma_pad * sigma_pad);
    }
    double dT = 0.;
    if (sigma_thr > 0)
    {
      TRandom3 tr((UInt_t) ((V[i].col & 0xFFFFFFFFFFULL) ^ 0x9E3779B9ULL));  // independent per-pad draw
      dT = tr.Gaus(0., sigma_thr);
    }
    for (size_t k = i; k < j; ++k)
    {
      double a = V[k].q * gpad * gaincal * MV_PER_E * ADU_PER_MV + PEDESTAL + rng.Gaus(0., NOISE_ADU);
      if (a > 1023)
      {
        a = 1023;
      }
      adu[k - i] = std::round(a - PEDESTAL);
    }
    const double t2 = ((thr2 >= 0) ? thr2 : thr) + dT;
    const double t1 = thr + dT;
    for (size_t k = 0; k < adu.size(); ++k)
    {
      if (adu[k] >= t1)
      {
        keep[k] = 1;
        for (int d = 1; d <= rpre; ++d)
        {
          if (k >= (size_t) d && V[i + k].tb - V[i + k - d].tb == d)
          {
            if (adu[k - d] >= t2 && (p_keep >= 1.0 || rng.Uniform() < p_keep))
            {
              keep[k - d] = 1;
            }
            else if (adu[k - d] >= 1 && p2 > 0 && rng.Uniform() < p2)
            {
              keep[k - d] = 1;
            }
          }
        }
        for (int d = 1; d <= rpost; ++d)
        {
          if (k + d < adu.size() && V[i + k + d].tb - V[i + k].tb == d)
          {
            if (adu[k + d] >= t2 && (p_keep >= 1.0 || rng.Uniform() < p_keep))
            {
              keep[k + d] = 1;
            }
            else if (adu[k + d] >= 1 && p2 > 0 && rng.Uniform() < p2)
            {
              keep[k + d] = 1;
            }
          }
        }
      }
    }
    // stats + run lengths over KEPT pixels (consecutive tbins)
    int rl = 0;
    int prev_tb = -99;
    for (size_t k = 0; k < adu.size(); ++k)
    {
      if (!keep[k] || adu[k] <= 0)
      {
        continue;
      }
      double a = adu[k];
      sum += a;
      n++;
      hz.Fill(a);
      if (a > 30)
      {
        sum30 += a;
        n30++;
      }
      int tb = V[i + k].tb;
      if (tb == prev_tb + 1)
      {
        rl++;
      }
      else
      {
        if (rl > 0)
        {
          runsum += rl;
          runn++;
          if (rl <= 3)
          {
            runle3++;
          }
        }
        rl = 1;
      }
      prev_tb = tb;
    }
    if (rl > 0)
    {
      runsum += rl;
      runn++;
      if (rl <= 3)
      {
        runle3++;
      }
    }
    i = j;
  }
  // near-threshold shape chi2 vs real (normalized, 10-120 ADU)
  double chi2 = 0;
  int ndf = 0;
  if (hz.Integral(hz.FindBin(8), hz.GetNbinsX()) > 0)
  {
    TH1D *rz = (TH1D *) r_adcz->Clone("rzc");
    int blo = hz.FindBin(8), bhi = hz.FindBin(30);
    double ri = rz->Integral(blo, bhi), si = hz.Integral(blo, bhi);
    for (int b = blo; b <= bhi; ++b)
    {
      double x = hz.GetBinContent(b) / si, y = rz->GetBinContent(b) / ri;
      double e2 = x / (si) + y / (ri);
      if (y > 0 && e2 > 0)
      {
        chi2 += (x - y) * (x - y) / e2;
        ndf++;
      }
    }
    delete rz;
  }
  if (hout)
  {
    *hout = (TH1D *) hz.Clone(Form("zs_%d", (int) (thr * 100 + thr2)));
    (*hout)->SetDirectory(nullptr);
  }
  double band = hz.Integral(hz.FindBin(11), hz.FindBin(19));
  double main = hz.Integral(hz.FindBin(20), hz.FindBin(40));
  double lo13 = hz.Integral(hz.FindBin(1), hz.FindBin(3));
  double lo810 = hz.Integral(hz.FindBin(8), hz.FindBin(10));
  double sub10 = n > 0 ? hz.Integral(hz.FindBin(1), hz.FindBin(10)) / n : 0;
  return {n ? sum / n : 0, n30 ? sum30 / n30 : 0, runn ? runsum / runn : 0,
          runn ? runle3 / runn : 0, n, ndf ? chi2 / ndf : 1e9, main > 0 ? band / main : -1,
          sub10, lo810 > 0 ? lo13 / lo810 : -1};
}
}  // namespace SCN
using namespace SCN;

void scan_day2(const char *rawin = "raw_pix_2ev.root")
{
  TFile *fr = TFile::Open("ref_real.root");
  TH1D *r_adc = (TH1D *) fr->Get("r_adc");
  TH1D *r_adcz = (TH1D *) fr->Get("r_adcz");
  TH1D *r_run = (TH1D *) fr->Get("r_run");
  double real_mean = r_adc->GetMean();
  TH1D *tmp = (TH1D *) r_adc->Clone("t30");
  tmp->GetXaxis()->SetRangeUser(30, 1100);
  double real_mean30 = tmp->GetMean();
  double real_runmean = r_run->GetMean();
  double real_runle3 = r_run->Integral(1, 3) / r_run->Integral();
  printf("REAL anchors: pixmean=%.2f mean(>30)=%.2f runmean=%.2f runle3=%.3f\n\n",
         real_mean, real_mean30, real_runmean, real_runle3);

  loadRaw(rawin);

  printf("== phase 1: gain scan (thr=30, ret 0/0; anchor mean(>30)) ==\n");
  printf("%8s %10s %10s\n", "gaincal", "mean(>30)", "target");
  double best_g = 1.0, best_d = 1e9;
  for (double g = 0.50; g <= 1.101; g += 0.05)
  {
    Res r = eval(g, 30., 0, 0, r_adcz);
    double d = std::fabs(r.mean30 - real_mean30);
    printf("%8.2f %10.2f %10.2f%s\n", g, r.mean30, real_mean30, d < best_d ? "  <-" : "");
    if (d < best_d)
    {
      best_d = d;
      best_g = g;
    }
  }
  printf("--> gaincal* = %.2f\n\n", best_g);

  printf("== phase 2: threshold x retention scan at gaincal*=%.2f ==\n", best_g);
  printf("%6s %5s/%-4s %9s %9s %9s %9s %11s %8s\n",
         "thr", "pre", "post", "pixmean", "runmean", "runle3", "chi2z/ndf", "pixels", "score");
  double bs = 1e9;
  double bthr = 15;
  int bpre = 0, bpost = 0;
  for (double thr : {5., 8., 10., 12., 15., 20.})
  {
    for (auto rp : {std::pair<int, int>{0, 0}, {1, 1}, {1, 2}, {2, 2}})
    {
      Res r = eval(best_g, thr, rp.first, rp.second, r_adcz);
      double score = std::fabs(r.runle3 - real_runle3) * 10 + std::fabs(r.runmean - real_runmean) +
                     std::min(r.chi2z, 50.) / 50.;
      printf("%6.0f %5d/%-4d %9.2f %9.2f %9.3f %9.2f %11.0f %8.3f%s\n",
             thr, rp.first, rp.second, r.pixmean, r.runmean, r.runle3, r.chi2z, r.nkept, score,
             score < bs ? "  <-" : "");
      if (score < bs)
      {
        bs = score;
        bthr = thr;
        bpre = rp.first;
        bpost = rp.second;
      }
    }
  }
  printf("\n--> phase-2 best: gaincal=%.2f thr=%.0f ret=%d/%d (score %.3f)\n", best_g, bthr, bpre, bpost, bs);

  printf("\n== phase 3: two-tier ZS (T1 x T2 x retention) — target: 10-20 shoulder + runs ==\n");
  printf("%5s %5s %4s/%-4s %9s %9s %9s %10s %11s %8s\n",
         "T1", "T2", "pre", "post", "pixmean", "runmean", "runle3", "chi2sh/ndf", "pixels", "score");
  double bs3 = 1e9, bT1 = 20, bT2 = 8;
  int b3pre = 1, b3post = 1;
  for (double T1 : {15., 18., 20., 22.})
  {
    for (double T2 : {6., 8., 10.})
    {
      for (auto rp : {std::pair<int, int>{1, 1}, {1, 2}, {2, 2}})
      {
        Res r = eval(best_g, T1, rp.first, rp.second, r_adcz, T2);
        double score = std::fabs(r.runle3 - real_runle3) * 10 + std::fabs(r.runmean - real_runmean) +
                       std::min(r.chi2z, 100.) / 25.;
        printf("%5.0f %5.0f %4d/%-4d %9.2f %9.2f %9.3f %10.2f %11.0f %8.3f%s\n",
               T1, T2, rp.first, rp.second, r.pixmean, r.runmean, r.runle3, r.chi2z, r.nkept, score,
               score < bs3 ? "  <-" : "");
        if (score < bs3)
        {
          bs3 = score;
          bT1 = T1;
          bT2 = T2;
          b3pre = rp.first;
          b3post = rp.second;
        }
      }
    }
  }
  printf("\n--> BEST two-tier: gaincal=%.2f T1=%.0f T2=%.0f ret=%d/%d (score %.3f vs phase-2 %.3f)\n",
         best_g, bT1, bT2, b3pre, b3post, bs3, bs);
}


// joint gain x ZS-model shoot-out on full statistics; writes zs_shapes.root
void zs_pick(const char *rawin = "raw_pix_20ev.root")
{
  TFile *fr = TFile::Open("ref_real.root");
  TH1D *r_adcz = (TH1D *) fr->Get("r_adcz");
  loadRaw(rawin);
  printf("%8s %-28s %9s %9s %9s\n", "gaincal", "config", "pixmean", "runmean", "runle3");
  for (double g : {0.70, 0.75, 0.80, 0.85, 0.90})
  {
    Res r1 = eval(g, 12, 0, 0, r_adcz);
    Res r2 = eval(g, 20, 1, 1, r_adcz, 10);
    printf("%8.2f %-28s %9.2f %9.2f %9.3f\n", g, "plain thr12", r1.pixmean, r1.runmean, r1.runle3);
    printf("%8.2f %-28s %9.2f %9.2f %9.3f\n", g, "two-tier T1=20 T2=10 r1/1", r2.pixmean, r2.runmean, r2.runle3);
  }
  printf("(real anchors: pixmean 99.15, runmean 2.70, runle3 0.816)\n");

  printf("\n== selective-retention corner (target shoulder ratio 0.23) ==\n");
  printf("%5s %5s %4s/%-4s %9s %9s %9s %9s\n","T1","T2","pre","post","pixmean","runmean","runle3","sh/main");
  for (double T2 : {10., 12., 14., 16.})
  {
    for (auto rp : {std::pair<int,int>{0,1},{1,1}})
    {
      TH1D *hh = nullptr;
      Res r = eval(0.80, 20, rp.first, rp.second, r_adcz, T2, &hh);
      double sh = hh->Integral(hh->FindBin(10), hh->FindBin(19));
      double mn = hh->Integral(hh->FindBin(20), hh->FindBin(40));
      printf("%5.0f %5.0f %4d/%-4d %9.2f %9.2f %9.3f %9.2f\n",
             20., T2, rp.first, rp.second, r.pixmean, r.runmean, r.runle3, mn > 0 ? sh / mn : -1);
      delete hh;
    }
  }
  // dump shoulder shapes at the gain that best matches pixmean for each model
  TH1D *h1 = nullptr, *h2 = nullptr;
  eval(0.80, 12, 0, 0, r_adcz, -1, &h1);   // config A
  eval(0.70, 20, 1, 1, r_adcz, 16, &h2);   // config B
  TFile *fo = new TFile("zs_shapes.root", "RECREATE");
  h1->Write("zs_plain");
  h2->Write("zs_twotier");
  r_adcz->Write("zs_real");
  fo->Close();
  printf("wrote zs_shapes.root\n");

  printf("\n== final joint re-anchor ==\n");
  printf("%8s %5s %5s %4s/%-4s %9s %9s %9s %9s\n","gaincal","T1","T2","pre","post","pixmean","runmean","runle3","sh/main");
  for (double g : {0.66, 0.70, 0.74, 0.78})
  {
    for (auto cfg : {std::tuple<double,int,int>{14.,0,1},{16.,1,1}})
    {
      TH1D *hh=nullptr;
      Res r = eval(g, 20, std::get<1>(cfg), std::get<2>(cfg), r_adcz, std::get<0>(cfg), &hh);
      double sh=hh->Integral(hh->FindBin(10),hh->FindBin(19)), mn=hh->Integral(hh->FindBin(20),hh->FindBin(40));
      printf("%8.2f %5.0f %5.0f %4d/%-4d %9.2f %9.2f %9.3f %9.2f\n",
             g, 20., std::get<0>(cfg), std::get<1>(cfg), std::get<2>(cfg),
             r.pixmean, r.runmean, r.runle3, mn>0?sh/mn:-1);
      delete hh;
    }
  }
}


void zs_bridge(const char *rawin = "raw_pix_20ev.root")
{
  TFile *fr = TFile::Open("ref_real.root");
  TH1D *r_adcz = (TH1D *) fr->Get("r_adcz");
  double rb = r_adcz->Integral(r_adcz->FindBin(11), r_adcz->FindBin(19)) /
              r_adcz->Integral(r_adcz->FindBin(20), r_adcz->FindBin(40));
  loadRaw(rawin);
  printf("real band(11-19)/main(20-40) = %.3f ; anchors pixmean 99.15 runmean 2.70 runle3 0.816\n\n", rb);
  printf("%6s %6s %5s %9s %9s %9s %9s\n", "sigpad", "p", "T2", "pixmean", "runmean", "runle3", "band");
  for (double sp : {0., 0.25, 0.35})
  {
    for (double p : {1.0, 0.5, 0.35})
    {
      Res r = eval(0.70, 20, 1, 1, r_adcz, 11, nullptr, sp, p);
      printf("%6.2f %6.2f %5.0f %9.2f %9.2f %9.3f %9.3f\n", sp, p, 11., r.pixmean, r.runmean, r.runle3, r.band);
    }
  }
}

void zs_refine(const char *rawin = "raw_pix_20ev.root")
{
  TFile *fr = TFile::Open("ref_real.root");
  TH1D *rz = (TH1D *) fr->Get("r_adcz");
  loadRaw(rawin);
  printf("%6s %6s %9s %9s %9s %9s\n", "g", "p", "pixmean", "runmean", "runle3", "band");
  for (double g : {0.68, 0.70, 0.72})
  {
    for (double p : {0.38, 0.40, 0.42})
    {
      Res r = eval(g, 20, 1, 1, rz, 11, nullptr, 0., p);
      printf("%6.2f %6.2f %9.2f %9.2f %9.3f %9.3f\n", g, p, r.pixmean, r.runmean, r.runle3, r.band);
    }
  }
}


void zs_b3(const char *rawin = "raw_pix_20ev.root")
{
  TFile *fr = TFile::Open("ref_real.root");
  TH1D *rz = (TH1D *) fr->Get("r_adcz");
  loadRaw(rawin);
  printf("targets: sub10 2.27e-4, flat 1.83, band 0.233, pixmean 99.15, run 2.70/0.816\n\n");
  printf("%6s %7s %10s %6s %7s %9s %9s %9s\n", "sigT", "p2", "sub10", "flat", "band", "pixmean", "runmean", "runle3");
  for (double sT : {0., 2., 3., 4.})
  {
    for (double p2 : {0., 0.002, 0.005, 0.01})
    {
      Res r = eval(0.70, 20, 1, 1, rz, 11, nullptr, 0., 0.39, sT, p2);
      printf("%6.1f %7.3f %10.2e %6.2f %7.3f %9.2f %9.2f %9.3f\n",
             sT, p2, r.sub10, r.flat, r.band, r.pixmean, r.runmean, r.runle3);
    }
  }
}

void zs_b3fine(const char *rawin = "raw_pix_20ev.root")
{
  TFile *fr = TFile::Open("ref_real.root");
  TH1D *rz = (TH1D *) fr->Get("r_adcz");
  loadRaw(rawin);
  printf("%8s %10s %6s %7s %9s %9s %9s\n", "p2", "sub10", "flat", "band", "pixmean", "runmean", "runle3");
  for (double p2 : {0.0004, 0.0005, 0.0006})
  {
    Res r = eval(0.70, 20, 1, 1, rz, 11, nullptr, 0., 0.39, 0., p2);
    printf("%8.4f %10.2e %6.2f %7.3f %9.2f %9.2f %9.3f\n",
           p2, r.sub10, r.flat, r.band, r.pixmean, r.runmean, r.runle3);
  }
}
