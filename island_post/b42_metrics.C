// b42_metrics.C — one-line metric summary of a digi file for the B4.2 tail scan.
// Prints: kept/frame, pixel mean, sub-10 fraction, shoulder(10-19)/hi(20-40),
// run-length mean and tail P(run>=5/10/20/30)  [run = consecutive tbins per
// (side,layer,pad), event-bucketed sort, same definition as the real-side probe].
#include <TFile.h>
#include <TH1D.h>
#include <TNtuple.h>
#include <algorithm>
#include <cstdio>
#include <vector>

void b42_metrics(const char *fn = "b42_tmp.root", const char *tag = "cfg")
{
  TFile *f = TFile::Open(fn);
  TNtuple *t = (TNtuple *) f->Get("ntp_hit");
  float ev, lay, pb, tb, ze, adc;
  t->SetBranchStatus("*", 0);
  for (auto b : {"event", "layer", "phibin", "tbin", "zelem", "adc"}) t->SetBranchStatus(b, 1);
  t->SetBranchAddress("event", &ev);
  t->SetBranchAddress("layer", &lay);
  t->SetBranchAddress("phibin", &pb);
  t->SetBranchAddress("tbin", &tb);
  t->SetBranchAddress("zelem", &ze);
  t->SetBranchAddress("adc", &adc);
  TH1D hr("hr", "", 400, 0.5, 400.5), ha("ha", "", 1101, -0.5, 1100.5);
  double rn[3] = {0, 0, 0}, radc[3] = {0, 0, 0};
  std::vector<std::pair<long long, int>> buf;
  buf.reserve(3000000);
  int curev = -1;
  int evmin = 1 << 30, evmax = -(1 << 30);
  auto flush = [&]() {
    if (buf.empty()) return;
    std::sort(buf.begin(), buf.end());
    long long pk = -1;
    int pt = -1000000, run = 0;
    for (auto &p : buf)
    {
      if (p.first == pk && p.second == pt) continue;
      if (p.first == pk && p.second == pt + 1) { ++run; }
      else
      {
        if (run > 0) hr.Fill(run);
        run = 1;
      }
      pk = p.first;
      pt = p.second;
    }
    if (run > 0) hr.Fill(run);
    buf.clear();
  };
  const Long64_t N = t->GetEntries();
  for (Long64_t i = 0; i < N; ++i)
  {
    t->GetEntry(i);
    ha.Fill(adc);
    int rg = lay < 23 ? 0 : (lay < 39 ? 1 : 2);
    rn[rg] += 1;
    radc[rg] += adc;
    if ((int) ev != curev) { flush(); curev = (int) ev; }
    if ((int) ev < evmin) evmin = (int) ev;
    if ((int) ev > evmax) evmax = (int) ev;
    long long key = ((long long) (ze > 0.5)) * 100000000LL + (long long) lay * 1000000LL + (long long) pb;
    buf.emplace_back(key, (int) tb);
  }
  flush();
  double rt = hr.Integral(), at = ha.Integral();
  double nev = evmax - evmin + 1;  // appended blocks (minis) make transition-counting wrong
  double rtot = rn[0] + rn[1] + rn[2];
  printf("B43REGION %-14s shares R1 %.3f R2 %.3f R3 %.3f | pixm R1 %.1f R2 %.1f R3 %.1f\n",
         tag, rn[0] / rtot, rn[1] / rtot, rn[2] / rtot,
         radc[0] / rn[0], radc[1] / rn[1], radc[2] / rn[2]);
  printf("B42METRIC %-14s kept/fr %7.0f | pixmean %6.2f | sub10 %.2e | sh/hi %.3f | run mean %.3f P5 %.2e P10 %.2e P20 %.2e P30 %.2e\n",
         tag, at / nev, ha.GetMean(),
         ha.Integral(ha.FindBin(1), ha.FindBin(10)) / at,
         ha.Integral(ha.FindBin(10), ha.FindBin(19)) / ha.Integral(ha.FindBin(20), ha.FindBin(40)),
         hr.GetMean(), hr.Integral(hr.FindBin(5), 400) / rt, hr.Integral(hr.FindBin(10), 400) / rt,
         hr.Integral(hr.FindBin(20), 400) / rt, hr.Integral(hr.FindBin(30), 400) / rt);
  f->Close();
}
