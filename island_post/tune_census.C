// tune_census.C — per-candidate census for the v5.3 tune check: fired-mean
// and uniform-mean kept px per generator event, fired flags from the gen's
// 4-column FIRED file (1-based event ids).
#include <TFile.h>
#include <TTree.h>
#include <cstdio>
#include <map>

void tune_census(const char *censfile, const char *firedfile, const char *tag)
{
  std::map<int, int> fired;
  {
    FILE *f = fopen(firedfile, "r");
    char lab[16];
    int ev, n, s, fl;
    char line[256];
    while (f && fgets(line, 256, f))
      if (sscanf(line, "%15s %d %d %d %d", lab, &ev, &n, &s, &fl) == 5)
        fired[ev - 1] = fl;
    if (f) fclose(f);
  }
  TFile *f = TFile::Open(censfile);
  TTree *t = (TTree *) f->Get("ntp_hit");
  float ev;
  t->SetBranchStatus("*", 0);
  t->SetBranchStatus("event", 1);
  t->SetBranchAddress("event", &ev);
  std::map<int, long> n;
  for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); n[(int) ev]++; }
  double su = 0, sf = 0;
  long nu = 0, nf = 0;
  for (auto &kv : fired)
  {
    long c = n.count(kv.first) ? n[kv.first] : 0;
    su += c; nu++;
    if (kv.second) { sf += c; nf++; }
  }
  printf("TUNECHK %s: fired-mean %.0f (anchor 10897) | uniform-mean %.0f | eps %.4f | events %ld\n",
         tag, sf / nf, su / nu, (double) nf / nu, nu);
}
