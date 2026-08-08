// gen_scan_cfg.cc — probe-1 variant of gen_scan (2026-07-30): identical
// meters, but the generator is configured from a Pythia .cfg FILE (readFile)
// instead of the hardcoded settings — so the OFFICIAL sPHENIX pp MB config
// (calibrations: Generators/HeavyFlavor_TG/phpythia8_minBias_MDC2.cfg) runs
// verbatim. Seed/quiet are applied AFTER the file so ours win.
//   usage: gen_scan_cfg <nev> <seed> <cfgfile> <label>
#include "Pythia8/Pythia.h"
#include <cmath>
#include <cstdio>
using namespace Pythia8;
int main(int argc, char **argv)
{
  int nev = argc > 1 ? atoi(argv[1]) : 4000;
  int seed = argc > 2 ? atoi(argv[2]) : 20260850;
  const char *cfg = argc > 3 ? argv[3] : "official_pp_mb_mdc2.cfg";
  const char *lab = argc > 4 ? argv[4] : "cfg";
  Pythia py;
  if (!py.readFile(cfg)) { fprintf(stderr, "readFile %s failed\n", cfg); return 1; }
  py.readString("Random:setSeed = on");
  py.readString(("Random:seed = " + std::to_string(seed)).c_str());
  py.readString("Print:quiet = on");
  if (!py.init()) { fprintf(stderr, "init failed\n"); return 1; }
  double np[4] = {0, 0, 0, 0}, nch = 0, nfwd = 0, nfin = 0;
  long fired = 0;
  for (int e = 0; e < nev;)
  {
    if (!py.next()) continue;
    int mbdN = 0, mbdS = 0;
    for (int i = 0; i < py.event.size(); ++i)
    {
      const Particle &p = py.event[i];
      if (!p.isFinal()) continue;
      nfin++;
      double pt = p.pT(), eta = p.eta();
      int ip = pt < 0.2 ? 0 : (pt < 0.5 ? 1 : (pt < 1.0 ? 2 : 3));
      np[ip]++;
      if (std::fabs(eta) > 1.1) nfwd++;
      if (p.isCharged())
      {
        nch++;
        if (eta > 3.51 && eta < 4.61) mbdN++;
        if (eta < -3.51 && eta > -4.61) mbdS++;
      }
    }
    fired += (mbdN > 0 && mbdS > 0);
    e++;
  }
  printf("GENSCAN %-14s | Nfin %.2f Nch %.2f fwd-frac %.3f | classN %.3f %.3f %.3f %.3f | eps %.4f\n",
         lab, nfin / nev, nch / nev, nfwd / nfin,
         np[0] / nev, np[1] / nev, np[2] / nev, np[3] / nev, (double) fired / nev);
  return 0;
}
