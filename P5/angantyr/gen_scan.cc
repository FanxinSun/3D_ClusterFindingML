// gen_scan.cc — generator-level pre-scan for the v5.3 soft-production tune
// (2026-07-25). No HepMC output: per-point statistics only.
//   usage: gen_scan <nev> <seed> <pT0Ref>
// Reports per-event means of ALL-final counts in the attribution pT classes
// (<0.2 / 0.2-0.5 / 0.5-1.0 / >1.0 — the classes carrying 29.3/43.4/22.7/4.6%
// of TPC deposited energy in the v5.2 attribution), forward fraction, charged
// multiplicity, and the MBD-proxy fired fraction. The content proxy =
// sum_p share_p * N_p(point)/N_p(default) is computed offline from these.
#include "Pythia8/Pythia.h"
#include <cmath>
#include <cstdio>
using namespace Pythia8;
int main(int argc, char **argv)
{
  int nev = argc > 1 ? atoi(argv[1]) : 4000;
  int seed = argc > 2 ? atoi(argv[2]) : 20260800;
  double pt0 = argc > 3 ? atof(argv[3]) : 2.28;
  Pythia py;
  py.readString("Beams:idA = 2212");
  py.readString("Beams:idB = 2212");
  py.readString("Beams:eCM = 200.");
  py.readString("SoftQCD:inelastic = on");
  py.readString(("MultipartonInteractions:pT0Ref = " + std::to_string(pt0)).c_str());
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
  printf("GENSCAN pT0Ref %.3f | Nfin %.2f Nch %.2f fwd-frac %.3f | classN %.3f %.3f %.3f %.3f | eps %.4f\n",
         pt0, nfin / nev, nch / nev, nfwd / nfin,
         np[0] / nev, np[1] / nev, np[2] / nev, np[3] / nev, (double) fired / nev);
  return 0;
}
