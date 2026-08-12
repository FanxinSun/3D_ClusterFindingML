// gen_scan2.cc — 2-knob world-constrained retune scan (2026-08-08):
// (pT0Ref, StringPT:sigma) -> prints dNch/deta(|eta|<0.5) + the attribution
// pT-class counts (TPC-charge proxy inputs) + eps proxy, per point.
#include "Pythia8/Pythia.h"
#include <cmath>
#include <cstdio>
using namespace Pythia8;
int main(int argc, char **argv)
{
  int nev = argc > 1 ? atoi(argv[1]) : 6000;
  int seed = argc > 2 ? atoi(argv[2]) : 20260860;
  double pt0 = argc > 3 ? atof(argv[3]) : 2.28;
  double sps = argc > 4 ? atof(argv[4]) : 0.335;
  Pythia py;
  py.readString("Beams:idA = 2212");
  py.readString("Beams:idB = 2212");
  py.readString("Beams:eCM = 200.");
  py.readString("SoftQCD:inelastic = on");
  py.readString(("MultipartonInteractions:pT0Ref = " + std::to_string(pt0)).c_str());
  py.readString(("StringPT:sigma = " + std::to_string(sps)).c_str());
  py.readString("Random:setSeed = on");
  py.readString(("Random:seed = " + std::to_string(seed)).c_str());
  py.readString("Print:quiet = on");
  if (!py.init()) { fprintf(stderr, "init failed\n"); return 1; }
  double np[4] = {0, 0, 0, 0};
  long n05 = 0, fired = 0;
  for (int e = 0; e < nev;)
  {
    if (!py.next()) continue;
    int mbdN = 0, mbdS = 0;
    for (int i = 0; i < py.event.size(); ++i)
    {
      const Particle &p = py.event[i];
      if (!p.isFinal()) continue;
      double pt = p.pT(), eta = p.eta();
      np[pt < 0.2 ? 0 : (pt < 0.5 ? 1 : (pt < 1.0 ? 2 : 3))]++;
      if (p.isCharged())
      {
        if (std::fabs(eta) < 0.5) n05++;
        if (eta > 3.51 && eta < 4.61) mbdN++;
        if (eta < -3.51 && eta > -4.61) mbdS++;
      }
    }
    fired += (mbdN > 0 && mbdS > 0);
    e++;
  }
  printf("SCAN2 pt0 %.2f sig %.3f | dnde %.3f | classN %.3f %.3f %.3f %.3f | eps %.4f\n",
         pt0, sps, n05 / (1.0 * nev), np[0] / (double) nev, np[1] / (double) nev,
         np[2] / (double) nev, np[3] / (double) nev, (double) fired / nev);
  return 0;
}
