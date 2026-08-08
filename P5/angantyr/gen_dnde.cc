// gen_dnde.cc — external cross-check vs world data (2026-08-08): mid-rapidity
// charged multiplicity density dNch/deta for INELASTIC pp 200 GeV, to compare
// against PHOBOS 2.29 +- 0.08 (and the RHIC world band). Modes:
//   gen_dnde <nev> <seed> pt0 <value>   — our line (Monash + pT0Ref)
//   gen_dnde <nev> <seed> cfg <file>    — official config verbatim
#include "Pythia8/Pythia.h"
#include <cmath>
#include <cstdio>
#include <cstring>
using namespace Pythia8;
int main(int argc, char **argv)
{
  int nev = argc > 1 ? atoi(argv[1]) : 4000;
  int seed = argc > 2 ? atoi(argv[2]) : 20260850;
  const char *mode = argc > 3 ? argv[3] : "pt0";
  const char *arg = argc > 4 ? argv[4] : "2.28";
  Pythia py;
  if (!strcmp(mode, "cfg"))
  {
    if (!py.readFile(arg)) { fprintf(stderr, "readFile failed\n"); return 1; }
  }
  else
  {
    py.readString("Beams:idA = 2212");
    py.readString("Beams:idB = 2212");
    py.readString("Beams:eCM = 200.");
    py.readString("SoftQCD:inelastic = on");
    py.readString(("MultipartonInteractions:pT0Ref = " + std::string(arg)).c_str());
  }
  py.readString("Random:setSeed = on");
  py.readString(("Random:seed = " + std::to_string(seed)).c_str());
  py.readString("Print:quiet = on");
  if (!py.init()) { fprintf(stderr, "init failed\n"); return 1; }
  long n05 = 0, n10 = 0;
  for (int e = 0; e < nev;)
  {
    if (!py.next()) continue;
    for (int i = 0; i < py.event.size(); ++i)
    {
      const Particle &p = py.event[i];
      if (!p.isFinal() || !p.isCharged()) continue;
      double ae = std::fabs(p.eta());
      if (ae < 0.5) n05++;
      if (ae < 1.0) n10++;
    }
    e++;
  }
  printf("DNDE %-4s %-28s : dNch/deta(|eta|<0.5) %.3f | (|eta|<1.0) %.3f  [inelastic, %d ev]\n",
         mode, arg, n05 / (1.0 * nev), n10 / (2.0 * nev), nev);
  return 0;
}
