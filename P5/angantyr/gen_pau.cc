// Angantyr pAu 200 GeV MB -> minimal HepMC2-ASCII (exactly what the
// standalone_tpc reader consumes: E separators, one V(0,0,0,0), P lines
// "P bc pdg px py pz e m status"). No HepMC library needed.
#include "Pythia8/Pythia.h"
#include <cstdio>
using namespace Pythia8;
int main(int argc, char **argv)
{
  int nev = argc > 1 ? atoi(argv[1]) : 60;
  const char *out = argc > 2 ? argv[2] : "angantyr_pau200.dat";
  int seed = argc > 3 ? atoi(argv[3]) : 20260715;
  Pythia py;
  py.readString("Beams:idA = 2212");
  py.readString("Beams:idB = 1000791970");
  py.readString("Beams:eCM = 200.");
  py.readString("Beams:frameType = 1");
  py.readString("Random:setSeed = on");
  py.readString(("Random:seed = " + std::to_string(seed)).c_str());
  py.readString("Print:quiet = on");
  if (!py.init()) { fprintf(stderr, "init failed\n"); return 1; }
  FILE *f = fopen(out, "w");
  fprintf(f, "HepMC::Version 2.06.11\nHepMC::IO_GenEvent-START_EVENT_LISTING\n");
  int done = 0;
  while (done < nev)
  {
    if (!py.next()) continue;
    int nfin = 0, mbdN = 0, mbdS = 0;
    for (int i = 0; i < py.event.size(); ++i)
      if (py.event[i].isFinal())
      {
        nfin++;
        if (py.event[i].isCharged())
        {
          double eta = py.event[i].eta();
          if (eta > 3.51 && eta < 4.61) mbdN++;
          if (eta < -3.51 && eta > -4.61) mbdS++;
        }
      }
    fprintf(f, "E %d -1 -1 -1 -1 0 0 1 0 0 0 0\nU GEV MM\n", done + 1);
    fprintf(f, "V -1 0 0 0 0 0 0 %d 0\n", nfin);
    int bc = 1;
    for (int i = 0; i < py.event.size(); ++i)
    {
      const Particle &p = py.event[i];
      if (!p.isFinal()) continue;
      fprintf(f, "P %d %d %.8g %.8g %.8g %.8g %.8g 1 0 0 0 0\n",
              bc++, p.id(), p.px(), p.py(), p.pz(), p.e(), p.m());
    }
    fprintf(stderr, "FIRED %d %d\n", done + 1, (mbdN > 0 && mbdS > 0) ? 1 : 0);
    done++;
    if (done % 10 == 0) printf("GEN %d/%d\n", done, nev);
  }
  fprintf(f, "HepMC::IO_GenEvent-END_EVENT_LISTING\n");
  fclose(f);
  printf("GEN-DONE %d events -> %s\n", nev, out);
  return 0;
}
