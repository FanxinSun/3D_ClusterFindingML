// Native Pythia8 pp 200 GeV MB (SoftQCD:inelastic = ND+SD+DD, diffraction ON —
// eps_MBD needs it) -> minimal HepMC2-ASCII, exactly what the standalone_tpc
// reader consumes (same writer as gen_pau.cc; species is the only physics change).
// MBD-fired proxy identical to v4.0: charged finals in 3.51<|eta|<4.61, both arms.
// stderr per event: "FIRED <ev> <nNorth> <nSouth> <fired>"  (4-column — carries
// the arm counts, unlike the old 2-column fired_*.txt of the Angantyr era).
#include "Pythia8/Pythia.h"
#include <cstdio>
#include <random>
using namespace Pythia8;
int main(int argc, char **argv)
{
  int nev = argc > 1 ? atoi(argv[1]) : 100;
  const char *out = argc > 2 ? argv[2] : "pp_run.dat";
  int seed = argc > 3 ? atoi(argv[3]) : 20260722;
  // optional 5th arg (v5.4 probe): vertex z-spread sigma in MM (0 = all
  // collisions at origin, byte-identical output to the pre-probe binary).
  // Vertices draw from an INDEPENDENT RNG so the Pythia stream — and thus
  // the particle content and FIRED flags — is unchanged at equal seed (CRN).
  double vzsig = argc > 5 ? atof(argv[5]) : 0.;
  std::mt19937 vrng((unsigned) seed ^ 0x5A5A5A5AU);
  std::normal_distribution<double> vgaus(0., 1.);
  if (vzsig > 0) fprintf(stderr, "VZSPREAD sigma %.1f mm\n", vzsig);
  Pythia py;
  py.readString("Beams:idA = 2212");
  py.readString("Beams:idB = 2212");
  py.readString("Beams:eCM = 200.");
  py.readString("Beams:frameType = 1");
  py.readString("SoftQCD:inelastic = on");
  if (argc > 4)  // v5.3 soft-production tune knob (omitted = Monash default,
  {              // preserving byte-level reproducibility of earlier productions)
    py.readString(("MultipartonInteractions:pT0Ref = " + std::string(argv[4])).c_str());
    fprintf(stderr, "TUNE pT0Ref %s\n", argv[4]);
  }
  py.readString("Random:setSeed = on");
  py.readString(("Random:seed = " + std::to_string(seed)).c_str());
  py.readString("Print:quiet = on");
  if (!py.init()) { fprintf(stderr, "init failed\n"); return 1; }
  FILE *f = fopen(out, "w");
  fprintf(f, "HepMC::Version 2.06.11\nHepMC::IO_GenEvent-START_EVENT_LISTING\n");
  int done = 0;
  long nfired = 0;
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
    double vz = vzsig > 0 ? vzsig * vgaus(vrng) : 0.;
    fprintf(f, "V -1 0 0 0 %.6g 0 0 %d 0\n", vz, nfin);
    int bc = 1;
    for (int i = 0; i < py.event.size(); ++i)
    {
      const Particle &p = py.event[i];
      if (!p.isFinal()) continue;
      fprintf(f, "P %d %d %.8g %.8g %.8g %.8g %.8g 1 0 0 0 0\n",
              bc++, p.id(), p.px(), p.py(), p.pz(), p.e(), p.m());
    }
    int fired = (mbdN > 0 && mbdS > 0) ? 1 : 0;
    nfired += fired;
    fprintf(stderr, "FIRED %d %d %d %d\n", done + 1, mbdN, mbdS, fired);
    done++;
    if (done % 500 == 0) printf("GEN %d/%d\n", done, nev);
  }
  fprintf(f, "HepMC::IO_GenEvent-END_EVENT_LISTING\n");
  fclose(f);
  printf("GEN-DONE %d events -> %s | MBD-fired %ld (%.3f)\n", nev, out, nfired,
         (double) nfired / nev);
  return 0;
}
