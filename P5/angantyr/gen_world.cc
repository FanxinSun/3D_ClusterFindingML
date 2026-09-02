// gen_world.cc — generator-level confrontation with world pp-200 spectra
// (reshape campaign, 2026-09-02). Pythia 8 pp 200 GeV, SoftQCD:inelastic.
//   usage: gen_world <nev> <seed> <label> pt0 <pT0Ref>      (Monash + pT0Ref)
//          gen_world <nev> <seed> <label> cfg <file.cfg>    (readFile, e.g. official MDC2 = Detroit)
// Output <label>_world.txt: summary line + fine-binned (0.05 GeV, 0-5) invariant
// yields (1/Nev) d2N/(2 pi pT dpT dy) per species for two selections:
//   sel0: |y|<0.35, all inelastic events        (PHENIX PRC 83 064903 convention, sigma_inel 42 mb)
//   sel1: |y|<0.5,  NSD events only (no SD)     (STAR PLB 616 convention)
// p/pbar exclude hyperon feed-down (PHENIX/STAR feed-down corrected); pi/K inclusive.
// Also: dNch/deta(|eta|<0.5), MBD-fired fraction (3.51<|eta|<4.61 both arms),
// gen_scan2 attribution pT-class counts, charged |eta|<0.5 spectrum, <pT>.
#include "Pythia8/Pythia.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
using namespace Pythia8;
static const int NB = 100; static const double DPT = 0.05;
static bool hyperonAncestor(const Event &ev, int i, int depth = 0)
{
  if (depth > 30) return false;
  for (int m : ev[i].motherList())
  {
    if (m <= 0 || m >= ev.size()) continue;
    int a = std::abs(ev[m].id());
    if (a == 3122 || a == 3212 || a == 3222 || a == 3112 || a == 3312 || a == 3322 || a == 3334) return true;
    if (hyperonAncestor(ev, m, depth + 1)) return true;
  }
  return false;
}
int main(int argc, char **argv)
{
  if (argc < 6) { fprintf(stderr, "usage: gen_world <nev> <seed> <label> pt0|cfg <val>\n"); return 1; }
  int nev = atoi(argv[1]); int seed = atoi(argv[2]); std::string label = argv[3], mode = argv[4], val = argv[5];
  Pythia py;
  if (mode == "cfg" || mode == "cfgset")
  { // cfgset: "<file>|key=val;key=val" — official cfg base with overrides applied after
    std::string f = val, ex; size_t bar = val.find('|'); if (bar != std::string::npos) { f = val.substr(0, bar); ex = val.substr(bar + 1); }
    py.readFile(f);
    size_t a = 0; while (a < ex.size()) { size_t b = ex.find(';', a); if (b == std::string::npos) b = ex.size();
      std::string kv = ex.substr(a, b - a); size_t eq = kv.find('='); if (eq != std::string::npos) py.readString(kv.substr(0, eq) + " = " + kv.substr(eq + 1)); a = b + 1; }
  }
  else { py.readString("Beams:idA = 2212"); py.readString("Beams:idB = 2212"); py.readString("Beams:eCM = 200.");
         py.readString("SoftQCD:inelastic = on");
         if (mode == "pt0") py.readString("MultipartonInteractions:pT0Ref = " + val);
         else { size_t a = 0; while (a < val.size()) { size_t b = val.find(';', a); if (b == std::string::npos) b = val.size();
                std::string kv = val.substr(a, b - a); size_t eq = kv.find('='); if (eq != std::string::npos) py.readString(kv.substr(0, eq) + " = " + kv.substr(eq + 1)); a = b + 1; } } }
  py.readString("Random:setSeed = on"); py.readString("Random:seed = " + std::to_string(seed));
  py.readString("Print:quiet = on"); py.readString("Init:showChangedSettings = off"); py.readString("Next:numberShowInfo = 0");
  py.readString("Next:numberShowEvent = 0");
  if (!py.init()) { fprintf(stderr, "init failed\n"); return 1; }
  const int ids[8] = {211, -211, 321, -321, 2212, -2212, 2212, -2212};  // 6,7 = inclusive p/pbar (no feed-down cut)
  static double H[8][2][NB]; memset(H, 0, sizeof(H));
  static double Hch[NB]; memset(Hch, 0, sizeof(Hch)); static double Heta[24]; memset(Heta, 0, sizeof(Heta)); long ntpc = 0, ntpc_soft = 0;
  double np[4] = {0, 0, 0, 0}; long n05 = 0, fired = 0, nnsd = 0; double sumpt = 0; long nch = 0;
  double sptS[8][2] = {{0}}; long nS[8][2] = {{0}};
  for (int e = 0; e < nev;)
  {
    if (!py.next()) continue;
    int code = py.info.code(); bool sd = (code == 103 || code == 104); if (!sd) nnsd++;
    int mbdN = 0, mbdS = 0;
    for (int i = 0; i < py.event.size(); ++i)
    {
      const Particle &p = py.event[i]; if (!p.isFinal()) continue;
      double pt = p.pT(), eta = p.eta();
      np[pt < 0.2 ? 0 : (pt < 0.5 ? 1 : (pt < 1.0 ? 2 : 3))]++;
      if (p.isCharged())
      {
        if (std::fabs(eta) < 0.5) { n05++; sumpt += pt; nch++; int b = (int) (pt / DPT); if (b < NB) Hch[b]++; }
        if (std::fabs(eta) < 3.0) Heta[(int) (std::fabs(eta) / 0.25)]++;
        if (std::fabs(eta) < 1.1 && pt > 0.1) { ntpc++; if (pt < 0.3) ntpc_soft++; }
        if (eta > 3.51 && eta < 4.61) mbdN++;
        if (eta < -3.51 && eta > -4.61) mbdS++;
      }
      int s0 = -1; for (int k = 0; k < 6; ++k) if (p.id() == ids[k]) { s0 = k; break; }
      if (s0 < 0) continue;
      double y = p.y(); int b = (int) (pt / DPT); if (b >= NB) continue;
      bool fd = (s0 >= 4) && hyperonAncestor(py.event, i);
      for (int pass = 0; pass < 2; ++pass)
      {
        int s = s0; if (pass == 1) { if (s0 < 4) break; s = s0 + 2; } else if (fd) continue;  // pass1 = inclusive p/pbar
        if (std::fabs(y) < 0.35) { H[s][0][b]++; sptS[s][0] += pt; nS[s][0]++; }
        if (!sd && std::fabs(y) < 0.5) { H[s][1][b]++; sptS[s][1] += pt; nS[s][1]++; }
      }
    }
    fired += (mbdN > 0 && mbdS > 0); e++;
  }
  FILE *fo = fopen((label + "_world.txt").c_str(), "w");
  fprintf(fo, "# gen_world %s mode %s %s nev %d seed %d\n", label.c_str(), mode.c_str(), val.c_str(), nev, seed);
  fprintf(fo, "SUMMARY nev %d nnsd %ld dnde %.4f eps %.4f classN %.4f %.4f %.4f %.4f meanpt_ch05 %.4f\n",
          nev, nnsd, n05 / (double) nev, fired / (double) nev, np[0] / nev, np[1] / nev, np[2] / nev, np[3] / nev, sumpt / std::max(1L, nch));
  const char *nm[8] = {"pip", "pim", "kp", "km", "p", "pbar", "pinc", "pbarinc"}; const double DY[2] = {0.7, 1.0};
  for (int s = 0; s < 8; ++s) for (int sel = 0; sel < 2; ++sel)
  {
    double norm = (sel == 0 ? nev : nnsd);
    fprintf(fo, "SPEC %s sel%d meanpt %.4f yield %.5f", nm[s], sel, sptS[s][sel] / std::max(1L, nS[s][sel]), nS[s][sel] / norm);
    for (int b = 0; b < NB; ++b) { double pt = (b + 0.5) * DPT; fprintf(fo, " %.6e", H[s][sel][b] / (norm * 2 * M_PI * pt * DPT * DY[sel])); }
    fprintf(fo, "\n");
  }
  fprintf(fo, "ETA dNch/deta |eta| bins 0.25 (0-3):"); for (int b = 0; b < 12; ++b) fprintf(fo, " %.4f", Heta[b] / (nev * 2 * 0.25)); fprintf(fo, " | TPC-acc charged(|eta|<1.1,pT>0.1)/ev %.3f soft(pT<0.3) %.3f\n", ntpc / (double) nev, ntpc_soft / (double) nev);
  fprintf(fo, "SPEC chpm sel0eta05 meanpt %.4f yield %.5f", sumpt / std::max(1L, nch), nch / (double) nev);
  for (int b = 0; b < NB; ++b) { double pt = (b + 0.5) * DPT; fprintf(fo, " %.6e", Hch[b] / (nev * 2 * M_PI * pt * DPT * 1.0)); }
  fprintf(fo, "\n"); fclose(fo);
  printf("gen_world %s: nev %d nnsd %ld dnde %.3f eps %.3f -> %s_world.txt\n", label.c_str(), nev, nnsd, n05 / (double) nev, fired / (double) nev, label.c_str());
  return 0;
}
