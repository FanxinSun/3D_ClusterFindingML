// dl_gen — run the collaboration's PHG4TpcDirectLaser VERBATIM (host stubs), one aim,
// dump hits in tpc_transport's ntp_g4hit format (layer from radius via geom table).
#include "PHG4TpcDirectLaser.h"
#include <phool/getClass.h>
#include <TFile.h>
#include <TNtuple.h>
#include <cmath>
#include <cstdio>
#include <vector>
static std::vector<std::pair<int,double>> LAYR;
static void loadGeom(const char* p)
{
  FILE* g = fopen(p, "r"); char line[256];
  while (fgets(line, 256, g))
  { int L, nb; double r, sl, p0, p1;
    if (line[0] == '#') continue;
    if (sscanf(line, "%d %d %lf %lf %lf %lf", &L, &nb, &r, &sl, &p0, &p1) == 6) LAYR.push_back({L, r}); }
  fclose(g);
}
static int layerFromR(double r)
{
  int best = -1; double bd = 1e9;
  for (auto& lr : LAYR) { double d = std::fabs(lr.second - r); if (d < bd) { bd = d; best = lr.first; } }
  return (bd < 0.7) ? best : -1;
}
PHG4HitContainer* g_stub_container = nullptr;
int main(int argc, char** argv)
{
  const double theta = (argc > 1) ? atof(argv[1]) : 30.;   // degrees
  const double phi = (argc > 2) ? atof(argv[2]) : 0.;
  loadGeom("/home/rog/sPHENIX/3D_ClusterFindingML/island_post/tpc_geom_table.txt");
  g_stub_container = new PHG4HitContainer();
  PHG4TpcDirectLaser dl("DL");
  dl.Detector("TPC");
  dl.set_double_param("arbitrary_theta", theta);
  dl.set_double_param("arbitrary_phi", phi);
  dl.InitRun(nullptr);
  dl.process_event(nullptr);
  printf("dl_gen: %zu hits (theta %.0f deg, phi %.0f deg, all 8 lasers)\n",
         g_stub_container->hits.size(), theta, phi);
  int nskip = 0;
  TFile f("laser_dl_g4hits.root", "RECREATE");
  TNtuple nt("ntp_g4hit", "DL tracks", "event:glayer:gx:gy:gz:gt:gpl:gpx:gpy:gpz:gedep:gtrackID");
  for (auto* h : g_stub_container->hits)
  {
    double dx = h->x[1] - h->x[0], dy = h->y[1] - h->y[0], dz = h->z[1] - h->z[0];
    double L = std::sqrt(dx * dx + dy * dy + dz * dz);
    double mx = 0.5 * (h->x[0] + h->x[1]), my = 0.5 * (h->y[0] + h->y[1]);
    int lay = layerFromR(std::sqrt(mx * mx + my * my));
    if (lay < 7) { nskip++; continue; }
    nt.Fill(0., (float) lay, mx, my, 0.5 * (h->z[0] + h->z[1]),
            h->t[0], L, L > 0 ? dx / L : 1., L > 0 ? dy / L : 0., L > 0 ? dz / L : 0.,
            h->edep, -1.);
  }
  nt.Write(); f.Close();
  printf("dl_gen: wrote laser_dl_g4hits.root (%d off-row steps skipped)\n", nskip);
  return 0;
}
