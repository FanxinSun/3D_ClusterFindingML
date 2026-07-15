// native-side probe: MapField (standalone app implementation) at the same points
#define main standalone_main_unused
#include "standalone_tpc.cc"
#undef main
int main()
{
  MapField f("/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/FIELDMAP_GAP/65/a9/65a930ed6de9c0e049cd0f3ef226e6b4_sphenix3dbigmapxyz_gap_rebuild_v2.root");
  double pts[9][3] = {{0,0,0},{40,0,0},{0,60,50},{33.7,-21.3,47.9},{-55.1,12.8,-88.6},
                      {70,0,100},{25,25,-25},{0,0,300},{150,150,0}};
  for (auto &q : pts)
  {
    double p[4] = {q[0] * cm, q[1] * cm, q[2] * cm, 0};
    double B[3] = {0, 0, 0};
    f.GetFieldValue(p, B);
    printf("FPROBE %7.1f %7.1f %7.1f : %+9.5f %+9.5f %+9.5f\n",
           q[0], q[1], q[2], B[0] / tesla, B[1] / tesla, B[2] / tesla);
  }
  return 0;
}
