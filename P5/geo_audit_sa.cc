// native-side geometry audit: walk the GDML-loaded store, same record format
#define main standalone_main_unused
#include "standalone_tpc.cc"
#undef main
#include <G4LogicalVolumeStore.hh>
int main()
{
  G4GDMLParser p;
  p.Read("sphenix_p5.gdml", false);
  p.GetWorldVolume();
  for (auto *lv : *G4LogicalVolumeStore::GetInstance())
  {
    std::string n = lv->GetName();
    auto sp = n.find("0x");
    if (sp != std::string::npos) n = n.substr(0, sp);
    double cv = -1;
    try { cv = lv->GetSolid()->GetCubicVolume() / CLHEP::cm3; } catch (...) {}
    printf("GEOAUD %s | %s | %.6g\n", n.c_str(),
           lv->GetMaterial() ? lv->GetMaterial()->GetName().c_str() : "none", cv);
  }
  printf("GEOAUD-END %zu volumes\n", G4LogicalVolumeStore::GetInstance()->size());
  return 0;
}
