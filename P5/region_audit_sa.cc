// native region audit: after full app init (physics built), walk regions
#define main standalone_main_unused
#include "standalone_tpc.cc"
#undef main
#include <G4RegionStore.hh>
#include <G4ProductionCuts.hh>
int main()
{
  auto *rm = new G4RunManager;
  auto *sd = new TpcSD("region_audit_tmp.root");
  rm->SetUserInitialization(new Det("sphenix_p5.gdml",
    "/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/FIELDMAP_GAP/65/a9/65a930ed6de9c0e049cd0f3ef226e6b4_sphenix3dbigmapxyz_gap_rebuild_v2.root", sd));
  auto *phys = new FTFP_BERT(0);
  phys->RegisterPhysics(new G4StepLimiterPhysics);
  rm->SetUserInitialization(phys);
  rm->SetUserAction(new PrimGen(nullptr));
  rm->Initialize();
  for (auto *rg : *G4RegionStore::GetInstance())
  {
    G4ProductionCuts *pc = rg->GetProductionCuts();
    if (pc)
      printf("REGAUD %s : gamma %.3f e- %.3f e+ %.3f p %.3f mm (%zu rootLV)\n",
             rg->GetName().c_str(), pc->GetProductionCut(0) / CLHEP::mm,
             pc->GetProductionCut(1) / CLHEP::mm, pc->GetProductionCut(2) / CLHEP::mm,
             pc->GetProductionCut(3) / CLHEP::mm, rg->GetNumberOfRootVolumes());
    else
      printf("REGAUD %s : NO CUTS OBJECT (%zu rootLV)\n", rg->GetName().c_str(),
             rg->GetNumberOfRootVolumes());
  }
  return 0;
}
