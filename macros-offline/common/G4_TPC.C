#ifndef MACRO_G4TPC_C
#define MACRO_G4TPC_C

#include <GlobalVariables.C>

#include <G4_ActsGeom.C>
#include <G4_Intt.C>
#include <G4_Mvtx.C>
#include <QA.C>

#include <g4tpc/PHG4TpcCentralMembrane.h>
#include <g4tpc/PHG4TpcDigitizer.h>
#include <g4tpc/PHG4TpcDirectLaser.h>
#include <g4tpc/PHG4TpcDistortion.h>
#include <g4tpc/PHG4TpcElectronDrift.h>
#include <g4tpc/PHG4TpcEndCapSubsystem.h>
#include <g4tpc/PHG4TpcPadPlane.h>
#include <g4tpc/PHG4TpcPadPlaneReadout.h>
#include <g4tpc/PHG4TpcSubsystem.h>

#include <g4main/PHG4Reco.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wundefined-internal"
#include <tpc/TpcClusterizer.h>
#include <tpc/TpcSimpleClusterizer.h>
#include <tpc/TpcClusterZCrossingCorrection.h>
#pragma GCC diagnostic pop

#include <tpc/TpcClusterCleaner.h>

#include <tpccalib/PHTpcCentralMembraneClusterizer.h>
#include <tpccalib/PHTpcCentralMembraneMatcher.h>
#include <tpccalib/TpcDirectLaserReconstruction.h>

#include <qa_modules/QAG4SimulationTpc.h>

#include <fun4all/Fun4AllServer.h>

R__LOAD_LIBRARY(libg4tpc.so)
R__LOAD_LIBRARY(libtpc.so)
R__LOAD_LIBRARY(libtpccalib.so)
R__LOAD_LIBRARY(libqa_modules.so)

namespace Enable
{
  bool TPC = false;
  bool TPC_ABSORBER = false;
  bool TPC_OVERLAPCHECK = false;
  bool TPC_CELL = false;
  bool TPC_CLUSTER = false;
  bool TPC_QA = false;

  bool TPC_ENDCAP = true;

  int TPC_VERBOSITY = 0;
}  // namespace Enable

namespace G4TPC
{
  int n_tpc_layer_inner = 16;
  int tpc_layer_rphi_count_inner = 1128; // 94 * 12
  int n_tpc_layer_mid = 16;
  int n_tpc_layer_outer = 16;
  int n_gas_layer = n_tpc_layer_inner + n_tpc_layer_mid + n_tpc_layer_outer;
  double tpc_outer_radius = 77. + 2.;

  // drift velocity is set here for all relevant modules
  double tpc_drift_velocity_sim= 0.007550;  // cm/ns   // this is the ArCF4Isobutane version of the gas
//  double tpc_drift_velocity_reco now set in GlobalVariables.C
//  double tpc_drift_velocity_reco= 0.007550;  // cm/ns   // this is the ArCF4Isobutane version of the gas

  // use simple clusterizer
  bool USE_SIMPLE_CLUSTERIZER = false;

  // distortions
  // static_distortion_filename is the map fed to PHG4TpcDistortion during electron drift
  // (SIMULATION). It is the static-only map MERGED with the TPC module-edge residual map
  // (see merge_module_edge.C) so the module-edge (RPhi) distortion is baked into the
  // simulated TRKR_HITSET/clusters rather than applied as a reco correction.
  // To revert to static-only, point this back at:
  //   .../5a/3d/5a3d0b9b5268b8bc6921ddd16e801c8f_static_only.distortion_map.hist.root
  bool ENABLE_STATIC_DISTORTIONS = true;
  auto static_distortion_filename = std::string("/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/TPC_STATIC_DISTORTION_with_ModuleEdge.root");

  bool ENABLE_TIME_ORDERED_DISTORTIONS = true;
  std::string time_ordered_distortion_filename = std::string("/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/TPC_TIMEORDERED_DISTORTION/51/72/51725df52f40ba2f3b5208fdd1bd39c2_TimeOrderedDistortions.root");

  // distortion corrections (applied during tracking/reco, not during G4 hit generation)
  bool ENABLE_CORRECTIONS = true;
  auto correction_filename = std::string("/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/TPC_STATIC_CORRECTION_MODEL/ec/7b/ec7bd756f9fc7274af6b479ee39580e3_static_only_inverted_10-new.root");

  // module-edge distortion corrections (loaded into slot 1 of TpcLoadDistortionCorrection)
  // DISABLED: module-edge is now applied during SIMULATION (merged into
  // static_distortion_filename above), not as a reco correction. Re-enabling this would
  // correct the module-edge distortion back out at tracking time, defeating the purpose.
  bool ENABLE_MODULE_EDGE_CORRECTIONS = false;
  auto module_edge_correction_filename = std::string("/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/TPC_Module_Edge/90/c9/90c95346a3e25a7a0753d71c46becbc6_residual_map_field_layers10iter_2D_acts_static_enabled.hist.root");

  // enable central membrane g4hits generation
  bool ENABLE_CENTRAL_MEMBRANE_HITS = false;
  
  // enable direct laser g4hits generation
  bool ENABLE_DIRECT_LASER_HITS = false;

  // save histograms
  bool DIRECT_LASER_SAVEHISTOGRAMS = false;

  // do cluster <-> hit association
  bool DO_HIT_ASSOCIATION = true;
  
  // space charge calibration output file
  std::string DIRECT_LASER_ROOTOUTPUT_FILENAME = "TpcSpaceChargeMatrices.root";
  std::string DIRECT_LASER_HISTOGRAMOUTPUT_FILENAME = "TpcDirectLaserReconstruction.root";

  // ---- 2026 gas mixture: ArCF4Isobutane (75:20:05) ----
  // Diffusion constants for ArCF4Isobutane
  double diffusion_long  = 0.014596;  // cm/sqrt(cm)
  double diffusion_trans = 0.005313;  // cm/sqrt(cm)

  // Physical TPC drift geometry (2026 measured)
  double maxDriftLength = 102.325;    // cm, CM face to top of GEM stack
  double CM_halfwidth   = 0.28;       // cm

}  // namespace G4TPC

void TPCInit()
{
  BlackHoleGeometry::max_radius = std::max(BlackHoleGeometry::max_radius, G4TPC::tpc_outer_radius);

  if (Enable::TPC_ENDCAP)
  {
    BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, 130.);
    BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, -130.);
  }
  else
  {
    BlackHoleGeometry::max_z = std::max(BlackHoleGeometry::max_z, 211. / 2.);
    BlackHoleGeometry::min_z = std::min(BlackHoleGeometry::min_z, -211. / 2.);
  }

  // the mvtx is not called if disabled but the default number of layers is set to 3,
  // so we need to set it to zero
  if (!Enable::MVTX)
  {
    G4MVTX::n_maps_layer = 0;
  }
  // same for the INTT
  if (!Enable::INTT)
  {
    G4INTT::n_intt_layer = 0;
  }

  // Set the (static) drift velocity in the cluster Z crossing correction module
  TpcClusterZCrossingCorrection::_vdrift = G4TPC::tpc_drift_velocity_reco;
}

//! TPC end cap, wagon wheel, electronics
void TPC_Endcaps(PHG4Reco* g4Reco)
{
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::TPC_OVERLAPCHECK;
  bool AbsorberActive = Enable::ABSORBER || Enable::TPC_ABSORBER;

  PHG4TpcEndCapSubsystem* tpc_endcap = new PHG4TpcEndCapSubsystem("TPC_ENDCAP");
  tpc_endcap->SuperDetector("TPC_ENDCAP");

  if (AbsorberActive) tpc_endcap->SetActive();
  tpc_endcap->OverlapCheck(OverlapCheck);

  //  tpc_endcap->set_int_param("construction_verbosity", 2);

  g4Reco->registerSubsystem(tpc_endcap);

  return;
}

double TPC(PHG4Reco* g4Reco,
           double radius)
{
  bool OverlapCheck = Enable::OVERLAPCHECK || Enable::TPC_OVERLAPCHECK;
  bool AbsorberActive = Enable::ABSORBER || Enable::TPC_ABSORBER;

  PHG4TpcSubsystem* tpc = new PHG4TpcSubsystem("TPC");
  tpc->SetActive();
  tpc->SuperDetector("TPC");
  tpc->set_double_param("steplimits", 1);  // 1cm steps

  // 2026 drift geometry: pass maxDriftLength and CM_halfwidth to the subsystem
  tpc->set_double_param("tpc_length", G4TPC::maxDriftLength * 2 + G4TPC::CM_halfwidth * 2);
  tpc->set_double_param("maxdriftlength", G4TPC::maxDriftLength);

  // Pass inner-layer binning to subsystem (modern macros do this at subsystem level)
  tpc->set_int_param("tpc_minlayer_inner", G4MVTX::n_maps_layer + G4INTT::n_intt_layer);
  tpc->set_int_param("ntpc_layers_inner", G4TPC::n_tpc_layer_inner);
  tpc->set_int_param("ntpc_phibins_inner", G4TPC::tpc_layer_rphi_count_inner);

  if (AbsorberActive)
    {
      tpc->SetAbsorberActive();
    }
  tpc->OverlapCheck(OverlapCheck);
  
  g4Reco->registerSubsystem(tpc);
  
  if (Enable::TPC_ENDCAP)
    {
      TPC_Endcaps(g4Reco);
    }
  
  radius = G4TPC::tpc_outer_radius;
  
  radius += no_overlapp;
  
  return radius;
}

void TPC_Cells()
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::TPC_VERBOSITY);
  auto se = Fun4AllServer::instance();

  // central membrane G4Hit generation
  if( G4TPC::ENABLE_CENTRAL_MEMBRANE_HITS )
  {
    auto centralMembrane = new PHG4TpcCentralMembrane;
    centralMembrane->setCentralMembraneDelay(0);
    se->registerSubsystem(centralMembrane);
  }

  // direct laser G4Hit generation
  if( G4TPC::ENABLE_DIRECT_LASER_HITS )
  {
    auto directLaser = new PHG4TpcDirectLaser;

    // setup phi and theta steps
    /* use 5deg steps */
    static constexpr double deg_to_rad = M_PI/180.;
    directLaser->SetPhiStepping( 72, 0*deg_to_rad, 360*deg_to_rad );
    directLaser->SetThetaStepping( 17, 5*deg_to_rad, 90*deg_to_rad );
    directLaser->SetDirectLaserAuto( true );
    directLaser->set_double_param("drift_velocity", G4TPC::tpc_drift_velocity_sim);
    se->registerSubsystem(directLaser);
  }

  //=========================
  // setup Tpc readout for filling cells
  // g4tpc/PHG4TpcElectronDrift uses
  // g4tpc/PHG4TpcPadPlaneReadout
  //=========================

  auto padplane = new PHG4TpcPadPlaneReadout;
  padplane->Verbosity(verbosity);

  auto edrift = new PHG4TpcElectronDrift;
  edrift->Detector("TPC");
  edrift->Verbosity(verbosity);
  if( G4TPC::ENABLE_STATIC_DISTORTIONS || G4TPC::ENABLE_TIME_ORDERED_DISTORTIONS )
  {
    auto distortionMap = new PHG4TpcDistortion;
    distortionMap->set_do_static_distortions( G4TPC::ENABLE_STATIC_DISTORTIONS );
    distortionMap->set_static_distortion_filename( G4TPC::static_distortion_filename );

    distortionMap->set_do_time_ordered_distortions( G4TPC::ENABLE_TIME_ORDERED_DISTORTIONS );
    distortionMap->set_time_ordered_distortion_filename( G4TPC::time_ordered_distortion_filename );

    distortionMap->Init();
    edrift->setTpcDistortion( distortionMap );
  }

  // 2026: use maxDriftLength instead of hardcoded 105.5 cm
  double tpc_readout_time = G4TPC::maxDriftLength / G4TPC::tpc_drift_velocity_sim;  // ns
  double extended_readout_time = 0.0;
  if(TRACKING::pp_mode) extended_readout_time = TRACKING::pp_extended_readout_time;
  if (gSystem->Getenv("TPC_EXTENDED_READOUT_NS"))
  {
    // Match the real streaming-style readout depth (run-79507 window ~ 51 us total = full
    // drift 13.97 us + ~37 us). Folded entirely into max_time; edrift's internal
    // extended_readout_time (default 7000 ns) is zeroed to avoid double-counting, and the
    // padplane gets the same value so its time axis (MaxT/NTBins) records the late hits.
    extended_readout_time = atof(gSystem->Getenv("TPC_EXTENDED_READOUT_NS"));
    edrift->set_double_param("extended_readout_time", 0);
    // ana.331 padplane has no 'extended_readout_time' param; its recording axis is
    // MaxT = 2*maxdriftlength/v_drift. Inflate its maxdriftlength (pure readout
    // bookkeeping -- true drift stays bounded by the G4 gas volume and edrift's
    // max_time) so late pileup hits get time bins: MaxT -> 2*drift + extended.
    padplane->set_double_param("maxdriftlength",
                               G4TPC::maxDriftLength + 0.5 * extended_readout_time * G4TPC::tpc_drift_velocity_sim);
  }
  edrift->set_double_param("max_time", tpc_readout_time + extended_readout_time);
  std::cout << "PHG4TpcElectronDrift readout window is from 0 to " <<  tpc_readout_time + extended_readout_time << std::endl;

  // override the default drift velocity parameter specification
  edrift->set_double_param("drift_velocity", G4TPC::tpc_drift_velocity_sim);
  padplane->SetDriftVelocity(G4TPC::tpc_drift_velocity_sim);

  // 2026 ArCF4Isobutane diffusion constants
  edrift->set_double_param("diffusion_long", G4TPC::diffusion_long);
  edrift->set_double_param("diffusion_trans", G4TPC::diffusion_trans);

  edrift->registerPadPlane(padplane);
  se->registerSubsystem(edrift);

  // The pad plane readout default is set in PHG4TpcPadPlaneReadout

  // We may want to change the number of inner layers, and can do that here
  padplane->set_int_param("tpc_minlayer_inner", G4MVTX::n_maps_layer + G4INTT::n_intt_layer);  // sPHENIX layer number of first Tpc readout layer
  padplane->set_int_param("ntpc_layers_inner", G4TPC::n_tpc_layer_inner);
  padplane->set_int_param("ntpc_phibins_inner", G4TPC::tpc_layer_rphi_count_inner);

  // Tpc digitizer
  //=========
  PHG4TpcDigitizer* digitpc = new PHG4TpcDigitizer();
  digitpc->SetTpcMinLayer(G4MVTX::n_maps_layer + G4INTT::n_intt_layer);
  double ENC = 670.0;  // standard
  digitpc->SetENC(ENC);
  double ADC_threshold = 4.0 * ENC;
  digitpc->SetADCThreshold(ADC_threshold);  // 4 * ENC seems OK
  digitpc->Verbosity(verbosity);
  cout << " Tpc digitizer: Setting ENC to " << ENC << " ADC threshold to " << ADC_threshold
       << " maps+Intt layers set to " << G4MVTX::n_maps_layer + G4INTT::n_intt_layer << endl;

  se->registerSubsystem(digitpc);

}

void TPC_Clustering()
{
  int verbosity = std::max(Enable::VERBOSITY, Enable::TPC_VERBOSITY);
  ACTSGEOM::ActsGeomInit();
  Fun4AllServer* se = Fun4AllServer::instance();

  //-------------
  // Cluster Hits
  //-------------

  // --- optional: mask dead TPC channels (real-data TPC_DEADCHANNELMAP payload from
  //     CDB_offline) BEFORE clustering, so sim occupancy sees the same dead regions.
  //     Env-gated: TPC_DEADMAP=<payload.root> + TRKRNTUP_SO + TPC_DEADMASK_H.
  if (gSystem->Getenv("TPC_DEADMAP") && gSystem->Getenv("TRKRNTUP_SO") && gSystem->Getenv("TPC_DEADMASK_H"))
  {
    gSystem->Load(gSystem->Getenv("TRKRNTUP_SO"));
    gROOT->ProcessLine(Form("#include \"%s\"", gSystem->Getenv("TPC_DEADMASK_H")));
    gROOT->ProcessLine(Form("{ TpcDeadChannelMasker* _dm = new TpcDeadChannelMasker(\"%s\");"
                            " Fun4AllServer::instance()->registerSubsystem(_dm); }",
                            gSystem->Getenv("TPC_DEADMAP")));
    std::cout << "[TpcDeadChannelMasker] registered, payload = " << gSystem->Getenv("TPC_DEADMAP") << std::endl;
  }

  // For the Tpc
  //==========
  if( G4TPC::USE_SIMPLE_CLUSTERIZER )
  {

    auto tpcclusterizer = new TpcSimpleClusterizer;
    tpcclusterizer->Verbosity(verbosity);
    se->registerSubsystem(tpcclusterizer);

  } else {

    auto tpcclusterizer = new TpcClusterizer;
    tpcclusterizer->Verbosity(verbosity);
    tpcclusterizer->set_cluster_version(G4TRACKING::cluster_version);
    tpcclusterizer->set_do_hit_association( G4TPC::DO_HIT_ASSOCIATION );
    {
      // Island-like size windows. ana.331 defaults (half 3/5 -> max 7 pads / 11-12 tbins)
      // guillotine exactly the merged/looper clusters dense AuAu events produce; the real
      // 2025 island reco is uncapped (phisize up to 97, zsize up to 126).
      int hphi = gSystem->Getenv("TPC_CLUS_HALF_PHI") ? atoi(gSystem->Getenv("TPC_CLUS_HALF_PHI")) : 16;
      int hz   = gSystem->Getenv("TPC_CLUS_HALF_Z")   ? atoi(gSystem->Getenv("TPC_CLUS_HALF_Z"))   : 40;
      tpcclusterizer->set_max_cluster_half_size_phi(hphi);
      tpcclusterizer->set_max_cluster_half_size_z(hz);
      std::cout << "TpcClusterizer max cluster half sizes: phi " << hphi << " z " << hz << std::endl;
    }
    se->registerSubsystem(tpcclusterizer);

  }
  
  if( !G4TPC::ENABLE_DIRECT_LASER_HITS )
  {
    auto tpcclustercleaner = new TpcClusterCleaner;
    tpcclustercleaner->Verbosity(verbosity);
    tpcclustercleaner->set_cluster_version(G4TRACKING::cluster_version);
    se->registerSubsystem(tpcclustercleaner);
  }

  // direct laser reconstruction
  if( G4TPC::ENABLE_DIRECT_LASER_HITS )
  { 
    auto directLaserReconstruction = new TpcDirectLaserReconstruction;
    directLaserReconstruction->set_outputfile( G4TPC::DIRECT_LASER_ROOTOUTPUT_FILENAME );
    directLaserReconstruction->set_savehistograms( G4TPC::DIRECT_LASER_SAVEHISTOGRAMS );
    directLaserReconstruction->set_histogram_outputfile( G4TPC::DIRECT_LASER_HISTOGRAMOUTPUT_FILENAME );
    se->registerSubsystem(directLaserReconstruction); 
  }
  
  // central membrane reconstruction
  if( G4TPC::ENABLE_CENTRAL_MEMBRANE_HITS )
  {
    // central membrane clusterizer
    se->registerSubsystem(new PHTpcCentralMembraneClusterizer);
    
    // match central membrane clusters to pads and generate distortion correction
    auto centralMembraneMatcher = new PHTpcCentralMembraneMatcher;
    se->registerSubsystem(centralMembraneMatcher);
  }
}


void TPC_QA()
{
  int verbosity = std::max(Enable::QA_VERBOSITY, Enable::TPC_VERBOSITY);

  Fun4AllServer* se = Fun4AllServer::instance();
  QAG4SimulationTpc * qa =  new QAG4SimulationTpc;
  qa->Verbosity(verbosity);
  qa->set_cluster_version(G4TRACKING::cluster_version);
  se->registerSubsystem(qa);
}


#endif
