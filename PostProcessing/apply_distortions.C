#include <iostream>
#include <TFile.h>
#include <TNtuple.h>
#include <TVector3.h>

#include "TpcDistortionCorrectionContainer.cc"
#include "TpcDistortionCorrection.cc"

void apply_distortions() {
    std::cout << "Loading Distortion Container..." << std::endl;
    TpcDistortionCorrectionContainer dcc;
    
    // Choose one of the static distortion maps you downloaded from the CDB
    std::string distortion_file = "/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/TPC_STATIC_DISTORTION/5a/3d/5a3d0b9b5268b8bc6921ddd16e801c8f_static_only.distortion_map.hist.root";
    
    std::cout << "Loading histograms from: " << distortion_file << std::endl;
    dcc.load_histograms(distortion_file);
    
    // Set parameters as they are set in the real gcc-12 code
    dcc.m_dimensions = 3;
    dcc.m_phi_hist_in_radians = true;
    dcc.m_interpolate_z = true;

    // Create the correction tool
    TpcDistortionCorrection tpc_distortion_correction;

    // Now let's open your idealized clusters NTuple
    std::string ntuple_file = "/home/rog/sPHENIX/3D_ClusterFindingML/clusters_seeds_island_79507-0.root_ntuplizer.root";
    TFile *f = TFile::Open(ntuple_file.c_str());
    if (!f || f->IsZombie()) {
        std::cout << "Error opening NTuple file" << std::endl;
        return;
    }
    
    TNtuple *ntp_cluster = (TNtuple*)f->Get("ntp_cluster");
    if (!ntp_cluster) {
        std::cout << "Error finding ntp_cluster tree" << std::endl;
        return;
    }
    
    float x, y, z;
    ntp_cluster->SetBranchAddress("x", &x);
    ntp_cluster->SetBranchAddress("y", &y);
    ntp_cluster->SetBranchAddress("z", &z);
    
    std::cout << "Testing distortion on 5 TPC hits..." << std::endl;
    int tpc_hits_found = 0;
    int i = 0;
    while (tpc_hits_found < 5 && i < ntp_cluster->GetEntries()) {
        ntp_cluster->GetEntry(i);
        
        TVector3 ideal_pos(x, y, z);
        
        // ONLY apply distortion if the hit is in the TPC (radius between ~30 and ~78 cm)
        double r = ideal_pos.Perp();
        if (r > 30.0 && r < 78.0) {
            // Apply the official sPHENIX distortion!
            TVector3 distorted_pos = tpc_distortion_correction.get_corrected_position(ideal_pos, &dcc);
            
            std::cout << "TPC Hit " << i << " (r=" << r << ")\n"
                      << "  Ideal:     (" << ideal_pos.X() << ", " << ideal_pos.Y() << ", " << ideal_pos.Z() << ")\n"
                      << "  Distorted: (" << distorted_pos.X() << ", " << distorted_pos.Y() << ", " << distorted_pos.Z() << ")\n"
                      << "  Delta:     (" << distorted_pos.X() - ideal_pos.X() << ", " << distorted_pos.Y() - ideal_pos.Y() << ", " << distorted_pos.Z() - ideal_pos.Z() << ")\n" << std::endl;
            tpc_hits_found++;
        }
        i++;
    }
    
    f->Close();
}
