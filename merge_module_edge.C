/// merge_module_edge.C
/// Merges TPC_Module_Edge distortion histograms into the static distortion map
/// so that PHG4TpcDistortion applies both during electron drift simulation.
///
/// The static map has TH3 histograms, the module-edge map has TH2 histograms.
/// For TH2 module-edge maps, the distortion is uniform in Z (same correction
/// at every drift distance), so we replicate the 2D values across all Z bins
/// of the 3D static histogram.
///
/// Usage: root -l -b -q merge_module_edge.C

void merge_module_edge()
{
  const char* static_file = "/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/TPC_STATIC_DISTORTION/5a/3d/5a3d0b9b5268b8bc6921ddd16e801c8f_static_only.distortion_map.hist.root";
  const char* module_edge_file = "/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/TPC_Module_Edge/90/c9/90c95346a3e25a7a0753d71c46becbc6_residual_map_field_layers10iter_2D_acts_static_enabled.hist.root";
  const char* output_file = "/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/TPC_STATIC_DISTORTION_with_ModuleEdge.root";

  // Open source files
  TFile* f_static = TFile::Open(static_file, "READ");
  TFile* f_edge   = TFile::Open(module_edge_file, "READ");

  if (!f_static || f_static->IsZombie()) { std::cerr << "Cannot open static file" << std::endl; return; }
  if (!f_edge || f_edge->IsZombie())     { std::cerr << "Cannot open module edge file" << std::endl; return; }

  // Create output file
  TFile* f_out = TFile::Open(output_file, "RECREATE");

  const char* suffixes[2] = {"_negz", "_posz"};
  const char* axes[3]     = {"R", "P", "Z"};

  for (int iz = 0; iz < 2; iz++)
  {
    for (int ia = 0; ia < 3; ia++)
    {
      TString hname = TString::Format("hIntDistortion%s%s", axes[ia], suffixes[iz]);

      // Get the 3D static histogram
      TH3* h3_static = dynamic_cast<TH3*>(f_static->Get(hname));
      if (!h3_static)
      {
        std::cerr << "WARNING: " << hname << " not found in static file, skipping" << std::endl;
        continue;
      }

      // Clone it to the output
      f_out->cd();
      TH3* h3_merged = dynamic_cast<TH3*>(h3_static->Clone(hname));
      h3_merged->SetDirectory(f_out);

      // Get the 2D module-edge histogram
      TH2* h2_edge = dynamic_cast<TH2*>(f_edge->Get(hname));
      if (!h2_edge)
      {
        std::cout << "INFO: " << hname << " not found in module-edge file (OK, using static only)" << std::endl;
        continue;
      }

      // The 3D static hist axes are (phi, r, z)
      // The 2D module-edge hist axes are (phi, r)
      // Add the 2D edge values to every Z slice of the 3D histogram
      int nbins_phi = h3_merged->GetNbinsX();
      int nbins_r   = h3_merged->GetNbinsY();
      int nbins_z   = h3_merged->GetNbinsZ();

      // Verify axis compatibility between 2D and 3D
      std::cout << hname << ": merging 2D module-edge ("
                << h2_edge->GetNbinsX() << "x" << h2_edge->GetNbinsY()
                << ") into 3D static ("
                << nbins_phi << "x" << nbins_r << "x" << nbins_z << ")" << std::endl;

      for (int iphi = 1; iphi <= nbins_phi; iphi++)
      {
        // Find the corresponding phi bin in the 2D hist
        double phi_center = h3_merged->GetXaxis()->GetBinCenter(iphi);
        int iphi_2d = h2_edge->GetXaxis()->FindBin(phi_center);

        for (int ir = 1; ir <= nbins_r; ir++)
        {
          double r_center = h3_merged->GetYaxis()->GetBinCenter(ir);
          int ir_2d = h2_edge->GetYaxis()->FindBin(r_center);

          // Get the module-edge distortion value at this (phi, r)
          double edge_val = h2_edge->GetBinContent(iphi_2d, ir_2d);

          // Add it to every Z bin
          for (int iz_bin = 1; iz_bin <= nbins_z; iz_bin++)
          {
            double old_val = h3_merged->GetBinContent(iphi, ir, iz_bin);
            h3_merged->SetBinContent(iphi, ir, iz_bin, old_val + edge_val);
          }
        }
      }

      std::cout << "  -> merged successfully" << std::endl;
    }
  }

  // Also copy any ReachesReadout histograms if they exist
  for (int iz = 0; iz < 2; iz++)
  {
    TString hname = TString::Format("hReachesReadout%s", suffixes[iz]);
    TH3* h_rr = dynamic_cast<TH3*>(f_static->Get(hname));
    if (h_rr)
    {
      f_out->cd();
      TH3* h_rr_clone = dynamic_cast<TH3*>(h_rr->Clone(hname));
      h_rr_clone->SetDirectory(f_out);
      std::cout << hname << ": copied from static file" << std::endl;
    }
  }

  f_out->Write();
  f_out->Close();
  f_static->Close();
  f_edge->Close();

  std::cout << "\nDone! Merged output written to:\n  " << output_file << std::endl;
}
