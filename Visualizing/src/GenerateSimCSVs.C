// GenerateSimCSVs.C — SIM v4.0 exported in the exact schemas of the
// real-data CSVs in ../csv (GenerateEvent74CSVs.C + PlotEvent74Hits3D.C),
// so the existing plotting scripts run on sim by switching filenames.
//   sim_event<N>_hits.csv   event,layer,phi,tbin,adc,zelem,x,y,z,plot_x,plot_y,plot_z
//     SINGLE sim event N - the number declares the one-event content, like the
//     real event74_* convention; N is the SIM's own event index chosen by
//     criterion (43 = median-occupancy event), NEVER mirrored from a real event.
//   sim_event<N>_ntp_cluster.csv event,x,y,z,adc,zelem,layer
//   sim_event<N>_ntp_clus_trk.csv same, TRUTH track-labeled clusters (cls==0) —
//                           the sim analog of clusters-on-tracks, no tracking needed.
// plot projection identical: plot_x=tbin, plot_y=layer*cos(phi), plot_z=layer*sin(phi).
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <fstream>

void GenerateSimCSVs(int targetFrame = 43)  // default = median-occupancy sim event
{
  TString scriptPath = gSystem->DirName(gInterpreter->GetCurrentMacroName());
  if (scriptPath.IsNull()) scriptPath = ".";
  TString base = scriptPath + "/../../island_post/";

  // ---- hits ----
  {
    TFile *f = TFile::Open(base + "digi_frames_production_v40b.root");
    TTree *t = (TTree *) f->Get("ntp_hit");
    Float_t event, layer, phi, tbin, adc, zelem, z;
    t->SetBranchStatus("*", 0);
    for (auto b : {"event", "layer", "phi", "tbin", "adc", "zelem", "z"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("event", &event);
    t->SetBranchAddress("layer", &layer);
    t->SetBranchAddress("phi", &phi);
    t->SetBranchAddress("tbin", &tbin);
    t->SetBranchAddress("adc", &adc);
    t->SetBranchAddress("zelem", &zelem);
    t->SetBranchAddress("z", &z);
    // radius per layer from the transport geometry
    // loadGeo() reads tpc_geom_table.txt RELATIVE -> must run from island_post
    TString oldcwd = gSystem->pwd();
    gSystem->ChangeDirectory(base);
    gROOT->ProcessLine(".L tpc_digitize.C+");
    gROOT->ProcessLine("TDG::loadGeo();");
    double R[55] = {0};
    for (int L = 7; L <= 54; ++L)
      R[L] = *(double *) gROOT->ProcessLineFast(Form("&TDG::GEO[%d].radius;", L));
    if (R[7] < 1) { printf("GEO LOAD FAILED (R[7]=%f)\n", R[7]); return; }
    gSystem->ChangeDirectory(oldcwd);
    std::ofstream csv((scriptPath + Form("/../csv_sim/sim_event%d_hits.csv", targetFrame)).Data());
    csv << "event,layer,phi,tbin,adc,zelem,x,y,z,plot_x,plot_y,plot_z\n";
    // dedupe (iso/mini appended-block collisions; see island_post/sim_to_csv.C):
    // buffer the frame, merge adc on duplicate (layer,zelem,phibin,tbin)
    t->SetBranchStatus("phibin", 1);
    Float_t pbv; t->SetBranchAddress("phibin", &pbv);
    struct Row { float layer, phi, tbin, adc, zelem, z; };
    std::map<ULong64_t, Row> buf;
    long ndup = 0;
    for (Long64_t i = 0; i < t->GetEntries(); ++i)
    {
      t->GetEntry(i);
      if ((int) event != targetFrame || adc <= 0) continue;
      ULong64_t k = ((ULong64_t)(int) layer << 25) | ((ULong64_t)(int) zelem << 24) |
                    ((ULong64_t)(int) pbv << 12) | (ULong64_t)(int) tbin;
      auto ins = buf.emplace(k, Row{layer, phi, tbin, adc, zelem, z});
      if (!ins.second) { ins.first->second.adc += adc; ndup++; }
    }
    long n = 0;
    for (auto &kv : buf)
    {
      Row &r = kv.second;
      int L = (int) r.layer;
      double x = R[L] * TMath::Cos(r.phi), y = R[L] * TMath::Sin(r.phi);
      csv << targetFrame << "," << r.layer << "," << r.phi << "," << r.tbin << ","
          << r.adc << "," << r.zelem << "," << x << "," << y << "," << r.z << ","
          << r.tbin << "," << r.layer * TMath::Cos(r.phi) << "," << r.layer * TMath::Sin(r.phi) << "\n";
      n++;
    }
    printf("sim_event%d_hits.csv: %ld rows (%ld dupes merged)\n", targetFrame, n, ndup);
  }

  // ---- clusters (all) + truth-track clusters ----
  {
    TFile *f = TFile::Open(base + "island91_frames_production_v40b.root");
    TTree *c = (TTree *) f->Get("ntp_cluster");
    TTree *u = (TTree *) f->Get("ntp_truth");
    Float_t event, x, y, z, adc, zelem, layer, cls;
    c->SetBranchStatus("*", 0);
    for (auto b : {"event", "x", "y", "z", "adc", "zelem", "layer"}) c->SetBranchStatus(b, 1);
    c->SetBranchAddress("event", &event);
    c->SetBranchAddress("x", &x);
    c->SetBranchAddress("y", &y);
    c->SetBranchAddress("z", &z);
    c->SetBranchAddress("adc", &adc);
    c->SetBranchAddress("zelem", &zelem);
    c->SetBranchAddress("layer", &layer);
    u->SetBranchStatus("*", 0);
    u->SetBranchStatus("cls", 1);
    u->SetBranchAddress("cls", &cls);
    std::ofstream ca((scriptPath + Form("/../csv_sim/sim_event%d_ntp_cluster.csv", targetFrame)).Data());
    std::ofstream ct((scriptPath + Form("/../csv_sim/sim_event%d_ntp_clus_trk.csv", targetFrame)).Data());
    ca << "event,x,y,z,adc,zelem,layer\n";
    ct << "event,x,y,z,adc,zelem,layer\n";
    long na = 0, nt = 0;
    for (Long64_t i = 0; i < c->GetEntries(); ++i)
    {
      c->GetEntry(i);
      if ((int) event != targetFrame) continue;
      u->GetEntry(i);  // row-aligned truth
      ca << targetFrame << "," << x << "," << y << "," << z << ","
         << adc << "," << zelem << "," << layer << "\n";
      na++;
      if ((int) cls == 0)
      {
        ct << targetFrame << "," << x << "," << y << "," << z << ","
           << adc << "," << zelem << "," << layer << "\n";
        nt++;
      }
    }
    printf("sim_event%d_ntp_cluster.csv: %ld rows | sim_event%d_ntp_clus_trk.csv (truth track): %ld rows\n",
           targetFrame, na, targetFrame, nt);
  }
}
