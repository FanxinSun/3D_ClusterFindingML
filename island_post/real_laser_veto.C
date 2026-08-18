// real_laser_veto.C — rebuild the REAL reference files with the GL1 laser event
// (event 44) vetoed, mirroring the official reco (2026-08-17). Backups keep the
// pre-veto files (*_prelaserveto.root). Products:
//   ref_real.root                 all-99 hist reference (make_ref.C logic + veto)
//   real_complete61_hits.root     slim ntp_hit of the 61 complete non-laser events
//   island_real.root              islands, event!=44 (all-99)
//   island_real_complete.root     islands, complete-61
//   island91_real.root            ntp_cluster (+ any other trees), event!=44
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TSystem.h>
#include <TEntryList.h>
#include <TDirectory.h>
#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>

static void filter_tree_file(const char *fn, const char *treename)
{
  TString bak = TString(fn).ReplaceAll(".root", "_prelaserveto.root");
  if (gSystem->AccessPathName(bak)) gSystem->CopyFile(fn, bak);
  TFile *fi = TFile::Open(bak);
  TTree *t = (TTree *) fi->Get(treename);
  TFile *fo = new TFile(fn, "RECREATE");
  TTree *tc = t->CopyTree("event!=44");
  printf("  %s/%s: %lld -> %lld rows (event!=44)\n", fn, treename, t->GetEntries(), tc->GetEntries());
  tc->Write();
  fo->Close();
  fi->Close();
}

void real_laser_veto()
{
  const char *RF = "/home/rog/sPHENIX/3D_ClusterFindingML/clusters_seeds_island_79507-0.root_ntuplizer.root";
  // ---- 1. slim complete-61 hits + all-99 histogram reference in one pass
  TFile *f = TFile::Open(RF);
  TTree *h = (TTree *) f->Get("ntp_hit");
  float event, layer, phibin, tb, adc, side, phi, z;
  h->SetBranchStatus("*", 0);
  for (const char *b : {"event", "layer", "phibin", "tbin", "adc", "zelem", "phi", "z"}) h->SetBranchStatus(b, 1);
  h->SetBranchAddress("event", &event); h->SetBranchAddress("layer", &layer);
  h->SetBranchAddress("phibin", &phibin); h->SetBranchAddress("tbin", &tb);
  h->SetBranchAddress("adc", &adc); h->SetBranchAddress("zelem", &side);
  h->SetBranchAddress("phi", &phi); h->SetBranchAddress("z", &z);
  // pass A: endpoints -> complete set
  std::map<int, std::vector<float>> tbs;
  Long64_t N = h->GetEntries();
  for (Long64_t i = 0; i < N; ++i)
  {
    h->GetEntry(i);
    if (layer < 7 || layer > 54 || adc <= 0) continue;
    tbs[(int) event].push_back(tb);
  }
  std::set<int> comp;
  for (auto &kv : tbs)
  {
    auto &v = kv.second; std::sort(v.begin(), v.end());
    if (v[(size_t) (0.999 * v.size())] > 950 && kv.first != 44) comp.insert(kv.first);
  }
  printf("real_laser_veto: complete non-laser events = %zu (was 62 with event 44)\n", comp.size());
  // pass B: write products
  if (gSystem->AccessPathName("ref_real_prelaserveto.root")) gSystem->CopyFile("ref_real.root", "ref_real_prelaserveto.root");
  TFile *o = new TFile("ref_real.root", "RECREATE");
  TH1D *r_adc = new TH1D("r_adc", "real per-pixel ADC;ADC", 1101, -0.5, 1100.5);
  TH1D *r_adcz = new TH1D("r_adcz", "real ADC near threshold;ADC", 121, -0.5, 120.5);
  TH1D *r_run = new TH1D("r_run", "real run length;consecutive tbins", 25, 0.5, 25.5);
  TH1D *r_lay = new TH1D("r_lay", "real hits/event/layer;layer", 48, 6.5, 54.5);
  TH1D *r_fold = new TH1D("r_fold", "real phi fold;phi mod 30deg", 120, 0, 0.5236);
  TH1D *r_tbin = new TH1D("r_tbin", "real tbin;tbin", 1000, 0, 1000);
  TFile *oc = new TFile("real_complete61_hits.root", "RECREATE");
  TTree *tc = new TTree("ntp_hit", "real complete-61 non-laser slim hits");
  float c_ev, c_lay, c_adc, c_ze, c_z, c_tb, c_phi;
  tc->Branch("event", &c_ev, "event/F"); tc->Branch("layer", &c_lay, "layer/F");
  tc->Branch("adc", &c_adc, "adc/F"); tc->Branch("zelem", &c_ze, "zelem/F");
  tc->Branch("z", &c_z, "z/F"); tc->Branch("tbin", &c_tb, "tbin/F"); tc->Branch("phi", &c_phi, "phi/F");
  std::map<uint64_t, std::vector<std::pair<int, int>>> cols;
  long nslim = 0;
  for (Long64_t i = 0; i < N; ++i)
  {
    h->GetEntry(i);
    int ev = (int) event;
    if (comp.count(ev))
    {
      c_ev = event; c_lay = layer; c_adc = adc; c_ze = side; c_z = z; c_tb = tb; c_phi = phi;
      tc->Fill(); nslim++;
    }
    if (ev == 44) continue;
    if (layer < 7 || layer > 54 || adc <= 0) continue;
    r_adc->Fill(adc); r_adcz->Fill(adc); r_lay->Fill(layer); r_tbin->Fill(tb);
    double ph = phi; while (ph < 0) ph += 2 * M_PI; r_fold->Fill(fmod(ph, 0.5235988));
    uint64_t k = ((uint64_t) (uint32_t) event << 24U) | ((uint64_t) (uint32_t) layer << 8U) | (uint64_t) (((int) side == 1) ? 1 : 0);
    cols[k].push_back({(int) phibin, (int) tb});
  }
  for (auto &g : cols)
  {
    auto &v = g.second; std::sort(v.begin(), v.end());
    int rl = 1;
    for (size_t i = 1; i < v.size(); ++i)
    {
      if (v[i].first == v[i - 1].first && v[i].second == v[i - 1].second + 1) rl++;
      else { r_run->Fill(rl); rl = 1; }
    }
    r_run->Fill(rl);
  }
  r_lay->Scale(1. / 99.);
  oc->cd(); tc->Write(); oc->Close();
  o->cd(); r_adc->Write(); r_adcz->Write(); r_run->Write(); r_lay->Write(); r_fold->Write(); r_tbin->Write(); o->Close();
  printf("real_laser_veto: ref_real.root rebuilt (99 events); real_complete61_hits.root %ld rows\n", nslim);
  f->Close();
  // ---- 2. island-level references
  filter_tree_file("island_real.root", "island");
  filter_tree_file("island_real_complete.root", "island");
  {  // island91_real: ntp_cluster + row-aligned ntp_truth -> same entry list on both
    const char *fn = "island91_real.root";
    TString bak = TString(fn).ReplaceAll(".root", "_prelaserveto.root");
    if (gSystem->AccessPathName(bak)) gSystem->CopyFile(fn, bak);
    TFile *fi = TFile::Open(bak);
    TTree *tc = (TTree *) fi->Get("ntp_cluster");
    TTree *tt = (TTree *) fi->Get("ntp_truth");
    tc->Draw(">>elist", "event!=44", "entrylist");
    TEntryList *el = (TEntryList *) gDirectory->Get("elist");
    TFile *fo = new TFile(fn, "RECREATE");
    tc->SetEntryList(el);
    TTree *tcc = tc->CopyTree("");
    TTree *ttc = nullptr;
    if (tt) { tt->SetEntryList(el); ttc = tt->CopyTree(""); }
    printf("  %s: ntp_cluster %lld -> %lld, ntp_truth %lld -> %lld (same entry list)\n", fn,
           tc->GetEntries(), tcc->GetEntries(), tt ? tt->GetEntries() : 0LL, ttc ? ttc->GetEntries() : 0LL);
    tcc->Write(); if (ttc) ttc->Write();
    fo->Close(); fi->Close();
  }
  printf("real_laser_veto: DONE\n");
}
