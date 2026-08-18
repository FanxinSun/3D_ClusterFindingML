// p3_cmp3.C — three-line production comparison: every p3_cmp panel with BOTH
// real references (dual-reference directive, same convention as
// arrivals_bothrefs): REAL all-100 as recorded / REAL complete-62 windows /
// SIM. Pixel curves recomputed from source trees (manual fills — TTree::Draw
// >>h resolves the name via gDirectory and can miss stack/other-directory
// histograms; see the arrivals3 empty-sim incident, 2026-07-25). Island
// curves from island91_real filtered by the complete-62 event set.
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <algorithm>
#include <set>

void p3_cmp3(const char *simfile = "digi_frames_production_v53.root",
             const char *suffix = "_v53_bothrefs", const char *tag = "SIM v5.3 pp",
             const char *simisl = "island91_frames_production_v53.root")
{
  gROOT->SetBatch(1); gStyle->SetOptStat(0);
  TH1D h1a("h1a", "", 121, -0.5, 120.5), h1c("h1c", "", 121, -0.5, 120.5), h1s("h1s", "", 121, -0.5, 120.5);
  TH1D h2a("h2a", "", 1000, 0, 1000), h2c("h2c", "", 1000, 0, 1000), h2s("h2s", "", 1000, 0, 1000);
  TH1D h3a("h3a", "", 100, 0.5, 100.5), h3c("h3c", "", 100, 0.5, 100.5), h3s("h3s", "", 100, 0.5, 100.5);
  TH1D h4a("h4a", "", 25, 0.5, 25.5), h4c("h4c", "", 25, 0.5, 25.5), h4s("h4s", "", 25, 0.5, 25.5);
  TH1D h5a("h5a", "", 40, 0.5, 40.5), h5c("h5c", "", 40, 0.5, 40.5), h5s("h5s", "", 40, 0.5, 40.5);
  TH1D h6a("h6a", "", 150, 0, 6000), h6c("h6c", "", 150, 0, 6000), h6s("h6s", "", 150, 0, 6000);

  std::set<int> evset;  // complete-62 event ids (island filter)
  {  // real all-99 pixels (as recorded, laser event 44 vetoed — 2026-08-19 fix:
     // this block had read the raw ntuplizer without the veto)
    TFile *f = TFile::Open("/home/rog/sPHENIX/3D_ClusterFindingML/clusters_seeds_island_79507-0.root_ntuplizer.root");
    TTree *t = (TTree *) f->Get("ntp_hit");
    float lay, adc, tb, evr;
    t->SetBranchStatus("*", 0);
    for (auto b : {"event", "layer", "adc", "tbin"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("event", &evr); t->SetBranchAddress("layer", &lay); t->SetBranchAddress("adc", &adc); t->SetBranchAddress("tbin", &tb);
    for (Long64_t i = 0; i < t->GetEntries(); ++i)
    {
      t->GetEntry(i);
      if (lay < 7 || lay > 54 || adc <= 0 || (int) evr == 44) continue;
      h1a.Fill(adc); h2a.Fill(tb);
    }
    f->Close();
  }
  {  // real complete-62 pixels + event set
    TFile *f = TFile::Open("real_complete61_hits.root");  // laser-vetoed complete set (2026-08-17)
    TTree *t = (TTree *) f->Get("ntp_hit");
    float ev, lay, adc, tb;
    t->SetBranchStatus("*", 0);
    for (auto b : {"event", "layer", "adc", "tbin"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("event", &ev); t->SetBranchAddress("layer", &lay);
    t->SetBranchAddress("adc", &adc); t->SetBranchAddress("tbin", &tb);
    for (Long64_t i = 0; i < t->GetEntries(); ++i)
    {
      t->GetEntry(i);
      evset.insert((int) ev);
      if (lay < 7 || lay > 54 || adc <= 0) continue;
      h1c.Fill(adc); h2c.Fill(tb);
    }
    f->Close();
  }
  {  // sim pixels
    TFile *f = TFile::Open(simfile);
    TTree *t = (TTree *) f->Get("ntp_hit");
    float adc, tb;
    t->SetBranchStatus("*", 0);
    for (auto b : {"adc", "tbin"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("adc", &adc); t->SetBranchAddress("tbin", &tb);
    for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); h1s.Fill(adc); h2s.Fill(tb); }
    f->Close();
  }
  {  // real islands (all + complete-filtered)
    TFile *f = TFile::Open("island91_real.root");
    TTree *t = (TTree *) f->Get("ntp_cluster");
    float ev, sz, ps, zs, ad;
    t->SetBranchStatus("*", 0);
    for (auto b : {"event", "size", "phisize", "zsize", "adc"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("event", &ev); t->SetBranchAddress("size", &sz);
    t->SetBranchAddress("phisize", &ps); t->SetBranchAddress("zsize", &zs); t->SetBranchAddress("adc", &ad);
    for (Long64_t i = 0; i < t->GetEntries(); ++i)
    {
      t->GetEntry(i);
      h3a.Fill(sz); h4a.Fill(ps); h5a.Fill(zs); h6a.Fill(ad);
      if (evset.count((int) ev)) { h3c.Fill(sz); h4c.Fill(ps); h5c.Fill(zs); h6c.Fill(ad); }
    }
    f->Close();
  }
  {  // sim islands
    TFile *f = TFile::Open(simisl);
    TTree *t = (TTree *) f->Get("ntp_cluster");
    float sz, ps, zs, ad;
    t->SetBranchStatus("*", 0);
    for (auto b : {"size", "phisize", "zsize", "adc"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("size", &sz); t->SetBranchAddress("phisize", &ps);
    t->SetBranchAddress("zsize", &zs); t->SetBranchAddress("adc", &ad);
    for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); h3s.Fill(sz); h4s.Fill(ps); h5s.Fill(zs); h6s.Fill(ad); }
    f->Close();
  }
  printf("p3_cmp3: complete event set %zu events (laser-vetoed)\n", evset.size());

  auto d3 = [&](TVirtualPad *p, TH1D *a, TH1D *c, TH1D *s, const char *ti, bool logy) {
    p->cd();
    for (TH1D *h : {a, c, s}) if (h->Integral() > 0) h->Scale(1. / h->Integral());
    a->SetLineColor(kBlue + 1); a->SetLineWidth(2);
    c->SetLineColor(kAzure + 6); c->SetLineWidth(2); c->SetLineStyle(2);
    s->SetLineColor(kMagenta + 1); s->SetLineWidth(2);
    for (TH1D *h : {a, c, s}) h->SetStats(0);
    a->SetTitle(ti);
    double mx = std::max({a->GetMaximum(), c->GetMaximum(), s->GetMaximum()});
    a->SetMaximum(logy ? mx * 2.5 : mx * 1.3);
    if (logy) { gPad->SetLogy(); a->SetMinimum(1e-7); } else a->SetMinimum(0);
    a->Draw("HIST"); c->Draw("HIST SAME"); s->Draw("HIST SAME");
    TLegend *L = new TLegend(0.40, 0.70, 0.89, 0.89);
    L->SetBorderSize(0); L->SetFillStyle(0);
    L->AddEntry(a, "REAL all 99 non-laser (as recorded)", "l");
    L->AddEntry(c, "REAL 61 complete non-laser", "l");
    L->AddEntry(s, tag, "l");
    L->Draw();
  };
  TCanvas cv("cv", "", 1600, 900); cv.Divide(3, 2);
  d3(cv.cd(1), &h1a, &h1c, &h1s, "per-pixel ADC (ZS region);ADC;norm", true);
  d3(cv.cd(2), &h2a, &h2c, &h2s, "arrivals (adc>0);tbin;norm", false);
  d3(cv.cd(3), &h3a, &h3c, &h3s, "island size;pixels;norm", true);
  d3(cv.cd(4), &h4a, &h4c, &h4s, "island #phi-size;pads;norm", true);
  d3(cv.cd(5), &h5a, &h5c, &h5s, "island z-size;tbins;norm", true);
  d3(cv.cd(6), &h6a, &h6c, &h6s, "island ADC;raw sum;norm", true);
  cv.SaveAs(Form("/home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/p3_production_cmp%s.png", suffix));
  printf("p3_production_cmp%s.png saved (three-line dual-reference)\n", suffix);
}
