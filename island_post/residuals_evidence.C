// residuals_evidence.C — ORIGINAL-PURPOSE producer for sim_validation_plots/residuals_evidence.png
//
// Day-2 origin (then an inline throwaway; persistent now per the provenance rule):
// evidence panels for the two NAMED residuals of that era, same variables/binning/axes:
//   (1) island #phi-size overlay, unit-normalized, log-y  — "THE width residual" (-14% then)
//   (2) islands/event per layer, per-frame norm, log-y    — "THE content residual" (flat x4.7 then)
// Replotting under a new config keeps the panels and lets the legends carry the
// re-measured numbers. Island trees are already canonical (islandize applies the
// TPC cut and the side-0 z fix on ingestion), so no extra cuts here — see canon.h.
//
// v33_b42 edition: REAL islands (island_real.root, 100 frames) vs SIM v3.3 B4.2
// (island_frames_v33.root, 250 composed frames at the P1-derived rate).
#include "canon.h"

void residuals_evidence(const char *realisl = "island_real.root",
                        const char *simisl = "island_frames_v40b.root",
                        const char *tag = "SIM v4.0", const char *suffix = "")
{
  gROOT->SetBatch(1);
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.05);

  TFile *fi = TFile::Open(realisl);
  TTree *tr = (TTree *) fi->Get("island");
  TFile *fd = TFile::Open(simisl);
  TTree *td = (TTree *) fd->Get("island");
  // distinct-event counts (range arithmetic breaks on non-contiguous subsets
  // such as the complete-62 real reference)
  auto nev = [](TTree *t) {
    float ev; t->SetBranchStatus("*", 0);
    t->SetBranchStatus("event", 1); t->SetBranchAddress("event", &ev);
    std::set<int> s;
    for (Long64_t i = 0; i < t->GetEntries(); ++i) { t->GetEntry(i); s.insert((int) ev); }
    t->SetBranchStatus("*", 1);
    return (double) s.size();
  };
  const double nR = nev(tr);
  const double nS = nev(td);

  TCanvas c("c", "", 1300, 520);
  c.Divide(2, 1);

  // panel 1: phi-size shape overlay (the width residual)
  c.cd(1);
  gPad->SetLogy();
  TH1D *p1 = new TH1D("p1", "island #phi-size: THE width residual;#phi-size [pads];norm", 25, 0.5, 25.5);
  TH1D *p2 = new TH1D("p2", "", 25, 0.5, 25.5);
  tr->Draw("phisize>>p1", "", "goff");
  td->Draw("phisize>>p2", "", "goff");
  p1->Scale(1. / p1->Integral());
  p2->Scale(1. / p2->Integral());
  p1->SetLineColor(kBlue + 1);
  p2->SetLineColor(kMagenta + 1);
  p1->SetLineWidth(2);
  p2->SetLineWidth(2);
  p1->SetMaximum(std::max(p1->GetMaximum(), p2->GetMaximum()) * 2);
  p1->SetMinimum(1e-6);
  p1->SetStats(0);
  p1->Draw("HIST");
  p2->Draw("HIST SAME");
  const double m1 = p1->GetMean(), m2 = p2->GetMean();
  TLegend *L1 = new TLegend(0.35, 0.72, 0.89, 0.89);
  L1->SetBorderSize(0);
  L1->SetFillStyle(0);
  L1->AddEntry(p1, Form("REAL  <#phi-size>=%.2f", m1), "l");
  L1->AddEntry(p2, Form("%s  <#phi-size>=%.2f (%+.1f%%)", tag, m2, 100. * (m2 - m1) / m1), "l");
  L1->Draw();

  // panel 2: per-frame island rate by layer (the content residual)
  c.cd(2);
  gPad->SetLogy();
  TH1D *q1 = new TH1D("q1", "islands/event per layer: THE content residual;TPC layer;islands/event", 48, 6.5, 54.5);
  TH1D *q2 = new TH1D("q2", "", 48, 6.5, 54.5);
  tr->Draw("layer>>q1", "", "goff");
  td->Draw("layer>>q2", "", "goff");
  q1->Scale(1. / nR);
  q2->Scale(1. / nS);
  q1->SetLineColor(kBlue + 1);
  q2->SetLineColor(kMagenta + 1);
  q1->SetLineWidth(2);
  q2->SetLineWidth(2);
  q1->SetStats(0);
  q1->SetMaximum(std::max(q1->GetMaximum(), q2->GetMaximum()) * 3);
  q1->SetMinimum(50);
  q1->Draw("HIST");
  q2->Draw("HIST SAME");
  const double ratio = q2->Integral() / q1->Integral();
  TLegend *L2 = new TLegend(0.35, 0.72, 0.89, 0.89);
  L2->SetBorderSize(0);
  L2->SetFillStyle(0);
  L2->AddEntry(q1, Form("REAL (per frame, %d)", (int) nR), "l");
  L2->AddEntry(q2, Form("%s (per frame, %d): flat #times%.2f", tag, (int) nS, ratio), "l");
  L2->Draw();

  c.SaveAs(Form("/home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/residuals_evidence%s.png", suffix));
  printf("islands/frame: real %.0f | sim %.0f (x%.2f); <phisize> real %.2f | sim %.2f (%+.1f%%)\n",
         q1->Integral(), q2->Integral(), ratio, m1, m2, 100. * (m2 - m1) / m1);
}
