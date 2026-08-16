// trunc_check.C — the run-79507 recording-truncation problem, printed and plotted.
// 38 of the 100 events stop recording before the end of the 965-tbin window;
// 62 run to completion (endpoint ~957-970). Discovered 2026-07-22 (endpoint
// census), verified independently by the side review; the COMPLETE-62 subset is
// what targets_v51.txt is measured on.
//
// Method (no truth, no sim): per event, the recording endpoint is the tbin below
// which 99.9% of that event's TPC hits lie. A complete event's endpoint sits at
// the window edge; a truncated one's sits wherever the DAQ stopped.
//
// usage:  root -l -b -q 'trunc_check.C+()'
//         root -l -b -q 'trunc_check.C+("<other file>", 100)'
// out:    stdout table + sim_validation_plots/trunc_check.png
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <cstdio>
#include <vector>

void trunc_check(const char *realf =
                     "/home/rog/sPHENIX/3D_ClusterFindingML/clusters_seeds_island_79507-0.root_ntuplizer.root",
                 int nev = 100, int ntb = 965, int fullcut = 900)
{
  gROOT->SetBatch(1);
  gStyle->SetOptStat(0);
  TFile *f = TFile::Open(realf);
  TTree *t = (TTree *) f->Get("ntp_hit");
  const char *CUT = "layer>=7&&layer<=54&&adc>0";  // canon.h TPC_CUT

  // ---- one pass: (event, tbin) map at full tbin resolution ----
  TH2F *m = new TH2F("m", "", nev, 0.5, nev + 0.5, ntb, 0, ntb);
  t->Draw("tbin:event>>m", CUT, "goff");

  // ---- per-event endpoint = tbin containing the 99.9% quantile ----
  std::vector<int> ep(nev + 1, 0);
  std::vector<double> tot(nev + 1, 0);
  std::vector<int> trunc;
  int nc = 0;
  printf("\n  event   endpoint   TPC hits\n  ------------------------------\n");
  for (int e = 1; e <= nev; ++e)
  {
    for (int b = 1; b <= ntb; ++b) tot[e] += m->GetBinContent(e, b);
    if (tot[e] < 1) { printf("  %5d      --         empty\n", e); continue; }
    double c = 0;
    for (int b = ntb; b >= 1; --b)
    {
      c += m->GetBinContent(e, b);
      if (c > 0.001 * tot[e]) { ep[e] = (int) m->GetYaxis()->GetBinCenter(b); break; }
    }
    bool full = ep[e] >= fullcut;
    if (full) nc++; else trunc.push_back(e);
    printf("  %5d   %6d   %10.0f %s\n", e, ep[e], tot[e], full ? "" : "  <== TRUNCATED");
  }
  printf("\n  COMPLETE %d   TRUNCATED %zu   (endpoint cut %d)\n", nc, trunc.size(), fullcut);
  printf("  truncated events: ");
  for (int e : trunc) printf("%d ", e);
  printf("\n\n");

  // ---- cut strings for the two classes ----
  TString trcut, cocut;
  for (int e = 1; e <= nev; ++e)
  {
    if (tot[e] < 1) continue;
    TString &s = (ep[e] >= fullcut) ? cocut : trcut;
    s += Form("%sevent==%d", s.Length() ? "||" : "", e);
  }

  // ---- figure ----
  TCanvas c("c", "", 1800, 1200);
  c.Divide(2, 2);

  // [1] recording-extent map: the one-glance view
  c.cd(1);
  TH2F *hm = new TH2F("hm", "recording extent per event;event;tbin", nev, 0.5, nev + 0.5, 193, 0, ntb);
  t->Draw("tbin:event>>hm", CUT, "goff");
  hm->Draw("COLZ");
  gPad->SetLogz();
  gPad->SetRightMargin(0.13);

  // [2] class-average arrival curves
  c.cd(2);
  TH1D *hco = new TH1D("hco", "arrival curve, class average;tbin;hits / event / 5 tbins", 193, 0, ntb);
  TH1D *htr = new TH1D("htr", "", 193, 0, ntb);
  t->Draw("tbin>>hco", Form("(%s)&&(%s)", CUT, cocut.Data()), "goff");
  t->Draw("tbin>>htr", Form("(%s)&&(%s)", CUT, trcut.Data()), "goff");
  hco->Scale(1. / nc);
  htr->Scale(1. / trunc.size());
  hco->SetLineColor(kBlue + 1);
  htr->SetLineColor(kRed + 1);
  hco->SetLineWidth(2);
  htr->SetLineWidth(2);
  hco->SetMinimum(0);
  hco->Draw("HIST");
  htr->Draw("HIST SAME");
  TLegend *L = new TLegend(0.45, 0.72, 0.89, 0.88);
  L->SetBorderSize(0);
  L->SetFillStyle(0);
  L->AddEntry(hco, Form("COMPLETE (%d events)", nc), "l");
  L->AddEntry(htr, Form("TRUNCATED (%zu events)", trunc.size()), "l");
  L->Draw();

  // [3] endpoint distribution
  c.cd(3);
  TH1D *he = new TH1D("he", "recording endpoint;endpoint tbin;events", 40, 0, ntb + 35);
  for (int e = 1; e <= nev; ++e)
    if (tot[e] >= 1) he->Fill(ep[e]);
  he->SetLineColor(kBlack);
  he->SetFillColorAlpha(kAzure - 9, 0.6);
  he->SetLineWidth(2);
  he->Draw("HIST");
  TLatex tx;
  tx.SetNDC();
  tx.SetTextSize(0.033);
  tx.DrawLatex(0.16, 0.82, Form("#color[600]{%d complete} pile up at the window edge", nc));
  tx.DrawLatex(0.16, 0.77, Form("#color[632]{%zu truncated} spread ~uniformly", trunc.size()));

  // [4] survival curve: fraction of events still recording at tbin t
  c.cd(4);
  TH1D *hs = new TH1D("hs", "fraction of events still recording;tbin;S(t)", 193, 0, ntb);
  int ngood = nc + (int) trunc.size();
  for (int b = 1; b <= 193; ++b)
  {
    double tb = hs->GetBinLowEdge(b);
    int alive = 0;
    for (int e = 1; e <= nev; ++e)
      if (tot[e] >= 1 && ep[e] >= tb) alive++;
    hs->SetBinContent(b, (double) alive / ngood);
  }
  hs->SetMinimum(0);
  hs->SetMaximum(1.05);
  hs->SetLineColor(kGreen + 3);
  hs->SetLineWidth(2);
  hs->Draw("HIST");
  tx.DrawLatex(0.18, 0.34, "a linear fall = stops distributed uniformly in time");
  tx.DrawLatex(0.18, 0.28, Form("S(t) plateau = %.2f = the COMPLETE fraction", (double) nc / ngood));
  tx.DrawLatex(0.18, 0.22, "the final cliff is the 965-tbin window edge itself");

  const char *out = "/home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/trunc_check.png";
  c.SaveAs(out);
  printf("  wrote %s\n", out);
}
