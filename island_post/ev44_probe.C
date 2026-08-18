// ev44_probe.C — what the bright dot at (event 44, tbin ~330) in the recording-extent map is.
//
// ANSWER (2026-08-17 probe): event 44 is the ONE laser-flash event of the 100 (a GL1
// laser-triggered readout: the sPHENIX diffuse laser fires ~4.5 us after the trigger; the
// flash of photoelectrons off the central-membrane stripes arrives one full drift later at
// tbin 329). 129.8k hits in [322,340]; the OTHER 99 events show NO excess at tbin 327-331
// (their sum is flat to <1%) -> the tbin-330 spike of the all-100 arrival curve is event 44
// alone; laser_frames.txt (16-tbin window minus a later control) was flagging arrival-curve
// slopes, not flashes. It is NOT a tbin=0 / trigger-start artifact (no tbin=0 spike in any
// event) — but the reco knows it: coresoftware LaserEventIdentifier (peak>1000 & peak/mean>=7
// over samples>=320) tags it and TpcClusterizer skips it -> ntp_info ntrk=ntpcseed=nclustpc=0
// for event 44 while nhittpcall=535,606 (largest of the run).
// Collateral inside event 44 (all measured here): prompt detector-wide spike at tbin 86-89 =
// the laser fire time (one full drift, 242 tbins, before the CM flash); 3,707 pads saturate
// (adc 940-963); the saturated pads keep firing over tbin 340-400 (0.2 hits/pad/tbin, 100x
// baseline) = the ~40k-hit tail; 15.8k adc==0 hits (43% of the file's) sit exactly at the
// flash peak (a burst-only unpacker artifact, dropped by islandize).
//
// usage: cd island_post && root -l -b -q 'ev44_probe.C+()'
// out:   sim_validation_plots/ev44_probe.png  (+ stdout tables)
#include "canon.h"
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2F.h>
#include <THStack.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TArrow.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <TGaxis.h>
#include <cstdio>
#include <cmath>
#include <map>
#include <set>
#include <vector>
#include <algorithm>

void ev44_probe(const char *realf =
                    "/home/rog/sPHENIX/3D_ClusterFindingML/clusters_seeds_island_79507-0.root_ntuplizer.root",
                int LASER = 44, int CONTRAST = 18, int nev = 100, int ntb = 965)
{
  gROOT->SetBatch(1);
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.05);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  TFile *f = TFile::Open(realf);
  TTree *t = (TTree *) f->Get("ntp_hit");
  const char *CUT = CANON::TPC_CUT;

  // ---------- pass 1: (event, tbin) map, integer-centred ----------
  TH2F *m = new TH2F("m", "", nev, 0.5, nev + 0.5, ntb, -0.5, ntb - 0.5);
  t->Draw("tbin:event>>m", CUT, "goff");

  // per-event narrow-spike excess at [328,331] vs immediate neighbours [324,327]+[332,335]
  TH1D *hex = new TH1D("hex", "narrow excess at tbin 328-331 per event;event;hits[328-331] - neighbours (log)", nev, 0.5, nev + 0.5);
  TH1D *hexf = (TH1D *) hex->Clone("hexf");  // the 43 events laser_frames.txt flagged (besides 44)
  TH1D *hex44 = (TH1D *) hex->Clone("hex44");  // the laser event itself
  int flagged[] = {2, 4, 10, 12, 13, 14, 17, 20, 22, 23, 25, 30, 31, 35, 36, 41, 42, 47, 49, 54, 56, 57, 58, 59, 60, 61, 63, 64, 66, 68, 69, 70, 76, 80, 82, 83, 87, 88, 90, 92, 94, 97, 100};
  std::set<int> F(flagged, flagged + 43);
  printf("\n  narrow-spike excess S[328,331]-N per event (N = mean of [324,327],[332,335]):\n");
  for (int e = 1; e <= nev; ++e)
  {
    double S = 0, N = 0;
    for (int tb = 328; tb <= 331; ++tb) S += m->GetBinContent(e, tb + 1);
    for (int tb = 324; tb <= 327; ++tb) N += m->GetBinContent(e, tb + 1);
    for (int tb = 332; tb <= 335; ++tb) N += m->GetBinContent(e, tb + 1);
    N /= 2;
    double x = S - N;
    (e == LASER ? hex44 : F.count(e) ? hexf : hex)->SetBinContent(e, std::max(x, 0.5));
    if (e == LASER || std::fabs(x) > 250) printf("    event %3d: %+8.0f\n", e, x);
  }
  // mean arrival curve of the other 99 events (per event) and event LASER's own
  TH1D *h44 = new TH1D("h44", "event 44 vs the mean of the other 99;tbin;TPC hits per tbin (log)", ntb, -0.5, ntb - 0.5);
  TH1D *hmean = new TH1D("hmean", "", ntb, -0.5, ntb - 0.5);
  for (int tb = 0; tb < ntb; ++tb)
  {
    h44->SetBinContent(tb + 1, m->GetBinContent(LASER, tb + 1));
    double s = 0;
    for (int e = 1; e <= nev; ++e)
      if (e != LASER) s += m->GetBinContent(e, tb + 1);
    hmean->SetBinContent(tb + 1, s / (nev - 1));
  }
  // sum of the 99 others around the flash: the "no residual" statement
  double s99w = 0, s99n = 0;
  for (int tb = 327; tb <= 331; ++tb) s99w += hmean->GetBinContent(tb + 1) * (nev - 1);
  for (int tb = 320; tb <= 324; ++tb) s99n += hmean->GetBinContent(tb + 1) * (nev - 1);
  for (int tb = 334; tb <= 338; ++tb) s99n += hmean->GetBinContent(tb + 1) * (nev - 1);
  printf("  sum of the other 99 events: [327,331] = %.0f vs neighbours %.0f -> residual %+.2f%%\n", s99w, s99n / 2, 100 * (s99w / (s99n / 2) - 1));

  // ---------- pass 2: per-pad loop over event LASER (+ side x sector for LASER and CONTRAST) ----------
  float event, layer, phielem, zelem, phibin, tbin, adc;
  t->SetBranchStatus("*", 0);
  for (const char *b : {"event", "layer", "phielem", "zelem", "phibin", "tbin", "adc"}) t->SetBranchStatus(b, 1);
  t->SetBranchAddress("event", &event);
  t->SetBranchAddress("layer", &layer);
  t->SetBranchAddress("phielem", &phielem);
  t->SetBranchAddress("zelem", &zelem);
  t->SetBranchAddress("phibin", &phibin);
  t->SetBranchAddress("tbin", &tbin);
  t->SetBranchAddress("adc", &adc);
  struct Hit { int tb; float adc; };
  std::map<long, std::vector<Hit>> pads;
  double SS[2][24] = {{0}}, NN[2][24] = {{0}};  // [0]=LASER, [1]=CONTRAST ; 24 = side*12+sector
  Long64_t N = t->GetEntries();
  for (Long64_t i = 0; i < N; ++i)
  {
    t->GetEntry(i);
    int e = (int) event;
    if (e != LASER && e != CONTRAST) continue;
    if (layer < 7 || layer > 54) continue;
    int tb = (int) tbin, cell = (int) zelem * 12 + (int) phielem, k = (e == LASER) ? 0 : 1;
    if (adc > 0)
    {
      if (tb >= 328 && tb <= 331) SS[k][cell]++;
      else if ((tb >= 324 && tb <= 327) || (tb >= 332 && tb <= 335)) NN[k][cell] += 0.5;
    }
    if (e == LASER)
      pads[((long) layer * 2 + (long) zelem) * 4096 + (long) phibin].push_back({tb, adc});
  }
  // pad classes by flash amplitude A = max adc in [327,332]
  TH1D *hzoom[4];
  const char *cn[4] = {"no flash hit / soft (A<500)", "hard (500<=A<930)", "SATURATED (A>=930)", "adc==0 hits (all pads)"};
  int col[4] = {kAzure - 9, kOrange - 3, kRed + 1, kBlack};
  for (int c = 0; c < 4; ++c) hzoom[c] = new TH1D(Form("hz%d", c), "", 130, 299.5, 429.5);
  int npad[3] = {0, 0, 0}; double tail[3] = {0, 0, 0}, sat_all = 0;
  for (auto &kv : pads)
  {
    auto &v = kv.second;
    float A = 0;
    for (auto &h : v) if (h.tb >= 327 && h.tb <= 332) A = std::max(A, h.adc);
    int c = (A < 500) ? 0 : (A < 930) ? 1 : 2;
    npad[c]++;
    for (auto &h : v)
    {
      if (h.adc == 0) { hzoom[3]->Fill(h.tb); continue; }
      hzoom[c]->Fill(h.tb);
      if (h.tb >= 340 && h.tb < 400) tail[c]++;
      if (h.adc >= 930 && h.tb >= 322 && h.tb <= 340) sat_all++;
    }
  }
  printf("  pads by flash amplitude: soft/none %d, hard %d, saturated %d | tail [340,400) hits/pad/tbin: %.4f %.4f %.4f | saturated flash hits %.0f\n",
         npad[0], npad[1], npad[2], tail[0] / npad[0] / 60, tail[1] / npad[1] / 60, tail[2] / npad[2] / 60, sat_all);
  // side x sector maps: fraction of the event's positive narrow excess carried by each cell (%)
  TH2D *hmap = new TH2D("hmap", Form("where the narrow tbin-328-331 excess sits: event %d vs event %d;sector;", LASER, CONTRAST), 12, -0.5, 11.5, 4, 0, 4);
  const char *rl[4] = {Form("ev%d side 0", CONTRAST), Form("ev%d side 1", CONTRAST), Form("ev%d side 0", LASER), Form("ev%d side 1", LASER)};
  for (int k = 0; k < 2; ++k)
  {
    double pos = 0; int n3 = 0;
    for (int c = 0; c < 24; ++c) { double x = SS[k][c] - NN[k][c]; if (x > 0) pos += x; if (x > 3 * sqrt(std::max(NN[k][c], 1.))) n3++; }
    for (int c = 0; c < 24; ++c)
    {
      double x = SS[k][c] - NN[k][c];
      int row = (k == 0 ? 2 : 0) + c / 12;  // LASER rows 2,3 ; CONTRAST rows 0,1
      hmap->SetBinContent(c % 12 + 1, row + 1, 100 * std::max(x, 0.) / std::max(pos, 1.));
    }
    printf("  event %d: positive narrow excess %.0f hits, cells(>3 sigma) %d/24\n", k == 0 ? LASER : CONTRAST, pos, n3);
  }
  for (int r = 0; r < 4; ++r) hmap->GetYaxis()->SetBinLabel(r + 1, rl[r]);

  // ---------- figure ----------
  TCanvas c("c", "", 1800, 1300);
  c.Divide(2, 2);
  TLatex tx; tx.SetNDC(); tx.SetTextSize(0.034);

  // [1] arrival curves, log y
  c.cd(1); gPad->SetLogy(); gPad->SetLeftMargin(0.11); gPad->SetRightMargin(0.03);
  h44->SetLineColor(kRed + 1); h44->SetLineWidth(2);
  hmean->SetLineColor(kBlue + 1); hmean->SetLineWidth(2);
  h44->SetMinimum(20); h44->SetMaximum(2e5);
  h44->GetXaxis()->SetTitleSize(0.045); h44->GetYaxis()->SetTitleSize(0.045);
  h44->Draw("HIST"); hmean->Draw("HIST SAME");
  TArrow ar; ar.SetLineColor(kGray + 2); ar.SetLineWidth(2); ar.SetArrowSize(0.012);
  ar.DrawArrow(87, 3.0e4, 329, 3.0e4, 0.012, "<|>");
  tx.SetTextColor(kGray + 3);
  tx.DrawLatex(0.30, 0.86, "one full drift (242 tbins)");
  tx.SetTextColor(kRed + 1);
  tx.DrawLatex(0.14, 0.62, "#splitline{prompt spike tbin 86-89}{= laser fire time}");
  tx.DrawLatex(0.42, 0.75, "#splitline{CM flash tbin 329}{129.8k hits, 3.7k pads saturate}");
  tx.DrawLatex(0.53, 0.55, "#splitline{tail 340-400}{= saturated-channel recovery}");
  TLegend *L1 = new TLegend(0.58, 0.12, 0.94, 0.23); L1->SetBorderSize(0); L1->SetFillStyle(0); L1->SetTextSize(0.033);
  L1->AddEntry(h44, Form("event %d (%s)", LASER, "TPC, adc>0"), "l");
  L1->AddEntry(hmean, "mean of the other 99 events", "l");
  L1->Draw();

  // [2] narrow excess per event, log y
  c.cd(2); gPad->SetLogy(); gPad->SetLeftMargin(0.11); gPad->SetRightMargin(0.03);
  hex->SetFillColor(kBlue - 9); hex->SetLineColor(kBlue + 1);
  hexf->SetFillColor(kOrange - 3); hexf->SetLineColor(kOrange + 7);
  hex->SetMinimum(0.7); hex->SetMaximum(5e5);
  hex->GetXaxis()->SetTitleSize(0.045); hex->GetYaxis()->SetTitleSize(0.045);
  hex44->SetFillColor(kRed + 1); hex44->SetLineColor(kRed + 1);
  hex->Draw("HIST"); hexf->Draw("HIST SAME"); hex44->Draw("HIST SAME");
  TLine ln; ln.SetLineStyle(2); ln.SetLineColor(kGray + 2); ln.DrawLine(0.5, 1000, nev + 0.5, 1000);
  tx.SetTextColor(kGray + 3); tx.DrawLatex(0.56, 0.555, "1000 hits = LaserEventIdentifier floor");
  tx.SetTextColor(kBlack);
  tx.DrawLatex(0.14, 0.86, Form("event %d: 97.6k  |  every other event < 500 (localised track/looper structure)", LASER));
  TLegend *L2 = new TLegend(0.50, 0.64, 0.96, 0.80); L2->SetBorderSize(0); L2->SetFillStyle(0); L2->SetTextSize(0.031);
  L2->AddEntry(hex44, Form("event %d (flagged too)", LASER), "f");
  L2->AddEntry(hexf, "43 others flagged by laser_frames.txt", "f");
  L2->AddEntry(hex, "unflagged", "f");
  L2->Draw();

  // [3] zoom on event 44 by pad class (stacked) + adc==0
  c.cd(3); gPad->SetLogy(); gPad->SetLeftMargin(0.11); gPad->SetRightMargin(0.03);
  THStack *st = new THStack("st", Form("event %d, hits per tbin by pad class (stacked);tbin;hits per tbin (log)", LASER));
  for (int cidx = 0; cidx < 3; ++cidx) { hzoom[cidx]->SetFillColor(col[cidx]); hzoom[cidx]->SetLineColor(col[cidx]); st->Add(hzoom[cidx]); }
  st->SetMinimum(20); st->SetMaximum(2e5);
  st->Draw("HIST");
  st->GetXaxis()->SetTitleSize(0.045); st->GetYaxis()->SetTitleSize(0.045);
  hzoom[3]->SetLineColor(kBlack); hzoom[3]->SetLineWidth(2); hzoom[3]->SetLineStyle(1); hzoom[3]->Draw("HIST SAME");
  TLegend *L3 = new TLegend(0.42, 0.62, 0.96, 0.88); L3->SetBorderSize(0); L3->SetFillStyle(0); L3->SetTextSize(0.031);
  L3->SetHeader("pad class by its flash amplitude A = max adc in 327-332");
  for (int cidx = 0; cidx < 3; ++cidx) L3->AddEntry(hzoom[cidx], Form("%s: %d pads", cn[cidx], npad[cidx]), "f");
  L3->AddEntry(hzoom[3], "adc==0 hits (15.8k, all at 328-330)", "l");
  L3->Draw();
  tx.SetTextColor(kRed + 1);
  tx.DrawLatex(0.42, 0.55, Form("tail 340-400: %.3f hits/pad/tbin on saturated pads", tail[2] / npad[2] / 60));
  tx.DrawLatex(0.42, 0.50, Form("vs %.4f on the rest (x%.0f)", tail[0] / npad[0] / 60, (tail[2] / npad[2]) / (tail[0] / npad[0])));

  // [4] side x sector maps
  c.cd(4); gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.12); gPad->SetBottomMargin(0.12); gPad->SetTopMargin(0.16);
  hmap->GetZaxis()->SetTitle("% of the event's positive excess");
  hmap->GetXaxis()->SetTitleSize(0.045); hmap->GetYaxis()->SetLabelSize(0.045);
  hmap->SetMarkerSize(1.3);
  gStyle->SetPaintTextFormat(".0f");
  hmap->Draw("COLZ TEXT");
  ln.SetLineColor(kBlack); ln.SetLineStyle(1); ln.SetLineWidth(2); ln.DrawLine(-0.5, 2, 11.5, 2);
  tx.SetTextColor(kBlack);
  tx.SetTextSize(0.031);
  tx.DrawLatex(0.14, 0.865, Form("event %d: all 24 cells lit (detector-wide flash)  |  event %d: one cell (a looper)", LASER, CONTRAST));

  const char *out = "/home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/ev44_probe.png";
  c.SaveAs(out);
  printf("  wrote %s\n", out);
}
