// funny_shapes.C — ORIGINAL-PURPOSE producer for sim_validation_plots/funny_shapes.png
//
// Reproduces the user's notebook view (Visualizing: "Event 74 layer 15 hits:
// zoomed phi-tbin view" — phi in [0,1] rad, tbin in [600,800], viridis, ADC color)
// and places the SIM counterpart beside it, to check whether the simulation
// reproduces the asymmetric cluster shapes seen in real data.
//
//   LEFT : REAL event 74, layer 15 (canonical cuts layer==15 && adc>0)
//   RIGHT: SIM v4.0 (exam6 library through P0-P3), layer 15 — the frame
//          whose hit count inside the window is closest to the real event's
//          (matched-activity choice, printed).
// Identical axes, binning (pad-pitch-scale phi bins x 1-tbin bins), palette and
// color scale (0..850, as in the notebook screenshot).

#include <TCanvas.h>
#include <TFile.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <cmath>
#include <cstdio>

void funny_shapes(const char *realpix =
                      "/home/rog/sPHENIX/3D_ClusterFindingML/clusters_seeds_island_79507-0.root_ntuplizer.root",
                  const char *simpix = "digi_frames_production_v40b.root",
                  int realevent = 74, int simlayer = 15, int simframe = -1, int reallayer = -1,
                  double philo = 0.0, double phihi = 1.0, int tlo = 600, int thi = 800,
                  const char *out = "/home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/funny_shapes.png",
                  const char *tag = "pAu v3.6")
{
  gROOT->SetBatch(1);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);
  if (reallayer < 0)
  {
    reallayer = simlayer;
  }
  const int NPHI = 180;  // ~ pad pitch of layer 15 (1128 pads / 2pi) over 1 rad

  // ---- real panel ----
  TFile *fr = TFile::Open(realpix);
  TTree *rh = (TTree *) fr->Get("ntp_hit");
  TH2D *hr = new TH2D("hr", Form("REAL event %d layer %d: zoomed #phi-tbin view;#phi [rad];tbin", realevent, reallayer),
                      NPHI, philo, phihi, thi - tlo, tlo, thi);
  rh->Draw(Form("tbin:phi>>hr"),
           Form("event==%d&&layer==%d&&adc>0&&phi>=%f&&phi<%f&&tbin>=%d&&tbin<%d",
                realevent, reallayer, philo, phihi, tlo, thi),
           "goff COLZ");
  // ADC-weighted fill (Draw above counts; refill weighted)
  hr->Reset();
  rh->Draw(Form("tbin:phi>>hr"),
           Form("adc*(event==%d&&layer==%d&&adc>0&&phi>=%f&&phi<%f&&tbin>=%d&&tbin<%d)",
                realevent, reallayer, philo, phihi, tlo, thi),
           "goff");
  double nreal = hr->GetEntries();

  // ---- sim panel: matched-activity frame ----
  TFile *fs = TFile::Open(simpix);
  TTree *sh = (TTree *) fs->Get("ntp_hit");
  int bestfr = -1;
  double bestd = 1e18, realcount = 0;
  {
    TH1D hc("hc", "", 1, 0, 1);
    rh->Draw("0.5>>hc", Form("event==%d&&layer==%d&&adc>0&&phi>=%f&&phi<%f&&tbin>=%d&&tbin<%d",
                             realevent, reallayer, philo, phihi, tlo, thi), "goff");
    realcount = hc.GetBinContent(1);
    if (simframe >= 0)
    {
      bestfr = simframe;
      bestd = -1;
    }
    else
    for (int f = 0; f < 20; ++f)
    {
      TH1D hs("hs", "", 1, 0, 1);
      // sim time axis: tbin (digi zbin==tbin exactly; hit69 exports have zbin=NaN)
      sh->Draw("0.5>>hs", Form("event==%d&&layer==%d&&phi>=%f&&phi<%f&&tbin>=%d&&tbin<%d",
                               f, simlayer, philo, phihi, tlo, thi), "goff");
      double d = std::fabs(hs.GetBinContent(1) - realcount);
      if (d < bestd)
      {
        bestd = d;
        bestfr = f;
      }
    }
  }
  printf("funny_shapes: real window hits %.0f; matched sim frame %d (delta %.0f)\n", realcount, bestfr, bestd);
  TH2D *hs2 = new TH2D("hs2", Form("SIM %s (event %d) layer %d: zoomed #phi-tbin view;#phi [rad];tbin", tag, bestfr, simlayer),
                       NPHI, philo, phihi, thi - tlo, tlo, thi);
  sh->Draw(Form("tbin:phi>>hs2"),
           Form("adc*(event==%d&&layer==%d&&phi>=%f&&phi<%f&&tbin>=%d&&tbin<%d)",
                bestfr, simlayer, philo, phihi, tlo, thi),
           "goff");

  // ---- draw, identical color scale like the notebook (0..850) ----
  TCanvas c("c", "", 1700, 750);
  c.Divide(2, 1, 0.005, 0.01);
  for (auto pr : {std::pair<int, TH2D *>{1, hr}, {2, hs2}})
  {
    c.cd(pr.first);
    gPad->SetRightMargin(0.14);
    pr.second->SetMinimum(0);
    pr.second->SetMaximum(850);
    pr.second->GetZaxis()->SetTitle("ADC");
    pr.second->Draw("COLZ");
  }
  c.SaveAs(out);
  printf("funny_shapes: original-purpose view regenerated (real evt %d vs %s frame %d)\n", realevent, tag, bestfr);
}
