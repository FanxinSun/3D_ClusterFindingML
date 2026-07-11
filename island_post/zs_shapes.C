// zs_shapes.C — ORIGINAL-PURPOSE producer for sim_validation_plots/zs_shapes.png
//
// Day-2 origin (then an inline throwaway; persistent now per the provenance rule):
// ZS-region per-pixel ADC SHAPE overlay — REAL vs our readout configs — integer-ADC
// bins, every curve normalized over the ADC 8-40 window, log-y, x shown 0-60.
// It is the figure that exposed the quantized-pedestal binning comb and drove the
// A/B/B3 shoot-out; the ZS scheme itself (T1=20/T2=11/p=0.39/p2=5e-4) is frozen since B3.
//
// v33_b42 edition: REAL anchor vs the v3.3 (B4.2) production, with v3.2 (B4.1,
// no saturation-triggered slow disturbance) kept as the reference curve.
#include "canon.h"

void zs_shapes(const char *vcur = "digi_frames_production_v34.root",
               const char *vref = "digi_frames_production_v33.root")
{
  gROOT->SetBatch(1);
  gStyle->SetOptStat(0);

  // REAL: consume the frozen anchors (canonical cut baked in). r_adc = full spectrum
  // (denominator for leakage fractions), r_adcz = near-threshold zoom (the drawn shape).
  TH1D *ra = CANON::anchor("r_adc");
  TH1D *hr = CANON::anchor("r_adcz");
  if (!ra || !hr) return;

  // SIM: full-spectrum integer-bin fills (identical bin width 1 -> shape-comparable).
  TFile *f1 = TFile::Open(vcur);
  TTree *t1 = (TTree *) f1->Get("ntp_hit");
  TH1D *h1 = new TH1D("h_v33", "", ra->GetNbinsX(), ra->GetXaxis()->GetXmin(), ra->GetXaxis()->GetXmax());
  t1->Draw("adc>>h_v33", "", "goff");
  TFile *f2 = TFile::Open(vref);
  TTree *t2 = (TTree *) f2->Get("ntp_hit");
  TH1D *h2 = new TH1D("h_vref", "", ra->GetNbinsX(), ra->GetXaxis()->GetXmin(), ra->GetXaxis()->GetXmax());
  t2->Draw("adc>>h_vref", "", "goff");

  // numbers BEFORE normalization (fractions of ALL kept pixels)
  auto frac = [](TH1D *h, double a, double b) {
    return h->Integral(h->FindBin(a), h->FindBin(b)) / h->Integral();
  };
  printf("sub-10 (1-10) fraction   : real %.2e | v3.3 %.2e | v3.2 %.2e\n",
         frac(ra, 1, 10), frac(h1, 1, 10), frac(h2, 1, 10));
  printf("shoulder(10-19) fraction : real %.3f | v3.3 %.3f | v3.2 %.3f\n",
         frac(ra, 10, 19), frac(h1, 10, 19), frac(h2, 10, 19));
  printf("hi(20-40) fraction       : real %.3f | v3.3 %.3f | v3.2 %.3f\n",
         frac(ra, 20, 40), frac(h1, 20, 40), frac(h2, 20, 40));
  printf("shoulder/hi ratio        : real %.3f | v3.3 %.3f | v3.2 %.3f\n",
         frac(ra, 10, 19) / frac(ra, 20, 40), frac(h1, 10, 19) / frac(h1, 20, 40),
         frac(h2, 10, 19) / frac(h2, 20, 40));

  // shape normalization over the 8-40 window (frozen convention of this figure)
  auto nrm = [](TH1D *h) {
    double s = h->Integral(h->FindBin(8), h->FindBin(40));
    if (s > 0) h->Scale(1. / s);
  };
  nrm(hr);
  nrm(h1);
  nrm(h2);

  hr->SetLineColor(kBlue + 1);
  h1->SetLineColor(kMagenta + 1);
  h2->SetLineColor(kGreen + 2);
  for (auto h : {hr, h1, h2})
  {
    h->SetLineWidth(2);
    h->SetStats(0);
  }
  TCanvas c("c", "", 850, 600);
  gPad->SetLogy();
  hr->SetTitle("ZS-region shape (integer-ADC bins);per-pixel ADC;norm (8-40)");
  hr->GetXaxis()->SetRangeUser(0, 60);
  hr->SetMaximum(std::max({hr->GetMaximum(), h1->GetMaximum(), h2->GetMaximum()}) * 2);
  hr->SetMinimum(2e-6);
  hr->Draw("HIST");
  h1->Draw("HIST SAME");
  h2->Draw("HIST SAME");
  // legend in the empty lower-right (curves occupy the upper band at all x)
  TLegend L(0.45, 0.16, 0.89, 0.38);
  L.SetBorderSize(0);
  L.SetFillStyle(0);
  L.AddEntry(hr, "REAL", "l");
  L.AddEntry(h1, "SIM v3.4 (B4.3: real region ZS)", "l");
  L.AddEntry(h2, "SIM v3.3 (B4.2, reference)", "l");
  L.Draw();
  c.SaveAs("/home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/zs_shapes.png");
  printf("saved zs_shapes.png (v3.4 vs v3.3 vs real)\n");
}
