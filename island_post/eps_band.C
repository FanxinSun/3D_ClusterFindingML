// eps_band.C — the eps_MBD scan band figure (2026-07-24). Reads
// eps_scan_results.txt: fit-family meters (solid, colored) with linear fits
// and zero-crossings; Class-A content family (gray dashed, never fit);
// <size> (flat = response-owned). Vertical bands: fitted eps 0.588+-0.033
// (green) and the geometric proxy 0.519 (line).
#include <TCanvas.h>
#include <TGraph.h>
#include <TF1.h>
#include <TH2D.h>
#include <TBox.h>
#include <TLine.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TROOT.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

void eps_band()
{
  gROOT->SetBatch(1);
  gStyle->SetOptStat(0);
  std::ifstream fi("eps_scan_results.txt");
  std::vector<double> e;
  std::vector<std::vector<double>> col(10);
  std::string L;
  while (std::getline(fi, L))
  {
    if (L.empty() || L[0] == '#') continue;
    for (auto &c : L) if (c == '|') c = ' ';
    std::stringstream ss(L);
    double v; ss >> v; e.push_back(v);
    for (int k = 0; k < 10; ++k) { ss >> v; col[k].push_back(v); }
  }
  const int N = e.size();
  // columns: 0 w1, 1 w2, 2 w3, 3 bump, 4 step, 5 CV, 6 px, 7 lowadc, 8 islands, 9 size
  const char *nm[10] = {"w1", "w2", "w3", "bump", "step", "CV",
                        "px/event", "low-adc", "islands", "island <size>"};
  int colr[10] = {kBlue + 1, kGray + 1, kAzure + 7, kViolet - 5, kRed + 1, kOrange + 7,
                  kGray + 2, kGray + 1, kGray + 3, kGreen + 3};
  TCanvas c("c", "", 1050, 760);
  c.SetLeftMargin(0.09); c.SetRightMargin(0.03); c.SetTopMargin(0.07);
  TH2D fr("fr", "#varepsilon_{MBD} scan: residuals vs COMPLETE-62;#varepsilon_{MBD};deviation [%]",
          10, 0.43, 0.67, 10, -34, 14);
  fr.Draw();
  TBox band(0.588 - 0.033, -34, 0.588 + 0.033, 14);
  band.SetFillColorAlpha(kGreen - 8, 0.35); band.Draw();
  TBox noise(0.43, -3, 0.67, 3);
  noise.SetFillColorAlpha(kGray, 0.25); noise.Draw();
  TLine z(0.43, 0, 0.67, 0); z.SetLineColor(kBlack); z.Draw();
  TLine px(0.519, -34, 0.519, 14); px.SetLineStyle(2); px.SetLineColor(kGray + 2); px.Draw();
  TLegend Lg(0.115, 0.12, 0.47, 0.42);
  Lg.SetBorderSize(0); Lg.SetFillStyle(0); Lg.SetTextSize(0.028);
  for (int k = 0; k < 10; ++k)
  {
    if (k == 1) continue;  // w2: trendless, omitted for legibility
    TGraph *g = new TGraph(N);
    for (int i = 0; i < N; ++i) g->SetPoint(i, e[i], col[k][i]);
    bool fitfam = (k == 0 || k == 2 || k == 3 || k == 4 || k == 5);
    g->SetLineColor(colr[k]); g->SetMarkerColor(colr[k]);
    g->SetLineWidth(fitfam ? 2 : 2);
    g->SetLineStyle(fitfam ? 1 : 7);
    g->SetMarkerStyle(fitfam ? 20 : 24);
    g->SetMarkerSize(1.0);
    g->Draw("PL SAME");
    Lg.AddEntry(g, Form("%s%s", nm[k], fitfam ? "" : (k == 9 ? "  [flat: response]" : "  [Class A: report only]")), "pl");
  }
  Lg.Draw();
  TLatex tx; tx.SetTextSize(0.030);
  tx.DrawLatex(0.556, 11.6, "#varepsilon = 0.588 #pm 0.033");
  tx.SetTextSize(0.026); tx.SetTextColor(kGray + 2);
  tx.DrawLatex(0.4955, 11.6, "proxy 0.519");
  c.SaveAs("/home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/eps_scan_band_v51.png");
  printf("saved eps_scan_band_v51.png\n");
}
