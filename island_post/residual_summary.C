// residual_summary.C — ONE-PLOT residual overview: every measured sim-vs-real
// deviation as a horizontal bar row, grouped by family, with BOTH references
// where they exist (filled blue = vs ALL real events; open red = vs the
// STEADY subset excluding end-truncated boundary windows). Shaded band =
// +-5% window-gate noise floor. Reads a residuals_<ver>.txt data file
// (group|label|devAll|devSteady, 999 = n/a) so each era regenerates by
// swapping the data file. Values beyond the axis are clipped with printed
// arrows + numbers.
#include <TCanvas.h>
#include <TBox.h>
#include <TLine.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TROOT.h>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

void residual_summary(const char *datafile = "residuals_v50.txt",
                      const char *title = "v5.0 pp pilot: all residuals (sim/real - 1)",
                      const char *out = "/home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/residual_summary_v50.png")
{
  gROOT->SetBatch(1);
  gStyle->SetOptStat(0);
  struct Row { std::string grp, lab; double a, s; };
  std::vector<Row> rows;
  std::ifstream fi(datafile);
  std::string L;
  while (std::getline(fi, L))
  {
    if (L.empty() || L[0] == '#') continue;
    std::stringstream ss(L);
    std::string g, l, a, s;
    std::getline(ss, g, '|'); std::getline(ss, l, '|');
    std::getline(ss, a, '|'); std::getline(ss, s, '|');
    rows.push_back({g, l, atof(a.c_str()), atof(s.c_str())});
  }
  const int N = rows.size();
  const double XMIN = -32, XMAX = 32;
  TCanvas c("c", "", 1150, 60 + 42 * N);
  c.SetLeftMargin(0.30); c.SetRightMargin(0.10);
  c.SetTopMargin(0.06); c.SetBottomMargin(0.07);
  TH2D fr("fr", Form("%s;deviation [%%];", title), 10, XMIN, XMAX, N, 0, N);
  fr.GetYaxis()->SetLabelSize(0);
  fr.GetYaxis()->SetTickLength(0);
  fr.Draw();
  TBox band(-5, 0, 5, N);
  band.SetFillColorAlpha(kGray, 0.30); band.Draw();
  TLine zero(0, 0, 0, N);
  zero.SetLineColor(kBlack); zero.SetLineWidth(2); zero.Draw();
  TLatex tx; tx.SetTextSize(0.014 * 22. / N * 1.5 > 0.03 ? 0.03 : 0.0145 * 22. / N + 0.008);
  tx.SetTextAlign(32);
  TLatex gx; gx.SetTextAlign(12); gx.SetTextFont(62);
  gx.SetTextSize(tx.GetTextSize());
  std::string lastg;
  for (int i = 0; i < N; ++i)
  {
    int y = N - 1 - i;   // first row on top
    Row &r = rows[i];
    if (r.grp != lastg)
    {
      lastg = r.grp;
      TLine *sep = new TLine(XMIN, y + 1, XMAX, y + 1);
      sep->SetLineStyle(3); sep->SetLineColor(kGray + 1);
      if (i) sep->Draw();
      gx.SetTextColor(kGray + 2);
      gx.DrawLatex(XMIN + 0.7, y + 0.62, r.grp.c_str());
    }
    tx.DrawLatex(XMIN - 0.8, y + 0.5, r.lab.c_str());
    auto mark = [&](double v, int col, int sty, double dy) {
      if (v > 900) return;
      double vc = std::max(XMIN + 1.2, std::min(XMAX - 1.2, v));
      TLine *bar = new TLine(0, y + 0.5 + dy, vc, y + 0.5 + dy);
      bar->SetLineColor(col); bar->SetLineWidth(3); bar->Draw();
      TGraph *g = new TGraph(1);
      g->SetPoint(0, vc, y + 0.5 + dy);
      g->SetMarkerColor(col); g->SetMarkerStyle(sty); g->SetMarkerSize(1.3);
      g->Draw("P SAME");
      if (v != vc)
      {
        TLatex *ov = new TLatex(vc + (v > 0 ? -0.4 : 0.4), y + 0.5 + dy + 0.16,
                                Form("%+.0f%%", v));
        ov->SetTextSize(tx.GetTextSize() * 0.85);
        ov->SetTextColor(col);
        ov->SetTextAlign(v > 0 ? 32 : 12);
        ov->Draw();
      }
    };
    mark(r.a, kBlue + 1, 20, +0.12);   // vs ALL real
    mark(r.s, kRed + 1, 24, -0.14);    // vs STEADY subset
  }
  TLegend *lg = new TLegend(0.315, 0.075, 0.52, 0.145);
  lg->SetBorderSize(0); lg->SetFillStyle(0);
  TGraph *ga = new TGraph(1); ga->SetMarkerColor(kBlue + 1); ga->SetMarkerStyle(20);
  TGraph *gs = new TGraph(1); gs->SetMarkerColor(kRed + 1); gs->SetMarkerStyle(24);
  lg->AddEntry(ga, "vs ALL real events as recorded (99, laser event vetoed)", "p");
  lg->AddEntry(gs, "vs COMPLETE subset (61 non-laser events, full windows)", "p");
  lg->Draw();
  c.SaveAs(out);
  printf("residual_summary: %d rows -> %s\n", N, out);
}
