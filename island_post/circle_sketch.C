// circle_sketch.C — geometry explainer for the truth-circle numbers:
// how pT=0.5 GeV + B=1.4 T + TPC radii (20..78 cm) chain into R=119 cm,
// chord 58.8 cm, sagitta 3.6 cm. Companion of truth_circle.C.
// Geometry: track circle of radius R through the vertex, center at (0,R);
// point after turning angle th sits at (R sin th, R(1-cos th)), at distance
// r = 2R sin(th/2) from the beamline.
#include <TCanvas.h>
#include <TH1F.h>
#include <TEllipse.h>
#include <TGraph.h>
#include <TLine.h>
#include <TArrow.h>
#include <TMarker.h>
#include <TLatex.h>
#include <TStyle.h>
#include <cmath>

void circle_sketch()
{
  const double R = 119.1;                    // = 0.5 GeV / (0.2998 * 1.4 T), in cm
  const double th_in = 2 * std::asin(20. / (2 * R));    // crosses r = 20 cm
  const double th_out = 2 * std::asin(78. / (2 * R));   // crosses r = 78 cm
  auto px = [&](double th) { return R * std::sin(th); };
  auto py = [&](double th) { return R * (1 - std::cos(th)); };

  gStyle->SetOptStat(0);
  TCanvas *cv = new TCanvas("cv", "circle geometry", 950, 720);
  TH1 *fr = gPad->DrawFrame(-8, -14, 100, 58,
                            "one chain: p_{T} #rightarrow R #rightarrow chord #rightarrow sagitta;x [cm];y [cm]");
  fr->GetYaxis()->SetTitleOffset(1.1);

  // TPC boundaries: circles AROUND THE BEAMLINE (position radius r)
  TEllipse *e1 = new TEllipse(0, 0, 20, 20);
  e1->SetFillStyle(0); e1->SetLineColor(kGray + 1); e1->SetLineWidth(2); e1->Draw();
  TEllipse *e2 = new TEllipse(0, 0, 78, 78);
  e2->SetFillStyle(0); e2->SetLineColor(kGray + 1); e2->SetLineWidth(2); e2->Draw();

  // the track's own circle (radius R, center NOT on the beamline): thin dashed
  const int N = 300;
  double xs[N], ys[N];
  for (int i = 0; i < N; ++i)
  {
    double th = -0.10 + 0.95 * i / (N - 1.);
    xs[i] = px(th); ys[i] = py(th);
  }
  TGraph *gall = new TGraph(N, xs, ys);
  gall->SetLineColor(kBlue - 9); gall->SetLineStyle(2); gall->Draw("L same");
  // the piece the TPC records (between r=20 and r=78): thick solid
  double xi[N], yi[N];
  for (int i = 0; i < N; ++i)
  {
    double th = th_in + (th_out - th_in) * i / (N - 1.);
    xi[i] = px(th); yi[i] = py(th);
  }
  TGraph *gin = new TGraph(N, xi, yi);
  gin->SetLineColor(kBlue + 1); gin->SetLineWidth(4); gin->Draw("L same");

  // chord between entry and exit
  TLine *ch = new TLine(px(th_in), py(th_in), px(th_out), py(th_out));
  ch->SetLineColor(kRed + 1); ch->SetLineStyle(7); ch->SetLineWidth(2); ch->Draw();

  // sagitta: chord midpoint -> arc midpoint
  double mx = 0.5 * (px(th_in) + px(th_out)), my = 0.5 * (py(th_in) + py(th_out));
  double ax = px(0.5 * (th_in + th_out)), ay = py(0.5 * (th_in + th_out));
  TArrow *sg = new TArrow(mx, my, ax, ay, 0.012, "<|>");
  sg->SetLineColor(kBlack); sg->SetLineWidth(2); sg->Draw();

  // radius pointer from arc midpoint toward the circle center (0, R) (off-plot)
  double dxn = (0 - ax) / R, dyn = (R - ay) / R, tlen = 40;
  TLine *rl = new TLine(ax, ay, ax + dxn * tlen, ay + dyn * tlen);
  rl->SetLineColor(kBlue + 1); rl->SetLineStyle(3); rl->Draw();

  TMarker *vx = new TMarker(0, 0, 29); vx->SetMarkerColor(kBlack); vx->SetMarkerSize(1.8); vx->Draw();

  TLatex t; t.SetTextSize(0.031);
  t.SetTextColor(kGray + 2);
  t.DrawLatex(3.5, 15.0, "TPC inner r = 20 cm");
  t.DrawLatex(63, 48.0, "TPC outer r = 78 cm");
  t.SetTextColor(kBlack);
  t.DrawLatex(1.5, -6.5, "vertex (beamline)");
  t.SetTextColor(kBlue + 1);
  t.DrawLatex(74, 30.5, "track, p_{T} = 0.5 GeV");
  t.DrawLatex(1.5, 51.5, "dotted line points at the track-circle center (0, 119)");
  t.SetTextColor(kRed + 1);
  t.SetTextAngle(24);
  t.DrawLatex(30, 12.6, "chord c = 58.8 cm");
  t.SetTextAngle(0);
  t.SetTextColor(kBlack);
  t.DrawLatex(50.5, 8.0, "sagitta s = 3.6 cm");

  // the numeric chain (each number computed from the previous one)
  TLatex n; n.SetTextSize(0.030);
  n.DrawLatex(38, -2.0, "1)  R = p_{T}/(0.3B) = 0.5/(0.3#times1.4) = 1.19 m = 119 cm");
  n.DrawLatex(38, -5.5, "2)  TPC sees the piece with 20 < r < 78 cm");
  n.DrawLatex(38, -9.0, "3)  chord of that piece: c = 58.8 cm (#approx 78#minus20)");
  n.DrawLatex(38, -12.5, "4)  sagitta s = c^{2}/(8R) = 58.8^{2}/(8#times119) = 3.6 cm");

  cv->SaveAs("../sim_validation_plots/circle_geometry_sketch.png");
  printf("wrote ../sim_validation_plots/circle_geometry_sketch.png\n");
}
