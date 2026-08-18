// residual_fingerprint.C — Q4 verification (2026-07-16): is the ~25% content
// gap a DISTINCT background population or a pure collision-rate scale?
// Method: real(x) - s * sim_v40b(x), s = expected/measured content rate
// (defaults 590/846 = 0.6974), per-frame normalized. If the residual's shape
// differs across observables from collision content -> background (structured
// fingerprint). If residual/real is ~constant everywhere -> scale-like (eps).
#include "canon.h"
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <algorithm>

namespace RFP
{
const char *SIMTAG = "SIM v4.0";
void panel(TVirtualPad *p, TH1D *hr, TH1D *hs, double s, const char *ti, bool logy)
{
  p->cd();
  TH1D *res = (TH1D *) hr->Clone(Form("%s_res", hr->GetName()));
  res->Add(hs, -s);
  hr->SetLineColor(kBlue + 1); hr->SetLineWidth(2);
  hs->SetLineColor(kMagenta + 1); hs->SetLineWidth(2); hs->Scale(s);
  res->SetLineColor(kGreen + 2); res->SetLineWidth(2); res->SetFillColorAlpha(kGreen - 8, 0.4);
  for (TH1D *h : {hr, hs, res}) h->SetStats(0);
  hr->SetTitle(ti);
  double mx = hr->GetMaximum();
  hr->SetMaximum(logy ? mx * 3 : mx * 1.3);
  if (logy) { gPad->SetLogy(); hr->SetMinimum(std::max(1e-1, res->GetMinimum())); }
  else hr->SetMinimum(std::min(0., res->GetMinimum() * 1.2));
  hr->Draw("HIST");
  hs->Draw("HIST SAME");
  res->Draw("HIST SAME");
  TLegend *L = new TLegend(0.45, 0.68, 0.89, 0.89);
  L->SetBorderSize(0); L->SetFillStyle(0);
  L->AddEntry(hr, "REAL (per frame)", "l");
  L->AddEntry(hs, Form("%s #times %.2f (collision part)", SIMTAG, s), "l");
  L->AddEntry(res, "RESIDUAL (unexplained)", "f");
  L->Draw();
}
}  // namespace RFP
using namespace RFP;

// simz: sim apparent-z expression — "z" for digi (stores apparent z directly);
// "z-105.5*(zelem==0)" for hit69 exports (real storage convention, same corr as real).
// realfile/NR: default = as-recorded 100 events; pass real_complete62_hits.root
// with NR=62 for the COMPLETE-window reference (dual-reference directive).
void residual_fingerprint(double s = 0.6974,   // NR default 99: the real as-recorded reference is laser-vetoed (event 44) since 2026-08-17
                          const char *simfile = "digi_frames_production_v40b.root",
                          const char *simz = "z", const char *suffix = "",
                          const char *tag = "SIM v4.0",
                          const char *realfile = "/home/rog/sPHENIX/3D_ClusterFindingML/clusters_seeds_island_79507-0.root_ntuplizer.root",
                          double NR = 99.)
{
  RFP::SIMTAG = tag;
  gROOT->SetBatch(1); gStyle->SetOptStat(0); gStyle->SetTitleFontSize(0.055);
  TFile *fr = TFile::Open(realfile);
  TTree *tr = (TTree *) fr->Get("ntp_hit");
  TFile *fs = TFile::Open(simfile);
  TTree *ts = (TTree *) fs->Get("ntp_hit");
  const double NS = 250.;  // sim frames

  auto mk = [&](TTree *t, const char *v, const char *hn, int nb, double lo, double hi,
                const char *cut, double nfr) {
    TH1D *h = new TH1D(hn, "", nb, lo, hi);
    t->Draw(Form("%s>>%s", v, hn), cut, "goff");
    h->Scale(1. / nfr);
    return h;
  };
  TCanvas c("c", "", 1500, 1250); c.Divide(2, 3);
  // 1. arrivals (envelope shape), flash window excluded via canon cut.
  // Axis ends at 960.5: real recording is live only to tbin ~961 (0.04% of
  // real hits above 960) while the sim window runs to 970 (0.96% of px) —
  // beyond 960 the panel would compare a live sim bin against a dead real
  // bin and report the window-length difference as "unexplained" residual.
  panel(c.cd(1), mk(tr, "tbin", "r1", 96, 0.5, 960.5, CANON::TPC_CUT_NOLASER, NR),
        mk(ts, "tbin", "s1", 96, 0.5, 960.5, CANON::SIM_NOLASER, NS), s,
        "arrival time;tbin;hits/frame", false);
  // 2. layer profile (radial structure)
  panel(c.cd(2), mk(tr, "layer", "r2", 48, 6.5, 54.5, CANON::TPC_CUT, NR),
        mk(ts, "layer", "s2", 48, 6.5, 54.5, "", NS), s,
        "layer profile;layer;hits/frame", false);
  // 3. per-hit adc spectrum
  panel(c.cd(3), mk(tr, "adc", "r3", 100, 0.5, 400.5, CANON::TPC_CUT, NR),
        mk(ts, "adc", "s3", 100, 0.5, 400.5, "", NS), s,
        "per-hit ADC;ADU;hits/frame", true);
  // 4. phi (sector structure)
  panel(c.cd(4), mk(tr, "phi", "r4", 90, -3.15, 3.15, CANON::TPC_CUT, NR),
        mk(ts, "phi", "s4", 90, -3.15, 3.15, "", NS), s,
        "azimuth;#phi;hits/frame", false);
  // 5. apparent z
  panel(c.cd(5), mk(tr, "z-105.5*(zelem==0)", "r5", 88, -110, 110, CANON::TPC_CUT, NR),
        mk(ts, simz, "s5", 88, -110, 110, "", NS), s,
        "apparent z;z [cm];hits/frame", false);
  // 6. residual fraction per region (the scalar fingerprint)
  c.cd(6);
  TH1D *fr6 = new TH1D("fr6", "residual fraction of REAL;;fraction", 5, 0, 5);
  const char *labs[5] = {"total", "R1", "R2", "R3", "adc<30"};
  const char *rc[5] = {CANON::TPC_CUT, "layer>=7&&layer<23&&adc>0", "layer>=23&&layer<39&&adc>0",
                       "layer>=39&&layer<=54&&adc>0", "layer>=7&&layer<=54&&adc>0&&adc<30"};
  const char *sc[5] = {"", "layer<23", "layer>=23&&layer<39", "layer>=39", "adc<30"};
  for (int i = 0; i < 5; ++i)
  {
    double R = tr->GetEntries(rc[i]) / NR;
    double S = ts->GetEntries(sc[i][0] ? sc[i] : "1") / NS * s;
    fr6->SetBinContent(i + 1, (R - S) / R);
    fr6->GetXaxis()->SetBinLabel(i + 1, labs[i]);
    printf("RESID %-6s real %.0f sim*s %.0f residual %.0f (%.1f%% of real)\n",
           labs[i], R, S, R - S, 100. * (R - S) / R);
  }
  fr6->SetStats(0); fr6->SetFillColorAlpha(kGreen - 8, 0.5); fr6->SetLineColor(kGreen + 2);
  fr6->SetMinimum(-0.5); fr6->SetMaximum(0.5);  // negative = sim above real
  fr6->Draw("HIST");
  TLine *z6 = new TLine(0, 0, 5, 0); z6->SetLineColor(kBlack); z6->Draw();
  c.SaveAs(Form("/home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/residual_fingerprint%s.png", suffix));
  printf("saved residual_fingerprint%s.png\n", suffix);
}
