// cluster_shapes_cmp.C — supervisor request (2026-07-15): "compare cluster shapes
// and cluster property distributions between data and simulation."
// Two canvases at the SEALED v3.6revF state:
//   [1] cluster_props_cmp.png  — PRODUCTION SEMANTICS (ported TpcClusterizer,
//       fixed_window=0): real ntp_cluster vs port x real pixels [gate] vs
//       port x SIM v4.0 — the apples-to-apples the collaboration reads.
//   [2] cluster_shapes_cmp.png — ISLAND SEMANTICS (8-connected, true pixel
//       lists): scalar props + SHAPE MOMENTS (asym, rho) + phisize x zsize
//       joint maps, real vs sim. Conventions from canon.h.
#include "canon.h"
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLatex.h>
#include <algorithm>

namespace CSC
{
const char *SIMTAG = "SIM v4.0";
TH1D *mk(TTree *t, const char *v, const char *hn, int nb, double lo, double hi, const char *cut)
{
  TH1D *h = new TH1D(hn, "", nb, lo, hi);
  t->Draw(Form("%s>>%s", v, hn), cut, "goff");
  if (h->Integral() > 0) h->Scale(1. / h->Integral());
  return h;
}
void d2(TVirtualPad *p, TH1D *a, TH1D *b, const char *ti, bool logy,
        TH1D *g = nullptr)
{
  p->cd();
  a->SetLineColor(kBlue + 1); a->SetLineWidth(2);
  b->SetLineColor(kMagenta + 1); b->SetLineWidth(2);
  if (g) { g->SetLineColor(kGray + 2); g->SetLineStyle(2); g->SetLineWidth(2); }
  a->SetStats(0); b->SetStats(0);
  a->SetTitle(ti);
  double mx = std::max(a->GetMaximum(), b->GetMaximum());
  if (g) mx = std::max(mx, g->GetMaximum());
  a->SetMaximum(logy ? mx * 2.5 : mx * 1.3);
  if (logy) { gPad->SetLogy(); a->SetMinimum(1e-7); }
  else a->SetMinimum(0);
  a->Draw("HIST");
  if (g) g->Draw("HIST SAME");
  b->Draw("HIST SAME");
  TLegend *L = new TLegend(0.42, 0.70, 0.89, 0.89);
  L->SetBorderSize(0); L->SetFillStyle(0);
  L->AddEntry(a, "REAL", "l");
  if (g) L->AddEntry(g, "port #times real pixels [gate]", "l");
  L->AddEntry(b, SIMTAG, "l");
  L->Draw();
}
}  // namespace CSC
using namespace CSC;

void cluster_shapes_cmp(const char *suffix = "",
                        const char *simprod = "prodclus_v40b.root",
                        const char *simisl = "island_frames_v40b.root",
                        const char *tag = "SIM v4.0")
{
  CSC::SIMTAG = tag;
  gROOT->SetBatch(1); gStyle->SetOptStat(0); gStyle->SetTitleFontSize(0.055);
  const char *REALNT = "/home/rog/sPHENIX/3D_ClusterFindingML/clusters_seeds_island_79507-0.root_ntuplizer.root";

  // ---------- [1] production semantics ----------
  TFile *fr = TFile::Open(REALNT);
  TTree *rc = (TTree *) fr->Get("ntp_cluster");
  TFile *fg = TFile::Open("prodclus_real_fw0.root"); TTree *gc = (TTree *) fg->Get("ntp_clus");
  TFile *fs = TFile::Open(simprod);                   TTree *sc = (TTree *) fs->Get("ntp_clus");
  const char *RC = "layer>=7&&layer<=54";
  TCanvas c1("c1", "", 1500, 950); c1.Divide(3, 2);
  d2(c1.cd(1), mk(rc, "zsize", "p1", 30, 0.5, 30.5, RC), mk(sc, "zsize", "p1s", 30, 0.5, 30.5, ""),
     "z-size (production semantics);tbins;norm", true, mk(gc, "zsize", "p1g", 30, 0.5, 30.5, ""));
  d2(c1.cd(2), mk(rc, "phisize", "p2", 16, 0.5, 16.5, RC), mk(sc, "phisize", "p2s", 16, 0.5, 16.5, ""),
     "#phi-size (production semantics);pads;norm", true, mk(gc, "phisize", "p2g", 16, 0.5, 16.5, ""));
  d2(c1.cd(3), mk(rc, "adc", "p3", 120, 0, 3600, RC), mk(sc, "adc", "p3s", 120, 0, 3600, ""),
     "cluster adc;ADU;norm", true, mk(gc, "adc", "p3g", 120, 0, 3600, ""));
  d2(c1.cd(4), mk(rc, "maxadc", "p4", 105, 0, 1050, RC), mk(sc, "maxadc", "p4s", 105, 0, 1050, ""),
     "cluster maxadc;ADU;norm", true, mk(gc, "maxadc", "p4g", 105, 0, 1050, ""));
  d2(c1.cd(5), mk(rc, "phisize*zsize", "p5", 60, 0.5, 120.5, RC), mk(sc, "phisize*zsize", "p5s", 60, 0.5, 120.5, ""),
     "bbox area #phi#times z (real `size` semantics);pads#times tbins;norm", true,
     mk(gc, "phisize*zsize", "p5g", 60, 0.5, 120.5, ""));
  d2(c1.cd(6), mk(rc, "adc/(phisize*zsize)", "p6", 80, 0, 400, RC), mk(sc, "adc/(phisize*zsize)", "p6s", 80, 0, 400, ""),
     "charge density adc/area;ADU per cell;norm", true, mk(gc, "adc/(phisize*zsize)", "p6g", 80, 0, 400, ""));
  c1.SaveAs(Form("/home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/cluster_props_cmp%s.png", suffix));
  printf("saved cluster_props_cmp%s.png\n", suffix);

  // ---------- [2] island semantics: shapes ----------
  TFile *fi = TFile::Open("island_real.root");      TTree *ir = (TTree *) fi->Get("island");
  TFile *fj = TFile::Open(simisl);                    TTree *is = (TTree *) fj->Get("island");
  TCanvas c2("c2", "", 1500, 1250); c2.Divide(3, 3);
  d2(c2.cd(1), mk(ir, "size", "s1", 60, 0.5, 60.5, ""), mk(is, "size", "s1s", 60, 0.5, 60.5, ""),
     "island size (TRUE pixel count);pixels;norm", true);
  d2(c2.cd(2), mk(ir, "phisize", "s2", 14, 0.5, 14.5, ""), mk(is, "phisize", "s2s", 14, 0.5, 14.5, ""),
     "island #phi-size;pads;norm", true);
  d2(c2.cd(3), mk(ir, "zsize", "s3", 30, 0.5, 30.5, ""), mk(is, "zsize", "s3s", 30, 0.5, 30.5, ""),
     "island z-size;tbins;norm", true);
  d2(c2.cd(4), mk(ir, "asym", "s4", 50, -1, 1, ""), mk(is, "asym", "s4s", 50, -1, 1, ""),
     "SHAPE: charge asymmetry (signed, along t);asym;norm", false);
  // NOTE: islandize `rho` IS the bbox fill (duplicate of pad 6) - replaced
  // with charge concentration, the observable behind the nbright/hot-core family.
  d2(c2.cd(5), mk(ir, "maxadc/adc", "s5", 50, 0, 1, ""), mk(is, "maxadc/adc", "s5s", 50, 0, 1, ""),
     "SHAPE: charge concentration maxadc/adc;fraction;norm", false);
  d2(c2.cd(6), mk(ir, "size/(phisize*zsize)", "s6", 40, 0, 1.01, ""), mk(is, "size/(phisize*zsize)", "s6s", 40, 0, 1.01, ""),
     "SHAPE: fill factor size/bbox;fill;norm", false);
  // joint phisize x zsize maps, column-normalized visual
  c2.cd(7);
  TH2D *j1 = new TH2D("j1", "REAL: #phi-size #times z-size;pads;tbins", 10, 0.5, 10.5, 20, 0.5, 20.5);
  ir->Draw("zsize:phisize>>j1", "", "goff");
  j1->Scale(1. / j1->Integral()); j1->SetStats(0); gPad->SetLogz(); j1->SetMinimum(1e-6); j1->Draw("COLZ");
  c2.cd(8);
  TH2D *j2 = new TH2D("j2", Form("%s: #phi-size #times z-size;pads;tbins", SIMTAG), 10, 0.5, 10.5, 20, 0.5, 20.5);
  is->Draw("zsize:phisize>>j2", "", "goff");
  j2->Scale(1. / j2->Integral()); j2->SetStats(0); gPad->SetLogz(); j2->SetMinimum(1e-6); j2->Draw("COLZ");
  c2.cd(9);
  TH2D *j3 = (TH2D *) j2->Clone("j3"); j3->Divide(j1);
  j3->SetTitle("ratio SIM/REAL;pads;tbins");
  j3->SetStats(0); j3->SetMinimum(0.2); j3->SetMaximum(5); gPad->SetLogz(); j3->Draw("COLZ");
  c2.SaveAs(Form("/home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/cluster_shapes_cmp%s.png", suffix));
  printf("saved cluster_shapes_cmp%s.png\n", suffix);

  // headline numbers
  for (auto pr : {std::make_pair(ir, "REAL"), std::make_pair(is, "SIM ")})
  {
    TTree *t = pr.first;
    TH1D ha("ha", "", 50, -1, 1), hb("hb", "", 50, 0, 1), hc("hc", "", 50, 0, 1);
    t->Draw("asym>>ha", "", "goff"); t->Draw("abs(asym)>>hb", "", "goff");
    t->Draw("maxadc/adc>>hc", "", "goff");
    printf("SHAPES %s: <asym signed> %+.3f | <|asym|> %.3f | <maxadc/adc> %.3f\n",
           pr.second, ha.GetMean(), hb.GetMean(), hc.GetMean());
  }
}
