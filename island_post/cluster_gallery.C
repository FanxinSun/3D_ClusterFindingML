// cluster_gallery.C — persistent producer for sim_validation_plots/cluster_gallery.png
// (original AuAu-era figure came from an inline throwaway; reconstructed under the
//  provenance rule with a DOCUMENTED selection, exam6/pAu-v2 sim through P0-P3.)
//
// "Funny" cluster gallery: pad x tbin ADC maps of asymmetric clusters,
//   selection: size in [12,80] && |asym| >= 0.25, max 2 per event (variety),
//   top 2 rows REAL (canonical cuts: layer 7-54 && adc>0),
//   bottom 2 rows SIM = composed v3.5revC (B3 readout, 390 kHz).
// Panels show the island bounding box +-1 pad/tbin (neighbourhood display).

#include <TCanvas.h>
#include <TROOT.h>
#include <TFile.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TTree.h>
#include <cmath>
#include <cstdio>
#include <map>
#include <vector>

namespace FS
{
struct Sel
{
  int event, layer, side, plo, phi_, tlo, thi;
  float size, asym;
};

std::vector<Sel> pick(const char *islfile, int nwant)
{
  TFile *f = TFile::Open(islfile);
  TTree *t = (TTree *) f->Get("island");
  float ev, lay, side, size, cpb, ctb, psz, zsz, asym;
  t->SetBranchAddress("event", &ev);
  t->SetBranchAddress("layer", &lay);
  t->SetBranchAddress("side", &side);
  t->SetBranchAddress("size", &size);
  t->SetBranchAddress("cphibin", &cpb);
  t->SetBranchAddress("ctbin", &ctb);
  t->SetBranchAddress("phisize", &psz);
  t->SetBranchAddress("zsize", &zsz);
  t->SetBranchAddress("asym", &asym);
  std::vector<Sel> out;
  std::map<int, int> perev;
  Long64_t N = t->GetEntries();
  for (Long64_t i = 0; i < N && (int) out.size() < nwant; ++i)
  {
    t->GetEntry(i);
    if (size < 12 || size > 80 || std::fabs(asym) < 0.25)
    {
      continue;
    }
    if (perev[(int) ev] >= 2)
    {
      continue;
    }
    perev[(int) ev]++;
    Sel s;
    s.event = (int) ev;
    s.layer = (int) lay;
    s.side = (int) side;
    s.plo = (int) std::lround(cpb - psz / 2) - 1;
    s.phi_ = (int) std::lround(cpb + psz / 2) + 1;
    s.tlo = (int) std::lround(ctb - zsz / 2) - 1;
    s.thi = (int) std::lround(ctb + zsz / 2) + 1;
    s.size = size;
    s.asym = asym;
    out.push_back(s);
  }
  f->Close();
  printf("cluster_gallery: %zu candidates selected from %s\n", out.size(), islfile);
  return out;
}

// one pass over the pixel tree fills all panels for this dataset
void fill(const char *pixfile, bool isSim, std::vector<Sel> &sel, std::vector<TH2D *> &hh)
{
  for (size_t i = 0; i < sel.size(); ++i)
  {
    const Sel &s = sel[i];
    hh.push_back(new TH2D(Form("h%s%zu", isSim ? "s" : "r", i), "",
                          s.phi_ - s.plo + 1, s.plo - 0.5, s.phi_ + 0.5,
                          s.thi - s.tlo + 1, s.tlo - 0.5, s.thi + 0.5));
  }
  TFile *f = TFile::Open(pixfile);
  TTree *h = (TTree *) f->Get("ntp_hit");
  float ev, lay, pb, tb, adc, side;
  h->SetBranchStatus("*", 0);
  for (const char *b : {"event", "layer", "phibin", "adc", "zelem"})
  {
    h->SetBranchStatus(b, 1);
  }
  h->SetBranchStatus(isSim ? "zbin" : "tbin", 1);
  h->SetBranchAddress("event", &ev);
  h->SetBranchAddress("layer", &lay);
  h->SetBranchAddress("phibin", &pb);
  h->SetBranchAddress(isSim ? "zbin" : "tbin", &tb);
  h->SetBranchAddress("adc", &adc);
  h->SetBranchAddress("zelem", &side);
  Long64_t N = h->GetEntries();
  for (Long64_t k = 0; k < N; ++k)
  {
    h->GetEntry(k);
    if (lay < 7 || lay > 54 || adc <= 0)  // canonical selection
    {
      continue;
    }
    for (size_t i = 0; i < sel.size(); ++i)
    {
      const Sel &s = sel[i];
      if ((int) ev != s.event || (int) lay != s.layer || ((int) side == 1 ? 1 : 0) != s.side)
      {
        continue;
      }
      if (pb < s.plo || pb > s.phi_ || tb < s.tlo || tb > s.thi)
      {
        continue;
      }
      hh[i]->Fill(pb, tb, adc);
    }
  }
  f->Close();
}
}  // namespace FS
using namespace FS;

void cluster_gallery(const char *realpix =
                      "/home/rog/sPHENIX/3D_ClusterFindingML/clusters_seeds_island_79507-0.root_ntuplizer.root",
                  const char *realisl = "island_real.root",
                  const char *simpix = "digi_frames_production_v35.root",
                  const char *simisl = "island_frames_v35.root")
{
  gROOT->SetBatch(1);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);
  const int NSHOW = 8;
  auto rsel = pick(realisl, NSHOW);
  auto ssel = pick(simisl, NSHOW);
  std::vector<TH2D *> rh, sh;
  fill(realpix, false, rsel, rh);
  fill(simpix, true, ssel, sh);

  TCanvas c("c", "", 1800, 950);
  c.Divide(4, 4, 0.004, 0.008);
  TLatex tx;
  tx.SetNDC(1);
  tx.SetTextSize(0.09);
  for (int i = 0; i < NSHOW && i < (int) rh.size(); ++i)
  {
    c.cd(i + 1);
    rh[i]->SetTitle(";pad;tbin");
    rh[i]->Draw("COLZ");
    tx.DrawLatex(0.14, 0.86, Form("REAL  L%d  n=%.0f  asym=%.2f", rsel[i].layer, rsel[i].size, rsel[i].asym));
  }
  for (int i = 0; i < NSHOW && i < (int) sh.size(); ++i)
  {
    c.cd(NSHOW + i + 1);
    sh[i]->SetTitle(";pad;tbin");
    sh[i]->Draw("COLZ");
    tx.DrawLatex(0.14, 0.86, Form("SIM pAu v3.5revC  L%d  n=%.0f  asym=%.2f", ssel[i].layer, ssel[i].size, ssel[i].asym));
  }
  c.SaveAs("/home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/cluster_gallery.png");
  printf("cluster_gallery.png regenerated (real vs v3.5revC)\n");
}
