// asym_showcase.C — 凸-shaped PURE-SINGLE cluster gallery (user spec, 2026-07-11):
//   * ntrks == 1 only (row-aligned truth), size 6-8 (typically 7)
//   * pixel-level 凸 test v4: 4 ORIENTATIONS + UNIFORM brightness (max/med<=3) —
//     one extreme-row pixel, non-corner, ON a >= 3-wide base (in rotated frame)
// Outputs: asym_showcase.png (gallery, <=16) + per-specimen FULL-LAYER context views
// with the specimen circled in red: showcase_ctx_<i>_f<ev>_L<lay>.png
#include <TCanvas.h>
#include <TEllipse.h>
#include <TFile.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <cmath>
#include <cstdio>
#include <map>
#include <vector>

namespace AS2
{
struct Cand
{
  int event, layer, side, plo, phi_, tlo, thi;
  float size, mx, cpb, ctb, cphi;
  TH2D *box = nullptr;
  bool keep = false;
};
}  // namespace AS2
using namespace AS2;

void asym_showcase(const char *isl91 = "island91_frames_production_v40b.root",
                   const char *pix = "digi_frames_production_v40b.root")
{
  gROOT->SetBatch(1);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kViridis);

  // ---- candidates: ntrks==1, size 6-8, wide base geometry ----
  TFile *f = TFile::Open(isl91);
  TTree *c = (TTree *) f->Get("ntp_cluster");
  TTree *t = (TTree *) f->Get("ntp_truth");
  float ev, lay, zel, size, mx, pb, tb, psz, zsz, ntk, cph;
  c->SetBranchStatus("*", 0);
  for (const char *b : {"event", "layer", "zelem", "size", "maxadc", "phibin", "tbin", "phisize", "zsize", "phi"})
  {
    c->SetBranchStatus(b, 1);
  }
  c->SetBranchAddress("event", &ev);
  c->SetBranchAddress("layer", &lay);
  c->SetBranchAddress("zelem", &zel);
  c->SetBranchAddress("size", &size);
  c->SetBranchAddress("maxadc", &mx);
  c->SetBranchAddress("phibin", &pb);
  c->SetBranchAddress("tbin", &tb);
  c->SetBranchAddress("phisize", &psz);
  c->SetBranchAddress("zsize", &zsz);
  c->SetBranchAddress("phi", &cph);
  t->SetBranchStatus("*", 0);
  t->SetBranchStatus("ntrks", 1);
  t->SetBranchAddress("ntrks", &ntk);
  std::vector<Cand> cand;
  std::map<int, int> perframe;
  Long64_t N = c->GetEntries();
  for (Long64_t i = 0; i < N && cand.size() < 60; ++i)
  {
    c->GetEntry(i);
    if (size < 6 || size > 8 || mx > 200 || psz < 3 || zsz < 2 || size > 0.85 * psz * zsz)
    {
      continue;
    }
    t->GetEntry(i);
    if ((int) ntk != 1)
    {
      continue;
    }
    if (perframe[(int) ev] >= 3)
    {
      continue;
    }
    perframe[(int) ev]++;
    Cand s;
    s.event = (int) ev;
    s.layer = (int) lay;
    s.side = ((int) zel == 1) ? 1 : 0;
    s.plo = (int) std::lround(pb - psz / 2) - 1;
    s.phi_ = (int) std::lround(pb + psz / 2) + 1;
    s.tlo = (int) std::lround(tb - zsz / 2) - 1;
    s.thi = (int) std::lround(tb + zsz / 2) + 1;
    s.size = size;
    s.mx = mx;
    s.cpb = pb;
    s.ctb = tb;
    s.cphi = cph;
    cand.push_back(s);
  }
  printf("asym_showcase: %zu ntrks=1 candidates (size 6-8, wide-base)\n", cand.size());

  // ---- pass 1: fill candidate bbox maps ----
  for (auto &s : cand)
  {
    s.box = new TH2D(Form("b%d_%d_%d", s.event, s.layer, (int) s.cpb), "",
                     s.phi_ - s.plo + 1, s.plo - 0.5, s.phi_ + 0.5,
                     s.thi - s.tlo + 1, s.tlo - 0.5, s.thi + 0.5);
    s.box->SetDirectory(nullptr);
  }
  TFile *fp = TFile::Open(pix);
  TTree *h = (TTree *) fp->Get("ntp_hit");
  float hev, hlay, hpb, htb, hadc, hsd;
  h->SetBranchStatus("*", 0);
  for (const char *b : {"event", "layer", "phibin", "zbin", "adc", "zelem", "phi"})
  {
    h->SetBranchStatus(b, 1);
  }
  h->SetBranchAddress("event", &hev);
  h->SetBranchAddress("layer", &hlay);
  h->SetBranchAddress("phibin", &hpb);
  h->SetBranchAddress("zbin", &htb);
  h->SetBranchAddress("adc", &hadc);
  h->SetBranchAddress("zelem", &hsd);
  float hphi;
  h->SetBranchAddress("phi", &hphi);
  Long64_t M = h->GetEntries();
  for (Long64_t k = 0; k < M; ++k)
  {
    h->GetEntry(k);
    for (auto &s : cand)
    {
      if ((int) hev != s.event || (int) hlay != s.layer || (((int) hsd == 1) ? 1 : 0) != s.side)
      {
        continue;
      }
      if (hpb < s.plo || hpb > s.phi_ || htb < s.tlo || htb > s.thi)
      {
        continue;
      }
      s.box->Fill(hpb, htb, hadc);
    }
  }

  // ---- 凸 test on the maps ----
  int kept = 0;
  for (auto &s : cand)
  {
    if (kept >= 16)
    {
      break;
    }
    TH2D *b = s.box;
    int nx = b->GetNbinsX(), ny = b->GetNbinsY();
    // 4-orientation tu test (user 2026-07-14): exact, 90-left/right, 180
    // rotations all in scope. Extract occupied bins, apply the upright test
    // in 4 coordinate transforms.
    std::vector<std::array<int, 2>> cp;
    std::vector<double> adcs;
    for (int ix = 1; ix <= nx; ++ix)
      for (int iy = 1; iy <= ny; ++iy)
        if (b->GetBinContent(ix, iy) > 0)
        {
          cp.push_back({ix, iy});
          adcs.push_back(b->GetBinContent(ix, iy));
        }
    // definition v4: brightness uniformity (max/median <= 3), hot cores read
    // as overlaps to the eye
    std::sort(adcs.begin(), adcs.end());
    if (adcs.empty() || adcs[adcs.size() / 2] <= 0 ||
        adcs.back() / adcs[adcs.size() / 2] > 3.0)
    {
      continue;
    }
    auto passes = [](const std::vector<std::array<int, 2>> &q) {
      int tmax = -1 << 30;
      for (auto &p : q) tmax = std::max(tmax, p[1]);
      int ntop = 0, xtop = -1 << 30, plo = 1 << 30, phi2 = -1 << 30, nbelow = 0;
      bool onbase = false;
      for (auto &p : q)
      {
        plo = std::min(plo, p[0]);
        phi2 = std::max(phi2, p[0]);
        if (p[1] == tmax) { ntop++; xtop = p[0]; }
        if (p[1] == tmax - 1) nbelow++;
      }
      for (auto &p : q)
        if (p[1] == tmax - 1 && p[0] == xtop) onbase = true;
      return ntop == 1 && nbelow >= 3 && onbase && xtop != plo && xtop != phi2;
    };
    bool istu = false;
    for (int o = 0; o < 4 && !istu; ++o)
    {
      std::vector<std::array<int, 2>> q;
      q.reserve(cp.size());
      for (auto &p : cp)
      {
        if (o == 0) q.push_back({p[0], p[1]});
        else if (o == 1) q.push_back({p[0], -p[1]});
        else if (o == 2) q.push_back({p[1], p[0]});
        else q.push_back({p[1], -p[0]});
      }
      istu = passes(q);
    }
    if (!istu)
    {
      continue;
    }
    s.keep = true;
    kept++;
  }
  printf("asym_showcase: %d 凸 specimens kept\n", kept);

  // ---- gallery ----
  TCanvas cg("cg", "", 1800, 950);
  cg.Divide(4, 4, 0.004, 0.008);
  TLatex tx;
  tx.SetNDC(1);
  tx.SetTextSize(0.085);
  int ip = 0;
  for (auto &s : cand)
  {
    if (!s.keep || ip >= 16)
    {
      continue;
    }
    cg.cd(++ip);
    s.box->SetTitle(";pad;tbin");
    s.box->SetMinimum(0);
    s.box->SetMaximum(200);
    s.box->Draw("COL");
    tx.DrawLatex(0.13, 0.86, Form("f%d L%d n=%.0f mx=%.0f ntrks=1", s.event, s.layer, s.size, s.mx));
  }
  cg.SaveAs("/home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/asym_showcase.png");

  // ---- pass 2: full-layer context views with red circles ----
  std::vector<Cand *> keeps;
  for (auto &s : cand)
  {
    if (s.keep)
    {
      keeps.push_back(&s);
    }
  }
  std::vector<TH2D *> ctx;
  for (size_t i = 0; i < keeps.size(); ++i)
  {
    int npads = (keeps[i]->layer <= 22) ? 1128 : (keeps[i]->layer <= 38 ? 1536 : 2304);
    ctx.push_back(new TH2D(Form("ctx%zu", i),
                           Form("SIM frame %d layer %d side %d: #phi-tbin view;#phi [rad];tbin",
                                keeps[i]->event, keeps[i]->layer, keeps[i]->side),
                           npads / 2, 0., 2. * M_PI, 486, -0.5, 971.5));
    ctx.back()->SetDirectory(nullptr);
  }
  for (Long64_t k = 0; k < M; ++k)
  {
    h->GetEntry(k);
    for (size_t i = 0; i < keeps.size(); ++i)
    {
      Cand *s = keeps[i];
      if ((int) hev != s->event || (int) hlay != s->layer || (((int) hsd == 1) ? 1 : 0) != s->side)
      {
        continue;
      }
      double fph = hphi < 0 ? hphi + 2. * M_PI : hphi;
      ctx[i]->Fill(fph, htb, hadc);
    }
  }
  for (size_t i = 0; i < keeps.size(); ++i)
  {
    Cand *s = keeps[i];
    TCanvas cc("cc", "", 1400, 700);
    cc.SetRightMargin(0.12);
    ctx[i]->SetMinimum(0);
    ctx[i]->SetMaximum(300);
    ctx[i]->Draw("COL");
    double fcp = s->cphi < 0 ? s->cphi + 2. * M_PI : s->cphi;
    TEllipse el(fcp, s->ctb, 0.14, 28);
    el.SetLineColor(kRed);
    el.SetLineWidth(3);
    el.SetFillStyle(0);
    el.Draw();
    cc.SaveAs(Form("/home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/showcase_ctx_%02zu_f%d_L%d.png",
                   i, s->event, s->layer));
  }
  printf("asym_showcase: gallery + %zu context views saved\n", keeps.size());
}
