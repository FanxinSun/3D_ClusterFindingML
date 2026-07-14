// tank_specimens.C — uniformity-calibration sheet (decision artifact, 2026-07-14).
// Every shape+isolation-passing tank in the windows of interest, drawn with its
// ADC numbers + max/median ratio and difference, PASS/FAIL under both candidate
// uniformity forms: ratio (max/med <= 3) and diff (max-med <= 200).
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TROOT.h>
#include <vector>
#include <array>
#include <map>
#include <algorithm>
#include <cmath>

namespace TSPEC
{
struct Cand { int win; std::vector<std::array<int,3>> px; double med, mx; };

void collect(const char *fn, bool isReal, int ev, int lay, const char *tag,
             std::vector<Cand> &out, std::vector<std::string> &labels)
{
  TFile *fd = TFile::Open(fn);
  TTree *t = (TTree *) fd->Get("ntp_hit");
  float e, l, pb, tb, adc, ph;
  t->SetBranchStatus("*", 0);
  for (auto bn : {"event", "layer", "phibin", "adc", "phi"}) t->SetBranchStatus(bn, 1);
  t->SetBranchStatus(isReal ? "tbin" : "zbin", 1);
  t->SetBranchAddress("event", &e);
  t->SetBranchAddress("layer", &l);
  t->SetBranchAddress("phibin", &pb);
  t->SetBranchAddress(isReal ? "tbin" : "zbin", &tb);
  t->SetBranchAddress("adc", &adc);
  t->SetBranchAddress("phi", &ph);
  std::vector<std::array<int,3>> px;
  for (Long64_t i = 0; i < t->GetEntries(); ++i)
  {
    t->GetEntry(i);
    if ((int) e != ev || (int) l != lay) continue;
    if (adc <= 0 || ph < 0 || ph >= 1 || tb < 600 || tb > 800) continue;
    px.push_back({(int) pb, (int) tb, (int) adc});
  }
  fd->Close();
  int n = (int) px.size();
  std::vector<int> comp(n, -1);
  int nc = 0;
  for (int i = 0; i < n; ++i)
  {
    if (comp[i] >= 0) continue;
    std::vector<int> st = {i};
    comp[i] = nc;
    while (!st.empty())
    {
      int j = st.back(); st.pop_back();
      for (int k = 0; k < n; ++k)
        if (comp[k] < 0 && std::abs(px[j][0]-px[k][0]) + std::abs(px[j][1]-px[k][1]) == 1)
        { comp[k] = nc; st.push_back(k); }
    }
    nc++;
  }
  auto passes = [](const std::vector<std::array<int,2>> &q) {
    int tmax = -1 << 30;
    for (auto &p : q) tmax = std::max(tmax, p[1]);
    int ntop = 0, xtop = -1 << 30, plo = 1 << 30, phi2 = -1 << 30, nbelow = 0;
    bool onbase = false;
    for (auto &p : q)
    {
      plo = std::min(plo, p[0]); phi2 = std::max(phi2, p[0]);
      if (p[1] == tmax) { ntop++; xtop = p[0]; }
      if (p[1] == tmax - 1) nbelow++;
    }
    for (auto &p : q)
      if (p[1] == tmax - 1 && p[0] == xtop) onbase = true;
    return ntop == 1 && nbelow >= 3 && onbase && xtop != plo && xtop != phi2;
  };
  for (int c = 0; c < nc; ++c)
  {
    std::vector<std::array<int,3>> cp;
    for (int i = 0; i < n; ++i) if (comp[i] == c) cp.push_back(px[i]);
    if (cp.size() < 6 || cp.size() > 8) continue;
    bool istu = false;
    for (int o = 0; o < 4 && !istu; ++o)
    {
      std::vector<std::array<int,2>> q;
      for (auto &p : cp)
      {
        if (o == 0) q.push_back({p[0], p[1]});
        else if (o == 1) q.push_back({p[0], -p[1]});
        else if (o == 2) q.push_back({p[1], p[0]});
        else q.push_back({p[1], -p[0]});
      }
      istu = passes(q);
    }
    if (!istu) continue;
    bool iso = true;
    for (int i = 0; i < n && iso; ++i)
    {
      if (comp[i] == c) continue;
      for (auto &p : cp)
        if (std::abs(px[i][0]-p[0]) <= 1 && std::abs(px[i][1]-p[1]) <= 1) { iso = false; break; }
    }
    if (!iso) continue;
    std::vector<int> adcs;
    for (auto &p : cp) adcs.push_back(p[2]);
    std::sort(adcs.begin(), adcs.end());
    Cand cd; cd.px = cp; cd.med = adcs[adcs.size()/2]; cd.mx = adcs.back();
    out.push_back(cd);
    labels.push_back(tag);
  }
}
}  // namespace TSPEC
using namespace TSPEC;

void tank_specimens()
{
  gROOT->SetBatch(1); gStyle->SetOptStat(0);
  const char *REAL = "/home/rog/sPHENIX/3D_ClusterFindingML/clusters_seeds_island_79507-0.root_ntuplizer.root";
  const char *SIM = "digi_frames_production_v36.root";
  std::vector<Cand> cands; std::vector<std::string> labs;
  collect(REAL, true, 74, 15, "REAL e74 L15", cands, labs);
  collect(SIM, false, 12, 15, "SIM f12 L15", cands, labs);
  collect(SIM, false, 82, 7, "SIM f82 L7", cands, labs);
  collect(SIM, false, 209, 12, "SIM f209 L12", cands, labs);
  collect(SIM, false, 245, 15, "SIM f245 L15", cands, labs);
  int nsp = (int) cands.size();
  printf("tank_specimens: %d candidates\n", nsp);
  int cols = 4, rows = (nsp + cols - 1) / cols;
  TCanvas c("c", "", 380 * cols, 340 * std::max(rows, 1));
  c.Divide(cols, std::max(rows, 1), 0.004, 0.004);
  for (int i = 0; i < nsp; ++i)
  {
    c.cd(i + 1);
    auto &cd = cands[i];
    int plo = 1 << 30, phi2 = -1 << 30, tlo = 1 << 30, thi = -1 << 30;
    for (auto &p : cd.px)
    { plo = std::min(plo, p[0]); phi2 = std::max(phi2, p[0]);
      tlo = std::min(tlo, p[1]); thi = std::max(thi, p[1]); }
    TH2D *h = new TH2D(Form("sp%d", i), "", phi2-plo+5, plo-2.5, phi2+2.5, thi-tlo+5, tlo-2.5, thi+2.5);
    for (auto &p : cd.px) h->Fill(p[0], p[1], p[2]);
    h->SetMaximum(850); h->SetMinimum(0);
    h->Draw("COL");
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.065);
    for (auto &p : cd.px)
    {
      TLatex t2; t2.SetTextSize(0.05); t2.SetTextAlign(22);
      t2.DrawLatex(p[0], p[1], Form("%d", p[2]));
    }
    double ratio = cd.mx / std::max(cd.med, 1.0), diff = cd.mx - cd.med;
    tx.SetTextColor(kBlack);
    tx.DrawLatex(0.08, 0.93, Form("%s", labs[i].c_str()));
    tx.DrawLatex(0.08, 0.02, Form("med %.0f max %.0f | r%.1f %s | d%.0f %s",
      cd.med, cd.mx, ratio, ratio <= 3.0 ? "PASS" : "FAIL",
      diff, diff <= 200 ? "PASS" : "FAIL"));
  }
  c.SaveAs("/home/rog/sPHENIX/3D_ClusterFindingML/island_post/tank_specimens.png");
  printf("saved tank_specimens.png\n");
}
