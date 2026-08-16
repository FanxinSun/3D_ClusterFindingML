// looper_view.C — SIM twin of the supervisor's real event-74 3D cluster view
// (Visualizing/images/event74_ntp_cluster_py.pdf, made by
// Visualizing/src/plot_event74_hits.py): all ntp_cluster rows with zelem==1
// and 7 < layer < 55, drawn as (z, x, y) points coloured by adc, view
// elev 24 / azim 62. Rendered in ROOT because no matplotlib is available on
// this host; the projection and cuts are the same, the styling is a close
// approximation (marker size, blue->yellow adc colouring).
//
// Two panels: REAL event 74 (left, from the canonical ntuplizer file) and a
// SIM event (right, island91 v5.5). Purpose: show that the sim clusters
// contain the same "small looper" objects — tight low-pT helices spanning
// the drift — that catch the eye in the real view.
//
// usage: root -l -b -q 'looper_view.C+(236)'   -> sim_validation_plots/looper_view_ev236_v55.png
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH3F.h>
#include <TPolyMarker3D.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TView.h>
#include <TLatex.h>
#include <TColor.h>
#include <cstdio>
#include <vector>
#include <cmath>
#include <algorithm>
#include <map>
#include <tuple>

namespace LV {
struct Pt { float x, y, z, adc; };

std::vector<Pt> load(const char* fn, int ev, bool needSimGuard) {
  std::vector<Pt> out;
  TFile* f = TFile::Open(fn);
  TTree* t = (TTree*) f->Get("ntp_cluster");
  float e, x, y, z, adc, layer, zelem;
  t->SetBranchStatus("*", 0);
  for (auto b : {"event", "x", "y", "z", "adc", "layer", "zelem"}) t->SetBranchStatus(b, 1);
  t->SetBranchAddress("event", &e);   t->SetBranchAddress("x", &x);   t->SetBranchAddress("y", &y);
  t->SetBranchAddress("z", &z);       t->SetBranchAddress("adc", &adc);
  t->SetBranchAddress("layer", &layer); t->SetBranchAddress("zelem", &zelem);
  for (Long64_t i = 0; i < t->GetEntries(); ++i) {
    t->GetEntry(i);
    if ((int) e != ev) continue;
    if ((int) zelem != 1 || !(layer > 7 && layer < 55)) continue;
    out.push_back({x, y, z, adc});
  }
  f->Close();
  return out;
}

// draw one panel: (z, x, y) axes like the reference, adc-coloured markers.
// Empty TH3F frame + TPolyMarker3D per adc tertile, pad theta/phi set BEFORE the
// markers are drawn (the verified batch-safe path).
void panel(TVirtualPad* p, std::vector<Pt>& v, const char* title, double zlo, double zhi) {
  p->cd();
  TH3F* fr = new TH3F(Form("fr_%s", title), Form("%s;z (cm);x (cm);y (cm)", title), 2, zlo, zhi, 2, -80, 80, 2, -80, 80);
  fr->SetStats(0);
  fr->GetXaxis()->SetTitleOffset(1.7); fr->GetYaxis()->SetTitleOffset(1.9); fr->GetZaxis()->SetTitleOffset(1.2);
  for (auto* ax : {fr->GetXaxis(), fr->GetYaxis(), fr->GetZaxis()}) ax->SetLabelSize(0.026);
  fr->Draw();
  p->SetTheta(24.); p->SetPhi(-118.);  // reproduces the notebook's view_init(elev=24, azim=62) orientation in ROOT's convention
  std::vector<double> adcs; for (auto& q : v) adcs.push_back(q.adc);
  std::sort(adcs.begin(), adcs.end());
  double a50 = adcs[adcs.size() * 0.50], a90 = adcs[adcs.size() * 0.90];
  int cols[3] = {TColor::GetColor("#3B5BA5"), TColor::GetColor("#4FA3A0"), TColor::GetColor("#E3C13A")};
  TPolyMarker3D* pm[3];
  for (int k = 0; k < 3; k++) { pm[k] = new TPolyMarker3D(); pm[k]->SetMarkerStyle(20); pm[k]->SetMarkerSize(0.30); pm[k]->SetMarkerColor(cols[k]); }
  for (auto& q : v) { int k = q.adc < a50 ? 0 : (q.adc < a90 ? 1 : 2); pm[k]->SetNextPoint(q.z, q.x, q.y); }
  for (int k = 0; k < 3; k++) pm[k]->Draw();
  TLatex L; L.SetNDC(); L.SetTextSize(0.024); L.SetTextColor(kGray + 2);
  L.DrawLatex(0.02, 0.02, Form("%zu clusters, zelem==1, 7<layer<55; colour = adc tertiles (blue < teal < yellow)", v.size()));
}

// ---- finder-free small-looper census (no truth on either side) ----
// Step 1: connect clusters into chains by 3D proximity (link if within LINK cm),
//         the finder-free analogue of "the eye follows a curl".
// Step 2: for each component with >= NMIN clusters, Kasa circle fit in xy;
//         accept as a small looper if R <= RMAX, xy rms <= TOL, and the
//         azimuth about the fitted centre unwraps to >= TURNS turns while z
//         progresses (z-span <= ZMAX). Applied identically to real and sim.
struct Curl { double cx, cy, R, rms, turns, zspan; std::vector<int> mem; };
std::vector<Curl> census(std::vector<Pt>& v, double RMAX = 12., double TOL = 0.6, int NMIN = 30, double TURNS = 2.5,
                         double LINK = 3.0, double ZMAX = 130.) {
  std::vector<Curl> out; const int n = v.size();
  // 3D grid for linking
  const double cell = LINK; std::map<std::tuple<int,int,int>, std::vector<int>> grid;
  auto key = [&](const Pt& p){ return std::make_tuple((int) floor(p.x/cell), (int) floor(p.y/cell), (int) floor(p.z/cell)); };
  for (int i = 0; i < n; i++) grid[key(v[i])].push_back(i);
  std::vector<int> comp(n, -1); int nc = 0; std::vector<int> stack;
  for (int i = 0; i < n; i++) {
    if (comp[i] >= 0) continue;
    comp[i] = nc; stack.assign(1, i);
    while (!stack.empty()) { int a = stack.back(); stack.pop_back(); auto k = key(v[a]);
      for (int dx=-1; dx<=1; dx++) for (int dy=-1; dy<=1; dy++) for (int dz=-1; dz<=1; dz++) {
        auto it = grid.find(std::make_tuple(std::get<0>(k)+dx, std::get<1>(k)+dy, std::get<2>(k)+dz)); if (it == grid.end()) continue;
        for (int b : it->second) { if (comp[b] >= 0) continue;
          double d = sqrt(pow(v[a].x-v[b].x,2)+pow(v[a].y-v[b].y,2)+pow(v[a].z-v[b].z,2));
          if (d <= LINK) { comp[b] = nc; stack.push_back(b); } } } }
    nc++;
  }
  std::vector<std::vector<int>> members(nc); for (int i = 0; i < n; i++) members[comp[i]].push_back(i);
  for (auto& mem : members) {
    if ((int) mem.size() < NMIN) continue;
    // Kasa fit
    double Sx=0,Sy=0,Sxx=0,Syy=0,Sxy=0,Sxz=0,Syz=0,Sz=0; int m = mem.size();
    for (int q : mem) { double X=v[q].x, Y=v[q].y, Z=X*X+Y*Y; Sx+=X; Sy+=Y; Sxx+=X*X; Syy+=Y*Y; Sxy+=X*Y; Sxz+=X*Z; Syz+=Y*Z; Sz+=Z; }
    double M[3][3]={{Sxx,Sxy,Sx},{Sxy,Syy,Sy},{Sx,Sy,(double)m}}, b[3]={-Sxz,-Syz,-Sz}; bool ok=true;
    for (int i=0;i<3&&ok;i++){ int piv=i; for(int k=i+1;k<3;k++) if(fabs(M[k][i])>fabs(M[piv][i])) piv=k;
      std::swap(M[i],M[piv]); std::swap(b[i],b[piv]); if(fabs(M[i][i])<1e-9){ok=false;break;}
      for(int k=i+1;k<3;k++){ double f=M[k][i]/M[i][i]; for(int j=i;j<3;j++) M[k][j]-=f*M[i][j]; b[k]-=f*b[i]; } }
    if (!ok) continue;
    double sol[3]; for(int i=2;i>=0;i--){ sol[i]=b[i]; for(int j=i+1;j<3;j++) sol[i]-=M[i][j]*sol[j]; sol[i]/=M[i][i]; }
    Curl c; c.cx=-sol[0]/2; c.cy=-sol[1]/2; c.R=sqrt(std::max(0.0,c.cx*c.cx+c.cy*c.cy-sol[2]));
    if (c.R > RMAX || c.R < 1.0) continue;
    double rms=0; for (int q : mem){ double d=hypot(v[q].x-c.cx,v[q].y-c.cy)-c.R; rms+=d*d; } c.rms=sqrt(rms/m);
    if (c.rms > TOL) continue;
    std::vector<std::pair<float,double>> zp; for (int q : mem) zp.push_back({v[q].z, atan2(v[q].y-c.cy, v[q].x-c.cx)});
    std::sort(zp.begin(), zp.end());
    double tot=0; for (size_t k=1;k<zp.size();k++){ double dd=zp[k].second-zp[k-1].second; while(dd>M_PI)dd-=2*M_PI; while(dd<-M_PI)dd+=2*M_PI; tot+=fabs(dd); }
    c.turns = tot/(2*M_PI); c.zspan = zp.back().first - zp.front().first;
    if (c.turns < TURNS || c.zspan > ZMAX) continue;
    c.mem = mem; out.push_back(c);
  }
  return out;
}
}  // namespace LV

void looper_view(int simev = 236,
                 const char* simf = "island91_frames_production_v55.root",
                 const char* realf = "../clusters_seeds_island_79507-0.root_ntuplizer.root",
                 int realev = 74, const char* ver = "v55")
{
  gROOT->SetBatch(1); gStyle->SetOptStat(0);
  auto R = LV::load(realf, realev, false);
  auto S = LV::load(simf, simev, true);
  printf("looper_view: real ev%d %zu clusters | sim ev%d %zu clusters (zelem==1, 7<layer<55)\n",
         realev, R.size(), simev, S.size());
  // isolate small loopers with the same finder-free census on both sides
  auto cR = LV::census(R), cS = LV::census(S);
  std::vector<LV::Pt> onlyR, onlyS; int tightR = 0, tightS = 0;
  for (auto& c : cR) { for (int q : c.mem) onlyR.push_back(R[q]); if (c.R <= 8 && c.turns >= 4) tightR++; }
  for (auto& c : cS) { for (int q : c.mem) onlyS.push_back(S[q]); if (c.R <= 8 && c.turns >= 4) tightS++; }
  printf("small loopers (linked chains, R<=12 cm, >=2.5 turns): real %zu (%zu clusters) | sim %zu (%zu clusters); tight (R<=8, >=4 turns): real %d | sim %d\n",
         cR.size(), onlyR.size(), cS.size(), onlyS.size(), tightR, tightS);
  TCanvas c("c", "", 2000, 2000);
  c.Divide(2, 2);
  LV::panel(c.cd(1), R, Form("REAL event %d (run 79507): all clusters", realev), -310, 110);
  LV::panel(c.cd(2), S, Form("SIM %s event %d: all clusters", ver, simev), -310, 110);
  LV::panel(c.cd(3), onlyR, Form("REAL event %d: small loopers only, R#leq12 cm, #geq2.5 turns (%zu found)", realev, cR.size()), -310, 110);
  LV::panel(c.cd(4), onlyS, Form("SIM %s event %d: small loopers only, R#leq12 cm, #geq2.5 turns (%zu found)", ver, simev, cS.size()), -310, 110);
  const char* out = Form("/home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/looper_view_ev%d_%s.png", simev, ver);
  c.SaveAs(out);
  printf("wrote %s\n", out);
}
