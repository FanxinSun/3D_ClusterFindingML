// pp_accept.C — permanent pixel-level acceptance meters for a pp production
// (2026-08-17). Definitions (all vs the COMPLETE-62 real reference):
//   WINDOWS w1/w2/w3 = fractions of hits in tbin 60-240 / 270-600 / 650-950 of
//     the 60-950 total, the flash window 318-345 EXCLUDED everywhere;
//   BUMP = mean hits/tbin over 60-240 divided by mean over 270-450 with the
//     flash window 318-345 EXCLUDED (flash-blind since v5.6 — the old
//     definition let event 44's flash and the sim's injected flash both sit in
//     the denominator; real target moved 1.1039 -> 1.1376);
//   STEP = hits in 230-246 over hits in 254-270 (drift-edge sharpness);
//   px/event, CV, low-adc/event; islands/event, <size>, <phisize>.
// Real targets: island_post/targets_v6.txt.
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <cmath>
#include <cstdio>
#include <map>
void pp_accept(const char *digi = "digi_frames_production_v6.root",
               const char *island = "island_frames_v6.root", const char *tag = "v6", double nfr = 250., int isreal = 0)
{
  TFile *fs = TFile::Open(digi);
  TTree *ts = (TTree *) fs->Get("ntp_hit");
  // canonical TPC selection (layer 7-54, adc>0): a no-op on sim digi files (which
  // hold only kept TPC pixels) but REQUIRED on real slim files (silicon layers,
  // adc==0 rows). Real-side laser veto (event 44) is applied by construction in
  // the laser-vetoed reference files.
  // isreal=1: also apply the laser veto (event 44) - needed only when reading the
  // raw ntuplizer (the slim reference files are vetoed by construction).
  const char *CUT = isreal ? "layer>=7&&layer<=54&&adc>0&&event!=44" : "layer>=7&&layer<=54&&adc>0";
  TH1D h("h", "", 971, -0.5, 970.5);
  ts->Draw("tbin>>h", CUT, "goff");
  double tot = 0, w1 = 0, w2 = 0, w3 = 0;
  for (int b = 60; b <= 950; ++b)
  {
    if (b >= 318 && b <= 345) continue;
    double v = h.GetBinContent(b + 1);
    tot += v;
    if (b <= 240) w1 += v;
    else if (b >= 270 && b <= 600) w2 += v;
    else if (b >= 650) w3 += v;
  }
  double pre = 0, post = 0, a = 0, b2 = 0;
  for (int x = 230; x <= 246; ++x) pre += h.GetBinContent(x + 1);
  for (int x = 254; x <= 270; ++x) post += h.GetBinContent(x + 1);
  for (int x = 60; x <= 240; ++x) a += h.GetBinContent(x + 1);
  for (int x = 270; x <= 450; ++x) if (!(x >= 318 && x <= 345)) b2 += h.GetBinContent(x + 1);
  float ev, adc, lay;
  ts->SetBranchStatus("*", 0);
  ts->SetBranchStatus("event", 1); ts->SetBranchAddress("event", &ev);
  ts->SetBranchStatus("adc", 1); ts->SetBranchAddress("adc", &adc);
  ts->SetBranchStatus("layer", 1); ts->SetBranchAddress("layer", &lay);
  std::map<int, long> n;
  long la = 0;
  for (Long64_t i = 0; i < ts->GetEntries(); ++i)
  {
    ts->GetEntry(i);
    if (lay < 7 || lay > 54 || adc <= 0) continue;
    if (isreal && (int) ev == 44) continue;
    n[(int) ev]++;
    if (adc < 30) la++;
  }
  double s = 0, s2 = 0;
  for (auto &kv : n) { s += kv.second; s2 += 1.0 * kv.second * kv.second; }
  double m = s / n.size(), sd = std::sqrt(s2 / n.size() - m * m);
  printf("ACC %s: WINDOWS %.4f %.4f %.4f | BUMP(flash-blind) %.4f | STEP %.4f | CV %.3f | px/ev %.0f | lowadc/ev %.0f\n",
         tag, w1 / tot, w2 / tot, w3 / tot, (a / 181.) / (b2 / 153.), pre / post, sd / m, m, (double) la / n.size());
  printf("ACC targets (complete-61, laser-vetoed): 0.2278 0.3431 0.3411 | 1.1396 | 1.1192 | 0.230 | 320586 | 90657\n");
  TFile *fi = TFile::Open(island);
  TTree *t = (TTree *) fi->Get("island");
  // FULL-RANGE island means (2026-08-17 fix): the old real target 10.79 was a
  // full-range mean while sim values were computed over 0-200 -> the "-10%
  // <size>" residual was largely a definition mismatch (real 0-200 = 9.82).
  TH1D h1("h1", "", 2000, 0, 2000), h2("h2", "", 100, 0, 100);
  t->Draw("size>>h1", "", "goff");
  t->Draw("phisize>>h2", "", "goff");
  printf("ACC %s islands: %.0f/ev | <size> %.2f | <phisize> %.2f  (complete-61 targets 29656 | 10.79 | 3.37)\n",
         tag, t->GetEntries() / nfr, h1.GetMean(), h2.GetMean());
}
