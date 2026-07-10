// p1_opener.C — the P1 opener measurements (2026-07-09), time-domain method.
// Decomposed the old "x2.8 occupancy residual". Findings (see PIPELINE.md, P1 section):
//   1. real events = 51.5 us streaming frames, NO trigger bump (flat arrivals);
//   2. rmbd GL1 scaler (mean ~240 kHz) confirmed as collision counter by the
//      frame-fluctuation test (sigma/mu = 0.45 ~ 12-15 collisions/frame);
//   3. real content ~18-30k px per MBD-counted collision; sim 491k px/collision,
//      of which 68% G4 secondaries — sim primaries alone (152k) = naive MIP expectation;
//   4. sim pileup bug: arrivals end at 20.8 us => pileup delivered only to +7.6 us (not +14).
//
// usage: root -l -b -q 'p1_opener.C("<real ntuplizer file>","digi_cal.root","island91_sim.root")'

void p1_opener(const char *realfile =
                   "/home/rog/sPHENIX/3D_ClusterFindingML/clusters_seeds_island_79507-0.root_ntuplizer.root",
               const char *simdigi = "digi_cal.root",
               const char *sim91 = "island91_sim.root")
{
  gROOT->SetBatch(1);
  gStyle->SetOptStat(0);

  // ---- (1) real arrivals: pedestal, drainage excess, readout extent ----
  TFile *fr = TFile::Open(realfile);
  TTree *rh = (TTree *) fr->Get("ntp_hit");
  TH1D *ht = new TH1D("ht", "real arrivals;tbin (53 ns);pixels/event/tbin", 1400, 0, 1400);
  rh->Draw("tbin>>ht", "layer>=7&&layer<=54", "goff");
  ht->Scale(1. / 100.);
  int lastbin = 0;
  for (int i = ht->GetNbinsX(); i >= 1; --i)
  {
    if (ht->GetBinContent(i) > 0.5)
    {
      lastbin = i;
      break;
    }
  }
  auto avg = [&](TH1D *h, int a, int b) {
    double s = 0;
    for (int i = a; i <= b; ++i)
    {
      s += h->GetBinContent(i);
    }
    return s / (b - a + 1);
  };
  const double drift_tb = 13200. / 53.;
  double P = avg(ht, 400, 900);
  double tot = 0;
  for (int i = 1; i <= (int) drift_tb; ++i)
  {
    tot += ht->GetBinContent(i);
  }
  double drain = tot - P * drift_tb;
  printf("REAL: readout to tbin %d (%.1f us) | pedestal %.0f px/us | pre-frame drainage %.0fk px\n",
         lastbin, lastbin * 0.053, P / 0.053, drain / 1e3);

  // ---- (2) GL1 scalers ----
  TTree *rc = (TTree *) fr->Get("ntp_cluster");
  for (const char *v : {"rmbd", "rawmbd", "livembd", "rzdc"})
  {
    TH1D hs("hs", "", 1, -1e18, 1e18);
    rc->Draw(Form("%s>>hs", v), "", "goff", 50000);
    printf("REAL scaler %-8s mean %.6g (min %.6g max %.6g)\n", v, hs.GetMean(), rc->GetMinimum(v), rc->GetMaximum(v));
  }

  // ---- (3) frame-fluctuation discriminator ----
  TH2D h2("h2", "", 102, 0, 102, 1, 0, 1e9);
  rh->Draw("0.5:event>>h2", "layer>=7&&layer<=54", "goff");
  std::vector<double> vv;
  for (int e = 0; e < 102; ++e)
  {
    double x = h2.GetBinContent(e + 1, 1);
    if (x > 0)
    {
      vv.push_back(x);
    }
  }
  double s = 0, s2 = 0;
  for (double x : vv)
  {
    s += x;
    s2 += x * x;
  }
  double m = s / vv.size(), sd = std::sqrt(s2 / vv.size() - m * m);
  printf("REAL frames: mean %.0fk sigma %.0fk  sigma/mu = %.2f  (0.35-0.45 => ~12-15 coll/frame; >1 => single)\n",
         m / 1e3, sd / 1e3, sd / m);
  double R = 240.;  // kHz, rmbd mean
  printf("=> Npix per MBD-counted collision = pedestal/R = %.1fk px\n", (P / 0.053) / R);

  // ---- (4) sim: bump, pedestal, pileup-window check ----
  TFile *fs = TFile::Open(simdigi);
  TTree *sh = (TTree *) fs->Get("ntp_hit");
  TH1D *hs2 = new TH1D("hs2", "sim arrivals;tbin;px/event/tbin", 1000, 0, 1000);
  sh->Draw("zbin>>hs2", "", "goff");
  hs2->Scale(1. / 20.);
  int slast = 0;
  for (int i = hs2->GetNbinsX(); i >= 1; --i)
  {
    if (hs2->GetBinContent(i) > 0.5)
    {
      slast = i;
      break;
    }
  }
  double Ps = avg(hs2, 300, 480), stot = 0;
  for (int i = 1; i <= 249; ++i)
  {
    stot += hs2->GetBinContent(i);
  }
  double B = stot - Ps * drift_tb;
  printf("SIM: arrivals to tbin %d (%.1f us => pileup delivered to +%.1f us; designed +14) | Npix/collision = %.0fk\n",
         slast, slast * 0.053, slast * 0.053 - 13.2, B / 1e3);

  // ---- (5) truth decomposition ----
  double totpx = sh->GetEntries();
  double secpx = sh->GetEntries("gtrackID<0");
  printf("SIM decomposition: G4 secondaries %.1f%% of pixels; primaries %.0fk px/event\n",
         100. * secpx / totpx, (totpx - secpx) / 20 / 1e3);
  TFile *f9 = TFile::Open(sim91);
  if (f9 && f9->Get("ntp_truth"))
  {
    TTree *t = (TTree *) f9->Get("ntp_truth");
    double n = t->GetEntries();
    printf("SIM islands: primaries %.1f%% | secondary loopers (pT<0.164) %.1f%% | secondary tracks %.1f%%\n",
           100. * t->GetEntries("gtrackID>=0") / n,
           100. * t->GetEntries("gtrackID<0&&cls==1") / n,
           100. * t->GetEntries("gtrackID<0&&cls==0") / n);
  }

  // ---- figure ----
  TCanvas c("c", "", 900, 550);
  gPad->SetLogy();
  ht->SetLineColor(kBlue + 1);
  ht->SetLineWidth(2);
  ht->Draw("HIST");
  hs2->SetLineColor(kMagenta + 1);
  hs2->SetLineWidth(2);
  hs2->Draw("HIST SAME");
  TLine l(drift_tb, 0.1, drift_tb, ht->GetMaximum());
  l.SetLineStyle(2);
  l.Draw();
  TLegend L(0.5, 0.75, 0.89, 0.89);
  L.SetBorderSize(0);
  L.AddEntry(ht, "REAL arrivals (flat = streaming frame)", "l");
  L.AddEntry(hs2, "SIM arrivals (bump + truncated pileup)", "l");
  L.Draw();
  c.SaveAs("/home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/p1_tbin.png");
}
