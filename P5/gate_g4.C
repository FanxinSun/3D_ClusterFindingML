// gate_g4.C — P5 detachment GATE 1: standalone app vs container reference,
// same 2 pau200 generator events. Robust observables (vertex smearing differs:
// reference had sPHENIX beam params, standalone uses generator vertices).
void gate_g4(const char *sa = "SA10_g4hit.root", const char *ref = "P5_g4hit_eval.root", int nev = 10)
{
  gROOT->SetBatch(1); gStyle->SetOptStat(0);
  TFile *fa = TFile::Open(sa);  TTree *ta = (TTree *) fa->Get("ntp_g4hit");
  TFile *fr = TFile::Open(ref); TTree *tr = (TTree *) fr->Get("ntp_g4hit");
  double tsa = 0, tref = 0, psa = 0, pref = 0, ssa = 0, sref = 0;
  for (int e = 0; e < nev; ++e)
  {
    auto q = [&](TTree *t, const char *c) { return (double) t->GetEntries(Form("event==%d%s", e, c)); };
    printf("GATE ev%d hits: SA %.0f REF %.0f (ratio %.3f) | prim SA %.0f REF %.0f | sec SA %.0f REF %.0f\n",
           e, q(ta, ""), q(tr, ""), q(ta, "") / q(tr, ""),
           q(ta, "&&gtrackID>0"), q(tr, "&&gtrackID>0"),
           q(ta, "&&gtrackID<0"), q(tr, "&&gtrackID<0"));
    TH1D hea("hea", "", 1, -1e9, 1e9), her("her", "", 1, -1e9, 1e9);
    ta->Draw("gedep>>hea", Form("event==%d", e), "goff");
    tr->Draw("gedep>>her", Form("event==%d", e), "goff");
    double ea = hea.GetBinContent(1) * hea.GetMean() > 0 ? 0 : 0;  // placeholder
    printf("GATE ev%d edep sum: SA %.4g REF %.4g (ratio %.3f)\n", e,
           hea.Integral() * hea.GetMean(), her.Integral() * her.GetMean(),
           (hea.Integral() * hea.GetMean()) / (her.Integral() * her.GetMean()));
    tsa += q(ta, ""); tref += q(tr, "");
    psa += q(ta, "&&gtrackID>0"); pref += q(tr, "&&gtrackID>0");
    ssa += q(ta, "&&gtrackID<0"); sref += q(tr, "&&gtrackID<0");
  }
  printf("GATE AGGREGATE (%d ev): hits SA/REF %.3f | prim %.3f | sec %.3f\n",
         nev, tsa / tref, psa / pref, ssa / sref);
  auto cmp1 = [&](const char *v, const char *hn, int nb, double lo, double hi, const char *ti, bool logy, TVirtualPad *p) {
    p->cd();
    TH1D *ha = new TH1D(Form("%sa", hn), "", nb, lo, hi);
    TH1D *hr2 = new TH1D(Form("%sr", hn), "", nb, lo, hi);
    ta->Draw(Form("%s>>%sa", v, hn), "", "goff");
    tr->Draw(Form("%s>>%sr", v, hn), "", "goff");
    if (ha->Integral() > 0) ha->Scale(1. / ha->Integral());
    if (hr2->Integral() > 0) hr2->Scale(1. / hr2->Integral());
    hr2->SetLineColor(kBlue + 1); hr2->SetLineWidth(2);
    ha->SetLineColor(kMagenta + 1); ha->SetLineWidth(2);
    hr2->SetTitle(ti); hr2->SetStats(0); ha->SetStats(0);
    double mx = std::max(ha->GetMaximum(), hr2->GetMaximum());
    hr2->SetMaximum(logy ? mx * 2.5 : mx * 1.3);
    if (logy) { gPad->SetLogy(); hr2->SetMinimum(1e-6); } else hr2->SetMinimum(0);
    hr2->Draw("HIST"); ha->Draw("HIST SAME");
    TLegend *L = new TLegend(0.55, 0.72, 0.89, 0.89);
    L->SetBorderSize(0); L->SetFillStyle(0);
    L->AddEntry(hr2, "REF (ana.558 ctr)", "l");
    L->AddEntry(ha, "STANDALONE (G4 11.4)", "l");
    L->Draw();
  };
  TCanvas c("c", "", 1500, 950); c.Divide(3, 2);
  cmp1("sqrt(gx*gx+gy*gy)", "g1", 60, 20, 80, "hit radius;r [cm];norm", false, c.cd(1));
  cmp1("gz", "g2", 60, -120, 120, "hit z;z [cm];norm", false, c.cd(2));
  cmp1("gpl", "g3", 60, 0, 1.5, "step length;gpl [cm];norm", true, c.cd(3));
  cmp1("gedep*1e6", "g4", 60, 0, 30, "step edep;keV;norm", true, c.cd(4));
  cmp1("gt", "g5", 60, 0, 120, "hit time;ns;norm", true, c.cd(5));
  cmp1("sqrt(gpx*gpx+gpy*gpy)", "g6", 60, 0, 2, "pre-step pT;GeV;norm", true, c.cd(6));
  c.SaveAs("gate_g4.png");
  printf("saved gate_g4.png\n");
}
