// z_bridge.C — replot of the earliest-era "are the z conventions bridgable?" demo,
// now with exam6/pAu v3.5revC frames: both datasets carry APPARENT z with the same mapping
// z = +-(105.5 - v_d * t). Panels: real z:tbin map | sim z:tbin map | signed-z overlay.
void z_bridge(){
  gROOT->SetBatch(1); gStyle->SetOptStat(0); gStyle->SetPalette(kViridis);
  TFile*fr=TFile::Open("/home/rog/sPHENIX/3D_ClusterFindingML/clusters_seeds_island_79507-0.root_ntuplizer.root");
  TTree*rh=(TTree*)fr->Get("ntp_hit");
  TFile*fs=TFile::Open("digi_frames_production_v35.root"); TTree*sh=(TTree*)fs->Get("ntp_hit");
  TCanvas c("c","",1800,600); c.Divide(3,1);
  c.cd(1); gPad->SetRightMargin(0.13); gPad->SetLeftMargin(0.12); gPad->SetLogz();
  TH2D*h1=new TH2D("h1","REAL: apparent z vs tbin;tbin;z [cm]",194,0,970,140,-320,320);
  rh->Draw("(z-105.5*(zelem==0)):tbin>>h1","layer>=7&&layer<=54&&adc>0","goff");  // canon: real side0 z quirk
  h1->GetYaxis()->SetTitleOffset(1.15); h1->Draw("COLZ");
  c.cd(2); gPad->SetRightMargin(0.13); gPad->SetLeftMargin(0.12); gPad->SetLogz();
  TH2D*h2=new TH2D("h2","SIM v3.5revC: apparent z vs tbin;tbin;z [cm]",194,0,970,140,-320,320);
  sh->Draw("z:zbin>>h2","","goff");
  h2->GetYaxis()->SetTitleOffset(1.15); h2->Draw("COLZ");
  c.cd(3); gPad->SetLogy();
  TH1D*z1=new TH1D("z1","signed apparent z;z [cm];norm",160,-320,320);
  TH1D*z2=new TH1D("z2","",160,-320,320);
  rh->Draw("(z-105.5*(zelem==0))>>z1","layer>=7&&layer<=54&&adc>0","goff"); sh->Draw("z>>z2","","goff");
  z1->Scale(1./z1->Integral()); z2->Scale(1./z2->Integral());
  z1->SetLineColor(kBlue+1); z2->SetLineColor(kMagenta+1); z1->SetLineWidth(2); z2->SetLineWidth(2);
  z1->SetMaximum(std::max(z1->GetMaximum(),z2->GetMaximum())*2.5); z1->SetMinimum(1e-6); z1->SetStats(0);
  z1->Draw("HIST"); z2->Draw("HIST SAME");
  TLegend L(0.35,0.14,0.89,0.27); L.SetBorderSize(0); L.SetFillStyle(0);
  L.AddEntry(z1,"REAL","l"); L.AddEntry(z2,"SIM v3.5revC","l"); L.Draw();
  c.SaveAs("/home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/z_bridge.png");
  printf("saved z_bridge.png (pAu v3.5revC)\n");
}
