// RETIRED 2026-08-17: the "REAL" flash curves in this macro are event 44 alone
// (laser_flash_remodel_request.md); the spike-ratio gate and the "44 flashes" logic
// no longer apply. Kept for provenance only.
// laser_assess.C — assessment figure: real flash vs collaboration models through the
// detached pipeline (CM injected in frames; DL isolated signature).
void laser_assess(){
  gROOT->SetBatch(1); gStyle->SetOptStat(0); gStyle->SetTitleFontSize(0.055);
  TFile*fr=TFile::Open("/home/rog/sPHENIX/3D_ClusterFindingML/clusters_seeds_island_79507-0.root_ntuplizer.root");
  TTree*rh=(TTree*)fr->Get("ntp_hit");
  TFile*fc=TFile::Open("digi_cmtest.root");  TTree*ch=(TTree*)fc->Get("ntp_hit");
  TFile*fd=TFile::Open("digi_dl.root");      TTree*dh=(TTree*)fd->Get("ntp_hit");
  TCanvas c("c","",1700,520); c.Divide(3,1);
  // (1) arrivals: real vs CM-injected frames vs DL (shape)
  c.cd(1);
  TH1D*a=new TH1D("a","arrivals;tbin;norm",971,0,971);
  TH1D*b=new TH1D("b","",971,0,971);
  TH1D*d=new TH1D("d","",971,0,971);
  rh->Draw("tbin>>a","layer>=7&&layer<=54&&adc>0","goff");
  ch->Draw("zbin>>b","","goff");
  dh->Draw("zbin>>d","","goff");
  a->Scale(1./a->Integral()); b->Scale(1./b->Integral()); d->Scale(1./d->Integral());
  a->SetLineColor(kBlue+1); b->SetLineColor(kMagenta+1);
  d->SetLineColor(kGray+2); d->SetLineStyle(2); d->SetLineWidth(1);
  a->SetLineWidth(2); b->SetLineWidth(2); a->SetStats(0);
  gPad->SetLogy(); a->SetMaximum(std::max(d->GetMaximum(),std::max(a->GetMaximum(),b->GetMaximum()))*3); a->SetMinimum(1e-5);
  a->Draw("HIST"); b->Draw("HIST SAME"); d->Draw("HIST SAME");
  TLegend*L1=new TLegend(0.32,0.68,0.89,0.89); L1->SetBorderSize(0); L1->SetFillStyle(0);
  L1->AddEntry(a,"REAL (flash @329)","l");
  L1->AddEntry(b,"SIM frames + CM model (PHG4TpcCentralMembrane)","l");
  L1->AddEntry(d,"DirectLaser model (isolated; endcap-peaked)","l"); L1->Draw();
  // (2) zoom on the spike
  c.cd(2);
  TH1D*a2=(TH1D*)a->Clone("a2"); TH1D*b2=(TH1D*)b->Clone("b2");
  a2->SetMaximum(-1111); b2->SetMaximum(-1111);  // clones inherit panel-1 SetMaximum
  a2->SetTitle("flash window zoom;tbin;norm");
  a2->GetXaxis()->SetRangeUser(300,360); gPad->SetLogy(0);
  double wmax=0; for(int i=301;i<=360;++i) wmax=std::max({wmax,a2->GetBinContent(i),b2->GetBinContent(i)});
  a2->SetMaximum(wmax*1.3); a2->SetMinimum(0);
  a2->Draw("HIST"); b2->Draw("HIST SAME");
  TLegend*L2=new TLegend(0.4,0.76,0.89,0.89); L2->SetBorderSize(0); L2->SetFillStyle(0);
  L2->AddEntry(a2,"REAL","l"); L2->AddEntry(b2,"SIM + CM model","l"); L2->Draw();
  // (3) phi-fold stripes in the window
  c.cd(3);
  TH1D*p1=new TH1D("p1","#phi-fold in flash window (R2);#phi mod 30#circ [rad];norm",300,0,0.5236);
  TH1D*p2=new TH1D("p2","",300,0,0.5236);
  rh->Draw("fmod(phi+2*pi,pi/6)>>p1","layer>=23&&layer<=38&&adc>0&&tbin>=322&&tbin<=340","goff");
  ch->Draw("fmod(phi+2*pi,pi/6)>>p2","layer>=23&&layer<=38&&zbin>=322&&zbin<=340","goff");
  p1->Scale(1./p1->Integral()); p2->Scale(1./p2->Integral());
  p1->SetLineColor(kBlue+1); p2->SetLineColor(kMagenta+1); p1->SetLineWidth(2); p2->SetLineWidth(2);
  p1->SetStats(0); p1->SetMaximum(std::max(p1->GetMaximum(),p2->GetMaximum())*1.3);
  p2->SetStats(0); p2->SetMaximum(std::max(p1->GetMaximum(),p2->GetMaximum())*1.3);
  p2->SetTitle("#phi-fold in flash window (R2);#phi mod 30#circ [rad];norm");
  p2->Draw("HIST"); p1->Draw("HIST SAME");
  TLegend*L3=new TLegend(0.35,0.76,0.89,0.89); L3->SetBorderSize(0); L3->SetFillStyle(0);
  L3->AddEntry(p1,"REAL (smeared stripes)","l"); L3->AddEntry(p2,"SIM + CM model (sharp stripes)","l"); L3->Draw();
  c.SaveAs("/home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/laser_assessment.png");
  printf("saved laser_assessment.png\n");
}
