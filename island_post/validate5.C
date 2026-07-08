void norm(TH1*h){ if(h->Integral()>0) h->Scale(1.0/h->Integral()); }
void validate5(){
  gROOT->SetBatch(1); gStyle->SetOptStat(0); gStyle->SetTitleFontSize(0.05);
  TFile*fp=TFile::Open("/home/rog/sPHENIX/3D_ClusterFindingML/clusters_seeds_island_79507-0.root_ntuplizer.root");
  TFile*fm=TFile::Open("island_real.root");
  TTree*p=(TTree*)fp->Get("ntp_cluster");   // BNL production island clusters
  TTree*m=(TTree*)fm->Get("island");        // my islandize.C on the same raw hits
  const char* tpc="layer>=7&&layer<=54";
  TCanvas c("c","",1300,500); c.Divide(2,1);
  auto d2=[&](TH1*a,TH1*b,const char*xt,const char*ti){
    a->SetLineColor(kBlue+1);a->SetLineWidth(2); b->SetLineColor(kRed+1);b->SetLineWidth(2);
    a->SetStats(0);b->SetStats(0); gPad->SetLogy();
    a->SetMaximum(std::max(a->GetMaximum(),b->GetMaximum())*2); a->SetMinimum(1e-7);
    a->SetTitle(ti); a->GetXaxis()->SetTitle(xt); a->GetYaxis()->SetTitle("norm.");
    a->Draw("HIST"); b->Draw("HIST SAME");
    TLegend*L=new TLegend(0.35,0.74,0.89,0.89);L->SetBorderSize(0);L->SetFillStyle(0);
    L->AddEntry(a,"BNL production ntp_cluster","l");
    L->AddEntry(b,"my islandize.C (same raw hits)","l"); L->Draw(); };
  c.cd(1);
  TH1D*a1=new TH1D("a1","",60,0.5,60.5); p->Draw("size>>a1",Form("%s&&size>0",tpc),"goff"); norm(a1);
  TH1D*b1=new TH1D("b1","",60,0.5,60.5); m->Draw("size>>b1","","goff"); norm(b1);
  d2(a1,b1,"pixels/cluster","REAL data: production vs my islander -- size");
  c.cd(2);
  TH1D*a2=new TH1D("a2","",25,0.5,25.5); p->Draw("zsize>>a2",tpc,"goff"); norm(a2);
  TH1D*b2=new TH1D("b2","",25,0.5,25.5); m->Draw("zsize>>b2","","goff"); norm(b2);
  d2(a2,b2,"z-size [tbins]","REAL data: production vs my islander -- z-size");
  c.SaveAs("/home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/islander_validation.png");
  printf("means: production size=%.2f zsize=%.2f | mine size=%.2f zsize=%.2f\n",
    a1->GetMean(),a2->GetMean(),b1->GetMean(),b2->GetMean());
}
