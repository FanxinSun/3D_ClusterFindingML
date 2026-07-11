// zsize_physics.C — replot of the z-size structure explainer:
//  (1) the production-clusterizer CAP artifact (real ntp_cluster zsize cliff; the old
//      "zsize 12-20 empty" observation) vs uncapped islands on the SAME real data;
//  (2) real vs sim islands zsize (uncapped physics tails);
//  (3) diffusion growth: mean island z-size vs |z| (drift length), real vs sim.
// Sim side follows the CURRENT production (v3.3 B4.2 via island_frames_v33.root).
void zsize_physics(){
  gROOT->SetBatch(1); gStyle->SetOptStat(0);
  TFile*fr=TFile::Open("/home/rog/sPHENIX/3D_ClusterFindingML/clusters_seeds_island_79507-0.root_ntuplizer.root");
  TTree*rc=(TTree*)fr->Get("ntp_cluster");
  TFile*fi=TFile::Open("island_real.root");      TTree*ri=(TTree*)fi->Get("island");
  TFile*fv=TFile::Open("island_frames_v33.root"); TTree*si=(TTree*)fv->Get("island");
  TCanvas c("c","",1800,600); c.Divide(3,1);
  auto mk=[&](TTree*t,const char*v,const char*hn,const char*cut){
    TH1D*h=new TH1D(hn,"",40,0.5,40.5); t->Draw(Form("%s>>%s",v,hn),cut,"goff");
    if(h->Integral()>0)h->Scale(1./h->Integral()); return h; };
  c.cd(1); gPad->SetLogy();
  TH1D*a=mk(rc,"zsize","a","layer>=7&&layer<=54"); TH1D*b=mk(ri,"zsize","b","");
  a->SetLineColor(kBlue+1); b->SetLineColor(kGreen+2); a->SetLineWidth(2); b->SetLineWidth(2);
  a->SetTitle("REAL: production clusterizer vs islands;z-size [tbins];norm"); a->SetStats(0);
  a->SetMaximum(std::max(a->GetMaximum(),b->GetMaximum())*2.5); a->SetMinimum(1e-7);
  a->Draw("HIST"); b->Draw("HIST SAME");
  TLegend*L1=new TLegend(0.35,0.72,0.89,0.89); L1->SetBorderSize(0); L1->SetFillStyle(0);
  L1->AddEntry(a,"production ntp_cluster (CAPPED)","l"); L1->AddEntry(b,"islandize (uncapped)","l"); L1->Draw();
  c.cd(2); gPad->SetLogy();
  TH1D*d=mk(ri,"zsize","d",""); TH1D*e=mk(si,"zsize","e","");
  d->SetLineColor(kBlue+1); e->SetLineColor(kMagenta+1); d->SetLineWidth(2); e->SetLineWidth(2);
  d->SetTitle("islands: real vs sim;z-size [tbins];norm"); d->SetStats(0);
  d->SetMaximum(std::max(d->GetMaximum(),e->GetMaximum())*2.5); d->SetMinimum(1e-7);
  d->Draw("HIST"); e->Draw("HIST SAME");
  TLegend*L2=new TLegend(0.42,0.76,0.89,0.89); L2->SetBorderSize(0); L2->SetFillStyle(0);
  L2->AddEntry(d,"REAL islands","l"); L2->AddEntry(e,"SIM v3.3 (B4.2)","l"); L2->Draw();
  c.cd(3);
  TProfile*p1=new TProfile("p1","<z-size> vs |apparent z|  [flat = correct: #sigma_{L} sub-tbin, t0-mixed];|apparent z| [cm];<z-size> [tbins]",30,0,300);
  TProfile*p2=new TProfile("p2","",30,0,300);
  ri->Draw("zsize:abs(cz)>>p1","","goff"); si->Draw("zsize:abs(cz)>>p2","","goff");
  p1->SetLineColor(kBlue+1); p2->SetLineColor(kMagenta+1); p1->SetLineWidth(2); p2->SetLineWidth(2);
  p1->SetStats(0); p1->SetMinimum(0); p1->SetMaximum(8);
  p1->Draw(); p2->Draw("SAME");
  TLegend*L3=new TLegend(0.42,0.76,0.89,0.89); L3->SetBorderSize(0); L3->SetFillStyle(0);
  L3->AddEntry(p1,"REAL islands","l"); L3->AddEntry(p2,"SIM v3.3 (B4.2)","l"); L3->Draw();
  c.SaveAs("/home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/zsize_physics.png");
  printf("saved zsize_physics.png (v3.3 B4.2)\n");
}
