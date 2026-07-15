// cmp_prodclus.C — P4 figure: PRODUCTION-SEMANTICS comparison (new file).
// Three curves per panel:
//   REAL ntp_cluster  (the production data itself)
//   port x REAL pixels (the GATE: same algorithm, same pixels -> must overlay real)
//   port x SIM v4.0 (the first apples-to-apples sim comparison)
void cmp_prodclus()
{
  gROOT->SetBatch(1); gStyle->SetOptStat(0);
  TFile*fr=TFile::Open("/home/rog/sPHENIX/3D_ClusterFindingML/clusters_seeds_island_79507-0.root_ntuplizer.root");
  TTree*rc=(TTree*)fr->Get("ntp_cluster");
  TFile*fg=TFile::Open("prodclus_real_fw0.root"); TTree*gc=(TTree*)fg->Get("ntp_clus");
  TFile*fs=TFile::Open("prodclus_v40b.root");      TTree*sc=(TTree*)fs->Get("ntp_clus");
  const char* RC="layer>=7&&layer<=54";
  auto mk=[&](TTree*t,const char*v,const char*hn,int nb,double lo,double hi,const char*cut){
    TH1D*h=new TH1D(hn,"",nb,lo,hi); t->Draw(Form("%s>>%s",v,hn),cut,"goff");
    if(h->Integral()>0)h->Scale(1./h->Integral()); return h; };
  auto d3=[&](TVirtualPad*p,TH1D*a,TH1D*g,TH1D*b,const char*ti,bool logy){
    p->cd();
    a->SetLineColor(kBlue+1); a->SetLineWidth(2);
    g->SetLineColor(kGray+2); g->SetLineStyle(2); g->SetLineWidth(2);
    b->SetLineColor(kMagenta+1); b->SetLineWidth(2);
    for(TH1D*x:{a,g,b}) x->SetStats(0);
    a->SetTitle(ti);
    double mx=std::max({a->GetMaximum(),g->GetMaximum(),b->GetMaximum()});
    a->SetMaximum(logy?mx*2.5:mx*1.3);
    if(logy){gPad->SetLogy(); a->SetMinimum(1e-7);} else a->SetMinimum(0);
    a->Draw("HIST"); g->Draw("HIST SAME"); b->Draw("HIST SAME");
    TLegend*L=new TLegend(0.35,0.68,0.89,0.89); L->SetBorderSize(0); L->SetFillStyle(0);
    L->AddEntry(a,"REAL ntp_cluster","l");
    L->AddEntry(g,"port #times real pixels  [gate]","l");
    L->AddEntry(b,"port #times SIM v4.0","l"); L->Draw(); };
  TCanvas c("c","",1500,950); c.Divide(2,2);
  d3(c.cd(1), mk(rc,"zsize","a1",45,0.5,45.5,RC), mk(gc,"zsize","g1",45,0.5,45.5,""),
     mk(sc,"zsize","b1",45,0.5,45.5,""), "cluster z-size (production semantics);tbins;norm", true);
  d3(c.cd(2), mk(rc,"phisize","a2",25,0.5,25.5,RC), mk(gc,"phisize","g2",25,0.5,25.5,""),
     mk(sc,"phisize","b2",25,0.5,25.5,""), "cluster #phi-size (production semantics);pads;norm", true);
  d3(c.cd(3), mk(rc,"adc","a3",150,0,6000,RC), mk(gc,"adc","g3",150,0,6000,""),
     mk(sc,"adc","b3",150,0,6000,""), "cluster adc;adc;norm", true);
  d3(c.cd(4), mk(rc,"maxadc","a4",120,0,1080,RC), mk(gc,"maxadc","g4",120,0,1080,""),
     mk(sc,"maxadc","b4",120,0,1080,""), "cluster maxadc;ADC;norm", true);
  c.SaveAs("/home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/cmp_prodclus.png");
  printf("saved cmp_prodclus.png\n");
}
