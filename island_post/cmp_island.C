void norm(TH1*h){ if(h->Integral()>0) h->Scale(1.0/h->Integral()); }
TH1D* H(TTree*t,const char*v,const char*cut,const char*hn,int nb,double lo,double hi){
  TH1D*h=new TH1D(hn,"",nb,lo,hi); t->Draw(Form("%s>>%s",v,hn),cut,"goff"); norm(h); return h; }
void d2(TH1*a,TH1*b,const char*xt,const char*ti,bool logy){
  a->SetLineColor(kBlue+1);a->SetLineWidth(2); b->SetLineColor(kRed+1);b->SetLineWidth(2);
  a->SetStats(0);b->SetStats(0);
  double mx=std::max(a->GetMaximum(),b->GetMaximum()); a->SetMaximum(logy?mx*2:mx*1.3);
  if(logy){gPad->SetLogy(); a->SetMinimum(1e-7);} else a->SetMinimum(0);
  a->SetTitle(ti); a->GetXaxis()->SetTitle(xt); a->GetYaxis()->SetTitle("norm.");
  a->Draw("HIST"); b->Draw("HIST SAME");
  TLegend*L=new TLegend(0.42,0.76,0.89,0.89);L->SetBorderSize(0);L->SetFillStyle(0);
  L->AddEntry(a,"REAL (island post)","l"); L->AddEntry(b,"SIM AuAu (island post)","l"); L->Draw(); }
void cmp_island(){
  gROOT->SetBatch(1); gStyle->SetOptStat(0); gStyle->SetTitleFontSize(0.05);
  TFile*fr=TFile::Open("island_real.root"); TFile*fs=TFile::Open("island_sim.root");
  TTree*r=(TTree*)fr->Get("island"); TTree*s=(TTree*)fs->Get("island");
  auto st=[&](TTree*t,const char*tag){
    auto m=[&](const char*v,double lo,double hi)->double{TH1D h("h","",4000,lo,hi);t->Draw(Form("%s>>h",v),"","goff");return h.GetMean();};
    long n=(long)t->GetEntries();
    printf("  %-4s N=%-8ld <size>=%.1f <phisize>=%.2f <zsize>=%.2f <adc>=%.0f <maxadc>=%.0f <asym>=%+.3f <rho>=%.2f\n",
      tag,n,m("size",0,3000),m("phisize",0,300),m("zsize",0,1000),m("adc",0,3e5),m("maxadc",0,66000),m("asym",-1,1),m("rho",0,1.01));
    printf("       tails: size>20 %.3f  phisize>7 %.4f  zsize>11 %.4f  (of all)\n",
      (double)t->GetEntries("size>20")/n,(double)t->GetEntries("phisize>7")/n,(double)t->GetEntries("zsize>11")/n);
  };
  printf("== island clusters, SAME algorithm both files ==\n"); st(r,"REAL"); st(s,"SIM");
  TCanvas c("c","",1500,950); c.Divide(3,2);
  c.cd(1); d2(H(r,"size","","a1",100,0.5,100.5),H(s,"size","","b1",100,0.5,100.5),"pixels/cluster","island size",true);
  c.cd(2); d2(H(r,"phisize","","a2",30,0.5,30.5),H(s,"phisize","","b2",30,0.5,30.5),"#phi-size [pads]","island #phi-size",true);
  c.cd(3); d2(H(r,"zsize","","a3",40,0.5,40.5),H(s,"zsize","","b3",40,0.5,40.5),"z-size [tbins]","island z-size",true);
  c.cd(4); d2(H(r,"adc","","a4",150,0,6000),H(s,"adc","","b4",150,0,6000),"raw ADC sum","island ADC (RAW sums both)",true);
  c.cd(5); d2(H(r,"asym","asym>-1&&asym<1","a5",60,-1,1),H(s,"asym","asym>-1&&asym<1","b5",60,-1,1),"(late-early)/sum","time-axis ADC asymmetry",false);
  c.cd(6); d2(H(r,"rho","","a6",50,0,1.001),H(s,"rho","","b6",50,0,1.001),"pixel fill of bbox","raggedness (1=solid block)",false);
  c.SaveAs("/home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/island_cmp.png");
  printf("saved island_cmp.png\n");
}
