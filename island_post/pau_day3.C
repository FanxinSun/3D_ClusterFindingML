void pau_day3(){
  gROOT->SetBatch(1); gStyle->SetOptStat(0); gStyle->SetTitleFontSize(0.055);
  TFile*fr=TFile::Open("island91_real.root"); TTree*r=(TTree*)fr->Get("ntp_cluster");
  TFile*fs=TFile::Open("island91_frames_production_v35.root");  TTree*s=(TTree*)fs->Get("ntp_cluster");
  TTree*t=(TTree*)fs->Get("ntp_truth"); s->AddFriend(t,"t");
  auto d2=[&](TH1*a,TH1*b,const char*xt,const char*ti,bool logy,bool perEv){
    if(perEv){ a->Scale(1./100.); b->Scale(1./250.); }
    else { if(a->Integral()>0)a->Scale(1./a->Integral()); if(b->Integral()>0)b->Scale(1./b->Integral()); }
    a->SetLineColor(kBlue+1);a->SetLineWidth(2); b->SetLineColor(kMagenta+1);b->SetLineWidth(2);
    a->SetStats(0);b->SetStats(0);
    double mx=std::max(a->GetMaximum(),b->GetMaximum()); a->SetMaximum(logy?mx*2.5:mx*1.35);
    if(logy){gPad->SetLogy(); a->SetMinimum(1e-7);} else a->SetMinimum(0);
    a->SetTitle(ti); a->GetYaxis()->SetTitleOffset(1.4); a->GetXaxis()->SetTitle(xt); a->GetYaxis()->SetTitle(perEv?"per event":"norm.");
    a->Draw("HIST"); b->Draw("HIST SAME");
    TLegend*L=new TLegend(0.38,0.14,0.89,0.27);L->SetBorderSize(0);L->SetFillStyle(0);
    L->AddEntry(a,"REAL island91","l"); L->AddEntry(b,"SIM v3.5revC","l"); L->Draw(); };
  auto H=[&](TTree*x,const char*v,const char*hn,int nb,double lo,double hi,const char*cut=""){
    TH1D*h=new TH1D(hn,"",nb,lo,hi); x->Draw(Form("%s>>%s",v,hn),cut,"goff"); return h; };

  TCanvas c("c","",1600,1350); c.Divide(3,3);
  c.cd(1); d2(H(r,"adc","a1",150,0,6000),H(s,"adc","b1",150,0,6000),"raw ADC sum","cluster ADC",true,false);
  c.cd(2); d2(H(r,"size","a2",100,0.5,100.5),H(s,"size","b2",100,0.5,100.5),"pixels","cluster size",true,false);
  c.cd(3); d2(H(r,"phisize","a3",25,0.5,25.5),H(s,"phisize","b3",25,0.5,25.5),"#phi-size [pads]","cluster #phi-size",true,false);
  c.cd(4); d2(H(r,"zsize","a4",40,0.5,40.5),H(s,"zsize","b4",40,0.5,40.5),"z-size [tbins]","cluster z-size",true,false);
  c.cd(5); d2(H(r,"r","a5",60,28,80),H(s,"r","b5",60,28,80),"r [cm]","cluster radius",false,false);
  c.cd(6); d2(H(r,"z","a6",100,-320,320),H(s,"z","b6",100,-320,320),"z [cm]","cluster z (frame span)",true,false);
  c.cd(7); d2(H(r,"eta","a7",60,-2.5,2.5),H(s,"eta","b7",60,-2.5,2.5),"#eta","cluster #eta",false,false);
  c.cd(8); d2(H(r,"phielem","a8",12,-0.5,11.5),H(s,"phielem","b8",12,-0.5,11.5),"TPC sector","sector occupancy [dead-map check]",false,false);
  c.cd(9); d2(H(r,"layer","a9",48,6.5,54.5),H(s,"layer","b9",48,6.5,54.5),"TPC layer","clusters/frame/layer",true,true);
  c.SaveAs("/home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/pau_day3.png");

  printf("\n== day-3 residual table (island91 real vs sim) ==\n");
  auto m=[&](TTree*x,const char*v,double lo,double hi){TH1D h("h","",4000,lo,hi);x->Draw(Form("%s>>h",v),"","goff");return h.GetMean();};
  const char* obs[7] = {"adc","size","phisize","zsize","r","eta","maxadc"};
  double hi[7] = {3e5,3000,300,1000,90,3,66000};
  for(int i=0;i<7;++i){
    double mr=m(r,obs[i],(i==5?-3:0),hi[i]), ms=m(s,obs[i],(i==5?-3:0),hi[i]);
    printf("  %-8s real %8.2f | sim %8.2f | %+5.1f%%\n",obs[i],mr,ms,100.*(ms-mr)/(mr!=0?mr:1));
  }
  printf("  %-8s real %8.0f | sim %8.0f | x%.2f  (per event)\n","islands",r->GetEntries()/100.,s->GetEntries()/250.,(s->GetEntries()/250.)/(r->GetEntries()/100.));
  printf("  z-span   real +-%.0f cm | sim +-%.0f cm  (frame span: sim 971 tbins; real apparent-z reaches further)\n",
    std::max(fabs(r->GetMinimum("z")),r->GetMaximum("z")), std::max(fabs(s->GetMinimum("z")),s->GetMaximum("z")));
  printf("  truth composition (sim): track %.1f%%  looper %.1f%%  merged(ntrks>=2) %.1f%%  <purity> %.3f\n",
    100.*t->GetEntries("cls==0")/t->GetEntries(), 100.*t->GetEntries("cls==1")/t->GetEntries(),
    100.*t->GetEntries("ntrks>=2")/t->GetEntries(),
    ({TH1D h("hp","",100,0,1.0001); t->Draw("purity>>hp","","goff"); h.GetMean();}));
}
