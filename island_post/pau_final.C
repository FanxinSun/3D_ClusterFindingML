void pau_final(){
  gROOT->SetBatch(1); gStyle->SetOptStat(0); gStyle->SetTitleFontSize(0.055);
  const char* RF="/home/rog/sPHENIX/3D_ClusterFindingML/clusters_seeds_island_79507-0.root_ntuplizer.root";
  TFile*fr=TFile::Open(RF); TTree*rh=(TTree*)fr->Get("ntp_hit");
  TFile*fs=TFile::Open("digi_frames_production_v35.root"); TTree*sh=(TTree*)fs->Get("ntp_hit");
  TFile*fri=TFile::Open("island91_real.root"); TTree*ri=(TTree*)fri->Get("ntp_cluster");
  TFile*fsi=TFile::Open("island91_frames_production_v35.root"); TTree*si=(TTree*)fsi->Get("ntp_cluster");
  TTree*st=(TTree*)fsi->Get("ntp_truth");
  auto d2=[&](TVirtualPad*p,TH1D*a,TH1D*b,const char*ti,bool logy){
    p->cd(); if(a->Integral()>0)a->Scale(1./a->Integral()); if(b->Integral()>0)b->Scale(1./b->Integral());
    a->SetLineColor(kBlue+1);a->SetLineWidth(2); b->SetLineColor(kMagenta+1);b->SetLineWidth(2);
    a->SetStats(0);b->SetStats(0); a->SetTitle(ti); a->GetYaxis()->SetTitleOffset(1.4);
    double mx=std::max(a->GetMaximum(),b->GetMaximum()); a->SetMaximum(logy?mx*2.5:mx*1.3);
    if(logy){gPad->SetLogy();a->SetMinimum(1e-7);} else a->SetMinimum(0);
    a->Draw("HIST"); b->Draw("HIST SAME");
    TLegend*L=new TLegend(0.42,0.75,0.89,0.89); L->SetBorderSize(0); L->SetFillStyle(0);
    L->AddEntry(a,"REAL","l"); L->AddEntry(b,"SIM v3.5revC","l"); L->Draw(); };
  auto H=[&](TTree*x,const char*v,const char*hn,int nb,double lo,double hi,const char*cut=""){
    TH1D*hh=new TH1D(hn,"",nb,lo,hi); x->Draw(Form("%s>>%s",v,hn),cut,"goff"); return hh; };

  // ---- pixel-level figure ----
  // real pixel-level curves = the FROZEN ANCHORS (ref_real.root: layer 7-54 && adc>0),
  // so they are identical to the day-2 era figures BY CONSTRUCTION.
  TFile*fref=TFile::Open("ref_real.root");
  TH1D*RA=(TH1D*)((TH1D*)fref->Get("r_adc"))->Clone("RA");   RA->SetDirectory(nullptr);
  TH1D*RZ=(TH1D*)((TH1D*)fref->Get("r_adcz"))->Clone("RZ");  RZ->SetDirectory(nullptr);
  TH1D*RR=(TH1D*)((TH1D*)fref->Get("r_run"))->Clone("RR");   RR->SetDirectory(nullptr);
  TH1D*RT=(TH1D*)((TH1D*)fref->Get("r_tbin"))->Clone("RT");  RT->SetDirectory(nullptr);
  TH1D*RL=(TH1D*)((TH1D*)fref->Get("r_lay"))->Clone("RL");   RL->SetDirectory(nullptr);
  TCanvas c1("c1","",1600,900); c1.Divide(3,2);
  d2(c1.cd(1),RA,H(sh,"adc","b1",1101,-0.5,1100.5),"per-pixel ADC;ADC;norm",true);
  d2(c1.cd(2),RZ,H(sh,"adc","b2",121,-0.5,120.5),"near threshold (ZS);ADC;norm",true);
  d2(c1.cd(3),RT,H(sh,"zbin","b3",1000,0,1000),"arrivals (adc>0);tbin;norm",false);
  // run lengths (consecutive tbins per pad)
  auto runs=[&](TTree*t,const char*tb,const char*hn){
    TH1D*hr=new TH1D(hn,"",25,0.5,25.5);
    float ev,l,pb,z,sd; t->SetBranchStatus("*",0);
    for(const char*b:{"event","layer","phibin","zelem"}) t->SetBranchStatus(b,1);
    t->SetBranchStatus(tb,1);
    t->SetBranchAddress("event",&ev); t->SetBranchAddress("layer",&l);
    t->SetBranchAddress("phibin",&pb); t->SetBranchAddress(tb,&z); t->SetBranchAddress("zelem",&sd);
    Long64_t N=t->GetEntries();
    std::vector<std::pair<uint64_t,int>> v; v.reserve(N);
    for(Long64_t i=0;i<N;++i){t->GetEntry(i);
      if(l<7||l>54) continue;
      v.push_back({((uint64_t)(uint32_t)ev<<40)|((uint64_t)(uint32_t)l<<32)|((uint64_t)(uint32_t)sd<<24)|((uint64_t)(uint32_t)pb<<8),(int)z});}
    std::sort(v.begin(),v.end(),[](auto&a,auto&b){return a.first==b.first?a.second<b.second:a.first<b.first;});
    int run=0;
    for(size_t i=0;i<v.size();++i){
      if(i&&v[i].first==v[i-1].first&&v[i].second==v[i-1].second+1) run++;
      else{ if(run)hr->Fill(run); run=1;}
    }
    if(run)hr->Fill(run);
    t->SetBranchStatus("*",1);
    return hr; };
  d2(c1.cd(4),RR,runs(sh,"zbin","b4"),"run length;consecutive tbins/pad;norm",true);
  // per-frame pixel count
  auto perev=[&](TTree*t,const char*hn,int nev,const char*cut){
    TH2D htmp("htmp","",1000,0,1000,1,0,1e9); t->Draw(Form("0.5:event>>htmp"),cut,"goff");
    TH1D*hh=new TH1D(hn,"",30,0,600);
    for(int e=0;e<1000;++e){double x=htmp.GetBinContent(e+1,1); if(x>0)hh->Fill(x/1e3);}
    return hh; };
  d2(c1.cd(5),perev(rh,"a5",100,"layer>=7&&layer<=54&&adc>0"),perev(sh,"b5",250,""),"pixels per frame;k pixels;frames (norm)",false);
  d2(c1.cd(6),RL,H(sh,"layer","b6",48,6.5,54.5),"pixels/layer;TPC layer;norm",false);
  c1.SaveAs("/home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/pau_day2.png");

  // ---- island-level figure ----
  TCanvas c2("c2","",1600,900); c2.Divide(3,2);
  d2(c2.cd(1),H(ri,"size","i1",100,0.5,100.5),H(si,"size","j1",100,0.5,100.5),"island size;pixels;norm",true);
  d2(c2.cd(2),H(ri,"phisize","i2",25,0.5,25.5),H(si,"phisize","j2",25,0.5,25.5),"island #phi-size;pads;norm",true);
  d2(c2.cd(3),H(ri,"zsize","i3",40,0.5,40.5),H(si,"zsize","j3",40,0.5,40.5),"island z-size;tbins;norm",true);
  d2(c2.cd(4),H(ri,"adc","i4",150,0,6000),H(si,"adc","j4",150,0,6000),"island ADC;raw sum;norm",true);
  d2(c2.cd(5),H(ri,"eta","i5",60,-2.5,2.5),H(si,"eta","j5",60,-2.5,2.5),"island #eta (apparent);#eta;norm",false);
  TH1D*i6=H(ri,"layer","i6",48,6.5,54.5); i6->Scale(1./100.);
  TH1D*j6=H(si,"layer","j6",48,6.5,54.5); j6->Scale(1./250.);
  c2.cd(6); gPad->SetLogy(); i6->SetLineColor(kBlue+1); j6->SetLineColor(kMagenta+1);
  i6->SetLineWidth(2); j6->SetLineWidth(2); i6->SetStats(0);
  i6->SetTitle("islands/frame/layer;TPC layer;per frame"); i6->SetMinimum(50);
  i6->SetMaximum(std::max(i6->GetMaximum(),j6->GetMaximum())*3);
  i6->Draw("HIST"); j6->Draw("HIST SAME");
  TLegend L6(0.42,0.75,0.89,0.89); L6.SetBorderSize(0); L6.SetFillStyle(0);
  L6.AddEntry(i6,"REAL","l"); L6.AddEntry(j6,"SIM v3.5revC","l"); L6.Draw();
  c2.SaveAs("/home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/pau_islands.png");

  // ---- residual table ----
  printf("\n== CORRECTED residual table (real vs v3.5revC @ 390 kHz) ==\n");
  TH1D ht("ht","",971,0,971); sh->Draw("zbin>>ht","","goff"); ht.Scale(1./250.);
  auto avg=[&](TH1D&hh,int a,int b){double s=0;for(int i=a;i<=b;++i)s+=hh.GetBinContent(i);return s/(b-a+1);};
  TH2D h2("h2","",250,0,250,1,0,1e9); sh->Draw("0.5:event>>h2","","goff");
  double s=0,s2=0;int n=0; for(int e=0;e<250;++e){double x=h2.GetBinContent(e+1,1); if(x>0){s+=x;s2+=x*x;n++;}}
  double m=s/n, sd=sqrt(s2/n-m*m);
  printf("px/frame:   real 257k | sim %.0fk\n", m/1e3);
  printf("pedestal:   real 4429 | sim %.0f px/us\n", avg(ht,400,900)/0.053);
  printf("early/late: real ~1.3 | sim %.2f\n", avg(ht,20,240)/avg(ht,400,900));
  printf("sigma/mu:   real 0.45 | sim %.2f\n", sd/m);
  auto mm=[&](TTree*t,const char*v,double hi){TH1D h("h","",4000,0,hi);t->Draw(Form("%s>>h",v),"","goff");return h.GetMean();};
  printf("pixmean:    real %.2f (anchor) | sim %.2f\n", RA->GetMean(), mm(sh,"adc",1101));
  printf("island size/phi/z/adc: real %.1f/%.2f/%.2f/%.0f | sim %.1f/%.2f/%.2f/%.0f\n",
    mm(ri,"size",3000),mm(ri,"phisize",300),mm(ri,"zsize",1000),mm(ri,"adc",3e5),
    mm(si,"size",3000),mm(si,"phisize",300),mm(si,"zsize",1000),mm(si,"adc",3e5));
  printf("islands/frame: real %.1fk | sim %.1fk\n", ri->GetEntries()/100./1e3, si->GetEntries()/250./1e3);
  double N=st->GetEntries();
  TH1D hp("hp","",100,0,1.0001); st->Draw("purity>>hp","","goff");
  printf("DATASET v2: %.2fM clusters | track %.1f%% looper %.1f%% | purity %.3f | merged %.1f%%\n",
    N/1e6, 100.*st->GetEntries("cls==0")/N, 100.*st->GetEntries("cls==1")/N, hp.GetMean(), 100.*st->GetEntries("ntrks>=2")/N);
}
