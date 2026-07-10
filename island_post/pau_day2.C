#include <map>
#include <vector>
#include <algorithm>
void pau_day2(){
  gROOT->SetBatch(1); gStyle->SetOptStat(0); gStyle->SetTitleFontSize(0.05);
  TFile*fr=TFile::Open("ref_real.root");
  TH1D*r_adc=(TH1D*)fr->Get("r_adc"); TH1D*r_adcz=(TH1D*)fr->Get("r_adcz");
  TH1D*r_run=(TH1D*)fr->Get("r_run"); TH1D*r_lay=(TH1D*)fr->Get("r_lay");
  // digi side (single pass like make_ref)
  TFile*f=TFile::Open("digi_frames_production_v2.root"); TTree*h=(TTree*)f->Get("ntp_hit");
  float event,layer,phibin,tb,adc,side;
  h->SetBranchStatus("*",0);
  for(const char* b : {"event","layer","phibin","zbin","adc","zelem"}) h->SetBranchStatus(b,1);
  h->SetBranchAddress("event",&event); h->SetBranchAddress("layer",&layer);
  h->SetBranchAddress("phibin",&phibin); h->SetBranchAddress("zbin",&tb);
  h->SetBranchAddress("adc",&adc); h->SetBranchAddress("zelem",&side);
  TH1D* d_adc=new TH1D("d_adc","",1101,-0.5,1100.5); TH1D* d_adcz=new TH1D("d_adcz","",121,-0.5,120.5);
  TH1D* d_run=new TH1D("d_run","",25,0.5,25.5); TH1D* d_lay=new TH1D("d_lay","",48,6.5,54.5);
  std::map<uint64_t,std::vector<std::pair<int,int>>> cols;
  for(Long64_t i=0;i<h->GetEntries();++i){ h->GetEntry(i);
    d_adc->Fill(adc); d_adcz->Fill(adc); d_lay->Fill(layer);
    uint64_t k=((uint64_t)(uint32_t)event<<24U)|((uint64_t)(uint32_t)layer<<8U)|(uint64_t)(((int)side==1)?1:0);
    cols[k].push_back({(int)phibin,(int)tb}); }
  for(auto&g:cols){ auto&v=g.second; std::sort(v.begin(),v.end());
    int rl=1;
    for(size_t i=1;i<=v.size();++i){
      if(i<v.size() && v[i].first==v[i-1].first && v[i].second==v[i-1].second+1) rl++;
      else { d_run->Fill(rl); rl=1; } } }
  d_lay->Scale(1./250.);
  // islands
  TFile*ir=TFile::Open("island_real.root"); TTree*tir=(TTree*)ir->Get("island");
  TFile*id=TFile::Open("island_frames_v2.root"); TTree*tid=(TTree*)id->Get("island");
  printf("== pau_day2 (day-2 harness, frozen anchors) ==\n");
  printf("pixel:  mean ADC  real %.2f | digi %.2f\n", r_adc->GetMean(), d_adc->GetMean());
  printf("run:    mean %.2f/%.2f  <=3 frac %.3f/%.3f (real/digi)\n",
    r_run->GetMean(), d_run->GetMean(),
    r_run->Integral(1,3)/r_run->Integral(), d_run->Integral(1,3)/d_run->Integral());
  auto m=[&](TTree*t,const char*v,double lo,double hi){TH1D hh("hh","",3000,lo,hi);t->Draw(Form("%s>>hh",v),"","goff");return hh.GetMean();};
  printf("island: size %.1f/%.1f  phisize %.2f/%.2f  zsize %.2f/%.2f  asym %+.3f/%+.3f (real/digi)\n",
    m(tir,"size",0,3000),m(tid,"size",0,3000), m(tir,"phisize",0,300),m(tid,"phisize",0,300),
    m(tir,"zsize",0,1000),m(tid,"zsize",0,1000), m(tir,"asym",-1,1),m(tid,"asym",-1,1));
  printf("occupancy: pixels/event real 257k | digi %.0fk (x%.1f)\n",
    d_adc->GetEntries()/250/1000., d_adc->GetEntries()/250./257079.);

  auto d2=[&](TH1*a,TH1*b,const char*xt,const char*ti,bool logy,bool perEv){
    if(!perEv){ if(a->Integral()>0)a->Scale(1./a->Integral()); if(b->Integral()>0)b->Scale(1./b->Integral()); }
    a->SetLineColor(kBlue+1);a->SetLineWidth(2); b->SetLineColor(kMagenta+1);b->SetLineWidth(2);
    a->SetStats(0);b->SetStats(0);
    double mx=std::max(a->GetMaximum(),b->GetMaximum()); a->SetMaximum(logy?mx*2:mx*1.3);
    if(logy){gPad->SetLogy(); a->SetMinimum(1e-7);} else a->SetMinimum(0);
    a->SetTitle(ti); a->GetXaxis()->SetTitle(xt); a->GetYaxis()->SetTitle(perEv?"per event":"norm.");
    a->Draw("HIST"); b->Draw("HIST SAME");
    TLegend*L=new TLegend(0.38,0.75,0.89,0.89);L->SetBorderSize(0);L->SetFillStyle(0);
    L->AddEntry(a,"REAL","l"); L->AddEntry(b,"SIM pAu frames v2 (B3)","l"); L->Draw(); };
  TCanvas c("c","",1500,950); c.Divide(3,2);
  c.cd(1); d2((TH1*)r_adc->Clone(),(TH1*)d_adc->Clone(),"per-pixel ADC","pixel ADC spectrum",true,false);
  c.cd(2); d2((TH1*)r_adcz->Clone(),(TH1*)d_adcz->Clone(),"per-pixel ADC","near threshold (ZS region)",true,false);
  c.cd(3); d2((TH1*)r_run->Clone(),(TH1*)d_run->Clone(),"consecutive tbins/pad","run length",true,false);
  c.cd(4); d2((TH1*)r_lay->Clone(),(TH1*)d_lay->Clone(),"TPC layer","hits/event/layer",true,true);
  // island panels
  auto H=[&](TTree*t,const char*v,const char*hn,int nb,double lo,double hi){TH1D*x=new TH1D(hn,"",nb,lo,hi);t->Draw(Form("%s>>%s",v,hn),"","goff");return x;};
  c.cd(5); d2(H(tir,"size","i1",100,0.5,100.5),H(tid,"size","i2",100,0.5,100.5),"pixels/cluster","island size (same algo)",true,false);
  c.cd(6); d2(H(tir,"zsize","i3",40,0.5,40.5),H(tid,"zsize","i4",40,0.5,40.5),"z-size [tbins]","island z-size (same algo)",true,false);
  c.SaveAs("/home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/pau_day2.png");
  printf("saved pau_day2.png\n");
}
