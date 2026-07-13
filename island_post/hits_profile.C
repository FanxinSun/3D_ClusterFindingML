// hit-level profiler: the ntp_hit deep-dive (source data for detached clustering)
#include <map>
#include <vector>
#include <algorithm>
struct HP { TH1D *run, *adc, *adcz, *lay, *fold; double nev; };
HP profile(const char* in, bool isSim, const char* tag){
  TFile*f=TFile::Open(in); TTree*h=(TTree*)f->Get("ntp_hit");
  float event,layer,phibin,tb,adc,side,phi;
  h->SetBranchStatus("*",0);
  for(const char* b : {"event","layer","phibin","adc","zelem","phi"}) h->SetBranchStatus(b,1);
  h->SetBranchStatus(isSim?"zbin":"tbin",1);
  h->SetBranchAddress("event",&event); h->SetBranchAddress("layer",&layer);
  h->SetBranchAddress("phibin",&phibin); h->SetBranchAddress(isSim?"zbin":"tbin",&tb);
  h->SetBranchAddress("adc",&adc); h->SetBranchAddress("zelem",&side); h->SetBranchAddress("phi",&phi);
  HP o; o.run=new TH1D(Form("run_%s",tag),"",25,0.5,25.5); o.adc=new TH1D(Form("adc_%s",tag),"",220,0,1100);
  o.adcz=new TH1D(Form("adcz_%s",tag),"",120,0,120); o.lay=new TH1D(Form("lay_%s",tag),"",48,6.5,54.5);
  o.fold=new TH1D(Form("fold_%s",tag),"",60,0,0.5236); // phi mod 30deg sector period
  std::map<uint64_t,std::vector<std::pair<int,int>>> cols; // (evt,layer,side) -> (pad,tbin)
  double evmin=1e9,evmax=-1e9;
  Long64_t N=h->GetEntries();
  for(Long64_t i=0;i<N;++i){ h->GetEntry(i);
    if(layer<7||layer>54||adc<=0) continue;
    evmin=std::min(evmin,(double)event); evmax=std::max(evmax,(double)event);
    o.adc->Fill(adc); o.adcz->Fill(adc); o.lay->Fill(layer);
    double ph=phi; while(ph<0)ph+=2*M_PI; o.fold->Fill(fmod(ph,0.5235988));
    int sd=((int)side==1)?1:0;
    uint64_t k=((uint64_t)(uint32_t)event<<24U)|((uint64_t)(uint32_t)layer<<8U)|(uint64_t)sd;
    cols[k].push_back({(int)phibin,(int)tb});
  }
  for(auto&g:cols){ auto&v=g.second; std::sort(v.begin(),v.end());
    int rl=1;
    for(size_t i=1;i<=v.size();++i){
      if(i<v.size() && v[i].first==v[i-1].first && v[i].second==v[i-1].second+1) rl++;
      else { o.run->Fill(rl); rl=1; }
    } }
  o.nev=evmax-evmin+1;
  printf("%s: %lld raw rows, %.0f events, <runlen>=%.2f  run==5 frac=%.3f  run<=3 frac=%.3f\n",
    tag,N,o.nev,o.run->GetMean(),o.run->GetBinContent(5)/std::max(1.,o.run->Integral()),
    o.run->Integral(1,3)/std::max(1.,o.run->Integral()));
  return o;
}
void hits_profile(){
  gROOT->SetBatch(1); gStyle->SetOptStat(0); gStyle->SetTitleFontSize(0.05);
  HP r=profile("/home/rog/sPHENIX/3D_ClusterFindingML/clusters_seeds_island_79507-0.root_ntuplizer.root",false,"real");
  HP s=profile("/home/rog/sPHENIX/3D_ClusterFindingML/macros-offline/detectors/sPHENIX/exam5_g4svtx_eval.root",true,"sim");
  auto d2=[&](TH1*a,TH1*b,const char*xt,const char*ti,bool logy,bool scaleEv){
    if(scaleEv){a->Scale(1./r.nev); b->Scale(1./s.nev);} else {if(a->Integral()>0)a->Scale(1./a->Integral()); if(b->Integral()>0)b->Scale(1./b->Integral());}
    a->SetLineColor(kBlue+1);a->SetLineWidth(2); b->SetLineColor(kRed+1);b->SetLineWidth(2);
    double mx=std::max(a->GetMaximum(),b->GetMaximum()); a->SetMaximum(logy?mx*2:mx*1.3);
    if(logy){gPad->SetLogy(); a->SetMinimum(1e-7);} else a->SetMinimum(0);
    a->SetTitle(ti); a->GetXaxis()->SetTitle(xt); a->GetYaxis()->SetTitle(scaleEv?"per event":"norm.");
    a->Draw("HIST"); b->Draw("HIST SAME");
    TLegend*L=new TLegend(0.45,0.76,0.89,0.89);L->SetBorderSize(0);L->SetFillStyle(0);
    L->AddEntry(a,"REAL ntp_hit","l"); L->AddEntry(b,"SIM ntp_hit (exam5)","l"); L->Draw(); };
  TCanvas c("c","",1500,950); c.Divide(3,2);
  c.cd(1); d2((TH1*)r.run->Clone(),(TH1*)s.run->Clone(),"consecutive tbins per pad","time-column run length  [KEY: 5-bin comb]",true,false);
  c.cd(2); d2(r.run,s.run,"consecutive tbins per pad","run length (linear zoom)",false,false); r.run->GetXaxis()->SetRangeUser(0.5,10.5);
  c.cd(3); d2(r.adc,s.adc,"per-hit ADC","per-hit ADC spectrum",true,false);
  c.cd(4); d2(r.adcz,s.adcz,"per-hit ADC","ADC near threshold (ZS region)",true,false);
  c.cd(5); d2(r.lay,s.lay,"TPC layer","hits/event per layer  [occupancy]",true,true);
  c.cd(6); d2(r.fold,s.fold,"phi mod 30#circ [rad]","sector-gap structure",false,false);
  c.SaveAs("/home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/hit_profile.png");
  printf("saved hit_profile.png\n");
}

// exam6/pAu replot (2026-07-10): same 6 panels, sim = composed v3.6 (P0-P3 bridge).
// The original hits_profile() above is kept as the AuAu/exam5-era record (5-bin comb discovery).
void hits_profile_pau(){
  gROOT->SetBatch(1); gStyle->SetOptStat(0); gStyle->SetTitleFontSize(0.05);
  HP r=profile("/home/rog/sPHENIX/3D_ClusterFindingML/clusters_seeds_island_79507-0.root_ntuplizer.root",false,"real");
  HP s4=profile("digi_frames_production_v34.root",true,"v34");
  HP s=profile("digi_frames_production_v36.root",true,"v36");
  auto d2=[&](TH1*a,TH1*g,TH1*b,const char*xt,const char*ti,bool logy,bool scaleEv){
    if(scaleEv){a->Scale(1./r.nev); g->Scale(1./s4.nev); b->Scale(1./s.nev);}
    else {for(TH1*x:{a,g,b}) if(x->Integral()>0)x->Scale(1./x->Integral());}
    a->SetLineColor(kBlue+1);a->SetLineWidth(2);
    g->SetLineColor(kGreen+2);g->SetLineWidth(2);
    b->SetLineColor(kMagenta+1);b->SetLineWidth(2);
    double mx=std::max({a->GetMaximum(),g->GetMaximum(),b->GetMaximum()}); a->SetMaximum(logy?mx*2:mx*1.3);
    if(logy){gPad->SetLogy(); a->SetMinimum(1e-7);} else a->SetMinimum(0);
    a->SetTitle(ti); a->GetXaxis()->SetTitle(xt); a->GetYaxis()->SetTitle(scaleEv?"per frame":"norm.");
    a->Draw("HIST"); g->Draw("HIST SAME"); b->Draw("HIST SAME");
    TLegend*L=new TLegend(0.42,0.70,0.89,0.89);L->SetBorderSize(0);L->SetFillStyle(0);
    L->AddEntry(a,"REAL ntp_hit","l"); L->AddEntry(g,"SIM v3.4 (B4.3)","l");
    L->AddEntry(b,"SIM v3.6","l"); L->Draw(); };
  TCanvas c("c","",1500,950); c.Divide(3,2);
  c.cd(1); d2((TH1*)r.run->Clone(),(TH1*)s4.run->Clone(),(TH1*)s.run->Clone(),"consecutive tbins per pad","time-column run length",true,false);
  c.cd(2); d2(r.run,s4.run,s.run,"consecutive tbins per pad","run length (linear zoom)",false,false); r.run->GetXaxis()->SetRangeUser(0.5,10.5);
  c.cd(3); d2(r.adc,s4.adc,s.adc,"per-hit ADC","per-hit ADC spectrum",true,false);
  c.cd(4); d2(r.adcz,s4.adcz,s.adcz,"per-hit ADC","ADC near threshold (ZS region)",true,false);
  c.cd(5); d2(r.lay,s4.lay,s.lay,"TPC layer","hits/frame per layer  [occupancy]",true,true);
  c.cd(6); d2(r.fold,s4.fold,s.fold,"phi mod 30#circ [rad]","sector-gap structure",false,false);
  c.SaveAs("/home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/hit_profile.png");
  printf("saved hit_profile.png (v3.4 vs v3.6)\n");
}
