#include <map>
#include <vector>
#include <algorithm>
void digi_check(){
  gROOT->SetBatch(1); gStyle->SetOptStat(0);
  // --- run-length + per-pixel ADC of the DIGITIZED pixels vs REAL ---
  auto prof=[&](const char* in, bool isSim, TH1D*&run, TH1D*&adc){
    TFile*f=TFile::Open(in); TTree*h=(TTree*)f->Get("ntp_hit");
    float event,layer,phibin,tb,a,side;
    h->SetBranchStatus("*",0);
    for(const char* b : {"event","layer","phibin","adc","zelem"}) h->SetBranchStatus(b,1);
    h->SetBranchStatus(isSim?"zbin":"tbin",1);
    h->SetBranchAddress("event",&event); h->SetBranchAddress("layer",&layer);
    h->SetBranchAddress("phibin",&phibin); h->SetBranchAddress(isSim?"zbin":"tbin",&tb);
    h->SetBranchAddress("adc",&a); h->SetBranchAddress("zelem",&side);
    run=new TH1D(Form("run_%d",(int)isSim),"",25,0.5,25.5);
    adc=new TH1D(Form("adc_%d",(int)isSim),"",220,0,1100);
    std::map<uint64_t,std::vector<std::pair<int,int>>> cols;
    Long64_t N=h->GetEntries();
    for(Long64_t i=0;i<N;++i){ h->GetEntry(i);
      if(layer<7||layer>54||a<=0) continue;
      adc->Fill(a);
      uint64_t k=((uint64_t)(uint32_t)event<<24U)|((uint64_t)(uint32_t)layer<<8U)|(uint64_t)(((int)side==1)?1:0);
      cols[k].push_back({(int)phibin,(int)tb}); }
    for(auto&g:cols){ auto&v=g.second; std::sort(v.begin(),v.end());
      int rl=1;
      for(size_t i=1;i<=v.size();++i){
        if(i<v.size() && v[i].first==v[i-1].first && v[i].second==v[i-1].second+1) rl++;
        else { run->Fill(rl); rl=1; } } }
    printf("%s: <runlen>=%.2f run<=3 frac=%.3f run==5 frac=%.3f | pix adc mean=%.1f\n",
      in, run->GetMean(), run->Integral(1,3)/std::max(1.,run->Integral()),
      run->GetBinContent(5)/std::max(1.,run->Integral()), adc->GetMean());
  };
  TH1D *rr,*ra,*dr,*da;
  prof("/home/rog/sPHENIX/3D_ClusterFindingML/clusters_seeds_island_79507-0.root_ntuplizer.root",false,rr,ra);
  prof("digi_sim.root",true,dr,da);
  TCanvas c("c","",1300,500); c.Divide(2,1);
  auto d2=[&](TH1*x,TH1*y,const char*xt,const char*ti,bool logy){
    if(x->Integral()>0)x->Scale(1./x->Integral()); if(y->Integral()>0)y->Scale(1./y->Integral());
    x->SetLineColor(kBlue+1);x->SetLineWidth(2); y->SetLineColor(kMagenta+1);y->SetLineWidth(2);
    x->SetStats(0); double mx=std::max(x->GetMaximum(),y->GetMaximum());
    x->SetMaximum(logy?mx*2:mx*1.3); if(logy){gPad->SetLogy();x->SetMinimum(1e-6);} else x->SetMinimum(0);
    x->SetTitle(ti); x->GetXaxis()->SetTitle(xt); x->Draw("HIST"); y->Draw("HIST SAME");
    TLegend*L=new TLegend(0.4,0.75,0.89,0.89);L->SetBorderSize(0);L->SetFillStyle(0);
    L->AddEntry(x,"REAL ntp_hit","l"); L->AddEntry(y,"DETACHED DIGITIZER v1","l"); L->Draw(); };
  c.cd(1); d2(rr,dr,"consecutive tbins per pad","run length: THE acceptance test",true);
  c.cd(2); d2(ra,da,"per-pixel ADC","per-pixel ADC (gain uncalibrated)",true);
  c.SaveAs("/home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/digi_day1.png");
  printf("saved digi_day1.png\n");
}
