// one-pass real-data reference histograms for day-2 calibration (built ONCE)
#include <map>
#include <vector>
#include <algorithm>
void make_ref(){
  TFile*f=TFile::Open("/home/rog/sPHENIX/3D_ClusterFindingML/clusters_seeds_island_79507-0.root_ntuplizer.root");
  TTree*h=(TTree*)f->Get("ntp_hit");
  float event,layer,phibin,tb,adc,side,phi;
  h->SetBranchStatus("*",0);
  for(const char* b : {"event","layer","phibin","tbin","adc","zelem","phi"}) h->SetBranchStatus(b,1);
  h->SetBranchAddress("event",&event); h->SetBranchAddress("layer",&layer);
  h->SetBranchAddress("phibin",&phibin); h->SetBranchAddress("tbin",&tb);
  h->SetBranchAddress("adc",&adc); h->SetBranchAddress("zelem",&side); h->SetBranchAddress("phi",&phi);
  TFile*o=new TFile("ref_real.root","RECREATE");
  TH1D* r_adc  = new TH1D("r_adc","real per-pixel ADC;ADC",1100,0,1100);
  TH1D* r_adcz = new TH1D("r_adcz","real ADC near threshold;ADC",240,0,120);
  TH1D* r_run  = new TH1D("r_run","real run length;consecutive tbins",25,0.5,25.5);
  TH1D* r_lay  = new TH1D("r_lay","real hits/event/layer;layer",48,6.5,54.5);
  TH1D* r_fold = new TH1D("r_fold","real phi fold;phi mod 30deg",120,0,0.5236);
  TH1D* r_tbin = new TH1D("r_tbin","real tbin;tbin",1000,0,1000);
  std::map<uint64_t,std::vector<std::pair<int,int>>> cols;
  Long64_t N=h->GetEntries();
  for(Long64_t i=0;i<N;++i){ h->GetEntry(i);
    if(layer<7||layer>54||adc<=0) continue;
    r_adc->Fill(adc); r_adcz->Fill(adc); r_lay->Fill(layer); r_tbin->Fill(tb);
    double ph=phi; while(ph<0)ph+=2*M_PI; r_fold->Fill(fmod(ph,0.5235988));
    uint64_t k=((uint64_t)(uint32_t)event<<24U)|((uint64_t)(uint32_t)layer<<8U)|(uint64_t)(((int)side==1)?1:0);
    cols[k].push_back({(int)phibin,(int)tb}); }
  for(auto&g:cols){ auto&v=g.second; std::sort(v.begin(),v.end());
    int rl=1;
    for(size_t i=1;i<=v.size();++i){
      if(i<v.size() && v[i].first==v[i-1].first && v[i].second==v[i-1].second+1) rl++;
      else { r_run->Fill(rl); rl=1; } } }
  r_lay->Scale(1./100.);  // per event (100 events)
  o->Write();
  printf("ref_real.root: adc mean %.2f | run mean %.2f | anchors: pixmean 99.1, maxadc-sat ~949, runle3 0.816\n",
    r_adc->GetMean(), r_run->GetMean());
  o->Close();
}
