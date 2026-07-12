// v34_vs_v35.C — DEDICATED version-decision figure (new file per the new-discovery
// rule; frozen-purpose figures are not distorted with extra curves).
// Four island-level panels, three curves each: REAL vs v3.4 (B4.3, the rewind
// baseline) vs v3.5revC (composition pass + near-threshold fixes):
//   (1) island pixel-count composition 1-30 (log)   — the singles/minis story
//   (2) cluster adc 0-300 (log)                     — the "drips" region
//   (3) cluster maxadc 900-1000                     — saturation bump vs delta
//   (4) islands/frame per layer                     — content + region shape
void v34_vs_v35()
{
  gROOT->SetBatch(1); gStyle->SetOptStat(0);
  TFile*fr=TFile::Open("island_real.root");      TTree*rr=(TTree*)fr->Get("island");
  TFile*f4=TFile::Open("island_frames_v34.root"); TTree*s4=(TTree*)f4->Get("island");
  TFile*f5=TFile::Open("island_frames_v35.root"); TTree*s5=(TTree*)f5->Get("island");
  const double NR=100., N4=250., N5=250.;
  auto mk=[&](TTree*t,const char*v,const char*hn,int nb,double lo,double hi,const char*cut){
    TH1D*h=new TH1D(hn,"",nb,lo,hi); t->Draw(Form("%s>>%s",v,hn),cut,"goff"); return h; };
  auto d3=[&](TVirtualPad*p,TH1D*a,TH1D*b,TH1D*c,const char*ti,bool logy,bool norm){
    p->cd();
    if(norm) for(TH1D*x:{a,b,c}) if(x->Integral()>0) x->Scale(1./x->Integral());
    a->SetLineColor(kBlue+1); b->SetLineColor(kGreen+2); c->SetLineColor(kMagenta+1);
    a->SetLineWidth(2); b->SetLineWidth(2); c->SetLineWidth(2);
    for(TH1D*x:{a,b,c}) x->SetStats(0);
    a->SetTitle(ti);
    double mx=std::max({a->GetMaximum(),b->GetMaximum(),c->GetMaximum()});
    a->SetMaximum(logy?mx*2.5:mx*1.3);
    if(logy){gPad->SetLogy(); a->SetMinimum(norm?1e-6:0.5);} else a->SetMinimum(0);
    a->Draw("HIST"); b->Draw("HIST SAME"); c->Draw("HIST SAME");
    TLegend*L=new TLegend(0.40,0.70,0.89,0.89); L->SetBorderSize(0); L->SetFillStyle(0);
    L->AddEntry(a,"REAL","l");
    L->AddEntry(b,"SIM v3.4 (B4.3)  [rewind baseline]","l");
    L->AddEntry(c,"SIM v3.5revC","l"); L->Draw(); };
  TCanvas c("c","",1500,1000); c.Divide(2,2);
  d3(c.cd(1), mk(rr,"size","a1",30,0.5,30.5,""), mk(s4,"size","b1",30,0.5,30.5,""),
     mk(s5,"size","c1",30,0.5,30.5,""), "island pixel count (composition);pixels;norm", true, true);
  d3(c.cd(2), mk(rr,"adc","a2",75,0,300,""), mk(s4,"adc","b2",75,0,300,""),
     mk(s5,"adc","c2",75,0,300,""), "cluster adc (low end);adc;norm", true, true);
  d3(c.cd(3), mk(rr,"maxadc","a3",50,900.5,1000.5,""), mk(s4,"maxadc","b3",50,900.5,1000.5,""),
     mk(s5,"maxadc","c3",50,900.5,1000.5,""), "cluster maxadc (saturation);ADC;norm", false, true);
  TH1D*a4=mk(rr,"layer","a4",48,6.5,54.5,""); a4->Scale(1./NR);
  TH1D*b4=mk(s4,"layer","b4",48,6.5,54.5,""); b4->Scale(1./N4);
  TH1D*c4=mk(s5,"layer","c4",48,6.5,54.5,""); c4->Scale(1./N5);
  d3(c.cd(4), a4, b4, c4, "islands per frame per layer;TPC layer;islands/frame", false, false);
  c.SaveAs("/home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/v34_vs_v35revC.png");
  printf("saved v34_vs_v35revC.png\n");
}
