// p3_cmp.C — production comparison, ANCHOR-CONSISTENT (follows the CURRENT
// production since the 2026-07-13 full-suite replot; was the v1/240kHz record —
// that version lives in git history at 81cd06d):
// real pixel curves from ref_real.root (layer 7-54 && adc>0); islands from island91_real.
void p3_cmp(const char *simfile = "digi_frames_production_v40b.root",
            const char *suffix = "", const char *tag = "SIM v4.0",
            const char *simisl = "island91_frames_production_v40b.root"){
  gROOT->SetBatch(1); gStyle->SetOptStat(0);
  TFile*fref=TFile::Open("ref_real.root");
  TH1D*RZ=(TH1D*)((TH1D*)fref->Get("r_adcz"))->Clone("RZ"); RZ->SetDirectory(nullptr);
  TH1D*RT=(TH1D*)((TH1D*)fref->Get("r_tbin"))->Clone("RT"); RT->SetDirectory(nullptr);
  TFile*fd=TFile::Open(simfile); TTree*h=(TTree*)fd->Get("ntp_hit");
  TFile*fri=TFile::Open("island91_real.root"); TTree*ri=(TTree*)fri->Get("ntp_cluster");
  TFile*fsi=TFile::Open(simisl); TTree*si=(TTree*)fsi->Get("ntp_cluster");
  auto d2=[&](TVirtualPad*p,TH1D*a,TH1D*b,const char*ti,bool logy){
    p->cd(); if(a->Integral()>0)a->Scale(1./a->Integral()); if(b->Integral()>0)b->Scale(1./b->Integral());
    a->SetLineColor(kBlue+1);a->SetLineWidth(2); b->SetLineColor(kMagenta+1);b->SetLineWidth(2);
    a->SetStats(0);b->SetStats(0); a->SetTitle(ti);
    double mx=std::max(a->GetMaximum(),b->GetMaximum()); a->SetMaximum(logy?mx*2.5:mx*1.3);
    if(logy){gPad->SetLogy();a->SetMinimum(1e-7);} else a->SetMinimum(0);
    a->Draw("HIST"); b->Draw("HIST SAME");
    TLegend*L=new TLegend(0.45,0.75,0.89,0.89); L->SetBorderSize(0); L->SetFillStyle(0);
    L->AddEntry(a,"REAL","l"); L->AddEntry(b,tag,"l"); L->Draw(); };
  auto H=[&](TTree*x,const char*v,const char*hn,int nb,double lo,double hi){
    TH1D*hh=new TH1D(hn,"",nb,lo,hi); x->Draw(Form("%s>>%s",v,hn),"","goff"); return hh; };
  TCanvas c("c","",1600,900); c.Divide(3,2);
  d2(c.cd(1),RZ,H(h,"adc","b1",121,-0.5,120.5),"per-pixel ADC (ZS region);ADC;norm",true);
  d2(c.cd(2),RT,H(h,"tbin","b2",1000,0,1000),"arrivals (adc>0);tbin;norm",false);
  d2(c.cd(3),H(ri,"size","a3",100,0.5,100.5),H(si,"size","b3",100,0.5,100.5),"island size;pixels;norm",true);
  d2(c.cd(4),H(ri,"phisize","a4",25,0.5,25.5),H(si,"phisize","b4",25,0.5,25.5),"island #phi-size;pads;norm",true);
  d2(c.cd(5),H(ri,"zsize","a5",40,0.5,40.5),H(si,"zsize","b5",40,0.5,40.5),"island z-size;tbins;norm",true);
  d2(c.cd(6),H(ri,"adc","a6",150,0,6000),H(si,"adc","b6",150,0,6000),"island ADC;raw sum;norm",true);
  c.SaveAs(Form("/home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/p3_production_cmp%s.png",suffix));
  printf("p3_production_cmp%s.png regenerated anchor-consistent\n",suffix);
}
