// cmp_pau.C — exam6 CONTAINER-SIDE comparisons (successor of the exam2/5-era
// cmp_hit / cmp_cluster / cmp_auau_* figures): the in-container chain (SvtxEvaluator
// digitized hits + ported TrkrNtuplizer 91-branch clusters) run on the pAu single-collision
// library, vs real. Shape-level (single collisions vs streaming frames — occupancies differ
// by construction; all curves area-normalized).
void cmp_pau(){
  gROOT->SetBatch(1); gStyle->SetOptStat(0);
  const char* M="/home/rog/sPHENIX/3D_ClusterFindingML/macros-offline/detectors/sPHENIX";
  TFile*fr=TFile::Open("/home/rog/sPHENIX/3D_ClusterFindingML/clusters_seeds_island_79507-0.root_ntuplizer.root");
  TTree*rh=(TTree*)fr->Get("ntp_hit"); TTree*rc=(TTree*)fr->Get("ntp_cluster");
  TChain*sh=new TChain("ntp_hit");     // container-digitized hits (evaluator)
  TChain*sc=new TChain("ntp_cluster"); // ported 91-branch on container clusters
  for(const char* t : {"pau_a_eval","pauL_a_eval","pauL_b_eval","pauL_c_eval","pauL_d_eval"}){
    sh->Add(Form("%s/%s_g4svtx_eval.root",M,t));
    sc->Add(Form("%s/%s_trkrntuple.root",M,t));
  }
  auto d2=[&](TVirtualPad*p,TH1D*a,TH1D*b,const char*ti,bool logy,const char*lb){
    p->cd(); if(a->Integral()>0)a->Scale(1./a->Integral()); if(b->Integral()>0)b->Scale(1./b->Integral());
    a->SetLineColor(kBlue+1);a->SetLineWidth(2); b->SetLineColor(kOrange+7);b->SetLineWidth(2);
    a->SetStats(0);b->SetStats(0); a->SetTitle(ti);
    double mx=std::max(a->GetMaximum(),b->GetMaximum()); a->SetMaximum(logy?mx*2.5:mx*1.3);
    if(logy){gPad->SetLogy();a->SetMinimum(1e-7);} else a->SetMinimum(0);
    a->Draw("HIST"); b->Draw("HIST SAME");
    TLegend*L=new TLegend(0.38,0.76,0.89,0.89); L->SetBorderSize(0); L->SetFillStyle(0);
    L->AddEntry(a,"REAL","l"); L->AddEntry(b,lb,"l"); L->Draw(); };
  auto H=[&](TTree*t,const char*v,const char*hn,int nb,double lo,double hi,const char*cut){
    TH1D*h=new TH1D(hn,"",nb,lo,hi); t->Draw(Form("%s>>%s",v,hn),cut,"goff"); return h; };
  const char* RC="layer>=7&&layer<=54&&adc>0";
  const char* SC="layer>=7&&layer<=54&&adc>0";
  // ---- hits: THREE-WAY — real vs detached frames (CURRENT production v3.3) vs container ----
  TFile*fb=TFile::Open("digi_frames_production_v36.root"); TTree*bh=(TTree*)fb->Get("ntp_hit");
  auto d3=[&](TVirtualPad*p,TH1D*a,TH1D*b,TH1D*g,const char*ti,bool logy,const char*slab="SIM detached pipeline (v3.6)"){
    p->cd(); for(TH1D*x:{a,b,g}) if(x->Integral()>0)x->Scale(1./x->Integral());
    a->SetLineColor(kBlue+1); b->SetLineColor(kMagenta+1);
    g->SetLineColor(kGray+2); g->SetLineStyle(2); g->SetLineWidth(1);
    a->SetLineWidth(2); b->SetLineWidth(2);
    for(TH1D*x:{a,b,g}) x->SetStats(0);
    a->SetTitle(ti);
    double mx=std::max({a->GetMaximum(),b->GetMaximum(),g->GetMaximum()});
    a->SetMaximum(logy?mx*2.5:mx*1.3);
    if(logy){gPad->SetLogy();a->SetMinimum(1e-7);} else a->SetMinimum(0);
    a->Draw("HIST"); b->Draw("HIST SAME"); g->Draw("HIST SAME");
    TLegend*L=new TLegend(0.35,0.70,0.89,0.89); L->SetBorderSize(0); L->SetFillStyle(0);
    L->AddEntry(a,"REAL","l");
    L->AddEntry(b,slab,"l");
    L->AddEntry(g,"stock ana.331 (reference)","l"); L->Draw(); };
  TCanvas c1("c1","",1500,500); c1.Divide(3,1);
  d3(c1.cd(1),H(rh,"adc","a1",220,-0.5,1099.5,RC),H(bh,"adc","m1",220,-0.5,1099.5,""),
     H(sh,"adc-74.4","b1",220,-0.5,1099.5,SC),"per-hit ADC;ADC;norm",true);
  d3(c1.cd(2),H(rh,"adc","a2",121,-0.5,120.5,RC),H(bh,"adc","m2",121,-0.5,120.5,""),
     H(sh,"adc-74.4","b2",121,-0.5,120.5,SC),"near threshold;ADC;norm",true);
  d3(c1.cd(3),H(rh,"layer","a3",48,6.5,54.5,RC),H(bh,"layer","m3",48,6.5,54.5,""),
     H(sh,"layer","b3",48,6.5,54.5,SC),"layer profile (shape);TPC layer;norm",false);
  c1.SaveAs("/home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/cmp_pau_hit.png");
  // ---- clusters: three-way (REAL production | detached islandize91 v2 | stock ghost) ----
  // P4: sim clusters now come from the PORTED production clusterizer on sim digi
  // (prodclus_v36.root) — same algorithm as the real ntp_cluster curve, so the
  // cluster canvas is apples-to-apples (the old island91 source mixed algorithms;
  // named residual (a) resolved here).
  TFile*fb2=TFile::Open("prodclus_v36.root"); TTree*bc=(TTree*)fb2->Get("ntp_clus");
  TCanvas c2("c2","",1500,900); c2.Divide(3,2);
  const char* RCC="layer>=7&&layer<=54";
  const char* PL="SIM v3.6 #times ported prod. clusterizer";
  d3(c2.cd(1),H(rc,"adc","d1",150,0,6000,RCC),H(bc,"adc","f1",150,0,6000,""),H(sc,"adc","e1",150,0,6000,""),
     "cluster ADC;adc;norm",true,PL);
  d3(c2.cd(2),H(rc,"size","d2",60,0.5,60.5,RCC),H(bc,"phisize*zsize","f2",60,0.5,60.5,""),H(sc,"phisize*zsize","e2",60,0.5,60.5,""),
     "cluster bbox AREA (= real size semantics: 99.5% size==phisize#timeszsize);phisize #times zsize;norm",true,PL);
  d3(c2.cd(3),H(rc,"phisize","d3h",20,0.5,20.5,RCC),H(bc,"phisize","f3",20,0.5,20.5,""),H(sc,"phisize","e3",20,0.5,20.5,""),
     "cluster #phi-size  [ALL production-capped now];pads;norm",true,PL);
  d3(c2.cd(4),H(rc,"zsize","d4",20,0.5,20.5,RCC),H(bc,"zsize","f4",20,0.5,20.5,""),H(sc,"zsize","e4",20,0.5,20.5,""),
     "cluster z-size  [ALL production-capped now];tbins;norm",true,PL);
  d3(c2.cd(5),H(rc,"maxadc","d5",120,0,1080,RCC),H(bc,"maxadc","f5",120,0,1080,""),H(sc,"maxadc-74.4","e5",120,0,1080,""),
     "cluster maxadc  [stock: -74.4];ADC;norm",true,PL);
  d3(c2.cd(6),H(rc,"e","d6",150,0,6000,RCC),H(bc,"adc","f6",150,0,6000,""),H(sc,"e","e6",150,0,6000,""),
     "cluster e  [port: e#equiv adc, verified];e;norm",true,PL);
  c2.SaveAs("/home/rog/sPHENIX/3D_ClusterFindingML/sim_validation_plots/cmp_pau_cluster.png");
  printf("saved cmp_pau_hit.png + cmp_pau_cluster.png (container-side exam6)\n");
}
