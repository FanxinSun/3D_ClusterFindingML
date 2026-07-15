void v40_buildlibs()
{
  std::map<int,std::vector<std::pair<int,int>>> acc;
  { std::ifstream fi("v40_accept_list.txt"); int f,e,j,fl;
    while(fi>>f>>e>>j>>fl) acc[f].push_back({e,j}); }
  for(int L=0;L<10;++L){
    TFile*fr=TFile::Open(Form("raw_v40_%d.root",L));
    TNtuple*r=(TNtuple*)fr->Get("raw_pix");
    TFile*fe=TFile::Open(Form("/home/rog/sPHENIX/3D_ClusterFindingML/P5/ANG_g4hit_%d.root",L));
    TTree*g=(TTree*)fe->Get("ntp_g4hit");
    std::map<int,int> rm; for(auto&pe:acc[L]) rm[pe.first]=pe.second;
    TFile*orw=new TFile(Form("raw_lib_v40_%d.root",L),"RECREATE");
    TNtuple*nr=new TNtuple("raw_pix","v40 thin","event:layer:side:pad:tbin:q:trk");
    float rv[7]; TObjArray*br=r->GetListOfBranches();
    for(int b=0;b<7;++b) r->SetBranchAddress(br->At(b)->GetName(),&rv[b]);
    for(Long64_t i=0;i<r->GetEntries();++i){ r->GetEntry(i);
      auto it=rm.find((int)rv[0]); if(it==rm.end()) continue;
      float o[7]={(float)it->second,rv[1],rv[2],rv[3],rv[4],rv[5],rv[6]}; nr->Fill(o);}
    nr->Write(); orw->Close();
    TFile*oev=new TFile(Form("eval_v40_%d.root",L),"RECREATE");
    TTree*cg=g->CloneTree(0); cg->SetName("ntp_g4hit");
    float gev; g->SetBranchAddress("event",&gev);
    for(Long64_t i=0;i<g->GetEntries();++i){ g->GetEntry(i);
      auto it=rm.find((int)gev); if(it==rm.end()) continue;
      gev=(float)it->second; cg->Fill(); }
    cg->Write(); oev->Close();
    printf("V40LIB %d built\n", L);
  }
}
