// data-driven TPC pad geometry table from existing sim hits:
// per (layer, side): linear fit phi = phi0 + slope*(phibin+0.5); plus mean radius per layer
void geomfit(){
  TFile*f=TFile::Open("/home/rog/sPHENIX/3D_ClusterFindingML/macros-offline/detectors/sPHENIX/exam5_g4svtx_eval.root");
  TTree*h=(TTree*)f->Get("ntp_hit");
  TFile*fc=TFile::Open("/home/rog/sPHENIX/3D_ClusterFindingML/macros-offline/detectors/sPHENIX/exam5_trkrntuple.root");
  TTree*c=(TTree*)fc->Get("ntp_cluster");
  FILE*out=fopen("tpc_geom_table.txt","w");
  fprintf(out,"# layer nbins radius slope phi0_side0 phi0_side1\n");
  h->SetEstimate(3000000);
  for(int L=7;L<55;++L){
    // radius from clusters
    TH1D hr("hr","",4000,20,90); c->Draw("r>>hr",Form("layer==%d",L),"goff");
    double rad=hr.GetMean();
    double slope=0, phi0[2]={0,0};
    for(int sd=0;sd<2;++sd){
      long n=h->Draw("phi:phibin",Form("layer==%d&&zelem==%d&&event<3",L,sd),"goff");
      if(n<50){ phi0[sd]=-9999; continue; }
      double*py=h->GetVal(0), *px=h->GetVal(1);
      // robust linear fit: slope known ~2pi/nbins; fit intercept via median of phi - slope*(bin+0.5)
      // first estimate slope by least squares
      double sx=0,sy=0,sxx=0,sxy=0;
      for(long i=0;i<n;++i){ sx+=px[i]; sy+=py[i]; sxx+=px[i]*px[i]; sxy+=px[i]*py[i]; }
      double sl=(n*sxy-sx*sy)/(n*sxx-sx*sx);
      double in=(sy-sl*sx)/n;
      if(sd==0||slope==0) slope=sl;
      phi0[sd]=in;
    }
    int nb=(int)std::lround(2*M_PI/std::fabs(slope));
    fprintf(out,"%d %d %.4f %.10f %.6f %.6f\n",L,nb,rad,slope,phi0[0],phi0[1]);
  }
  fclose(out);
  printf("wrote tpc_geom_table.txt\n");
  // show a sample
  gSystem->Exec("head -4 tpc_geom_table.txt; grep '^23 ' tpc_geom_table.txt; grep '^40 ' tpc_geom_table.txt");
}
