// pilot_compmeter.C — v6.1 width re-balance PILOT composition meters
// (width_rebalance_request.md ADJUDICATED AMENDMENT, 2026-08-23).
// Sim-only, chunk-0 pilot digi: truth-grouped (event,trk) raw-pixel tracks,
// fc/nf bar (>=12 rows, span>=15 cm, R_fit>=45 cm), global circle fit,
// 3x(1.4826xMAD) clip x3 (gtail discipline), then per-row mean residuals ->
//   C(d)/C(0) row-lag autocorrelation, d = 1,2,4,8 (per-track de-meaned,
//     pooled)  [probe production values: real 0.30/0.18/0.08/-0.06;
//     sim post-twist C(1) 0.57 — pilot values RANK candidates, the
//     probe battery on the production is the acceptance instrument]
//   cell coherence = RMS over (side,sector,layer) cells of the per-cell
//     mean row residual (n>=20 rows/cell)  [real 307 um, sim 145].
#include <TFile.h>
#include <TNtuple.h>
#include <cmath>
#include <cstdio>
#include <map>
#include <vector>
#include <algorithm>

namespace PCM
{
struct Fit { double a=0,b=0,R=0; int n=0; bool ok=false; };
Fit fitC(const std::vector<double>&X,const std::vector<double>&Y)
{
  Fit F; F.n=(int)X.size(); if(F.n<3) return F;
  double mx=0,my=0; for(int i=0;i<F.n;++i){mx+=X[i];my+=Y[i];} mx/=F.n; my/=F.n;
  double suu=0,suv=0,svv=0,suuu=0,svvv=0,suvv=0,svuu=0;
  for(int i=0;i<F.n;++i){double u=X[i]-mx,v=Y[i]-my;
    suu+=u*u;suv+=u*v;svv+=v*v;suuu+=u*u*u;svvv+=v*v*v;suvv+=u*v*v;svuu+=v*u*u;}
  double det=suu*svv-suv*suv; if(std::fabs(det)<1e-12) return F;
  double uc=(svv*(suuu+suvv)-suv*(svvv+svuu))/(2*det);
  double vc=(suu*(svvv+svuu)-suv*(suuu+suvv))/(2*det);
  F.a=uc+mx; F.b=vc+my; F.R=std::sqrt(uc*uc+vc*vc+(suu+svv)/F.n);
  for(int it=0;it<6;++it){
    double JtJ00=0,JtJ01=0,JtJ02=0,JtJ11=0,JtJ12=0,JtJ22=0,Jr0=0,Jr1=0,Jr2=0;
    for(int i=0;i<F.n;++i){double dx=X[i]-F.a,dy=Y[i]-F.b,ri=std::hypot(dx,dy);
      if(ri<1e-9)continue; double e=ri-F.R,j0=-dx/ri,j1=-dy/ri;
      JtJ00+=j0*j0;JtJ01+=j0*j1;JtJ02+=-j0;JtJ11+=j1*j1;JtJ12+=-j1;JtJ22+=1;
      Jr0+=j0*e;Jr1+=j1*e;Jr2+=-e;}
    double A[3][4]={{JtJ00,JtJ01,JtJ02,-Jr0},{JtJ01,JtJ11,JtJ12,-Jr1},{JtJ02,JtJ12,JtJ22,-Jr2}};
    for(int c=0;c<3;++c){int p=c; for(int r2=c+1;r2<3;++r2) if(std::fabs(A[r2][c])>std::fabs(A[p][c]))p=r2;
      if(std::fabs(A[p][c])<1e-12){F.ok=false;return F;} std::swap(A[p],A[c]);
      for(int r2=0;r2<3;++r2){ if(r2==c)continue; double f=A[r2][c]/A[c][c];
        for(int cc=c;cc<4;++cc)A[r2][cc]-=f*A[c][cc]; } }
    F.a+=A[0][3]/A[0][0]; F.b+=A[1][3]/A[1][1]; F.R+=A[2][3]/A[2][2];
  }
  F.ok=true; return F;
}
}

void pilot_compmeter(const char *digi = "v61pilot_base6t_digi.root", const char *tag = "base6t", double keepfrac = 1.0)
{
  using namespace PCM;
  double rowR[55]={0}; { FILE*g=fopen("tpc_geom_table.txt","r"); char l[256];
    while(fgets(l,256,g)){int L,nb;double r,sl,p0,p1; if(l[0]=='#')continue;
      if(sscanf(l,"%d %d %lf %lf %lf %lf",&L,&nb,&r,&sl,&p0,&p1)==6) rowR[L]=r;} fclose(g);}
  TFile*f=TFile::Open(digi); TNtuple*t=(TNtuple*)f->Get("ntp_hit");
  float ev,lay,pad,tb,adc,sd,phi,z,trk;
  t->SetBranchStatus("*",0);
  for(auto b:{"event","layer","phibin","tbin","adc","zelem","phi","z","gtrackID"}) t->SetBranchStatus(b,1);
  t->SetBranchAddress("event",&ev); t->SetBranchAddress("layer",&lay); t->SetBranchAddress("phibin",&pad);
  t->SetBranchAddress("tbin",&tb); t->SetBranchAddress("adc",&adc); t->SetBranchAddress("zelem",&sd);
  t->SetBranchAddress("phi",&phi); t->SetBranchAddress("z",&z); t->SetBranchAddress("gtrackID",&trk);
  struct Px{float x,y;int L,sct,sd;};
  std::map<std::pair<int,int>,std::vector<Px>> G;
  for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    int L=(int)lay; if(L<7||L>54||trk<=0) continue;
    if(keepfrac<1.0){  // deterministic px-density degrade (real 3.25 vs sim 6.67 px/row)
      unsigned h=2166136261u; for(int vv:{(int)ev,(int)trk,L,(int)pad,(int)tb}){h^=(unsigned)vv; h*=16777619u;}
      if((h%10000)>=(unsigned)(keepfrac*10000)) continue; }
    double r=rowR[L]; double phw=phi<0?phi+2*M_PI:phi;
    G[{(int)ev,(int)trk}].push_back({(float)(r*cos(phi)),(float)(r*sin(phi)),L,std::min(11,(int)(phw/(M_PI/6.))),(int)sd}); }
  double sC[9]={0}; long nC[9]={0};
  std::map<long,std::pair<double,long>> cell; // (side*1000+sct)*100+L -> (sum m, n)
  long ntr=0;
  for(auto&kv:G){ auto&v=kv.second;
    std::vector<int> rows; for(auto&p:v) rows.push_back(p.L);
    std::sort(rows.begin(),rows.end()); rows.erase(std::unique(rows.begin(),rows.end()),rows.end());
    if((int)rows.size()<12) continue;
    if(rowR[rows.back()]-rowR[rows.front()]<15.) continue;
    std::vector<double> X,Y; for(auto&p:v){X.push_back(p.x);Y.push_back(p.y);}
    Fit F=fitC(X,Y); if(!F.ok||F.R<45||F.R>2e4) continue;
    std::vector<double> res(v.size());
    for(size_t i=0;i<v.size();++i) res[i]=std::hypot(v[i].x-F.a,v[i].y-F.b)-F.R;
    std::vector<char> keep(v.size(),1);
    for(int it=0;it<3;++it){ std::vector<double> rr2; for(size_t i=0;i<v.size();++i) if(keep[i]) rr2.push_back(res[i]);
      if(rr2.size()<5)break; std::nth_element(rr2.begin(),rr2.begin()+rr2.size()/2,rr2.end());
      double med=rr2[rr2.size()/2]; std::vector<double> ad; for(double q:rr2) ad.push_back(std::fabs(q-med));
      std::nth_element(ad.begin(),ad.begin()+ad.size()/2,ad.end());
      double clip=std::max(0.05, 3*1.4826*ad[ad.size()/2]);
      for(size_t i=0;i<v.size();++i) keep[i]=std::fabs(res[i]-med)<=clip; }
    std::map<int,std::pair<double,int>> rm; std::map<int,std::pair<int,int>> rmc;
    for(size_t i=0;i<v.size();++i){ if(!keep[i])continue;
      rm[v[i].L].first+=res[i]; rm[v[i].L].second++;
      rmc[v[i].L]={v[i].sct,v[i].sd}; }
    if((int)rm.size()<8) continue;
    std::map<int,double> M; double mb=0; int nm=0;
    for(auto&r2:rm){ if(r2.second.second<2)continue; M[r2.first]=r2.second.first/r2.second.second; mb+=M[r2.first]; nm++; }
    if(nm<8) continue; mb/=nm; ntr++;
    for(auto&m:M){ double a=m.second-mb;
      long ck=((long)(rmc[m.first].second*1000+rmc[m.first].first))*100+m.first;
      cell[ck].first+=m.second; cell[ck].second++;
      for(int d:{0,1,2,4,8}){ auto it2=M.find(m.first+d); if(it2==M.end())continue;
        sC[d]+=a*(it2->second-mb); nC[d]++; } }
  }
  double C0=nC[0]>0?sC[0]/nC[0]:1;
  double ss=0; long ncell=0;
  for(auto&c:cell){ if(c.second.second<20)continue; double m=c.second.first/c.second.second; ss+=m*m; ncell++; }
  printf("PCOMP %-14s tracks %ld | C(d)/C(0): d1 %.3f d2 %.3f d4 %.3f d8 %.3f | cellcoh %.0f um (%ld cells) | rowRMS %.0f um\n",
         tag,ntr,sC[1]/nC[1]/C0,sC[2]/nC[2]/C0,sC[4]/nC[4]/C0,sC[8]/nC[8]/C0,
         std::sqrt(ss/std::max(1L,ncell))*1e4,ncell,std::sqrt(C0)*1e4);
  f->Close();
}
