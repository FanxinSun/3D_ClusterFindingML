#!/usr/bin/env bash
set -o pipefail
cd /home/rog/sPHENIX/3D_ClusterFindingML/island_post
DM=/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/TPC_DEADCHANNELMAP/ff/c3/ffc3f6498934c5a8ba31065292c6ebcc_TPCDeadMap_79471.root
root -l -b -q -e 'gROOT->ProcessLine(".L tpc_digitize.C+");
tpc_transport("/home/rog/sPHENIX/3D_ClusterFindingML/P5/ANG30_g4hit.root","raw_ang30.root",30,0.05,0.00953);' 2>&1 | grep tpc_transport:
root -l -b -q -e "gROOT->ProcessLine(\".L tpc_digitize.C+\");
tpc_readout(\"raw_ang30.root\",\"digi_ang30.root\",0.86,20.0,1,1,4711,\"$DM\",11.0,0.39,0.55,0.80,0.021,7.0,36.0,70.0,11.0,940.0,2,0.29,10.0,1.24,1.06,-1.0,5.0,0.0);" 2>&1 | grep "tpc_readout: raw_ang30"
root -l -b -q -e '
TFile*f=TFile::Open("digi_ang30.root"); TTree*t=(TTree*)f->Get("ntp_hit");
double s=0,s2=0; int n=30;
for(int e=0;e<n;++e){ double k=t->GetEntries(Form("event==%d",e)); s+=k; s2+=k*k; }
double m=s/n, r=sqrt(s2/n-m*m);
printf("ANGCENSUS 30 collisions | kept mean %.0f | rms/mean %.2f\n", m, r/m);
printf("VERDICT-TABLE: HIJING 17936 | ANGANTYR %.0f | REAL-uniform 7212 | REAL-fired 8200\n", m);' 2>&1 | grep -E "ANGCENSUS|VERDICT"
echo ANG-CENSUS-DONE
