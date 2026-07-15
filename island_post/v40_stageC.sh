#!/usr/bin/env bash
# v4.0 stage C: full-set census -> thin refit -> libs -> production -> acceptance
set -o pipefail
cd /home/rog/sPHENIX/3D_ClusterFindingML/island_post
DM=/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/TPC_DEADCHANNELMAP/ff/c3/ffc3f6498934c5a8ba31065292c6ebcc_TPCDeadMap_79471.root
NEV=(125 125 125 125 250 250 250 250 250 250)
echo "--- C1: transport + census readout, 10 files ---"
for i in 0 1 2 3 4 5 6 7 8 9; do
root -l -b -q -e "gROOT->ProcessLine(\".L tpc_digitize.C+\");
tpc_transport(\"/home/rog/sPHENIX/3D_ClusterFindingML/P5/ANG_g4hit_$i.root\",\"raw_v40_$i.root\",${NEV[$i]},0.040,0.0070);
tpc_readout(\"raw_v40_$i.root\",\"cens_v40_$i.root\",1.03,20.0,1,1,4711,\"$DM\",11.0,0.39,0.55,0.70,0.021,7.0,36.0,70.0,11.0,940.0,2,0.29,10.0,1.24,1.06,-1.0,5.0,0.0);" 2>&1 | grep -cE "transport|readout" > /dev/null
echo "C1 file $i done"
done
root -l -b -q -e '
FILE*fo=fopen("v40_census.txt","w");
int NEV[10]={125,125,125,125,250,250,250,250,250,250};
for(int L=0;L<10;++L){
  TFile*f=TFile::Open(Form("cens_v40_%d.root",L));
  TTree*t=(TTree*)f->Get("ntp_hit");
  for(int e=0;e<NEV[L];++e) fprintf(fo,"%d %d %.0f\n",L,e,(double)t->GetEntries(Form("event==%d",e)));
  f->Close();}
fclose(fo); printf("census done\n");' 2>&1 | grep census
echo "C1-DONE"
python3 v40_thin_fit.py
echo "C2-DONE (thin fit + accept list + mbd file)"
root -l -b -q v40_buildlibs.C 2>&1 | grep -E "V40LIB|Error" | head -12
echo "C3-DONE (thinned libs)"
./v40_production.sh "5.80:1.3542,11.10:1.2654,16.40:1.1901,21.70:0.9930,27.00:1.1215,32.30:1.0999,37.60:0.8274,42.90:0.7611,47.96:0.8312"
echo "C4-DONE (production)"
./v40_accept.sh 2>&1 | grep -E "WINDOWS|BUMP|STEP|B42METRIC|B43REGION|ISLANDS"
echo "STAGE-C-COMPLETE"
