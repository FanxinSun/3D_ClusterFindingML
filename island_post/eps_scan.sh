#!/bin/bash
# eps_scan.sh — eps_MBD scan (user protocol 2026-07-23): 6 points, per-point
# envelope re-derivation (env_solve.py), CRN seeds (flat 20260731, pilot
# 20260732), frozen v5.1 electronics, battery vs COMPLETE-62 targets.
# FIT FAMILY = Class B (w1/w3, step, bump, CV). Class A (content family)
# reported but NEVER fit (laundering trap). Class C = untouched holdouts.
# Idempotent: point skipped if eps_done_<eps> marker exists.
set -eo pipefail
cd /home/rog/sPHENIX/3D_ClusterFindingML/island_post
DM=/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/TPC_DEADCHANNELMAP/ff/c3/ffc3f6498934c5a8ba31065292c6ebcc_TPCDeadMap_79471.root
LIBS=""; EVALS=""
for i in 0 1 2 3 4 5 6 7 8 9; do LIBS+="raw_lib_pp_$i.root,"; EVALS+="eval_pp_$i.root,"; done
LIBS=${LIBS%,}; EVALS=${EVALS%,}
MBD="pp_mbd.txt|pp_run_0.dat,pp_run_1.dat,pp_run_2.dat,pp_run_3.dat,pp_run_4.dat,pp_run_5.dat,pp_run_6.dat,pp_run_7.dat,pp_run_8.dat,pp_run_9.dat"
SPEC="0.008:1,0.009:1,0.011:1,0.012:1,0.013:1,0.014:1,0.018:1,0.021:1,0.027:1,0.037:1,2.2:0.25"
ENVFLAT="5.80:1.0,11.10:1.0,16.40:1.0,21.70:1.0,27.00:1.0,32.30:1.0,37.60:1.0,42.90:1.0,47.96:1.0"
REALC="345.6 371.7 339.5 322.4 325.9 334.1 324.9 328.0 323.7 304.3"
TRIGC="27.96 31.85 18.74 0.20 0.03 0 0 0 0 0"
RO='1.005,20.0,1,1,4711,"'$DM'",11.0,0.39,0.55,0.81,0.016,7.0,36.0,70.0,11.0,940.0,2,0.29,10.0,1.24,1.06,-1.0,5.0'

coarse () {  # file -> 10 per-tbin/frame values on stdout
  root -l -b -q -e '
  TFile*f=TFile::Open("'$1'"); TTree*t=(TTree*)f->Get("ntp_hit");
  TH1D h("h","",971,-0.5,970.5); t->Draw("tbin>>h","","goff");
  for(int k=0;k<10;++k){ double s=0; int nb=0;
    for(int b=k*97;b<(k+1)*97&&b<=970;++b){ if(b>=318&&b<=345) continue; s+=h.GetBinContent(b+1); nb++; }
    printf("%.2f ", s/nb/'$2'.); } printf("\n");' 2>/dev/null | tail -1 || true
}

for EPS in 0.45 0.50 0.519 0.55 0.60 0.65; do
  [ -f eps_done_$EPS ] && { echo "EPSSCAN $EPS: marker exists, skipping"; continue; }
  echo "=== EPSSCAN point eps=$EPS ==="
  RS=$(python3 -c "
vals=[float(v.split(':')[0]) for v in open('rspec99_v51.txt').readlines()[-1].split(',')]
print(','.join(f'{v*0.519/$EPS:.0f}:1' for v in vals))")
  root -l -b -q -e "
  gROOT->ProcessLine(\".L frame_composer.C+\");
  frame_composer(\"$LIBS\",\"fFLx.root\",100,600.,20260731,0,\"$EVALS\",\"raw_lib_cmflash_w.root\",0.44,1.0,\"$SPEC\",1.5,2.0,1.1,0.75,\"$RS\",0.27,\"$MBD\",1.0,0.0,\"$ENVFLAT\");
  gROOT->ProcessLine(\".L tpc_digitize.C+\");
  tpc_readout(\"fFLx.root\",\"dFLx.root\",$RO,1.0);" >> pp_logs/eps_scan_steps.log 2>&1
  FLATC=$(coarse dFLx.root 100)
  ENVX=$(python3 env_solve.py "$REALC" "$TRIGC" "$FLATC")
  echo "EPSSCAN $EPS flatcoarse: $FLATC"
  echo "EPSSCAN $EPS env: $ENVX"
  root -l -b -q -e "
  gROOT->ProcessLine(\".L frame_composer.C+\");
  frame_composer(\"$LIBS\",\"fEPx.root\",100,600.,20260732,0,\"$EVALS\",\"raw_lib_cmflash_w.root\",0.44,1.0,\"$SPEC\",1.5,2.0,1.1,0.75,\"$RS\",0.27,\"$MBD\",1.0,0.0,\"$ENVX\");
  gROOT->ProcessLine(\".L tpc_digitize.C+\");
  tpc_readout(\"fEPx.root\",\"dEPx.root\",$RO,1.0);" >> pp_logs/eps_scan_steps.log 2>&1
  root -l -b -q -e '
  TFile*f=TFile::Open("dEPx.root"); TTree*t=(TTree*)f->Get("ntp_hit");
  TH1D h("h","",971,-0.5,970.5); t->Draw("tbin>>h","","goff");
  double tot=0,w1=0,w2=0,w3=0;
  for(int b=60;b<=950;++b){ if(b>=318&&b<=345) continue; double v=h.GetBinContent(b+1); tot+=v;
    if(b<=240) w1+=v; else if(b>=270&&b<=600) w2+=v; else if(b>=650) w3+=v; }
  double pre=0,post=0,a=0,b2=0;
  for(int x=230;x<=246;++x) pre+=h.GetBinContent(x+1);
  for(int x=254;x<=270;++x) post+=h.GetBinContent(x+1);
  for(int x=60;x<=240;++x) a+=h.GetBinContent(x+1);
  for(int x=270;x<=450;++x) b2+=h.GetBinContent(x+1);
  float ev,lay,adc; t->SetBranchStatus("*",0);
  t->SetBranchStatus("event",1); t->SetBranchAddress("event",&ev);
  t->SetBranchStatus("adc",1); t->SetBranchAddress("adc",&adc);
  std::map<int,long> n; long la=0;
  for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i); n[(int)ev]++; if(adc<30) la++; }
  double s=0,s2=0; for(auto&kv:n){ s+=kv.second; s2+=1.0*kv.second*kv.second; }
  double m=s/n.size(), sd=sqrt(s2/n.size()-m*m);
  printf("EPSSCAN_B: WINDOWS %.4f %.4f %.4f | BUMP %.4f | STEP %.4f | CV %.3f\n",
    w1/tot,w2/tot,w3/tot,(a/181.)/(b2/181.),pre/post,sd/m);
  printf("EPSSCAN_A: px/ev %.0f | lowadc/ev %.0f\n", m, (double)la/n.size());' 2>&1 | grep EPSSCAN || true
  root -l -b -q -e 'gROOT->ProcessLine(".L islandize.C+"); islandize("dEPx.root","iEPx.root",1);' >> pp_logs/eps_scan_steps.log 2>&1
  root -l -b -q -e '
  TFile*f=TFile::Open("iEPx.root"); TTree*t=(TTree*)f->Get("island");
  TH1D h1("h1","",200,0,200); t->Draw("size>>h1","","goff");
  printf("EPSSCAN_A2: islands/ev %.0f | <size> %.2f\n", t->GetEntries()/100., h1.GetMean());' 2>&1 | grep EPSSCAN || true
  rm -f fFLx.root dFLx.root fEPx.root dEPx.root iEPx.root
  touch eps_done_$EPS
  echo "EPSSCAN $EPS COMPLETE"
done
echo "(targets: WINDOWS 0.2284 0.3427 0.3409 | BUMP 1.1039 | STEP 1.1253 | CV 0.239 | px/ev 323799 | lowadc 92025 | islands 30006 | size 10.79)"
echo "EPS SCAN ALL COMPLETE"
