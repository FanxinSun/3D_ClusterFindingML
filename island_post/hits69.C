// hits69.C — export the sim digi ntp_hit in the REAL ntuplizer's 69-branch
// ntp_hit layout (clusters_seeds_island_79507-0.root_ntuplizer.root), the hit-
// tree companion of islandize91: real-data scripts run on sim unmodified.
// Conventions measured on the real file (2026-07-20) and mirrored:
//   ecell = adc (exact in real TPC rows)      cellID = 0 (real: all 0)
//   zbin  = NaN (real TPC rows carry NaN)     run/seg/job = 0, seed = 666 (real)
//   z     = apparent z + 105.5*(zelem==0)     (real storage convention; the
//           canon.h real-side correction inverts this exactly)
//   phielem = 12 sectors, phi=0 -> 5, centers (5-n)*pi/6 (sim has no spoke
//           gaps, so gap-phi rows get the nearest sector)
//   hitID = opaque per-pixel key (real's is a detector key; ours packs
//           layer/zelem/phibin/tbin)          e = 0 (real e is truth-matched
//           energy; sim keeps no per-pixel truth energy)
//   per-event context: nhittpcall = TPC rows, nclusall = nclustpc = island91
//           clusters, nclustpcpos/neg = clusters by side (real pos/neg is
//           side-based: ev1 check 2981/3620 vs z-sign 2976/3625);
//           everything DAQ/tracking (scalers, bco, occ*, ntrk, nhittpcin/...)
//           = 0 — no source in the detached sim (real semantics of
//           nhittpcin/mid/out are NOT region hit counts; do not fake).
// Duplicate pixels (iso/mini appended-block collisions, 808/67.3M in v40b)
// are adc-merged in place, consistent with the CSV exporters.
// Output is EVENT-CONTIGUOUS (real ntp_hit is; the digi order interleaves
// appended background blocks). The sort is native here as of 2026-07-22 —
// in the 2026-07-20 artifact it was a post-hoc step (user review finding);
// content identical, so the sealed hit69 md5 stands. Note: regeneration is
// content-identical but never byte-identical (ROOT embeds timestamps).
#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TMath.h>
#include <TROOT.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <unordered_map>
#include <vector>

void hits69(const char *digi = "digi_frames_production_v40b.root",
            const char *islands = "island91_frames_production_v40b.root",
            const char *out = "hit69_frames_production_v40b.root")
{
  gROOT->ProcessLine(".L tpc_digitize.C+");
  gROOT->ProcessLine("TDG::loadGeo();");
  double R[55] = {0};
  for (int L = 7; L <= 54; ++L)
    R[L] = *(double *) gROOT->ProcessLineFast(Form("&TDG::GEO[%d].radius;", L));
  if (R[7] < 1) { printf("GEO LOAD FAILED\n"); return; }

  // per-event cluster counts from island91 (nclusall/nclustpc/pos/neg)
  std::map<int, long> nclus, ncpos, ncneg;
  {
    TFile *fi = TFile::Open(islands);
    TTree *c = (TTree *) fi->Get("ntp_cluster");
    float cev, czel;
    c->SetBranchStatus("*", 0);
    c->SetBranchStatus("event", 1); c->SetBranchAddress("event", &cev);
    c->SetBranchStatus("zelem", 1); c->SetBranchAddress("zelem", &czel);
    for (Long64_t i = 0; i < c->GetEntries(); ++i)
    {
      c->GetEntry(i);
      nclus[(int) cev]++;
      if ((int) czel == 1) ncpos[(int) cev]++; else ncneg[(int) cev]++;
    }
    fi->Close();
  }

  TFile *fd = TFile::Open(digi);
  TTree *t = (TTree *) fd->Get("ntp_hit");
  float ev, lay, phi, tb, adc, zel, z, pbv;
  t->SetBranchStatus("*", 0);
  for (auto b : {"event", "layer", "phi", "tbin", "adc", "zelem", "z", "phibin"})
    t->SetBranchStatus(b, 1);
  t->SetBranchAddress("event", &ev);
  t->SetBranchAddress("layer", &lay);
  t->SetBranchAddress("phi", &phi);
  t->SetBranchAddress("tbin", &tb);
  t->SetBranchAddress("adc", &adc);
  t->SetBranchAddress("zelem", &zel);
  t->SetBranchAddress("z", &z);
  t->SetBranchAddress("phibin", &pbv);
  const Long64_t NE = t->GetEntries();
  auto key = [](float lay, float zel, float pbv, float tb) {
    return ((ULong64_t)(int) lay << 25) | ((ULong64_t)(int) zel << 24) |
           ((ULong64_t)(int) pbv << 12) | (ULong64_t)(int) tb;
  };

  // pass 1: per-event (key,adc) lists -> duplicate keys with merged adc;
  // per-event deduped row count (nhittpcall); per-event ENTRY index so pass 2
  // can write event-contiguous output (events ascending, original order within)
  std::map<int, std::vector<std::pair<ULong64_t, float>>> ka;
  std::map<int, std::vector<Long64_t>> ent;
  for (Long64_t i = 0; i < NE; ++i)
  {
    t->GetEntry(i);
    if (adc <= 0) continue;
    ka[(int) ev].emplace_back(key(lay, zel, pbv, tb), adc);
    ent[(int) ev].push_back(i);
  }
  std::map<int, std::unordered_map<ULong64_t, float>> dupsum;  // ev -> key -> merged adc
  std::map<int, long> nhit;
  long ndup = 0;
  for (auto &kv : ka)
  {
    auto &v = kv.second;
    std::sort(v.begin(), v.end());
    long uniq = 0;
    for (size_t i = 0; i < v.size();)
    {
      size_t j = i + 1;
      float s = v[i].second;
      while (j < v.size() && v[j].first == v[i].first) { s += v[j].second; ++j; }
      if (j - i > 1) { dupsum[kv.first][v[i].first] = s; ndup += j - i - 1; }
      uniq++; i = j;
    }
    nhit[kv.first] = uniq;
  }
  ka.clear();
  printf("pass1: %lld rows, %ld duplicate pixels to merge\n", NE, ndup);

  // pass 2: stream, merge dups at first occurrence, fill the 69-branch layout
  const char *LAYOUT =
      "event:seed:run:seg:job:hitID:e:adc:layer:phielem:zelem:cellID:ecell:"
      "phibin:zbin:tbin:phi:r:x:y:z:occ11:occ116:occ21:occ216:occ31:occ316:"
      "rawzdc:livezdc:scaledzdc:rawmbd:livembd:scaledmbd:rawmbdv10:livembdv10:"
      "scaledmbdv10:rawzdc1:livezdc1:scaledzdc1:rawmbd1:livembd1:scaledmbd1:"
      "rawmbdv101:livembdv101:scaledmbdv101:rzdc:rmbd:rmbdv10:bco1:bco:bcotr:"
      "bcotr1:ntrk:ntpcseed:nsiseed:nhitmvtx:nhitintt:nhittpot:nhittpcall:"
      "nhittpcin:nhittpcmid:nhittpcout:nclusall:nclustpc:nclustpcpos:"
      "nclustpcneg:nclusintt:nclusmaps:nclusmms";
  TFile *fo = TFile::Open(out, "RECREATE");
  fo->SetCompressionLevel(5);
  TNtuple *nt = new TNtuple("ntp_hit", "sim hits, real-data 69-branch layout", LAYOUT);
  const float NaN = std::numeric_limits<float>::quiet_NaN();
  std::map<int, std::unordered_map<ULong64_t, bool>> written;  // dup keys already emitted
  float b[69];
  long nout = 0;
  for (auto &epair : ent)
  for (Long64_t i : epair.second)
  {
    t->GetEntry(i);
    if (adc <= 0) continue;
    int e = (int) ev, L = (int) lay;
    ULong64_t k = key(lay, zel, pbv, tb);
    float a = adc;
    auto de = dupsum.find(e);
    if (de != dupsum.end())
    {
      auto dk = de->second.find(k);
      if (dk != de->second.end())
      {
        if (written[e].count(k)) continue;  // later duplicate: already merged
        written[e][k] = true;
        a = dk->second;
      }
    }
    int sec = 5 - (int) TMath::Nint(phi * 6.0 / TMath::Pi());
    if (sec < 0) sec += 12;
    if (sec > 11) sec -= 12;
    for (int m = 0; m < 69; ++m) b[m] = 0;
    b[0] = ev;                       // event
    b[1] = 666;                      // seed (real sentinel)
    b[5] = (float) k;                // hitID (opaque pixel key)
    b[7] = a;                        // adc (dup-merged)
    b[8] = lay;
    b[9] = sec;                      // phielem
    b[10] = zel;
    b[12] = a;                       // ecell == adc (real convention)
    b[13] = pbv;                     // phibin
    b[14] = NaN;                     // zbin (real TPC rows: NaN)
    b[15] = tb;
    b[16] = phi;
    b[17] = R[L];                    // r
    b[18] = R[L] * TMath::Cos(phi);  // x
    b[19] = R[L] * TMath::Sin(phi);  // y
    b[20] = z + 105.5f * ((int) zel == 0);  // real storage convention
    b[58] = nhit[e];                 // nhittpcall
    b[62] = nclus[e] ? nclus[e] : 0; // nclusall
    b[63] = b[62];                   // nclustpc
    b[64] = ncpos[e];
    b[65] = ncneg[e];
    nt->Fill(b);
    nout++;
  }
  nt->Write();
  fo->Close();
  printf("hits69: %ld rows (%ld merged) -> %s [%d branches]\n",
         nout, ndup, out, 69);
}
