// sim_to_csv.C — export the v4.0 sim in the EXACT labeled_hits.csv schema
// (event,hitID,layer,phi,tbin,x,y,z,adc,label,trans_x,trans_y) + the
// per-event summary schema (event,total_active_hits,raw_islands,
// valid_classical_tracks). Labels are AUTOMATIC truth (vs the manual 82-hit
// real sample): track / looper (gpt<0.164) / noise (iso, minis, gpt<0).
// trans_x/y = (x/r, y/r) * layer — reverse-engineered from labeled_hits.csv
// rows (matches to 5 decimals).
// x,y from layer radius (transport GEO) x phi; z from the sim z branch.
#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>
#include <cstdio>
#include <map>
#include <set>
#include <cmath>

void sim_to_csv(int f_lo = 0, int f_hi = 9,
                const char *digi = "digi_frames_production_v40b.root",
                const char *frames = "frames_pau_production_v40b.root",
                const char *islands = "island_frames_v40b.root",
                const char *outhits = "sim_labeled_hits_v40.csv",
                const char *outsum = "sim_event_summary_v40.csv")
{
  // layer radii from the transport geometry (single source of truth)
  gROOT->ProcessLine(".L tpc_digitize.C+");
  gROOT->ProcessLine("TDG::loadGeo()");
  double R[55];
  for (int L = 7; L <= 54; ++L)
    R[L] = *(double *) gROOT->ProcessLineFast(Form("&TDG::GEO[%d].radius", L));

  // truth: (frame, id) -> gpt
  TFile *ff = TFile::Open(frames);
  TTree *ft = (TTree *) ff->Get("frame_truth");
  float tev, tid, tgpt;
  ft->SetBranchAddress("event", &tev);
  ft->SetBranchAddress("newid", &tid);
  ft->SetBranchAddress("gpt", &tgpt);
  std::map<std::pair<int, int>, float> gpt;
  for (Long64_t i = 0; i < ft->GetEntries(); ++i)
  {
    ft->GetEntry(i);
    if ((int) tev < f_lo || (int) tev > f_hi) continue;
    gpt[{(int) tev, (int) tid}] = tgpt;
  }

  TFile *fd = TFile::Open(digi);
  TTree *t = (TTree *) fd->Get("ntp_hit");
  float ev, lay, phi, tb, adc, z, trk;
  t->SetBranchStatus("*", 0);
  for (auto b : {"event", "layer", "phi", "tbin", "adc", "z", "gtrackID"}) t->SetBranchStatus(b, 1);
  t->SetBranchAddress("event", &ev);
  t->SetBranchAddress("layer", &lay);
  t->SetBranchAddress("phi", &phi);
  t->SetBranchAddress("tbin", &tb);
  t->SetBranchAddress("adc", &adc);
  t->SetBranchAddress("z", &z);
  t->SetBranchAddress("gtrackID", &trk);

  // dedupe: iso/mini injection can append a second row on an occupied
  // (pad,tbin) (808/67M in v40b) - merge adc, keep first row identity.
  t->SetBranchStatus("phibin", 1);
  t->SetBranchStatus("zelem", 1);
  float pbv, zev; t->SetBranchAddress("phibin", &pbv);
  t->SetBranchAddress("zelem", &zev);
  FILE *fo = fopen(outhits, "w");
  fprintf(fo, "event,hitID,layer,phi,tbin,x,y,z,adc,label,trans_x,trans_y\n");
  long hid = 0, ndup = 0;
  std::map<int, long> nhits;
  std::map<int, std::set<int>> trktracks;
  // NOTE: iso/mini blocks are APPENDED with event numbers revisited (readout
  // layout quirk) - so buffers must span the WHOLE tree per frame, not flush
  // on event transitions. Memory ~13 MB/frame: export big ranges in chunks.
  struct Row { double phi, tb, z, adc; int L; float trk; };
  std::map<int, std::map<ULong64_t, Row>> allbuf;
  auto flush = [&](int fr, std::map<ULong64_t, Row> &buf) {
    for (auto &kv : buf)
    {
      Row &r = kv.second;
      double x = R[r.L] * cos(r.phi), y = R[r.L] * sin(r.phi);
      const char *label = "noise";
      auto it = gpt.find({fr, (int) r.trk});
      if (it != gpt.end() && it->second >= 0)
      {
        label = (it->second < 0.164) ? "looper" : "track";
        if (it->second >= 0.164) trktracks[fr].insert((int) r.trk);
      }
      fprintf(fo, "%.1f,%ld.0,%.1f,%f,%.1f,%f,%f,%f,%.1f,%s,%f,%f\n",
              (double) (fr + 1), hid++, (double) r.L, r.phi, r.tb, x, y, r.z, r.adc,
              label, x / R[r.L] * r.L, y / R[r.L] * r.L);
      nhits[fr]++;
    }
    buf.clear();
  };
  for (Long64_t i = 0; i < t->GetEntries(); ++i)
  {
    t->GetEntry(i);
    int fr = (int) ev;
    if (fr < f_lo || fr > f_hi) continue;
    ULong64_t k = ((ULong64_t)(int) lay << 25) | ((ULong64_t)(int) zev << 24) |
                  ((ULong64_t)(int) pbv << 12) | (ULong64_t)(int) tb;
    auto ins = allbuf[fr].emplace(k, Row{phi, tb, z, adc, (int) lay, trk});
    if (!ins.second) { ins.first->second.adc += adc; ndup++; }
  }
  for (auto &fb : allbuf) flush(fb.first, fb.second);
  fclose(fo);
  printf("dedupe: %ld duplicate pixels merged\n", ndup);

  // per-event summary (islands from the island tree)
  TFile *fi = TFile::Open(islands);
  TTree *ti = (TTree *) fi->Get("island");
  FILE *fs = fopen(outsum, "w");
  fprintf(fs, "event,total_active_hits,raw_islands,valid_classical_tracks\n");
  for (int fr = f_lo; fr <= f_hi; ++fr)
    fprintf(fs, "%.1f,%ld,%lld,%zu\n", (double) (fr + 1), nhits[fr],
            ti->GetEntries(Form("event==%d", fr)), trktracks[fr].size());
  fclose(fs);
  printf("SIMCSV frames %d-%d: %ld hit rows -> %s | summary -> %s\n",
         f_lo, f_hi, hid, outhits, outsum);
}
