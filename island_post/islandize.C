// islandize.C — standalone "island" clusterizer post-processor (option D, v1).
// One algorithm for BOTH real and sim hit ntuples -> clustering drops out of any
// real-vs-sim comparison by construction.
//
//   island(in, out, isSim):
//     in : real  clusters_seeds_island_*_ntuplizer.root (ntp_hit: phibin/tbin, zelem=side)
//          sim   *_g4svtx_eval.root                     (ntp_hit: phibin/zbin,  zelem=side)
//     out: TTree "island" — one row per island cluster:
//          event, layer, side, nhits(size), adc(RAW sum), maxadc,
//          cphibin, ctbin (ADC-weighted centroids), phisize, zsize (bbox),
//          cphi, cz (ADC-weighted means of per-hit phi / z as stored in the input),
//          asym (ADC asymmetry along tbin: (late-early)/sum), rho (pixel fill of bbox)
//
// Pixels are grouped per (event, layer, side) with 8-connectivity. TPC layers only.
// Run: root -b -q 'islandize.C("in.root","out.root",1)'

#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>
#include <cmath>
#include <cstdint>
#include <map>
#include <unordered_map>
#include <vector>

namespace
{
struct Pix
{
  int pad, tb;
  float adc, phi, z;
};
}  // namespace

void island(const char *in, const char *out, bool isSim)
{
  TFile *fi = TFile::Open(in);
  TTree *h = (TTree *) fi->Get("ntp_hit");
  // use highest-cycle tree (files carry two cycles; Get returns highest)
  float event, layer, phibin, tb, adc, side, phi, z;
  h->SetBranchStatus("*", 0);
  for (const char *b : {"event", "layer", "phibin", "adc", "zelem", "phi", "z"})
  {
    h->SetBranchStatus(b, 1);
  }
  h->SetBranchStatus(isSim ? "zbin" : "tbin", 1);
  h->SetBranchAddress("event", &event);
  h->SetBranchAddress("layer", &layer);
  h->SetBranchAddress("phibin", &phibin);
  h->SetBranchAddress(isSim ? "zbin" : "tbin", &tb);
  h->SetBranchAddress("adc", &adc);
  h->SetBranchAddress("zelem", &side);
  h->SetBranchAddress("phi", &phi);
  h->SetBranchAddress("z", &z);

  TFile *fo = new TFile(out, "RECREATE");
  TNtuple *t = new TNtuple("island", "island clusters (post-processed)",
                           "event:layer:side:size:adc:maxadc:cphibin:ctbin:phisize:zsize:cphi:cz:asym:rho");

  // group pixels by (event, layer, side)
  std::map<uint64_t, std::vector<Pix> > groups;
  const Long64_t N = h->GetEntries();
  for (Long64_t i = 0; i < N; ++i)
  {
    h->GetEntry(i);
    if (layer < 7 || layer > 54 || adc <= 0)
    {
      continue;
    }
    int sd = ((int) side == 1) ? 1 : 0;
    uint64_t key = ((uint64_t) (uint32_t) event << 24U) | ((uint64_t) (uint32_t) layer << 8U) | (uint64_t) sd;
    // real-ntuple convention quirk (canon.h): side0 z = v*t, offset +105.5 vs physical
    float zc = (!isSim && ((int) side) == 0) ? z - 105.5f : z;
    groups[key].push_back({(int) phibin, (int) tb, adc, phi, zc});
    if (i % 5000000 == 0)
    {
      printf("  read %lld/%lld hits (%zu groups)\n", i, N, groups.size());
    }
  }
  printf("read done: %lld hits, %zu (event,layer,side) groups\n", N, groups.size());

  long nclus = 0;
  for (auto &g : groups)
  {
    auto &v = g.second;
    const float ev = (float) (g.first >> 24U);
    const float ly = (float) ((g.first >> 8U) & 0xFFFF);
    const float sd = (float) (g.first & 0xFF);
    // pixel lookup for this group
    std::unordered_map<uint64_t, int> where;
    where.reserve(v.size() * 2);
    auto pk = [](int p, int tbn) { return ((uint64_t) (uint32_t) p << 20U) | (uint32_t) tbn; };
    for (size_t i = 0; i < v.size(); ++i)
    {
      where[pk(v[i].pad, v[i].tb)] = (int) i;
    }
    std::vector<char> used(v.size(), 0);
    std::vector<int> stack;
    for (size_t s = 0; s < v.size(); ++s)
    {
      if (used[s])
      {
        continue;
      }
      // flood fill (8-connectivity)
      stack.assign(1, (int) s);
      used[s] = 1;
      std::vector<int> members;
      while (!stack.empty())
      {
        int c = stack.back();
        stack.pop_back();
        members.push_back(c);
        for (int dp = -1; dp <= 1; ++dp)
        {
          for (int dt = -1; dt <= 1; ++dt)
          {
            if (!dp && !dt)
            {
              continue;
            }
            auto it = where.find(pk(v[c].pad + dp, v[c].tb + dt));
            if (it != where.end() && !used[it->second])
            {
              used[it->second] = 1;
              stack.push_back(it->second);
            }
          }
        }
      }
      // summarize
      double W = 0, cp = 0, ct = 0, cphi = 0, cz = 0, maxa = 0;
      int plo = 1 << 30, phe = -(1 << 30), tlo = 1 << 30, the = -(1 << 30);
      for (int m : members)
      {
        const Pix &p = v[m];
        W += p.adc;
        cp += p.adc * p.pad;
        ct += p.adc * p.tb;
        cphi += p.adc * p.phi;
        cz += p.adc * p.z;
        maxa = std::max(maxa, (double) p.adc);
        plo = std::min(plo, p.pad);
        phe = std::max(phe, p.pad);
        tlo = std::min(tlo, p.tb);
        the = std::max(the, p.tb);
      }
      if (W <= 0)
      {
        continue;
      }
      cp /= W;
      ct /= W;
      // ADC asymmetry along time: (late - early)/(total), early/late split at centroid
      double early = 0, late = 0;
      for (int m : members)
      {
        (v[m].tb < ct ? early : late) += v[m].adc;
      }
      const float psz = phe - plo + 1, zsz = the - tlo + 1;
      float row[14] = {ev, ly, sd, (float) members.size(), (float) W, (float) maxa,
                       (float) cp, (float) ct, psz, zsz,
                       (float) (cphi / W), (float) (cz / W),
                       (float) ((late - early) / W),
                       (float) (members.size() / (psz * zsz))};
      t->Fill(row);
      nclus++;
    }
  }
  printf("islandize: %s -> %s : %ld clusters\n", in, out, nclus);
  fo->cd();
  t->Write();
  fo->Close();
}

void islandize(const char *in, const char *out, int isSim)
{
  island(in, out, isSim != 0);
}
