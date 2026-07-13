// prodclus.C — P4: detached port of the PRODUCTION TpcClusterizer (coresoftware
// offline/packages/tpc/TpcClusterizer.cc, ana.331-era local clone), faithful to the
// set_rawdata_reco() preset semantics:
//   - per (event, layer, sector, side) hitset; clusters cannot cross sector edges
//   - adcval grid keeps pixels with adc > edge_threshold(3); seeds are adc > seed(5),
//     processed HIGHEST-ADC-FIRST; consumed pixels marked USHRT_MAX ("touch")
//   - remove_singles: seed skipped if its 8-neighbourhood adc sum is zero
//   - seeded rectangular growth: t-range at the seed pad, then phi-range per t row;
//     FixedWindow(3) short-circuits ranges to +-3 with edge truncation
//   - ESCAPE HATCH (verbatim): if the windowed cluster fills >= half the (2W+1)^2 box
//     or spans the full width/height, re-cluster with the window OFF -> half-size
//     caps (MaxClusterHalfSizePhi 10 / MaxClusterHalfSizeT 20) -> zsize cap spike 40,
//     phisize edge ~20-21
//   - size<=1 discarded; adc_sum < min_adc_sum(5) discarded
//   - phisize/zsize = bbox spans, adc = sum, maxadc = max (rounded), centroids
//     adc-weighted; min_err cut passes trivially (min_err_squared=0)
// Pedestal = 0: our inputs (real ntuplizer ntp_hit, sim digi ntp_hit) are already
// pedestal-subtracted, matching the online-ZS input the raw preset assumes.
//
// GATE (self-validating): prodclus on REAL pixels must reproduce real ntp_cluster
// (28,065 clusters/frame; zsize cap structure; phisize comb; size semantics).
#include <TFile.h>
#include <TNtuple.h>
#include <TTree.h>
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdio>
#include <map>
#include <vector>

namespace PC
{
struct Cfg
{
  double seed_thr = 5.0;
  double edge_thr = 3.0;
  double min_adc_sum = 5.0;
  int fixed_window = 3;
  int max_half_phi = 10;
  int max_half_t = 20;
  bool do_singles = true;   // remove isolated seeds
  bool do_split = false;
  int ntb = 971;
};

struct LayG
{
  int nbins = 0;
  double slope = 0, phi0 = 0;
};
LayG GEO[55];
bool geo_ok = false;
void loadGeo()
{
  if (geo_ok) return;
  FILE *g = fopen("tpc_geom_table.txt", "r");
  if (!g)
  {
    printf("ERROR: no tpc_geom_table.txt\n");
    return;
  }
  char line[256];
  while (fgets(line, 256, g))
  {
    int L, nb;
    double r, sl, p0, p1;
    if (line[0] == '#') continue;
    if (sscanf(line, "%d %d %lf %lf %lf %lf", &L, &nb, &r, &sl, &p0, &p1) == 6)
    {
      GEO[L].nbins = nb;
      GEO[L].slope = sl;
      GEO[L].phi0 = p0;
    }
  }
  fclose(g);
  geo_ok = true;
}

typedef std::vector<std::vector<unsigned short>> Grid;
struct IH
{
  int iphi, it;
  double adc;
};
// per-hitset (pad,tbin) -> trk lookup for sim truth attachment
typedef std::map<int, int> TrkMap;  // key = iphi*4096 + it

bool isolated(int iphi, int it, int NP, int NT, const Grid &a)
{
  int lo_p = std::max(iphi - 1, 0), hi_p = std::min(iphi + 1, NP - 1);
  int lo_t = std::max(it - 1, 0), hi_t = std::min(it + 1, NT - 1);
  long sum = 0;
  for (int p = lo_p; p <= hi_p; ++p)
    for (int t = lo_t; t <= hi_t; ++t)
    {
      if (p == iphi && t == it) continue;
      sum += a[p][t];
    }
  return sum == 0;
}

void find_t_range(int phibin, int tbin, const Cfg &c, int fixw, const Grid &a,
                  int &tdown, int &tup, int &touch, int &edge)
{
  const int NT = c.ntb;
  tup = 0;
  tdown = 0;
  if (fixw != 0)
  {
    tup = fixw;
    tdown = fixw;
    if (tbin + tup >= NT)
    {
      tup = NT - tbin - 1;
      edge++;
    }
    if ((tbin - tdown) <= 0)
    {
      tdown = tbin;
      edge++;
    }
    return;
  }
  for (int it = 0; it < c.max_half_t; it++)
  {
    int ct = tbin + it;
    if (ct <= 0 || ct >= NT)
    {
      edge++;
      break;
    }
    if (a[phibin][ct] <= 0) break;
    if (a[phibin][ct] == USHRT_MAX)
    {
      touch++;
      break;
    }
    tup = it;
  }
  for (int it = 0; it < c.max_half_t; it++)
  {
    int ct = tbin - it;
    if (ct <= 0 || ct >= NT)
    {
      edge++;
      break;
    }
    if (a[phibin][ct] <= 0) break;
    if (a[phibin][ct] == USHRT_MAX)
    {
      touch++;
      break;
    }
    tdown = it;
  }
}

void find_phi_range(int phibin, int tbin, const Cfg &c, int fixw, int NP, const Grid &a,
                    int &phidown, int &phiup, int &touch, int &edge)
{
  phidown = 0;
  phiup = 0;
  if (fixw != 0)
  {
    phiup = fixw;
    phidown = fixw;
    if (phibin + phiup >= NP)
    {
      phiup = NP - phibin - 1;
      edge++;
    }
    if (phibin - phidown <= 0)
    {
      phidown = phibin;
      edge++;
    }
    return;
  }
  for (int ip = 0; ip < c.max_half_phi; ip++)
  {
    int cp = phibin + ip;
    if (cp < 0 || cp >= NP)
    {
      edge++;
      break;
    }
    if (a[cp][tbin] <= 0) break;
    if (a[cp][tbin] == USHRT_MAX)
    {
      touch++;
      break;
    }
    phiup = ip;
  }
  for (int ip = 0; ip < c.max_half_phi; ip++)
  {
    int cp = phibin - ip;
    if (cp < 0 || cp >= NP)
    {
      edge++;
      break;
    }
    if (a[cp][tbin] <= 0) break;
    if (a[cp][tbin] == USHRT_MAX)
    {
      touch++;
      break;
    }
    phidown = ip;
  }
}

void get_cluster(int phibin, int tbin, const Cfg &c, int fixw, int NP, const Grid &a,
                 std::vector<IH> &out, int &touch, int &edge)
{
  int tup = 0, tdown = 0;
  find_t_range(phibin, tbin, c, fixw, a, tdown, tup, touch, edge);
  for (int it = tbin - tdown; it <= tbin + tup; it++)
  {
    int pu = 0, pd = 0;
    find_phi_range(phibin, it, c, fixw, NP, a, pd, pu, touch, edge);
    for (int ip = phibin - pd; ip <= phibin + pu; ip++)
    {
      if (a[ip][it] > 0 && a[ip][it] != USHRT_MAX)
      {
        if (c.do_singles && isolated(ip, it, NP, c.ntb, a)) continue;
        out.push_back({ip, it, (double) a[ip][it]});
      }
    }
  }
}
}  // namespace PC
using namespace PC;

// in: pixel tree (ntp_hit). isSim: 1 = sim digi (zbin=tbin already; no extra cut),
//                            0 = real ntuplizer (canonical layer/adc cut applied).
void prodclus(const char *in = "/home/rog/sPHENIX/3D_ClusterFindingML/clusters_seeds_island_79507-0.root_ntuplizer.root",
              const char *out = "prodclus_real.root", int isSim = 0,
              double seed_thr = 5.0, double edge_thr = 3.0, double min_adc_sum = 5.0,
              int fixed_window = 3, int max_half_phi = 10, int max_half_t = 20,
              int do_singles = 1,
              const char *maskfile = "", const char *truthfile = "")
{
  loadGeo();
  if (!geo_ok) return;
  Cfg cfg;
  cfg.seed_thr = seed_thr;
  cfg.edge_thr = edge_thr;
  cfg.min_adc_sum = min_adc_sum;
  cfg.fixed_window = fixed_window;
  cfg.max_half_phi = max_half_phi;
  cfg.max_half_t = max_half_t;
  cfg.do_singles = do_singles;

  // optional channel veto (production is_pad_masked semantics at ingestion):
  // deadmap-format payload (Multiple tree; Fx/Fy -> phi -> global pad, as in
  // tpc_digitize's loader). Production masked these in reco even where the
  // ntuplizer still carries their hits.
  std::map<long, char> mask;
  if (maskfile && maskfile[0])
  {
    TFile *fm = TFile::Open(maskfile);
    TTree *tm = fm ? (TTree *) fm->Get("Multiple") : nullptr;
    if (tm)
    {
      Int_t dl, dside;
      Float_t fx, fy;
      tm->SetBranchAddress("Ilayer", &dl);
      tm->SetBranchAddress("Iside", &dside);
      tm->SetBranchAddress("Fx", &fx);
      tm->SetBranchAddress("Fy", &fy);
      for (long k = 0; k < tm->GetEntries(); ++k)
      {
        tm->GetEntry(k);
        if (dl < 7 || dl > 54) continue;
        double phi = std::atan2((double) fy, (double) fx);
        int pad = (int) std::lround((phi - GEO[dl].phi0) / GEO[dl].slope);
        if (pad < 0) pad += GEO[dl].nbins;
        if (pad >= GEO[dl].nbins) pad -= GEO[dl].nbins;
        mask[((long) dl << 20) | ((long) dside << 16) | pad] = 1;
      }
      printf("prodclus: mask %s -> %zu channels\n", maskfile, mask.size());
      fm->Close();
    }
  }
  // optional truth sidecar (sim): frame_truth (event,newid) -> (pt,flavor,embed,primary)
  std::map<long long, std::array<float, 4>> truth;
  if (truthfile && truthfile[0])
  {
    TFile *ft = TFile::Open(truthfile);
    TTree *tt = ft ? (TTree *) ft->Get("frame_truth") : nullptr;
    if (tt)
    {
      float tev, tid, tpt, tfl, tem, tpr;
      tt->SetBranchAddress("event", &tev);
      tt->SetBranchAddress("newid", &tid);
      tt->SetBranchAddress("gpt", &tpt);
      tt->SetBranchAddress("gflavor", &tfl);
      tt->SetBranchAddress("gembed", &tem);
      tt->SetBranchAddress("gprimary", &tpr);
      for (Long64_t k = 0; k < tt->GetEntries(); ++k)
      {
        tt->GetEntry(k);
        truth[((long long) (int) tev << 32) | (uint32_t) (int) tid] = {tpt, tfl, tem, tpr};
      }
      printf("prodclus: truth sidecar %s -> %zu ids\n", truthfile, truth.size());
      ft->Close();
    }
  }
  TFile *fi = TFile::Open(in);
  TTree *t = (TTree *) fi->Get("ntp_hit");
  float ev, lay, pb, tb, adc, ze, trk = 0;
  t->SetBranchStatus("*", 0);
  for (auto b : {"event", "layer", "phibin", "adc", "zelem"}) t->SetBranchStatus(b, 1);
  if (isSim)
  {
    t->SetBranchStatus("gtrackID", 1);
    t->SetBranchAddress("gtrackID", &trk);
  }
  t->SetBranchStatus(isSim ? "zbin" : "tbin", 1);
  t->SetBranchAddress("event", &ev);
  t->SetBranchAddress("layer", &lay);
  t->SetBranchAddress("phibin", &pb);
  t->SetBranchAddress(isSim ? "zbin" : "tbin", &tb);
  t->SetBranchAddress("adc", &adc);
  t->SetBranchAddress("zelem", &ze);

  // bucket pixels per (event, layer, side, sector)
  std::map<long long, std::vector<std::array<int, 4>>> buckets;  // (localpad, tbin, adc, trk)
  Long64_t N = t->GetEntries();
  long skipped = 0;
  for (Long64_t i = 0; i < N; ++i)
  {
    t->GetEntry(i);
    int L = (int) lay;
    if (L < 7 || L > 54 || adc <= 0) continue;
    int nper = GEO[L].nbins / 12;
    int pad = (int) pb;
    if (pad < 0 || pad >= GEO[L].nbins)
    {
      skipped++;
      continue;
    }
    int sec = pad / nper;
    int side = (((int) ze) == 1) ? 1 : 0;
    if (!mask.empty() && mask.count(((long) L << 20) | ((long) side << 16) | pad))
    {
      skipped++;
      continue;
    }
    long long key = ((long long) (int) ev << 20) | ((long long) L << 12) | (side << 8) | sec;
    buckets[key].push_back({pad - sec * nper, (int) tb, (int) adc, isSim ? (int) trk : 0});
  }
  printf("prodclus: %lld rows -> %zu hitsets (%ld skipped)\n", N, buckets.size(), skipped);

  TFile *fo = new TFile(out, "RECREATE");
  TNtuple *o = new TNtuple("ntp_clus", "ported production clusters",
                           "event:layer:side:sector:phisize:zsize:adc:maxadc:npix:cphi:ctbin:nedge:ntouch");
  TNtuple *ot = (isSim && truthfile && truthfile[0])
                    ? new TNtuple("ntp_truth", "per-cluster truth (row-aligned)",
                                  "event:iclus:gtrackID:purity:gpt:gflavor:gembed:gprimary:cls:ntrks")
                    : nullptr;
  long nclus = 0, ntrk_lab = 0, nloop_lab = 0, nnoise_lab = 0;
  for (auto &kv : buckets)
  {
    int sec = (int) (kv.first & 0xFF);
    int side = (int) ((kv.first >> 8) & 0xF);
    int L = (int) ((kv.first >> 12) & 0xFF);
    int evt = (int) (kv.first >> 20);
    int NP = GEO[L].nbins / 12;
    Grid a(NP, std::vector<unsigned short>(cfg.ntb, 0));
    std::multimap<unsigned short, std::pair<int, int>> seeds;
    TrkMap trkmap;
    for (auto &px : kv.second)
    {
      if (px[0] < 0 || px[0] >= NP || px[1] < 0 || px[1] >= cfg.ntb) continue;
      unsigned short v = (unsigned short) px[2];
      if (v > cfg.seed_thr) seeds.insert({v, {px[0], px[1]}});
      if (v > cfg.edge_thr)
      {
        a[px[0]][px[1]] = v;
        if (ot) trkmap[px[0] * 4096 + px[1]] = px[3];
      }
    }
    auto remove_hit = [&](double adcv, int ip, int it) {
      a[ip][it] = USHRT_MAX;
      auto rng = seeds.equal_range((unsigned short) adcv);
      for (auto it2 = rng.first; it2 != rng.second; ++it2)
      {
        if (it2->second.first == ip && it2->second.second == it)
        {
          seeds.erase(it2);
          break;
        }
      }
    };
    while (!seeds.empty())
    {
      auto iter = std::prev(seeds.end());  // highest adc
      int ip = iter->second.first, it = iter->second.second;
      double sadc = iter->first;
      if (cfg.do_singles && isolated(ip, it, NP, cfg.ntb, a))
      {
        remove_hit(sadc, ip, it);
        continue;
      }
      std::vector<IH> lst;
      int ntouch = 0, nedge = 0;
      get_cluster(ip, it, cfg, cfg.fixed_window, NP, a, lst, ntouch, nedge);
      if (cfg.fixed_window > 0)
      {
        int plo = 1 << 30, phi2 = -(1 << 30), tlo = 1 << 30, thi = -(1 << 30);
        for (auto &h : lst)
        {
          plo = std::min(plo, h.iphi);
          phi2 = std::max(phi2, h.iphi);
          tlo = std::min(tlo, h.it);
          thi = std::max(thi, h.it);
        }
        int W = 2 * cfg.fixed_window + 1;
        if ((int) lst.size() > (int) (0.5 * W * W) ||
            (phi2 - plo + 1) >= W || (thi - tlo + 1) >= W)
        {
          lst.clear();
          get_cluster(ip, it, cfg, 0, NP, a, lst, ntouch, nedge);
        }
      }
      if (lst.size() <= 1)
      {
        for (auto &h : lst) remove_hit(h.adc, h.iphi, h.it);
        remove_hit(sadc, ip, it);
        continue;
      }
      // parameters
      double asum = 0, ipsum = 0, itsum = 0;
      int mx = 0, plo = 1 << 30, phi2 = -(1 << 30), tlo = 1 << 30, thi = -(1 << 30);
      for (auto &h : lst)
      {
        asum += h.adc;
        ipsum += h.iphi * h.adc;
        itsum += h.it * h.adc;
        mx = std::max(mx, (int) std::lround(h.adc));
        plo = std::min(plo, h.iphi);
        phi2 = std::max(phi2, h.iphi);
        tlo = std::min(tlo, h.it);
        thi = std::max(thi, h.it);
      }
      if (asum >= cfg.min_adc_sum)
      {
        double cphi = GEO[L].phi0 + GEO[L].slope * (ipsum / asum + sec * NP);
        float row[13] = {(float) evt, (float) L, (float) side, (float) sec,
                         (float) (phi2 - plo + 1), (float) (thi - tlo + 1),
                         (float) asum, (float) mx, (float) lst.size(),
                         (float) cphi, (float) (itsum / asum), (float) nedge, (float) ntouch};
        o->Fill(row);
        if (ot)
        {
          // pixel-majority truth: dominant contributor by adc; islandize91 label rules
          std::map<int, double> peradc;
          for (auto &h : lst)
          {
            auto itk = trkmap.find(h.iphi * 4096 + h.it);
            peradc[(itk != trkmap.end()) ? itk->second : -999999] += h.adc;
          }
          int dom = -999999, ntrks = 0;
          double domadc = -1;
          for (auto &pa : peradc)
          {
            ntrks++;
            if (pa.second > domadc)
            {
              domadc = pa.second;
              dom = pa.first;
            }
          }
          float gpt = -1, gfl = 0, gem = 0, gpr = 0;
          auto itt = truth.find(((long long) evt << 32) | (uint32_t) dom);
          if (itt != truth.end())
          {
            gpt = itt->second[0];
            gfl = itt->second[1];
            gem = itt->second[2];
            gpr = itt->second[3];
          }
          int cls = (gpt < 0) ? 2 : (gpt < 0.164 ? 1 : 0);
          if (cls == 0) ntrk_lab++;
          else if (cls == 1) nloop_lab++;
          else nnoise_lab++;
          float trow[10] = {(float) evt, (float) nclus, (float) dom, (float) (domadc / asum),
                            gpt, gfl, gem, gpr, (float) cls, (float) ntrks};
          ot->Fill(trow);
        }
        nclus++;
      }
      for (auto &h : lst) remove_hit(h.adc, h.iphi, h.it);
      remove_hit(sadc, ip, it);
    }
  }
  printf("prodclus: %s -> %s : %ld clusters", in, out, nclus);
  if (ot)
  {
    printf("  [labels: track %ld, looper %ld, noise %ld]", ntrk_lab, nloop_lab, nnoise_lab);
  }
  printf("\n");
  fo->cd();
  o->Write();
  if (ot) ot->Write();
  fo->Close();
}
