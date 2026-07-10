// islandize91.C — Day-3 (P2): labeled, real-format island clusters.
//
//   islandize91(pixels.root, out.root, isSim, truthsrc.root)
//     pixels : digi_cal.root (detached digitizer, has per-pixel gtrackID)
//              or the real ntuplizer file (no truth)
//     out    : TWO trees, row-aligned 1:1
//       ntp_cluster — EXACT real-data 91-branch layout (event/cluster/info blocks,
//                     same names+order as production & the container port)
//       ntp_truth   — event:iclus:gtrackID:purity:gpt:gflavor:gembed:gprimary:cls:ntrks
//                     cls: 0=track  1=looper(truth pT<0.164 GeV/c: curls before r_outer)
//                          2=noise/unmatched   -1=real data (no truth exists)
//     truthsrc (sim only): eval file with ntp_g4hit -> trackID -> (pT,flavor,embed,primary)
//
// Same 8-connected flood-fill as islandize.C (the validated comparison tool, left intact).
// Columns without meaning in the detached frame are NaN/0 by design and documented:
//   locx/locy (Acts-surface locals), fee/chan/sampa (electronics IDs), ex/ey/pez/pephi,
//   thick/afac/bfac, trackID/niter (reco tracks), GL1/occ scaler block.
// ez/ephi use the v5-style base (r^2*var_phi/(adc*0.14)) for semantic consistency
// with the container port.
//
// usage:
//   root -b -q 'islandize91.C+("digi_cal.root","island91_sim.root",1,
//                "/home/rog/sPHENIX/3D_ClusterFindingML/macros-offline/detectors/sPHENIX/exam5_g4svtx_eval.root")'
//   root -b -q 'islandize91.C+("<real ntuplizer file>","island91_real.root",0,"")'

#include <TFile.h>
#include <TNtuple.h>
#include <TTree.h>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <map>
#include <unordered_map>
#include <vector>

namespace I91
{
const double CLOCK = 53.0, VD = 0.0080, PT_LOOP = 0.164;

struct Lay
{
  int nbins = 0;
  double radius = 0, slope = 0, phi0 = 0;
};
Lay GEO[55];
bool loadGeo()
{
  FILE *g = fopen("tpc_geom_table.txt", "r");
  if (!g)
  {
    printf("ERROR: tpc_geom_table.txt missing (run geomfit.C)\n");
    return false;
  }
  char line[256];
  while (fgets(line, 256, g))
  {
    int L, nb;
    double r, sl, p0, p1;
    if (line[0] == '#')
    {
      continue;
    }
    if (sscanf(line, "%d %d %lf %lf %lf %lf", &L, &nb, &r, &sl, &p0, &p1) == 6)
    {
      GEO[L] = {nb, r, sl, p0};
    }
  }
  fclose(g);
  return true;
}

struct Pix
{
  int pad, tb;
  float adc, phi, z;
  int trk;
};

struct TruthRec
{
  float pt = -1, flavor = 0, embed = 0, primary = 0;
};

// 91-branch layout (identical strings to TrkrNtuplizerMin.cc / real production)
const char *VARLIST_CLUSTER =
    "event:seed:run:seg:job:"
    "locx:locy:x:y:z:r:phi:eta:theta:phibin:tbin:fee:chan:sampa:ex:ey:ez:ephi:pez:pephi:"
    "e:adc:maxadc:thick:afac:bfac:dcal:layer:phielem:zelem:size:phisize:zsize:pedge:redge:"
    "ovlp:trackID:niter:"
    "occ11:occ116:occ21:occ216:occ31:occ316:rawzdc:livezdc:scaledzdc:rawmbd:livembd:scaledmbd:"
    "rawmbdv10:livembdv10:scaledmbdv10:rawzdc1:livezdc1:scaledzdc1:rawmbd1:livembd1:scaledmbd1:"
    "rawmbdv101:livembdv101:scaledmbdv101:rzdc:rmbd:rmbdv10:bco1:bco:bcotr:bcotr1:"
    "ntrk:ntpcseed:nsiseed:nhitmvtx:nhitintt:nhittpot:nhittpcall:nhittpcin:nhittpcmid:nhittpcout:"
    "nclusall:nclustpc:nclustpcpos:nclustpcneg:nclusintt:nclusmaps:nclusmms";
enum
{
  iEVT = 0, iSEED, iRUN, iSEG, iJOB,
  iLOCX, iLOCY, iX, iY, iZ, iR, iPHI, iETA, iTHETA, iPHIBIN, iTBIN, iFEE, iCHAN, iSAMPA,
  iEX, iEY, iEZ, iEPHI, iPEZ, iPEPHI, iE, iADC, iMAXADC, iTHICK, iAFAC, iBFAC, iDCAL,
  iLAYER, iPHIELEM, iZELEM, iSIZE, iPHISIZE, iZSIZE, iPEDGE, iREDGE, iOVLP, iTRACKID, iNITER,
  iINFO0,  // first info column (occ11)
  NCOL = 91
};
// info offsets (relative to iINFO0)
enum
{
  jNTRK = 31, jNTPCSEED, jNSISEED, jNHITMVTX, jNHITINTT, jNHITTPOT,
  jNHITTPCALL, jNHITTPCIN, jNHITTPCMID, jNHITTPCOUT,
  jNCLUSALL, jNCLUSTPC, jNCLUSTPCPOS, jNCLUSTPCNEG, jNCLUSINTT, jNCLUSMAPS, jNCLUSMMS
};
}  // namespace I91
using namespace I91;

void islandize91(const char *pixels, const char *out, int isSim, const char *truthsrc = "", const char *sidecar = "")
{
  if (!loadGeo())
  {
    return;
  }

  // --- truth-track table (sim only): trackID -> (pT, flavor, embed, primary)
  std::unordered_map<int, TruthRec> ttab;
  if (isSim && truthsrc && truthsrc[0])
  {
    TFile *ft = TFile::Open(truthsrc);
    TTree *g = ft ? (TTree *) ft->Get("ntp_g4hit") : nullptr;
    if (g)
    {
      float gid, gpx, gpy, gfl, gem, gpr;
      g->SetBranchStatus("*", 0);
      for (const char *b : {"gtrackID", "gpx", "gpy", "gflavor", "gembed", "gprimary"})
      {
        g->SetBranchStatus(b, 1);
      }
      g->SetBranchAddress("gtrackID", &gid);
      g->SetBranchAddress("gpx", &gpx);
      g->SetBranchAddress("gpy", &gpy);
      g->SetBranchAddress("gflavor", &gfl);
      g->SetBranchAddress("gembed", &gem);
      g->SetBranchAddress("gprimary", &gpr);
      Long64_t N = g->GetEntries();
      for (Long64_t i = 0; i < N; ++i)
      {
        g->GetEntry(i);
        if (!std::isfinite(gid))
        {
          continue;
        }
        int id = (int) gid;
        if (ttab.count(id))
        {
          continue;
        }
        TruthRec t;
        t.pt = std::isfinite(gpx) && std::isfinite(gpy) ? std::sqrt(gpx * gpx + gpy * gpy) : -1;
        t.flavor = gfl;
        t.embed = gem;
        t.primary = gpr;
        ttab[id] = t;
      }
      printf("islandize91: truth table %zu tracks from %s\n", ttab.size(), truthsrc);
      ft->Close();
    }
  }

  // --- v2 sidecar: frame-scoped truth from frame_composer (key = event<<32 | newid)
  const bool SIDE = sidecar && sidecar[0];
  std::unordered_map<uint64_t, TruthRec> stab;
  if (SIDE)
  {
    TFile *fs = TFile::Open(sidecar);
    TTree *t = fs ? (TTree *) fs->Get("frame_truth") : nullptr;
    if (!t)
    {
      printf("islandize91: no frame_truth in %s\n", sidecar);
      return;
    }
    float ev, nid, gpt, gfl, gem, gpr;
    t->SetBranchAddress("event", &ev);
    t->SetBranchAddress("newid", &nid);
    t->SetBranchAddress("gpt", &gpt);
    t->SetBranchAddress("gflavor", &gfl);
    t->SetBranchAddress("gembed", &gem);
    t->SetBranchAddress("gprimary", &gpr);
    Long64_t N = t->GetEntries();
    for (Long64_t i = 0; i < N; ++i)
    {
      t->GetEntry(i);
      TruthRec r;
      r.pt = gpt;
      r.flavor = gfl;
      r.embed = gem;
      r.primary = gpr;
      stab[(((uint64_t) (uint32_t) ev) << 32U) | (uint32_t) (int32_t) nid] = r;
    }
    printf("islandize91: sidecar truth %zu (frame,id) entries from %s\n", stab.size(), sidecar);
    fs->Close();
  }

  // --- read pixels
  TFile *fi = TFile::Open(pixels);
  TTree *h = (TTree *) fi->Get("ntp_hit");
  float event, layer, phibin, tb, adc, side, phi, z, gtrk = -9999;
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
  if (isSim)
  {
    h->SetBranchStatus("gtrackID", 1);
    h->SetBranchAddress("gtrackID", &gtrk);
  }

  std::map<uint64_t, std::vector<Pix> > groups;  // (event,layer,side) — event-major order
  double nhit_tpc[4] = {0, 0, 0, 0};             // per-file placeholders (filled per event below)
  std::map<int, std::array<double, 4> > evhits;  // event -> {all,in,mid,out}
  Long64_t N = h->GetEntries();
  for (Long64_t i = 0; i < N; ++i)
  {
    h->GetEntry(i);
    if (layer < 7 || layer > 54 || adc <= 0)
    {
      continue;
    }
    if (!isSim && ((int) side) == 0)
    {
      z -= 105.5f;  // real-ntuple side0 z convention quirk (canon.h)
    }
    int sd = ((int) side == 1) ? 1 : 0;
    uint64_t key = ((uint64_t) (uint32_t) event << 24U) | ((uint64_t) (uint32_t) layer << 8U) | (uint64_t) sd;
    groups[key].push_back({(int) phibin, (int) tb, adc, phi, z, isSim && std::isfinite(gtrk) ? (int) gtrk : -9999});
    auto &eh = evhits[(int) event];
    eh[0]++;
    if ((int) layer == 7) eh[1]++;
    if ((int) layer == 30) eh[2]++;
    if ((int) layer == 54) eh[3]++;
  }
  printf("islandize91: %lld rows -> %zu groups, %zu events\n", N, groups.size(), evhits.size());

  TFile *fo = new TFile(out, "RECREATE");
  TNtuple *oc = new TNtuple("ntp_cluster", "island clusters, real-data 91-branch layout", VARLIST_CLUSTER);
  TNtuple *ot = new TNtuple("ntp_truth", "per-cluster truth labels (row-aligned with ntp_cluster)",
                            "event:iclus:gtrackID:purity:gpt:gflavor:gembed:gprimary:cls:ntrks");

  // per-event buffering so the info block carries complete per-event tallies
  struct BufRow
  {
    float c[NCOL];
    float t[10];
  };
  std::vector<BufRow> buf;
  int curev = -1;
  double nclus_side[2] = {0, 0};
  long ncl_tot = 0, cls_cnt[3] = {0, 0, 0};

  auto flush = [&](int ev) {
    auto &eh = evhits[ev];
    for (auto &b : buf)
    {
      b.c[iINFO0 + jNHITTPCALL] = eh[0];
      b.c[iINFO0 + jNHITTPCIN] = eh[1];
      b.c[iINFO0 + jNHITTPCMID] = eh[2];
      b.c[iINFO0 + jNHITTPCOUT] = eh[3];
      b.c[iINFO0 + jNCLUSALL] = (float) buf.size();
      b.c[iINFO0 + jNCLUSTPC] = (float) buf.size();
      b.c[iINFO0 + jNCLUSTPCPOS] = (float) nclus_side[1];
      b.c[iINFO0 + jNCLUSTPCNEG] = (float) nclus_side[0];
      oc->Fill(b.c);
      ot->Fill(b.t);
    }
    buf.clear();
    nclus_side[0] = nclus_side[1] = 0;
  };

  std::vector<int> stack, members;
  const float NaN = std::numeric_limits<float>::quiet_NaN();
  for (auto &g : groups)
  {
    int ev = (int) (g.first >> 24U);
    int L = (int) ((g.first >> 8U) & 0xFFFF);
    int sd = (int) (g.first & 0xFF);
    if (ev != curev)
    {
      if (curev >= 0)
      {
        flush(curev);
      }
      curev = ev;
    }
    auto &v = g.second;
    std::unordered_map<uint64_t, int> where;
    where.reserve(v.size() * 2);
    auto pk = [](int p, int t) { return ((uint64_t) (uint32_t) p << 20U) | (uint32_t) t; };
    for (size_t i = 0; i < v.size(); ++i)
    {
      where[pk(v[i].pad, v[i].tb)] = (int) i;
    }
    std::vector<char> used(v.size(), 0);
    for (size_t s = 0; s < v.size(); ++s)
    {
      if (used[s])
      {
        continue;
      }
      stack.assign(1, (int) s);
      used[s] = 1;
      members.clear();
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
      // --- summarize island
      double W = 0, cp = 0, ct = 0, cphi = 0, cz = 0, maxa = 0, vp = 0, vt = 0;
      int plo = 1 << 30, phe = -(1 << 30), tlo = 1 << 30, the = -(1 << 30);
      std::map<int, double> trkq;
      for (int m : members)
      {
        const Pix &p = v[m];
        W += p.adc;
        cp += p.adc * p.pad;
        ct += p.adc * p.tb;
        vp += p.adc * (double) p.pad * p.pad;
        vt += p.adc * (double) p.tb * p.tb;
        cphi += p.adc * p.phi;
        cz += p.adc * p.z;
        maxa = std::max(maxa, (double) p.adc);
        plo = std::min(plo, p.pad);
        phe = std::max(phe, p.pad);
        tlo = std::min(tlo, p.tb);
        the = std::max(the, p.tb);
        if (p.trk != -9999 && !(SIDE && p.trk == 0))
        {
          trkq[p.trk] += p.adc;
        }
      }
      if (W <= 0)
      {
        continue;
      }
      cp /= W;
      ct /= W;
      cphi /= W;
      cz /= W;
      vp = std::max(0., vp / W - cp * cp);
      vt = std::max(0., vt / W - ct * ct);
      const Lay &geo = GEO[L];
      double r = geo.radius;
      double x = r * std::cos(cphi), y = r * std::sin(cphi);
      double theta = std::atan2(r, (double) cz);
      double eta = -std::log(std::tan(0.5 * theta));
      // v5-style base errors for semantic consistency with the container port
      double pitch_rphi = geo.slope * r;
      double ephi = (phe == plo) ? std::sqrt(9. * pitch_rphi * pitch_rphi / 12.)
                                 : std::sqrt(r * r * (vp * geo.slope * geo.slope) / (W * 0.14));
      double ez = (the == tlo) ? std::sqrt(9. * (CLOCK * VD) * (CLOCK * VD) / 12.)
                               : std::sqrt(vt * (CLOCK * VD) * (CLOCK * VD) / (W * 0.14));
      int psz = phe - plo + 1, zsz = the - tlo + 1;
      int nper = geo.nbins / 12;
      int localpad = ((int) std::lround(cp)) % std::max(1, nper);
      float pedge = (localpad < psz / 2. || localpad > nper - 1 - psz / 2.) ? 1.f : 0.f;
      float redge = (L == 7 || L == 22 || L == 23 || L == 28 || L == 39 || L == 54) ? 1.f : 0.f;

      BufRow b;
      for (int i = 0; i < NCOL; ++i)
      {
        b.c[i] = 0;
      }
      b.c[iEVT] = (float) ev;
      b.c[iSEED] = 0;
      b.c[iLOCX] = NaN;
      b.c[iLOCY] = NaN;
      b.c[iX] = (float) x;
      b.c[iY] = (float) y;
      b.c[iZ] = (float) cz;
      b.c[iR] = (float) r;
      b.c[iPHI] = (float) cphi;
      b.c[iETA] = (float) eta;
      b.c[iTHETA] = (float) theta;
      b.c[iPHIBIN] = (float) cp;
      b.c[iTBIN] = (float) ct;
      b.c[iFEE] = NaN;
      b.c[iCHAN] = NaN;
      b.c[iSAMPA] = NaN;
      b.c[iEX] = NaN;
      b.c[iEY] = NaN;
      b.c[iEZ] = (float) ez;
      b.c[iEPHI] = (float) ephi;
      b.c[iPEZ] = NaN;
      b.c[iPEPHI] = NaN;
      b.c[iE] = (float) W;
      b.c[iADC] = (float) W;
      b.c[iMAXADC] = (float) maxa;
      b.c[iTHICK] = NaN;
      b.c[iAFAC] = NaN;
      b.c[iBFAC] = NaN;
      b.c[iDCAL] = 1;
      b.c[iLAYER] = (float) L;
      b.c[iPHIELEM] = (float) (((int) std::lround(cp)) / std::max(1, nper));  // TPC sector
      b.c[iZELEM] = (float) sd;
      b.c[iSIZE] = (float) members.size();
      b.c[iPHISIZE] = (float) psz;
      b.c[iZSIZE] = (float) zsz;
      b.c[iPEDGE] = pedge;
      b.c[iREDGE] = redge;
      b.c[iOVLP] = 3;
      b.c[iTRACKID] = NaN;
      b.c[iNITER] = 0;

      // --- truth
      int bestid = -9999;
      double bestq = 0, sumq = 0;
      for (auto &kv : trkq)
      {
        sumq += kv.second;
        if (kv.second > bestq)
        {
          bestq = kv.second;
          bestid = kv.first;
        }
      }
      float purity = (W > 0) ? (float) (bestq / W) : 0.f;
      float gpt = -1, gfl = 0, gem = 0, gpr = 0;
      int cls = isSim ? 2 : -1;  // default: noise/unmatched (sim), no-truth (real)
      if (isSim && bestid != -9999)
      {
        const TruthRec *tr = nullptr;
        if (SIDE)
        {
          auto it = stab.find((((uint64_t) (uint32_t) ev) << 32U) | (uint32_t) (int32_t) bestid);
          if (it != stab.end())
          {
            tr = &it->second;
          }
        }
        else
        {
          auto it = ttab.find(bestid);
          if (it != ttab.end())
          {
            tr = &it->second;
          }
        }
        if (tr)
        {
          gpt = tr->pt;
          gfl = tr->flavor;
          gem = tr->embed;
          gpr = tr->primary;
          // pt<0 = sentinel (no kinematics known, e.g. injected calibration flash) -> noise
          cls = (gpt < 0) ? 2 : ((gpt < PT_LOOP) ? 1 : 0);
        }
      }
      if (isSim && cls >= 0 && cls <= 2)
      {
        cls_cnt[cls]++;
      }
      float t[10] = {(float) ev, (float) buf.size(), (float) bestid, purity, gpt, gfl, gem, gpr,
                     (float) cls, (float) trkq.size()};
      for (int i = 0; i < 10; ++i)
      {
        b.t[i] = t[i];
      }
      nclus_side[sd]++;
      buf.push_back(b);
      ncl_tot++;
    }
  }
  if (curev >= 0)
  {
    flush(curev);
  }
  printf("islandize91: %s -> %s : %ld clusters", pixels, out, ncl_tot);
  if (isSim)
  {
    printf("  [labels: track %ld, looper %ld, noise %ld]", cls_cnt[0], cls_cnt[1], cls_cnt[2]);
  }
  printf("\n");
  fo->cd();
  oc->Write();
  ot->Write();
  fo->Close();
}
