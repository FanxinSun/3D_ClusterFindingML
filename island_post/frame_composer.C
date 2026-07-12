// frame_composer.C — P3: compose real-structure streaming frames from a library of
// single-collision transported pixels (tpc_transport output), at ANY collision rate,
// without new Geant4. See PIPELINE.md "P3 LAUNCHED" section.
//
//   frame_composer(libs, out, nframes, rate_khz, seed, ev0, evals)
//     libs   : comma-separated raw_pix library files (one collision per `event`)
//     out    : raw_pix-format output (event = frame id) -> feed tpc_readout unchanged
//     ev0    : frame-id offset (batch dense frames: readout <=10 at a time, then hadd)
//     evals  : v2 TRUTH MODE — comma-separated source eval files (ntp_g4hit), SAME ORDER
//              as libs. When given: per-frame remap (source,trk) -> sequential newid>=1;
//              output trk = newid (dominant contributor per merged pixel), plus sidecar
//              tree `frame_truth` (event:newid:gpt:gflavor:gembed:gprimary) for
//              islandize91's sidecar mode. Empty evals -> v1 behaviour (trk=0, no sidecar).
//
// Frame model (P1-derived): 971 tbins x 53 ns = 51.46 us live window; collisions sampled
// Poisson over [-60 us, +51.46 us] — the 60 us pre-window carries drift drainage AND
// looper-tail history (13.2 us inverted the early/late ratio: 0.58 -> 0.89+).
// Collisions drawn uniformly WITH replacement; pixels shift by round(t0/53ns) tbins;
// duplicates on (layer,side,pad,tbin) merge by charge sum.

#include <TFile.h>
#include <TNtuple.h>
#include <TRandom3.h>
#include <TString.h>
#include <TTree.h>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <unordered_map>
#include <vector>

namespace FC
{
struct Pix
{
  uint8_t layer, side;
  uint16_t pad, tb;
  float q;
  int32_t trk;  // collision-local majority track id (from tpc_transport)
};
struct Contrib
{
  float q = 0, qbest = -1;
  int32_t best = 0;  // dominant newid
};
struct TRec
{
  float pt = -1, flavor = 0, embed = 0, primary = 0;
};
}  // namespace FC
using namespace FC;

void frame_composer(const char *libs = "raw_lib_a.root,raw_lib_b.root",
                    const char *out = "frames_raw.root",
                    int nframes = 100, double rate_khz = 240., int seed = 20260710,
                    int ev0 = 0, const char *evals = "",
                    const char *flashlib = "", double flash_prob = 0., double flash_scale = 1.,
                    const char *flash_spec = "", double flash_jitter_tb = 0., double flash_pixdisp = 0., double flash_stripedisp = 0., double flash_blur = 0.,
                    const char *rate_spec = "", double iso_scale = 0.0)
{
  const int NTB = 971;       // real frame: 971 tbins x 53 ns
  const double CLK = 0.053;  // us per tbin
  const double PRE = 60.0;   // us pre-window: drainage + looper-tail history
  const double SPAN = PRE + NTB * CLK;

  // ---- load pixel library ----
  std::vector<std::vector<Pix> > lib;
  std::vector<int> libsrc;  // source-file index per collision
  TString s(libs);
  auto *tok = s.Tokenize(",");
  for (int i = 0; i < tok->GetEntries(); ++i)
  {
    TFile *f = TFile::Open(tok->At(i)->GetName());
    if (!f || f->IsZombie())
    {
      printf("frame_composer: cannot open %s\n", tok->At(i)->GetName());
      continue;
    }
    TTree *t = (TTree *) f->Get("raw_pix");
    float ev, layer, side, pad, tb, q, trk;
    t->SetBranchAddress("event", &ev);
    t->SetBranchAddress("layer", &layer);
    t->SetBranchAddress("side", &side);
    t->SetBranchAddress("pad", &pad);
    t->SetBranchAddress("tbin", &tb);
    t->SetBranchAddress("q", &q);
    t->SetBranchAddress("trk", &trk);
    int cur = -1;
    size_t base = lib.size();
    Long64_t N = t->GetEntries();
    for (Long64_t k = 0; k < N; ++k)
    {
      t->GetEntry(k);
      if ((int) ev != cur)
      {
        cur = (int) ev;
        lib.push_back({});
        lib.back().reserve(500000);
        libsrc.push_back(i);
      }
      lib.back().push_back({(uint8_t) layer, (uint8_t) side, (uint16_t) pad, (uint16_t) tb, q, (int32_t) trk});
    }
    printf("frame_composer: %s -> %zu collisions (running total %zu)\n",
           tok->At(i)->GetName(), lib.size() - base, lib.size());
    f->Close();
  }
  if (lib.empty())
  {
    printf("frame_composer: empty library\n");
    return;
  }

  // ---- optional calibration-flash injection (e.g. CM laser): merged at fixed t0=0
  // (arrival phase baked into the flash lib via its g4hit times); trk from lib (-1 ->
  // no truth-table entry -> islandize91 labels these cls=2 automatically).
  std::vector<Pix> flash;
  if (flashlib && flashlib[0])
  {
    TFile *ff = TFile::Open(flashlib);
    TTree *t = ff ? (TTree *) ff->Get("raw_pix") : nullptr;
    if (t)
    {
      float ev, layer, side, pad, tb, q, trk;
      t->SetBranchAddress("event", &ev);
      t->SetBranchAddress("layer", &layer);
      t->SetBranchAddress("side", &side);
      t->SetBranchAddress("pad", &pad);
      t->SetBranchAddress("tbin", &tb);
      t->SetBranchAddress("q", &q);
      t->SetBranchAddress("trk", &trk);
      Long64_t N = t->GetEntries();
      for (Long64_t k = 0; k < N; ++k)
      {
        t->GetEntry(k);
        flash.push_back({(uint8_t) layer, (uint8_t) side, (uint16_t) pad, (uint16_t) tb,
                         (float) q, (int32_t) trk});  // scale applied at merge (per-flash sampling)
      }
      printf("frame_composer: flash lib %s -> %zu pixels, prob %.2f, scale %.2f\n",
             flashlib, flash.size(), flash_prob, flash_scale);
      ff->Close();
    }
  }
  // per-flash intensity sampling: "s1:w1,s2:w2,..." (weights normalized); empty -> flash_scale
  std::vector<std::pair<double, double> > fspec;
  if (flash_spec && flash_spec[0])
  {
    TString fs(flash_spec);
    auto *ft = fs.Tokenize(",");
    double wsum = 0;
    for (int i = 0; i < ft->GetEntries(); ++i)
    {
      TString p(ft->At(i)->GetName());
      auto *sw = p.Tokenize(":");
      fspec.push_back({atof(sw->At(0)->GetName()), atof(sw->At(1)->GetName())});
      wsum += fspec.back().second;
    }
    for (auto &q : fspec) q.second /= wsum;
    printf("frame_composer: flash intensity spec: %zu buckets\n", fspec.size());
  }

  // ---- v2: truth tables per source eval ----
  const bool TRUTH = evals && evals[0];
  std::vector<std::unordered_map<int, TRec> > ttab;
  if (TRUTH)
  {
    TString es(evals);
    auto *etok = es.Tokenize(",");
    for (int i = 0; i < etok->GetEntries(); ++i)
    {
      ttab.push_back({});
      TFile *f = TFile::Open(etok->At(i)->GetName());
      TTree *g = f ? (TTree *) f->Get("ntp_g4hit") : nullptr;
      if (!g)
      {
        printf("frame_composer: no ntp_g4hit in %s\n", etok->At(i)->GetName());
        continue;
      }
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
      for (Long64_t k = 0; k < N; ++k)
      {
        g->GetEntry(k);
        if (!std::isfinite(gid))
        {
          continue;
        }
        int id = (int) gid;
        if (ttab.back().count(id))
        {
          continue;
        }
        TRec r;
        r.pt = (std::isfinite(gpx) && std::isfinite(gpy)) ? std::sqrt(gpx * gpx + gpy * gpy) : -1;
        r.flavor = gfl;
        r.embed = gem;
        r.primary = gpr;
        ttab.back()[id] = r;
      }
      printf("frame_composer: truth table %d: %zu tracks (%s)\n", i, ttab.back().size(), etok->At(i)->GetName());
      f->Close();
    }
  }

  // ---- compose ----
  TFile *fo = new TFile(out, "RECREATE");
  TNtuple *op = new TNtuple("raw_pix", "composed frames", "event:layer:side:pad:tbin:q:trk");
  TNtuple *ot = TRUTH ? new TNtuple("frame_truth", "frame-scoped truth sidecar",
                                    "event:newid:gpt:gflavor:gembed:gprimary")
                      : nullptr;
  TRandom3 rng(seed);
  const double mu = rate_khz * 1e-3 * SPAN;
  // per-frame rate jitter: "khz1:w1,khz2:w2,..." sampled each frame (real rmbd varies
  // 196-362 kHz across frames; fixed-rate frames under-disperse px/frame: sigma/mu 0.29 vs 0.45)
  std::vector<std::pair<double, double> > rspec;
  if (rate_spec && rate_spec[0])
  {
    TString rs(rate_spec);
    auto *rt = rs.Tokenize(",");
    double wsum = 0;
    for (int i = 0; i < rt->GetEntries(); ++i)
    {
      TString p(rt->At(i)->GetName());
      auto *sw = p.Tokenize(":");
      rspec.push_back({atof(sw->At(0)->GetName()), atof(sw->At(1)->GetName())});
      wsum += rspec.back().second;
    }
    for (auto &q : rspec) q.second /= wsum;
    printf("frame_composer: rate spec: %zu buckets\n", rspec.size());
  }
  std::unordered_map<uint64_t, Contrib> merge;
  merge.reserve(1 << 22);
  std::unordered_map<uint64_t, int32_t> remap;  // (src<<32 | trk) -> newid
  std::vector<std::pair<uint64_t, int32_t> > remap_list;
  int32_t isoid = 0;  // isolated-background truth ids (src 0xFE)
  double totpx = 0, totcoll = 0;
  for (int fr = 0; fr < nframes; ++fr)
  {
    merge.clear();
    remap.clear();
    remap_list.clear();
    int32_t nextid = 1;
    double mu_f = mu;
    if (!rspec.empty())
    {
      double u = rng.Uniform(), acc = 0;
      for (auto &q : rspec) { acc += q.second; if (u <= acc) { mu_f = q.first * 1e-3 * SPAN; break; } }
    }
    int nc = rng.Poisson(mu_f);
    totcoll += nc;
    if (!flash.empty() && rng.Uniform() < flash_prob)
    {
      double fsc = flash_scale;
      if (!fspec.empty())
      {
        double u = rng.Uniform(), acc = 0;
        for (auto &q : fspec) { acc += q.second; if (u <= acc) { fsc = q.first; break; } }
      }
      int fdtb = (flash_jitter_tb > 0) ? (int) std::lround(rng.Gaus(0., flash_jitter_tb)) : 0;
      std::unordered_map<int32_t, double> sdisp;  // per-stripe factor, resampled per flash
      for (const Pix &p : flash)
      {
        int32_t nid = 0;
        if (TRUTH && p.trk != -9999)
        {
          uint64_t rk = (0xFFULL << 32U) | (uint32_t) p.trk;  // src=0xFF = flash
          auto it = remap.find(rk);
          if (it == remap.end())
          {
            it = remap.emplace(rk, nextid++).first;
            remap_list.push_back({rk, it->second});
          }
          nid = it->second;
        }
        int ftb = (int) p.tb + fdtb;
        if (ftb < 0 || ftb >= NTB)
        {
          continue;
        }
        double disp = 1.0;
        if (flash_stripedisp > 0)
        {  // per-stripe illumination nonuniformity (correlated within a stripe)
          auto it = sdisp.find(p.trk);
          if (it == sdisp.end())
          {
            it = sdisp.emplace(p.trk, std::exp(rng.Gaus(0., flash_stripedisp) -
                                               0.5 * flash_stripedisp * flash_stripedisp)).first;
          }
          disp *= it->second;
        }
        if (flash_pixdisp > 0)
        {  // sub-stripe optical speckle, per pixel
          disp *= std::exp(rng.Gaus(0., flash_pixdisp) - 0.5 * flash_pixdisp * flash_pixdisp);
        }
        double qq0 = p.q * fsc * disp;
        // optical halo: spread flash charge over neighbours (gaussian kernel, sigma=flash_blur
        // in pad AND tbin units) — real diffuse-laser light spills beyond the 1 mm stripes
        const int BW = (flash_blur > 0) ? 2 : 0;
        double ksum = 0;
        double kern[5][5];
        for (int dp = -BW; dp <= BW; ++dp)
        {
          for (int dt = -BW; dt <= BW; ++dt)
          {
            double w = (flash_blur > 0) ? std::exp(-0.5 * (dp * dp + dt * dt) / (flash_blur * flash_blur)) : 1.0;
            kern[dp + 2][dt + 2] = w;
            ksum += w;
          }
        }
        for (int dp = -BW; dp <= BW; ++dp)
        {
          for (int dt = -BW; dt <= BW; ++dt)
          {
            int tb2 = ftb + dt;
            if (tb2 < 0 || tb2 >= NTB)
            {
              continue;
            }
            int pd2 = (int) p.pad + dp;
            if (pd2 < 0)
            {
              continue;
            }
            double qq = qq0 * kern[dp + 2][dt + 2] / ksum;
            uint64_t key = ((uint64_t) p.layer << 40U) | ((uint64_t) p.side << 32U) |
                           ((uint64_t) (uint16_t) pd2 << 16U) | (uint64_t) tb2;
            Contrib &c = merge[key];
            c.q += qq;
            if (qq > c.qbest)
            {
              c.qbest = qq;
              c.best = nid;
            }
          }
        }
      }
    }
    for (int ic = 0; ic < nc; ++ic)
    {
      size_t ci = rng.Integer(lib.size());
      const auto &coll = lib[ci];
      int src = libsrc[ci];
      double t0 = -PRE + rng.Uniform() * SPAN;
      int dtb = (int) std::lround(t0 / CLK);
      for (const Pix &p : coll)
      {
        int ntb = (int) p.tb + dtb;
        if (ntb < 0 || ntb >= NTB)
        {
          continue;
        }
        int32_t nid = 0;
        if (TRUTH && p.trk != -9999)
        {
          uint64_t rk = ((uint64_t) (uint32_t) src << 32U) | (uint32_t) p.trk;
          auto it = remap.find(rk);
          if (it == remap.end())
          {
            it = remap.emplace(rk, nextid++).first;
            remap_list.push_back({rk, it->second});
          }
          nid = it->second;
        }
        uint64_t key = ((uint64_t) p.layer << 40U) | ((uint64_t) p.side << 32U) |
                       ((uint64_t) p.pad << 16U) | (uint64_t) ntb;
        Contrib &c = merge[key];
        c.q += p.q;
        if (p.q > c.qbest)
        {
          c.qbest = p.q;
          c.best = nid;
        }
      }
    }
    // ---- isolated-hit background injection (v3.5): real frames carry ~2.5k/frame
    // SINGLE-PIXEL islands (83% diffuse over ~93k pads, 17% on ~100 hot pads; adc just
    // above the region threshold, ~exponential falloff) that the G4 event content lacks
    // (low-energy background + hot channels). Injected at the composer so the truth
    // sidecar labels them: src 0xFE ids stay unmatched -> pt<0 -> cls=2, like the flash.
    // Rates/spectra measured on run 79507 island singles (2026-07-12):
    //   R1 789/fr thr 11 tau 4.5 | R2 588/fr thr 21 tau 6.5 | R3 1143/fr thr 21 tau 7.0
    // adu -> raw q via the NOMINAL v3.4 response chain (0.93 * greg * MV/e * ADU/mV);
    // readout's per-pad gain spread and noise add the natural variation on top.
    if (iso_scale > 0)
    {
      static const double IRATE[3] = {974., 460., 1491.};  // injected rates: measured real
      // singles (789/588/1143) x pilot-2 yield correction (some injected hits merge
      // into clusters and stop being singles: effective yield 0.81/1.11/0.77)
      static const double ITHR[3] = {11., 21., 21.};
      static const double ITAU[3] = {4.5, 6.5, 7.0};
      static const double IGREG[3] = {1.24, 1.0, 1.06};
      static const int NPADS[3] = {1128, 1536, 2304};
      const double Q_PER_ADU = 1.0 / (0.93 * 7.68e-3 * (1024. / 2200.));
      static std::vector<std::array<int, 3>> hotpads;  // (layer, side, pad), hash-static
      if (hotpads.empty())
      {
        TRandom3 hr(20260712);
        for (int i = 0; i < 100; ++i)
        {
          int rg = (int) (hr.Uniform() * 3);
          hotpads.push_back({7 + rg * 16 + (int) hr.Integer(16), (int) hr.Integer(2), (int) hr.Integer(NPADS[rg])});
        }
      }
      auto inject = [&](int lay, int side, int pad, int rg) {
        double adu = ITHR[rg] + rng.Exp(ITAU[rg]);
        double q = adu * Q_PER_ADU / IGREG[rg];
        int tb = (int) rng.Integer(NTB);
        int32_t nid = 0;
        if (TRUTH)
        {
          uint64_t rk = (0xFEULL << 32U) | (uint32_t) isoid++;
          auto it = remap.emplace(rk, nextid++).first;
          remap_list.push_back({rk, it->second});
          nid = it->second;
        }
        uint64_t key = ((uint64_t) lay << 40U) | ((uint64_t) side << 32U) |
                       ((uint64_t) (uint16_t) pad << 16U) | (uint64_t) tb;
        Contrib &c = merge[key];
        c.q += q;
        if (q > c.qbest)
        {
          c.qbest = q;
          c.best = nid;
        }
      };
      for (int rg = 0; rg < 3; ++rg)
      {
        int nd = rng.Poisson(0.83 * IRATE[rg] * iso_scale);
        for (int i = 0; i < nd; ++i)
        {
          inject(7 + rg * 16 + (int) rng.Integer(16), (int) rng.Integer(2), (int) rng.Integer(NPADS[rg]), rg);
        }
      }
      int nh = rng.Poisson(0.17 * (IRATE[0] + IRATE[1] + IRATE[2]) * iso_scale);
      for (int i = 0; i < nh; ++i)
      {
        auto &hp = hotpads[rng.Integer(100)];
        int rg = (hp[0] - 7) / 16;
        inject(hp[0], hp[1], hp[2], rg);
      }
    }
    for (const auto &kv : merge)
    {
      op->Fill((float) (fr + ev0), (float) ((kv.first >> 40U) & 0xFF), (float) ((kv.first >> 32U) & 0xFF),
               (float) ((kv.first >> 16U) & 0xFFFF), (float) (kv.first & 0xFFFF),
               kv.second.q, (float) kv.second.best);
    }
    if (TRUTH)
    {
      for (const auto &rl : remap_list)
      {
        int src = (int) (rl.first >> 32U);
        int32_t trk = (int32_t) (uint32_t) (rl.first & 0xFFFFFFFFULL);
        TRec r;
        if (src < (int) ttab.size())
        {
          auto it = ttab[src].find(trk);
          if (it != ttab[src].end())
          {
            r = it->second;
          }
        }
        ot->Fill((float) (fr + ev0), (float) rl.second, r.pt, r.flavor, r.embed, r.primary);
      }
    }
    totpx += merge.size();
  }
  printf("frame_composer: %d frames @ %.1f kHz -> <%.1f> collisions/frame, <%.2fM> raw px/frame%s -> %s\n",
         nframes, rate_khz, totcoll / nframes, totpx / nframes / 1e6,
         TRUTH ? " [TRUTH sidecar]" : "", out);
  fo->cd();
  op->Write();
  if (ot)
  {
    ot->Write();
  }
  fo->Close();
}
