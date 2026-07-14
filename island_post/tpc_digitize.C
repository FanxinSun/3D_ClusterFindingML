// tpc_digitize.C — detached TPC response digitizer v2 (two-stage, calibration-ready).
//
//  STAGE A  tpc_transport(in, raw.root, NEV)
//    ntp_g4hit -> ionization -> drift+diffusion -> GEM gain -> zigzag pads ->
//    SAMPA Gamma4 time shares -> RAW pixel charge tree "raw_pix" (no electronics).
//    Expensive (~90 s/AuAu event) but runs ONCE.
//
//  STAGE B  tpc_readout(raw.root, out.root, gaincal, thr_adu, ret_pre, ret_post)
//    gain scale (exact: exponential draws scale linearly) -> mV -> ADU + ENC noise ->
//    clamp 1023 -> PEDESTAL SUBTRACT (74.4; real saturation bump at ~950 = 1023-ped) ->
//    ZS threshold + SAMPA-DSP pre/post-sample retention -> "ntp_hit" pixel tree
//    (islandize.C / hits_profile.C compatible). Seconds per iteration -> scan-friendly.
//
//  tpc_digitize(in, out, NEV) = A + B with day-1-compatible defaults.
//
// Physics sources: coresoftware master (see PIPELINE.md component table).
// usage examples:
//   root -b -q 'tpc_digitize.C+("eval.root","digi.root",2)'
//   root -l -b -q 'tpc_digitize.C+' -e 'tpc_transport("eval.root","raw2.root",2)'
//   root -l -b -q 'tpc_digitize.C+' -e 'tpc_readout("raw2.root","digi.root",0.84,8,1,2)'

#include <TFile.h>
#include <TNtuple.h>
#include <TRandom3.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <unordered_map>
#include <vector>

namespace CFG
{
const double VDRIFT = 0.0080;    // cm/ns (real calibration)
const double HALFZ = 105.5;      // cm
const double CLOCK = 53.0;       // ns
const int NTBIN = 965;           // ~51.1 us window
const double TS = 45.0;          // ns EFFECTIVE peaking time — day-2 calibrated vs real
                                 // run-length + near-threshold spectrum (master default 55;
                                 // 45 shortens the Gamma4 tail to match real runs at low thr)
const double EPG = 28.43e6;      // electrons per GeV
const double GAIN = 1400.;       // GEM mean gain (calibrate via gaincal in stage B)
const double CLOUD = 0.12;
const double DIFF_T = 0.005313;  // cm/sqrt(cm)
const double DIFF_L = 0.014596;  // cm/sqrt(cm)
const double N_SIGMA = 5;
const double MV_PER_E = 7.68e-3;          // mV per electron (peak, incl x2.4)
const double ADU_PER_MV = 1024. / 2200.;  // 10-bit over 2.2 V
const double NOISE_ADU = 670. * (7.68e-3 / 2.4) * (1024. / 2200.);  // ENC in ADU (~1.0)
const double PEDESTAL = 74.4;    // ADU (TpcClusterizer.h; real saturation bump 1023-ped~949)
const bool GAPS = true;
const double GAPPHASE = M_PI / 12.;
const double ACTIVE[3] = {0.5024, 0.5087, 0.5097};
}  // namespace CFG

namespace TDG
{
struct Lay
{
  int nbins;
  double radius, slope, phi0;
};
Lay GEO[55];
bool geo_ok = false;
// v3.5: transverse charge-cloud sigma, settable per tpc_transport call (CLOUD
// recalibration under the corrected regional ZS); defaults to the frozen CFG value.
double CLOUD_EFF = CFG::CLOUD;

void loadGeo()
{
  if (geo_ok)
  {
    return;
  }
  FILE *g = fopen("tpc_geom_table.txt", "r");
  if (!g)
  {
    printf("ERROR: no tpc_geom_table.txt — run geomfit.C first\n");
    return;
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
  geo_ok = true;
}

struct Cell
{
  float q = 0;
  float qbest = 0;
  int trk = -9999;
};

double g4pulse(double dt)
{
  if (dt <= 0)
  {
    return 0;
  }
  double x = dt / CFG::TS;
  return std::exp(-4 * x) * std::pow(x, 4.0);
}

const int NLUT = 212;
float LUT[NLUT][8];
void buildLUT()
{
  for (int iu = 0; iu < NLUT; ++iu)
  {
    double u = (iu + 0.5) * CFG::CLOCK / NLUT;
    double tot = 0;
    double tfe = CFG::CLOCK - u;
    LUT[iu][0] = 0.5 * g4pulse(tfe) * tfe;
    tot += LUT[iu][0];
    for (int k = 1; k < 8; ++k)
    {
      double s = 0;
      for (int is = 0; is < 6; ++is)
      {
        s += g4pulse(k * CFG::CLOCK - u + (is + 0.5) * CFG::CLOCK / 6.) * (CFG::CLOCK / 6.);
      }
      LUT[iu][k] = s;
      tot += s;
    }
    for (int k = 0; k < 8; ++k)
    {
      LUT[iu][k] = (tot > 0) ? LUT[iu][k] / tot : 0;
    }
  }
}

inline double gaus(double x, double s) { return std::exp(-0.5 * x * x / (s * s)) / (s * std::sqrt(2 * M_PI)); }
inline int region(int L) { return (L < 23) ? 0 : (L < 39 ? 1 : 2); }

bool inGap(double phi, int reg)
{
  if (!CFG::GAPS)
  {
    return false;
  }
  const double SEC = M_PI / 6.;
  double half_gap = 0.5 * (SEC - CFG::ACTIVE[reg]);
  double m = std::fmod(phi - CFG::GAPPHASE + 4 * M_PI, SEC);
  return (m < half_gap || m > SEC - half_gap);
}

void zigzag(int L, double phi, double sig, std::vector<int> &pads, std::vector<double> &share)
{
  pads.clear();
  share.clear();
  const Lay &g = GEO[L];
  const double radius = g.radius;
  const double step = g.slope;
  const double rphi = phi * radius;
  int blo = (int) std::floor((phi - (CFG::N_SIGMA * sig / radius) - step - g.phi0) / step);
  int bhi = (int) std::floor((phi + (CFG::N_SIGMA * sig / radius) + step - g.phi0) / step);
  int npads = bhi - blo;
  if (npads < 0 || npads > 9)
  {
    npads = 9;
  }
  const double pitch = step * radius;
  for (int ip = 0; ip <= npads; ++ip)
  {
    int pad = blo + ip;
    if (pad >= g.nbins)
    {
      pad -= g.nbins;
    }
    if (pad < 0)
    {
      pad += g.nbins;
    }
    double x = (g.phi0 + g.slope * pad) * radius - rphi;
    if (x > M_PI * radius)
    {
      x -= 2 * M_PI * radius;
    }
    if (x < -M_PI * radius)
    {
      x += 2 * M_PI * radius;
    }
    const double s = sig;
    double ov = (pitch - x) * (std::erf(x / (M_SQRT2 * s)) - std::erf((x - pitch) / (M_SQRT2 * s))) / (pitch * 2) +
                (pitch + x) * (std::erf((x + pitch) / (M_SQRT2 * s)) - std::erf(x / (M_SQRT2 * s))) / (pitch * 2) +
                (gaus(x - pitch, s) - gaus(x, s)) * s * s / pitch +
                (gaus(x + pitch, s) - gaus(x, s)) * s * s / pitch;
    if (ov > 1e-6)
    {
      pads.push_back(pad);
      share.push_back(ov);
    }
  }
}
}  // namespace TDG
using namespace TDG;

// ---------------- STAGE A: transport ----------------
void tpc_transport(const char *in, const char *rawout, int NEV = 2, double cloud = -1.0, double cloud_k = 0.0)
{
  CLOUD_EFF = (cloud > 0) ? cloud : CFG::CLOUD;
  loadGeo();
  if (!geo_ok)
  {
    return;
  }
  buildLUT();
  TRandom3 rng(20260708);

  TFile *fi = TFile::Open(in);
  TTree *t = (TTree *) fi->Get("ntp_g4hit");
  float event, glayer, gx, gy, gz, gt, gpl, gpx, gpy, gpz, gedep, gtrackID;
  t->SetBranchStatus("*", 0);
  for (const char *b : {"event", "glayer", "gx", "gy", "gz", "gt", "gpl", "gpx", "gpy", "gpz", "gedep", "gtrackID"})
  {
    t->SetBranchStatus(b, 1);
  }
  t->SetBranchAddress("event", &event);
  t->SetBranchAddress("glayer", &glayer);
  t->SetBranchAddress("gx", &gx);
  t->SetBranchAddress("gy", &gy);
  t->SetBranchAddress("gz", &gz);
  t->SetBranchAddress("gt", &gt);
  t->SetBranchAddress("gpl", &gpl);
  t->SetBranchAddress("gpx", &gpx);
  t->SetBranchAddress("gpy", &gpy);
  t->SetBranchAddress("gpz", &gpz);
  t->SetBranchAddress("gedep", &gedep);
  t->SetBranchAddress("gtrackID", &gtrackID);

  TFile *fo = new TFile(rawout, "RECREATE");
  TNtuple *o = new TNtuple("raw_pix", "raw pixel charge (pre-electronics)",
                           "event:layer:side:pad:tbin:q:trk");

  std::unordered_map<uint64_t, Cell> pix;
  std::vector<int> pads;
  std::vector<double> pshare;
  long npix = 0, ne = 0;
  int curev = -1;

  auto flush = [&](int ev) {
    for (auto &kv : pix)
    {
      uint64_t k = kv.first;
      float row[7] = {(float) ev, (float) (k >> 40U), (float) ((k >> 32U) & 0xFF),
                      (float) ((k >> 16U) & 0xFFFF), (float) (k & 0xFFFF),
                      kv.second.q, (float) kv.second.trk};
      o->Fill(row);
      npix++;
    }
    pix.clear();
  };

  Long64_t N = t->GetEntries();
  for (Long64_t i = 0; i < N; ++i)
  {
    t->GetEntry(i);
    int ev = (int) event;
    if (ev >= NEV)
    {
      break;
    }
    if (ev != curev)
    {
      if (curev >= 0)
      {
        flush(curev);
        printf("  transport: event %d done (%ld raw pixels, %ld electrons)\n", curev, npix, ne);
      }
      curev = ev;
    }
    int L = (int) glayer;
    if (L < 7 || L > 54)
    {
      // P5 dumper writes glayer=99 (PHG4Hit carries no offline layer);
      // derive from hit radius via the geometry table (nearest ring, 0.7 cm gate)
      double rr = std::sqrt((double) gx * gx + (double) gy * gy);
      double dbest = 0.7;
      int Lb = -1;
      for (int q = 7; q <= 54; ++q)
      {
        double d = std::fabs(rr - GEO[q].radius);
        if (d < dbest) { dbest = d; Lb = q; }
      }
      L = Lb;
    }
    if (L < 7 || L > 54 || gedep <= 0)
    {
      continue;
    }
    int nel = rng.Poisson((double) gedep * CFG::EPG);
    if (nel <= 0)
    {
      continue;
    }
    ne += nel;
    double pl = std::isfinite(gpl) ? (double) gpl : 0.;  // evaluator writes gpl=NaN for TPC
    double pn = std::sqrt((double) gpx * gpx + (double) gpy * gpy + (double) gpz * gpz);
    double dx = pn > 0 ? gpx / pn : 0, dy = pn > 0 ? gpy / pn : 0, dz = pn > 0 ? gpz / pn : 0;
    for (int ie = 0; ie < nel; ++ie)
    {
      double u = rng.Uniform() - 0.5;
      double x = gx + dx * pl * u, y = gy + dy * pl * u, z = gz + dz * pl * u;
      if (!std::isfinite(x) || !std::isfinite(z) || std::fabs(z) > CFG::HALFZ)
      {
        continue;
      }
      int sd = (z >= 0) ? 1 : 0;
      double Ld = CFG::HALFZ - std::fabs(z);
      double sT = CFG::DIFF_T * std::sqrt(Ld);
      x += rng.Gaus(0., sT);
      y += rng.Gaus(0., sT);
      double tarr = gt + Ld / CFG::VDRIFT + rng.Gaus(0., CFG::DIFF_L * std::sqrt(Ld)) / CFG::VDRIFT;
      if (!(tarr >= 0 && tarr < CFG::NTBIN * CFG::CLOCK))
      {
        continue;
      }
      double phi = std::atan2(y, x);
      if (inGap(phi, region(L)))
      {
        continue;
      }
      double q = rng.Exp(CFG::GAIN);
      double sprf = (cloud_k > 0) ? std::sqrt(CLOUD_EFF * CLOUD_EFF + cloud_k * cloud_k * Ld) : CLOUD_EFF;
      zigzag(L, phi, sprf, pads, pshare);
      if (pads.empty())
      {
        continue;
      }
      int tb0 = (int) (tarr / CFG::CLOCK);
      int iu = (int) ((tarr - tb0 * CFG::CLOCK) / CFG::CLOCK * NLUT);
      iu = std::max(0, std::min(NLUT - 1, iu));
      for (size_t ip = 0; ip < pads.size(); ++ip)
      {
        for (int k = 0; k < 8; ++k)
        {
          int tb = tb0 + k;
          if (tb >= CFG::NTBIN)
          {
            break;
          }
          float dq = (float) (q * pshare[ip] * LUT[iu][k]);
          if (dq <= 0)
          {
            continue;
          }
          uint64_t key = ((uint64_t) L << 40U) | ((uint64_t) sd << 32U) |
                         ((uint64_t) (uint32_t) pads[ip] << 16U) | (uint64_t) tb;
          Cell &c = pix[key];
          c.q += dq;
          if (dq > c.qbest)
          {
            c.qbest = dq;
            c.trk = (int) gtrackID;
          }
        }
      }
    }
  }
  if (curev >= 0)
  {
    flush(curev);
  }
  printf("tpc_transport: %s -> %s : %ld raw pixels, %ld electrons (%d events)\n", in, rawout, npix, ne, NEV);
  fo->cd();
  o->Write();
  fo->Close();
}

// ---------------- STAGE B: electronics/readout ----------------
void tpc_readout(const char *rawin, const char *out,
                 double gaincal = 1.0, double thr_adu = 15.0,
                 int ret_pre = 0, int ret_post = 0, int seed = 4711,
                 const char *deadmap = "", double thr2_adu = -1.0, double p_keep = 1.0, double sigma_pad = 0.0, double p2 = 0.0,
                 double tail_frac = 0.0, double tail_tau = 4.0,
                 double tail2_frac = 0.0, double tail2_tau = 60.0, double tail2_floor = 8.0, double tail2_q0 = 0.0, int tail2_emitonly = 0, double tail2_ptrig = 1.0,
                 double thr_r1 = -1.0, double gain_r1 = -1.0, double gain_r3 = -1.0, double sigma_pad_r1 = -1.0, double sigma_ped = 0.0, double mini_scale = 0.0)
{
  loadGeo();
  if (!geo_ok)
  {
    return;
  }
  TRandom3 rng(seed);

  // optional dead-channel mask: TPC_DEADCHANNELMAP payload (Multiple tree, LOCAL pad
  // within sector). Key = (layer, side, sector, localpad); sector assumed = globalpad/nper.
  std::unordered_map<uint32_t, char> dead;
  if (deadmap && deadmap[0])
  {
    TFile *fd = TFile::Open(deadmap);
    TTree *td = fd ? (TTree *) fd->Get("Multiple") : nullptr;
    if (td)
    {
      // payload Ipad convention is FEE-relative (values exceed pads/sector); use the
      // unambiguous channel POSITION (Fx,Fy in mm) -> phi -> global pad via our geometry
      Int_t dl, dside;
      Float_t fx, fy;
      td->SetBranchAddress("Ilayer", &dl);
      td->SetBranchAddress("Iside", &dside);
      td->SetBranchAddress("Fx", &fx);
      td->SetBranchAddress("Fy", &fy);
      for (long k = 0; k < td->GetEntries(); ++k)
      {
        td->GetEntry(k);
        if (dl < 7 || dl > 54)
        {
          continue;
        }
        double phi = std::atan2((double) fy, (double) fx);
        int pad = (int) std::lround((phi - GEO[dl].phi0) / GEO[dl].slope - 0.0);
        if (pad < 0)
        {
          pad += GEO[dl].nbins;
        }
        if (pad >= GEO[dl].nbins)
        {
          pad -= GEO[dl].nbins;
        }
        dead[((uint32_t) dl << 20U) | ((uint32_t) dside << 16U) | (uint32_t) pad] = 1;
      }
      printf("tpc_readout: dead map %s -> %zu (layer,side,globalpad) channels\n", deadmap, dead.size());
      fd->Close();
    }
  }
  long nmasked = 0;
  TFile *fi = TFile::Open(rawin);
  TNtuple *r = (TNtuple *) fi->Get("raw_pix");
  float ev, L, sd, pad, tb, q, trk;
  r->SetBranchAddress("event", &ev);
  r->SetBranchAddress("layer", &L);
  r->SetBranchAddress("side", &sd);
  r->SetBranchAddress("pad", &pad);
  r->SetBranchAddress("tbin", &tb);
  r->SetBranchAddress("q", &q);
  r->SetBranchAddress("trk", &trk);

  struct Row
  {
    uint64_t col;  // (event,layer,side,pad)
    int tb;
    float q, trk;
  };
  std::vector<Row> v;
  v.reserve((size_t) r->GetEntries());
  for (Long64_t i = 0; i < r->GetEntries(); ++i)
  {
    r->GetEntry(i);
    uint64_t col = ((uint64_t) (uint32_t) ev << 40U) | ((uint64_t) (uint32_t) L << 32U) |
                   ((uint64_t) (((int) sd) ? 1 : 0) << 24U) | (uint64_t) (uint32_t) pad;
    v.push_back({col, (int) tb, q, trk});
  }
  std::sort(v.begin(), v.end(), [](const Row &a, const Row &b) {
    return a.col == b.col ? a.tb < b.tb : a.col < b.col;
  });

  TFile *fo = new TFile(out, "RECREATE");
  TNtuple *o = new TNtuple("ntp_hit", "digitized pixels",
                           "event:layer:phibin:zbin:tbin:adc:zelem:phi:z:gtrackID");
  long nkept = 0;
  size_t i = 0;
  std::vector<double> adu;
  std::vector<char> keep;
  while (i < v.size())
  {
    size_t j = i;
    while (j < v.size() && v[j].col == v[i].col)
    {
      ++j;
    }
    // electronics for this pad column
    adu.assign(j - i, 0.);
    // B4.3: REAL region-dependent deterministic ZS (CDB payload TpcADUThresholds10R1_20R23;
    // confirmed by run 79507 per-region cliffs: R1 keeps adc>=11, R2/R3 keep adc>=21).
    // thr_r1 >= 0 enables it: thr_adu is then the R2/R3 payload value (20) and the whole
    // band/retention machinery (thr2/p_keep/ret_pre/ret_post) is DISABLED — the global
    // "band 11-19 + step at 20" of B3 was purely R1(thr10) stacked on R23(thr20).
    // p2 becomes a STANDALONE sub-threshold trace leak (real level ~1e-4).
    const int colL = (int) ((v[i].col >> 32U) & 0xFF);
    const bool regZS = (thr_r1 >= 0);
    const double thr_eff = regZS ? ((colL < 23 ? thr_r1 : thr_adu) + 1.0) : thr_adu;
    // B4.3 region gain factors (phenomenological; coarse stand-in for the CAEN
    // per-module corrections = B4.4): the real per-region pixel means demand
    // R1 ~1.2x, R3 ~1.06x relative to the R2-anchored global gaincal.
    double greg = 1.0;
    if (gain_r1 > 0 && colL < 23) greg = gain_r1;
    if (gain_r3 > 0 && colL >= 39) greg = gain_r3;
    const double spad_eff = (colL < 23 && sigma_pad_r1 > 0) ? sigma_pad_r1 : sigma_pad;
    // per-pad static pedestal spread (2026-07-12): real saturated maxadc is a BUMP
    // (931-963, sigma ~5: saturated adc = 1023 - ped_pad) where sim had a delta at
    // 949. The offset cancels in subtracted values everywhere EXCEPT at the clamp,
    // so ONLY the saturation ceiling spreads. Same static hash family as the gain.
    double dped = 0.0;
    if (sigma_ped > 0)
    {
      UInt_t pseed = ((UInt_t) (v[i].col & 0xFFFFFFFFFFULL)) ^ 0x5DEECE66U;
      if (pseed == 0) pseed = 0xBADC0FFEU;  // TRandom3(0) is UUID-random
      TRandom3 pr(pseed);
      dped = pr.Gaus(0., sigma_ped);
    }
    double gpad = 1.0;
    if (spad_eff > 0)
    {
      // NB seed hash: (UInt_t) truncation makes the gain static per (side,pad) shared
      // across layers (documented quirk, frozen for v3.2 reproducibility), and seed 0
      // (side 0, pad 0) would be UUID-seeded = NONDETERMINISTIC in TRandom3 - remap it.
      UInt_t gseed = (UInt_t) (v[i].col & 0xFFFFFFFFFFULL);
      if (gseed == 0) gseed = 0x9E3779B9;
      TRandom3 gr(gseed);
      gpad = std::exp(gr.Gaus(0., spad_eff) - 0.5 * spad_eff * spad_eff);
    }
    for (size_t k = i; k < j; ++k)
    {
      double a = v[k].q * gpad * greg * gaincal * CFG::MV_PER_E * CFG::ADU_PER_MV + (CFG::PEDESTAL + dped) + rng.Gaus(0., CFG::NOISE_ADU);
      if (a > 1023)
      {
        a = 1023;  // 10-bit saturation BEFORE pedestal subtraction
      }
      // real ADC is an integer ADU count; quantize like the hardware, then subtract pedestal
      adu[k - i] = std::round(a - (CFG::PEDESTAL + dped));
    }
    // ZS + SAMPA-DSP retention: keep above-threshold samples plus ret_pre/ret_post
    // CONSECUTIVE-tbin neighbours of any above-threshold sample
    const double t2 = (thr2_adu >= 0) ? thr2_adu : thr_adu;  // two-tier ZS: neighbours kept if >= thr2
    // B4: SAMPA/ion tail — big pulses induce a decaying tail that re-crosses threshold,
    // elongating high-charge clusters in TIME (real zsize grows with adc; sim did not).
    // B4.2 adds a SLOW second exponential (ion tail proper): recursive one-pole
    // accumulator per component, the same architecture ALICE's GEM-TPC online ion-tail
    // filter corrects for (arXiv:2304.03881: exponential tail, ~0.7% of peak but ~9% of
    // signal integral; recursive Q_corr weighted by k2=e^-slope). The stock sPHENIX sim
    // truncates the shaper at 400 ns (PHG4TpcPadPlaneReadout, 8 clocks) so it has NO
    // long-time response at all — this block is the detached-stage stand-in, amplitudes
    // tuned to run 79507's per-pad run-length tail (P(run>=20), P(run>=30)).
    // The tail rides on EVERY clock sample, not just tbins where transported charge
    // exists — so with the slow component on, synthetic samples are EMITTED into empty
    // tbins while the state is still retainable (>= EMIT_FLOOR, just below T2=11).
    // Synthetic samples get baseline noise, carry the source pixel's track id (the tail
    // belongs to its source cluster for truth purposes), respect the ADC ceiling and the
    // 971-tbin frame, and do NOT feed the state back (tail-of-tail is second order; this
    // also guarantees stability for any tail2_frac). With tail2_frac=0 this reduces
    // bit-for-bit to the B4.1 behavior.
    const double CEIL_ADU = 1023.0 - (CFG::PEDESTAL + dped);  // per-pad ceiling (sigma_ped)
    struct Px
    {
      int tb;
      double a;
      float trk;
    };
    std::vector<Px> col;
    col.reserve(adu.size() + 64);
    if (tail_frac > 0 || tail2_frac > 0)
    {
      // emission floor: real late tails sit JUST above T1 (measured run 79507 long-run
      // profile: 22-44 ADU plateau); emitting below T1 inflates medium runs through the
      // retention band (scan round 1), so the physical floor is ~T1.
      const double EMIT_FLOOR = regZS ? thr_eff : tail2_floor;  // chains die at the region threshold
      // SLOW COMPONENT (B4.2, scan round 4+): FIXED-DURATION LINEAR ramp, per ALICE's
      // measured GEM slow ion component ("nearly linear — ions uniformly produced in the
      // induction gap", arXiv:2304.03881): each source pixel with a > tail2_q0 launches
      // a tail of amplitude tail2_frac*(a-q0) that ramps linearly to zero over tail2_tau
      // (reinterpreted: DURATION in tbins). Fixed duration decouples trigger RATE from
      // run LENGTH — an exponential accumulator provably cannot fit both (rounds 1-3:
      // P10/P20 x2 overshoot at any setting reaching P30).
      const double T2LEN = tail2_tau;
      double s1 = 0.;
      struct Src
      {
        int tb;
        double amp, T;
      };
      std::vector<Src> src;  // active slow sources
      auto slow_at = [&](int t) {
        double s = 0.;
        size_t w = 0;
        for (size_t m = 0; m < src.size(); ++m)
        {
          double dt = t - src[m].tb;
          if (dt > src[m].T)
          {
            continue;  // expired source: compact away
          }
          src[w++] = src[m];
          if (dt > 0)
          {
            s += src[m].amp * (1.0 - dt / src[m].T);
          }
        }
        src.resize(w);
        return s;
      };
      int prevtb = -1000000;
      float lasttrk = -1.f;
      for (size_t k = 0; k < adu.size(); ++k)
      {
        int tb = v[i + k].tb;
        if (prevtb > -1000000)
        {
          int at = prevtb;  // time the fast state is currently valid at
          for (int gt = prevtb + 1; gt < tb && tail2_frac > 0; ++gt)
          {
            s1 *= std::exp(-1.0 / tail_tau);
            at = gt;
            double val = s1 + slow_at(gt);
            if (val < EMIT_FLOOR || gt > 970)
            {
              break;  // slow part only decreases inside a gap -> safe to stop
            }
            double a = std::round(val + rng.Gaus(0., CFG::NOISE_ADU));
            if (a > CEIL_ADU)
            {
              a = CEIL_ADU;
            }
            if (a >= 1)
            {
              col.push_back({gt, a, lasttrk});
            }
          }
          if (at < tb)
          {
            s1 *= std::exp(-(double) (tb - at) / tail_tau);
          }
        }
        // emit-only mode: the slow tail materializes ONLY as synthetic samples in
        // empty tbins (run-out chains) and never boosts real charged samples — scans
        // 1-4 showed any additive slow tail GLUES the dense sub-threshold raw
        // occupancy around medium clusters into 10-20 tbin runs real data lacks.
        double a = adu[k] + s1 + (tail2_emitonly ? 0. : slow_at(tb));
        // ADC saturates AFTER the induced tail adds at the input: re-clamp to the
        // hardware ceiling (1023 raw = 1023-PEDESTAL post-subtraction) — without this,
        // tail-augmented pixels exceeded the physical maximum (caught by user zoom).
        if (a > CEIL_ADU)
        {
          a = CEIL_ADU;
        }
        s1 += a * tail_frac;
        if (tail2_frac > 0)
        {
          if (tail2_emitonly == 2)
          {
            // mode 2: BINARY saturation-triggered disturbance — fixed amplitude
            // (tail2_frac, in ADU), retriggerable (one active disturbance per column).
            // Real run-length data: additions over baseline are ~equal at >=10/20/30,
            // i.e. tails either don't fire or run 20+ tbins — amplitude-proportional
            // models can't do that with a steeply falling charge spectrum (rounds 1-5).
            if (a >= tail2_q0 && (tail2_ptrig >= 1.0 || rng.Uniform() < tail2_ptrig))
            {
              // per-trigger dispersion (mean-preserving lognormal amp, gaussian
              // duration): fixed amp/duration produced a BIMODAL island zsize tail
              // (characteristic chain length ~ T*(1-T1/A)) where real is smooth.
              src.clear();
              double amp = tail2_frac * std::exp(rng.Gaus(0., 0.35) - 0.5 * 0.35 * 0.35);
              double dur = std::max(10., tail2_tau * (1.0 + rng.Gaus(0., 0.30)));
              src.push_back({tb, amp, dur});
            }
          }
          else if (a > tail2_q0)
          {
            src.push_back({tb, (a - tail2_q0) * tail2_frac, T2LEN});
          }
        }
        col.push_back({tb, a, v[i + k].trk});
        lasttrk = v[i + k].trk;
        prevtb = tb;
      }
      // tail run-out beyond the column's last charged sample
      if (tail2_frac > 0 && prevtb > -1000000)
      {
        for (int gt = prevtb + 1; gt <= 970; ++gt)
        {
          s1 *= std::exp(-1.0 / tail_tau);
          double val = s1 + slow_at(gt);
          if (val < EMIT_FLOOR)
          {
            break;
          }
          double a = std::round(val + rng.Gaus(0., CFG::NOISE_ADU));
          if (a > CEIL_ADU)
          {
            a = CEIL_ADU;
          }
          if (a >= 1)
          {
            col.push_back({gt, a, lasttrk});
          }
        }
      }
    }
    else
    {
      for (size_t k = 0; k < adu.size(); ++k)
      {
        col.push_back({v[i + k].tb, adu[k], v[i + k].trk});
      }
    }
    keep.assign(col.size(), 0);
    for (size_t k = 0; k < col.size(); ++k)
    {
      if (col[k].a >= thr_eff)
      {
        keep[k] = 1;
        if (regZS)
        {
          continue;  // deterministic region ZS: no neighbour retention, no band
        }
        for (int d = 1; d <= ret_pre; ++d)
        {
          if (k >= (size_t) d && col[k].tb - col[k - d].tb == d)
          {
            if (col[k - d].a >= t2 && (p_keep >= 1.0 || rng.Uniform() < p_keep))
            {
              keep[k - d] = 1;
            }
            else if (col[k - d].a >= 1 && p2 > 0 && rng.Uniform() < p2)
            {
              keep[k - d] = 1;  // sub-floor retention: real ZS leaks a flat trace (B3)
            }
          }
        }
        for (int d = 1; d <= ret_post; ++d)
        {
          if (k + d < col.size() && col[k + d].tb - col[k].tb == d)
          {
            if (col[k + d].a >= t2 && (p_keep >= 1.0 || rng.Uniform() < p_keep))
            {
              keep[k + d] = 1;
            }
            else if (col[k + d].a >= 1 && p2 > 0 && rng.Uniform() < p2)
            {
              keep[k + d] = 1;  // sub-floor retention: real ZS leaks a flat trace (B3)
            }
          }
        }
      }
      else if (regZS && col[k].a >= 1 && p2 > 0)
      {
        // SHAPED sub-threshold trace (2026-07-12): keep-probability per ADU count,
        // fitted as (real per-region trace)/(no-ZS raw spectrum) so the KEPT trace
        // reproduces run 79507 in SHAPE, not just rate (flat p2 inherited the steep
        // raw falloff: x2.5 hot at adc 1, x3 cold at 9-10). p2 = SCALE (nominal 1.0).
        static const double P2R1[11] = {0, 1.50e-05, 2.66e-05, 3.69e-05, 4.17e-05, 7.73e-05, 1.13e-04, 1.06e-04, 1.38e-04, 1.93e-04, 1.56e-04};
        static const double P2R2[21] = {0, 2.84e-05, 2.96e-05, 6.14e-05, 9.70e-05, 1.29e-04, 1.64e-04, 2.46e-04, 2.51e-04, 3.07e-04, 3.00e-04, 2.79e-04, 3.47e-04, 3.97e-04, 3.39e-04, 3.34e-04, 4.04e-04, 4.34e-04, 5.04e-04, 4.83e-04, 4.78e-04};
        static const double P2R3[21] = {0, 7.79e-05, 9.63e-05, 1.66e-04, 2.15e-04, 2.84e-04, 2.64e-04, 3.11e-04, 3.55e-04, 4.06e-04, 4.13e-04, 3.40e-04, 4.30e-04, 3.66e-04, 3.75e-04, 3.38e-04, 4.62e-04, 3.24e-04, 3.04e-04, 3.31e-04, 4.70e-04};
        int ai = (int) std::lround(col[k].a);
        double pk = 0;
        if (colL < 23 && ai >= 1 && ai <= 10) pk = P2R1[ai];
        else if (colL >= 23 && ai >= 1 && ai <= 20) pk = (colL < 39 ? P2R2[ai] : P2R3[ai]);
        if (pk > 0 && rng.Uniform() < pk * p2)
        {
          keep[k] = 1;
        }
      }
    }
    {
      int LL = (int) ((v[i].col >> 32U) & 0xFF);
      int ss = (int) ((v[i].col >> 24U) & 0x1);
      int pp = (int) (v[i].col & 0xFFFFFF);
      float evf = (float) (v[i].col >> 40U);
      bool masked = false;
      if (!dead.empty())
      {
        uint32_t dkey = ((uint32_t) LL << 20U) | ((uint32_t) ss << 16U) | (uint32_t) pp;
        masked = dead.count(dkey) > 0;
      }
      double phi = GEO[LL].phi0 + GEO[LL].slope * pp;
      for (size_t k = 0; k < col.size(); ++k)
      {
        if (!keep[k] || col[k].a <= 0)
        {
          continue;
        }
        if (masked)
        {
          nmasked++;
          continue;
        }
        double zapp = (ss == 1 ? 1 : -1) * (CFG::HALFZ - col[k].tb * CFG::CLOCK * CFG::VDRIFT);
        float row[10] = {evf, (float) LL, (float) pp, (float) col[k].tb, (float) col[k].tb,
                         (float) col[k].a, (float) ss, (float) phi, (float) zapp, col[k].trk};
        o->Fill(row);
        nkept++;
      }
    }
    i = j;
  }
  // mini-cluster background injection (2026-07-12): real low-adc (<80) 2px and
  // 3-4px clusters fill the gap between the sub-threshold trace spike and the main
  // population; deficits per frame 2px {178,55,101}, 3-4px {141,20,22}. Injected at
  // the OUTPUT stage (post-gain: exact adu spectra, yield ~1 — a composer-level
  // attempt died to the per-pad lognormal gain knocking near-threshold pixels under
  // the cut). Pixel adu ~ thr+1+Exp(7) (real R1 2px sums peak 24-36). trk sentinel
  // -7 stays unmatched in the truth sidecar -> islandize91 labels cls=2.
  if (mini_scale > 0)
  {
    static const double MTHR[3] = {11., 21., 21.};
    static const int NPADS_R[3] = {1128, 1536, 2304};
    static const double MR2[3] = {209., 65., 119.};   // deficits x1.17 merge-loss comp (pilot 6)
    static const double MR34[3] = {166., 30., 32.};  // R2/R3: 3px-only (4px spills adc<80)
    static const double MTAU34[3] = {5.0, 2.5, 2.5};  // 3-4px pixel Exp tail: sums must stay
                                                      // inside adc<80 (R2/R3: 3x(21+Exp(7))>80
                                                      // killed the class entirely, pilot 6)
    int evlo = (int) (v.empty() ? 0 : (v.front().col >> 40U));
    int evhi = (int) (v.empty() ? 0 : (v.back().col >> 40U));
    long nmini = 0;
    for (int ev = evlo; ev <= evhi; ++ev)
    {
      for (int rg = 0; rg < 3; ++rg)
      {
        for (int cls2 = 0; cls2 < 2; ++cls2)
        {
          int n = rng.Poisson((cls2 ? MR34[rg] : MR2[rg]) * mini_scale);
          for (int im = 0; im < n; ++im)
          {
            int lay = 7 + rg * 16 + (int) rng.Integer(16);
            int side = (int) rng.Integer(2);
            int pad = (int) rng.Integer(NPADS_R[rg] - 1);
            int tb0 = (int) rng.Integer(969);
            int npx = cls2 ? (rg == 0 ? 3 + (int) rng.Integer(2) : 3) : 2;
            bool horiz = rng.Uniform() < 0.5;
            if (!dead.empty())
            {
              uint32_t dk1 = ((uint32_t) lay << 20U) | ((uint32_t) side << 16U) | (uint32_t) pad;
              uint32_t dk2 = ((uint32_t) lay << 20U) | ((uint32_t) side << 16U) | (uint32_t) (pad + 1);
              if (dead.count(dk1) || dead.count(dk2))
              {
                continue;
              }
            }
            double phi = GEO[lay].phi0 + GEO[lay].slope * pad;
            for (int p = 0; p < npx; ++p)
            {
              int dp = (npx <= 2) ? (horiz ? p : 0) : (p & 1);
              int dt = (npx <= 2) ? (horiz ? 0 : p) : (p >> 1);
              double a = MTHR[rg] + 1.0 + rng.Exp(npx <= 2 ? 7.0 : MTAU34[rg]);
              float row[10] = {(float) ev, (float) lay, (float) (pad + dp), (float) (tb0 + dt), (float) (tb0 + dt),
                               (float) std::round(a), (float) side,
                               (float) (phi + GEO[lay].slope * dp),
                               (float) ((side == 1 ? 1 : -1) * (CFG::HALFZ - (tb0 + dt) * CFG::CLOCK * CFG::VDRIFT)),
                               -7.f};
              o->Fill(row);
              nkept++;
              nmini++;
            }
          }
        }
      }
    }
    printf("tpc_readout: mini background: %ld pixels injected\n", nmini);
  }
  printf("tpc_readout: %s -> %s : %ld pixels kept, %ld dead-masked (gaincal=%.3f thr=%.1f ret=%d/%d p=%.2f)\n",
         rawin, out, nkept, nmasked, gaincal, thr_adu, ret_pre, ret_post, p_keep);
  fo->cd();
  o->Write();
  fo->Close();
}

// ---------------- wrapper (day-1 compatible defaults) ----------------
void tpc_digitize(const char *in = "/home/rog/sPHENIX/3D_ClusterFindingML/macros-offline/detectors/sPHENIX/exam5_g4svtx_eval.root",
                  const char *out = "digi_sim.root", int NEV = 2)
{
  tpc_transport(in, "raw_pix.root", NEV);
  tpc_readout("raw_pix.root", out, 1.0, 15.0, 0, 0);
}
