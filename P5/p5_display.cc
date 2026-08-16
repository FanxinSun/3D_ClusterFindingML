// p5_display.cc — P5 interactive Geant4 event display for the v5.5 pp chain.
//
// Replaces the CVMFS/singularity-era sphenix_display.sh path: that script drove
// a container Fun4All macro and knows nothing about the P5 standalone chain.
// This one runs on the SAME host toolchain the batch pipeline runs on (local
// Geant4 11.4.2 + local ROOT), loads the SAME sphenix_p5.gdml the batch g4 stage
// loads, and overlays the batch final products on top of it.
//
// What it shows (frame mode, the default):
//   pixels   island_post/digi_frames_production_<VER>.root  ntp_hit
//            (x,y from tpc_geom_table.txt radius + pad phi; z = APPARENT z from
//             drift time, exactly as the digitizer wrote it)
//   clusters island_post/island91_frames_production_<VER>.root  ntp_cluster
//            (x,y,z taken straight from the 91-branch export)
//   tracks   truth cluster chains: ntp_truth gtrackID groups the clusters of one
//            truth particle; the polyline is those cluster centroids in radius
//            order. This chain has NO reconstructed tracker — "track" here means
//            the truth trajectory the clusters were labeled from.
//
// What it shows (library mode, /p5/g4event):
//   g4hits   P5/PP_g4hit_<i>.root ntp_g4hit for ONE library collision, plus the
//            per-G4-track polylines. Library events are the composer's INPUT;
//            they do not correspond 1:1 to production frames (a frame is a
//            pileup draw over many library collisions and the composer records
//            only the remapped newid, not the source event). So the two modes
//            are mutually exclusive by construction, not by preference.
//
// Geant4 itself is live: /run/beamOn works in the GUI with the pipeline's own
// HepMC input (--hepmc) or with the built-in gun, so you can watch G4 tracking
// on the same geometry that produced the batch g4hits.
//
// Build: ./build_display.sh     Drive: ../p5_display.sh

#include <G4Circle.hh>
#include <G4Colour.hh>
#include <G4Event.hh>
#include <G4FieldManager.hh>
#include <G4GDMLParser.hh>
#include <G4IonTable.hh>
#include <G4LogicalVolume.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4MagneticField.hh>
#include <G4NistManager.hh>
#include <G4ParticleGun.hh>
#include <G4ParticleTable.hh>
#include <G4Polyline.hh>
#include <G4Polymarker.hh>
#include <G4PrimaryParticle.hh>
#include <G4PrimaryVertex.hh>
#include <G4RunManager.hh>
#include <G4Square.hh>
#include <G4StepLimiterPhysics.hh>
#include <G4SystemOfUnits.hh>
#include <G4Text.hh>
#include <G4TransportationManager.hh>
#include <G4UIExecutive.hh>
#include <G4UIcmdWithABool.hh>
#include <G4UIcmdWithADouble.hh>
#include <G4UIcmdWithAString.hh>
#include <G4UIcmdWithAnInteger.hh>
#include <G4UIcmdWithoutParameter.hh>
#include <G4UIcommand.hh>
#include <G4UIdirectory.hh>
#include <G4UImanager.hh>
#include <G4UImessenger.hh>
#include <G4UIparameter.hh>
#include <G4UserLimits.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4VUserVisAction.hh>
#include <G4VVisManager.hh>
#include <G4VisAttributes.hh>
#include <G4VisExecutive.hh>
#include <G4VisExtent.hh>
#include <G4Scene.hh>
#include <G4VisManager.hh>
#include <FTFP_BERT.hh>

#include <TFile.h>
#include <TTree.h>

#include <sys/stat.h>
#include <sys/types.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// =====================================================================
//  palette
// =====================================================================
// Three data layers are on screen at once, so the categorical assignment is
// held to the three slots that clear the all-pairs CVD/normal-vision gates on a
// dark surface (validated: worst pair CVD dE 9.4, normal-vision 20.9 against
// #1a1a19). A fourth series-coloured layer does not exist here: the g4hit
// overlay is a separate MODE that re-uses the same three slots, not a fourth
// colour. Anything past the budget renders in chrome grey ("Other").
namespace PAL
{
// dark-surface categorical slots 1-3
const G4Colour BLUE(0.224, 0.529, 0.898);    // #3987e5  slot 1  pixels (hue)
const G4Colour ORANGE(0.851, 0.349, 0.149);  // #d95926  slot 2  clusters
const G4Colour AQUA(0.098, 0.620, 0.439);    // #199e70  slot 3  tracks
// chrome: not series colours
const G4Colour MUTED(0.537, 0.529, 0.506);   // #898781  "Other" / noise / labels
const G4Colour OTHER(0.290, 0.286, 0.275);   // #4a4a46  background marks in spotlight mode
const G4Colour GRID(0.173, 0.173, 0.165);    // #2c2c2a  dimmed data, hairlines
const G4Colour AXIS(0.220, 0.220, 0.208);    // #383835  geometry chrome
const G4Colour INK(1.000, 1.000, 1.000);     // #ffffff  primary ink (legend text)
const G4Colour SURFACE(0.102, 0.102, 0.098); // #1a1a19  viewer background

// Blue sequential ramp for magnitude, dark -> light because the surface is dark.
// It stops at step 400 (#3987e5) ON PURPOSE: the pixel layer is the BACKGROUND,
// and the full 100..700 scale ran to #cde2fb at luminance 0.87 — brighter than
// the cluster accent (0.478) and the chain accent (0.443), so high-ADC pixels
// owned the brightest marks on screen and the hierarchy read inverted. Capped
// here, every stop clears both accents by 20.9-34.3 normal-vision dE (floor 15)
// and none of them outshines an accent.
const int NRAMP = 5;
const G4Colour RAMP[NRAMP] = {
    G4Colour(0.094, 0.310, 0.584),  // #184f95  600
    G4Colour(0.110, 0.361, 0.671),  // #1c5cab  550
    G4Colour(0.145, 0.416, 0.749),  // #256abf  500
    G4Colour(0.165, 0.471, 0.839),  // #2a78d6  450
    G4Colour(0.224, 0.529, 0.898)}; // #3987e5  400

// the 8 documented dark categorical slots, in the documented order, used ONLY
// by the track-spotlight mode. Slots past 3 do not clear the all-pairs floors —
// /p5/trk/top says so out loud when you ask for more.
const int NSLOT = 8;
const G4Colour SLOT[NSLOT] = {
    G4Colour(0.224, 0.529, 0.898),  // #3987e5 blue
    G4Colour(0.851, 0.349, 0.149),  // #d95926 orange
    G4Colour(0.098, 0.620, 0.439),  // #199e70 aqua
    G4Colour(0.788, 0.522, 0.000),  // #c98500 yellow
    G4Colour(0.835, 0.318, 0.506),  // #d55181 magenta
    G4Colour(0.000, 0.514, 0.000),  // #008300 green
    G4Colour(0.565, 0.522, 0.914),  // #9085e9 violet
    G4Colour(0.902, 0.404, 0.404)}; // #e66767 red

G4Colour ramp(double t)
{
  if (!(t > 0)) t = 0;
  if (t > 1) t = 1;
  double f = t * (NRAMP - 1);
  int i = (int) f;
  if (i >= NRAMP - 1) return RAMP[NRAMP - 1];
  double u = f - i;
  return G4Colour(RAMP[i].GetRed() + u * (RAMP[i + 1].GetRed() - RAMP[i].GetRed()),
                  RAMP[i].GetGreen() + u * (RAMP[i + 1].GetGreen() - RAMP[i].GetGreen()),
                  RAMP[i].GetBlue() + u * (RAMP[i + 1].GetBlue() - RAMP[i].GetBlue()));
}
}  // namespace PAL

// =====================================================================
//  configuration
// =====================================================================
struct Cfg
{
  std::string repo = "/home/rog/sPHENIX/3D_ClusterFindingML";
  std::string ver = "v55";
  std::string gdml, digi, isl, geomtab, rowdr, fieldmap, hepmc, g4pat, g4out;
  std::string session = "qt", visdriver = "OGL", vissize = "1400x900", macro;
  // UI commands queued from the shell and replayed once the viewer exists —
  // anything that needs a live viewer (styles, views, /vis/...) goes here
  // rather than becoming a bespoke flag.
  std::vector<std::string> cmds;
  bool batch = false;

  int frame = 0;         // production frame (0..NB*PER-1)
  int g4chunk = 0;       // PP_g4hit_<chunk>.root
  int g4event = -1;      // first library event; <0 = frame mode
  std::vector<int> g4evList;  // every library event to superimpose
  int layerLo = 7, layerHi = 54;
  int sector = -1;       // TPC 30-degree wedge 0..11; -1 = all
  double zLo = -105.5, zHi = 105.5;  // physical TPC by default; see /p5/zrange
  double adcMin = 0, adcMax = 1023;
  long pixMax = 250000;  // stride-sample above this, and say so
  double pixScreen = 0;  // 0 = dots (1 px); >0 = filled circles of that size
  bool pixScreenSet = false;  // explicit user choice beats the per-mode default
  double clusScreen = 2.2;  // base screen size of a cluster mark
  int trkMinClus = 6;
  double trkPtMin = 0.2;  // above the 0.164 GeV looper wall by default
  bool trkPrimaryOnly = false;
  int trkTop = 3;        // spotlight slots in colour mode "track"
  int trkSelect = 0;     // spotlight one specific truth id
  bool showPix = true, showClus = true, showTrk = true, showG4 = true;
  bool showLegend = true;
  std::string colour = "level";   // level | class | track
  std::string geometry = "tpc";   // startup geometry preset
  bool applyRowDr = true;
} C;

// geometry chrome lives far below (it needs the G4 stores) but the data loaders
// and the overlay legend both reach into it
static void applyChromeColours();
static bool viewerLive();
static std::vector<int> parseEvSpec(const std::string &spec);
static std::string evSpecLabel();
static void setG4Events(const std::string &spec);

// =====================================================================
//  TPC row geometry (the pipeline's own tables)
// =====================================================================
struct Lay
{
  int nbins = 0;
  double radius = 0, slope = 0, phi0 = 0;
};
static Lay GEO[55];
static double DROW[55] = {0};
static bool GEO_OK = false;

// TPC sector convention, identical to islandize91: wrap phi into [0,2pi) and
// divide by 30 degrees.
static bool inSector(double phi)
{
  if (C.sector < 0) return true;
  double p = phi < 0 ? phi + 2 * M_PI : phi;
  int s = (int) (p / (M_PI / 6.));
  return std::min(11, std::max(0, s)) == C.sector;
}

static void loadGeom()
{
  FILE *g = fopen(C.geomtab.c_str(), "r");
  if (!g)
  {
    printf("p5_display: WARNING tpc_geom_table.txt not found at %s — pixel radii unavailable\n",
           C.geomtab.c_str());
    return;
  }
  char line[256];
  int n = 0;
  while (fgets(line, sizeof line, g))
  {
    int L, nb;
    double r, sl, p0, p1;
    if (line[0] == '#') continue;
    if (sscanf(line, "%d %d %lf %lf %lf %lf", &L, &nb, &r, &sl, &p0, &p1) == 6 && L >= 0 && L < 55)
    {
      GEO[L] = {nb, r, sl, p0};
      n++;
    }
  }
  fclose(g);
  GEO_OK = n > 0;
  printf("p5_display: row geometry %d layers from %s\n", n, C.geomtab.c_str());

  if (!C.rowdr.empty())
  {
    FILE *f = fopen(C.rowdr.c_str(), "r");
    if (!f)
    {
      printf("p5_display: WARNING row-dr table %s missing — pixels drawn at GDML nominal radius\n",
             C.rowdr.c_str());
      return;
    }
    int m = 0;
    while (fgets(line, sizeof line, f))
    {
      int L;
      double dr;
      if (line[0] == '#') continue;
      if (sscanf(line, "%d %lf", &L, &dr) >= 2 && L >= 0 && L < 55)
      {
        DROW[L] = dr;
        m++;
      }
    }
    fclose(f);
    printf("p5_display: row-dr offsets %d layers from %s (matches the islandize91 export)\n",
           m, C.rowdr.c_str());
  }
}

static double layerR(int L)
{
  if (L < 0 || L >= 55) return 0;
  return GEO[L].radius + (C.applyRowDr ? DROW[L] : 0.0);
}

// =====================================================================
//  per-event entry index (cached)
// =====================================================================
// The production trees are ~10^8 rows and are NOT globally sorted by event
// (the digitizer appends its injected mini-cluster background after the main
// per-frame block), but they ARE piecewise contiguous: 2.00 runs/event in
// ntp_hit, 1.00 in ntp_cluster. So the whole index is a few hundred ranges —
// scan the event branch once, cache it, and every later frame switch is free.
struct Range
{
  Long64_t first, n;
};
using EvtMap = std::map<int, std::vector<Range>>;

static std::string cacheDir(const std::string &file)
{
  const char *env = getenv("P5_DISP_CACHE");
  if (env && *env) return env;
  size_t s = file.find_last_of('/');
  return (s == std::string::npos ? std::string(".") : file.substr(0, s)) + "/.p5disp";
}

static std::string cachePath(const std::string &file, const char *tree)
{
  std::string dir = cacheDir(file);
  ::mkdir(dir.c_str(), 0755);
  size_t s = file.find_last_of('/');
  std::string base = (s == std::string::npos) ? file : file.substr(s + 1);
  return dir + "/" + base + "." + tree + ".idx";
}

static bool idxLoad(const std::string &path, Long64_t nent, EvtMap &m)
{
  std::ifstream in(path);
  if (!in) return false;
  std::string tag;
  Long64_t got = -1;
  in >> tag >> got;
  if (tag != "p5idx1" || got != nent) return false;
  int nev = 0;
  in >> nev;
  for (int i = 0; i < nev; ++i)
  {
    int ev, nr;
    in >> ev >> nr;
    std::vector<Range> v(nr);
    for (int k = 0; k < nr; ++k) in >> v[k].first >> v[k].n;
    m[ev] = v;
  }
  return in.good() || in.eof();
}

static void idxSave(const std::string &path, Long64_t nent, const EvtMap &m)
{
  std::ofstream out(path);
  if (!out) return;
  out << "p5idx1 " << nent << "\n" << m.size() << "\n";
  for (auto &kv : m)
  {
    out << kv.first << " " << kv.second.size();
    for (auto &r : kv.second) out << " " << r.first << " " << r.n;
    out << "\n";
  }
}

// scan the 'event' branch and collapse into contiguous ranges
static EvtMap idxBuild(TTree *t)
{
  EvtMap m;
  t->SetBranchStatus("*", 0);
  t->SetBranchStatus("event", 1);
  float ev = 0;
  t->SetBranchAddress("event", &ev);
  Long64_t N = t->GetEntries(), start = 0;
  int cur = -2147483647;
  for (Long64_t i = 0; i < N; ++i)
  {
    t->GetEntry(i);
    int e = (int) ev;
    if (e != cur)
    {
      if (i > start) m[cur].push_back({start, i - start});
      cur = e;
      start = i;
    }
  }
  if (N > start) m[cur].push_back({start, N - start});
  t->SetBranchStatus("*", 1);
  t->ResetBranchAddresses();
  return m;
}

static const EvtMap &getIdx(const std::string &file, const char *tree, TTree *t)
{
  static std::map<std::string, EvtMap> cache;
  std::string key = file + "#" + tree;
  auto it = cache.find(key);
  if (it != cache.end()) return it->second;
  EvtMap m;
  std::string p = cachePath(file, tree);
  if (!idxLoad(p, t->GetEntries(), m))
  {
    printf("p5_display: indexing %s:%s (%lld rows, one time) ...\n", file.c_str(), tree,
           t->GetEntries());
    fflush(stdout);
    m = idxBuild(t);
    idxSave(p, t->GetEntries(), m);
    printf("p5_display: index cached -> %s (%zu events)\n", p.c_str(), m.size());
  }
  return cache.emplace(key, std::move(m)).first->second;
}

// =====================================================================
//  data model
// =====================================================================
struct Pix
{
  float x, y, z, adc;
  int layer;
};
struct Clus
{
  float x, y, z, r, adc, size, phisize, zsize, gpt, purity;
  int layer, trk, cls, prim;
};
struct G4H
{
  float x, y, z, edep, t;
  int ev, trk, prim;
};
struct Trk
{
  int id = 0, prim = 0;
  float pt = -1;
  std::vector<const Clus *> cl;
};

static std::vector<Pix> g_pix;
static std::vector<Clus> g_clus;
static std::vector<G4H> g_g4;
static std::vector<Trk> g_trk;
static std::map<long long, std::vector<G4H *>> g_g4trk;
static long g_pixTotal = 0, g_pixDrop = 0, g_clusTotal = 0;
static long g_pixOutZ = 0, g_clusOutZ = 0;
static int g_pixStride = 1;
static std::string g_status = "no data loaded";

// ---------------------------------------------------------------- pixels
static void loadPixels()
{
  g_pix.clear();
  g_pixTotal = g_pixDrop = g_pixOutZ = 0;
  g_pixStride = 1;
  if (C.digi.empty()) return;
  TFile *f = TFile::Open(C.digi.c_str());
  if (!f || f->IsZombie())
  {
    printf("p5_display: cannot open %s\n", C.digi.c_str());
    return;
  }
  TTree *t = (TTree *) f->Get("ntp_hit");
  if (!t)
  {
    printf("p5_display: no ntp_hit in %s\n", C.digi.c_str());
    f->Close();
    return;
  }
  if (!GEO_OK)
    printf("p5_display: WARNING no row geometry — pixel radii are 0, nothing will be placed\n");
  const EvtMap &idx = getIdx(C.digi, "ntp_hit", t);
  auto it = idx.find(C.frame);
  if (it == idx.end())
  {
    printf("p5_display: frame %d not in %s\n", C.frame, C.digi.c_str());
    f->Close();
    return;
  }
  float lay, adc, phi, z;
  t->SetBranchStatus("*", 0);
  for (const char *b : {"layer", "adc", "phi", "z"}) t->SetBranchStatus(b, 1);
  t->SetBranchAddress("layer", &lay);
  t->SetBranchAddress("adc", &adc);
  t->SetBranchAddress("phi", &phi);
  t->SetBranchAddress("z", &z);

  // first pass over the frame's ranges: count survivors, then stride to the cap
  for (int pass = 0; pass < 2; ++pass)
  {
    long kept = 0;
    for (auto &r : it->second)
      for (Long64_t i = r.first; i < r.first + r.n; ++i)
      {
        t->GetEntry(i);
        int L = (int) lay;
        if (pass == 0) g_pixTotal++;
        if (adc < C.adcMin) { if (pass == 0) g_pixDrop++; continue; }
        if (z < C.zLo || z > C.zHi) { if (pass == 0) g_pixOutZ++; continue; }
        if (L < C.layerLo || L > C.layerHi) { if (pass == 0) g_pixDrop++; continue; }
        if (!inSector(phi)) { if (pass == 0) g_pixDrop++; continue; }
        if (pass == 0) { kept++; continue; }
        if (g_pixStride > 1 && (kept++ % g_pixStride)) continue;
        double R = layerR(L);
        g_pix.push_back({(float) (R * std::cos(phi)), (float) (R * std::sin(phi)), z, adc, L});
      }
    if (pass == 0)
    {
      g_pixStride = (C.pixMax > 0 && kept > C.pixMax) ? (int) ((kept + C.pixMax - 1) / C.pixMax) : 1;
      g_pix.reserve(kept / g_pixStride + 8);
      if (g_pixStride > 1)
        printf("p5_display: %ld pixels pass the cuts, drawing every %dth (cap /p5/pix/max %ld)\n",
               kept, g_pixStride, C.pixMax);
    }
  }
  t->ResetBranchAddresses();
  f->Close();
}

// ------------------------------------------------------- clusters + truth
static void loadClusters()
{
  g_clus.clear();
  g_trk.clear();
  g_clusTotal = g_clusOutZ = 0;
  if (C.isl.empty()) return;
  TFile *f = TFile::Open(C.isl.c_str());
  if (!f || f->IsZombie())
  {
    printf("p5_display: cannot open %s\n", C.isl.c_str());
    return;
  }
  TTree *c = (TTree *) f->Get("ntp_cluster");
  TTree *u = (TTree *) f->Get("ntp_truth");
  if (!c)
  {
    printf("p5_display: no ntp_cluster in %s\n", C.isl.c_str());
    f->Close();
    return;
  }
  const EvtMap &idx = getIdx(C.isl, "ntp_cluster", c);
  auto it = idx.find(C.frame);
  if (it == idx.end())
  {
    printf("p5_display: frame %d not in %s\n", C.frame, C.isl.c_str());
    f->Close();
    return;
  }
  float x, y, z, r, adc, lay, sz, psz, zsz;
  c->SetBranchStatus("*", 0);
  for (const char *b : {"x", "y", "z", "r", "adc", "layer", "size", "phisize", "zsize"})
    c->SetBranchStatus(b, 1);
  c->SetBranchAddress("x", &x);
  c->SetBranchAddress("y", &y);
  c->SetBranchAddress("z", &z);
  c->SetBranchAddress("r", &r);
  c->SetBranchAddress("adc", &adc);
  c->SetBranchAddress("layer", &lay);
  c->SetBranchAddress("size", &sz);
  c->SetBranchAddress("phisize", &psz);
  c->SetBranchAddress("zsize", &zsz);

  float gid = 0, gpt = -1, gpr = 0, cls = -1, pur = 0;
  if (u)
  {
    u->SetBranchStatus("*", 0);
    for (const char *b : {"gtrackID", "gpt", "gprimary", "cls", "purity"}) u->SetBranchStatus(b, 1);
    u->SetBranchAddress("gtrackID", &gid);
    u->SetBranchAddress("gpt", &gpt);
    u->SetBranchAddress("gprimary", &gpr);
    u->SetBranchAddress("cls", &cls);
    u->SetBranchAddress("purity", &pur);
  }
  for (auto &rg : it->second)
    for (Long64_t i = rg.first; i < rg.first + rg.n; ++i)
    {
      c->GetEntry(i);
      g_clusTotal++;
      int L = (int) lay;
      if (z < C.zLo || z > C.zHi) { g_clusOutZ++; continue; }
      if (L < C.layerLo || L > C.layerHi) continue;
      if (!inSector(std::atan2((double) y, (double) x))) continue;
      // ntp_truth is written row-for-row with ntp_cluster — same entry index
      if (u) u->GetEntry(i);
      g_clus.push_back({x, y, z, r, adc, sz, psz, zsz, u ? gpt : -1.f, u ? pur : 0.f, L,
                        u ? (int) gid : 0, u ? (int) cls : -1, u ? (int) gpr : 0});
    }
  c->ResetBranchAddresses();
  if (u) u->ResetBranchAddresses();
  f->Close();

  // ---- truth cluster chains = "tracks" -------------------------------
  std::unordered_map<int, std::vector<const Clus *>> by;
  for (auto &cl : g_clus)
    if (cl.trk != 0 && cl.trk != -9999) by[cl.trk].push_back(&cl);
  for (auto &kv : by)
  {
    if ((int) kv.second.size() < C.trkMinClus) continue;
    const Clus *s = kv.second.front();
    if (C.trkPrimaryOnly && s->prim <= 0) continue;
    if (s->gpt < C.trkPtMin) continue;
    Trk tk;
    tk.id = kv.first;
    tk.pt = s->gpt;
    tk.prim = s->prim;
    tk.cl = kv.second;
    std::sort(tk.cl.begin(), tk.cl.end(),
              [](const Clus *a, const Clus *b) { return a->r < b->r; });
    g_trk.push_back(std::move(tk));
  }
  std::sort(g_trk.begin(), g_trk.end(), [](const Trk &a, const Trk &b) { return a.pt > b.pt; });
}

// --------------------------------------------------------------- g4 hits
static void loadG4()
{
  g_g4.clear();
  g_g4trk.clear();
  if (C.g4event < 0 || C.g4pat.empty()) return;
  char buf[1024];
  snprintf(buf, sizeof buf, C.g4pat.c_str(), C.g4chunk);
  const std::string fn(buf);
  TFile *f = TFile::Open(fn.c_str());
  if (!f || f->IsZombie())
  {
    printf("p5_display: cannot open %s\n", fn.c_str());
    return;
  }
  TTree *t = (TTree *) f->Get("ntp_g4hit");
  if (!t)
  {
    printf("p5_display: no ntp_g4hit in %s\n", fn.c_str());
    f->Close();
    return;
  }
  const EvtMap &idx = getIdx(fn, "ntp_g4hit", t);
  float gx, gy, gz, ge, gt, gid, gpr;
  t->SetBranchStatus("*", 0);
  for (const char *b : {"gx", "gy", "gz", "gedep", "gt", "gtrackID", "gprimary"})
    t->SetBranchStatus(b, 1);
  t->SetBranchAddress("gx", &gx);
  t->SetBranchAddress("gy", &gy);
  t->SetBranchAddress("gz", &gz);
  t->SetBranchAddress("gedep", &ge);
  t->SetBranchAddress("gt", &gt);
  t->SetBranchAddress("gtrackID", &gid);
  t->SetBranchAddress("gprimary", &gpr);
  int nfound = 0, nempty = 0;
  for (int ev : C.g4evList)
  {
    auto it = idx.find(ev);
    if (it == idx.end()) { nempty++; continue; }
    nfound++;
    for (auto &r : it->second)
      for (Long64_t i = r.first; i < r.first + r.n; ++i)
      {
        t->GetEntry(i);
        g_g4.push_back({gx, gy, gz, ge, gt, ev, (int) gid, (int) gpr});
      }
  }
  t->ResetBranchAddresses();
  f->Close();
  if (nempty)
    printf("p5_display: %d of %zu requested events have no TPC hits\n", nempty, C.g4evList.size());
  // key on (event, trackID): each library event numbers its tracks from 1, so a
  // bare trackID would weld unrelated particles from different collisions into
  // one polyline
  for (auto &h : g_g4) g_g4trk[((long long) h.ev << 32) | (uint32_t) h.trk].push_back(&h);
  // order by STEP TIME, not radius: a G4 track is a trajectory, and low-energy
  // secondaries spiral inward and outward again — radius order would connect
  // those points into a scribble
  for (auto &kv : g_g4trk)
    std::sort(kv.second.begin(), kv.second.end(),
              [](const G4H *a, const G4H *b) { return a->t < b->t; });
  printf("p5_display: library %s event %s -> %zu g4 hits, %zu G4 tracks\n", fn.c_str(),
         evSpecLabel().c_str(), g_g4.size(), g_g4trk.size());
}

// "7" | "0-9" | "0,3,7-9" | "-1" (frame mode). A leading '-' is a sign, not a
// range separator, so the search for the dash starts at index 1.
static std::vector<int> parseEvSpec(const std::string &spec)
{
  std::vector<int> out;
  std::stringstream ss(spec);
  std::string tok;
  while (std::getline(ss, tok, ','))
  {
    while (!tok.empty() && isspace((unsigned char) tok.front())) tok.erase(tok.begin());
    while (!tok.empty() && isspace((unsigned char) tok.back())) tok.pop_back();
    if (tok.empty()) continue;
    size_t d = tok.find('-', 1);
    try
    {
      if (d == std::string::npos) out.push_back(std::stoi(tok));
      else
      {
        int a = std::stoi(tok.substr(0, d)), b = std::stoi(tok.substr(d + 1));
        if (b < a) std::swap(a, b);
        for (int i = a; i <= b; ++i) out.push_back(i);
      }
    }
    catch (...) { printf("p5_display: cannot parse event spec '%s'\n", tok.c_str()); }
  }
  return out;
}

// a runaway range (a typo like 0-19999) would pull millions of rows, so cap the
// event COUNT and say so rather than silently truncating
static void setG4Events(const std::string &spec)
{
  std::vector<int> v = parseEvSpec(spec);
  if (v.empty()) { printf("p5_display: /p5/g4event <n> | <lo>-<hi> | <list>; -1 = frame\n"); return; }
  if (v.size() == 1 && v[0] < 0) { C.g4event = -1; C.g4evList.clear(); return; }
  const size_t CAP = 500;
  if (v.size() > CAP)
  {
    printf("p5_display: %zu events requested, drawing the first %zu\n", v.size(), CAP);
    v.resize(CAP);
  }
  C.g4evList = v;
  C.g4event = v.front();
}

static std::string evSpecLabel()
{
  if (C.g4evList.size() <= 1) return std::to_string(C.g4event);
  bool contiguous = true;
  for (size_t i = 1; i < C.g4evList.size(); ++i)
    if (C.g4evList[i] != C.g4evList[i - 1] + 1) { contiguous = false; break; }
  char b[96];
  if (contiguous)
  {
    snprintf(b, sizeof b, "%d-%d (%zu)", C.g4evList.front(), C.g4evList.back(), C.g4evList.size());
    return b;
  }
  if (C.g4evList.size() <= 4)
  {
    std::string j;
    for (size_t i = 0; i < C.g4evList.size(); ++i)
      j += (i ? "," : "") + std::to_string(C.g4evList[i]);
    return j;
  }
  snprintf(b, sizeof b, "%zu picked", C.g4evList.size());
  return b;
}

// A minimum-bias pp collision is genuinely sparse — the median library event
// puts ~3k G4 steps in the TPC gas and the spread is an order of magnitude — so
// "which event is worth looking at" is a real question. The cached range index
// already knows the answer without touching the ROOT payload.
static void listG4Events(int n)
{
  char buf[1024];
  snprintf(buf, sizeof buf, C.g4pat.c_str(), C.g4chunk);
  const std::string fn(buf);
  TFile *f = TFile::Open(fn.c_str());
  if (!f || f->IsZombie()) { printf("p5_display: cannot open %s\n", fn.c_str()); return; }
  TTree *t = (TTree *) f->Get("ntp_g4hit");
  if (!t) { f->Close(); return; }
  const EvtMap &idx = getIdx(fn, "ntp_g4hit", t);
  std::vector<std::pair<long, int>> v;
  for (auto &kv : idx)
  {
    long s = 0;
    for (auto &r : kv.second) s += r.n;
    v.push_back({s, kv.first});
  }
  f->Close();
  if (v.empty()) return;
  std::sort(v.rbegin(), v.rend());
  long med = v[v.size() / 2].first;
  printf("  chunk %d: %zu events with TPC hits, median %ld g4 hits/event\n", C.g4chunk, v.size(),
         med);
  printf("  %-10s %-10s %s\n", "event", "g4 hits", "vs median");
  for (int i = 0; i < n && i < (int) v.size(); ++i)
    printf("  %-10d %-10ld x%.1f\n", v[i].second, v[i].first, (double) v[i].first / std::max(1L, med));
  printf("  /p5/g4event <n> to load one\n");
}

// -------------------------------------------------------------- summary
static void makeStatus()
{
  char b[512];
  if (C.g4event >= 0)
    snprintf(b, sizeof b, "chunk %d event%s %s   %zu g4 hits   %zu G4 tracks", C.g4chunk,
             C.g4evList.size() > 1 ? "s" : "", evSpecLabel().c_str(), g_g4.size(), g_g4trk.size());
  else
    snprintf(b, sizeof b, "frame %d   %zu/%ld px   %zu/%ld clus   %zu truth chains", C.frame,
             g_pix.size(), g_pixTotal, g_clus.size(), g_clusTotal, g_trk.size());
  g_status = b;
}

static void reload()
{
  if (C.g4event >= 0)
  {
    g_pix.clear();
    g_clus.clear();
    g_trk.clear();
    loadG4();
  }
  else
  {
    g_g4.clear();
    g_g4trk.clear();
    loadPixels();
    loadClusters();
  }
  makeStatus();
  printf("p5_display: %s\n", g_status.c_str());
  if (C.g4event < 0 && (g_pixOutZ || g_clusOutZ))
    printf("p5_display: outside the z window [%.1f,%.1f] cm: %ld px, %ld clus "
           "(apparent z runs to +-305 cm over the 965-tbin readout — /p5/zrange full)\n",
           C.zLo, C.zHi, g_pixOutZ, g_clusOutZ);
  if (C.g4event < 0 && g_pixDrop)
    printf("p5_display: %ld px cut by the row/sector/ADC selection\n", g_pixDrop);
}

// forward declarations: the geometry system table lives below (it needs the
// G4 stores), but the overlay's legend names the systems that are on screen
struct GeoSys;
extern std::vector<GeoSys> G_SYS;
extern std::vector<char> g_sysOn;
extern bool g_chromeTint;
const char *sysName(size_t i);
bool sysTinted(size_t i);
G4Colour sysColour(size_t i);
size_t sysCount();

// =====================================================================
//  the overlay
// =====================================================================
class Overlay : public G4VUserVisAction
{
 public:
  void Draw() override
  {
    G4VVisManager *vis = G4VVisManager::GetConcreteInstance();
    if (!vis) return;
    if (C.g4event >= 0) drawLibrary(vis);
    else drawFrame(vis);
    if (C.showLegend) drawLegend(vis);
  }

 private:
  static G4Point3D P(double x, double y, double z)
  {
    return G4Point3D(x * cm, y * cm, z * cm);
  }

  // pixels: one polymarker per ramp bucket (a polymarker carries one colour)
  void drawPix(G4VVisManager *vis, const std::vector<std::pair<G4Point3D, double>> &pts, bool dim,
               double lo, double hi) const
  {
    if (pts.empty()) return;
    const int NB = PAL::NRAMP;
    std::vector<G4Polymarker> bucket((size_t) NB);
    for (int i = 0; i < NB; ++i)
    {
      // a frame has ~200k pixels (dots are right); a library event has ~5k G4
      // steps, which as 1-px dots are nearly invisible — so unless the user said
      // otherwise, library mode draws them as small circles
      double sz = C.pixScreenSet ? C.pixScreen : (C.g4event >= 0 ? 2.5 : 0.0);
      bucket[i].SetMarkerType(sz > 0 ? G4Polymarker::circles : G4Polymarker::dots);
      bucket[i].SetFillStyle(G4VMarker::filled);
      if (sz > 0) bucket[i].SetScreenSize(sz);
      G4VisAttributes va(dim ? PAL::OTHER : PAL::ramp((double) i / (NB - 1)));
      bucket[i].SetVisAttributes(va);
    }
    const double denom = std::log1p(std::max(hi, lo + 1) - lo);
    for (auto &p : pts)
    {
      double t = std::log1p(std::max(0.0, p.second - lo)) / denom;
      int b = (int) (t * (NB - 1) + 0.5);
      bucket[std::min(NB - 1, std::max(0, b))].push_back(p.first);
    }
    for (auto &b : bucket)
      if (!b.empty()) vis->Draw(b);
  }

  void drawFrame(G4VVisManager *vis) const
  {
    const bool dimPix = (C.colour != "level");

    if (C.showPix && !g_pix.empty())
    {
      std::vector<std::pair<G4Point3D, double>> pts;
      pts.reserve(g_pix.size());
      for (auto &p : g_pix) pts.push_back({P(p.x, p.y, p.z), p.adc});
      drawPix(vis, pts, dimPix, C.adcMin, C.adcMax);
    }

    if (C.showClus && !g_clus.empty())
    {
      // spotlight table for colour mode "track": the N highest-pT truth chains
      std::unordered_map<int, int> slot;
      if (C.colour == "track")
      {
        if (C.trkSelect != 0) slot[C.trkSelect] = 0;
        else
          for (int i = 0; i < (int) g_trk.size() && i < C.trkTop; ++i) slot[g_trk[i].id] = i;
      }
      // one polymarker per (colour, size bucket)
      struct Key
      {
        int col, sz;
        bool operator<(const Key &o) const { return col != o.col ? col < o.col : sz < o.sz; }
      };
      std::map<Key, G4Polymarker> pm;
      for (auto &cl : g_clus)
      {
        int col;
        if (C.colour == "class") col = (cl.cls == 0) ? 0 : (cl.cls == 1 ? 1 : -1);
        else if (C.colour == "track")
        {
          auto s = slot.find(cl.trk);
          col = (s == slot.end()) ? -1 : 100 + s->second;
        }
        else col = 200;  // "level": one accent for the whole cluster level
        int szb = std::min(3, (int) (cl.size <= 2 ? 0 : cl.size <= 5 ? 1 : cl.size <= 12 ? 2 : 3));
        Key k{col, szb};
        auto it = pm.find(k);
        if (it == pm.end())
        {
          G4Polymarker m;
          m.SetMarkerType(G4Polymarker::circles);
          m.SetFillStyle(G4VMarker::filled);
          m.SetScreenSize(C.clusScreen * (1.0 + 0.55 * szb));
          // "Other" is chrome, not a series: muted where it is a reported
          // category (class 2 = noise), darker where it is only background
          // behind a spotlight.
          G4Colour c = (C.colour == "track") ? PAL::OTHER : PAL::MUTED;
          if (col == 200) c = PAL::ORANGE;
          else if (col == 0) c = PAL::BLUE;
          else if (col == 1) c = PAL::ORANGE;
          else if (col >= 100) c = PAL::SLOT[(col - 100) % PAL::NSLOT];
          G4VisAttributes va(c);
          m.SetVisAttributes(va);
          it = pm.emplace(k, m).first;
        }
        it->second.push_back(P(cl.x, cl.y, cl.z));
      }
      for (auto &kv : pm) vis->Draw(kv.second);
    }

    if (C.showTrk && !g_trk.empty())
    {
      std::unordered_map<int, int> slot;
      if (C.colour == "track")
      {
        if (C.trkSelect != 0) slot[C.trkSelect] = 0;
        else
          for (int i = 0; i < (int) g_trk.size() && i < C.trkTop; ++i) slot[g_trk[i].id] = i;
      }
      for (auto &tk : g_trk)
      {
        G4Colour c = PAL::AQUA;
        double w = 1.5;
        if (C.colour == "track")
        {
          auto s = slot.find(tk.id);
          if (s == slot.end()) { c = PAL::GRID; w = 1.0; }
          else { c = PAL::SLOT[s->second % PAL::NSLOT]; w = 2.5; }
        }
        else if (C.colour == "class") { c = PAL::AQUA; }
        G4VisAttributes va(c);
        va.SetLineWidth(w);
        // break the chain where the r-ordering would draw a line that is not a
        // trajectory (looper wrap-around / a second truth segment on the far side)
        G4Polyline pl;
        pl.SetVisAttributes(va);
        const Clus *prev = nullptr;
        for (auto *cl : tk.cl)
        {
          if (prev)
          {
            double dz = std::fabs(cl->z - prev->z);
            double dphi = std::fabs(std::atan2(cl->y, cl->x) - std::atan2(prev->y, prev->x));
            if (dphi > M_PI) dphi = 2 * M_PI - dphi;
            if (dz > 20. || dphi > 0.35)
            {
              if (pl.size() > 1) vis->Draw(pl);
              pl.clear();
              pl.SetVisAttributes(va);
            }
          }
          pl.push_back(P(cl->x, cl->y, cl->z));
          prev = cl;
        }
        if (pl.size() > 1) vis->Draw(pl);
      }
    }
  }

  void drawLibrary(G4VVisManager *vis) const
  {
    if (C.showG4 && !g_g4.empty())
    {
      // edep spans decades and a handful of delta rays sit far above the bulk,
      // so the ramp tops out at the 95th percentile rather than the maximum
      std::vector<double> e;
      e.reserve(g_g4.size());
      for (auto &h : g_g4) e.push_back(1e6 * h.edep);  // GeV -> keV
      std::vector<double> s = e;
      std::nth_element(s.begin(), s.begin() + (size_t) (0.95 * (s.size() - 1)), s.end());
      double top = std::max(1.0, s[(size_t) (0.95 * (s.size() - 1))]);
      std::vector<std::pair<G4Point3D, double>> pts;
      pts.reserve(g_g4.size());
      for (size_t i = 0; i < g_g4.size(); ++i)
        pts.push_back({P(g_g4[i].x, g_g4[i].y, g_g4[i].z), e[i]});
      drawPix(vis, pts, false, 0, top);
    }
    if (C.showTrk)
      for (auto &kv : g_g4trk)
      {
        if (kv.second.size() < 3) continue;
        const bool prim = kv.second.front()->prim > 0;
        if (C.trkPrimaryOnly && !prim) continue;
        G4VisAttributes va(prim ? PAL::AQUA : PAL::OTHER);
        va.SetLineWidth(prim ? 1.8 : 1.0);
        G4Polyline pl;
        pl.SetVisAttributes(va);
        for (auto *h : kv.second) pl.push_back(P(h->x, h->y, h->z));
        if (pl.size() > 1) vis->Draw(pl);
      }
  }

  // geometry systems that carry a hue get named, so the tint is never the only
  // thing identifying them. Only the ones actually on screen, and only past the
  // default tpc preset (where everything is neutral anyway).
  template <class KeyFn>
  static void geomKeys(KeyFn &key, double &y)
  {
    if (!g_chromeTint) return;
    for (size_t i = 0; i < sysCount(); ++i)
      if (sysTinted(i) && i < g_sysOn.size() && g_sysOn[i]) key(sysColour(i), sysName(i));
    (void) y;
  }

  // 2D screen legend. Text always wears an ink token; the coloured mark beside
  // it is what carries identity — never the label itself.
  void drawLegend(G4VVisManager *vis) const
  {
    // Text is sized in points and the window is not, so the row pitch is set
    // for the smallest window worth using (~500 px tall) and simply reads airy
    // on a big one.
    const double X = -0.975, DY = 0.070;
    double y = 0.945;
    auto title = [&](const char *txt, const G4Colour &ink, double pt) {
      G4Text t(txt, G4Point3D(X, y, 0));
      t.SetScreenSize(pt);
      G4VisAttributes va(ink);
      t.SetVisAttributes(va);
      vis->Draw2D(t);
      y -= DY;
    };
    auto key = [&](const G4Colour &swatch, const char *txt) {
      G4VisAttributes vm(swatch);
      G4Square m(G4Point3D(X + 0.014, y + 0.011, 0));
      m.SetScreenSize(8);
      m.SetFillStyle(G4VMarker::filled);
      m.SetVisAttributes(vm);
      vis->Draw2D(m);
      G4Text t(txt, G4Point3D(X + 0.035, y, 0));
      t.SetScreenSize(11);
      G4VisAttributes va(PAL::INK);
      t.SetVisAttributes(va);
      vis->Draw2D(t);
      y -= DY;
    };
    char b[192];

    snprintf(b, sizeof b, "sPHENIX P5 %s   %s", C.ver.c_str(),
             C.g4event >= 0 ? "library event" : "production frame");
    title(b, PAL::INK, 13);
    title(g_status.c_str(), PAL::MUTED, 11);

    if (C.g4event >= 0)
    {
      if (C.showG4) key(PAL::RAMP[PAL::NRAMP - 1], "G4 step, edep (dark to light)");
      if (C.showTrk) key(PAL::AQUA, "G4 primary track (grey: secondary)");
      geomKeys(key, y);
    }
    else
    {
      if (C.showPix) key(PAL::RAMP[PAL::NRAMP - 1], "digitised pixel, ADC (dark to light)");
      if (C.showClus)
      {
        if (C.colour == "class")
        {
          key(PAL::BLUE, "cluster, truth class 0 (track-like)");
          key(PAL::ORANGE, "cluster, truth class 1 (looper)");
          key(PAL::MUTED, "cluster, truth class 2 (noise)");
        }
        else if (C.colour == "track")
        {
          int n = C.trkSelect ? 1 : std::min<int>(C.trkTop, (int) g_trk.size());
          for (int i = 0; i < n; ++i)
          {
            if (C.trkSelect) snprintf(b, sizeof b, "truth track %d", C.trkSelect);
            else snprintf(b, sizeof b, "truth track %d, pT %.2f GeV", g_trk[i].id, g_trk[i].pt);
            key(PAL::SLOT[i % PAL::NSLOT], b);
          }
          // no swatch: the dimmed grey is deliberately below swatch legibility,
          // so it is stated as a note rather than keyed as a series
          title("every other track and cluster dimmed", PAL::MUTED, 11);
        }
        else key(PAL::ORANGE, "island91 cluster (mark size = clus size)");
      }
      if (C.showTrk && C.colour != "track") key(PAL::AQUA, "truth cluster chain");

      geomKeys(key, y);
      y -= 0.012;
      snprintf(b, sizeof b, "z %.0f..%.0f cm   rows %d-%d   ADC >= %.0f", C.zLo, C.zHi, C.layerLo,
               C.layerHi, C.adcMin);
      title(b, PAL::MUTED, 11);
      snprintf(b, sizeof b, "sector %s   truth chain: >=%d clus, pT >= %.2f GeV%s",
               C.sector < 0 ? "all" : std::to_string(C.sector).c_str(), C.trkMinClus, C.trkPtMin,
               g_pixStride > 1 ? "   [pixels sub-sampled]" : "");
      title(b, PAL::MUTED, 11);
    }
  }
};
static Overlay *g_overlay = nullptr;

// =====================================================================
//  geometry presets
// =====================================================================
static std::vector<G4VisAttributes *> g_ownedVis;

static std::string baseName(const G4String &n)
{
  std::string s = n;
  size_t p = s.find("0x");
  return p == std::string::npos ? s : s.substr(0, p);
}

// Detector systems. Tints are held at 0.65x luminance and 0.95x chroma of the
// hues they name: MEASURED, because at full strength six of the fifteen
// geometry-vs-data pairs sat below the 15 normal-vision dE floor (worst was
// HCal vs the truth chains at 8.1) and the detector competed with the overlay.
// Luminance is the lever — desaturating alone only moved that 8.1 to 12.0,
// while dropping luminance clears every pair at >=16.9 with the hue intact.
// Uniform grey is unreadable once more than the TPC is on
// screen — recognition needs DIFFERENTIATION, not brightness — so each system
// gets its own tint and a legend row. The tints are deliberately dark and
// desaturated, and they avoid the overlay's blue/orange/aqua entirely, so the
// data still wins; identity is never on colour alone because the systems that
// carry a hue are named in the legend.
struct GeoSys
{
  const char *name;
  int tier;      // 0 tpc, 1 tracker, 2 calo, 3 beam line — matches the presets
  bool tinted;   // carries a hue and therefore a legend row
  G4Colour col;
  std::vector<std::string> pfx;   // empty = catch-all (must be last)
};
std::vector<GeoSys> G_SYS = {
    // the TPC is the subject, so it is the brightest thing and it gets a legend
    // row even though its swatch is neutral rather than a hue
    {"TPC gas", 0, true, G4Colour(0.493, 0.482, 0.445), {"tpc_gas"}},
    {"TPC cage", 0, false, G4Colour(0.296, 0.288, 0.259),
     {"tpc_cage_layer_", "tpc_window", "tpc_envelope"}},
    // the endcap/wagon-wheel assembly is ~2000 touchables and buries the
    // overlay at event-display zoom, so it joins only from 'tracker' outwards
    {"TPC endcap", 1, false, G4Colour(0.221, 0.214, 0.192),
     {"TPC_ENDCAP_", "tpc_hanger", "tpc_tie_rod"}},
    {"MVTX", 1, true, G4Colour(0.478, 0.183, 0.333),
     {"log_MVTX_", "MVTX", "CYSS", "ConeL", "EndWheel"}},
    {"INTT", 1, true, G4Colour(0.273, 0.229, 0.500),
     {"ladder_", "ladderext_", "stave_", "staveext_", "siactive", "siinactive", "si_glue",
      "fphx", "hdi_", "hdiext_", "rail_volume", "outer_skin_volume", "inner_skin_volume",
      "service_barrel_", "support_tube_volume", "endcap_AlPEEK_", "bus_extender_"}},
    {"TPOT", 1, true, G4Colour(0.429, 0.336, 0.026),
     {"MICROMEGAS_55_", "invisible_MICROMEGAS_55_", "micromegas_"}},
    {"beam pipe", 1, false, G4Colour(0.274, 0.263, 0.229),
     {"BE_PIPE", "VAC_BE_PIPE", "N_AL_PIPE", "S_AL_PIPE", "VAC_N_", "VAC_S_", "N_FLANGE_",
      "S_FLANGE_", "N_OUTER_PIPE_", "S_OUTER_PIPE_"}},
    {"CEMC", 2, true, G4Colour(0.472, 0.189, 0.138), {"CEMC_"}},
    {"HCal / cryostat", 2, true, G4Colour(0.166, 0.356, 0.170),
     {"HCAL_", "CRYOSTAT", "CRYOINT", "CONNECTOR", "BusbarD"}},
    {"EPD", 2, false, G4Colour(0.149, 0.146, 0.134), {"EPD_tile_"}},
    {"beam line", 3, false, G4Colour(0.142, 0.138, 0.123), {}},
};
static std::vector<G4VisAttributes *> g_sysVis;
std::vector<char> g_sysOn;                 // is this system currently on screen
static double g_chromeAlpha = 1.0;
static double g_chromeBright = 1.0;
bool g_chromeTint = true;                  // false = flat grey

// ONE colour map for both overlay modes: a frame and a library event show the
// same detector in the same colours. Density in frame mode is handled where it
// belongs — /p5/pix/adcmin, /p5/sector, a lighter /p5/geometry preset — not by
// silently repainting the geometry differently depending on what is loaded.
static double effBright() { return g_chromeBright; }

// the legend swatch must be the colour actually on screen, brightness and
// tint-toggle included, or the map lies about the render
G4Colour renderedSysColour(size_t i)
{
  G4Colour c = G_SYS[i].col;
  if (!g_chromeTint)
  {
    double y = 0.30 * c.GetRed() + 0.59 * c.GetGreen() + 0.11 * c.GetBlue();
    c = G4Colour(y, y, y);
  }
  auto lift = [](double v) { return std::min(1.0, v * effBright()); };
  return G4Colour(lift(c.GetRed()), lift(c.GetGreen()), lift(c.GetBlue()));
}

const char *sysName(size_t i) { return G_SYS[i].name; }
bool sysTinted(size_t i) { return G_SYS[i].tinted; }
G4Colour sysColour(size_t i) { return renderedSysColour(i); }
size_t sysCount() { return G_SYS.size(); }
static std::string g_geomStyle = "auto";
static G4VisAttributes *INVIS = nullptr;

static int sysOf(const std::string &base)
{
  for (size_t i = 0; i < G_SYS.size(); ++i)
    for (const auto &p : G_SYS[i].pfx)
      if (base.compare(0, p.size(), p) == 0) return (int) i;
  return (int) G_SYS.size() - 1;  // catch-all
}

static void ensureChrome()
{
  if (!g_sysVis.empty()) return;
  for (size_t i = 0; i < G_SYS.size(); ++i)
  {
    auto *va = new G4VisAttributes(G_SYS[i].col);
    va->SetVisibility(true);
    va->SetLineWidth(1.0);
    g_ownedVis.push_back(va);
    g_sysVis.push_back(va);
  }
  g_sysOn.assign(G_SYS.size(), 0);
  if (!INVIS) INVIS = new G4VisAttributes(false);
  applyChromeColours();
}

// The drawing style is NOT forced by default: a per-volume forced style
// overrides /vis/viewer/set/style. 'auto' hands that command back to the user.
static void applyChromeColours()
{
  for (size_t i = 0; i < G_SYS.size(); ++i)
  {
    G4Colour c = G_SYS[i].col;
    if (!g_chromeTint)  // flat grey: keep the luminance, drop the hue
    {
      double y = 0.30 * c.GetRed() + 0.59 * c.GetGreen() + 0.11 * c.GetBlue();
      c = G4Colour(y, y, y);
    }
    const double f = effBright();
    auto lift = [f](double v) { return std::min(1.0, v * f); };
    g_sysVis[i]->SetColour(
        G4Colour(lift(c.GetRed()), lift(c.GetGreen()), lift(c.GetBlue()), g_chromeAlpha));
    if (g_geomStyle == "wireframe") g_sysVis[i]->SetForceWireframe(true);
    else if (g_geomStyle == "surface") g_sysVis[i]->SetForceSolid(true);
    else g_sysVis[i]->SetForceWireframe(false);
  }
}

static bool viewerLive()
{
  auto *vm = G4VisManager::GetInstance();
  return vm && vm->GetCurrentViewer();
}

static void geomStyle(const std::string &style, double alpha)
{
  ensureChrome();
  g_geomStyle = style;
  g_chromeAlpha = alpha;
  applyChromeColours();
  auto *UI = G4UImanager::GetUIpointer();
  if (style == "surface") UI->ApplyCommand("/vis/viewer/set/style surface");
  else if (style == "wireframe") UI->ApplyCommand("/vis/viewer/set/style wireframe");
  printf("p5_display: geometry style '%s', alpha %.2f, brightness %.2f, tint %s%s\n",
         style.c_str(), alpha, effBright(), g_chromeTint ? "on" : "off",
         style == "surface" ? "  (alpha is what keeps the overlay visible)" : "");
  if (viewerLive()) UI->ApplyCommand("/vis/viewer/rebuild");
}

static void applyGeom(const std::string &preset)
{
  auto *store = G4LogicalVolumeStore::GetInstance();
  ensureChrome();
  int tier;
  if (preset == "none") tier = -1;
  else if (preset == "tpc") tier = 0;
  else if (preset == "tracker") tier = 1;
  else if (preset == "calo") tier = 2;
  else if (preset == "all") tier = 3;
  else
  {
    printf("p5_display: /p5/geometry tpc | tracker | calo | all | none\n");
    return;
  }
  std::fill(g_sysOn.begin(), g_sysOn.end(), 0);
  std::vector<int> n(G_SYS.size(), 0);
  int nv = 0;
  for (auto *lv : *store)
  {
    int k = sysOf(baseName(lv->GetName()));
    bool on = G_SYS[k].tier <= tier;
    lv->SetVisAttributes(on ? g_sysVis[k] : INVIS);
    if (on) { nv++; n[k]++; g_sysOn[k] = 1; }
  }
  printf("p5_display: geometry '%s' -> %d/%zu volumes%s\n", preset.c_str(), nv, store->size(),
         preset == "all" ? "  (66k touchables — interactive rotation gets slower)" : "");
  for (size_t i = 0; i < G_SYS.size(); ++i)
    if (n[i]) printf("            %-16s %4d\n", G_SYS[i].name, n[i]);
}

// free-form: any GDML volume by name substring, on top of the current preset.
// Shown volumes keep their own system's colour, so a pick stays recognisable.
static void geomSelect(const std::string &pattern, bool on)
{
  auto *store = G4LogicalVolumeStore::GetInstance();
  ensureChrome();
  int n = 0;
  for (auto *lv : *store)
  {
    std::string b = baseName(lv->GetName());
    if (b.find(pattern) == std::string::npos) continue;
    int k = sysOf(b);
    lv->SetVisAttributes(on ? g_sysVis[k] : INVIS);
    if (on) g_sysOn[k] = 1;
    n++;
  }
  printf("p5_display: geometry %s '%s' -> %d logical volume(s)%s\n", on ? "show" : "hide",
         pattern.c_str(), n, n ? "" : "  (nothing matched — try /p5/geom/list)");
}


static void geomList(const std::string &pattern)
{
  std::map<std::string, int> names;
  for (auto *lv : *G4LogicalVolumeStore::GetInstance())
  {
    std::string b = baseName(lv->GetName());
    if (pattern.empty() || b.find(pattern) != std::string::npos) names[b]++;
  }
  printf("  %zu distinct GDML volume names%s:\n", names.size(),
         pattern.empty() ? "" : (" matching '" + pattern + "'").c_str());
  int col = 0;
  for (auto &kv : names)
  {
    printf("  %-38s", kv.first.c_str());
    if (++col % 3 == 0) printf("\n");
  }
  if (col % 3) printf("\n");
  printf("  use /p5/geom/show <substring> to add any of them to the scene\n");
}

// =====================================================================
//  UI messenger
// =====================================================================
// /vis/viewer/zoomTo is relative to the scene's bounding sphere, and Geant4
// computes that sphere from the VISIBLE volumes — so the same zoom factor frames
// the TPC with preset 'tpc' and the whole beam line with preset 'all'. Every
// framing here is therefore expressed as "show this many cm of radius" and the
// zoom factor is derived from the live scene extent instead of hard-coded.
static double g_viewRadius = 150.;  // cm

static void applyZoom()
{
  auto *vm = G4VisManager::GetInstance();
  G4Scene *sc = vm ? vm->GetCurrentScene() : nullptr;
  if (!sc || g_viewRadius <= 0) return;
  double R = sc->GetExtent().GetExtentRadius() / cm;
  if (!(R > 0)) return;
  char b[96];
  snprintf(b, sizeof b, "/vis/viewer/zoomTo %.4f", R / g_viewRadius);
  G4UImanager::GetUIpointer()->ApplyCommand(b);
}

static void refresh(bool reloadData)
{
  if (reloadData) reload();
  else makeStatus();
  auto *UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand("/vis/scene/notifyHandlers");
  UI->ApplyCommand("/vis/viewer/flush");
}

class P5Messenger : public G4UImessenger
{
 public:
  P5Messenger()
  {
    dir = new G4UIdirectory("/p5/");
    dir->SetGuidance("P5 v5.5 batch-output event display");
    dirPix = new G4UIdirectory("/p5/pix/");
    dirTrk = new G4UIdirectory("/p5/trk/");

    cFrame = mk<G4UIcmdWithAnInteger>("/p5/frame", "production frame index (0..249 for 5x50)");
    cG4ev = mk<G4UIcmdWithAString>(
        "/p5/g4event",
        "library event(s): <n> | <lo>-<hi> | <a>,<b>,<c-d> superimposed; -1 = frame mode");
    cG4ch = mk<G4UIcmdWithAnInteger>("/p5/g4chunk", "PP_g4hit chunk index for /p5/g4event");
    cG4List = mk<G4UIcmdWithAnInteger>("/p5/g4list",
                                       "print the N busiest library events in the current chunk");
    cReload = mk<G4UIcmdWithoutParameter>("/p5/reload", "re-read the current selection");
    cPrint = mk<G4UIcmdWithoutParameter>("/p5/print", "print the current selection summary");
    dirGeom = new G4UIdirectory("/p5/geom/");
    cGeom = mk<G4UIcmdWithAString>("/p5/geometry",
                                   "tpc | tracker | calo | all | none");
    cGShow = mk<G4UIcmdWithAString>("/p5/geom/show", "add every GDML volume whose name contains this");
    cGHide = mk<G4UIcmdWithAString>("/p5/geom/hide", "remove every GDML volume whose name contains this");
    cGList = mk<G4UIcmdWithAString>("/p5/geom/list", "print GDML volume names (optional substring filter)");
    cGTint = mk<G4UIcmdWithABool>("/p5/geom/tint",
                                  "colour the geometry per detector system (false = flat grey)");
    cGBright = mk<G4UIcmdWithADouble>("/p5/geom/bright",
                                      "chrome brightness multiplier (1 = the tuned default)");
    cGStyle = new G4UIcommand("/p5/geom/style", this);
    cGStyle->SetGuidance("/p5/geom/style wireframe | surface [alpha] | auto");
    cGStyle->SetGuidance("  auto = let /vis/viewer/set/style decide (nothing forced per volume)");
    { auto *p = new G4UIparameter("style", 's', false); cGStyle->SetParameter(p); }
    // sentinel, not 0.12: the sensible alpha depends on the style (opaque lines in
    // wireframe, translucent shells in surface), so the handler resolves it
    { auto *p = new G4UIparameter("alpha", 'd', true); p->SetDefaultValue(-1.); cGStyle->SetParameter(p); }
    cColour = mk<G4UIcmdWithAString>("/p5/colour", "level | class | track");
    cLegend = mk<G4UIcmdWithABool>("/p5/legend", "draw the on-screen legend");
    cView = mk<G4UIcmdWithAString>("/p5/view", "3d | rphi | rz | iso  (viewpoint + zoom preset)");
    cSector = mk<G4UIcmdWithAnInteger>("/p5/sector", "TPC 30-degree sector 0..11; -1 = all");
    cZoom = mk<G4UIcmdWithADouble>("/p5/zoom", "frame this radius [cm] (TPC bounding radius ~131)");

    cShow = new G4UIcommand("/p5/show", this);
    cShow->SetGuidance("/p5/show <pix|clus|trk|g4|all> <true|false>");
    { auto *p = new G4UIparameter("layer", 's', false); cShow->SetParameter(p); }
    { auto *p = new G4UIparameter("on", 'b', true); p->SetDefaultValue(true); cShow->SetParameter(p); }

    cLayers = new G4UIcommand("/p5/layers", this);
    cLayers->SetGuidance("/p5/layers <lo> <hi>   TPC row range (7..54)");
    { auto *p = new G4UIparameter("lo", 'i', false); cLayers->SetParameter(p); }
    { auto *p = new G4UIparameter("hi", 'i', false); cLayers->SetParameter(p); }

    cZ = new G4UIcommand("/p5/zrange", this);
    cZ->SetGuidance("/p5/zrange <lo> <hi> [cm], or /p5/zrange full  (apparent z, +-305 cm)");
    { auto *p = new G4UIparameter("lo", 's', false); cZ->SetParameter(p); }
    { auto *p = new G4UIparameter("hi", 's', true); p->SetDefaultValue("0"); cZ->SetParameter(p); }

    cAdcMin = mk<G4UIcmdWithADouble>("/p5/pix/adcmin", "drop pixels below this ADC");
    cAdcMax = mk<G4UIcmdWithADouble>("/p5/pix/adcmax", "ADC at the top of the colour ramp");
    cPixMax = mk<G4UIcmdWithAnInteger>("/p5/pix/max", "cap on drawn pixels (stride-sampled above)");
    cPixSize = mk<G4UIcmdWithADouble>("/p5/pix/size", "0 = 1px dots; >0 = filled circles, screen px");
    cRowDr = mk<G4UIcmdWithABool>("/p5/pix/rowdr", "apply the v5.4 row-radius offsets to pixel r");
    cClusSize = mk<G4UIcmdWithADouble>("/p5/clussize", "base screen size of a cluster mark [px]");

    cMinClus = mk<G4UIcmdWithAnInteger>("/p5/trk/minclus", "minimum clusters for a truth chain");
    cPtMin = mk<G4UIcmdWithADouble>("/p5/trk/ptmin", "minimum truth pT [GeV] for a chain");
    cPrim = mk<G4UIcmdWithABool>("/p5/trk/primary", "keep only primary truth tracks");
    cTop = mk<G4UIcmdWithAnInteger>("/p5/trk/top", "spotlight slots in /p5/colour track");
    cSelect = mk<G4UIcmdWithAnInteger>("/p5/trk/select", "spotlight one truth id (0 = off)");
    cList = mk<G4UIcmdWithAnInteger>("/p5/trk/list", "print the N highest-pT truth chains");

    cField = mk<G4UIcmdWithAString>("/p5/field", "load the field map for /run/beamOn, or 'off'");
  }

  void SetNewValue(G4UIcommand *cmd, G4String val) override
  {
    if (cmd == cFrame) { C.frame = std::stoi(val); C.g4event = -1; C.g4evList.clear(); refresh(true); }
    else if (cmd == cG4ev) { setG4Events(val); refresh(true); }
    else if (cmd == cG4ch) { C.g4chunk = std::stoi(val); refresh(true); }
    else if (cmd == cReload) refresh(true);
    else if (cmd == cPrint) printSummary();
    else if (cmd == cGeom)
    {
      applyGeom(val);
      G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/rebuild");
      applyZoom();  // the visible set changed -> so did the scene extent
    }
    else if (cmd == cGShow || cmd == cGHide)
    {
      geomSelect(val, cmd == cGShow);
      G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/rebuild");
      applyZoom();
    }
    else if (cmd == cGList) geomList(val);
    else if (cmd == cGTint)
    {
      g_chromeTint = G4UIcmdWithABool::GetNewBoolValue(val);
      geomStyle(g_geomStyle, g_chromeAlpha);
      refresh(false);
    }
    else if (cmd == cGBright)
    {
      g_chromeBright = std::max(0.0, std::stod(val));
      geomStyle(g_geomStyle, g_chromeAlpha);  // re-tint at the current style
    }
    else if (cmd == cGStyle)
    {
      std::istringstream ss(val);
      std::string st;
      double al = -1;
      ss >> st;
      if (!(ss >> al)) al = -1;
      if (al < 0) al = (st == "surface") ? 0.12 : 1.0;
      if (st != "wireframe" && st != "surface" && st != "auto")
        printf("p5_display: /p5/geom/style wireframe | surface [alpha] | auto\n");
      else geomStyle(st, al);
    }
    else if (cmd == cColour)
    {
      std::string m = (val == "layers") ? "level" : val;  // pre-rename alias
      if (m != "level" && m != "class" && m != "track")
        printf("p5_display: unknown colour mode '%s' (level|class|track)\n", val.c_str());
      else { C.colour = m; refresh(false); }
    }
    else if (cmd == cLegend) { C.showLegend = G4UIcmdWithABool::GetNewBoolValue(val); refresh(false); }
    else if (cmd == cSector)
    {
      C.sector = std::stoi(val);
      if (C.sector > 11) C.sector = 11;
      refresh(true);
    }
    else if (cmd == cView) setView(val);
    else if (cmd == cZoom)
    {
      g_viewRadius = std::stod(val);
      applyZoom();
      G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/flush");
    }
    else if (cmd == cClusSize) { C.clusScreen = std::stod(val); refresh(false); }
    else if (cmd == cShow)
    {
      std::istringstream ss(val);
      std::string what, on = "true";
      ss >> what >> on;
      bool b = (on == "true" || on == "1" || on == "T" || on == "y");
      if (what == "all") C.showPix = C.showClus = C.showTrk = C.showG4 = b;
      else if (what == "pix") C.showPix = b;
      else if (what == "clus") C.showClus = b;
      else if (what == "trk") C.showTrk = b;
      else if (what == "g4") C.showG4 = b;
      else { printf("p5_display: /p5/show pix|clus|trk|g4|all <bool>\n"); return; }
      refresh(false);
    }
    else if (cmd == cLayers)
    {
      std::istringstream ss(val);
      ss >> C.layerLo >> C.layerHi;
      refresh(true);
    }
    else if (cmd == cZ)
    {
      std::istringstream ss(val);
      std::string a, b;
      ss >> a >> b;
      if (a == "full") { C.zLo = -310; C.zHi = 310; }
      else if (a == "tpc") { C.zLo = -105.5; C.zHi = 105.5; }
      else { C.zLo = std::stod(a); C.zHi = std::stod(b); }
      if (C.zHi - C.zLo > 260)
        printf("p5_display: the apparent-z volume is ~3x the TPC length — the view does not "
               "re-fit itself, use /vis/viewer/zoomTo 1.1 to see all of it\n");
      refresh(true);
    }
    else if (cmd == cAdcMin) { C.adcMin = std::stod(val); refresh(true); }
    else if (cmd == cAdcMax) { C.adcMax = std::stod(val); refresh(false); }
    else if (cmd == cPixMax) { C.pixMax = std::stol(val); refresh(true); }
    else if (cmd == cPixSize) { C.pixScreen = std::stod(val); C.pixScreenSet = true; refresh(false); }
    else if (cmd == cRowDr) { C.applyRowDr = G4UIcmdWithABool::GetNewBoolValue(val); refresh(true); }
    else if (cmd == cMinClus) { C.trkMinClus = std::stoi(val); refresh(true); }
    else if (cmd == cPtMin) { C.trkPtMin = std::stod(val); refresh(true); }
    else if (cmd == cPrim) { C.trkPrimaryOnly = G4UIcmdWithABool::GetNewBoolValue(val); refresh(true); }
    else if (cmd == cTop)
    {
      C.trkTop = std::max(1, std::stoi(val));
      if (C.trkTop > 3)
        printf("p5_display: NOTE slots beyond 3 are not all-pairs separable under CVD "
               "(yellow/orange collide) — use /p5/trk/select for a clean single spotlight\n");
      refresh(false);
    }
    else if (cmd == cSelect)
    {
      C.trkSelect = std::stoi(val);
      if (C.trkSelect) C.colour = "track";  // 0 = clear the spotlight, not enter it
      refresh(false);
    }
    else if (cmd == cList) listTracks(std::stoi(val));
    else if (cmd == cG4List) listG4Events(std::stoi(val));
    else if (cmd == cField) loadField(val);
  }

 private:
  template <class T>
  T *mk(const char *path, const char *guid)
  {
    T *c = new T(path, this);
    c->SetGuidance(guid);
    return c;
  }
  void printSummary() const
  {
    printf("---- p5_display -------------------------------------------------\n");
    printf("  version   %s\n", C.ver.c_str());
    printf("  mode      %s\n", C.g4event >= 0 ? "library (g4hit)" : "frame (production)");
    printf("  %s\n", g_status.c_str());
    printf("  layers %d-%d   sector %s   z [%.1f,%.1f] cm   adc >= %.0f (ramp top %.0f)\n",
           C.layerLo, C.layerHi, C.sector < 0 ? "all" : std::to_string(C.sector).c_str(), C.zLo,
           C.zHi, C.adcMin, C.adcMax);
    printf("  tracks: >=%d clusters, pT >= %.2f GeV%s\n", C.trkMinClus, C.trkPtMin,
           C.trkPrimaryOnly ? ", primary only" : "");
    printf("  colour    %s     pixel stride %d\n", C.colour.c_str(), g_pixStride);
    printf("  sources\n    %s\n    %s\n", C.digi.c_str(), C.isl.c_str());
    printf("-----------------------------------------------------------------\n");
  }
  void listTracks(int n) const
  {
    printf("  %-10s %-8s %-6s %-8s %s\n", "truth id", "pT[GeV]", "prim", "nclus", "r range [cm]");
    for (int i = 0; i < (int) g_trk.size() && i < n; ++i)
    {
      const Trk &t = g_trk[i];
      printf("  %-10d %-8.3f %-6d %-8zu %.1f .. %.1f\n", t.id, t.pt, t.prim, t.cl.size(),
             t.cl.front()->r, t.cl.back()->r);
    }
  }
  void loadField(const G4String &v);
  static void setView(const G4String &v)
  {
    auto *UI = G4UImanager::GetUIpointer();
    // the scene extent is the GDML world, so every preset re-targets and
    // re-zooms onto the TPC rather than trusting the default fit
    UI->ApplyCommand("/vis/viewer/set/targetPoint 0 0 0 cm");
    // viewpointThetaPhi is the direction FROM the target TO the camera:
    // (0,0) looks down -z = r-phi; (90,0) looks down -x = r-z with z across.
    // radius of the physical region each preset frames, in cm — the TPC is
    // r 31..76, |z| < 105.5, i.e. a bounding radius of ~131
    if (v == "rphi") { UI->ApplyCommand("/vis/viewer/set/viewpointThetaPhi 0 0"); g_viewRadius = 88.; }
    else if (v == "rz") { UI->ApplyCommand("/vis/viewer/set/viewpointThetaPhi 90 0"); g_viewRadius = 135.; }
    else if (v == "iso") { UI->ApplyCommand("/vis/viewer/set/viewpointThetaPhi 55 35"); g_viewRadius = 155.; }
    else if (v == "3d") { UI->ApplyCommand("/vis/viewer/set/viewpointThetaPhi 72 20"); g_viewRadius = 140.; }
    else { printf("p5_display: /p5/view 3d|rphi|rz|iso\n"); return; }
    UI->ApplyCommand("/vis/viewer/set/upVector 0 1 0");
    applyZoom();
    UI->ApplyCommand("/vis/viewer/flush");
  }

  G4UIdirectory *dir, *dirPix, *dirTrk, *dirGeom;
  G4UIcmdWithAnInteger *cFrame, *cG4ch, *cPixMax, *cMinClus, *cTop, *cSelect, *cList,
      *cSector, *cG4List;
  G4UIcmdWithoutParameter *cReload, *cPrint;
  G4UIcmdWithAString *cGeom, *cGShow, *cGHide, *cGList, *cColour, *cField, *cView, *cG4ev;
  G4UIcommand *cGStyle;
  G4UIcmdWithABool *cLegend, *cPrim, *cRowDr, *cGTint;
  G4UIcmdWithADouble *cAdcMin, *cAdcMax, *cPixSize, *cPtMin, *cClusSize, *cZoom, *cGBright;
  G4UIcommand *cShow, *cLayers, *cZ;
};

// =====================================================================
//  field (lazy: the map is 313 MB and a display does not need it)
// =====================================================================
class LazyField : public G4MagneticField
{
 public:
  bool load(const std::string &fn)
  {
    TFile *f = TFile::Open(fn.c_str());
    if (!f || f->IsZombie())
    {
      printf("p5_display: field map open failed: %s\n", fn.c_str());
      return false;
    }
    TTree *t = (TTree *) f->Get("fieldmap");
    if (!t)
    {
      printf("p5_display: no 'fieldmap' tree in %s\n", fn.c_str());
      f->Close();
      return false;
    }
    float x, y, z, bx, by, bz;
    t->SetBranchAddress("x", &x); t->SetBranchAddress("y", &y); t->SetBranchAddress("z", &z);
    t->SetBranchAddress("bx", &bx); t->SetBranchAddress("by", &by); t->SetBranchAddress("bz", &bz);
    m_bx.assign((size_t) NX * NY * NZ, 0.f);
    m_by.assign((size_t) NX * NY * NZ, 0.f);
    m_bz.assign((size_t) NX * NY * NZ, 0.f);
    Long64_t n = t->GetEntries();
    long used = 0;
    for (Long64_t i = 0; i < n; ++i)
    {
      t->GetEntry(i);
      int ix = (int) lrint((x - X0) / DX), iy = (int) lrint((y - X0) / DX),
          iz = (int) lrint((z - Z0) / DX);
      if (ix < 0 || ix >= NX || iy < 0 || iy >= NY || iz < 0 || iz >= NZ) continue;
      size_t k = idx(ix, iy, iz);
      m_bx[k] = bx; m_by[k] = by; m_bz[k] = bz;
      used++;
    }
    f->Close();
    m_on = true;
    printf("p5_display: field map %ld/%lld points (%.0f MB) — /run/beamOn now tracks in field\n",
           used, n, 3.0 * m_bx.size() * sizeof(float) / 1048576.);
    return true;
  }
  void unload()
  {
    m_on = false;
    m_bx.clear(); m_by.clear(); m_bz.clear();
    m_bx.shrink_to_fit(); m_by.shrink_to_fit(); m_bz.shrink_to_fit();
    printf("p5_display: field off\n");
  }
  void GetFieldValue(const G4double p[4], G4double *B) const override
  {
    B[0] = B[1] = B[2] = 0;
    if (!m_on) return;
    double xc = p[0] / cm, yc = p[1] / cm, zc = p[2] / cm;
    double fx = (xc - X0) / DX, fy = (yc - X0) / DX, fz = (zc - Z0) / DX;
    int ix = (int) floor(fx), iy = (int) floor(fy), iz = (int) floor(fz);
    if (ix < 0 || ix >= NX - 1 || iy < 0 || iy >= NY - 1 || iz < 0 || iz >= NZ - 1) return;
    double ux = fx - ix, uy = fy - iy, uz = fz - iz;
    double b[3] = {0, 0, 0};
    for (int dx = 0; dx < 2; ++dx)
      for (int dy = 0; dy < 2; ++dy)
        for (int dz = 0; dz < 2; ++dz)
        {
          double w = (dx ? ux : 1 - ux) * (dy ? uy : 1 - uy) * (dz ? uz : 1 - uz);
          size_t k = idx(ix + dx, iy + dy, iz + dz);
          b[0] += w * m_bx[k]; b[1] += w * m_by[k]; b[2] += w * m_bz[k];
        }
    B[0] = b[0] * tesla; B[1] = b[1] * tesla; B[2] = b[2] * tesla;
  }

 private:
  static constexpr int NX = 271, NY = 271, NZ = 356;
  static constexpr double X0 = -270., Z0 = -355., DX = 2.;
  static size_t idx(int ix, int iy, int iz) { return ((size_t) ix * NY + iy) * NZ + iz; }
  std::vector<float> m_bx, m_by, m_bz;
  bool m_on = false;
};
static LazyField *g_field = nullptr;

void P5Messenger::loadField(const G4String &v)
{
  if (!g_field) return;
  if (v == "off" || v == "0") { g_field->unload(); return; }
  g_field->load(v.empty() || v == "default" ? C.fieldmap : std::string(v));
}

// =====================================================================
//  detector: the pipeline's GDML, unchanged
// =====================================================================
class Det : public G4VUserDetectorConstruction
{
 public:
  G4VPhysicalVolume *Construct() override
  {
    printf("p5_display: reading %s ...\n", C.gdml.c_str());
    m_parser.Read(C.gdml, false);
    // same CEMC effective fill the batch g4 stage applies, so a GUI /run/beamOn
    // sees the same material budget the production g4hits were made with
    auto *nist = G4NistManager::Instance();
    auto *eff = new G4Material("CEMC_eff", 8.0 * g / cm3, 2);
    eff->AddMaterial(nist->FindOrBuildMaterial("G4_W"), 0.92);
    eff->AddMaterial(nist->FindOrBuildMaterial("G4_POLYSTYRENE"), 0.08);
    int nfill = 0;
    for (auto *lv : *G4LogicalVolumeStore::GetInstance())
      if (baseName(lv->GetName()) == "CEMC_0") { lv->SetMaterial(eff); nfill++; }
    printf("p5_display: CEMC_0 effective fill applied to %d volume(s)\n", nfill);
    return m_parser.GetWorldVolume();
  }
  void ConstructSDandField() override
  {
    int nvol = 0;
    for (auto *lv : *G4LogicalVolumeStore::GetInstance())
      if (strncmp(lv->GetName().c_str(), "tpc_gas", 7) == 0)
      {
        lv->SetUserLimits(new G4UserLimits(1.0 * cm));  // PHG4 TPC max step
        nvol++;
      }
    printf("p5_display: 1 cm step limit on %d tpc_gas volume(s)\n", nvol);
    g_field = new LazyField;
    auto *fm = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    fm->SetDetectorField(g_field);
    fm->CreateChordFinder(g_field);
    if (!C.fieldmap.empty() && getenv("P5_FIELD_PRELOAD")) g_field->load(C.fieldmap);
  }

 private:
  G4GDMLParser m_parser;
};

// =====================================================================
//  primaries: the pipeline's own HepMC chunks, or a gun
// =====================================================================
struct HepPart { int pdg; double px, py, pz, e, vx, vy, vz, vt; };
class HepMC2Reader
{
 public:
  explicit HepMC2Reader(const std::string &fn) : m_in(fn) {}
  bool ok() const { return (bool) m_in; }
  bool next(std::vector<HepPart> &out)
  {
    out.clear();
    std::string line;
    bool inEvent = false;
    double vx = 0, vy = 0, vz = 0, vt = 0;
    std::streampos lastpos = m_in.tellg();
    while (std::getline(m_in, line))
    {
      if (line.empty()) continue;
      char tag = line[0];
      if (tag == 'E')
      {
        if (inEvent) { m_in.seekg(lastpos); return true; }
        inEvent = true;
      }
      else if (tag == 'V' && inEvent)
      {
        std::istringstream ss(line.substr(1));
        long bc; int id;
        ss >> bc >> id >> vx >> vy >> vz >> vt;
      }
      else if (tag == 'P' && inEvent)
      {
        std::istringstream ss(line.substr(1));
        long bc; int pdg; double px, py, pz, e, m; int status;
        ss >> bc >> pdg >> px >> py >> pz >> e >> m >> status;
        if (status == 1) out.push_back({pdg, px, py, pz, e, vx, vy, vz, vt});
      }
      lastpos = m_in.tellg();
    }
    return inEvent;
  }

 private:
  std::ifstream m_in;
};

class PrimGen : public G4VUserPrimaryGeneratorAction
{
 public:
  PrimGen()
  {
    if (!C.hepmc.empty())
    {
      m_r = new HepMC2Reader(C.hepmc);
      if (!m_r->ok())
      {
        printf("p5_display: HepMC %s unreadable — falling back to the gun\n", C.hepmc.c_str());
        delete m_r;
        m_r = nullptr;
      }
      else printf("p5_display: primaries from %s (/run/beamOn steps through it)\n", C.hepmc.c_str());
    }
    if (!m_r)
    {
      m_gun = new G4ParticleGun(1);
      m_gun->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("pi-"));
      m_gun->SetParticleMomentum(G4ThreeVector(1.0 * GeV, 0, 0.4 * GeV));
    }
  }
  void GeneratePrimaries(G4Event *ev) override
  {
    if (!m_r) { m_gun->GeneratePrimaryVertex(ev); return; }
    std::vector<HepPart> parts;
    if (!m_r->next(parts)) { printf("p5_display: HepMC exhausted\n"); return; }
    long placed = 0;
    for (auto &p : parts)
    {
      G4ParticleDefinition *def = (p.pdg > 1000000000)
                                      ? G4IonTable::GetIonTable()->GetIon(p.pdg)
                                      : G4ParticleTable::GetParticleTable()->FindParticle(p.pdg);
      if (!def) continue;
      auto *vtx = new G4PrimaryVertex(p.vx * mm, p.vy * mm, p.vz * mm, p.vt * mm / CLHEP::c_light);
      vtx->SetPrimary(new G4PrimaryParticle(def, p.px * GeV, p.py * GeV, p.pz * GeV));
      ev->AddPrimaryVertex(vtx);
      placed++;
    }
    printf("p5_display: event %d, %ld primaries\n", ev->GetEventID(), placed);
  }

 private:
  HepMC2Reader *m_r = nullptr;
  G4ParticleGun *m_gun = nullptr;
};

// =====================================================================
//  main
// =====================================================================
static void usage(const char *p)
{
  printf(
      "usage: %s [options]\n"
      "  --repo <dir>       repository root (default %s)\n"
      "  --ver <tag>        production tag (default v55)\n"
      "  --frame <n>        production frame to load (default 0)\n"
      "  --g4event <spec>   library event(s) instead of a frame: n | lo-hi | a,b,c-d\n"
      "  --g4chunk <i>      PP_g4hit chunk for --g4event (default 0)\n"
      "  --digi <file>      override digi_frames_production_<ver>.root\n"
      "  --island <file>    override island91_frames_production_<ver>.root\n"
      "  --gdml <file>      override sphenix_p5.gdml\n"
      "  --hepmc <file>     HepMC2 chunk for /run/beamOn (default: no gun input)\n"
      "  --field <file|off> field map; loaded lazily unless P5_FIELD_PRELOAD is set\n"
      "  --vis <driver>     Geant4 vis driver (default OGL)\n"
      "  --size <WxH>       window size hint (default 1400x900)\n"
      "  --session <type>   qt | tcsh | csh (default qt)\n"
      "  --macro <file>     extra macro executed after the vis setup\n"
      "  --batch            run the macro and exit (no interactive session)\n"
      "  --zfull            start with the full +-310 cm apparent-z window\n"
      "  --sector <0..11>   start on one TPC 30-degree sector (-1 = all)\n"
      "  --layers <lo> <hi> start on a TPC row range\n"
      "  --adcmin <v>       start with a pixel ADC floor\n"
      "  --pixsize <px>     0 = 1px dots; >0 = circles (library mode defaults 2.5)\n"
      "  --clussize <px>    base screen size of a cluster mark\n"
      "  --colour <mode>    level | class | track\n"
      "  --view <preset>    3d | rphi | rz | iso\n"
      "  --geomstyle <s> [a]  wireframe | surface [alpha] | auto\n"
      "  --geombright <f>   chrome brightness multiplier (1 = tuned default)\n"
      "  --geomtint <bool>  colour geometry per detector system (default true)\n"
      "  --zoom <radius_cm> frame this radius instead of the --view preset's\n"
      "  --cmd \"<command>\"  any /p5/ or /vis/ command, after the viewer opens;\n"
      "                     repeatable, applied in the order given\n"
      "  --geometry <p>     tpc | tracker | calo | all | none\n",
      p, C.repo.c_str());
}

int main(int argc, char **argv)
{
  std::string over_digi, over_isl, over_gdml;
  for (int i = 1; i < argc; ++i)
  {
    std::string a = argv[i];
    auto nxt = [&]() { return (i + 1 < argc) ? std::string(argv[++i]) : std::string(); };
    if (a == "--repo") C.repo = nxt();
    else if (a == "--ver") C.ver = nxt();
    else if (a == "--frame") C.frame = std::stoi(nxt());
    else if (a == "--g4event") setG4Events(nxt());
    else if (a == "--g4chunk") C.g4chunk = std::stoi(nxt());
    else if (a == "--digi") over_digi = nxt();
    else if (a == "--island") over_isl = nxt();
    else if (a == "--gdml") over_gdml = nxt();
    else if (a == "--hepmc") C.hepmc = nxt();
    else if (a == "--field") C.fieldmap = nxt();
    else if (a == "--vis") C.visdriver = nxt();
    else if (a == "--size") C.vissize = nxt();
    else if (a == "--session") C.session = nxt();
    else if (a == "--macro") C.macro = nxt();
    else if (a == "--cmd") C.cmds.push_back(nxt());
    else if (a == "--view") C.cmds.push_back("/p5/view " + nxt());
    else if (a == "--geombright") C.cmds.push_back("/p5/geom/bright " + nxt());
    else if (a == "--geomtint") C.cmds.push_back("/p5/geom/tint " + nxt());
    else if (a == "--zoom") C.cmds.push_back("/p5/zoom " + nxt());
    else if (a == "--geomstyle")
    {
      std::string st = nxt(), al;
      // optional trailing alpha: a style is never numeric, so a numeric
      // next token belongs to this flag
      if (i + 1 < argc)
      {
        std::string nx = argv[i + 1];
        if (!nx.empty() && (isdigit((unsigned char) nx[0]) || nx[0] == '.')) al = argv[++i];
      }
      C.cmds.push_back("/p5/geom/style " + st + (al.empty() ? "" : " " + al));
    }
    else if (a == "--batch") C.batch = true;
    else if (a == "--zfull") { C.zLo = -310; C.zHi = 310; }
    else if (a == "--sector") C.sector = std::stoi(nxt());
    else if (a == "--colour" || a == "--color") { C.colour = nxt(); if (C.colour == "layers") C.colour = "level"; }
    else if (a == "--geometry") C.geometry = nxt();
    else if (a == "--adcmin") C.adcMin = std::stod(nxt());
    else if (a == "--pixsize") { C.pixScreen = std::stod(nxt()); C.pixScreenSet = true; }
    else if (a == "--clussize") C.clusScreen = std::stod(nxt());
    else if (a == "--layers") { C.layerLo = std::stoi(nxt()); C.layerHi = std::stoi(nxt()); }
    else if (a == "-h" || a == "--help") { usage(argv[0]); return 0; }
    else { printf("unknown option %s\n", a.c_str()); usage(argv[0]); return 1; }
  }
  const std::string IP = C.repo + "/island_post", P5 = C.repo + "/P5";
  C.gdml = over_gdml.empty() ? P5 + "/sphenix_p5.gdml" : over_gdml;
  C.digi = over_digi.empty() ? IP + "/digi_frames_production_" + C.ver + ".root" : over_digi;
  C.isl = over_isl.empty() ? IP + "/island91_frames_production_" + C.ver + ".root" : over_isl;
  C.geomtab = IP + "/tpc_geom_table.txt";
  C.rowdr = IP + "/real_row_radii_v54.txt";
  C.g4pat = P5 + "/PP_g4hit_%d.root";
  if (C.fieldmap == "off") C.fieldmap.clear();

  printf("=================================================================\n");
  printf("  P5 %s event display — GDML geometry + batch final output\n", C.ver.c_str());
  printf("=================================================================\n");
  loadGeom();

  // The UI session must exist BEFORE /vis/open: the TSG / OGL "Qt" drivers
  // check for a live Qt application and silently fall back to the Xt viewer
  // (no GUI panel, no scene tree) if it is not there yet.
  G4UIExecutive *ui = nullptr;
  if (!C.batch)
  {
    int fake_argc = 1;
    char *fake_argv[2] = {argv[0], nullptr};
    ui = new G4UIExecutive(fake_argc, fake_argv, C.session);
  }

  auto *rm = new G4RunManager;
  rm->SetUserInitialization(new Det);
  auto *phys = new FTFP_BERT(0);
  phys->RegisterPhysics(new G4StepLimiterPhysics);
  rm->SetUserInitialization(phys);
  rm->SetUserAction(new PrimGen);
  rm->Initialize();

  auto *vis = new G4VisExecutive("quiet");
  vis->Initialize();
  g_overlay = new Overlay;
  vis->RegisterRunDurationUserVisAction(
      "P5Overlay", g_overlay,
      G4VisExtent(-90 * cm, 90 * cm, -90 * cm, 90 * cm, -320 * cm, 320 * cm));
  new P5Messenger;

  applyGeom(C.geometry);
  reload();

  auto *UI = G4UImanager::GetUIpointer();
  setenv("P5_VIS_DRIVER", C.visdriver.c_str(), 1);
  setenv("P5_VIS_SIZE", C.vissize.c_str(), 1);
  UI->ApplyCommand("/control/macroPath " + P5);
  UI->ApplyCommand("/control/execute vis_p5.mac");
  if (!C.macro.empty()) UI->ApplyCommand("/control/execute " + C.macro);
  for (const auto &c : C.cmds) UI->ApplyCommand(c);

  if (ui)
  {
    printf("\np5_display: ready. /p5/print for the summary, help /p5/ for the commands.\n\n");
    ui->SessionStart();
    delete ui;
  }
  delete vis;
  delete rm;
  return 0;
}
