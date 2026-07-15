// standalone_tpc.cc — P5 standalone Geant4 application (detachment step 1).
// Local Geant4 11.4.2 + exported sphenix_p5.gdml + local field map + HepMC2
// ASCII events -> ntp_g4hit (identical 12-branch layout to the container
// reference dumper). No Fun4All, no container, no network.
//
// Conventions matched to the reference chain (P5_g4hit_eval.root):
//   - hit = one G4 step in the tpc_gas volume with edep > 0
//   - gx,gy,gz = step midpoint [cm]; gt = pre-step global time [ns]
//   - gpl = step length [cm]; gpx,gpy,gpz = pre-step momentum [GeV]
//   - gedep [GeV]; gtrackID: primary > 0, secondary < 0 (sign convention
//     consumed by tpc_transport / islandize91 labeling)
//   - glayer written as 99: layer derived downstream from radius (transport
//     patch of 2026-07-15), exactly like the container dumper path.
//
// Build: see build_standalone.sh. Usage:
//   ./standalone_tpc <gdml> <hepmc.dat> <out.root> <nevents> <fieldmap.root>

#include <G4RunManager.hh>
#include <G4GDMLParser.hh>
#include <G4LogicalVolume.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4VSensitiveDetector.hh>
#include <G4SDManager.hh>
#include <G4MagneticField.hh>
#include <G4FieldManager.hh>
#include <G4TransportationManager.hh>
#include <G4PrimaryVertex.hh>
#include <G4PrimaryParticle.hh>
#include <G4ParticleTable.hh>
#include <G4IonTable.hh>
#include <G4Event.hh>
#include <G4Step.hh>
#include <G4Track.hh>
#include <G4SystemOfUnits.hh>
#include <FTFP_BERT.hh>
#include <G4StepLimiterPhysics.hh>
#include <G4UserLimits.hh>
#include <G4NistManager.hh>
#include <G4LogicalVolumeStore.hh>

#include <TFile.h>
#include <TNtuple.h>
#include <TTree.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

// ---------------- field: 2-cm Cartesian grid, trilinear ----------------
class MapField : public G4MagneticField
{
 public:
  explicit MapField(const std::string &fn)
  {
    TFile *f = TFile::Open(fn.c_str());
    if (!f || f->IsZombie()) { fprintf(stderr, "field map open failed: %s\n", fn.c_str()); exit(1); }
    TTree *t = (TTree *) f->Get("fieldmap");
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
      int ix = (int) lrint((x - X0) / DX), iy = (int) lrint((y - X0) / DX), iz = (int) lrint((z - Z0) / DX);
      if (ix < 0 || ix >= NX || iy < 0 || iy >= NY || iz < 0 || iz >= NZ) continue;
      size_t k = idx(ix, iy, iz);
      m_bx[k] = bx; m_by[k] = by; m_bz[k] = bz;
      used++;
    }
    f->Close();
    printf("MapField: %ld/%lld grid points loaded (%.0f MB)\n", used, n,
           3.0 * m_bx.size() * sizeof(float) / 1048576.);
  }
  void GetFieldValue(const G4double p[4], G4double *B) const override
  {
    // p in mm -> grid coords in cm
    double xc = p[0] / cm, yc = p[1] / cm, zc = p[2] / cm;
    double fx = (xc - X0) / DX, fy = (yc - X0) / DX, fz = (zc - Z0) / DX;
    int ix = (int) floor(fx), iy = (int) floor(fy), iz = (int) floor(fz);
    if (ix < 0 || ix >= NX - 1 || iy < 0 || iy >= NY - 1 || iz < 0 || iz >= NZ - 1)
    { B[0] = B[1] = B[2] = 0; return; }
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
};

// ---------------- HepMC2 ASCII reader (E/V/P lines, U GEV MM) ----------------
struct HepPart { int pdg; double px, py, pz, e; double vx, vy, vz, vt; };
class HepMC2Reader
{
 public:
  explicit HepMC2Reader(const std::string &fn) : m_in(fn)
  {
    if (!m_in) { fprintf(stderr, "hepmc open failed: %s\n", fn.c_str()); exit(1); }
  }
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
        if (inEvent) { m_in.seekg(lastpos); return true; }  // next event begins
        inEvent = true;
      }
      else if (tag == 'V' && inEvent)
      {
        // V barcode id x y z ctau n_orphan n_out [weights...]
        std::istringstream ss(line.substr(1));
        long bc; int id;
        ss >> bc >> id >> vx >> vy >> vz >> vt;
      }
      else if (tag == 'P' && inEvent)
      {
        // P barcode pdg px py pz e m status ...
        std::istringstream ss(line.substr(1));
        long bc; int pdg; double px, py, pz, e, m; int status;
        ss >> bc >> pdg >> px >> py >> pz >> e >> m >> status;
        if (status == 1) out.push_back({pdg, px, py, pz, e, vx, vy, vz, vt});
      }
      lastpos = m_in.tellg();
    }
    return inEvent;  // last event of file
  }
 private:
  std::ifstream m_in;
};

// ---------------- primary generator ----------------
class PrimGen : public G4VUserPrimaryGeneratorAction
{
 public:
  explicit PrimGen(HepMC2Reader *r) : m_r(r) {}
  void GeneratePrimaries(G4Event *ev) override
  {
    std::vector<HepPart> parts;
    if (!m_r->next(parts)) { fprintf(stderr, "HepMC exhausted\n"); return; }
    long placed = 0, skipped = 0;
    for (auto &p : parts)
    {
      G4ParticleDefinition *def = nullptr;
      if (p.pdg > 1000000000)
        def = G4IonTable::GetIonTable()->GetIon(p.pdg);
      else
        def = G4ParticleTable::GetParticleTable()->FindParticle(p.pdg);
      if (!def) { skipped++; continue; }
      // U GEV MM: lengths mm, time mm/c
      auto *vtx = new G4PrimaryVertex(p.vx * mm, p.vy * mm, p.vz * mm, p.vt * mm / CLHEP::c_light);
      auto *pp = new G4PrimaryParticle(def, p.px * GeV, p.py * GeV, p.pz * GeV);
      vtx->SetPrimary(pp);
      ev->AddPrimaryVertex(vtx);
      placed++;
    }
    printf("PrimGen: event %d placed %ld primaries (%ld unknown-pdg skipped)\n",
           ev->GetEventID(), placed, skipped);
  }
 private:
  HepMC2Reader *m_r;
};

// forward decl: ancestry map (defined below, used by the SD)
class TrkMap;
extern TrkMap *g_trkmap;
int trkmap_ancestor(int tid);

// ---------------- TPC gas sensitive detector ----------------
class TpcSD : public G4VSensitiveDetector
{
 public:
  TpcSD(const std::string &out) : G4VSensitiveDetector("TpcSD")
  {
    m_f = new TFile(out.c_str(), "RECREATE");
    // 15 branches: the 12 transport consumes + gflavor/gembed/gprimary which
    // frame_composer reads PER HIT ROW for truth attachment (old-eval layout).
    // gpx/gpy/gpz = TRACK VERTEX momentum (old-eval semantics: constant per
    // track, used for the gpt looper cut AND by transport for spreading -
    // bug-compatible with the frozen pipeline currency).
    m_nt = new TNtuple("ntp_g4hit", "g4hit => truth",
                       "event:glayer:gx:gy:gz:gt:gpl:gpx:gpy:gpz:gedep:gtrackID:gflavor:gembed:gprimary:gancestor");
  }
  void setEvent(int e) { m_ev = e; }
  TNtuple *makePrimTable()
  {
    m_f->cd();
    m_pt = new TNtuple("ntp_prim", "primary table", "event:ptrk:ppt:peta:ppdg");
    return m_pt;
  }
  G4bool ProcessHits(G4Step *st, G4TouchableHistory *) override
  {
    double edep = st->GetTotalEnergyDeposit();
    if (edep <= 0) return false;
    auto *pre = st->GetPreStepPoint();
    auto *post = st->GetPostStepPoint();
    double gx = 0.5 * (pre->GetPosition().x() + post->GetPosition().x()) / cm;
    double gy = 0.5 * (pre->GetPosition().y() + post->GetPosition().y()) / cm;
    double gz = 0.5 * (pre->GetPosition().z() + post->GetPosition().z()) / cm;
    double gt = pre->GetGlobalTime() / ns;
    double gpl = st->GetStepLength() / cm;
    G4Track *trk = st->GetTrack();
    auto *def = trk->GetParticleDefinition();
    double ke = trk->GetVertexKineticEnergy();
    double m0 = def->GetPDGMass();
    double pmag = std::sqrt(ke * (ke + 2 * m0));
    auto dir = trk->GetVertexMomentumDirection();
    int tid = trk->GetTrackID();
    bool prim = (trk->GetParentID() == 0);
    if (!prim) tid = -tid;  // secondary: negative (sPHENIX convention)
    int anc = trkmap_ancestor(trk->GetTrackID());
    float row[16] = {(float) m_ev, 99, (float) gx, (float) gy, (float) gz, (float) gt, (float) gpl,
                     (float) (dir.x() * pmag / GeV), (float) (dir.y() * pmag / GeV),
                     (float) (dir.z() * pmag / GeV), (float) (edep / GeV), (float) tid,
                     (float) def->GetPDGEncoding(), 0.f, prim ? 1.f : 0.f, (float) anc};
    m_nt->Fill(row);
    m_n++;
    return true;
  }
  void finish()
  {
    m_f->cd(); m_nt->Write();
    if (m_pt) m_pt->Write();
    m_f->Close();
    printf("TpcSD: wrote %ld hits total\n", m_n);
  }
 private:
  TFile *m_f; TNtuple *m_nt;
  TNtuple *m_pt = nullptr;
  int m_ev = 0; long m_n = 0;
};

// ---------------- detector construction: GDML + SD + field ----------------
class Det : public G4VUserDetectorConstruction
{
 public:
  Det(const std::string &gdml, const std::string &fmap, TpcSD *sd) : m_gdml(gdml), m_fmap(fmap), m_sd(sd) {}
  G4VPhysicalVolume *Construct() override
  {
    m_parser.Read(m_gdml, false);
    // CEMC effective fill (2026-07-15): the Spacal detector code hardcodes
    // gdml_config->exclude_physical_vol for its towers/fibers, so the export
    // carries only the air-filled CEMC_0 envelope. Fill it with an effective
    // average W-epoxy (92/8 by mass, rho = 8.0 g/cm3, DECLARED approximate,
    // gate-tuned) to restore the calorimeter albedo physics (late-time
    // back-scatter band + part of the secondary budget).
    auto *nist = G4NistManager::Instance();
    auto *eff = new G4Material("CEMC_eff", 8.0 * g / cm3, 2);
    eff->AddMaterial(nist->FindOrBuildMaterial("G4_W"), 0.92);
    eff->AddMaterial(nist->FindOrBuildMaterial("G4_POLYSTYRENE"), 0.08);
    int nfill = 0;
    for (auto *lv : *G4LogicalVolumeStore::GetInstance())
    {
      std::string n = lv->GetName();
      auto sp = n.find("0x");
      if (sp != std::string::npos) n = n.substr(0, sp);
      if (n == "CEMC_0") { lv->SetMaterial(eff); nfill++; }
    }
    printf("Det: CEMC_0 effective fill applied to %d volume(s)\n", nfill);
    return m_parser.GetWorldVolume();
  }
  void ConstructSDandField() override
  {
    // attach SD to every logical volume named tpc_gas* (GDML pointer suffix)
    int nvol = 0;
    for (auto *lv : *G4LogicalVolumeStore::GetInstance())
      if (strncmp(lv->GetName().c_str(), "tpc_gas", 7) == 0)
      {
        G4SDManager::GetSDMpointer()->AddNewDetector(m_sd);
        lv->SetSensitiveDetector(m_sd);
        // PHG4 TPC gas max step = 1 cm (measured: reference gpl distribution
        // has a hard wall at exactly 1.0 cm) - required for hit granularity
        lv->SetUserLimits(new G4UserLimits(1.0 * cm));
        nvol++;
      }
    printf("Det: SD attached to %d tpc_gas volume(s)\n", nvol);
    auto *fld = new MapField(m_fmap);
    auto *fm = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    fm->SetDetectorField(fld);
    fm->CreateChordFinder(fld);
  }
 private:
  std::string m_gdml, m_fmap;
  TpcSD *m_sd;
  G4GDMLParser m_parser;
};

// ---------------- ancestry (attribution map, 2026-07-15) ----------------
#include <G4UserTrackingAction.hh>
#include <unordered_map>
class TrkMap : public G4UserTrackingAction
{
 public:
  void PreUserTrackingAction(const G4Track *t) override
  {
    m_parent[t->GetTrackID()] = t->GetParentID();
    if (t->GetParentID() == 0)
    {
      auto p = t->GetMomentum();
      m_prim[t->GetTrackID()] = {p.perp() / GeV, p.eta(),
                                 (double) t->GetParticleDefinition()->GetPDGEncoding()};
    }
  }
  int ancestor(int tid) const
  {
    int guard = 0;
    while (guard++ < 10000)
    {
      auto it = m_parent.find(tid);
      if (it == m_parent.end() || it->second == 0) return tid;
      tid = it->second;
    }
    return tid;
  }
  void clear() { m_parent.clear(); m_prim.clear(); }
  std::unordered_map<int, int> m_parent;
  std::unordered_map<int, std::array<double, 3>> m_prim;  // id -> pt, eta, pdg
};
TrkMap *g_trkmap = nullptr;
int trkmap_ancestor(int tid) { return g_trkmap ? g_trkmap->ancestor(tid) : 0; }

// ---------------- event bookkeeping ----------------
#include <G4UserEventAction.hh>
class EvtAct : public G4UserEventAction
{
 public:
  explicit EvtAct(TpcSD *sd) : m_sd(sd) {}
  void BeginOfEventAction(const G4Event *ev) override
  {
    m_sd->setEvent(ev->GetEventID());
    if (g_trkmap) g_trkmap->clear();
  }
  void EndOfEventAction(const G4Event *ev) override
  {
    if (g_trkmap && m_prim)
      for (auto &kv : g_trkmap->m_prim)
        m_prim->Fill(ev->GetEventID(), kv.first, kv.second[0], kv.second[1], kv.second[2]);
    printf("EvtAct: event %d done\n", ev->GetEventID());
  }
  TNtuple *m_prim = nullptr;
 private:
  TpcSD *m_sd;
};

int main(int argc, char **argv)
{
  if (argc < 6)
  {
    fprintf(stderr, "usage: %s <gdml> <hepmc.dat> <out.root> <nevents> <fieldmap.root>\n", argv[0]);
    return 1;
  }
  std::string gdml = argv[1], hepmc = argv[2], out = argv[3], fmap = argv[5];
  int nev = atoi(argv[4]);

  auto *rm = new G4RunManager;
  auto *sd = new TpcSD(out);
  rm->SetUserInitialization(new Det(gdml, fmap, sd));
  auto *phys = new FTFP_BERT(0);
  phys->RegisterPhysics(new G4StepLimiterPhysics);
  rm->SetUserInitialization(phys);
  auto *reader = new HepMC2Reader(hepmc);
  rm->SetUserAction(new PrimGen(reader));
  g_trkmap = new TrkMap;
  rm->SetUserAction(g_trkmap);
  auto *ea = new EvtAct(sd);
  ea->m_prim = sd->makePrimTable();
  rm->SetUserAction(ea);
  rm->Initialize();
  rm->BeamOn(nev);
  sd->finish();
  delete rm;
  return 0;
}
