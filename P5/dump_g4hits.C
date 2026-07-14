// P5: DST -> ntp_g4hit (the 12-branch format tpc_transport consumes) +
// per-event truth-particle census (the secondary-flood measurement).
// Runs in the alma9 env (needs coresoftware libs for PHG4 containers).
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/SubsysReco.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <phool/getClass.h>
#include <TFile.h>
#include <TNtuple.h>
#include <cmath>
R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libg4dst.so)

class G4HitDumper : public SubsysReco
{
 public:
  G4HitDumper(const std::string &out) : SubsysReco("G4HitDumper"), m_out(out) {}
  int Init(PHCompositeNode *) override
  {
    m_f = new TFile(m_out.c_str(), "RECREATE");
    m_nt = new TNtuple("ntp_g4hit", "g4hit => truth",
                       "event:glayer:gx:gy:gz:gt:gpl:gpx:gpy:gpz:gedep:gtrackID");
    return 0;
  }
  int process_event(PHCompositeNode *top) override
  {
    auto *hits = findNode::getClass<PHG4HitContainer>(top, "G4HIT_TPC");
    auto *truth = findNode::getClass<PHG4TruthInfoContainer>(top, "G4TruthInfo");
    long nh = 0;
    if (hits)
    {
      for (auto it = hits->getHits().first; it != hits->getHits().second; ++it)
      {
        PHG4Hit *h = it->second;
        int trk = h->get_trkid();
        float px = 0, py = 0, pz = 0;
        if (truth)
        {
          PHG4Particle *p = truth->GetParticle(trk);
          if (p) { px = p->get_px(); py = p->get_py(); pz = p->get_pz(); }
        }
        float x = 0.5f * (h->get_x(0) + h->get_x(1));
        float y = 0.5f * (h->get_y(0) + h->get_y(1));
        float z = 0.5f * (h->get_z(0) + h->get_z(1));
        float dx = h->get_x(1) - h->get_x(0), dy = h->get_y(1) - h->get_y(0), dz = h->get_z(1) - h->get_z(0);
        float pl = std::sqrt(dx * dx + dy * dy + dz * dz);
        m_nt->Fill(m_ev, h->get_layer(), x, y, z, h->get_t(0), pl, px, py, pz, h->get_edep(), trk);
        nh++;
      }
    }
    long npart = 0, nsec = 0, nseclo = 0;
    if (truth)
    {
      for (auto pit = truth->GetMap().begin(); pit != truth->GetMap().end(); ++pit)
      {
        npart++;
        if (pit->first < 0 || !truth->is_primary(pit->second)) nsec++;
      }
    }
    printf("P5DUMP event %d: g4hits_tpc %ld | truth particles %ld (sec %ld)\n", m_ev, nh, npart, nsec);
    m_ev++;
    return 0;
  }
  int End(PHCompositeNode *) override
  {
    m_f->cd(); m_nt->Write(); m_f->Close();
    printf("P5DUMP wrote %s\n", m_out.c_str());
    return 0;
  }
 private:
  std::string m_out;
  TFile *m_f = nullptr;
  TNtuple *m_nt = nullptr;
  int m_ev = 0;
};

void dump_g4hits(const char *dst, const char *out = "P5_g4hit_eval.root")
{
  Fun4AllServer *se = Fun4AllServer::instance();
  se->registerSubsystem(new G4HitDumper(out));
  Fun4AllInputManager *in = new Fun4AllDstInputManager("DSTIN");
  in->fileopen(dst);
  se->registerInputManager(in);
  se->run(0);
  se->End();
  delete se;
}
