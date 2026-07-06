#include "TpcDeadMasker.h"

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include <vector>

int TpcDeadChannelMasker::InitRun(PHCompositeNode * /*topNode*/)
{
  TFile *f = TFile::Open(m_payload.c_str());
  if (!f || f->IsZombie())
  {
    std::cout << PHWHERE << " ERROR: cannot open dead map payload " << m_payload << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  TTree *t = (TTree *) f->Get("Multiple");
  if (!t)
  {
    std::cout << PHWHERE << " ERROR: no 'Multiple' tree in " << m_payload << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  Int_t layer = 0, pad = 0, sector = 0, side = 0;
  t->SetBranchAddress("Ilayer", &layer);
  t->SetBranchAddress("Ipad", &pad);
  t->SetBranchAddress("Isector", &sector);
  t->SetBranchAddress("Iside", &side);
  long n = t->GetEntries();
  for (long i = 0; i < n; ++i)
  {
    t->GetEntry(i);
    m_dead[key(layer, sector, side)].insert(pad);
  }
  f->Close();
  std::cout << "TpcDeadChannelMasker: loaded " << n << " dead channels ("
            << m_dead.size() << " (layer,sector,side) groups) from " << m_payload << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcDeadChannelMasker::process_event(PHCompositeNode *topNode)
{
  TrkrHitSetContainer *hitsets = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!hitsets)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }
  auto range = hitsets->getHitSets(TrkrDefs::TrkrId::tpcId);
  for (auto hs = range.first; hs != range.second; ++hs)
  {
    unsigned int layer = TrkrDefs::getLayer(hs->first);
    unsigned int sector = TpcDefs::getSectorId(hs->first);
    unsigned int side = TpcDefs::getSide(hs->first);
    auto dead = m_dead.find(key(layer, sector, side));
    if (dead == m_dead.end())
    {
      continue;
    }
    std::vector<TrkrDefs::hitkey> to_remove;
    auto hits = hs->second->getHits();
    for (auto h = hits.first; h != hits.second; ++h)
    {
      m_nseen++;
      if (dead->second.count(TpcDefs::getPad(h->first)))
      {
        to_remove.push_back(h->first);
      }
    }
    for (auto k : to_remove)
    {
      hs->second->removeHit(k);
    }
    m_nremoved += to_remove.size();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcDeadChannelMasker::End(PHCompositeNode * /*topNode*/)
{
  std::cout << "TpcDeadChannelMasker: removed " << m_nremoved
            << " hits on dead channels (of " << m_nseen << " TPC hits inspected)" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
