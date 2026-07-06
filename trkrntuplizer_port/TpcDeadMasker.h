#ifndef TPCDEADMASKER_H
#define TPCDEADMASKER_H

// Masks dead TPC channels in simulation: removes TRKR_HITSET hits whose
// (layer, sector, side, pad) appears in a real-data TPC_DEADCHANNELMAP payload
// (CDB_offline; TTree "Multiple" with branches Ilayer/Isector/Iside/Ipad).
// Register BEFORE the clusterizer so sim clustering sees the same dead regions
// as real data. Fully offline - the payload is read straight from the file.

#include <fun4all/SubsysReco.h>

#include <cstdint>
#include <string>
#include <unordered_map>
#include <unordered_set>

class PHCompositeNode;

class TpcDeadChannelMasker : public SubsysReco
{
 public:
  explicit TpcDeadChannelMasker(const std::string &payload)
    : SubsysReco("TpcDeadChannelMasker")
    , m_payload(payload)
  {
  }
  ~TpcDeadChannelMasker() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

 private:
  static uint32_t key(int layer, int sector, int side)
  {
    return ((uint32_t) layer << 16U) | ((uint32_t) sector << 8U) | (uint32_t) side;
  }
  std::string m_payload;
  std::unordered_map<uint32_t, std::unordered_set<int> > m_dead;  // (layer,sector,side) -> dead pads
  long m_nremoved{0};
  long m_nseen{0};
};

#endif
