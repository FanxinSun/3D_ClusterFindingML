#ifndef TRKRNTUPLIZERMIN_H
#define TRKRNTUPLIZERMIN_H

// Cluster+info port of coresoftware TrkrNtuplizer (TrackingDiagnostics), back-ported
// to build & run against ana.331 (gcc-8.3) fully OFFLINE (no BNL network):
//  - ntp_cluster: real-data 91-branch layout. Physics columns + phibin/tbin (pad geom),
//    fee/chan/sampa (TPC_FEE_CHANNEL_MAP read directly from CDB_offline via CDBTTree),
//    ez/ephi (ana.331 ClusterErrorPara::get_simple_cluster_error), maxadc/size (exact,
//    from TRKR_CLUSTERHITASSOC hits), pedge (sector-edge estimate from pad position).
//  - ntp_info: real-data 53-branch layout; occupancy + nhit/nclus tallies + seed/track
//    container counts computed in-job. GL1/BCO/scaler columns stay 0 (real-data only).
// Track / clus_trk / hit / vertex paths of the upstream tool are omitted.
// Env: TRKRNTUP_FEEMAP overrides the FEE-map payload path.

#include <fun4all/SubsysReco.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase/TrkrDefs.h>

#include <TRandom3.h>

#include <cstdint>
#include <map>
#include <string>

class PHCompositeNode;
class TFile;
class TNtuple;
class TrkrClusterContainer;
class TrkrClusterHitAssoc;
class TrkrHitSetContainer;
class ActsGeometry;
class PHG4TpcCylinderGeomContainer;

class TrkrNtuplizer : public SubsysReco
{
 public:
  TrkrNtuplizer(const std::string &name = "TrkrNtuplizer",
                const std::string &filename = "trkrntuple.root",
                const std::string &trackmapname = "SvtxTrackMap",
                unsigned int nlayers_maps = 3,
                unsigned int nlayers_intt = 4,
                unsigned int nlayers_tpc = 48,
                unsigned int nlayers_mms = 2);
  ~TrkrNtuplizer() override;

  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void do_info_eval(bool b) { _do_info_eval = b; }
  void do_vertex_eval(bool b) { _do_vertex_eval = b; }
  void do_hit_eval(bool b) { _do_hit_eval = b; }
  void do_cluster_eval(bool b) { _do_cluster_eval = b; }
  void do_clus_trk_eval(bool b) { _do_clus_trk_eval = b; }
  void do_track_eval(bool b) { _do_track_eval = b; }
  void do_tpcseed_eval(bool b) { _do_tpcseed_eval = b; }
  void do_siseed_eval(bool b) { _do_siseed_eval = b; }
  void set_first_event(int value) { _ievent = value; }
  void segment(const int seg) { m_segment = seg; }
  void runnumber(const int run) { m_runnumber = run; }
  void job(const int j) { m_job = j; }

 private:
  void FillCluster(float *fXcluster, TrkrDefs::cluskey cluster_key);
  void FillInfo(PHCompositeNode *topNode, float *fx_info);
  void BuildFeeMap();
  static uint64_t feeKey(unsigned int layer, unsigned int sector,
                         unsigned int side, unsigned int phibin)
  {
    return ((uint64_t) layer << 40U) | ((uint64_t) sector << 32U) |
           ((uint64_t) side << 24U) | phibin;
  }

  bool _do_info_eval{true};
  bool _do_vertex_eval{false};
  bool _do_hit_eval{false};
  bool _do_cluster_eval{true};
  bool _do_clus_trk_eval{false};
  bool _do_track_eval{false};
  bool _do_tpcseed_eval{false};
  bool _do_siseed_eval{false};

  unsigned int _nlayers_maps{3};
  unsigned int _nlayers_intt{4};
  unsigned int _nlayers_tpc{48};
  unsigned int _nlayers_mms{2};

  int m_segment{0};
  int m_runnumber{0};
  int m_job{0};
  unsigned int _ievent{0};
  float m_fSeed{0};

  TNtuple *_ntp_cluster{nullptr};
  TNtuple *_ntp_info{nullptr};
  std::string _filename;
  std::string _trackmapname;
  ClusterErrorPara _ClusErrPara;
  TrkrClusterContainer *_cluster_map{nullptr};
  TrkrClusterHitAssoc *_cluster_hit_assoc{nullptr};
  TrkrHitSetContainer *_hitsets{nullptr};
  ActsGeometry *_tgeometry{nullptr};
  PHG4TpcCylinderGeomContainer *_geom_container{nullptr};
  TFile *_tfile{nullptr};

  struct fee_info
  {
    unsigned short fee{0};
    unsigned short channel{0};
    unsigned short sampa{0};
  };
  std::map<uint64_t, fee_info> _fee_map;
  bool _fee_map_ok{false};

  // GL1 emulation (sim has no trigger hardware): synthetic Poisson MB trigger sequence
  bool m_gl1_emul{true};
  double m_gl1_rate{50000.};  // Hz; TRKRNTUP_GL1RATE > PILEUPRATE > 50 kHz
  TRandom3 m_gl1_rng{20260705};
  uint64_t m_bco{0};
  uint64_t m_bco1{0};
  double m_rawmbd{0}, m_rawzdc{0}, m_rawmbdv10{0};
  double m_rawmbd1{0}, m_rawzdc1{0}, m_rawmbdv101{0};

  // upstream TrkrNtuplizer FEE cabling tables (coresoftware TrackingDiagnostics)
  int mc_sectors[12] = {5, 4, 3, 2, 1, 0, 11, 10, 9, 8, 7, 6};
  int FEE_map[26] = {4, 5, 0, 2, 1, 11, 9, 10, 8, 7, 6, 0, 1, 3, 7, 6, 5, 4, 3, 2, 0, 2, 1, 3, 5, 4};
  int FEE_R[26] = {2, 2, 1, 1, 1, 3, 3, 3, 3, 3, 3, 2, 2, 1, 2, 2, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3};
};

#endif
