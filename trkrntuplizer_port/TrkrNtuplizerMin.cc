#include "TrkrNtuplizerMin.h"

#include <trackbase/ActsGeometry.h>
#include <trackbase/InttDefs.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackSeedContainer.h>

#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <cdbobjects/CDBTTree.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>
#include <phool/recoConsts.h>

#include <TFile.h>
#include <TNtuple.h>
#include <TVector3.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

// ---- branch index enums copied from upstream TrkrNtuplizer ----
namespace
{
enum n_event
{
  evnev,
  evnseed,
  evnrun,
  evnseg,
  evnjob,
  evsize = evnjob + 1  // 5
};

enum n_info
{
  infonocc11,
  infonocc116,
  infonocc21,
  infonocc216,
  infonocc31,
  infonocc316,
  infonrawzdc,
  infonlivezdc,
  infonscaledzdc,
  infonrawmbd,
  infonlivembd,
  infonscaledmbd,
  infonrawmbdv10,
  infonlivembdv10,
  infonscaledmbdv10,
  infonrawzdc1,
  infonlivezdc1,
  infonscaledzdc1,
  infonrawmbd1,
  infonlivembd1,
  infonscaledmbd1,
  infonrawmbdv101,
  infonlivembdv101,
  infonscaledmbdv101,
  inforzdc,
  informbd,
  informbdv10,
  infonbco1,
  infonbco,
  infonbcotr,
  infonbcotr1,
  infontrk,
  infonntpcseed,
  infonnsiseed,
  infonhitmvtx,
  infonhitintt,
  infonhittpot,
  infonhittpcall,
  infonhittpcin,
  infonhittpcmid,
  infonhittpcout,
  infonclusall,
  infonclustpc,
  infonclustpcpos,
  infonclustpcneg,
  infonclusintt,
  infonclusmvtx,
  infonclustpot,
  infosize = infonclustpot + 1  // 48
};

enum n_cluster
{
  nclulocx,
  nclulocy,
  nclux,
  ncluy,
  ncluz,
  nclur,
  ncluphi,
  nclueta,
  nclutheta,
  ncluphibin,
  nclutbin,
  nclufee,
  ncluchan,
  nclusampa,
  ncluex,
  ncluey,
  ncluez,
  ncluephi,
  nclupez,
  nclupephi,
  nclue,
  ncluadc,
  nclumaxadc,
  ncluthick,
  ncluafac,
  nclubfac,
  ncludcal,
  nclulayer,
  ncluphielem,
  ncluzelem,
  nclusize,
  ncluphisize,
  ncluzsize,
  nclupedge,
  ncluredge,
  ncluovlp,
  nclutrackID,
  ncluniter,
  clusize = ncluniter + 1  // 38
};

const char *DEFAULT_FEEMAP =
    "/home/rog/sPHENIX/3D_ClusterFindingML/CDB_offline/TPC_FEE_CHANNEL_MAP/"
    "8d/c7/8dc723e570afae13e631329f092db3ec_TPCChannelMapV4.root";
}  // namespace

TrkrNtuplizer::TrkrNtuplizer(const std::string & /*name*/, const std::string &filename,
                             const std::string &trackmapname,
                             unsigned int nlayers_maps, unsigned int nlayers_intt,
                             unsigned int nlayers_tpc, unsigned int nlayers_mms)
  : SubsysReco("TrkrNtuplizer")
  , _nlayers_maps(nlayers_maps)
  , _nlayers_intt(nlayers_intt)
  , _nlayers_tpc(nlayers_tpc)
  , _nlayers_mms(nlayers_mms)
  , _filename(filename)
  , _trackmapname(trackmapname)
{
}

TrkrNtuplizer::~TrkrNtuplizer() = default;

int TrkrNtuplizer::Init(PHCompositeNode * /*topNode*/)
{
  _ievent = 0;
  _tfile = new TFile(_filename.c_str(), "RECREATE");
  _tfile->SetCompressionLevel(7);

  // identical variable lists to upstream TrkrNtuplizer (real-data file format)
  std::string str_event = {"event:seed:run:seg:job"};
  std::string str_cluster = {"locx:locy:x:y:z:r:phi:eta:theta:phibin:tbin:fee:chan:sampa:ex:ey:ez:ephi:pez:pephi:e:adc:maxadc:thick:afac:bfac:dcal:layer:phielem:zelem:size:phisize:zsize:pedge:redge:ovlp:trackID:niter"};
  std::string str_info = {"occ11:occ116:occ21:occ216:occ31:occ316:rawzdc:livezdc:scaledzdc:rawmbd:livembd:scaledmbd:rawmbdv10:livembdv10:scaledmbdv10:rawzdc1:livezdc1:scaledzdc1:rawmbd1:livembd1:scaledmbd1:rawmbdv101:livembdv101:scaledmbdv101:rzdc:rmbd:rmbdv10:bco1:bco:bcotr:bcotr1:ntrk:ntpcseed:nsiseed:nhitmvtx:nhitintt:nhittpot:nhittpcall:nhittpcin:nhittpcmid:nhittpcout:nclusall:nclustpc:nclustpcpos:nclustpcneg:nclusintt:nclusmaps:nclusmms"};

  if (_do_cluster_eval)
  {
    std::string varlist = str_event + ":" + str_cluster + ":" + str_info;
    _ntp_cluster = new TNtuple("ntp_cluster", "svtxcluster (TrkrNtuplizer port)", varlist.c_str());
  }
  if (_do_info_eval)
  {
    std::string varlist = str_event + ":" + str_info;
    _ntp_info = new TNtuple("ntp_info", "event info (TrkrNtuplizer port)", varlist.c_str());
  }

  // GL1 emulation config: TRKRNTUP_GL1_EMUL=0 disables (columns stay 0);
  // rate from TRKRNTUP_GL1RATE, else PILEUPRATE, else 50 kHz.
  const char *emul = std::getenv("TRKRNTUP_GL1_EMUL");
  m_gl1_emul = !(emul && std::string(emul) == "0");
  if (const char *r = std::getenv("TRKRNTUP_GL1RATE"))
  {
    m_gl1_rate = atof(r);
  }
  else if (const char *p = std::getenv("PILEUPRATE"))
  {
    m_gl1_rate = atof(p);
  }

  std::cout << "TrkrNtuplizer(min port): Init -> booked"
            << (_ntp_cluster ? " ntp_cluster(91br)" : "")
            << (_ntp_info ? " ntp_info(53br)" : "")
            << (m_gl1_emul ? Form(" | GL1 EMULATION at %.3g Hz", m_gl1_rate) : " | GL1 cols 0")
            << " | per-event AutoSave ON" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int TrkrNtuplizer::InitRun(PHCompositeNode *topNode)
{
  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "CORRECTED_TRKR_CLUSTER");
  if (!_cluster_map)
  {
    _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  }
  if (!_cluster_map)
  {
    std::cout << PHWHERE << " ERROR: TrkrNtuplizer(min port) cannot find TRKR_CLUSTER" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  _tgeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!_tgeometry)
  {
    std::cout << PHWHERE << " ERROR: TrkrNtuplizer(min port) cannot find ActsGeometry" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  _geom_container = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!_geom_container)
  {
    std::cout << PHWHERE << " WARNING: no CYLINDERCELLGEOM_SVTX -> phibin/tbin/fee/pedge/occ stay NaN/0" << std::endl;
  }
  _cluster_hit_assoc = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  _hitsets = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!_cluster_hit_assoc || !_hitsets)
  {
    std::cout << PHWHERE << " WARNING: no TRKR_CLUSTERHITASSOC/TRKR_HITSET -> maxadc/size stay NaN" << std::endl;
  }

  BuildFeeMap();
  return Fun4AllReturnCodes::EVENT_OK;
}

void TrkrNtuplizer::BuildFeeMap()
{
  // Offline replica of upstream TrkrNtuplizer::InitRun FEE mapping: reads the
  // TPC_FEE_CHANNEL_MAP payload straight from CDB_offline (no CDB/network),
  // then maps (layer, mc-sector, side, pad-bin) -> (fee, channel, sampa).
  if (!_geom_container)
  {
    return;
  }
  const char *env = std::getenv("TRKRNTUP_FEEMAP");
  std::string path = env ? env : DEFAULT_FEEMAP;
  CDBTTree *cdbttree = new CDBTTree(path);
  cdbttree->LoadCalibrations();

  long built = 0;
  for (unsigned int sector = 0; sector < 24; sector++)
  {
    int side = (sector > 11) ? 0 : 1;
    for (unsigned int fee = 0; fee <= 25; fee++)
    {
      for (unsigned int channel = 0; channel <= 255; channel++)
      {
        int feeM = FEE_map[fee];
        if (FEE_R[fee] == 2)
        {
          feeM += 6;
        }
        if (FEE_R[fee] == 3)
        {
          feeM += 14;
        }
        unsigned int key = (256 * feeM) + channel;
        int layer = cdbttree->GetIntValue(key, "layer", 0);
        if (layer <= 6)  // antenna pads sit in layer 0
        {
          continue;
        }
        double phi = ((side == 1 ? 1 : -1) * (cdbttree->GetDoubleValue(key, "phi", 0) - M_PI / 2.)) +
                     ((sector % 12) * M_PI / 6);
        PHG4TpcCylinderGeom *layergeom = _geom_container->GetLayerCellGeom(layer);
        if (!layergeom)
        {
          continue;
        }
        // ana.331 get_phibin has no side argument (newer builds use get_pad_float(phi, side))
        int phibin = layergeom->get_phibin(phi);
        if (phibin < 0)
        {
          continue;
        }
        fee_info fi;
        fi.fee = fee;
        fi.channel = channel;
        fi.sampa = channel % 32;
        // phibin is GLOBAL (encodes the sector already); keying additionally on the
        // mc-sector caused ~11/12 lookup misses when the ana.331 phi convention and
        // the upstream sector transform disagree -> key on (layer, side, phibin) only.
        _fee_map[feeKey(layer, 0, side, phibin)] = fi;
        built++;
      }
    }
  }
  delete cdbttree;
  _fee_map_ok = !_fee_map.empty();
  std::cout << "TrkrNtuplizer(min port): FEE map from " << path
            << " -> " << built << " pad entries (" << (_fee_map_ok ? "OK" : "EMPTY") << ")" << std::endl;
}

void TrkrNtuplizer::FillInfo(PHCompositeNode *topNode, float *fx_info)
{
  // per-layer hit tallies + occupancies (upstream fillOutputNtuples logic)
  float nhit[100];
  std::fill(nhit, nhit + 100, 0.f);
  if (_hitsets)
  {
    auto all_hitsets = _hitsets->getHitSets();
    for (auto hsiter = all_hitsets.first; hsiter != all_hitsets.second; ++hsiter)
    {
      unsigned int layer = TrkrDefs::getLayer(hsiter->first);
      auto hitrange = hsiter->second->getHits();
      for (auto h = hitrange.first; h != hitrange.second; ++h)
      {
        if (layer < 100)
        {
          nhit[layer]++;
        }
        if (layer < _nlayers_maps)
        {
          fx_info[infonhitmvtx]++;
        }
        if (layer >= _nlayers_maps && layer < _nlayers_maps + _nlayers_intt)
        {
          fx_info[infonhitintt]++;
        }
        if (layer >= _nlayers_maps + _nlayers_intt &&
            layer < _nlayers_maps + _nlayers_intt + _nlayers_tpc)
        {
          fx_info[infonhittpcall]++;
        }
        if (layer >= _nlayers_maps + _nlayers_intt + _nlayers_tpc)
        {
          fx_info[infonhittpot]++;
        }
        if (layer == _nlayers_maps + _nlayers_intt)
        {
          fx_info[infonhittpcin]++;
        }
        if (layer == _nlayers_maps + _nlayers_intt + _nlayers_tpc - 1)
        {
          fx_info[infonhittpcout]++;
        }
        if (layer == _nlayers_maps + _nlayers_intt + _nlayers_tpc / 2 - 1)
        {
          fx_info[infonhittpcmid]++;
        }
      }
    }
  }
  if (_geom_container)
  {
    const int occ_idx[6] = {infonocc11, infonocc116, infonocc21, infonocc216, infonocc31, infonocc316};
    const int occ_lay[6] = {0, 15, 16, 31, 32, 47};
    for (int i = 0; i < 6; ++i)
    {
      int layer = (int) (_nlayers_maps + _nlayers_intt) + occ_lay[i];
      PHG4TpcCylinderGeom *g = _geom_container->GetLayerCellGeom(layer);
      if (g && g->get_phibins() > 0 && g->get_zbins() > 0)
      {
        fx_info[occ_idx[i]] = nhit[layer] / ((float) g->get_phibins() * g->get_zbins());
      }
    }
  }
  // cluster tallies
  if (_cluster_map)
  {
    fx_info[infonclusall] = _cluster_map->size();
    for (const auto &hitsetkey : _cluster_map->getHitSetKeys())
    {
      auto range = _cluster_map->getClusters(hitsetkey);
      for (auto iter = range.first; iter != range.second; ++iter)
      {
        unsigned int layer = TrkrDefs::getLayer(iter->first);
        if (layer < _nlayers_maps)
        {
          fx_info[infonclusmvtx]++;
        }
        else if (layer < _nlayers_maps + _nlayers_intt)
        {
          fx_info[infonclusintt]++;
        }
        else if (layer < _nlayers_maps + _nlayers_intt + _nlayers_tpc)
        {
          fx_info[infonclustpc]++;
          if (TpcDefs::getSide(iter->first) == 0)
          {
            fx_info[infonclustpcneg]++;
          }
          else
          {
            fx_info[infonclustpcpos]++;
          }
        }
        else
        {
          fx_info[infonclustpot]++;
        }
      }
    }
  }
  // track / seed container counts (0 tracks locally is the true value)
  auto *trackmap = findNode::getClass<SvtxTrackMap>(topNode, _trackmapname);
  if (trackmap)
  {
    fx_info[infontrk] = (float) trackmap->size();
  }
  auto *tpcseedmap = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  if (tpcseedmap)
  {
    fx_info[infonntpcseed] = (float) tpcseedmap->size();
  }
  auto *siseedmap = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
  if (siseedmap)
  {
    fx_info[infonnsiseed] = (float) siseedmap->size();
  }
  // GL1 scalers / BCO / rates: real-data-only, remain 0 in simulation
}

int TrkrNtuplizer::process_event(PHCompositeNode *topNode)
{
  recoConsts *rc = recoConsts::instance();
  if (rc->FlagExist("RANDOMSEED"))
  {
    m_fSeed = rc->get_IntFlag("RANDOMSEED");
  }
  else
  {
    m_fSeed = std::numeric_limits<float>::quiet_NaN();
  }

  float fx_event[evsize] = {(float) _ievent, m_fSeed, (float) m_runnumber,
                            (float) m_segment, (float) m_job};
  float fx_info[infosize] = {0};
  FillInfo(topNode, fx_info);

  // --- GL1 EMULATION (documented: synthesized, not measured — sim has no trigger
  //     hardware). Models a Poisson min-bias trigger sequence at m_gl1_rate:
  //     BCO advances by Exp(1/rate) in 106 ns RHIC-clock ticks per triggered event;
  //     raw/live scalers count triggers (livetime 1.0, prescale off -> scaled=0);
  //     mbdv10 = |z_vtx|<10 cm subset for the sigma_z=7 cm vertex profile (P=0.847);
  //     rate columns carry the configured rate. Reproducible (fixed TRandom3 seed).
  if (m_gl1_emul)
  {
    double dt = m_gl1_rng.Exp(1.0 / m_gl1_rate);  // seconds since previous trigger
    m_bco += (uint64_t) (dt / 106e-9);
    m_rawmbd += 1;
    m_rawzdc += 1;  // AuAu: ZDC coincidence fires with ~every hadronic MB event
    if (std::fabs(m_gl1_rng.Gaus(0., 7.0)) < 10.)
    {
      m_rawmbdv10 += 1;
    }
    if (_ievent == 0)
    {
      m_bco1 = m_bco;
      m_rawmbd1 = m_rawmbd;
      m_rawzdc1 = m_rawzdc;
      m_rawmbdv101 = m_rawmbdv10;
    }
    uint64_t bcotr = (m_bco << 24U) >> 24U;
    uint64_t bcotr1 = (m_bco1 << 24U) >> 24U;
    fx_info[infonbco] = m_bco - m_bco1;
    fx_info[infonbcotr] = bcotr - bcotr1;
    fx_info[infonbco1] = m_bco1;
    fx_info[infonbcotr1] = bcotr1;
    fx_info[infonrawmbd] = m_rawmbd;
    fx_info[infonlivembd] = m_rawmbd;
    fx_info[infonrawzdc] = m_rawzdc;
    fx_info[infonlivezdc] = m_rawzdc;
    fx_info[infonrawmbdv10] = m_rawmbdv10;
    fx_info[infonlivembdv10] = m_rawmbdv10;
    fx_info[infonrawmbd1] = m_rawmbd1;
    fx_info[infonlivembd1] = m_rawmbd1;
    fx_info[infonrawzdc1] = m_rawzdc1;
    fx_info[infonlivezdc1] = m_rawzdc1;
    fx_info[infonrawmbdv101] = m_rawmbdv101;
    fx_info[infonlivembdv101] = m_rawmbdv101;
    fx_info[inforzdc] = m_gl1_rate;
    fx_info[informbd] = m_gl1_rate;
    fx_info[informbdv10] = m_gl1_rate * 0.847;
  }

  if (_ntp_info)
  {
    float row[evsize + infosize];
    std::copy(fx_event, fx_event + evsize, row);
    std::copy(fx_info, fx_info + infosize, row + evsize);
    _ntp_info->Fill(row);
  }

  if (_ntp_cluster && _cluster_map && _tgeometry)
  {
    const int NC = (int) evsize + (int) clusize + (int) infosize;
    float row[NC];
    for (const auto &hitsetkey : _cluster_map->getHitSetKeys())
    {
      auto range = _cluster_map->getClusters(hitsetkey);
      for (auto iter = range.first; iter != range.second; ++iter)
      {
        TrkrDefs::cluskey cluster_key = iter->first;
        float fxc[clusize] = {0};
        FillCluster(&fxc[0], cluster_key);
        std::copy(fx_event, fx_event + evsize, row);
        std::copy(fxc, fxc + (int) clusize, row + evsize);
        std::copy(fx_info, fx_info + infosize, row + evsize + (int) clusize);
        _ntp_cluster->Fill(row);
      }
    }
  }

  // kill-safe output: commit baskets + tree metadata every event, so a run killed at
  // any point keeps all previously completed events readable (no End() dependence)
  if (_tfile)
  {
    _tfile->cd();
    if (_ntp_cluster)
    {
      _ntp_cluster->AutoSave("SaveSelf,Overwrite");
    }
    if (_ntp_info)
    {
      _ntp_info->AutoSave("SaveSelf,Overwrite");
    }
  }

  ++_ievent;
  return Fun4AllReturnCodes::EVENT_OK;
}

void TrkrNtuplizer::FillCluster(float *fX, TrkrDefs::cluskey cluster_key)
{
  const float NaN = std::numeric_limits<float>::quiet_NaN();
  for (int i = 0; i < (int) clusize; ++i)
  {
    fX[i] = NaN;
  }

  unsigned int layer = TrkrDefs::getLayer(cluster_key);
  TrkrCluster *cluster = _cluster_map->findCluster(cluster_key);
  if (!cluster)
  {
    return;
  }

  Acts::Vector3 cglob = _tgeometry->getGlobalPosition(cluster_key, cluster);
  TVector3 pos(cglob(0), cglob(1), cglob(2));
  float r = pos.Perp();
  float phi = pos.Phi();
  float locx = cluster->getLocalX();
  float locy = cluster->getLocalY();

  fX[nclulocx] = locx;
  fX[nclulocy] = locy;
  fX[nclux] = cglob(0);
  fX[ncluy] = cglob(1);
  fX[ncluz] = cglob(2);
  fX[nclur] = r;
  fX[ncluphi] = phi;
  fX[nclueta] = pos.Eta();
  fX[nclutheta] = pos.Theta();
  fX[ncludcal] = 1;
  fX[nclue] = (float) cluster->getAdc();
  fX[ncluadc] = (float) cluster->getAdc();
  fX[nclulayer] = layer;
  fX[ncluphisize] = cluster->getPhiSize();
  fX[ncluzsize] = cluster->getZSize();
  fX[ncluovlp] = 3;
  fX[ncluniter] = 0;
  fX[ncluredge] = (layer == 7 || layer == 22 || layer == 23 ||
                   layer == 28 || layer == 39 || layer == 54)
                      ? 1.f
                      : 0.f;

  // cluster errors: default = ana.331-era parameterization; for TPC overwritten
  // below with a v5-style reconstruction (base error from charge-weighted hit RMS
  // as in the newer TpcClusterizer + the v5 'data' multipliers from upstream
  // ClusterErrorPara::get_clusterv5_modified_error).
  auto para_errors = _ClusErrPara.get_simple_cluster_error(cluster, r, cluster_key);
  fX[ncluez] = std::sqrt(para_errors.second);
  fX[ncluephi] = std::sqrt(para_errors.first);

  // maxadc + size (exact) and, for the TPC, charge-weighted pad/time statistics
  // for the v5-style error base — all from the cluster's own hits.
  double wsum = 0, pmean = 0, pvar = 0, tmean = 0, tvar = 0;
  int padlo = 1 << 30, padhi = -1, tlo = 1 << 30, thi = -1;
  bool isTpc = (layer >= _nlayers_maps + _nlayers_intt &&
                layer < _nlayers_maps + _nlayers_intt + _nlayers_tpc);
  if (_cluster_hit_assoc && _hitsets)
  {
    TrkrDefs::hitsetkey hskey = TrkrDefs::getHitSetKeyFromClusKey(cluster_key);
    TrkrHitSet *hitset = _hitsets->findHitSet(hskey);
    if (hitset)
    {
      auto hits = _cluster_hit_assoc->getHits(cluster_key);
      int nh = 0;
      float maxadc = 0;
      for (auto h = hits.first; h != hits.second; ++h)
      {
        TrkrHit *hit = hitset->getHit(h->second);
        if (hit)
        {
          nh++;
          float w = (float) hit->getAdc();
          maxadc = std::max(maxadc, w);
          if (isTpc)
          {
            int pad = TpcDefs::getPad(h->second);
            int tb = TpcDefs::getTBin(h->second);
            wsum += w;
            pmean += w * pad;
            pvar += w * pad * pad;
            tmean += w * tb;
            tvar += w * tb * tb;
            padlo = std::min(padlo, pad);
            padhi = std::max(padhi, pad);
            tlo = std::min(tlo, tb);
            thi = std::max(thi, tb);
          }
        }
      }
      if (nh > 0)
      {
        fX[nclusize] = nh;
        fX[nclumaxadc] = maxadc;
      }
    }
  }

  if (layer >= _nlayers_maps + _nlayers_intt &&
      layer < _nlayers_maps + _nlayers_intt + _nlayers_tpc)
  {
    // TPC
    unsigned int side = TpcDefs::getSide(cluster_key);
    unsigned int sector = TpcDefs::getSectorId(cluster_key);
    fX[ncluphielem] = sector;
    fX[ncluzelem] = side;

    if (_geom_container)
    {
      PHG4TpcCylinderGeom *g = _geom_container->GetLayerCellGeom(layer);
      if (g)
      {
        // --- v5-style ez/ephi: base error a la newer TpcClusterizer
        //     (phi_err^2 = r^2 * var_phi/(adc*0.14); single-pad -> 9*(r*pitch)^2/12),
        //     then the upstream get_clusterv5_modified_error 'data' multipliers.
        //     Edge/overlap multipliers are skipped (v4 clusters carry neither).
        if (wsum > 0)
        {
          double phistep = g->get_phistep();
          double clock = g->get_zstep();  // ns (ADC clock period in ana.331)
          double vd = _tgeometry->get_drift_velocity();
          pmean /= wsum;
          tmean /= wsum;
          double varp = pvar / wsum - pmean * pmean;
          double vart = tvar / wsum - tmean * tmean;
          double phierr2 = (padhi == padlo)
                               ? 9. * (r * phistep) * (r * phistep) / 12.
                               : r * r * (varp * phistep * phistep) / (wsum * 0.14);
          double terr2 = (thi == tlo)
                             ? 9. * clock * clock / 12.
                             : (vart * clock * clock) / (wsum * 0.14);
          double phierror = std::sqrt(phierr2);
          double zerror = std::sqrt(terr2) * vd;
          if (layer == 7 || layer == 22 || layer == 23 || layer == 38 || layer == 39 || layer == 54)
          {
            phierror *= 4;
            zerror *= 4;
          }
          int psz = (int) cluster->getPhiSize();
          int zsz = (int) cluster->getZSize();
          if (psz == 2) { phierror *= 3.15; }
          else if (psz == 3) { phierror *= 3.5; }
          else if (psz > 3) { phierror *= 4; }
          int reg = (layer < 7 + 16) ? 0 : (layer < 7 + 32 ? 1 : 2);
          static const double zf2[3] = {7., 4.5, 4.5};
          static const double zf34[3] = {7., 5., 5.};
          static const double zf5[3] = {20., 6., 7.};
          if (zsz == 2) { zerror *= zf2[reg]; }
          else if (zsz == 3 || zsz == 4) { zerror *= zf34[reg]; }
          else if (zsz >= 5) { zerror *= zf5[reg]; }
          auto pol2 = [](double x, double p0, double p1, double p2) { return p0 + p1 * x + p2 * x * x; };
          if (reg == 0) { phierror *= pol2(layer, 3.206, -0.252, 0.007); }
          if (reg == 1) { phierror *= pol2(layer, 4.48, -0.226, 0.00362); zerror *= pol2(layer, 5.593, -0.2458, 0.00333455); }
          if (reg == 2) { phierror *= pol2(layer, 14.8112, -0.577, 0.00605); zerror *= pol2(layer, 5.6964, -0.21338, 0.002502); }
          if (psz >= 5) { phierror *= 10; }
          fX[ncluephi] = phierror;
          fX[ncluez] = zerror;
        }

        int phibin = g->get_phibin(phi);
        fX[ncluphibin] = phibin;
        // tbin = drift time in ADC clock ticks. In ana.331 get_zstep() IS the ADC
        // clock period in ns (upstream uses it as AdcClockPeriod), so ticks =
        // drift_length_cm / (v_drift * clock). Sim t0=0; real data adds a trigger offset.
        double clock = g->get_zstep();
        double vd = _tgeometry->get_drift_velocity();
        if (clock > 0 && vd > 0)
        {
          // physical half drift length is 105.5 cm; do NOT derive it from get_zbins()
          // (the extended-readout config inflates the padplane time axis / zbins)
          const double halfz = 105.5;
          fX[nclutbin] = (halfz - std::fabs(cglob(2))) / (clock * vd);
        }
        // pedge: cluster pad range touches the sector boundary (upstream: clusterizer edge flag)
        int nper = g->get_phibins() / 12;
        if (nper > 0 && phibin >= 0)
        {
          int local = phibin % nper;
          float half = 0.5f * cluster->getPhiSize();
          fX[nclupedge] = (local < half || local > nper - 1 - half) ? 1.f : 0.f;
        }
        if (_fee_map_ok && phibin >= 0)
        {
          auto it = _fee_map.find(feeKey(layer, 0, side, phibin));
          if (it != _fee_map.end())
          {
            fX[nclufee] = it->second.fee;
            fX[ncluchan] = it->second.channel;
            fX[nclusampa] = it->second.sampa;
          }
          else
          {
            fX[nclufee] = 0;
            fX[ncluchan] = 0;
          }
        }
      }
    }
  }
  else
  {
    // non-TPC: upstream convention phibin=locx, tbin=locy
    fX[ncluphibin] = locx;
    fX[nclutbin] = locy;
    if (layer < 3)
    {
      fX[ncluphielem] = MvtxDefs::getStaveId(cluster_key);
      fX[ncluzelem] = MvtxDefs::getChipId(cluster_key);
    }
    else if (layer < 7)
    {
      fX[ncluphielem] = InttDefs::getLadderPhiId(cluster_key);
      fX[ncluzelem] = InttDefs::getLadderZId(cluster_key);
    }
  }
}

int TrkrNtuplizer::End(PHCompositeNode * /*topNode*/)
{
  if (_tfile)
  {
    _tfile->cd();
    _tfile->Write();
    _tfile->Close();
  }
  std::cout << "TrkrNtuplizer(min port): End -> wrote " << _filename
            << " (processed " << _ievent << " events)" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
