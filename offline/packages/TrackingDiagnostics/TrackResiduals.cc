
#include "TrackResiduals.h"

#include <trackbase/ClusterErrorPara.h>
#include <trackbase/InttDefs.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <globalvertex/MbdVertex.h>
#include <globalvertex/MbdVertexMap.h>

#include <intt/CylinderGeomIntt.h>
#include <intt/CylinderGeomInttHelper.h>

#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtHit.h>

#include <micromegas/CylinderGeomMicromegas.h>
#include <micromegas/MicromegasDefs.h>

#include <mvtx/CylinderGeom_Mvtx.h>
#include <mvtx/CylinderGeom_MvtxHelper.h>

#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/SvtxAlignmentState.h>
#include <trackbase_historic/SvtxAlignmentStateMap.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackAnalysisUtils.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeedHelper.h>

#include <tpc/LaserEventInfo.h>

#include <ffarawobjects/Gl1Packet.h>
#include <ffarawobjects/Gl1RawHit.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>

#include <TSystem.h>

#include <cmath>
#include <algorithm>
#include <limits>

namespace
{
  template <class T>
  inline T square(const T& t)
  {
    return t * t;
  }
  template <class T>
  inline T r(const T& x, const T& y)
  {
    return std::sqrt(square(x) + square(y));
  }

  std::vector<TrkrDefs::cluskey> get_cluster_keys(SvtxTrack* track)
  {
    std::vector<TrkrDefs::cluskey> out;

    if(!track)
      {
	return out;
      }
    
    for (const auto& seed : {track->get_silicon_seed(), track->get_tpc_seed()})
    {
      if (seed)
      {
        std::copy(seed->begin_cluster_keys(), seed->end_cluster_keys(), std::back_inserter(out));
      }
    }

    return out;
  }
}  // namespace

//____________________________________________________________________________..
TrackResiduals::TrackResiduals(const std::string& name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int TrackResiduals::InitRun(PHCompositeNode* topNode)
{
  m_outfile = new TFile(m_outfileName.c_str(), "RECREATE");
  createBranches();

  // global position wrapper
  m_globalPositionWrapper.loadNodes(topNode);
  m_globalPositionWrapper.set_suppressCrossing(m_convertSeeds);
  // clusterMover needs the correct radii of the TPC layers
  auto *tpccellgeo = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  m_clusterMover.initialize_geometry(tpccellgeo);
  m_clusterMover.set_verbosity(0);

  auto *se = Fun4AllServer::instance();
  m_runnumber = se->RunNumber();

  return Fun4AllReturnCodes::EVENT_OK;
}
void TrackResiduals::clearClusterStateVectors()
{
  m_cluskeys.clear();
  m_clusphisize.clear();
  m_cluszsize.clear();
  m_idealsurfcenterx.clear();
  m_idealsurfcentery.clear();
  m_idealsurfcenterz.clear();
  m_idealsurfnormx.clear();
  m_clsector.clear();
  m_clside.clear();
  m_idealsurfnormy.clear();
  m_idealsurfnormz.clear();
  m_missurfcenterx.clear();
  m_missurfcentery.clear();
  m_missurfcenterz.clear();
  m_missurfnormx.clear();
  m_missurfnormy.clear();
  m_missurfnormz.clear();
  m_clusgxideal.clear();
  m_clusgyideal.clear();
  m_clusgzideal.clear();
  m_missurfalpha.clear();
  m_missurfbeta.clear();
  m_missurfgamma.clear();
  m_idealsurfalpha.clear();
  m_idealsurfbeta.clear();
  m_idealsurfgamma.clear();

  m_statelxglobderivdx.clear();
  m_statelxglobderivdy.clear();
  m_statelxglobderivdz.clear();
  m_statelxglobderivdalpha.clear();
  m_statelxglobderivdbeta.clear();
  m_statelxglobderivdgamma.clear();

  m_statelxlocderivd0.clear();
  m_statelxlocderivz0.clear();
  m_statelxlocderivphi.clear();
  m_statelxlocderivtheta.clear();
  m_statelxlocderivqop.clear();

  m_statelzglobderivdx.clear();
  m_statelzglobderivdy.clear();
  m_statelzglobderivdz.clear();
  m_statelzglobderivdalpha.clear();
  m_statelzglobderivdbeta.clear();
  m_statelzglobderivdgamma.clear();

  m_statelzlocderivd0.clear();
  m_statelzlocderivz0.clear();
  m_statelzlocderivphi.clear();
  m_statelzlocderivtheta.clear();
  m_statelzlocderivqop.clear();

  m_clusedge.clear();
  m_clusoverlap.clear();
  m_cluslx.clear();
  m_cluslz.clear();
  m_cluselx.clear();
  m_cluselz.clear();
  m_clusgr.clear();
  m_clusgx.clear();
  m_clusgy.clear();
  m_clusgz.clear();
  m_clusgxunmoved.clear();
  m_clusgyunmoved.clear();
  m_clusgzunmoved.clear();
  m_clusAdc.clear();
  m_clusMaxAdc.clear();
  m_cluslayer.clear();

  m_statelx.clear();
  m_statelz.clear();
  m_stateelx.clear();
  m_stateelz.clear();
  m_stategx.clear();
  m_stategy.clear();
  m_stategz.clear();
  m_statepx.clear();
  m_statepy.clear();
  m_statepz.clear();
  m_statepl.clear();
}
//____________________________________________________________________________..
int TrackResiduals::process_event(PHCompositeNode* topNode)
{
  auto *trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
  auto *clustermap = findNode::getClass<TrkrClusterContainer>(topNode, m_clusterContainerName);
  auto *geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  auto *hitmap = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  auto *tpcGeom =
      findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  auto *mvtxGeom = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");
  auto *inttGeom = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
  auto *mmGeom = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MICROMEGAS_FULL");

  if (!mmGeom)
  {
    mmGeom = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MICROMEGAS");
  }
  if (!trackmap or !clustermap or !geometry or (!hitmap && m_doHits))
  {
    std::cout << "Missing node, can't continue" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  auto *gl1 = findNode::getClass<Gl1RawHit>(topNode, "GL1RAWHIT");
  if (gl1)
  {
    m_bco = gl1->get_bco();
    auto lbshift = m_bco << 24U;
    m_bcotr = lbshift >> 24U;
  }
  else
  {
    Gl1Packet* gl1PacketInfo = findNode::getClass<Gl1Packet>(topNode, "GL1RAWHIT");
    if (!gl1PacketInfo)
    {
      m_bco = std::numeric_limits<uint64_t>::quiet_NaN();
      m_bcotr = std::numeric_limits<uint64_t>::quiet_NaN();
    }
    m_firedTriggers.clear();

    if (gl1PacketInfo)
    {
      m_gl1BunchCrossing = gl1PacketInfo->getBunchNumber();
      uint64_t triggervec = gl1PacketInfo->getScaledVector();
      m_bco = gl1PacketInfo->getBCO();
      auto lbshift = m_bco << 24U;
      m_bcotr = lbshift >> 24U;
      for (int i = 0; i < 64; i++)
      {
        bool trig_decision = ((triggervec & 0x1U) == 0x1U);
        if (trig_decision)
        {
          m_firedTriggers.push_back(i);
        }
        triggervec = (triggervec >> 1U) & 0xffffffffU;
      }
    }
  }
  
  m_mbdvtxz = std::numeric_limits<float>::quiet_NaN();

  MbdVertexMap *mbdvertexmap = findNode::getClass<MbdVertexMap>(topNode, "MbdVertexMap");
  if(mbdvertexmap)
  {
    for (auto & it : *mbdvertexmap)
    {
      MbdVertex* mbdvertex = it.second;
      if (mbdvertex)
      {
        m_mbdvtxz = mbdvertex->get_z();
      }
    }
  }
  MbdPmtContainer* bbcpmts = findNode::getClass<MbdPmtContainer>(topNode, "MbdPmtContainer");
  m_totalmbd = 0;
  if (bbcpmts)
  {
    int nPMTs = bbcpmts->get_npmt();
    for (int i = 0; i < nPMTs; i++)
    {
      m_totalmbd += bbcpmts->get_pmt(i)->get_q();
    }
  }
  else
  {
    if (Verbosity() > 0)
    {
      std::cout << "TrackResiduals::process_event: Could not find MbdPmtContainer," << std::endl;
    }
  }

  m_ntpcclus = 0;
  if (Verbosity() > 1)
  {
    std::cout << "Track map size is " << trackmap->size() << std::endl;
  }

  if (m_doHits)
  {
    fillHitTree(hitmap, geometry, tpcGeom, mvtxGeom, inttGeom, mmGeom);
  }

  if (m_doClusters)
  {
    fillClusterTree(clustermap, geometry);
  }

  if (m_convertSeeds)
  {
    fillResidualTreeSeeds(topNode);
  }
  else
  {
    fillResidualTreeKF(topNode);
  }

  if (m_doVertex)
  {
    fillVertexTree(topNode);
  }
  if (m_doEventTree)
  {
    fillEventTree(topNode);
  }
  m_event++;
  clearClusterStateVectors();
  return Fun4AllReturnCodes::EVENT_OK;
}

float TrackResiduals::calc_dedx(TrackSeed* tpcseed, TrkrClusterContainer* clustermap, PHG4TpcCylinderGeomContainer* tpcGeom)
{
  std::vector<TrkrDefs::cluskey> clusterKeys;
  clusterKeys.insert(clusterKeys.end(), tpcseed->begin_cluster_keys(),
                     tpcseed->end_cluster_keys());

  std::vector<float> dedxlist;
  for (unsigned long cluster_key : clusterKeys)
  {
    auto detid = TrkrDefs::getTrkrId(cluster_key);
    if (detid != TrkrDefs::TrkrId::tpcId)
    {
      continue;  // the micromegas clusters are added to the TPC seeds
    }
    unsigned int layer_local = TrkrDefs::getLayer(cluster_key);
    TrkrCluster* cluster = clustermap->findCluster(cluster_key);
    float adc = cluster->getAdc();
    PHG4TpcCylinderGeom* GeoLayer_local = tpcGeom->GetLayerCellGeom(layer_local);
    float thick = GeoLayer_local->get_thickness();
    float r = GeoLayer_local->get_radius();
    float alpha = (r * r) / (2 * r * TMath::Abs(1.0 / tpcseed->get_qOverR()));
    float beta = std::atan(tpcseed->get_slope());
    float alphacorr = std::cos(alpha);
    if (alphacorr < 0 || alphacorr > 4)
    {
      alphacorr = 4;
    }
    float betacorr = std::cos(beta);
    if (betacorr < 0 || betacorr > 4)
    {
      betacorr = 4;
    }
    adc /= thick;
    adc *= alphacorr;
    adc *= betacorr;
    dedxlist.push_back(adc);
    sort(dedxlist.begin(), dedxlist.end());
  }
  int trunc_min = 0;
  if (dedxlist.empty())
  {
    return std::numeric_limits<float>::quiet_NaN();
  }
  int trunc_max = (int) dedxlist.size() * 0.7;
  float sumdedx = 0;
  int ndedx = 0;
  for (int j = trunc_min; j <= trunc_max; j++)
  {
    sumdedx += dedxlist.at(j);
    ndedx++;
  }
  sumdedx /= ndedx;
  return sumdedx;
}

void TrackResiduals::fillFailedSeedTree(PHCompositeNode* topNode, std::set<unsigned int>& tpc_seed_ids)
{
  auto *tpcseedmap = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  auto *trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
  auto *clustermap = findNode::getClass<TrkrClusterContainer>(topNode, m_clusterContainerName);
  auto *geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  auto *silseedmap = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
  auto *svtxseedmap = findNode::getClass<TrackSeedContainer>(topNode, "SvtxTrackSeedContainer");
  auto *tpcGeo = findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");

  if (!tpcseedmap or !trackmap or !clustermap or !silseedmap or !svtxseedmap or !geometry)
  {
    std::cout << "Missing node, can't continue" << std::endl;
    return;
  }

  for (const auto& seed : *svtxseedmap)
  {
    if (!seed)
    {
      continue;
    }
    m_trackid = svtxseedmap->find(seed);
    auto tpcseedindex = seed->get_tpc_seed_index();
    if (tpc_seed_ids.find(tpcseedindex) != tpc_seed_ids.end())
    {
      continue;
    }
    auto siliconseedindex = seed->get_silicon_seed_index();
    auto *tpcseed = tpcseedmap->get(tpcseedindex);
    auto *silseed = silseedmap->get(siliconseedindex);

    int crossing = SHRT_MAX;
    if (silseed)
    {
      const auto si_pos = TrackSeedHelper::get_xyz(silseed);
      m_silseedx = si_pos.x();
      m_silseedy = si_pos.y();
      m_silseedz = si_pos.z();
      crossing = silseed->get_crossing();
    }
    else
    {
      const auto tpc_pos = TrackSeedHelper::get_xyz(tpcseed);
      m_tpcseedx = tpc_pos.x();
      m_tpcseedy = tpc_pos.y();
      m_tpcseedz = tpc_pos.z();
    }

    if (m_zeroField)
    {
      float pt = fabs(1. / tpcseed->get_qOverR()) * (0.3 / 100) * 0.01;
      float phi = tpcseed->get_phi();
      m_tpcseedpx = pt * std::cos(phi);
      m_tpcseedpy = pt * std::sin(phi);
      m_tpcseedpz = pt * std::cosh(tpcseed->get_eta()) * std::cos(tpcseed->get_theta());
    }
    else
    {
      m_tpcseedpx = tpcseed->get_px();
      m_tpcseedpy = tpcseed->get_py();
      m_tpcseedpz = tpcseed->get_pz();
    }
    m_tpcseedcharge = tpcseed->get_qOverR() > 0 ? 1 : -1;
    m_dedx = calc_dedx(tpcseed, clustermap, tpcGeo);
    m_nmaps = 0;
    m_nmapsstate = 0;
    m_nintt = 0;
    m_ninttstate = 0;
    m_ntpc = 0;
    m_ntpcstate = 0;
    m_nmms = 0;
    m_nmmsstate = 0;
    clearClusterStateVectors();
    for (auto *tseed : {silseed, tpcseed})
    {
      if (!tseed)
      {
        continue;
      }
      for (auto it = tseed->begin_cluster_keys(); it != tseed->end_cluster_keys(); ++it)
      {
        auto ckey = *it;
        auto *cluster = clustermap->findCluster(ckey);
        const Acts::Vector3 global = m_globalPositionWrapper.getGlobalPositionDistortionCorrected(ckey, cluster, crossing);
        const auto local = geometry->getLocalCoords(ckey, cluster);
        m_cluslx.push_back(local.x());
        m_cluslz.push_back(local.y());
        m_clusgx.push_back(global.x());
        m_clusgy.push_back(global.y());
        m_clusgz.push_back(global.z());
        float cr = r(global.x(), global.y());
        if (global.y() < 0)
        {
          cr = -cr;
        }
        m_clusgr.push_back(cr);
        auto detid = TrkrDefs::getTrkrId(ckey);
        if (detid == TrkrDefs::TrkrId::mvtxId)
        {
          m_nmaps++;
        }
        else if (detid == TrkrDefs::TrkrId::inttId)
        {
          m_nintt++;
        }
        else if (detid == TrkrDefs::TrkrId::tpcId)
        {
          m_ntpc++;
        }
        else if (detid == TrkrDefs::TrkrId::micromegasId)
        {
          m_nmms++;
        }
      }
    }
    m_failedfits->Fill();
  }
}
void TrackResiduals::fillVertexTree(PHCompositeNode* topNode)
{
  auto *svtxvertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  auto *trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
  auto *clustermap = findNode::getClass<TrkrClusterContainer>(topNode, m_clusterContainerName);
  if (svtxvertexmap)
  {
    m_nvertices = svtxvertexmap->size();
    clearClusterStateVectors();

    for (const auto& [key, vertex] : *svtxvertexmap)
    {
      m_vertexid = key;
      m_vertex_crossing = vertex->get_beam_crossing();
      m_vx = vertex->get_x();
      m_vy = vertex->get_y();
      m_vz = vertex->get_z();
      m_ntracks = vertex->size_tracks();

      for (auto it = vertex->begin_tracks(); it != vertex->end_tracks(); ++it)
      {
        auto id = *it;
        auto *track = trackmap->find(id)->second;
        if (!track)
        {
          continue;
        }
        for (const auto& ckey : get_cluster_keys(track))
        {
          TrkrCluster* cluster = clustermap->findCluster(ckey);

          Acts::Vector3 clusglob = m_globalPositionWrapper.getGlobalPositionDistortionCorrected(key, cluster, track->get_crossing());

          m_clusgx.push_back(clusglob.x());
          m_clusgy.push_back(clusglob.y());
          m_clusgz.push_back(clusglob.z());
          float clusr = r(clusglob.x(), clusglob.y());
          if (clusglob.y() < 0)
          {
            clusr = -clusr;
          }
          m_clusgr.push_back(clusr);
        }
      }

      m_vertextree->Fill();
    }
  }
}

float TrackResiduals::convertTimeToZ(ActsGeometry* geometry, TrkrDefs::cluskey cluster_key, TrkrCluster* cluster)
{
  // must convert local Y from cluster average time of arival to local cluster z position
  double drift_velocity = geometry->get_drift_velocity();
  double zdriftlength = cluster->getLocalY() * drift_velocity;
  double surfCenterZ = 52.89;                // 52.89 is where G4 thinks the surface center is
  double zloc = surfCenterZ - zdriftlength;  // converts z drift length to local z position in the TPC in north
  unsigned int side = TpcDefs::getSide(cluster_key);
  if (side == 0)
  {
    zloc = -zloc;
  }
  float z = zloc;  // in cm

  return z;
}
void TrackResiduals::circleFitClusters(
    std::vector<TrkrDefs::cluskey>& keys,
    TrkrClusterContainer* clusters,
    const short int& crossing)
{
  std::vector<Acts::Vector3> clusPos;
  std::vector<Acts::Vector3> global_vec;
  for (auto& key : keys)
  {
    auto *cluster = clusters->findCluster(key);
    const Acts::Vector3 pos = m_globalPositionWrapper.getGlobalPositionDistortionCorrected(key, cluster, crossing);
    clusPos.push_back(pos);
  }
  TrackFitUtils::position_vector_t yzpoints;
  TrackFitUtils::position_vector_t xypoints;

  for (auto& pos : clusPos)
  {
    float clusr = r(pos.x(), pos.y());
    // exclude silicon and tpot clusters for now
    if (std::fabs(clusr) > 80 || (m_linefitTPCOnly && std::fabs(clusr) < 20.))
    {
      continue;
    }
    xypoints.emplace_back(pos.x(), pos.y());
    yzpoints.emplace_back(pos.z(), pos.y());
    global_vec.push_back(pos);
  }

  auto xyparams = TrackFitUtils::line_fit(xypoints);
  auto yzLineParams = TrackFitUtils::line_fit(yzpoints);
  auto fitpars = TrackFitUtils::fitClusters(global_vec, keys, false);
  // auto fitpars = TrackFitUtils::fitClusters(global_vec, keys, !m_linefitTPCOnly);
  m_xyint = std::get<1>(xyparams);
  m_xyslope = std::get<0>(xyparams);
  m_yzint = std::get<1>(yzLineParams);
  m_yzslope = std::get<0>(yzLineParams);
  if (!fitpars.empty())
  {
    m_R = fitpars[0];
    m_X0 = fitpars[1];
    m_Y0 = fitpars[2];
    m_rzslope = fitpars[3];
    m_rzint = fitpars[4];
  }
  else
  {
    m_R = std::numeric_limits<float>::quiet_NaN();
    m_X0 = std::numeric_limits<float>::quiet_NaN();
    m_Y0 = std::numeric_limits<float>::quiet_NaN();
    m_rzslope = std::numeric_limits<float>::quiet_NaN();
    m_rzint = std::numeric_limits<float>::quiet_NaN();
  }
}

void TrackResiduals::lineFitClusters(std::vector<TrkrDefs::cluskey>& keys,
                                     TrkrClusterContainer* clusters,
                                     const short int& crossing)
{
  std::vector<Acts::Vector3> clusPos;
  for (auto& key : keys)
  {
    auto *cluster = clusters->findCluster(key);
    const Acts::Vector3 pos = m_globalPositionWrapper.getGlobalPositionDistortionCorrected(key, cluster, crossing);
    clusPos.push_back(pos);
  }
  TrackFitUtils::position_vector_t xypoints;
  TrackFitUtils::position_vector_t rzpoints;
  TrackFitUtils::position_vector_t yzpoints;
  for (auto& pos : clusPos)
  {
    float clusr = r(pos.x(), pos.y());
    if (pos.y() < 0)
    {
      clusr *= -1;
    }

    // exclude 1d tpot clusters for now

    if (std::fabs(clusr) > 80 || (m_linefitTPCOnly && std::fabs(clusr) < 20.))
    {
      continue;
    }

    rzpoints.emplace_back(pos.z(), clusr);
    xypoints.emplace_back(pos.x(), pos.y());
    yzpoints.emplace_back(pos.z(), pos.y());
  }

  auto xyparams = TrackFitUtils::line_fit(xypoints);
  auto rzparams = TrackFitUtils::line_fit(rzpoints);
  auto yzparams = TrackFitUtils::line_fit(yzpoints);
  m_xyint = std::get<1>(xyparams);
  m_xyslope = std::get<0>(xyparams);
  m_rzint = std::get<1>(rzparams);
  m_rzslope = std::get<0>(rzparams);
  m_yzint = std::get<1>(yzparams);
  m_yzslope = std::get<0>(yzparams);
}

void TrackResiduals::fillClusterTree(TrkrClusterContainer* clusters,
                                     ActsGeometry* geometry)
{
  if (clusters->size() < m_min_cluster_size)
  {
    return;
  }
  for (const auto& det : {TrkrDefs::TrkrId::mvtxId, TrkrDefs::TrkrId::inttId,
                    TrkrDefs::TrkrId::tpcId, TrkrDefs::TrkrId::micromegasId})
  {
    for (const auto& hitsetkey : clusters->getHitSetKeys(det))
    {
      m_scluslayer = TrkrDefs::getLayer(hitsetkey);
      auto range = clusters->getClusters(hitsetkey);
      for (auto iter = range.first; iter != range.second; ++iter)
      {
        auto key = iter->first;
        auto *cluster = clusters->findCluster(key);
        Acts::Vector3 glob;
        // NOT IMPLEMENTED YET
        // if (TrkrDefs::getTrkrId(key) == TrkrDefs::tpcId)
        // {
        //   glob = geometry->getGlobalPosition(key, cluster);  // corrections make no sense if crossing is not known
        // }
        // else
        // {
        //   glob = geometry->getGlobalPosition(key, cluster);
        // }
        glob = geometry->getGlobalPosition(key, cluster);
        m_sclusgx = glob.x();
        m_sclusgy = glob.y();
        m_sclusgz = glob.z();
        m_sclusgr = r(m_sclusgx, m_sclusgy);
        m_sclusphi = atan2(glob.y(), glob.x());
        m_scluseta = acos(glob.z() / std::sqrt(square(glob.x()) + square(glob.y()) + square(glob.z())));
        m_adc = cluster->getAdc();
        m_clusmaxadc = cluster->getMaxAdc();
        m_scluslx = cluster->getLocalX();
        m_scluslz = cluster->getLocalY();
        auto para_errors = m_clusErrPara.get_clusterv5_modified_error(cluster, m_sclusgr, key);
        m_phisize = cluster->getPhiSize();
        m_zsize = cluster->getZSize();
        m_scluselx = std::sqrt(para_errors.first);
        m_scluselz = std::sqrt(para_errors.second);

        //! Fill relevant geom info that is specific to subsystem
        switch (det)
        {
        case TrkrDefs::TrkrId::mvtxId:
          m_staveid = MvtxDefs::getStaveId(key);
          m_chipid = MvtxDefs::getChipId(key);
          m_strobeid = MvtxDefs::getStrobeId(key);

          m_ladderzid = std::numeric_limits<int>::quiet_NaN();
          m_ladderphiid = std::numeric_limits<int>::quiet_NaN();
          m_timebucket = std::numeric_limits<int>::quiet_NaN();
          m_clussector = std::numeric_limits<int>::quiet_NaN();
          m_side = std::numeric_limits<int>::quiet_NaN();
          m_segtype = std::numeric_limits<int>::quiet_NaN();
          m_tileid = std::numeric_limits<int>::quiet_NaN();
          break;
        case TrkrDefs::TrkrId::inttId:
          m_ladderzid = InttDefs::getLadderZId(key);
          m_ladderphiid = InttDefs::getLadderPhiId(key);
          m_timebucket = InttDefs::getTimeBucketId(key);

          m_staveid = std::numeric_limits<int>::quiet_NaN();
          m_chipid = std::numeric_limits<int>::quiet_NaN();
          m_strobeid = std::numeric_limits<int>::quiet_NaN();
          m_clussector = std::numeric_limits<int>::quiet_NaN();
          m_side = std::numeric_limits<int>::quiet_NaN();
          m_segtype = std::numeric_limits<int>::quiet_NaN();
          m_tileid = std::numeric_limits<int>::quiet_NaN();
          break;
        case TrkrDefs::TrkrId::tpcId:
          m_clussector = TpcDefs::getSectorId(key);
          m_side = TpcDefs::getSide(key);

          m_staveid = std::numeric_limits<int>::quiet_NaN();
          m_chipid = std::numeric_limits<int>::quiet_NaN();
          m_strobeid = std::numeric_limits<int>::quiet_NaN();
          m_ladderzid = std::numeric_limits<int>::quiet_NaN();
          m_ladderphiid = std::numeric_limits<int>::quiet_NaN();
          m_timebucket = std::numeric_limits<int>::quiet_NaN();
          m_segtype = std::numeric_limits<int>::quiet_NaN();
          m_tileid = std::numeric_limits<int>::quiet_NaN();
          break;
        case TrkrDefs::TrkrId::micromegasId:
          m_segtype = (int) MicromegasDefs::getSegmentationType(key);
          m_tileid = MicromegasDefs::getTileId(key);

          m_staveid = std::numeric_limits<int>::quiet_NaN();
          m_chipid = std::numeric_limits<int>::quiet_NaN();
          m_strobeid = std::numeric_limits<int>::quiet_NaN();
          m_ladderzid = std::numeric_limits<int>::quiet_NaN();
          m_ladderphiid = std::numeric_limits<int>::quiet_NaN();
          m_timebucket = std::numeric_limits<int>::quiet_NaN();
          m_clussector = std::numeric_limits<int>::quiet_NaN();
          m_side = std::numeric_limits<int>::quiet_NaN();
          break;
        default:
          break;
        }

        m_clustree->Fill();
      }
    }
  }
}

//____________________________________________________________________________..
int TrackResiduals::End(PHCompositeNode* /*unused*/)
{
  m_outfile->cd();
  m_tree->Write();
  if (m_doClusters)
  {
    m_clustree->Write();
  }
  if (m_doHits)
  {
    m_hittree->Write();
  }
  if (m_doVertex)
  {
    m_vertextree->Write();
  }
  if (m_doFailedSeeds)
  {
    m_failedfits->Write();
  }
  if (m_doEventTree)
  {
    m_eventtree->Write();
  }
  m_outfile->Close();

  return Fun4AllReturnCodes::EVENT_OK;
}
void TrackResiduals::fillHitTree(TrkrHitSetContainer* hitmap,
                                 ActsGeometry* geometry,
                                 PHG4TpcCylinderGeomContainer* tpcGeom,
                                 PHG4CylinderGeomContainer* mvtxGeom,
                                 PHG4CylinderGeomContainer* inttGeom,
                                 PHG4CylinderGeomContainer* mmGeom)
{
  if (!tpcGeom or !mvtxGeom or !inttGeom or !mmGeom)
  {
    std::cout << PHWHERE << "missing hit map, can't continue with hit tree"
              << std::endl;
    return;
  }
  TrkrHitSetContainer::ConstRange all_hitsets = hitmap->getHitSets();
  for (TrkrHitSetContainer::ConstIterator hitsetiter = all_hitsets.first;
       hitsetiter != all_hitsets.second;
       ++hitsetiter)
  {
    m_hitsetkey = hitsetiter->first;
    TrkrHitSet* hitset = hitsetiter->second;

    m_hitlayer = TrkrDefs::getLayer(m_hitsetkey);
    auto det = TrkrDefs::getTrkrId(m_hitsetkey);
    //! Fill relevant geom info that is specific to subsystem
    switch (det)
    {
    case TrkrDefs::TrkrId::mvtxId:
    {
      m_staveid = MvtxDefs::getStaveId(m_hitsetkey);
      m_chipid = MvtxDefs::getChipId(m_hitsetkey);
      m_strobeid = MvtxDefs::getStrobeId(m_hitsetkey);

      m_ladderzid = std::numeric_limits<int>::quiet_NaN();
      m_ladderphiid = std::numeric_limits<int>::quiet_NaN();
      m_timebucket = std::numeric_limits<int>::quiet_NaN();
      m_sector = std::numeric_limits<int>::quiet_NaN();
      m_side = std::numeric_limits<int>::quiet_NaN();
      m_segtype = std::numeric_limits<int>::quiet_NaN();
      m_tileid = std::numeric_limits<int>::quiet_NaN();
      break;
    }
    case TrkrDefs::TrkrId::inttId:
    {
      m_ladderzid = InttDefs::getLadderZId(m_hitsetkey);
      m_ladderphiid = InttDefs::getLadderPhiId(m_hitsetkey);
      m_timebucket = InttDefs::getTimeBucketId(m_hitsetkey);

      m_staveid = std::numeric_limits<int>::quiet_NaN();
      m_chipid = std::numeric_limits<int>::quiet_NaN();
      m_strobeid = std::numeric_limits<int>::quiet_NaN();
      m_sector = std::numeric_limits<int>::quiet_NaN();
      m_side = std::numeric_limits<int>::quiet_NaN();
      m_segtype = std::numeric_limits<int>::quiet_NaN();
      m_tileid = std::numeric_limits<int>::quiet_NaN();
      break;
    }
    case TrkrDefs::TrkrId::tpcId:
    {
      m_sector = TpcDefs::getSectorId(m_hitsetkey);
      m_side = TpcDefs::getSide(m_hitsetkey);

      m_staveid = std::numeric_limits<int>::quiet_NaN();
      m_chipid = std::numeric_limits<int>::quiet_NaN();
      m_strobeid = std::numeric_limits<int>::quiet_NaN();
      m_ladderzid = std::numeric_limits<int>::quiet_NaN();
      m_ladderphiid = std::numeric_limits<int>::quiet_NaN();
      m_timebucket = std::numeric_limits<int>::quiet_NaN();
      m_segtype = std::numeric_limits<int>::quiet_NaN();
      m_tileid = std::numeric_limits<int>::quiet_NaN();

      break;
    }
    case TrkrDefs::TrkrId::micromegasId:
    {
      m_segtype = (int) MicromegasDefs::getSegmentationType(m_hitsetkey);
      m_tileid = MicromegasDefs::getTileId(m_hitsetkey);

      m_staveid = std::numeric_limits<int>::quiet_NaN();
      m_chipid = std::numeric_limits<int>::quiet_NaN();
      m_strobeid = std::numeric_limits<int>::quiet_NaN();
      m_ladderzid = std::numeric_limits<int>::quiet_NaN();
      m_ladderphiid = std::numeric_limits<int>::quiet_NaN();
      m_timebucket = std::numeric_limits<int>::quiet_NaN();
      m_sector = std::numeric_limits<int>::quiet_NaN();
      m_side = std::numeric_limits<int>::quiet_NaN();
      break;
    }
    default:
      break;
    }

    // Got all stave/ladder/sector/tile info, now get the actual hit info
    auto hitrangei = hitset->getHits();
    for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
         hitr != hitrangei.second;
         ++hitr)
    {
      auto hitkey = hitr->first;
      auto *hit = hitr->second;
      m_adc = hit->getAdc();

      switch (det)
      {
      case TrkrDefs::TrkrId::mvtxId:
      {
        m_row = MvtxDefs::getRow(hitkey);
        m_col = MvtxDefs::getCol(hitkey);
        auto *layergeom = dynamic_cast<CylinderGeom_Mvtx*>(mvtxGeom->GetLayerGeom(m_hitlayer));
        auto local_coords = layergeom->get_local_coords_from_pixel(m_row, m_col);
        TVector2 local;
        local.SetX(local_coords.X());
        local.SetY(local_coords.Z());
        auto surf = geometry->maps().getSiliconSurface(m_hitsetkey);
        auto glob = CylinderGeom_MvtxHelper::get_world_from_local_coords(surf, geometry, local);
        m_hitgx = glob.X();
        m_hitgy = glob.Y();
        m_hitgz = glob.Z();

        m_segtype = std::numeric_limits<int>::quiet_NaN();
        m_tileid = std::numeric_limits<int>::quiet_NaN();
        m_strip = std::numeric_limits<int>::quiet_NaN();
        m_hitpad = std::numeric_limits<int>::quiet_NaN();
        m_hittbin = std::numeric_limits<int>::quiet_NaN();

        m_zdriftlength = std::numeric_limits<float>::quiet_NaN();
        break;
      }
      case TrkrDefs::TrkrId::inttId:
      {
        m_row = InttDefs::getRow(hitkey);
        m_col = InttDefs::getCol(hitkey);
        auto *geom = dynamic_cast<CylinderGeomIntt*>(inttGeom->GetLayerGeom(m_hitlayer));
        double local_hit_loc[3] = {0, 0, 0};
        geom->find_strip_center_localcoords(m_ladderzid, m_row, m_col, local_hit_loc);
        auto surf = geometry->maps().getSiliconSurface(m_hitsetkey);
        TVector2 local;
        local.SetX(local_hit_loc[1]);
        local.SetY(local_hit_loc[2]);
        auto glob = CylinderGeomInttHelper::get_world_from_local_coords(surf, geometry, local);

        m_hitgx = glob.X();
        m_hitgy = glob.Y();
        m_hitgz = glob.Z();
        m_segtype = std::numeric_limits<int>::quiet_NaN();
        m_tileid = std::numeric_limits<int>::quiet_NaN();
        m_strip = std::numeric_limits<int>::quiet_NaN();
        m_hitpad = std::numeric_limits<int>::quiet_NaN();
        m_hittbin = std::numeric_limits<int>::quiet_NaN();
        m_zdriftlength = std::numeric_limits<float>::quiet_NaN();
        break;
      }
      case TrkrDefs::TrkrId::tpcId:
      {
        m_row = std::numeric_limits<int>::quiet_NaN();
        m_col = std::numeric_limits<int>::quiet_NaN();
        m_segtype = std::numeric_limits<int>::quiet_NaN();
        m_tileid = std::numeric_limits<int>::quiet_NaN();
        m_strip = std::numeric_limits<int>::quiet_NaN();

        m_hitpad = TpcDefs::getPad(hitkey);
        m_hittbin = TpcDefs::getTBin(hitkey);

        auto *geoLayer = tpcGeom->GetLayerCellGeom(m_hitlayer);
        auto phi = geoLayer->get_phicenter(m_hitpad, m_side);
        auto radius = geoLayer->get_radius();
        float AdcClockPeriod = geoLayer->get_zstep();
        auto glob = geometry->getGlobalPositionTpc(m_hitsetkey, hitkey, phi, radius, AdcClockPeriod);
        m_hitgx = glob.x();
        m_hitgy = glob.y();
        m_hitgz = glob.z();

        break;
      }
      case TrkrDefs::TrkrId::micromegasId:
      {
        auto *const layergeom = dynamic_cast<CylinderGeomMicromegas*>(mmGeom->GetLayerGeom(m_hitlayer));
        m_strip = MicromegasDefs::getStrip(hitkey);
        const auto global_coord = layergeom->get_world_coordinates(m_tileid, geometry, m_strip);
        m_hitgx = global_coord.X();
        m_hitgy = global_coord.Y();
        m_hitgz = global_coord.Z();
        m_row = std::numeric_limits<int>::quiet_NaN();
        m_col = std::numeric_limits<int>::quiet_NaN();
        m_segtype = std::numeric_limits<int>::quiet_NaN();
        m_tileid = std::numeric_limits<int>::quiet_NaN();
        m_hitpad = std::numeric_limits<int>::quiet_NaN();
        m_hittbin = std::numeric_limits<int>::quiet_NaN();

        m_zdriftlength = std::numeric_limits<float>::quiet_NaN();
      }
      default:
        break;
      }

      m_hittree->Fill();
    }
  }
}

void TrackResiduals::fillClusterBranchesKF(TrkrDefs::cluskey ckey, SvtxTrack* track,
                                           const std::vector<std::pair<TrkrDefs::cluskey, Acts::Vector3>>& global,
                                           PHCompositeNode* topNode)
{
  auto *clustermap = findNode::getClass<TrkrClusterContainer>(topNode, m_clusterContainerName);
  auto *geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");

  auto global_moved = global;  // if use transient transforms for distortion correction
  if (m_use_clustermover)
  {
    // move the corrected cluster positions back to the original readout surface
    global_moved = m_clusterMover.processTrack(global);
  }

  ActsTransformations transformer;
  TrkrCluster* cluster = clustermap->findCluster(ckey);

  // loop over global vectors and get this cluster
  Acts::Vector3 clusglob(0, 0, 0);
  for (const auto& pair : global)
  {
    auto thiskey = pair.first;
    clusglob = pair.second;
    if (thiskey == ckey)
    {
      break;
    }
  }

  Acts::Vector3 clusglob_moved(0, 0, 0);
  for (const auto& pair : global_moved)
  {
    auto thiskey = pair.first;
    clusglob_moved = pair.second;
    if (thiskey == ckey)
    {
      break;
    }
  }

  unsigned int layer = TrkrDefs::getLayer(ckey);

  if (Verbosity() > 1)
  {
    std::cout << "Called :fillClusterBranchesKF for ckey " << ckey << " layer " << layer << " trackid " << track->get_id() << " clusglob x " << clusglob(0) << " y " << clusglob(1) << " z " << clusglob(2) << std::endl;
  }

  switch (TrkrDefs::getTrkrId(ckey))
  {
  case TrkrDefs::mvtxId:
    m_nmaps++;
    break;
  case TrkrDefs::inttId:
    m_nintt++;
    break;
  case TrkrDefs::tpcId:
    m_ntpc++;
    m_clsector.push_back(TpcDefs::getSectorId(ckey));
    m_clside.push_back(TpcDefs::getSide(ckey));
    break;
  case TrkrDefs::micromegasId:
    m_nmms++;
    m_tileid = MicromegasDefs::getTileId(ckey);
    break;
  default:
    std::cout << PHWHERE << " unknown key " << ckey << std::endl;
    gSystem->Exit(1);
    exit(1);
  }

  SvtxTrackState* state = nullptr;

  // the track states from the Acts fit are fitted to fully corrected clusters, and are on the surface
  for (auto state_iter = track->begin_states();
       state_iter != track->end_states();
       ++state_iter)
  {
    SvtxTrackState* tstate = state_iter->second;
    auto stateckey = tstate->get_cluskey();
    if (stateckey == ckey)
    {
      state = tstate;
      break;
    }
  }

  if (!state)
  {
    if (Verbosity() > 1)
    {
      std::cout << "   no state for cluster " << ckey << "  in layer " << layer << std::endl;
    }
  }
  else
  {
    switch (TrkrDefs::getTrkrId(ckey))
    {
    case TrkrDefs::mvtxId:
      m_nmapsstate++;
      break;
    case TrkrDefs::inttId:
      m_ninttstate++;
      break;
    case TrkrDefs::tpcId:
      m_ntpcstate++;
      break;
    case TrkrDefs::micromegasId:
      m_nmmsstate++;
      break;
    default:
      std::cout << PHWHERE << " unknown key " << ckey << std::endl;
      gSystem->Exit(1);
      exit(1);
    }
  }

  m_cluskeys.push_back(ckey);

  //! have cluster and state, fill vectors
  m_clusedge.push_back(cluster->getEdge());
  m_clusoverlap.push_back(cluster->getOverlap());

  // get new local coords from moved cluster
  Surface surf = geometry->maps().getSurface(ckey, cluster);
  Surface surf_ideal = geometry->maps().getSurface(ckey, cluster);  // Unchanged by distortion corrections
  // if this is a TPC cluster, the crossing correction may have moved it across the central membrane, check the surface
  auto trkrid = TrkrDefs::getTrkrId(ckey);
  if (trkrid == TrkrDefs::tpcId)
  {
    TrkrDefs::hitsetkey hitsetkey = TrkrDefs::getHitSetKeyFromClusKey(ckey);
    TrkrDefs::subsurfkey new_subsurfkey = 0;
    surf = geometry->get_tpc_surface_from_coords(hitsetkey, clusglob_moved, new_subsurfkey);
  }
  if (!surf)
  {
    if (Verbosity() > 2)
    {
      std::cout << " Failed to find surface for cluskey " << ckey << std::endl;
    }
    return;
  }

  // get local coordinates
  Acts::Vector2 loc;
  loc = geometry->getLocalCoords(ckey, cluster, m_crossing);
  if (m_use_clustermover)
  {
    // in this case we get local coords from transform of corrected global coords
    clusglob_moved *= Acts::UnitConstants::cm;  // we want mm for transformations
    Acts::Vector3 normal = surf->normal(geometry->geometry().getGeoContext(),
                                        Acts::Vector3(1, 1, 1), Acts::Vector3(1, 1, 1));
    auto local = surf->globalToLocal(geometry->geometry().getGeoContext(),
                                     clusglob_moved, normal);
    if (local.ok())
    {
      loc = local.value() / Acts::UnitConstants::cm;
    }
    else
    {
      // otherwise take the manual calculation for the TPC
      // doing it this way just avoids the bounds check that occurs in the surface class method
      Acts::Vector3 loct = surf->transform(geometry->geometry().getGeoContext()).inverse() * clusglob_moved;  // global is in mm
      loct /= Acts::UnitConstants::cm;

      loc(0) = loct(0);
      loc(1) = loct(1);
    }
    clusglob_moved /= Acts::UnitConstants::cm;  // we want cm for the tree
  }

  m_cluslx.push_back(loc.x());
  m_cluslz.push_back(loc.y());

  if (Verbosity() > 2)
  {
    std::cout << "Trackresiduals cluster (cm): localX " << loc.x() << " localY " << loc.y() << std::endl
              << " global.x " << clusglob_moved(0) << " global.y " << clusglob_moved(1) << " global.z " << clusglob_moved(2) << std::endl;
  }

  float clusr = r(clusglob_moved.x(), clusglob_moved.y());
  auto para_errors = m_clusErrPara.get_clusterv5_modified_error(cluster,
                                                                clusr, ckey);
  m_cluselx.push_back(sqrt(para_errors.first));
  m_cluselz.push_back(sqrt(para_errors.second));
  m_clusgx.push_back(clusglob_moved.x());
  m_clusgy.push_back(clusglob_moved.y());
  m_clusgr.push_back(clusglob_moved.y() > 0 ? clusr : -1 * clusr);
  m_clusgz.push_back(clusglob_moved.z());
  m_clusgxunmoved.push_back(clusglob.x());
  m_clusgyunmoved.push_back(clusglob.y());
  m_clusgzunmoved.push_back(clusglob.z());
  m_clusAdc.push_back(cluster->getAdc());
  m_clusMaxAdc.push_back(cluster->getMaxAdc());
  m_cluslayer.push_back(TrkrDefs::getLayer(ckey));
  m_clusphisize.push_back(cluster->getPhiSize());
  m_cluszsize.push_back(cluster->getZSize());

  auto misaligncenter = surf->center(geometry->geometry().getGeoContext());
  auto misalignnorm = -1 * surf->normal(geometry->geometry().getGeoContext(), Acts::Vector3(1, 1, 1), Acts::Vector3(1, 1, 1));
  auto misrot = surf->transform(geometry->geometry().getGeoContext()).rotation();

  float mgamma = atan2(-misrot(1, 0), misrot(0, 0));
  float mbeta = -asin(misrot(0, 1));
  float malpha = atan2(misrot(1, 1), misrot(2, 1));

  //! Switch to get ideal transforms
  alignmentTransformationContainer::use_alignment = false;
  auto idealcenter = surf_ideal->center(geometry->geometry().getGeoContext());
  auto idealnorm = -1 * surf_ideal->normal(geometry->geometry().getGeoContext(), Acts::Vector3(1, 1, 1), Acts::Vector3(1, 1, 1));

  // replace the corrected moved cluster local position with the readout position from ideal geometry for now
  // This allows us to see the distortion corrections by subtracting this uncorrected position
  // revisit this when looking at the alignment case
  //  Acts::Vector3 ideal_local(loc.x(), loc.y(), 0.0);
  auto nominal_loc = geometry->getLocalCoords(ckey, cluster);
  Acts::Vector3 ideal_local(nominal_loc.x(), nominal_loc.y(), 0.0);
  Acts::Vector3 ideal_glob = surf_ideal->transform(geometry->geometry().getGeoContext()) * (ideal_local * Acts::UnitConstants::cm);
  auto idealrot = surf_ideal->transform(geometry->geometry().getGeoContext()).rotation();

  //! These calculations are taken from the wikipedia page for Euler angles,
  //! under the Tait-Bryan angle explanation. Formulas for the angles
  //! calculated from the rotation matrices depending on what order the
  //! rotation matrix is constructed are given
  //! They need to be modified to conform to the Acts basis of (x,z,y), for
  //! which the wiki page expects (x,y,z). This includes swapping the sign
  //! of some elements to account for the permutation
  //! https://en.wikipedia.org/wiki/Euler_angles#Conversion_to_other_orientation_representations
  float igamma = atan2(-idealrot(1, 0), idealrot(0, 0));
  float ibeta = -asin(idealrot(0, 1));
  float ialpha = atan2(idealrot(1, 1), idealrot(2, 1));

  alignmentTransformationContainer::use_alignment = true;

  idealcenter /= Acts::UnitConstants::cm;
  misaligncenter /= Acts::UnitConstants::cm;
  ideal_glob /= Acts::UnitConstants::cm;

  m_idealsurfalpha.push_back(ialpha);
  m_idealsurfbeta.push_back(ibeta);
  m_idealsurfgamma.push_back(igamma);
  m_missurfalpha.push_back(malpha);
  m_missurfbeta.push_back(mbeta);
  m_missurfgamma.push_back(mgamma);

  m_idealsurfcenterx.push_back(idealcenter.x());
  m_idealsurfcentery.push_back(idealcenter.y());
  m_idealsurfcenterz.push_back(idealcenter.z());
  m_idealsurfnormx.push_back(idealnorm.x());
  m_idealsurfnormy.push_back(idealnorm.y());
  m_idealsurfnormz.push_back(idealnorm.z());
  m_missurfcenterx.push_back(misaligncenter.x());
  m_missurfcentery.push_back(misaligncenter.y());
  m_missurfcenterz.push_back(misaligncenter.z());
  m_missurfnormx.push_back(misalignnorm.x());
  m_missurfnormy.push_back(misalignnorm.y());
  m_missurfnormz.push_back(misalignnorm.z());
  m_clusgxideal.push_back(ideal_glob.x());
  m_clusgyideal.push_back(ideal_glob.y());
  m_clusgzideal.push_back(ideal_glob.z());

  if (state)
  {
    Acts::Vector3 stateglob(state->get_x(), state->get_y(), state->get_z());
    Acts::Vector2 stateloc(state->get_localX(), state->get_localY());

    if (Verbosity() > 2)
    {
      std::cout << "Trackresiduals state (cm): localX " << stateloc(0) << " localY " << stateloc(1) << std::endl
                << " stateglobx " << stateglob(0) << " stategloby " << stateglob(1) << " stateglobz " << stateglob(2) << std::endl;
    }

    const auto actscov =
        transformer.rotateSvtxTrackCovToActs(state);

    m_statelx.push_back(stateloc(0));
    m_statelz.push_back(stateloc(1));
    m_stateelx.push_back(std::sqrt(actscov(Acts::eBoundLoc0, Acts::eBoundLoc0)) / Acts::UnitConstants::cm);
    m_stateelz.push_back(std::sqrt(actscov(Acts::eBoundLoc1, Acts::eBoundLoc1)) / Acts::UnitConstants::cm);
    m_stategx.push_back(state->get_x());
    m_stategy.push_back(state->get_y());
    m_stategz.push_back(state->get_z());
    m_statepx.push_back(state->get_px());
    m_statepy.push_back(state->get_py());
    m_statepz.push_back(state->get_pz());
    m_statepl.push_back(state->get_pathlength());
  }
  else
  {
    // cluster has no corresponding state, set state variables to NaNs
    m_statelx.push_back(std::numeric_limits<float>::quiet_NaN());
    m_statelz.push_back(std::numeric_limits<float>::quiet_NaN());
    m_stategx.push_back(std::numeric_limits<float>::quiet_NaN());
    m_stategy.push_back(std::numeric_limits<float>::quiet_NaN());
    m_stategz.push_back(std::numeric_limits<float>::quiet_NaN());
    m_statepx.push_back(std::numeric_limits<float>::quiet_NaN());
    m_statepy.push_back(std::numeric_limits<float>::quiet_NaN());
    m_statepz.push_back(std::numeric_limits<float>::quiet_NaN());
    m_statepl.push_back(std::numeric_limits<float>::quiet_NaN());
  }

  if (Verbosity() > 2)
  {
    if (ideal_glob(2) > 0)
    {
      double xideal = ideal_glob.x();
      double yideal = ideal_glob.y();
      double zideal = ideal_glob.z();

      double xunmoved = clusglob.x();
      double yunmoved = clusglob.y();
      double zunmoved = clusglob.z();

      double xmoved = clusglob_moved.x();
      double ymoved = clusglob_moved.y();
      double zmoved = clusglob_moved.z();

      double this_radius_ideal = sqrt((xideal * xideal) + (yideal * yideal));
      double this_radius_unmoved = sqrt((xunmoved * xunmoved) + (yunmoved * yunmoved));

      double this_phi_unmoved = atan2(yunmoved, xunmoved);
      double this_phi_ideal = atan2(yideal, xideal);
      double this_phi_moved = atan2(ymoved, xmoved);

      std::cout << " global:  unmoved " << xunmoved << "  " << yunmoved << "  " << zunmoved << " this_phi " << this_phi_unmoved << " radius_unmoved " << this_radius_unmoved
                << " r*phi " << this_radius_ideal * this_phi_unmoved << std::endl;
      std::cout << "              global: ideal " << xideal << "  " << yideal << "  " << zideal << " this_phi_ideal " << this_phi_ideal << " radius_ideal " << this_radius_ideal
                << " r*phi " << this_radius_ideal * this_phi_ideal << std::endl;
      std::cout << "              global: moved " << xmoved << "  " << ymoved << "  " << zmoved << " phi_moved " << this_phi_moved
                << " r*phi " << this_radius_ideal * this_phi_moved << std::endl;
      std::cout << "                        d_radius " << this_radius_unmoved - this_radius_ideal << " clusgz " << zideal << " r*dphi " << this_radius_ideal * (this_phi_unmoved - this_phi_ideal) << std::endl;
    }
  }
}

void TrackResiduals::fillClusterBranchesSeeds(TrkrDefs::cluskey ckey,  // SvtxTrack* track,
                                              const std::vector<std::pair<TrkrDefs::cluskey, Acts::Vector3>>& global,
                                              PHCompositeNode* topNode)
{
  // The input map global contains the corrected cluster positions - NOT moved back to the surfacer.
  // When filling the residualtree:
  //    clusgx etc are the corrected - but not moved back to the surface - cluster positions
  //    clugxideal etc are the completely uncorrected cluster positions - they do not even have crossing corrections
  // CircleFitClusters is called in this method. It applies TOF, crossing, and all distortion corrections before fitting
  //    stategx etc are at the intersection point of the helical fit with the cluster surface

  auto *clustermap = findNode::getClass<TrkrClusterContainer>(topNode, m_clusterContainerName);
  auto *geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");

  auto global_moved = global;
  if (m_use_clustermover)
  {
    // move the corrected cluster positions back to the original readout surface
    global_moved = m_clusterMover.processTrack(global);
  }

  TrkrCluster* cluster = clustermap->findCluster(ckey);

  // loop over global vectors and get this cluster
  Acts::Vector3 clusglob(0, 0, 0);
  for (const auto& pair : global)
  {
    auto thiskey = pair.first;
    clusglob = pair.second;
    if (thiskey == ckey)
    {
      break;
    }
  }
  Acts::Vector3 clusglob_moved(0, 0, 0);
  for (const auto& pair : global_moved)
  {
    auto thiskey = pair.first;
    clusglob_moved = pair.second;
    if (thiskey == ckey)
    {
      break;
    }
  }

  switch (TrkrDefs::getTrkrId(ckey))
  {
  case TrkrDefs::mvtxId:
    m_nmaps++;
    m_clsector.push_back(-1);
    m_clside.push_back(-1);
    break;
  case TrkrDefs::inttId:
    m_clsector.push_back(-1);
    m_clside.push_back(-1);
    m_nintt++;
    break;
  case TrkrDefs::tpcId:
    m_ntpc++;
    m_clsector.push_back(TpcDefs::getSectorId(ckey));
    m_clside.push_back(TpcDefs::getSide(ckey));
    break;
  case TrkrDefs::micromegasId:
    m_nmms++;
    m_tileid = MicromegasDefs::getTileId(ckey);
    m_clsector.push_back(-1);
    m_clside.push_back(-1);
    break;
  default:
    std::cout << PHWHERE << " unknown key " << ckey << std::endl;
    gSystem->Exit(1);
    exit(1);
  }

  m_cluskeys.push_back(ckey);

  if (Verbosity() > 1)
  {
    std::cout << "Called fillClusterBranchesSeeds for ckey " << ckey << " clusglob x " << clusglob(0) << " y " << clusglob(1) << " z " << clusglob(2) << std::endl;
  }

  //! have cluster and state, fill vectors
  m_clusedge.push_back(cluster->getEdge());
  m_clusoverlap.push_back(cluster->getOverlap());

  // This is the nominal position of the cluster in local coords, completely uncorrected - is that what we want?
  auto loc = geometry->getLocalCoords(ckey, cluster);

  m_cluslx.push_back(loc.x());
  m_cluslz.push_back(loc.y());

  float clusr = r(clusglob_moved.x(), clusglob_moved.y());
  auto para_errors = m_clusErrPara.get_clusterv5_modified_error(cluster,
                                                                clusr, ckey);
  m_cluselx.push_back(sqrt(para_errors.first));
  m_cluselz.push_back(sqrt(para_errors.second));
  m_clusgx.push_back(clusglob_moved.x());
  m_clusgy.push_back(clusglob_moved.y());
  m_clusgr.push_back(clusglob_moved.y() > 0 ? clusr : -1 * clusr);
  m_clusgz.push_back(clusglob_moved.z());
  m_clusgxunmoved.push_back(clusglob.x());
  m_clusgyunmoved.push_back(clusglob.y());
  m_clusgzunmoved.push_back(clusglob.z());
  m_clusAdc.push_back(cluster->getAdc());
  m_clusMaxAdc.push_back(cluster->getMaxAdc());
  m_cluslayer.push_back(TrkrDefs::getLayer(ckey));
  m_clusphisize.push_back(cluster->getPhiSize());
  m_cluszsize.push_back(cluster->getZSize());

  if (Verbosity() > 1)
  {
    std::cout << "Track state/clus in layer "
              << (unsigned int) TrkrDefs::getLayer(ckey) << " with pos "
              << clusglob.transpose() << std::endl;
  }

  auto surf = geometry->maps().getSurface(ckey, cluster);

  auto misaligncenter = surf->center(geometry->geometry().getGeoContext());
  auto misalignnorm = -1 * surf->normal(geometry->geometry().getGeoContext(), Acts::Vector3(1, 1, 1), Acts::Vector3(1, 1, 1));
  auto misrot = surf->transform(geometry->geometry().getGeoContext()).rotation();

  float mgamma = atan2(-misrot(1, 0), misrot(0, 0));
  float mbeta = -asin(misrot(0, 1));
  float malpha = atan2(misrot(1, 1), misrot(2, 1));

  //! Switch to get ideal transforms
  alignmentTransformationContainer::use_alignment = false;
  auto idealcenter = surf->center(geometry->geometry().getGeoContext());
  auto idealnorm = -1 * surf->normal(geometry->geometry().getGeoContext(), Acts::Vector3(1, 1, 1), Acts::Vector3(1, 1, 1));
  Acts::Vector3 ideal_local(loc.x(), loc.y(), 0.0);
  Acts::Vector3 ideal_glob = surf->transform(geometry->geometry().getGeoContext()) * (ideal_local * Acts::UnitConstants::cm);
  auto idealrot = surf->transform(geometry->geometry().getGeoContext()).rotation();

  //! These calculations are taken from the wikipedia page for Euler angles,
  //! under the Tait-Bryan angle explanation. Formulas for the angles
  //! calculated from the rotation matrices depending on what order the
  //! rotation matrix is constructed are given
  //! They need to be modified to conform to the Acts basis of (x,z,y), for
  //! which the wiki page expects (x,y,z). This includes swapping the sign
  //! of some elements to account for the permutation
  //! https://en.wikipedia.org/wiki/Euler_angles#Conversion_to_other_orientation_representations
  float igamma = atan2(-idealrot(1, 0), idealrot(0, 0));
  float ibeta = -asin(idealrot(0, 1));
  float ialpha = atan2(idealrot(1, 1), idealrot(2, 1));

  alignmentTransformationContainer::use_alignment = true;

  idealcenter /= Acts::UnitConstants::cm;
  misaligncenter /= Acts::UnitConstants::cm;
  ideal_glob /= Acts::UnitConstants::cm;

  m_idealsurfalpha.push_back(ialpha);
  m_idealsurfbeta.push_back(ibeta);
  m_idealsurfgamma.push_back(igamma);
  m_missurfalpha.push_back(malpha);
  m_missurfbeta.push_back(mbeta);
  m_missurfgamma.push_back(mgamma);

  m_idealsurfcenterx.push_back(idealcenter.x());
  m_idealsurfcentery.push_back(idealcenter.y());
  m_idealsurfcenterz.push_back(idealcenter.z());
  m_idealsurfnormx.push_back(idealnorm.x());
  m_idealsurfnormy.push_back(idealnorm.y());
  m_idealsurfnormz.push_back(idealnorm.z());
  m_missurfcenterx.push_back(misaligncenter.x());
  m_missurfcentery.push_back(misaligncenter.y());
  m_missurfcenterz.push_back(misaligncenter.z());
  m_missurfnormx.push_back(misalignnorm.x());
  m_missurfnormy.push_back(misalignnorm.y());
  m_missurfnormz.push_back(misalignnorm.z());
  m_clusgxideal.push_back(ideal_glob.x());
  m_clusgyideal.push_back(ideal_glob.y());
  m_clusgzideal.push_back(ideal_glob.z());

  if (m_zeroField)
  {
    fillStatesWithLineFit(ckey, cluster, geometry);
  }
  else
  {
    fillStatesWithCircleFit(ckey, cluster, clusglob, geometry);
  }

  //! skip filling the state information if a state is not there
  //! or we just ran the seeding. Fill with Nans to maintain the
  //! 1-to-1 mapping between cluster/state vectors
  m_statepx.push_back(std::numeric_limits<float>::quiet_NaN());
  m_statepy.push_back(std::numeric_limits<float>::quiet_NaN());
  m_statepz.push_back(std::numeric_limits<float>::quiet_NaN());
  m_statepl.push_back(std::numeric_limits<float>::quiet_NaN());
  return;
}

void TrackResiduals::fillStatesWithCircleFit(const TrkrDefs::cluskey& key,
                                             TrkrCluster* cluster, Acts::Vector3& glob, ActsGeometry* geometry)
{
  auto surf = geometry->maps().getSurface(key, cluster);
  std::vector<float> fitpars;
  fitpars.push_back(m_R);
  fitpars.push_back(m_X0);
  fitpars.push_back(m_Y0);
  fitpars.push_back(m_rzslope);
  fitpars.push_back(m_rzint);

  auto intersection = TrackFitUtils::get_helix_surface_intersection(surf, fitpars, glob, geometry);
  m_stategx.push_back(intersection.x());
  m_stategy.push_back(intersection.y());
  m_stategz.push_back(intersection.z());

  auto result = surf->globalToLocal(geometry->geometry().getGeoContext(), intersection, Acts::Vector3(1, 1, 1));
  if (result.ok())
  {
    auto loc = result.value() / Acts::UnitConstants::cm;
    m_statelx.push_back(loc.x());
    m_statelz.push_back(loc.y());
  }
  else
  {
    auto local = (surf->transform(geometry->geometry().getGeoContext())).inverse() * (intersection * Acts::UnitConstants::cm);
    local /= Acts::UnitConstants::cm;
    m_statelx.push_back(local.x());
    m_statelz.push_back(local.y());
  }
}
void TrackResiduals::fillStatesWithLineFit(const TrkrDefs::cluskey& key,
                                           TrkrCluster* cluster, ActsGeometry* geometry)
{
  auto intersection = TrackFitUtils::surface_3Dline_intersection(key, cluster, geometry, m_xyslope,
                                                                 m_xyint, m_yzslope, m_yzint);

  auto surf = geometry->maps().getSurface(key, cluster);
  Acts::Vector3 surfnorm = surf->normal(geometry->geometry().getGeoContext(), Acts::Vector3(1, 1, 1), Acts::Vector3(1, 1, 1));
  if (!std::isnan(intersection.x()))
  {
    auto locstateres = surf->globalToLocal(geometry->geometry().getGeoContext(),
                                           intersection * Acts::UnitConstants::cm,
                                           surfnorm);
    if (locstateres.ok())
    {
      Acts::Vector2 loc = locstateres.value() / Acts::UnitConstants::cm;
      m_statelx.push_back(loc(0));
      m_statelz.push_back(loc(1));
    }
    else
    {
      Acts::Vector3 loct = surf->transform(geometry->geometry().getGeoContext()).inverse() * (intersection * Acts::UnitConstants::cm);
      loct /= Acts::UnitConstants::cm;
      m_statelx.push_back(loct(0));
      m_statelz.push_back(loct(1));
    }
    m_stategx.push_back(intersection.x());
    m_stategy.push_back(intersection.y());
    m_stategz.push_back(intersection.z());
  }
  else
  {
    //! otherwise the line is parallel to the surface, should not happen if
    //! we have a cluster on the surface but just fill the state vecs with nan
    m_statelx.push_back(std::numeric_limits<float>::quiet_NaN());
    m_statelz.push_back(std::numeric_limits<float>::quiet_NaN());
    m_stategx.push_back(std::numeric_limits<float>::quiet_NaN());
    m_stategy.push_back(std::numeric_limits<float>::quiet_NaN());
    m_stategz.push_back(std::numeric_limits<float>::quiet_NaN());
  }
}
void TrackResiduals::createBranches()
{
  if (m_doEventTree)
  {
    m_eventtree = new TTree("eventtree", "A tree with all hits");
    m_eventtree->Branch("run", &m_runnumber, "m_runnumber/I");
    m_eventtree->Branch("segment", &m_segment, "m_segment/I");
    m_eventtree->Branch("event", &m_event, "m_event/I");
    m_eventtree->Branch("gl1bco", &m_bco, "m_bco/I");
    m_eventtree->Branch("nmvtx", &m_nmvtx_all, "m_nmvtx_all/I");
    m_eventtree->Branch("nintt", &m_nintt_all, "m_nintt_all/I");
    m_eventtree->Branch("nhittpc0", &m_ntpc_hits0, "m_ntpc_hits0/I");
    m_eventtree->Branch("nhittpc1", &m_ntpc_hits1, "m_ntpc_hits1/I");
    m_eventtree->Branch("nclustpc0", &m_ntpc_clus0, "m_ntpc_clus0/I");
    m_eventtree->Branch("nclustpc1", &m_ntpc_clus1, "m_ntpc_clus1/I");
    m_eventtree->Branch("nmms", &m_nmms_all, "m_nmms_all/I");
    m_eventtree->Branch("nsiseed", &m_nsiseed, "m_nsiseed/I");
    m_eventtree->Branch("ntpcseed", &m_ntpcseed, "m_ntpcseed/I");
    m_eventtree->Branch("ntracks", &m_ntracks_all, "m_ntracks_all/I");
    m_eventtree->Branch("mbdcharge",&m_totalmbd, "m_totalmbd/F");
    m_eventtree->Branch("ntpcClusSector", &m_ntpc_clus_sector);
  }

  m_failedfits = new TTree("failedfits", "tree with seeds from failed Acts fits");
  m_failedfits->Branch("run", &m_runnumber, "m_runnumber/I");
  m_failedfits->Branch("segment", &m_segment, "m_segment/I");
  m_failedfits->Branch("trackid", &m_trackid, "m_trackid/I");
  m_failedfits->Branch("event", &m_event, "m_event/I");
  m_failedfits->Branch("silseedx", &m_silseedx, "m_silseedx/F");
  m_failedfits->Branch("silseedy", &m_silseedy, "m_silseedy/F");
  m_failedfits->Branch("silseedz", &m_silseedz, "m_silseedz/F");
  m_failedfits->Branch("tpcseedx", &m_tpcseedx, "m_tpcseedx/F");
  m_failedfits->Branch("tpcseedy", &m_tpcseedy, "m_tpcseedy/F");
  m_failedfits->Branch("tpcseedz", &m_tpcseedz, "m_tpcseedz/F");
  m_failedfits->Branch("tpcseedpx", &m_tpcseedpx, "m_tpcseedpx/F");
  m_failedfits->Branch("tpcseedpy", &m_tpcseedpy, "m_tpcseedpy/F");
  m_failedfits->Branch("tpcseedpz", &m_tpcseedpz, "m_tpcseedpz/F");
  m_failedfits->Branch("tpcseedcharge", &m_tpcseedcharge, "m_tpcseedcharge/I");
  m_failedfits->Branch("dedx", &m_dedx, "m_dedx/F");
  m_failedfits->Branch("nmaps", &m_nmaps, "m_nmaps/I");
  m_failedfits->Branch("nintt", &m_nintt, "m_nintt/I");
  m_failedfits->Branch("ntpc", &m_ntpc, "m_ntpc/I");
  m_failedfits->Branch("nmms", &m_nmms, "m_nmms/I");
  m_failedfits->Branch("gx", &m_clusgx);
  m_failedfits->Branch("gy", &m_clusgy);
  m_failedfits->Branch("gz", &m_clusgz);
  m_failedfits->Branch("gr", &m_clusgr);
  m_failedfits->Branch("lx", &m_cluslx);
  m_failedfits->Branch("lz", &m_cluslz);

  m_vertextree = new TTree("vertextree", "tree with vertices");
  m_vertextree->Branch("run", &m_runnumber, "m_runnumber/I");
  m_vertextree->Branch("segment", &m_segment, "m_segment/I");
  m_vertextree->Branch("event", &m_event, "m_event/I");
  m_vertextree->Branch("firedTriggers", &m_firedTriggers);
  m_vertextree->Branch("gl1BunchCrossing", &m_gl1BunchCrossing, "m_gl1BunchCrossing/l");
  m_vertextree->Branch("gl1bco", &m_bco, "m_bco/l");
  m_vertextree->Branch("trbco", &m_bcotr, "m_bcotr/l");
  m_vertextree->Branch("vertexid", &m_vertexid);
  m_vertextree->Branch("vertex_crossing", &m_vertex_crossing, "m_vertex_crossing/I");
  m_vertextree->Branch("mbdzvtx", &m_mbdvtxz, "m_mbdvtxz/F");
  m_vertextree->Branch("vx", &m_vx, "m_vx/F");
  m_vertextree->Branch("vy", &m_vy, "m_vy/F");
  m_vertextree->Branch("vz", &m_vz, "m_vz/F");
  m_vertextree->Branch("ntracks", &m_ntracks, "m_ntracks/I");
  m_vertextree->Branch("nvertices", &m_nvertices, "m_nvertices/I");
  m_vertextree->Branch("gx", &m_clusgx);
  m_vertextree->Branch("gy", &m_clusgy);
  m_vertextree->Branch("gz", &m_clusgz);
  m_vertextree->Branch("gr", &m_clusgr);
  m_vertextree->Branch("mbdcharge", &m_totalmbd, "m_totalmbd/F");

  m_hittree = new TTree("hittree", "A tree with all hits");
  m_hittree->Branch("run", &m_runnumber, "m_runnumber/I");
  m_hittree->Branch("segment", &m_segment, "m_segment/I");
  m_hittree->Branch("event", &m_event, "m_event/I");
  m_hittree->Branch("gl1bco", &m_bco, "m_bco/l");
  m_hittree->Branch("hitsetkey", &m_hitsetkey, "m_hitsetkey/i");
  m_hittree->Branch("gx", &m_hitgx, "m_hitgx/F");
  m_hittree->Branch("gy", &m_hitgy, "m_hitgy/F");
  m_hittree->Branch("gz", &m_hitgz, "m_hitgz/F");
  m_hittree->Branch("layer", &m_hitlayer, "m_hitlayer/I");
  m_hittree->Branch("sector", &m_sector, "m_sector/I");
  m_hittree->Branch("side", &m_side, "m_side/I");
  m_hittree->Branch("stave", &m_staveid, "m_staveid/I");
  m_hittree->Branch("chip", &m_chipid, "m_chipid/I");
  m_hittree->Branch("strobe", &m_strobeid, "m_strobeid/I");
  m_hittree->Branch("ladderz", &m_ladderzid, "m_ladderzid/I");
  m_hittree->Branch("ladderphi", &m_ladderphiid, "m_ladderphiid/I");
  m_hittree->Branch("timebucket", &m_timebucket, "m_timebucket/I");
  m_hittree->Branch("pad", &m_hitpad, "m_hitpad/I");
  m_hittree->Branch("tbin", &m_hittbin, "m_hittbin/I");
  m_hittree->Branch("col", &m_col, "m_col/I");
  m_hittree->Branch("row", &m_row, "m_row/I");
  m_hittree->Branch("segtype", &m_segtype, "m_segtype/I");
  m_hittree->Branch("tile", &m_tileid, "m_tileid/I");
  m_hittree->Branch("strip", &m_strip, "m_strip/I");
  m_hittree->Branch("adc", &m_adc, "m_adc/F");
  m_hittree->Branch("zdriftlength", &m_zdriftlength, "m_zdriftlength/F");
  m_hittree->Branch("mbdcharge",&m_totalmbd, "m_totalmbd/F");

  m_clustree = new TTree("clustertree", "A tree with all clusters");
  m_clustree->Branch("run", &m_runnumber, "m_runnumber/I");
  m_clustree->Branch("segment", &m_segment, "m_segment/I");
  m_clustree->Branch("event", &m_event, "m_event/I");
  m_clustree->Branch("gl1bco", &m_bco, "m_bco/l");
  m_clustree->Branch("lx", &m_scluslx, "m_scluslx/F");
  m_clustree->Branch("lz", &m_scluslz, "m_scluslz/F");
  m_clustree->Branch("gx", &m_sclusgx, "m_sclusgx/F");
  m_clustree->Branch("gy", &m_sclusgy, "m_sclusgy/F");
  m_clustree->Branch("gz", &m_sclusgz, "m_sclusgz/F");
  m_clustree->Branch("phi", &m_sclusphi, "m_sclusphi/F");
  m_clustree->Branch("eta", &m_scluseta, "m_scluseta/F");
  m_clustree->Branch("adc", &m_adc, "m_adc/F");
  m_clustree->Branch("phisize", &m_phisize, "m_phisize/I");
  m_clustree->Branch("zsize", &m_zsize, "m_zsize/I");
  m_clustree->Branch("erphi", &m_scluselx, "m_scluselx/F");
  m_clustree->Branch("ez", &m_scluselz, "m_scluselz/F");
  m_clustree->Branch("maxadc", &m_clusmaxadc, "m_clusmaxadc/F");
  m_clustree->Branch("sector", &m_clussector, "m_clussector/I");
  m_clustree->Branch("side", &m_side, "m_side/I");
  m_clustree->Branch("stave", &m_staveid, "m_staveid/I");
  m_clustree->Branch("chip", &m_chipid, "m_chipid/I");
  m_clustree->Branch("strobe", &m_strobeid, "m_strobeid/I");
  m_clustree->Branch("ladderz", &m_ladderzid, "m_ladderzid/I");
  m_clustree->Branch("ladderphi", &m_ladderphiid, "m_ladderphiid/I");
  m_clustree->Branch("timebucket", &m_timebucket, "m_timebucket/I");
  m_clustree->Branch("segtype", &m_segtype, "m_segtype/I");
  m_clustree->Branch("tile", &m_tileid, "m_tileid/I");

  m_tree = new TTree("residualtree", "A tree with track, cluster, and state info");
  m_tree->Branch("run", &m_runnumber, "m_runnumber/I");
  m_tree->Branch("segment", &m_segment, "m_segment/I");
  m_tree->Branch("event", &m_event, "m_event/I");
  m_tree->Branch("mbdcharge",&m_totalmbd, "m_totalmbd/F");
  m_tree->Branch("mbdzvtx", &m_mbdvtxz, "m_mbdvtxz/F");
  m_tree->Branch("firedTriggers", &m_firedTriggers);
  m_tree->Branch("gl1BunchCrossing", &m_gl1BunchCrossing, "m_gl1BunchCrossing/l");
  m_tree->Branch("trackid", &m_trackid, "m_trackid/I");
  m_tree->Branch("tpcid", &m_tpcid, "m_tpcid/I");
  m_tree->Branch("silid", &m_silid, "m_silid/I");
  m_tree->Branch("gl1bco", &m_bco, "m_bco/l");
  m_tree->Branch("crossing", &m_crossing, "m_crossing/I");
  m_tree->Branch("crossing_estimate", &m_crossing_estimate, "m_crossing_estimate/I");
  m_tree->Branch("silseedx", &m_silseedx, "m_silseedx/F");
  m_tree->Branch("silseedy", &m_silseedy, "m_silseedy/F");
  m_tree->Branch("silseedz", &m_silseedz, "m_silseedz/F");
  m_tree->Branch("silseedpx", &m_silseedpx, "m_silseedpx/F");
  m_tree->Branch("silseedpy", &m_silseedpy, "m_silseedpy/F");
  m_tree->Branch("silseedpz", &m_silseedpz, "m_silseedpz/F");
  m_tree->Branch("silseedphi", &m_silseedphi, "m_silseedphi/F");
  m_tree->Branch("silseedeta", &m_silseedeta, "m_silseedeta/F");
  m_tree->Branch("silseedcharge", &m_silseedcharge, "m_silseedcharge/I");
  m_tree->Branch("tpcseedx", &m_tpcseedx, "m_tpcseedx/F");
  m_tree->Branch("tpcseedy", &m_tpcseedy, "m_tpcseedy/F");
  m_tree->Branch("tpcseedz", &m_tpcseedz, "m_tpcseedz/F");
  m_tree->Branch("tpcseedpx", &m_tpcseedpx, "m_tpcseedpx/F");
  m_tree->Branch("tpcseedpy", &m_tpcseedpy, "m_tpcseedpy/F");
  m_tree->Branch("tpcseedpz", &m_tpcseedpz, "m_tpcseedpz/F");
  m_tree->Branch("tpcseedphi", &m_tpcseedphi, "m_tpcseedphi/F");
  m_tree->Branch("tpcseedeta", &m_tpcseedeta, "m_tpcseedeta/F");
  m_tree->Branch("tpcseedcharge", &m_tpcseedcharge, "m_tpcseedcharge/I");
  m_tree->Branch("dedx", &m_dedx, "m_dedx/F");
  m_tree->Branch("tracklength", &m_tracklength, "m_tracklength/F");
  m_tree->Branch("px", &m_px, "m_px/F");
  m_tree->Branch("py", &m_py, "m_py/F");
  m_tree->Branch("pz", &m_pz, "m_pz/F");
  m_tree->Branch("pt", &m_pt, "m_pt/F");
  m_tree->Branch("eta", &m_eta, "m_eta/F");
  m_tree->Branch("phi", &m_phi, "m_phi/F");
  m_tree->Branch("deltapt", &m_deltapt, "m_deltapt/F");
  m_tree->Branch("charge", &m_charge, "m_charge/I");
  m_tree->Branch("quality", &m_quality, "m_quality/F");
  m_tree->Branch("ndf", &m_ndf, "m_ndf/F");
  m_tree->Branch("nhits", &m_nhits, "m_nhits/I");
  m_tree->Branch("nmaps", &m_nmaps, "m_nmaps/I");
  m_tree->Branch("nmapsstate", &m_nmapsstate, "m_nmapsstate/I");
  m_tree->Branch("nintt", &m_nintt, "m_nintt/I");
  m_tree->Branch("ninttstate", &m_ninttstate, "m_ninttstate/I");
  m_tree->Branch("ntpc", &m_ntpc, "m_ntpc/I");
  m_tree->Branch("ntpcstate", &m_ntpcstate, "m_ntpcstate/I");
  m_tree->Branch("nmms", &m_nmms, "m_nmms/I");
  m_tree->Branch("nmmsstate", &m_nmmsstate, "m_nmmsstate/I");
  m_tree->Branch("tile", &m_tileid, "m_tileid/I");
  m_tree->Branch("vertexid", &m_vertexid);
  m_tree->Branch("vertex_crossing", &m_vertex_crossing, "m_vertex_crossing/I");
  m_tree->Branch("vx", &m_vx, "m_vx/F");
  m_tree->Branch("vy", &m_vy, "m_vy/F");
  m_tree->Branch("vz", &m_vz, "m_vz/F");
  m_tree->Branch("vertex_ntracks", &m_vertex_ntracks, "m_vertex_ntracks/I");
  m_tree->Branch("pcax", &m_pcax, "m_pcax/F");
  m_tree->Branch("pcay", &m_pcay, "m_pcay/F");
  m_tree->Branch("pcaz", &m_pcaz, "m_pcaz/F");
  m_tree->Branch("rzslope", &m_rzslope, "m_rzslope/F");
  m_tree->Branch("xyslope", &m_xyslope, "m_xyslope/F");
  m_tree->Branch("yzslope", &m_yzslope, "m_yzslope/F");
  m_tree->Branch("rzint", &m_rzint, "m_rzint/F");
  m_tree->Branch("xyint", &m_xyint, "m_xyint/F");
  m_tree->Branch("yzint", &m_yzint, "m_yzint/F");
  m_tree->Branch("R", &m_R, "m_R/F");
  m_tree->Branch("X0", &m_X0, "m_X0/F");
  m_tree->Branch("Y0", &m_Y0, "m_Y0/F");
  m_tree->Branch("dcaxy", &m_dcaxy, "m_dcaxy/F");
  m_tree->Branch("dcaz", &m_dcaz, "m_dcaz/F");

  m_tree->Branch("clussector", &m_clsector);
  m_tree->Branch("cluslayer", &m_cluslayer);
  m_tree->Branch("clusside", &m_clside);
  m_tree->Branch("cluskeys", &m_cluskeys);
  m_tree->Branch("clusedge", &m_clusedge);
  m_tree->Branch("clusoverlap", &m_clusoverlap);
  m_tree->Branch("cluslx", &m_cluslx);
  m_tree->Branch("cluslz", &m_cluslz);
  m_tree->Branch("cluselx", &m_cluselx);
  m_tree->Branch("cluselz", &m_cluselz);
  m_tree->Branch("clusgx", &m_clusgx);
  m_tree->Branch("clusgy", &m_clusgy);
  m_tree->Branch("clusgz", &m_clusgz);
  m_tree->Branch("clusgr", &m_clusgr);
  if (m_doAlignment)
  {
    m_tree->Branch("clusgxunmoved", &m_clusgxunmoved);
    m_tree->Branch("clusgyunmoved", &m_clusgyunmoved);
    m_tree->Branch("clusgzunmoved", &m_clusgzunmoved);
  }
  m_tree->Branch("clusAdc", &m_clusAdc);
  m_tree->Branch("clusMaxAdc", &m_clusMaxAdc);
  m_tree->Branch("clusphisize", &m_clusphisize);
  m_tree->Branch("cluszsize", &m_cluszsize);

  if (m_doAlignment)
  {
    m_tree->Branch("idealsurfcenterx", &m_idealsurfcenterx);
    m_tree->Branch("idealsurfcentery", &m_idealsurfcentery);
    m_tree->Branch("idealsurfcenterz", &m_idealsurfcenterz);
    m_tree->Branch("idealsurfnormx", &m_idealsurfnormx);
    m_tree->Branch("idealsurfnormy", &m_idealsurfnormy);
    m_tree->Branch("idealsurfnormz", &m_idealsurfnormz);
    m_tree->Branch("missurfcenterx", &m_missurfcenterx);
    m_tree->Branch("missurfcentery", &m_missurfcentery);
    m_tree->Branch("missurfcenterz", &m_missurfcenterz);
    m_tree->Branch("missurfnormx", &m_missurfnormx);
    m_tree->Branch("missurfnormy", &m_missurfnormy);
    m_tree->Branch("missurfnormz", &m_missurfnormz);
    m_tree->Branch("clusgxideal", &m_clusgxideal);
    m_tree->Branch("clusgyideal", &m_clusgyideal);
    m_tree->Branch("clusgzideal", &m_clusgzideal);
    m_tree->Branch("missurfalpha", &m_missurfalpha);
    m_tree->Branch("missurfbeta", &m_missurfbeta);
    m_tree->Branch("missurfgamma", &m_missurfgamma);
    m_tree->Branch("idealsurfalpha", &m_idealsurfalpha);
    m_tree->Branch("idealsurfbeta", &m_idealsurfbeta);
    m_tree->Branch("idealsurfgamma", &m_idealsurfgamma);
  }

  m_tree->Branch("statelx", &m_statelx);
  m_tree->Branch("statelz", &m_statelz);
  m_tree->Branch("stateelx", &m_stateelx);
  m_tree->Branch("stateelz", &m_stateelz);
  m_tree->Branch("stategx", &m_stategx);
  m_tree->Branch("stategy", &m_stategy);
  m_tree->Branch("stategz", &m_stategz);
  m_tree->Branch("statepx", &m_statepx);
  m_tree->Branch("statepy", &m_statepy);
  m_tree->Branch("statepz", &m_statepz);
  m_tree->Branch("statepl", &m_statepl);

  if (m_doAlignment)
  {
    m_tree->Branch("statelxglobderivdx", &m_statelxglobderivdx);
    m_tree->Branch("statelxglobderivdy", &m_statelxglobderivdy);
    m_tree->Branch("statelxglobderivdz", &m_statelxglobderivdz);
    m_tree->Branch("statelxglobderivdalpha", &m_statelxglobderivdalpha);
    m_tree->Branch("statelxglobderivdbeta", &m_statelxglobderivdbeta);
    m_tree->Branch("statelxglobderivdgamma", &m_statelxglobderivdgamma);

    m_tree->Branch("statelxlocderivd0", &m_statelxlocderivd0);
    m_tree->Branch("statelxlocderivz0", &m_statelxlocderivz0);
    m_tree->Branch("statelxlocderivphi", &m_statelxlocderivphi);
    m_tree->Branch("statelxlocderivtheta", &m_statelxlocderivtheta);
    m_tree->Branch("statelxlocderivqop", &m_statelxlocderivqop);

    m_tree->Branch("statelzglobderivdx", &m_statelzglobderivdx);
    m_tree->Branch("statelzglobderivdy", &m_statelzglobderivdy);
    m_tree->Branch("statelzglobderivdz", &m_statelzglobderivdz);
    m_tree->Branch("statelzglobderivdalpha", &m_statelzglobderivdalpha);
    m_tree->Branch("statelzglobderivdbeta", &m_statelzglobderivdbeta);
    m_tree->Branch("statelzglobderivdgamma", &m_statelzglobderivdgamma);

    m_tree->Branch("statelzlocderivd0", &m_statelzlocderivd0);
    m_tree->Branch("statelzlocderivz0", &m_statelzlocderivz0);
    m_tree->Branch("statelzlocderivphi", &m_statelzlocderivphi);
    m_tree->Branch("statelzlocderivtheta", &m_statelzlocderivtheta);
    m_tree->Branch("statelzlocderivqop", &m_statelzlocderivqop);
  }
}

void TrackResiduals::fillResidualTreeKF(PHCompositeNode* topNode)
{
  auto *silseedmap = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
  auto *tpcseedmap = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  auto *tpcGeom =
      findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  auto *trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
  auto *clustermap = findNode::getClass<TrkrClusterContainer>(topNode, m_clusterContainerName);
  auto *vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  auto *alignmentmap = findNode::getClass<SvtxAlignmentStateMap>(topNode, m_alignmentMapName);

  std::set<unsigned int> tpc_seed_ids;
  for (const auto& [key, track] : *trackmap)
  {
    if (!track)
    {
      continue;
    }
    m_trackid = track->get_id();

    m_crossing = track->get_crossing();
    m_crossing_estimate = SHRT_MAX;
    m_px = track->get_px();
    m_py = track->get_py();
    m_pz = track->get_pz();

    m_pt = std::sqrt(square(m_px) + square(m_py));
    m_eta = std::atanh(m_pz / std::sqrt(square(m_pt) + square(m_pz)));
    m_phi = std::atan2(m_py, m_px);
    float CVxx = track->get_error(3, 3);
    float CVxy = track->get_error(3, 4);
    float CVyy = track->get_error(4, 4);
    m_deltapt = std::sqrt((CVxx * square(m_px) + 2 * CVxy * m_px * m_py + CVyy * square(m_py)) / (square(m_px) + square(m_py)));

    m_charge = track->get_charge();
    m_quality = track->get_quality();
    m_chisq = track->get_chisq();
    m_ndf = track->get_ndf();

    m_nmaps = 0;
    m_nmapsstate = 0;
    m_nintt = 0;
    m_ninttstate = 0;
    m_ntpc = 0;
    m_ntpcstate = 0;
    m_nmms = 0;
    m_nmmsstate = 0;
    m_silid = std::numeric_limits<unsigned int>::quiet_NaN();
    m_tpcid = std::numeric_limits<unsigned int>::quiet_NaN();
    m_silseedx = std::numeric_limits<float>::quiet_NaN();
    m_silseedy = std::numeric_limits<float>::quiet_NaN();
    m_silseedz = std::numeric_limits<float>::quiet_NaN();
    m_silseedpx = std::numeric_limits<float>::quiet_NaN();
    m_silseedpy = std::numeric_limits<float>::quiet_NaN();
    m_silseedpz = std::numeric_limits<float>::quiet_NaN();
    m_silseedphi = std::numeric_limits<float>::quiet_NaN();
    m_silseedeta = std::numeric_limits<float>::quiet_NaN();
    m_silseedcharge = std::numeric_limits<int>::quiet_NaN();
    m_tpcseedx = std::numeric_limits<float>::quiet_NaN();
    m_tpcseedy = std::numeric_limits<float>::quiet_NaN();
    m_tpcseedz = std::numeric_limits<float>::quiet_NaN();
    m_tpcseedpx = std::numeric_limits<float>::quiet_NaN();
    m_tpcseedpy = std::numeric_limits<float>::quiet_NaN();
    m_tpcseedpz = std::numeric_limits<float>::quiet_NaN();
    m_tpcseedphi = std::numeric_limits<float>::quiet_NaN();
    m_tpcseedeta = std::numeric_limits<float>::quiet_NaN();
    m_tpcseedcharge = std::numeric_limits<int>::quiet_NaN();

    m_pcax = track->get_x();
    m_pcay = track->get_y();
    m_pcaz = track->get_z();
    m_dcaxy = std::numeric_limits<float>::quiet_NaN();
    m_dcaz = std::numeric_limits<float>::quiet_NaN();
    m_vertexid = track->get_vertex_id();
    if (vertexmap)
    {
      auto vertexit = vertexmap->find(m_vertexid);
      if (vertexit != vertexmap->end())
      {
        auto *vertex = vertexit->second;
        m_vx = vertex->get_x();
        m_vy = vertex->get_y();
        m_vz = vertex->get_z();
        m_vertex_ntracks = vertex->size_tracks();
        Acts::Vector3 v(m_vx, m_vy, m_vz);
        auto dcapair = TrackAnalysisUtils::get_dca(track, v);
        m_dcaxy = dcapair.first.first;
        m_dcaz = dcapair.second.first;
      }
    }

    auto *tpcseed = track->get_tpc_seed();
    if (tpcseed)
    {
      m_tpcid = tpcseedmap->find(tpcseed);
      tpc_seed_ids.insert(tpcseedmap->find(tpcseed));
    }
    auto *silseed = track->get_silicon_seed();
    if (silseed)
    {
      m_silid = silseedmap->find(silseed);

      const auto si_pos = TrackSeedHelper::get_xyz(silseed);
      m_silseedx = si_pos.x();
      m_silseedy = si_pos.y();
      m_silseedz = si_pos.z();
      m_silseedpx = silseed->get_px();
      m_silseedpy = silseed->get_py();
      m_silseedpz = silseed->get_pz();
      m_silseedphi = silseed->get_phi();
      m_silseedeta = silseed->get_eta();
      m_silseedcharge = silseed->get_qOverR() > 0 ? 1 : -1;
    }
    if (tpcseed)
    {
      const auto tpc_pos = TrackSeedHelper::get_xyz(tpcseed);
      m_tpcseedx = tpc_pos.x();
      m_tpcseedy = tpc_pos.y();
      m_tpcseedz = tpc_pos.z();
      m_tpcseedpx = tpcseed->get_px();
      m_tpcseedpy = tpcseed->get_py();
      m_tpcseedpz = tpcseed->get_pz();
      m_tpcseedphi = tpcseed->get_phi();
      m_tpcseedeta = tpcseed->get_eta();
      m_tpcseedcharge = tpcseed->get_qOverR() > 0 ? 1 : -1;
    }
    if (tpcseed)
    {
      m_dedx = calc_dedx(tpcseed, clustermap, tpcGeom);
    }

    if (tpcseed)
    {
      m_tpcseedpx = tpcseed->get_px();
      m_tpcseedpy = tpcseed->get_py();
      m_tpcseedpz = tpcseed->get_pz();
    }
    else if (silseed)
    {
      m_silseedpx = silseed->get_px();
      m_silseedpy = silseed->get_py();
      m_silseedpz = silseed->get_pz();
    }

    clearClusterStateVectors();
    if (Verbosity() > 1)
    {
      std::cout << "Track " << key << " has cluster/states"
                << std::endl;
    }

    // get the fully corrected cluster global positions
    std::vector<std::pair<TrkrDefs::cluskey, Acts::Vector3>> global_raw;
    for (const auto& ckey : get_cluster_keys(track))
    {
      auto *cluster = clustermap->findCluster(ckey);

      // Fully correct the cluster positions for the crossing and all distortions
      Acts::Vector3 global = m_globalPositionWrapper.getGlobalPositionDistortionCorrected(ckey, cluster, m_crossing);

      // add the global positions to a vector to give to the cluster mover
      global_raw.emplace_back(ckey, global);
    }

    // move the cluster positions back to the original readout surface in the fillClusterBranchesKF method

    if (!m_doAlignment)
    {
      for (const auto& ckey : get_cluster_keys(track))
      {
        fillClusterBranchesKF(ckey, track, global_raw, topNode);
      }
    }

    m_nhits = m_nmaps + m_nintt + m_ntpc + m_nmms;

    if (m_doAlignment)
    {
      /// repopulate with info that is going into alignment
      clearClusterStateVectors();

      if (alignmentmap and alignmentmap->find(key) != alignmentmap->end())
      {
        auto& statevec = alignmentmap->find(key)->second;

        for (auto& state : statevec)
        {
          auto ckey = state->get_cluster_key();

          fillClusterBranchesKF(ckey, track, global_raw, topNode);

          const auto& globderivs = state->get_global_derivative_matrix();
          const auto& locderivs = state->get_local_derivative_matrix();

          m_statelxglobderivdalpha.push_back(globderivs(0, 0));
          m_statelxglobderivdbeta.push_back(globderivs(0, 1));
          m_statelxglobderivdgamma.push_back(globderivs(0, 2));
          m_statelxglobderivdx.push_back(globderivs(0, 3));
          m_statelxglobderivdy.push_back(globderivs(0, 4));
          m_statelxglobderivdz.push_back(globderivs(0, 5));

          m_statelzglobderivdalpha.push_back(globderivs(1, 0));
          m_statelzglobderivdbeta.push_back(globderivs(1, 1));
          m_statelzglobderivdgamma.push_back(globderivs(1, 2));
          m_statelzglobderivdx.push_back(globderivs(1, 3));
          m_statelzglobderivdy.push_back(globderivs(1, 4));
          m_statelzglobderivdz.push_back(globderivs(1, 5));

          m_statelxlocderivd0.push_back(locderivs(0, 0));
          m_statelxlocderivz0.push_back(locderivs(0, 1));
          m_statelxlocderivphi.push_back(locderivs(0, 2));
          m_statelxlocderivtheta.push_back(locderivs(0, 3));
          m_statelxlocderivqop.push_back(locderivs(0, 4));

          m_statelzlocderivd0.push_back(locderivs(1, 0));
          m_statelzlocderivz0.push_back(locderivs(1, 1));
          m_statelzlocderivphi.push_back(locderivs(1, 2));
          m_statelzlocderivtheta.push_back(locderivs(1, 3));
          m_statelzlocderivqop.push_back(locderivs(1, 4));
        }
      }
    }

    if (m_nmms > 0 || !m_doMicromegasOnly)
    {
      m_tree->Fill();
    }

  }  // end loop over tracks

  if (m_doFailedSeeds)
  {
    fillFailedSeedTree(topNode, tpc_seed_ids);
  }
}
void TrackResiduals::fillEventTree(PHCompositeNode* topNode)
{
  auto *silseedmap = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
  auto *tpcseedmap = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  auto *trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
  auto *clustermap = findNode::getClass<TrkrClusterContainer>(topNode, m_clusterContainerName);
  auto *hitmap = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");

  m_ntpc_hits0 = 0;
  m_ntpc_hits1 = 0;
  m_ntpc_clus0 = 0;
  m_ntpc_clus1 = 0;
  m_nmvtx_all = 0;
  m_nintt_all = 0;
  m_nmms_all = 0;
  m_nsiseed = 0;
  m_ntpcseed = 0;
  m_ntpc_clus_sector.resize(24, 0);
  m_nsiseed = silseedmap->size();
  m_ntpcseed = tpcseedmap->size();
  m_ntracks_all = trackmap->size();

  // Hits
  if (m_doHits)
  {
    TrkrHitSetContainer::ConstRange tpc_hitsets = hitmap->getHitSets(TrkrDefs::TrkrId::tpcId);
    for (TrkrHitSetContainer::ConstIterator hitsetiter = tpc_hitsets.first;
         hitsetiter != tpc_hitsets.second;
         ++hitsetiter)
    {
      TrkrDefs::hitsetkey hitsetkey = hitsetiter->first;
      TrkrHitSet* hitset = hitsetiter->second;

      int thisside = TpcDefs::getSide(hitsetkey);
      if (thisside == 0)
      {
        m_ntpc_hits0 += hitset->size();
      }
      if (thisside == 1)
      {
        m_ntpc_hits1 += hitset->size();
      }
    }
  }
  for (const auto& det : {TrkrDefs::TrkrId::mvtxId, TrkrDefs::TrkrId::inttId,
                    TrkrDefs::TrkrId::tpcId, TrkrDefs::TrkrId::micromegasId})
  {
    for (const auto& hitsetkey : clustermap->getHitSetKeys(det))
    {
      auto range = clustermap->getClusters(hitsetkey);
      int nclus = std::distance(range.first, range.second);
      int tpcside = TrkrDefs::getZElement(hitsetkey);
      int sector = TpcDefs::getSectorId(hitsetkey);

      switch (det)
      {
      case TrkrDefs::TrkrId::mvtxId:
        m_nmvtx_all += nclus;
        break;
      case TrkrDefs::TrkrId::inttId:
        m_nintt_all += nclus;
        break;
      case TrkrDefs::TrkrId::tpcId:
        if(tpcside == 1)
        {
          sector += 12;
        }
        m_ntpc_clus_sector[sector] += nclus;
        if (tpcside == 0)
        {
          m_ntpc_clus0 += nclus;
        }
        if (tpcside == 1)
        {
          m_ntpc_clus1 += nclus;
        }
        break;
      case TrkrDefs::TrkrId::micromegasId:
        m_nmms_all += nclus;
        break;
      default:
        break;
      }
    }
  }
  if (m_doEventTree)
  {
    if (Verbosity() > 1)
    {
      std::cout << " m_event:" << m_event << std::endl;
      std::cout << " m_ntpc_clus0:" << m_ntpc_clus0 << std::endl;
      std::cout << " m_ntpc_clus1: " << m_ntpc_clus1 << std::endl;
      std::cout << " m_nmvtx_all:" << m_nmvtx_all << std::endl;
      std::cout << " m_nintt_all: " << m_nintt_all << std::endl;
      std::cout << " m_nmms_all: " << m_nmms_all << std::endl;
    }
    m_eventtree->Fill();
  }
}

void TrackResiduals::fillResidualTreeSeeds(PHCompositeNode* topNode)
{
  auto *silseedmap = findNode::getClass<TrackSeedContainer>(topNode, "SiliconTrackSeedContainer");
  auto *tpcseedmap = findNode::getClass<TrackSeedContainer>(topNode, "TpcTrackSeedContainer");
  auto *tpcGeom =
      findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  auto *trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
  auto *clustermap = findNode::getClass<TrkrClusterContainer>(topNode, m_clusterContainerName);
  auto *vertexmap = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  auto *alignmentmap = findNode::getClass<SvtxAlignmentStateMap>(topNode, m_alignmentMapName);

  std::set<unsigned int> tpc_seed_ids;

  for (const auto& [key, track] : *trackmap)
  {
    if (!track)
    {
      continue;
    }
    m_trackid = track->get_id();

    m_crossing = track->get_crossing();
    m_crossing_estimate = SHRT_MAX;
    m_px = track->get_px();
    m_py = track->get_py();
    m_pz = track->get_pz();

    m_pt = std::sqrt(square(m_px) + square(m_py));
    m_eta = std::atanh(m_pz / std::sqrt(square(m_pt) + square(m_pz)));
    m_phi = std::atan2(m_py, m_px);
    float CVxx = track->get_error(3, 3);
    float CVxy = track->get_error(3, 4);
    float CVyy = track->get_error(4, 4);
    m_deltapt = std::sqrt((CVxx * square(m_px) + 2 * CVxy * m_px * m_py + CVyy * square(m_py)) / (square(m_px) + square(m_py)));

    m_charge = track->get_charge();
    m_quality = track->get_quality();
    m_chisq = track->get_chisq();
    m_ndf = track->get_ndf();

    if (Verbosity() > 1)
    {
      std::cout << "fillResidualTreeSeeds:  track " << m_trackid << " m_crossing " << m_crossing << " m_crossing_estimate " << m_crossing_estimate << " m_pt " << m_pt << std::endl;
    }

    m_nmaps = 0;
    m_nintt = 0;
    m_ntpc = 0;
    m_nmms = 0;
    m_silid = std::numeric_limits<unsigned int>::quiet_NaN();
    m_tpcid = std::numeric_limits<unsigned int>::quiet_NaN();
    m_silseedx = std::numeric_limits<float>::quiet_NaN();
    m_silseedy = std::numeric_limits<float>::quiet_NaN();
    m_silseedz = std::numeric_limits<float>::quiet_NaN();
    m_silseedpx = std::numeric_limits<float>::quiet_NaN();
    m_silseedpy = std::numeric_limits<float>::quiet_NaN();
    m_silseedpz = std::numeric_limits<float>::quiet_NaN();
    m_silseedphi = std::numeric_limits<float>::quiet_NaN();
    m_silseedeta = std::numeric_limits<float>::quiet_NaN();
    m_silseedcharge = std::numeric_limits<int>::quiet_NaN();
    m_tpcseedx = std::numeric_limits<float>::quiet_NaN();
    m_tpcseedy = std::numeric_limits<float>::quiet_NaN();
    m_tpcseedz = std::numeric_limits<float>::quiet_NaN();
    m_tpcseedpx = std::numeric_limits<float>::quiet_NaN();
    m_tpcseedpy = std::numeric_limits<float>::quiet_NaN();
    m_tpcseedpz = std::numeric_limits<float>::quiet_NaN();
    m_tpcseedphi = std::numeric_limits<float>::quiet_NaN();
    m_tpcseedeta = std::numeric_limits<float>::quiet_NaN();
    m_tpcseedcharge = std::numeric_limits<int>::quiet_NaN();

    m_vertexid = track->get_vertex_id();

    m_pcax = track->get_x();
    m_pcay = track->get_y();
    m_pcaz = track->get_z();
    m_dcaxy = std::numeric_limits<float>::quiet_NaN();
    m_dcaz = std::numeric_limits<float>::quiet_NaN();
    if (vertexmap)
    {
      auto vertexit = vertexmap->find(m_vertexid);
      if (vertexit != vertexmap->end())
      {
        auto *vertex = vertexit->second;
        m_vx = vertex->get_x();
        m_vy = vertex->get_y();
        m_vz = vertex->get_z();
        Acts::Vector3 vert(m_vx, m_vy, m_vz);
        auto dcapair = TrackAnalysisUtils::get_dca(track, vert);
        m_dcaxy = dcapair.first.first;
        m_dcaz = dcapair.second.first;
      }
    }

    auto *tpcseed = track->get_tpc_seed();
    if (tpcseed)
    {
      m_tpcid = tpcseedmap->find(tpcseed);
      tpc_seed_ids.insert(tpcseedmap->find(tpcseed));
    }
    auto *silseed = track->get_silicon_seed();
    if (silseed)
    {
      m_silid = silseedmap->find(silseed);
      const auto si_pos = TrackSeedHelper::get_xyz(silseed);
      m_silseedx = si_pos.x();
      m_silseedy = si_pos.y();
      m_silseedz = si_pos.z();
      m_silseedpx = silseed->get_px();
      m_silseedpy = silseed->get_py();
      m_silseedpz = silseed->get_pz();
      m_silseedphi = silseed->get_phi();
      m_silseedeta = silseed->get_eta();
      m_silseedcharge = silseed->get_qOverR() > 0 ? 1 : -1;
    }
    if (tpcseed)
    {
      const auto tpc_pos = TrackSeedHelper::get_xyz(tpcseed);
      m_tpcseedx = tpc_pos.x();
      m_tpcseedy = tpc_pos.y();
      m_tpcseedz = tpc_pos.z();
      m_tpcseedpx = tpcseed->get_px();
      m_tpcseedpy = tpcseed->get_py();
      m_tpcseedpz = tpcseed->get_pz();
      m_tpcseedphi = tpcseed->get_phi();
      m_tpcseedeta = tpcseed->get_eta();
      m_tpcseedcharge = tpcseed->get_qOverR() > 0 ? 1 : -1;
    }
    if (tpcseed)
    {
      m_dedx = calc_dedx(tpcseed, clustermap, tpcGeom);
    }
    if (m_zeroField)
    {
      float qor = std::numeric_limits<float>::quiet_NaN();
      float phi = std::numeric_limits<float>::quiet_NaN();
      float theta = std::numeric_limits<float>::quiet_NaN();
      float eta = std::numeric_limits<float>::quiet_NaN();
      if (tpcseed)
      {
        qor = tpcseed->get_qOverR();
        phi = tpcseed->get_phi();
        eta = tpcseed->get_eta();
        theta = tpcseed->get_theta();
      }
      else if (silseed)
      {
        qor = silseed->get_qOverR();
        phi = silseed->get_phi();
        eta = silseed->get_eta();
        theta = silseed->get_theta();
      }
      float pt = fabs(1. / qor) * (0.3 / 100) * 0.01;
      m_tpcseedpx = pt * std::cos(phi);
      m_tpcseedpy = pt * std::sin(phi);
      m_tpcseedpz = pt * std::cosh(eta) * std::cos(theta);
    }
    else
    {
      if (tpcseed)
      {
        m_tpcseedpx = tpcseed->get_px();
        m_tpcseedpy = tpcseed->get_py();
        m_tpcseedpz = tpcseed->get_pz();
      }
      else if (silseed)
      {
        m_silseedpx = silseed->get_px();
        m_silseedpy = silseed->get_py();
        m_silseedpz = silseed->get_pz();
      }
    }
    clearClusterStateVectors();
    if (Verbosity() > 1)
    {
      std::cout << "Track " << key << " has cluster/states"
                << std::endl;
    }

    // get the fully corrected cluster global positions
    std::vector<std::pair<TrkrDefs::cluskey, Acts::Vector3>> global_raw;
    float minR = std::numeric_limits<float>::max();
    float maxR = 0;
    for (const auto& ckey : get_cluster_keys(track))
    {
      auto *cluster = clustermap->findCluster(ckey);

      // Fully correct the cluster positions for the crossing and all distortions
      Acts::Vector3 global = m_globalPositionWrapper.getGlobalPositionDistortionCorrected(ckey, cluster, m_crossing);

      // add the global positions to a vector to give to the cluster mover
      global_raw.emplace_back(ckey, global);
      minR = std::min<double>(r(global.x(), global.y()), minR);
      maxR = std::max<double>(r(global.x(), global.y()), maxR);
    }
    m_tracklength = maxR - minR;
    // ---- we move the global positions back to the surface in fillClusterBranchesSeeds

    if (!m_doAlignment)
    {
      std::vector<TrkrDefs::cluskey> keys;
      for (const auto& ckey : get_cluster_keys(track))
      {
        keys.push_back(ckey);
      }
      if (m_zeroField)
      {
        lineFitClusters(keys, clustermap, m_crossing);
      }
      else
      {
        // this corrects the cluster positions and fits them, to fill the helical fit parameters
        //  that are used to calculate the "state" positions
        circleFitClusters(keys, clustermap, m_crossing);
      }

      for (const auto& ckey : get_cluster_keys(track))
      {
        fillClusterBranchesSeeds(ckey, global_raw, topNode);
      }
    }

    m_nhits = m_nmaps + m_nintt + m_ntpc + m_nmms;

    if (m_doAlignment)
    {
      /// repopulate with info that is going into alignment
      clearClusterStateVectors();

      if (alignmentmap and alignmentmap->find(key) != alignmentmap->end())
      {
        auto& statevec = alignmentmap->find(key)->second;

        for (auto& state : statevec)
        {
          auto ckey = state->get_cluster_key();

          fillClusterBranchesSeeds(ckey, global_raw, topNode);

          const auto& globderivs = state->get_global_derivative_matrix();
          const auto& locderivs = state->get_local_derivative_matrix();

          m_statelxglobderivdalpha.push_back(globderivs(0, 0));
          m_statelxglobderivdbeta.push_back(globderivs(0, 1));
          m_statelxglobderivdgamma.push_back(globderivs(0, 2));
          m_statelxglobderivdx.push_back(globderivs(0, 3));
          m_statelxglobderivdy.push_back(globderivs(0, 4));
          m_statelxglobderivdz.push_back(globderivs(0, 5));

          m_statelzglobderivdalpha.push_back(globderivs(1, 0));
          m_statelzglobderivdbeta.push_back(globderivs(1, 1));
          m_statelzglobderivdgamma.push_back(globderivs(1, 2));
          m_statelzglobderivdx.push_back(globderivs(1, 3));
          m_statelzglobderivdy.push_back(globderivs(1, 4));
          m_statelzglobderivdz.push_back(globderivs(1, 5));

          m_statelxlocderivd0.push_back(locderivs(0, 0));
          m_statelxlocderivz0.push_back(locderivs(0, 1));
          m_statelxlocderivphi.push_back(locderivs(0, 2));
          m_statelxlocderivtheta.push_back(locderivs(0, 3));
          m_statelxlocderivqop.push_back(locderivs(0, 4));

          m_statelzlocderivd0.push_back(locderivs(1, 0));
          m_statelzlocderivz0.push_back(locderivs(1, 1));
          m_statelzlocderivphi.push_back(locderivs(1, 2));
          m_statelzlocderivtheta.push_back(locderivs(1, 3));
          m_statelzlocderivqop.push_back(locderivs(1, 4));
        }
      }
    }
    if (!m_doMatchedOnly)
    {
      m_tree->Fill();
    }
    else
    {
      if (m_nmaps >= 3 && m_nintt >= 1 && m_ntpc > 32 && abs(m_crossing) < 5 && m_pt > 0.3)
      {
        m_tree->Fill();
      }
    }
  }
}
