
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

#include <intt/CylinderGeomIntt.h>

#include <micromegas/CylinderGeomMicromegas.h>
#include <micromegas/MicromegasDefs.h>

#include <mvtx/CylinderGeom_Mvtx.h>

#include <trackbase_historic/ActsTransformations.h>
#include <trackbase_historic/SvtxAlignmentState.h>
#include <trackbase_historic/SvtxAlignmentStateMap.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackSeed.h>

#include <ffarawobjects/Gl1RawHit.h>

#include <globalvertex/GlobalVertex.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>

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
TrackResiduals::~TrackResiduals()
{
}

//____________________________________________________________________________..
int TrackResiduals::Init(PHCompositeNode*)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int TrackResiduals::InitRun(PHCompositeNode*)
{
  m_outfile = new TFile(m_outfileName.c_str(), "RECREATE");
  createBranches();

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
  m_clusgx.clear();
  m_clusgy.clear();
  m_clusgz.clear();
  m_cluslayer.clear();
  m_clussize.clear();
  m_clushitsetkey.clear();

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
}
//____________________________________________________________________________..
int TrackResiduals::process_event(PHCompositeNode* topNode)
{
  auto trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
  auto clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  auto geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  auto vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  auto alignmentmap = findNode::getClass<SvtxAlignmentStateMap>(topNode, m_alignmentMapName);
  auto hitmap = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  auto tpcGeom =
      findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  auto mvtxGeom = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");
  auto inttGeom = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
  auto mmGeom = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MICROMEGAS_FULL");
  if (!mmGeom)
  {
    mmGeom = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MICROMEGAS");
  }
  if (!trackmap or !clustermap or !geometry or !hitmap)
  {
    std::cout << "Missing node, can't continue" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  auto gl1 = findNode::getClass<Gl1RawHit>(topNode, "GL1RAWHIT");
  if (gl1)
  {
    m_bco = gl1->get_bco();
    auto lbshift = m_bco << 24;
    m_bcotr = lbshift >> 24;
  }
  else
  {
    m_bco = std::numeric_limits<uint64_t>::quiet_NaN();
    m_bcotr = std::numeric_limits<uint64_t>::quiet_NaN();
  }
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

  for (const auto& [key, track] : *trackmap)
  {
    if (!track)
    {
      continue;
    }

    m_trackid = key;
    m_crossing = track->get_crossing();
    m_px = track->get_px();
    m_py = track->get_py();
    m_pz = track->get_pz();
    m_pt = std::sqrt(square(m_px) + square(m_py));
    m_eta = atanh(m_pz / std::sqrt(square(m_pt) + square(m_pz)));
    m_phi = atan2(m_py, m_px);
    float CVxx = track->get_error(3, 3);
    float CVxy = track->get_error(3, 4);
    float CVyy = track->get_error(4, 4);
    m_deltapt = std::sqrt((CVxx * square(m_px) + 2 * CVxy * m_px * m_py + CVyy * square(m_py)) / (square(m_px) + square(m_py)));

    m_charge = track->get_charge();
    m_quality = track->get_quality();
    m_chisq = track->get_chisq();
    m_ndf = track->get_ndf();
    m_nmaps = 0;
    m_nintt = 0;
    m_ntpc = 0;
    m_nmms = 0;
    m_vertexid = track->get_vertex_id();
    if (vertexmap)
    {
      auto vertexit = vertexmap->find(m_vertexid);
      if (vertexit != vertexmap->end())
      {
        auto vertex = vertexit->second;
        m_vx = vertex->get_x();
        m_vy = vertex->get_y();
        m_vz = vertex->get_z();
      }
    }
    m_pcax = track->get_x();
    m_pcay = track->get_y();
    m_pcaz = track->get_z();

    clearClusterStateVectors();
    if (Verbosity() > 1)
    {
      std::cout << "Track " << key << " has cluster/states"
                << std::endl;
    }

    if (!m_doAlignment)
    {
      std::vector<TrkrDefs::cluskey> keys;
      for (const auto& ckey : get_cluster_keys(track))
      {
        if (TrkrDefs::getTrkrId(ckey) == TrkrDefs::TrkrId::tpcId)
        {
          keys.push_back(ckey);
        }
      }
      if (m_zeroField)
      {
        lineFitClusters(keys, geometry, clustermap);
      }
      else
      {
        circleFitClusters(keys, geometry, clustermap);
      }
      for (const auto& ckey : get_cluster_keys(track))
      {
        fillClusterBranches(ckey, track, topNode);
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

          fillClusterBranches(ckey, track, topNode);

          auto& globderivs = state->get_global_derivative_matrix();
          auto& locderivs = state->get_local_derivative_matrix();

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
    m_tree->Fill();
  }

  m_event++;
  return Fun4AllReturnCodes::EVENT_OK;
}

float TrackResiduals::convertTimeToZ(ActsGeometry* geometry, TrkrDefs::cluskey cluster_key, TrkrCluster* cluster)
{
  // must convert local Y from cluster average time of arival to local cluster z position
  double drift_velocity = geometry->get_drift_velocity();
  double zdriftlength = cluster->getLocalY() * drift_velocity;
  double surfCenterZ = 52.89;                // 52.89 is where G4 thinks the surface center is
  double zloc = surfCenterZ - zdriftlength;  // converts z drift length to local z position in the TPC in north
  unsigned int side = TpcDefs::getSide(cluster_key);
  if (side == 0) zloc = -zloc;
  float z = zloc;  // in cm

  return z;
}
void TrackResiduals::circleFitClusters(std::vector<TrkrDefs::cluskey>& keys,
                                       ActsGeometry* geometry,
                                       TrkrClusterContainer* clusters)
{
  std::vector<Acts::Vector3> clusPos, global_vec;
  TrackFitUtils::getTrackletClusters(geometry, clusters,
                                     clusPos, keys);

  for (auto& pos : clusPos)
  {
    float clusr = r(pos.x(), pos.y());
    if (pos.y() < 0) clusr *= -1;

    // exclude silicon and tpot clusters for now
    if (fabs(clusr) > 80 || fabs(clusr) < 30)
    {
      continue;
    }
    global_vec.push_back(pos);
  }

  auto fitpars = TrackFitUtils::fitClusters(global_vec, keys);

  m_xyint = std::numeric_limits<float>::quiet_NaN();
  m_xyslope = std::numeric_limits<float>::quiet_NaN();
  m_R = fitpars[0];
  m_X0 = fitpars[1];
  m_Y0 = fitpars[2];
  m_rzslope = fitpars[3];
  m_rzint = fitpars[4];
}

void TrackResiduals::lineFitClusters(std::vector<TrkrDefs::cluskey>& keys,
                                     ActsGeometry* geometry,
                                     TrkrClusterContainer* clusters)
{
  std::vector<Acts::Vector3> clusPos;
  TrackFitUtils::getTrackletClusters(geometry, clusters,
                                     clusPos, keys);
  TrackFitUtils::position_vector_t xypoints, rzpoints;
  for (auto& pos : clusPos)
  {
    float clusr = r(pos.x(), pos.y());
    if (pos.y() < 0) clusr *= -1;

    // exclude silicon and tpot clusters for now
    if (fabs(clusr) > 80 || fabs(clusr) < 30)
    {
      continue;
    }
    rzpoints.push_back(std::make_pair(pos.z(), clusr));
    xypoints.push_back(std::make_pair(pos.x(), pos.y()));
  }

  auto xyparams = TrackFitUtils::line_fit(xypoints);
  auto rzparams = TrackFitUtils::line_fit(rzpoints);
  m_xyint = std::get<1>(xyparams);
  m_xyslope = std::get<0>(xyparams);
  m_rzint = std::get<1>(rzparams);
  m_rzslope = std::get<0>(rzparams);
}

void TrackResiduals::fillClusterTree(TrkrClusterContainer* clusters,
                                     ActsGeometry* geometry)
{
  for (auto& det : {TrkrDefs::TrkrId::mvtxId, TrkrDefs::TrkrId::inttId,
                    TrkrDefs::TrkrId::tpcId, TrkrDefs::TrkrId::micromegasId})
  {
    for (const auto& hitsetkey : clusters->getHitSetKeys(det))
    {
      m_scluslayer = TrkrDefs::getLayer(hitsetkey);
      auto range = clusters->getClusters(hitsetkey);
      for (auto iter = range.first; iter != range.second; ++iter)
      {
        auto key = iter->first;
        auto cluster = clusters->findCluster(key);
        auto glob = geometry->getGlobalPosition(key, cluster);
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
int TrackResiduals::End(PHCompositeNode*)
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
      auto hit = hitr->second;
      m_adc = hit->getAdc();

      switch (det)
      {
      case TrkrDefs::TrkrId::mvtxId:
      {
        m_row = MvtxDefs::getRow(hitkey);
        m_col = MvtxDefs::getCol(hitkey);
        auto layergeom = dynamic_cast<CylinderGeom_Mvtx*>(mvtxGeom->GetLayerGeom(m_hitlayer));
        auto local_coords = layergeom->get_local_coords_from_pixel(m_row, m_col);
        TVector2 local;
        local.SetX(local_coords.X());
        local.SetY(local_coords.Z());
        auto surf = geometry->maps().getSiliconSurface(m_hitsetkey);
        auto glob = layergeom->get_world_from_local_coords(surf, geometry, local);
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
        auto geom = dynamic_cast<CylinderGeomIntt*>(inttGeom->GetLayerGeom(m_hitlayer));
        double local_hit_loc[3] = {0, 0, 0};
        geom->find_strip_center_localcoords(m_ladderzid, m_row, m_col, local_hit_loc);
        auto surf = geometry->maps().getSiliconSurface(m_hitsetkey);
        TVector2 local;
        local.SetX(local_hit_loc[1]);
        local.SetY(local_hit_loc[2]);
        auto glob = geom->get_world_from_local_coords(surf, geometry, local);

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

        auto geoLayer = tpcGeom->GetLayerCellGeom(m_hitlayer);
        auto phi = geoLayer->get_phicenter(m_hitpad);
        auto radius = geoLayer->get_radius();

        float AdcClockPeriod = geoLayer->get_zstep();
        m_zdriftlength = m_hittbin * geometry->get_drift_velocity() * AdcClockPeriod;
	double NZBinsSide = 249;  // physical z bins per TPC side
	double tdriftmax = AdcClockPeriod * NZBinsSide;
        m_hitgz = (tdriftmax * geometry->get_drift_velocity()) - m_zdriftlength;
        if (m_side == 0)
        {
          m_hitgz *= -1;
        }
        m_hitgx = radius * std::cos(phi);
        m_hitgy = radius * std::sin(phi);
        break;
      }
      case TrkrDefs::TrkrId::micromegasId:
      {
        const auto layergeom = dynamic_cast<CylinderGeomMicromegas*>(mmGeom->GetLayerGeom(m_hitlayer));
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

void TrackResiduals::fillClusterBranches(TrkrDefs::cluskey ckey, SvtxTrack* track,
                                         PHCompositeNode* topNode)
{
  auto clustermap = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  auto geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");

  ActsTransformations transformer;
  TrkrCluster* cluster = clustermap->findCluster(ckey);
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
    break;
  case TrkrDefs::micromegasId:
    m_nmms++;
    break;
  }

  Acts::Vector3 clusglob = geometry->getGlobalPosition(ckey, cluster);

  SvtxTrackState* state = nullptr;

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

  m_cluskeys.push_back(ckey);

  //! have cluster and state, fill vectors
  m_cluslx.push_back(cluster->getLocalX());
  m_clusedge.push_back(cluster->getEdge());
  m_clusoverlap.push_back(cluster->getOverlap());
  float clusz = cluster->getLocalY();

  if (TrkrDefs::getTrkrId(ckey) == TrkrDefs::TrkrId::tpcId)
  {
    float rawclusz = convertTimeToZ(geometry, ckey, cluster);
    
    int crossing = track->get_crossing();
    unsigned int side = TpcDefs::getSide(ckey);
    clusz = m_clusterCrossingCorrection.correctZ(rawclusz, side, crossing);
    if(!m_ppmode)
    {
      clusz = rawclusz;
    }
  }

  m_cluslz.push_back(clusz);
  float clusr = r(clusglob.x(), clusglob.y());
  auto para_errors = m_clusErrPara.get_clusterv5_modified_error(cluster,
                                                                clusr, ckey);
  m_cluselx.push_back(sqrt(para_errors.first));
  m_cluselz.push_back(sqrt(para_errors.second));
  m_clusgx.push_back(clusglob.x());
  m_clusgy.push_back(clusglob.y());
  m_clusgz.push_back(clusglob.z());
  m_cluslayer.push_back(TrkrDefs::getLayer(ckey));
  m_clusphisize.push_back(cluster->getPhiSize());
  m_cluszsize.push_back(cluster->getZSize());
  m_clussize.push_back(cluster->getPhiSize() * cluster->getZSize());
  m_clushitsetkey.push_back(TrkrDefs::getHitSetKeyFromClusKey(ckey));

  if (Verbosity() > 1)
  {
    std::cout << "Track state/clus in layer "
              << (unsigned int) TrkrDefs::getLayer(ckey) << " with pos "
              << clusglob.transpose() << std::endl;
  }
  if (!state)
  {
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
    m_idealsurfalpha.push_back(NAN);
    m_idealsurfbeta.push_back(NAN);
    m_idealsurfgamma.push_back(NAN);
    m_missurfalpha.push_back(NAN);
    m_missurfbeta.push_back(NAN);
    m_missurfgamma.push_back(NAN);
    m_idealsurfcenterx.push_back(NAN);
    m_idealsurfcentery.push_back(NAN);
    m_idealsurfcenterz.push_back(NAN);
    m_idealsurfnormx.push_back(NAN);
    m_idealsurfnormy.push_back(NAN);
    m_idealsurfnormz.push_back(NAN);
    m_missurfcenterx.push_back(NAN);
    m_missurfcentery.push_back(NAN);
    m_missurfcenterz.push_back(NAN);
    m_missurfnormx.push_back(NAN);
    m_missurfnormy.push_back(NAN);
    m_missurfnormz.push_back(NAN);
    m_clusgxideal.push_back(NAN);
    m_clusgyideal.push_back(NAN);
    m_clusgzideal.push_back(NAN);
    m_statepx.push_back(NAN);
    m_statepy.push_back(NAN);
    m_statepz.push_back(NAN);
    m_statepl.push_back(NAN);
    return;
  }

  auto surf = geometry->maps().getSurface(ckey, cluster);
  Acts::Vector3 stateglob(state->get_x(), state->get_y(), state->get_z());
  Acts::Vector2 stateloc;
  auto misaligncenter = surf->center(geometry->geometry().getGeoContext());
  auto misalignnorm = -1 * surf->normal(geometry->geometry().getGeoContext());
  auto misrot = surf->transform(geometry->geometry().getGeoContext()).rotation();
  auto result = surf->globalToLocal(geometry->geometry().getGeoContext(),
                                    stateglob * Acts::UnitConstants::cm,
                                    misalignnorm);

  float mgamma = atan2(-misrot(1, 0), misrot(0, 0));
  float mbeta = -asin(misrot(0, 1));
  float malpha = atan2(misrot(1, 1), misrot(2, 1));

  //! Switch to get ideal transforms
  alignmentTransformationContainer::use_alignment = false;
  auto idealcenter = surf->center(geometry->geometry().getGeoContext());
  auto idealnorm = -1 * surf->normal(geometry->geometry().getGeoContext());
  Acts::Vector3 ideal_local(cluster->getLocalX(), clusz, 0.0);
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

  if (result.ok())
  {
    stateloc = result.value() / Acts::UnitConstants::cm;
  }
  else
  {
    //! manual transform for tpc
    Acts::Vector3 loct = surf->transform(geometry->geometry().getGeoContext()).inverse() * (stateglob * Acts::UnitConstants::cm);
    loct /= Acts::UnitConstants::cm;
    stateloc(0) = loct(0);
    stateloc(1) = loct(1);
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
  m_xyint, m_rzslope, m_rzint);
  
  auto surf = geometry->maps().getSurface(key, cluster);
  Acts::Vector3 surfnorm = surf->normal(geometry->geometry().getGeoContext());

  if(!std::isnan(intersection.x()))
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
    m_statelx.push_back(NAN);
    m_statelz.push_back(NAN);
    m_stategx.push_back(NAN);
    m_stategy.push_back(NAN);
    m_stategz.push_back(NAN);
  }
}
void TrackResiduals::createBranches()
{
  m_hittree = new TTree("hittree", "A tree with all hits");
  m_hittree->Branch("run", &m_runnumber, "m_runnumber/I");
  m_hittree->Branch("segment", &m_segment, "m_segment/I");
  m_hittree->Branch("event", &m_event, "m_event/I");
  m_hittree->Branch("gl1bco", &m_bco, "m_bco/l");
  m_hittree->Branch("trbco", &m_bcotr, "m_bcotr/l");
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

  m_clustree = new TTree("clustertree", "A tree with all clusters");
  m_clustree->Branch("run", &m_runnumber, "m_runnumber/I");
  m_clustree->Branch("segment", &m_segment, "m_segment/I");
  m_clustree->Branch("event", &m_event, "m_event/I");
  m_clustree->Branch("gl1bco", &m_bco, "m_bco/l");
  m_clustree->Branch("trbco", &m_bcotr, "m_bcotr/l");
  m_clustree->Branch("lx", &m_scluslx, "m_scluslx/F");
  m_clustree->Branch("lz", &m_scluslz, "m_scluslz/F");
  m_clustree->Branch("gx", &m_sclusgx, "m_sclusgx/F");
  m_clustree->Branch("gy", &m_sclusgy, "m_sclusgy/F");
  m_clustree->Branch("gz", &m_sclusgz, "m_sclusgz/F");
  m_clustree->Branch("r", &m_sclusgr, "m_sclusgr/F");
  m_clustree->Branch("phi", &m_sclusphi, "m_sclusphi/F");
  m_clustree->Branch("eta", &m_scluseta, "m_scluseta/F");
  m_clustree->Branch("adc", &m_adc, "m_adc/F");
  m_clustree->Branch("phisize", &m_phisize, "m_phisize/I");
  m_clustree->Branch("zsize", &m_zsize, "m_zsize/I");
  m_clustree->Branch("layer", &m_scluslayer, "m_scluslayer/I");
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
  m_tree->Branch("trackid", &m_trackid, "m_trackid/I");
  m_tree->Branch("gl1bco", &m_bco, "m_bco/l");
  m_tree->Branch("trbco", &m_bcotr, "m_bcotr/l");
  m_tree->Branch("crossing", &m_crossing, "m_crossing/I");
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
  m_tree->Branch("nintt", &m_nintt, "m_nintt/I");
  m_tree->Branch("ntpc", &m_ntpc, "m_ntpc/I");
  m_tree->Branch("nmms", &m_nmms, "m_nmms/I");
  m_tree->Branch("vertexid", &m_vertexid, "m_vertexid/I");
  m_tree->Branch("vx", &m_vx, "m_vx/F");
  m_tree->Branch("vy", &m_vy, "m_vy/F");
  m_tree->Branch("vz", &m_vz, "m_vz/F");
  m_tree->Branch("pcax", &m_pcax, "m_pcax/F");
  m_tree->Branch("pcay", &m_pcay, "m_pcay/F");
  m_tree->Branch("pcaz", &m_pcaz, "m_pcaz/F");
  m_tree->Branch("rzslope", &m_rzslope, "m_rzslope/F");
  m_tree->Branch("xyslope", &m_xyslope, "m_xyslope/F");
  m_tree->Branch("rzint", &m_rzint, "m_rzint/F");
  m_tree->Branch("xyint", &m_xyint, "m_xyint/F");
  m_tree->Branch("R", &m_R, "m_R/F");
  m_tree->Branch("X0", &m_X0, "m_X0/F");
  m_tree->Branch("Y0", &m_Y0, "m_Y0/F");

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
  m_tree->Branch("cluslayer", &m_cluslayer);
  m_tree->Branch("clussize", &m_clussize);
  m_tree->Branch("clusphisize",&m_clusphisize);
  m_tree->Branch("cluszsize",&m_cluszsize);
  m_tree->Branch("clushitsetkey", &m_clushitsetkey);
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
