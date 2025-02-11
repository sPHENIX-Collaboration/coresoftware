#include "g4evaltools.h"

#include <mvtx/CylinderGeom_Mvtx.h>

#include <intt/CylinderGeomIntt.h>

#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <g4tracking/EmbRecoMatchContainer.h>
#include <g4tracking/TrkrTruthTrack.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <TFile.h>
#include <TObjString.h>

#include <iostream>
#include <numeric> // for std::accumulate

namespace G4Eval
{
  std::vector<int> unmatchedSvtxTrkIds(EmbRecoMatchContainer* matches, SvtxTrackMap* m_SvtxTrackMap)
  {
    std::set<unsigned int> ids_matched{};
    for (auto trkid : matches->ids_RecoMatched())
    {
      ids_matched.insert(trkid);
    }

    std::set<int> ids_unmatched;
    for (auto& reco : *m_SvtxTrackMap)
    {
      auto trkid = reco.first;
      if (ids_matched.count(trkid) == 0)
      {
        ids_unmatched.insert(trkid);
      }
    }
    std::vector<int> ids_vec;
    ids_vec.reserve(ids_unmatched.size());
    for (auto id : ids_unmatched)
    {
      ids_vec.push_back((int) id);
    }
    std::sort(ids_vec.begin(), ids_vec.end());
    return ids_vec;
  }

  // function implementation mostly from
  // https://root-forum.cern.ch/t/is-it-possible-to-save-plain-text-in-a-root-file/27674/4
  // Note that there is also the code there to read it back from a TFile
  void write_StringToTFile(const std::string& msg_name, const std::string& msg)
  {
    // the string is either written to the current file, or to a new file
    // that is named f_outname
    TFile* s_current = gDirectory->GetFile();
    if (s_current == nullptr)
    {
      std::cout << PHWHERE << " Error no TFile open to which to wrote the "
                << std::endl
                << " TObjString mesaged." << std::endl;
      return;
    }
    TObjString obj(msg.c_str());
    s_current->WriteObject(&obj, msg_name.c_str());
    return;
  }

  // return the layer of the hit: 0:Mvtx 1:Intt 2:Tpc 3:Tpot
  int trklayer_0123(TrkrDefs::hitsetkey key)
  {
    auto layer = TrkrDefs::getLayer(key);
    if (layer < 3)
    {
      return 0;
    }
    if (layer < 7)
    {
      return 1;
    }
    if (layer < 55)
    {
      return 2;
    }
    return 3;
  }

  // Implementation of Cluster comparator
  TrkrClusterComparer::TrkrClusterComparer(float _nphi_widths, float _nz_widths)
    : m_nphi_widths{_nphi_widths}
    , m_nz_widths{_nz_widths} {};
  int TrkrClusterComparer::init(PHCompositeNode* topNode,
                                const std::string& name_phg4_clusters,
                                const std::string& name_reco_clusters)
  {
    // fill bin/pixel sizes
    // ------ MVTX data ------
    PHG4CylinderGeomContainer* geom_container_mvtx = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");
    if (!geom_container_mvtx)
    {
      std::cout << PHWHERE << " Could not locate CYLINDERGEOM_MVTX " << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    for (int this_layer = 0; this_layer < 3; ++this_layer)
    {
      auto layergeom = dynamic_cast<CylinderGeom_Mvtx*>(geom_container_mvtx->GetLayerGeom(this_layer));
      const double pitch = layergeom->get_pixel_x();
      const double length = layergeom->get_pixel_z();
      m_phistep[this_layer] = pitch;
      if (this_layer == 0)
      {
        m_zstep_mvtx = length;
      }
    }

    // ------ INTT data ------
    PHG4CylinderGeomContainer* geom_container_intt = findNode::getClass<PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
    if (!geom_container_intt)
    {
      std::cout << PHWHERE << " Could not locate CYLINDERGEOM_INTT " << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    // get phi and Z steps for intt
    for (int this_layer = 3; this_layer < 7; ++this_layer)
    {
      CylinderGeomIntt* geom =
          dynamic_cast<CylinderGeomIntt*>(geom_container_intt->GetLayerGeom(this_layer));
      float pitch = geom->get_strip_y_spacing();
      m_phistep[this_layer] = pitch;
    }

    // ------ TPC data ------
    auto geom_tpc =
        findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
    if (!geom_tpc)
    {
      std::cout << PHWHERE << " Could not locate CYLINDERCELLGEOM_SVTX node " << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    for (int this_layer = 7; this_layer < 55; ++this_layer)
    {
      PHG4TpcCylinderGeom* layergeom = geom_tpc->GetLayerCellGeom(this_layer);
      if (this_layer == 7)
      {
        m_zstep_tpc = layergeom->get_zstep();
      }
      m_phistep[this_layer] = layergeom->get_phistep();
    }

    m_TruthClusters =
        findNode::getClass<TrkrClusterContainer>(topNode, name_phg4_clusters.c_str());
    if (!m_TruthClusters)
    {
      std::cout << PHWHERE << " Could not locate " << name_phg4_clusters << " node" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    m_RecoClusters =
        findNode::getClass<TrkrClusterContainer>(topNode, name_reco_clusters.c_str());
    if (!m_TruthClusters)
    {
      std::cout << PHWHERE << " Could not locate " << name_reco_clusters << " node" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    m_ActsGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
    if (!m_ActsGeometry)
    {
      std::cout << PHWHERE << " Could not locate ActsGeometry node" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

    return Fun4AllReturnCodes::EVENT_OK;
  }

  // ok return tuple<bool is_match, float stat,
  // float delta_phi (in pixels), float half_phi_width (in pixels),
  // float half_eta_width (in pixels)
  // float delta_z (in bins), flaot
  //

  std::pair<bool, float> TrkrClusterComparer::operator()(TrkrDefs::cluskey key_T, TrkrDefs::cluskey key_R)
  {
    // note: can use returned values, or just pull from these values
    layer = TrkrDefs::getLayer(key_T);
    if (layer > 55)
    {
      std::cout << " Error! Trying to compar cluster in layer > 55, "
                << "which is not programmed yet!" << std::endl;
      return {false, 0.};
    }

    in_mvtx = (layer < 3);
    in_intt = (layer > 2 && layer < 7);
    in_tpc = (layer > 6 && layer < 55);

    float phi_step = m_phistep[layer];
    float z_step = in_mvtx ? m_zstep_mvtx : m_zstep_tpc;

    clus_T = m_TruthClusters->findCluster(key_T);
    clus_R = m_RecoClusters->findCluster(key_R);

    phi_T = clus_T->getPosition(0);
    phi_R = clus_R->getPosition(0);
    phisize_R = clus_R->getPhiSize() * m_nphi_widths;  // * phi_step;
    phisize_T = clus_T->getPhiSize() * m_nphi_widths;  // * phi_step; // only for user to get, if they want

    z_T = clus_T->getPosition(1);
    z_R = clus_R->getPosition(1);

    if (!in_intt)
    {
      zsize_R = clus_R->getZSize() * m_nz_widths;  // * z_step;
      zsize_T = clus_R->getZSize() * m_nz_widths;  // * z_step;
    }

    phi_delta = fabs(phi_T - phi_R);
    while (phi_delta > M_PI)
    {
      phi_delta = fabs(phi_delta - 2 * M_PI);
    }
    phi_delta /= phi_step;

    z_delta = fabs(z_T - z_R) / z_step;

    /* float phi_stat = (m_nphi_widths * phisize_R ); */

    float phi_stat = std::max(phisize_T, phisize_R);
    float fit_statistic = (phi_delta / phi_stat);
    is_match = (fit_statistic <= 1.);

    float z_stat = 0;
    if (!in_intt)
    {
      z_stat = std::max(zsize_T, zsize_R);

      is_match = is_match && (z_delta < z_stat);
      fit_statistic += z_delta / z_stat;
    }
    return {is_match, fit_statistic};
  }

  ClusLoc TrkrClusterComparer::clusloc_PHG4(
      std::pair<TrkrDefs::hitsetkey, TrkrDefs::cluskey> input)
  {
    auto cluster = m_TruthClusters->findCluster(input.second);
    Eigen::Vector3d gloc = m_ActsGeometry->getGlobalPosition(input.second, cluster);
    return {TrkrDefs::getLayer(input.first), gloc,
            (int) cluster->getPhiSize(), (int) cluster->getZSize()};
  }

  ClusLoc TrkrClusterComparer::clusloc_SVTX(
      std::pair<TrkrDefs::hitsetkey, TrkrDefs::cluskey> input)
  {
    auto cluster = m_RecoClusters->findCluster(input.second);
    Eigen::Vector3d gloc = m_ActsGeometry->getGlobalPosition(input.second, cluster);
    return {TrkrDefs::getLayer(input.first), gloc,
            (int) cluster->getPhiSize(), (int) cluster->getZSize()};
  }

  // Implementation of the iterable struct to get cluster keys from
  // a SvtxTrack. It is used like:
  // for (auto& cluskey : ClusKeyIter(svtx_track)) {
  //    ... // do things with cluster keys
  // }
  ClusKeyIter::ClusKeyIter(SvtxTrack* _track)
    : track{_track}
    , in_silicon{_track->get_silicon_seed() != nullptr}
    , has_tpc{_track->get_tpc_seed() != nullptr}
    , no_data{!in_silicon && !has_tpc}
  {
  }

  ClusKeyIter ClusKeyIter::begin()
  {
    ClusKeyIter iter0{track};
    if (iter0.no_data)
    {
      return iter0;
    }
    if (iter0.in_silicon)
    {
      iter0.iter = track->get_silicon_seed()->begin_cluster_keys();
      iter0.iter_end_silicon = track->get_silicon_seed()->end_cluster_keys();
    }
    else if (has_tpc)
    {
      iter0.iter = track->get_tpc_seed()->begin_cluster_keys();
    }
    return iter0;
  }

  ClusKeyIter ClusKeyIter::end()
  {
    ClusKeyIter iter0{track};
    if (iter0.no_data)
    {
      return iter0;
    }
    if (has_tpc)
    {
      iter0.iter = track->get_tpc_seed()->end_cluster_keys();
    }
    else if (in_silicon)
    {
      iter0.iter = track->get_silicon_seed()->end_cluster_keys();
    }
    return iter0;
  }

  void ClusKeyIter::operator++()
  {
    if (no_data)
    {
      return;
    }
    ++iter;
    if (in_silicon && has_tpc && iter == iter_end_silicon)
    {
      in_silicon = false;
      iter = track->get_tpc_seed()->begin_cluster_keys();
    }
  }

  bool ClusKeyIter::operator!=(const ClusKeyIter& rhs)
  {
    if (no_data)
    {
      return false;
    }
    return iter != rhs.iter;
  }

  TrkrDefs::cluskey ClusKeyIter::operator*()
  {
    return *iter;
  }

  TrkrClusterContainer* ClusCntr::get_PHG4_clusters()
  {
    if (comp == nullptr)
    {
      return nullptr;
    }
    else
    {
      return comp->m_TruthClusters;
    }
  }

  TrkrClusterContainer* ClusCntr::get_SVTX_clusters()
  {
    if (comp == nullptr)
    {
      return nullptr;
    }
    else
    {
      return comp->m_RecoClusters;
    }
  }

  std::array<int, 5> ClusCntr::cntclus(Vector& keys)
  {
    std::array<int, 5> cnt{0, 0, 0, 0, 0};
    for (auto& it : keys)
    {
      cnt[trklayer_0123(it.first)] += 1;
    }
    for (int i = 0; i < 4; ++i)
    {
      cnt[4] += cnt[i];
    }
    return cnt;
  }

  std::array<int, 5> ClusCntr::cnt_matchedclus(Vector& keys, std::vector<bool>& matches)
  {
    std::array<int, 5> cnt{0, 0, 0, 0, 0};
    if (keys.size() != matches.size())
    {
      std::cout << PHWHERE << " matching and key vector not the same size. "
                << std::endl
                << " run find_matches() first." << std::endl;
      return cnt;
    }
    for (unsigned int i = 0; i < keys.size(); ++i)
    {
      if (matches[i])
      {
        cnt[trklayer_0123(keys[i].first)] += 1;
      }
    }
    for (int i = 0; i < 4; ++i)
    {
      cnt[4] += cnt[i];
    }
    return cnt;
  }

  int ClusCntr::addClusKeys(SvtxTrack* track)
  {
    svtx_keys.clear();
    for (auto ckey : ClusKeyIter(track))
    {
      svtx_keys.push_back({TrkrDefs::getHitSetKeyFromClusKey(ckey), ckey});
    }
    std::sort(svtx_keys.begin(), svtx_keys.end());
    return svtx_keys.size();
  }

  int ClusCntr::addClusKeys(TrkrTruthTrack* track)
  {
    phg4_keys.clear();
    for (auto ckey : track->getClusters())
    {
      phg4_keys.push_back({TrkrDefs::getHitSetKeyFromClusKey(ckey), ckey});
    }
    std::sort(phg4_keys.begin(), phg4_keys.end());
    return phg4_keys.size();
  }

  void ClusCntr::reset()
  {
    phg4_keys.clear();
    phg4_matches.clear();
    svtx_keys.clear();
    svtx_matches.clear();
  }

  std::array<int, 3> ClusCntr::find_matches()
  {
    if (comp == nullptr)
    {
      std::cout << PHWHERE
                << " Won't compare tracks because of missing TrkrClusterComparer" << std::endl;
      return {0, 0, 0};
    }
    // find the matches between the svtx_keys and phg4_keys
    // also keep track of the sum of the comparison between then

    // ---------------------------------
    // set aliases for notation cleaness
    // use A for PHG4 and B for SVTX
    auto& vA = phg4_keys;
    auto& vB = svtx_keys;

    auto& matchesA = phg4_matches;
    auto& matchesB = svtx_matches;

    match_stat = 0.;

    // matches will say, cluster by cluster, which clusters are matched
    matchesA = std::vector<bool>(vA.size(), false);
    matchesB = std::vector<bool>(vB.size(), false);

    // user iterators to access the vectors
    auto iA0 = vA.begin();
    auto iA1 = vA.end();

    auto iB0 = vB.begin();
    auto iB1 = vB.end();

    auto iA = iA0;
    auto iB = iB0;

    int n_match{0};

    while (iA != iA1 && iB != iB1)
    {
      if (iA->first == iB->first)
      {
        auto hitset = iA->first;

        // must compare ALL sets of iA and iB with this same hitset
        auto sAend = iA + 1;  // search A end
        while (sAend != iA1 && sAend->first == hitset)
        {
          ++sAend;
        }

        auto sBend = iB + 1;  // search B end
        while (sBend != iB1 && sBend->first == hitset)
        {
          ++sBend;
        }

        for (auto A = iA; A != sAend; ++A)
        {
          for (auto B = iB; B != sBend; ++B)
          {
            auto comp_val = comp->operator()(A->second, B->second);
            if (comp_val.first)
            {
              matchesA[A - iA0] = true;
              matchesB[B - iB0] = true;
              match_stat += comp_val.second;
              ++n_match;
            }
          }
        }
        iA = sAend;
        iB = sBend;
      }
      else if (iA->first < iB->first)
      {
        ++iA;
      }
      else
      {
        ++iB;
      }
    }
    return {n_match, (int) phg4_keys.size(), (int) svtx_keys.size()};
  }

  std::array<int, 3> ClusCntr::find_matches(TrkrTruthTrack* g4_track, SvtxTrack* sv_track)
  {
    addClusKeys(sv_track);
    addClusKeys(g4_track);
    return find_matches();
  }

  int ClusCntr::phg4_n_matched()
  {
    return std::accumulate(phg4_matches.begin(), phg4_matches.end(), 0);
  }

  int ClusCntr::svtx_n_matched()
  {
    return std::accumulate(svtx_matches.begin(), svtx_matches.end(), 0);
  }

  std::vector<ClusLoc> ClusCntr::phg4_clusloc_all()
  {
    std::vector<ClusLoc> vec{};
    for (auto& cluspair : phg4_keys)
    {
      vec.push_back(comp->clusloc_PHG4(cluspair));
    }
    return vec;
  }

  std::vector<ClusLoc> ClusCntr::phg4_clusloc_unmatched()
  {
    std::vector<ClusLoc> vec{};
    auto cnt = phg4_keys.size();
    for (unsigned int i = 0; i < cnt; ++i)
    {
      if (!phg4_matches[i])
      {
        vec.push_back(comp->clusloc_PHG4(phg4_keys[i]));
      }
    }
    return vec;
  }

  std::vector<ClusLoc> ClusCntr::svtx_clusloc_all()
  {
    std::vector<ClusLoc> vec{};
    for (auto& cluspair : svtx_keys)
    {
      vec.push_back(comp->clusloc_SVTX(cluspair));
    }
    return vec;
  }

  std::vector<ClusLoc> ClusCntr::svtx_clusloc_unmatched()
  {
    std::vector<ClusLoc> vec{};
    auto cnt = svtx_keys.size();
    for (unsigned int i = 0; i < cnt; ++i)
    {
      if (!svtx_matches[i])
      {
        vec.push_back(comp->clusloc_SVTX(svtx_keys[i]));
      }
    }
    return vec;
  }

  std::vector<ClusLoc> ClusCntr::clusloc_matched()
  {
    std::vector<ClusLoc> vec{};
    auto cnt = phg4_keys.size();
    for (unsigned int i = 0; i < cnt; ++i)
    {
      if (phg4_matches[i])
      {
        vec.push_back(comp->clusloc_PHG4(phg4_keys[i]));
      }
    }
    return vec;
  }

  /* ClusCntr::layer_xyzLoc ClusCntr::xyzLoc(std::pair<TrkrDefs::hitsetkey,TrkrDefs::cluskey) { */
  /*   if (geom == nullptr) { */
  /*     std::cout << PHWHERE << " fatal: geom, type ActsGeometry*, must be set to call xyzLoc!" << std::endl; */
  /*     return {-1,{FLT_MAX,FLT_MAX,FLT_MAX}}; */
  /*   } */
  /*   Eigen::Vector3d gloc =    m_ActsGeometry->getGlobalPosition(reco_ckey, cluster); */
  /* } */

}  // namespace G4Eval
