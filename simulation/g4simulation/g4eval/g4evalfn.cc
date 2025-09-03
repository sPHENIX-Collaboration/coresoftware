#include <g4tracking/EmbRecoMatchContainer.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include "TrkrClusLoc.h"
#include "TrkrClusterIsMatcher.h"
#include "g4evalfn.h"

#include <algorithm>
#include <cmath>
#include <set>

namespace g4evalfn
{

  int trklayer_det(int layer)
  {
    if (layer < 3)
    {
      return DET::MVTX;
    }
    if (layer < 7)
    {
      return DET::INTT;
    }
    if (layer < 55)
    {
      return DET::TPC;
    }
    return DET::TPOT;
  }

  int trklayer_det(TrkrDefs::hitsetkey key)
  {
    return trklayer_det(TrkrDefs::getLayer(key));
  }

  int trklayer_det(TrkrDefs::cluskey key)
  {
    return trklayer_det(TrkrDefs::getLayer(key));
  }

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

  TrkrClusLoc clusloc_PHG4(TrkrClusterIsMatcher* ismatcher, TrkrDefs::cluskey key)
  {
    auto cluster = ismatcher->m_TruthClusters->findCluster(key);
    Eigen::Vector3d gloc = ismatcher->m_ActsGeometry->getGlobalPosition(key, cluster);
    return {TrkrDefs::getLayer(key), gloc,
            cluster->getPosition(0), cluster->getPhiSize(), cluster->getPosition(1), cluster->getZSize(), key};
  }

  TrkrClusLoc clusloc_SVTX(TrkrClusterIsMatcher* ismatcher, TrkrDefs::cluskey key)
  {
    auto cluster = ismatcher->m_RecoClusters->findCluster(key);
    Eigen::Vector3d gloc = ismatcher->m_ActsGeometry->getGlobalPosition(key, cluster);
    return {TrkrDefs::getLayer(key), gloc,
            cluster->getPosition(0), cluster->getPhiSize(), cluster->getPosition(1), cluster->getZSize(), key};
  }

  float calc_match_statistic(TrkrClusterIsMatcher* ismatcher, TrkrDefs::cluskey key_A, TrkrDefs::cluskey key_B)
  {
    int layer = TrkrDefs::getLayer(key_A);
    auto det_0123 = trklayer_det(layer);
    auto clus_T = ismatcher->m_TruthClusters->findCluster(key_A);
    auto clus_R = ismatcher->m_RecoClusters->findCluster(key_B);
    auto dphi = abs_dphi(clus_T->getPosition(0), clus_R->getPosition(0));

    float stat = 0.;
    switch (det_0123)
    {
    case 0:
    {  // MVTX
      if (ismatcher->single_pixel_phi_MVTX)
      {
        stat += dphi / ismatcher->tol_pitch_phi[layer];
      }
      else
      {
        stat += dphi / (ismatcher->tol_pitch_phi[layer] * std::max(clus_T->getPhiSize(), clus_R->getPhiSize()));
      }
      const float delta_z = std::abs(clus_T->getPosition(1) - clus_R->getPosition(1));
      if (ismatcher->single_pixel_z_MVTX)
      {
        stat += delta_z / ismatcher->tol_pitch_z_MVTX;
      }
      else
      {
        stat += delta_z / (ismatcher->tol_pitch_z_MVTX * std::max(clus_T->getZSize(), clus_R->getZSize()));
      }
    }
    break;
    case 1:
    {
      if (ismatcher->single_pixel_phi_INTT)
      {
        stat += dphi / ismatcher->tol_pitch_phi[layer];
      }
      else
      {
        stat += dphi / (ismatcher->tol_pitch_phi[layer] * std::max(clus_T->getPhiSize(), clus_R->getPhiSize()));
      }
    }
    break;
    case 2:
    {  // TPC
      if (ismatcher->single_bin_phi_TPC)
      {
        stat += dphi / ismatcher->tol_pitch_phi[layer];
      }
      else
      {
        stat += dphi / (ismatcher->tol_pitch_phi[layer] * std::max(clus_T->getPhiSize(), clus_R->getPhiSize()));
      }
      const float delta_t = std::abs(clus_T->getPosition(1) - clus_R->getPosition(1));
      if (ismatcher->single_bin_t_TPC)
      {
        stat += delta_t / ismatcher->tol_step_t_TPC;
      }
      else
      {
        stat += delta_t / (ismatcher->tol_step_t_TPC * std::max(clus_T->getZSize(), clus_R->getZSize()));
      }
    }
    break;
    // case 3:  // TPOT // empty case triggers clang-tidy warning
    //   break;
    default:
      break;
    }
    return stat;
  }

}  // namespace g4evalfn
