/*!
 *  \file PHCASiliconSeeding.cc
 *  \brief Silicon track seeding using ALICE-style "cellular automaton" (CA) algorithm
 *  \detail
 *  \author Michael Peters
 */

#include "PHCASiliconSeeding.h"

// sPHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

// trackbase_historic includes
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>  // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrClusterIterationMapv1.h>
#include <trackbase/TrkrDefs.h>  // for getLayer, clu...
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeed_v2.h>
#include <trackbase_historic/TrackSeedHelper.h>

// BOOST for combi seeding
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/policies/compare.hpp>

#include <algorithm>
#include <iostream>
#include <memory>
#include <numeric>
#include <unordered_set>
#include <utility>  // for pair, make_pair
#include <vector>

// anonymous namespace for local functions
namespace
{
  // square
  template <class T>
  inline constexpr T square(const T& x)
  {
    return x * x;
  }

  /// phi angle of Acts::Vector3
  inline double get_phi(const Acts::Vector3& position)
  {
    double phi = std::atan2(position.y(), position.x());
    if (phi < 0)
    {
      phi += 2. * M_PI;
    }
    return phi;
  }

  // note: assumes that a and b are in same range of phi;
  // this will fail if a\in[-2 pi,0] and b\in[0,2 pi] 
  // in this case is ok, as all are atan2 which [-pi,pi]
  // inline float wrap_dphi(float a, float b) {
  //   float   _dphi = b-a;
  //   return (_dphi < -M_PI) ? _dphi += 2*M_PI
  //        : (_dphi >  M_PI) ? _dphi -= 2*M_PI
  //        : _dphi;
  // }

  /// pseudo rapidity of Acts::Vector3
  /* inline double get_eta(const Acts::Vector3& position) */
  /* { */
  /*   const double norm = std::sqrt(square(position.x()) + square(position.y()) + square(position.z())); */
  /*   return std::log((norm + position.z()) / (norm - position.z())) / 2; */
  /* } */

}  // namespace

// using namespace ROOT::Minuit2;
namespace bgi = boost::geometry::index;

PHCASiliconSeeding::PHCASiliconSeeding(
    const std::string& name,
    unsigned int start_layer,
    unsigned int end_layer,
    unsigned int min_clusters_per_track,
    float neighbor_phi_width,
    float neighbor_z_width)
  : PHTrackSeeding(name)
  , _start_layer(start_layer)
  , _end_layer(end_layer)
  , _min_clusters_per_track(min_clusters_per_track)
  , _neighbor_phi_width(neighbor_phi_width)
  , _neighbor_z_width(neighbor_z_width)
{
}

int PHCASiliconSeeding::InitializeGeometry(PHCompositeNode* topNode)
{
  // geometry
  m_tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_tGeometry)
  {
    std::cout << PHWHERE << "No acts tracking geometry, can't proceed" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  // cluster container
  m_clusterMap = findNode::getClass<TrkrClusterContainer>(topNode,"TRKR_CLUSTER");
  if (!m_clusterMap)
  {
    std::cout << PHWHERE << "No cluster map, can't proceed" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  // cluster crossing associations
  m_clusterCrossingMap = findNode::getClass<TrkrClusterCrossingAssoc>(topNode,"TRKR_CLUSTERCROSSINGASSOC");
  if (!m_clusterCrossingMap)
  {
    std::cout << PHWHERE << "No cluster crossing association map, can't proceed" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

Acts::Vector3 PHCASiliconSeeding::getGlobalPosition(TrkrDefs::cluskey key, TrkrCluster* cluster) const
{
  return m_tGeometry->getGlobalPosition(key, cluster);
}

void PHCASiliconSeeding::QueryTree(const bgi::rtree<PHCASiliconSeeding::pointKey, bgi::quadratic<16>>& rtree, double phimin, double z_min, double phimax, double z_max, std::vector<pointKey>& returned_values) const
{
  bool query_both_ends = false;
  if (phimin < 0)
  {
    query_both_ends = true;
    phimin += 2 * M_PI;
  }
  if (phimax > 2 * M_PI)
  {
    query_both_ends = true;
    phimax -= 2 * M_PI;
  }
  if (query_both_ends)
  {
    rtree.query(bgi::intersects(box(point(phimin, z_min), point(2 * M_PI, z_max))), std::back_inserter(returned_values));
    rtree.query(bgi::intersects(box(point(0., z_min), point(phimax, z_max))), std::back_inserter(returned_values));
  }
  else
  {
    rtree.query(bgi::intersects(box(point(phimin, z_min), point(phimax, z_max))), std::back_inserter(returned_values));
  }
}

std::pair<PHCASiliconSeeding::PositionMap, PHCASiliconSeeding::keyListPerLayer> PHCASiliconSeeding::FillGlobalPositions()
{
  keyListPerLayer ckeys;
  PositionMap cachedPositions;
  cachedPositions.reserve(_cluster_map->size());  // avoid resizing mid-execution

  std::vector<TrkrDefs::hitsetkey> hskeys = _cluster_map->getHitSetKeys(TrkrDefs::mvtxId);
  std::vector<TrkrDefs::hitsetkey> intt_hskeys = _cluster_map->getHitSetKeys(TrkrDefs::inttId);
  hskeys.insert(hskeys.end(),intt_hskeys.begin(),intt_hskeys.end());

  for (const auto& hitsetkey : hskeys)
  {
    // filter for strobes with INTT clusters likely to be in them
    if(TrkrDefs::getTrkrId(hitsetkey) == TrkrDefs::mvtxId)
    {
      int strobeId = MvtxDefs::getStrobeId(hitsetkey);
      if(strobeId < _lowest_allowed_strobeid || strobeId > _highest_allowed_strobeid)
      {
        continue;
      }
    }
    if(TrkrDefs::getLayer(hitsetkey) < _start_layer || TrkrDefs::getLayer(hitsetkey) > _end_layer)
    {
      if (Verbosity() > 2)
      {
        std::cout << "skipping layer: " << TrkrDefs::getLayer(hitsetkey) << std::endl;
      }
      continue;
    }
    auto range = _cluster_map->getClusters(hitsetkey);
    for (auto clusIter = range.first; clusIter != range.second; ++clusIter)
    {
      TrkrDefs::cluskey ckey = clusIter->first;
      TrkrCluster* cluster = clusIter->second;
      unsigned int layer = TrkrDefs::getLayer(ckey);
      if (_iteration_map != nullptr && _n_iteration > 0)
      {
        if (_iteration_map->getIteration(ckey) > 0)
        {
          continue;  // skip hits used in a previous iteration
        }
      }

      // get global position, convert to Acts::Vector3 and store in map
      const Acts::Vector3 globalpos_d = getGlobalPosition(ckey, cluster);
      const Acts::Vector3 globalpos = {globalpos_d.x(), globalpos_d.y(), globalpos_d.z()};
      cachedPositions.insert(std::make_pair(ckey, globalpos));
      ckeys[layer].push_back(ckey);
    }
  }
  return std::make_pair(cachedPositions, ckeys);
}

std::vector<PHCASiliconSeeding::coordKey> PHCASiliconSeeding::FillTree(bgi::rtree<PHCASiliconSeeding::pointKey, bgi::quadratic<16>>& _rtree, const PHCASiliconSeeding::keyList& ckeys, const PHCASiliconSeeding::PositionMap& globalPositions, const int layer)
{
  // Fill _rtree with the clusters in ckeys; remove duplicates, and return a vector of the coordKeys
  // Note that layer is only used for a cout statement
  int n_dupli = 0;
  std::vector<coordKey> coords;
  // _rtree.clear(); DO NOT clear rtree if we're filling it multiple times (which for silicon, we are)
  /* _rtree.reserve(ckeys.size()); */
  for (const auto& ckey : ckeys)
  {
    const auto& globalpos_d = globalPositions.at(ckey);
    const double clus_phi = get_phi(globalpos_d);
    const double clus_z = globalpos_d.z();
    if (Verbosity() > 5)
    {
      /* int layer = TrkrDefs::getLayer(ckey); */
      std::cout << "Found cluster " << ckey << " in layer " << layer << std::endl;
    }
    //std::vector<pointKey> testduplicate;
    //QueryTree(_rtree, clus_phi - 0.00001, clus_z - 0.00001, clus_phi + 0.00001, clus_z + 0.00001, testduplicate);
    //if (!testduplicate.empty())
    //{
    //  ++n_dupli;
    //  continue;
    //}
    coords.push_back({{static_cast<float>(clus_phi), static_cast<float>(clus_z)}, ckey});
    _rtree.insert(std::make_pair(point(clus_phi, globalpos_d.z()), ckey));
  }
  if (Verbosity() > 5)
  {
    std::cout << "nhits in layer(" << layer << "): " << coords.size() << std::endl;
  }
  if (Verbosity() > 3)
  {
    std::cout << "number of duplicates : " << n_dupli << std::endl;
  }
  return coords;
}

int PHCASiliconSeeding::Process(PHCompositeNode* /*topNode*/)
{
  if (Verbosity() > 3)
  {
    std::cout << " Process...  " << std::endl;
  }
  if (_n_iteration > 0)
  {
    if (!_iteration_map)
    {
      std::cerr << PHWHERE << "Cluster Iteration Map missing, aborting." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  PositionMap globalPositions;
  keyListPerLayer ckeys;
  std::tie(globalPositions, ckeys) = FillGlobalPositions();

  //  int numberofseeds = 0;
  // numberofseeds += FindSeeds(globalPositions, ckeys);

  for(auto& rtree : _rtrees)
  {
    rtree.clear();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHCASiliconSeeding::FindSeeds(const PHCASiliconSeeding::PositionMap& globalPositions, const PHCASiliconSeeding::keyListPerLayer& ckeys)
{
  std::vector<std::vector<Triplet>> triplets = CreateLinks(globalPositions, ckeys);
  keyLists trackSeedKeyLists = FollowLinks(triplets);
  
  std::vector<TrackSeed_v2> seeds = FitSeeds(trackSeedKeyLists, globalPositions);
  HelixPropagate(seeds, globalPositions);
  HelixPropagate(seeds, globalPositions); // each call extends seed by up to one cluster
  
  publishSeeds(seeds);
  return seeds.size();
}

bool PHCASiliconSeeding::ClusterTimesAreCompatible(const TrkrDefs::cluskey clus_a, const TrkrDefs::cluskey clus_b) const
{
  const int time_index = GetClusterTimeIndex(clus_a);
  return ClusterTimesAreCompatible(TrkrDefs::getTrkrId(clus_a),time_index,clus_b);
}

bool PHCASiliconSeeding::ClusterTimesAreCompatible(const uint8_t trkr_id, const int time_index, const TrkrDefs::cluskey ckey) const
{
  if(TrkrDefs::getTrkrId(ckey) != trkr_id)
  {
    return true; // clusters can only be compared within the same detector, so compatibility cut is not applicable
  }
  else if(trkr_id == TrkrDefs::mvtxId)
  {
    if(Verbosity()>3)
    {
      std::cout << "strobe a " << time_index << " strobe b " << MvtxDefs::getStrobeId(ckey) << std::endl;
    }
    return (time_index == MvtxDefs::getStrobeId(ckey)); // cut on same MVTX strobe
  }
  else if(trkr_id == TrkrDefs::inttId)
  {
    short crossing = GetCleanINTTClusterCrossing(ckey);
    if(Verbosity()>3)
    {
      std::cout << "crossing a " << time_index << " crossing b " << crossing << std::endl;
    }
    return (crossing != SHRT_MAX && time_index != SHRT_MAX) && (abs(crossing - time_index) <= 1);
  }
  else
  {
    return true; // other detectors don't carry interpretable time info
  }
}

int PHCASiliconSeeding::GetClusterTimeIndex(const TrkrDefs::cluskey ckey) const
{
  if(TrkrDefs::getTrkrId(ckey) == TrkrDefs::mvtxId)
  {
    return MvtxDefs::getStrobeId(ckey);
  }
  else if(TrkrDefs::getTrkrId(ckey) == TrkrDefs::inttId)
  {
    return GetCleanINTTClusterCrossing(ckey);
  }
  else
  {
    return SHRT_MAX; // other detectors don't carry interpretable time info
  }
}

short PHCASiliconSeeding::GetCleanINTTClusterCrossing(const TrkrDefs::cluskey ckey) const
{
  std::set<short> crossings = GetINTTClusterCrossings(ckey);
  if(crossings.size()==0 || crossings.size() > 2)
  {
    if(Verbosity()>3 && crossings.size()>2)
    {
      std::cout << "more than two INTT crossings within cluster: ";
      for(short cross : crossings) std::cout << cross << " ";
      std::cout << std::endl;
    }
    return SHRT_MAX;
  }
  else if(crossings.size()==1)
  {
    return *(crossings.begin());
  }
  else // crossings.size()==2
  {
    // allow for crossing to be off by 1,
    // take lower crossing to be true value
    std::vector<short> crossings_vec;
    std::copy(crossings.begin(),crossings.end(),std::back_inserter(crossings_vec));

    if(abs(crossings_vec[1] - crossings_vec[0]) == 1)
    {
      std::cout << "INTT: resolving off-by-one between " << crossings_vec[0] << " and " << crossings_vec[1] << std::endl;
      return std::min(crossings_vec[0],crossings_vec[1]);
    }
    else
    {
      return SHRT_MAX;
    }
  }
}

std::set<short> PHCASiliconSeeding::GetINTTClusterCrossings(const TrkrDefs::cluskey ckey) const
{
  if(TrkrDefs::getTrkrId(ckey) != TrkrDefs::inttId)
  {
    return std::set<short>{}; // only INTT clusters have crossing info
  }
  else
  {
    std::set<short> crossings;
    TrkrCluster* clus = m_clusterMap->findCluster(ckey);
    if(!clus)
    {
      return std::set<short>{}; // a cluster that doesn't exist has no cluster crossing info
    }
    auto crossingrange = m_clusterCrossingMap->getCrossings(ckey);
    for(auto iter = crossingrange.first; iter != crossingrange.second; ++iter)
    {
      crossings.insert(iter->second);
    }
    return crossings;
  }
}

std::vector<std::vector<PHCASiliconSeeding::Triplet>> PHCASiliconSeeding::CreateLinks(const PHCASiliconSeeding::PositionMap& globalPositions, const PHCASiliconSeeding::keyListPerLayer& ckeys)
{
  std::vector<std::vector<Triplet>> triplets; // stored by layer

  std::vector<std::vector<coordKey>> clusterdata;

  for(size_t l = _start_layer; l<=_end_layer; l++)
  {
    const size_t l_index = l - _start_layer;
    _rtrees.push_back(bgi::rtree<pointKey, bgi::quadratic<16>>());
    triplets.push_back(std::vector<Triplet>());
    if(l==4)
    {
      const std::vector<coordKey> l4clusters = FillTree(_rtrees[3-_start_layer],ckeys[4],globalPositions,3);
      clusterdata[3-_start_layer].insert(clusterdata[3-_start_layer].end(),l4clusters.begin(),l4clusters.end());
    }
    else if(l==5)
    {
      clusterdata.push_back(FillTree(_rtrees[4-_start_layer],ckeys[5],globalPositions,4));
    }
    else if(l==6)
    {
      const std::vector<coordKey> l6clusters = FillTree(_rtrees[4-_start_layer],ckeys[6],globalPositions,4);
      clusterdata[4-_start_layer].insert(clusterdata[4-_start_layer].end(),l6clusters.begin(),l6clusters.end());
    }
    else
    {
      clusterdata.push_back(FillTree(_rtrees[l_index],ckeys[l],globalPositions,l));
    }
  }

  // form triplets in interior layers

  for(size_t l = _start_layer+1; l<=_end_layer-1; l++)
  {
    if(Verbosity()>2)
    {
      std::cout << "layer " << l << std::endl;
    }
    if(l>=4 && l<=6)
    {
      continue; // skip INTT double layers that have been combined into layers 3 and 4
    }
    const size_t l_index = l - _start_layer;

    for(coordKey centerCluster : clusterdata[l_index])
    {
      std::vector<pointKey> clustersBelow;
      std::vector<pointKey> clustersAbove;

      const float centerPhi = centerCluster.first[0];
      const float centerZ = centerCluster.first[1];

      const float dphiwindow_below = std::max(dphi_per_layer[l],dphi_per_layer[l-1]);
      const float dphiwindow_above = std::max(dphi_per_layer[l],dphi_per_layer[l+1]);
      const float dZwindow_below = std::max(dZ_per_layer[l],dZ_per_layer[l-1]);
      const float dZwindow_above = std::max(dZ_per_layer[l],dZ_per_layer[l+1]);

      QueryTree(_rtrees[l_index-1],
                centerPhi - dphiwindow_below,
                centerZ - dZwindow_below,
                centerPhi + dphiwindow_below,
                centerZ + dZwindow_below,
                clustersBelow);

      if(l_index>1)
      {
        const float dphiwindow_2below = std::max(dphi_per_layer[l],dphi_per_layer[l-2]);
        const float dZwindow_2below = std::max(dZ_per_layer[l],dZ_per_layer[l-2]);

        QueryTree(_rtrees[l_index-2],
                  centerPhi - dphiwindow_below - dphiwindow_2below,
                  centerZ - dZwindow_below - dZwindow_2below,
                  centerPhi + dphiwindow_below + dphiwindow_2below,
                  centerZ + dZwindow_below + dZwindow_2below,
                  clustersBelow);
      }

      QueryTree(_rtrees[l_index+1],
                centerPhi - dphiwindow_above,
                centerZ - dZwindow_above,
                centerPhi + dphiwindow_above,
                centerZ + dZwindow_above,
                clustersAbove);

      if(l_index<3)
      {
        const float dphiwindow_2above = std::max(dphi_per_layer[l],dphi_per_layer[l+2]);
        const float dZwindow_2above = std::max(dZ_per_layer[l],dZ_per_layer[l+2]);

        QueryTree(_rtrees[l_index+2],
                  centerPhi - dphiwindow_above - dphiwindow_2above,
                  centerZ - dZwindow_above - dZwindow_2above,
                  centerPhi + dphiwindow_above + dphiwindow_2above,
                  centerZ + dZwindow_above + dZwindow_2above,
                  clustersAbove);
      }

      if(Verbosity()>3)
      {
        std::cout << std::endl;
        std::cout << "found " << clustersBelow.size() << " clusters below, " << clustersAbove.size() << " clusters above" << std::endl;
      }

      float best_cos_angle = 1e9;
      TrkrDefs::cluskey best_below_ckey = 0;
      TrkrDefs::cluskey best_above_ckey = 0;

      std::vector<TrkrDefs::cluskey> passing_below_ckeys;
      std::vector<TrkrDefs::cluskey> passing_above_ckeys;

      const TrkrDefs::cluskey center_ckey = centerCluster.second;
      const Acts::Vector3 gpos_center = globalPositions.at(center_ckey);
      const int time_center = GetClusterTimeIndex(center_ckey);
      const uint8_t trkrid_center = TrkrDefs::getTrkrId(center_ckey);

      for(const pointKey& cbelow : clustersBelow)
      {
        if(!ClusterTimesAreCompatible(trkrid_center,time_center,cbelow.second))
        {
          if(Verbosity()>3)
          {
            std::cout << "below candidate has incompatible time" << std::endl;
          }
          continue;
        }
        const Acts::Vector3 gpos_below = globalPositions.at(cbelow.second);
        const Acts::Vector3 delta_below = gpos_below - gpos_center;

        for(const pointKey& cabove : clustersAbove)
        {
          if(!ClusterTimesAreCompatible(trkrid_center,time_center,cabove.second))
          {
            if(Verbosity()>3)
            {
              std::cout << "above candidate has incompatible time" << std::endl;
            }
            continue;
          }
          const Acts::Vector3 gpos_above = globalPositions.at(cabove.second);
          const Acts::Vector3 delta_above = gpos_above - gpos_center;

          float cos_angle;
          float dot_product;
          float mag2_below;
          float mag2_above;

          if(Verbosity()>3)
          {
            std::cout << "candidate triplet: " << std::endl;
            std::cout << "layer " << (int)TrkrDefs::getLayer(cbelow.second) << ": " << gpos_below.x() << ", " << gpos_below.y() << ", " << gpos_below.z() << std::endl;
            std::cout << "layer " << (int)TrkrDefs::getLayer(centerCluster.second) << ": " << gpos_center.x() << ", " << gpos_center.y() << ", " << gpos_center.z() << std::endl;
            std::cout << "layer " << (int)TrkrDefs::getLayer(cabove.second) << ": " << gpos_above.x() << ", " << gpos_above.y() << ", " << gpos_above.z() << std::endl;
          }

          if(l>=2 && l<=6) // use xy breaking angle only for any triplets that include INTT clusters
          {
            mag2_below = delta_below.x()*delta_below.x() + delta_below.y()*delta_below.y();
            mag2_above = delta_above.x()*delta_above.x() + delta_above.y()*delta_above.y();
            dot_product = delta_below.x()*delta_above.x() + delta_below.y()*delta_above.y();
          }
          else
          {
            mag2_below = delta_below.x()*delta_below.x() + delta_below.y()*delta_below.y() + delta_below.z()*delta_below.z();
	    mag2_above = delta_above.x()*delta_above.x() + delta_above.y()*delta_above.y() + delta_above.z()*delta_above.z();
	    dot_product = delta_below.x()*delta_above.x() + delta_below.y()*delta_above.y() + delta_below.z()*delta_above.z();
          }

          cos_angle = dot_product/sqrt(mag2_below*mag2_above);

          if(Verbosity()>3)
          {
            std::cout << "delta_below: " << std::endl;
            std::cout << delta_below.x() << ", " << delta_below.y() << ", " << delta_below.z() << " (magnitude " << sqrt(mag2_below) << ")" << std::endl;
            std::cout << "delta_above: " << std::endl;
            std::cout << delta_above.x() << ", " << delta_above.y() << ", " << delta_above.z() << " (magnitude " << sqrt(mag2_above) << ")" << std::endl;
            std::cout << "dot product: " << dot_product << std::endl;
            std::cout << "cos(breaking angle): " << cos_angle << std::endl;
          }

          if(cos_angle < _max_cos_angle)
          {
            if(_use_best)
            {
              if(cos_angle < best_cos_angle)
              {
                if(Verbosity()>3)
                {
                  std::cout << "beats best cos(angle) of " << best_cos_angle << std::endl;
                }
                best_cos_angle = cos_angle;
	        best_below_ckey = cbelow.second;
	        best_above_ckey = cabove.second;
              }
            }
	    else
	    {
              if(Verbosity()>3)
              {
                std::cout << "passes straightness criterion" << std::endl;
              }
              passing_below_ckeys.push_back(cbelow.second);
	      passing_above_ckeys.push_back(cabove.second);
	    }
          }
        }
      }

      if(_use_best && best_cos_angle < 1.)
      {
        if(Verbosity()>3)
        {
          std::cout << "adding triplet" << std::endl;
        }
        triplets[l_index].push_back({best_below_ckey,centerCluster.second,best_above_ckey});
      }
      else
      {
        for(size_t i=0;i<passing_below_ckeys.size(); i++)
        {
          triplets[l_index].push_back({passing_below_ckeys[i],centerCluster.second,passing_above_ckeys[i]});
        }
      }
    }
    if(Verbosity() > 1)
    {
      std::cout << "layer: " << l << " formed " << triplets[l_index].size() << " triplets" << std::endl;
    }
  }
  return triplets;
}

std::vector<PHCASiliconSeeding::keyList> PHCASiliconSeeding::FollowLinks(const std::vector<std::vector<PHCASiliconSeeding::Triplet>>& triplets)
{
  std::vector<keyList> finishedSeeds;
  std::vector<keyList> growingSeeds;

  for(const Triplet& start_triplet : triplets[1])
  {
    growingSeeds.push_back({start_triplet.bottom,start_triplet.center,start_triplet.top});
  }
  if(Verbosity() > 1)
  {
    std::cout << "Started with " << growingSeeds.size() << " stubs" << std::endl;
  }
  
  for(size_t l=_start_layer+1; l<=_end_layer-1; l++)
  {
    if(Verbosity() > 1)
    {
      std::cout << "layer " << l << std::endl;
      std::cout << growingSeeds.size() << " still-growing seeds" << std::endl;
      std::cout << finishedSeeds.size() << " finished seeds" << std::endl;
    }
    std::vector<keyList> tempSeeds;
    const size_t l_index = l - _start_layer;
    // grow existing seeds
    if(Verbosity()>3)
    {
      std::cout << "growing current seeds" << std::endl;
    }
    for(const keyList& seed : growingSeeds)
    {
      if(Verbosity()>3)
      {
        std::cout << "current keys: ";
        for(const TrkrDefs::cluskey& key : seed) std::cout << (uint64_t)key << ", ";
        std::cout << std::endl;
      }
      const TrkrDefs::cluskey currentTop = seed.back();
      const TrkrDefs::cluskey currentCenter = seed.crbegin()[1];
      bool finished = true;
      for(const Triplet& candidate_triplet : triplets[l_index+1])
      {
        if(candidate_triplet.center == currentTop && candidate_triplet.bottom == currentCenter)
        {
          if(Verbosity()>3)
          {
            std::cout << "found next candidate -- keys are " << (uint64_t)candidate_triplet.bottom << ", " << (uint64_t)candidate_triplet.center << ", " << (uint64_t)candidate_triplet.top << std::endl;
          }
          finished = false;
          keyList tempSeed = seed;
          tempSeed.push_back(candidate_triplet.top);
          tempSeeds.push_back(tempSeed);
        }
      }
      if(finished) finishedSeeds.push_back(seed);
    }
    // find starts of new seeds
    int new_seed_count = 0;
    if(Verbosity()>3)
    {
      std::cout << "starting new seeds" << std::endl;
    }
    for(const Triplet& triplet : triplets[l_index+1])
    {
      if(Verbosity()>3)
      {
        std::cout << "candidate triplet: " << (uint64_t)triplet.bottom << ", " << (uint64_t)triplet.center << ", " << (uint64_t)triplet.top << std::endl;
      }
      bool has_existing_seed = false;
      for(const keyList& seed : tempSeeds)
      {
        if(seed.back()==triplet.top && seed.crbegin()[1]==triplet.center && seed.crbegin()[2]==triplet.bottom)
        {
          if(Verbosity()>3)
          {
            std::cout << "has existing seed with keys ";
            for(const TrkrDefs::cluskey& key : seed) std::cout << (uint64_t)key << ", ";
            std::cout << std::endl;
          }
          has_existing_seed = true;
        }
      }
      if(!has_existing_seed)
      {
        if(Verbosity()>3)
        {
          std::cout << "did not find existing seed" << std::endl;
        }
        new_seed_count++;
        tempSeeds.push_back({triplet.bottom, triplet.center, triplet.top});
      }
    }
    if(Verbosity() > 1)
    {
      std::cout << "started " << new_seed_count << " new seeds this layer" << std::endl;
    }
    growingSeeds = tempSeeds;
  }

  finishedSeeds.insert(finishedSeeds.end(),growingSeeds.begin(),growingSeeds.end());

  return finishedSeeds;
}

float PHCASiliconSeeding::getSeedQuality(const TrackSeed_v2& seed, const PHCASiliconSeeding::PositionMap& globalPositions) const
{
  std::vector<std::pair<double,double>> xy_pts;
  std::vector<std::pair<double,double>> rz_pts;
  std::vector<float> xyerr;
  std::vector<float> zerr;
  for(auto iter = seed.begin_cluster_keys(); iter != seed.end_cluster_keys(); ++iter)
  {
    Acts::Vector3 pos = globalPositions.at(*iter);
    TrkrCluster* c = m_clusterMap->findCluster(*iter);
    xy_pts.push_back(std::make_pair(pos.x(),pos.y()));
    rz_pts.push_back(std::make_pair(sqrt(pos.x()*pos.x()+pos.y()*pos.y()),pos.z()));
    xyerr.push_back(c->getRPhiError());
    zerr.push_back(c->getZError());
  }
  std::vector<double> circle_residuals = TrackFitUtils::getCircleClusterResiduals(xy_pts,fabs(1./seed.get_qOverR()),seed.get_X0(),seed.get_Y0());
  std::vector<double> line_residuals = TrackFitUtils::getLineClusterResiduals(rz_pts,seed.get_slope(),seed.get_Z0());

  float chi2 = 0;
  for(size_t i=0; i<circle_residuals.size(); i++)
  {
    const float total_resid2 = circle_residuals[i]*circle_residuals[i]+line_residuals[i]*line_residuals[i];
    const float total_err2 = xyerr[i]*xyerr[i]+zerr[i]*zerr[i];
    chi2 += total_resid2/total_err2;
  }
  const int ndf = 2*seed.size_cluster_keys()-5;
  return chi2/ndf;
}

void PHCASiliconSeeding::HelixPropagate(std::vector<TrackSeed_v2>& seeds, const PHCASiliconSeeding::PositionMap& globalPositions) const
{
  for(TrackSeed_v2& seed : seeds)
  {
    if(Verbosity()>3)
    {
      std::cout << std::endl << std::endl;
      std::cout << "================================" << std::endl;
      seed.identify();
    }

    std::set<size_t> layers;
    for(size_t layer = _start_layer; layer <= _end_layer; layer++)
    {
      if(layer==5 || layer==6)
      {
        continue;
      }
      layers.insert(layer);
    }

    std::vector<Acts::Vector3> clusterpos;
    std::vector<TrkrDefs::cluskey> clusters;
    std::set<int> seed_strobes;

    for(auto iter = seed.begin_cluster_keys(); iter != seed.end_cluster_keys(); ++iter)
    {
      const TrkrDefs::cluskey ckey = *iter;
      clusters.push_back(ckey);
      clusterpos.push_back(globalPositions.at(ckey));
      if(TrkrDefs::getTrkrId(ckey) == TrkrDefs::mvtxId)
      {
        seed_strobes.insert(MvtxDefs::getStrobeId(ckey));
      }

      // remove any layers in which we already have clusters
      // handle INTT double-layer situation
      if(TrkrDefs::getLayer(ckey)==3 || TrkrDefs::getLayer(ckey)==4)
      {
        layers.erase(3);
      }
      else if(TrkrDefs::getLayer(ckey)==5 || TrkrDefs::getLayer(ckey)==6)
      {
        layers.erase(4);
      }
      else
      {
        layers.erase(TrkrDefs::getLayer(ckey));
      }
    }

    if(Verbosity()>1 && seed_strobes.size()>1)
    {
      std::cout << "WARNING: seed MVTX strobes not all equal, something probably went wrong earlier! Strobe cut will be ignored when propagating this seed." << std::endl;
    }

    if(Verbosity()>3)
    {
      std::cout << "layers already covered: ";
      for(auto iter = seed.begin_cluster_keys(); iter != seed.end_cluster_keys(); ++iter) std::cout << (size_t)TrkrDefs::getLayer(*iter) << ", ";
      std::cout << std::endl;

      std::cout << "layers to propagate: ";
      for(const auto& layer : layers) std::cout << layer << ", ";
      std::cout << std::endl;
    }

    std::vector<pointKey> closeClusters;

    // get fit parameters in TrkFitUtils vector format
    std::vector<float> fitpars;
    fitpars.push_back(1./fabs(seed.get_qOverR())); // radius of curvature
    fitpars.push_back(seed.get_X0());
    fitpars.push_back(seed.get_Y0());
    fitpars.push_back(seed.get_slope());
    fitpars.push_back(seed.get_Z0());
    
    // define width of initial search box based on max deviation of track at a radius of 12cm
    constexpr float maxRadius = 14.;
    const TrackFitUtils::circle_circle_intersection_output_t outer_xypt = TrackFitUtils::circle_circle_intersection(maxRadius,fitpars[0],fitpars[1],fitpars[2]);
    const double& xplus = std::get<0>(outer_xypt);
    const double& yplus = std::get<1>(outer_xypt);
    const double& xminus = std::get<2>(outer_xypt);
    const double& yminus = std::get<3>(outer_xypt);

    if(std::isnan(xplus) || std::isnan(yplus) || std::isnan(xminus) || std::isnan(yminus))
    {
      if(Verbosity()>3)
      {
        std::cout << "projection to outer radius failed, skipping layer" << std::endl;
      }
      continue;
    }

    const float z_outer = seed.get_Z0() + seed.get_slope()*maxRadius;

    if(Verbosity()>3)
    {
      std::cout << "circle intersection solutions: (" << xplus << ", " << yplus << "), (" << xminus << ", " << yminus << ")" << std::endl;
      std::cout << "projected Z: " << z_outer << std::endl;
    }

    // circle-circle intersection solution closest to outermost cluster is taken as the correct one
    const Acts::Vector3& last_clusterpos = clusterpos.back();
    const float dist_plus = sqrt((xplus-last_clusterpos.x())*(xplus-last_clusterpos.x()) + (yplus-last_clusterpos.y())*(yplus-last_clusterpos.y()));
    const float dist_minus = sqrt((xminus-last_clusterpos.x())*(xminus-last_clusterpos.x()) + (yminus-last_clusterpos.y())*(yminus-last_clusterpos.y()));

    float x_outer;
    float y_outer;
    if(dist_plus < dist_minus)
    {
      x_outer = xplus;
      y_outer = yplus;
    }
    else
    {
      x_outer = xminus;
      y_outer = yminus;
    }

    const float phi_outer = atan2(y_outer,x_outer);

    const float phi_min = std::min(seed.get_phi(),phi_outer);
    const float phi_max = std::max(seed.get_phi(),phi_outer);
    const float z_min = std::min(seed.get_Z0(),z_outer);
    const float z_max = std::max(seed.get_Z0(),z_outer);

    if(Verbosity()>3)
    {
      std::cout << "outer track position: (" << x_outer << ", " << y_outer << ", " << z_outer << ")" << std::endl;
      std::cout << "phi range: " << phi_min << " - " << phi_max << std::endl;
      std::cout << "z range: " << z_min << " - " << z_max << std::endl;
    }

    for(const size_t layer : layers)
    {
      if(Verbosity()>3)
      {
        std::cout << "------------------------------------" << std::endl;
        std::cout << "layer " << layer << ": " << std::endl;
        std::cout << "current seed:" << std::endl;
        seed.identify();
      }
      const size_t l_index = layer - _start_layer;

      QueryTree(_rtrees[l_index],
                phi_min-dphi_per_layer[layer],
                z_min-dZ_per_layer[layer],
                phi_max+dphi_per_layer[layer],
                z_max+dZ_per_layer[layer],
                closeClusters);

      if(Verbosity()>3)
      {
        std::cout << "found " << closeClusters.size() << " close clusters" << std::endl;
      }

    }
  
    // for all close clusters, find the closest one by DCA, within max-cut bounds

    TrkrDefs::cluskey best_added = 0;
    float best_dca3d = 1e9;

    for(const pointKey& closeCluster : closeClusters)
    {

      if(TrkrDefs::getTrkrId(closeCluster.second) == TrkrDefs::mvtxId && seed_strobes.size()==1
         && !ClusterTimesAreCompatible(TrkrDefs::mvtxId,*(seed_strobes.begin()),closeCluster.second))
      {
        continue;
      }
      
      if(_require_INTT_consistency && seed.get_crossing() != SHRT_MAX && TrkrDefs::getTrkrId(closeCluster.second) == TrkrDefs::inttId
         && !ClusterTimesAreCompatible(TrkrDefs::inttId,seed.get_crossing(),closeCluster.second))
      {
        continue;
      }
      
      const Acts::Vector3 closeclusterpos = globalPositions.at(closeCluster.second);
      const Acts::Vector3 pca = TrackFitUtils::get_helix_pca(fitpars,closeclusterpos);
      const Acts::Vector3 residual = closeclusterpos - pca;
      const float dca_xy = sqrt(residual.x()*residual.x() + residual.y()*residual.y());
      const float dca_z = fabs(residual.z());
      const float dca_3d = sqrt(residual.x()*residual.x() + residual.y()*residual.y() + residual.z()*residual.z());
      
      if(Verbosity()>3)
      {
        std::cout << "cluster key: " << (size_t)closeCluster.second << std::endl;
        std::cout << "close cluster position: (" << closeclusterpos.x() << ", " << closeclusterpos.y() << ", " << closeclusterpos.z() << ")" << std::endl;
        std::cout << "helix PCA: (" << pca.x() << ", " << pca.y() << ", " << pca.z() << ")" << std::endl;
        std::cout << "residual: (" << residual.x() << ", " << residual.y() << ", " << residual.z() << ")" << std::endl;
        std::cout << "DCA: " << dca_xy << " (xy) " << dca_z << " (z) " << dca_3d << " (3D)" << std::endl;
      }
      const size_t layer = TrkrDefs::getLayer(closeCluster.second);
      if(dca_xy <= max_dcaxy_perlayer[layer] && dca_z <= max_dcaz_perlayer[layer] && dca_3d <= best_dca3d)
      {
        if(Verbosity()>3)
        {
          std::cout << "passed cuts" << std::endl;
        }
        best_added = closeCluster.second;
        best_dca3d = dca_3d;
      }
    }
    if(best_added != 0)
    {
      if(Verbosity()>3)
      {
        std::cout << "adding clusterkey: " << (size_t)best_added << " with 3d dca " << best_dca3d << std::endl;
      }
      seed.insert_cluster_key(best_added);
      FitSeed(seed,globalPositions);
    }
  }
}

std::vector<TrackSeed_v2> PHCASiliconSeeding::FitSeeds(const std::vector<PHCASiliconSeeding::keyList>& chains, const PHCASiliconSeeding::PositionMap& globalPositions) const
{
  std::vector<TrackSeed_v2> clean_chains;

  for (const auto& chain : chains)
  {
    if (chain.size() < 3)
    {
      continue;
    }
    if (Verbosity() > 3)
    {
      std::cout << "chain size: " << chain.size() << std::endl;
    }

    TrackSeed_v2 trackseed;
    for (const auto& key : chain)
    {
      trackseed.insert_cluster_key(key);
    }

    TrackSeedHelper::position_map_t positions;
    std::set<short> crossings;
    size_t nintt = 0;
    for (const auto& cluskey : chain)
    {
      const auto& global = globalPositions.at(cluskey);
      positions.insert(std::make_pair(cluskey,global));

      if(TrkrDefs::getTrkrId(cluskey) == TrkrDefs::inttId)
      {
        nintt++;
        crossings.insert(GetClusterTimeIndex(cluskey));
      }
    }

    TrackSeedHelper::circleFitByTaubin(&trackseed,positions,_start_layer,_end_layer);
    if(!(positions.size()==3 && nintt==2))
    {
      TrackSeedHelper::lineFit(&trackseed,positions,_start_layer,2); // don't use INTT for line fit
    }
    else
    {
      TrackSeedHelper::lineFit(&trackseed,positions,_start_layer,_end_layer); // unless we absolutely have to
    }

    trackseed.set_phi(TrackSeedHelper::get_phi(&trackseed,positions));

    if(crossings.size()==1)
    {
      trackseed.set_crossing(*crossings.begin());
    }
    else if(crossings.size()==2)
    {
      std::vector<short> crossings_vec;
      std::copy(crossings.begin(),crossings.end(),std::back_inserter(crossings_vec));
      if(abs(crossings_vec[1] - crossings_vec[0]) == 1)
      {
        trackseed.set_crossing(std::min(crossings_vec[0],crossings_vec[1]));
      }
      else
      {
        if(Verbosity()>1)
        {
          std::cout << "Warning: seed has multiple crossings within INTT clusters, setting to SHRTMAX" << std::endl;
        }
        trackseed.set_crossing(SHRT_MAX);
      }
    }
    else
    {
      if(crossings.size()>2 && Verbosity()>1)
      {
        std::cout << "Warning: seed has multiple crossings within INTT clusters, setting to SHRTMAX" << std::endl;
      }
      trackseed.set_crossing(SHRT_MAX);
    }

    clean_chains.push_back(trackseed);
    if (Verbosity() > 2)
    {
      std::cout << "pushed clean chain with " << trackseed.size_cluster_keys() << " clusters" << std::endl;
    }
  }

  return clean_chains;
}

void PHCASiliconSeeding::FitSeed(TrackSeed_v2& seed, const PositionMap& globalPositions) const
{
  if(Verbosity()>3)
  {
    std::cout << "fitting seed:" << std::endl;
    seed.identify();
  }
  TrackSeedHelper::position_map_t positions;
  std::set<short> crossings;
  size_t nintt = 0;
  for(auto clusiter = seed.begin_cluster_keys(); clusiter != seed.end_cluster_keys(); ++clusiter) 
  {
    TrkrDefs::cluskey key = *clusiter;
    const auto& global = globalPositions.at(key);
    positions.insert(std::make_pair(key,global));

    if(TrkrDefs::getTrkrId(key) == TrkrDefs::inttId)
    {
      nintt++;
      crossings.insert(GetClusterTimeIndex(key));
    }
  }

  TrackSeedHelper::circleFitByTaubin(&seed,positions,_start_layer,_end_layer);
  if(!(positions.size()==3 && nintt==2))
  {
    TrackSeedHelper::lineFit(&seed,positions,_start_layer,2); // don't use INTT for line fit
  }
  else
  {
    TrackSeedHelper::lineFit(&seed,positions,_start_layer,_end_layer); // unless we absolutely have to
  }

  seed.set_phi(TrackSeedHelper::get_phi(&seed,positions));

  if(crossings.size()==1)
  {
    seed.set_crossing(*crossings.begin());
  }
  else if(crossings.size()==2)
  {
    std::vector<short> crossings_vec;
    std::copy(crossings.begin(),crossings.end(),std::back_inserter(crossings_vec));
    if(abs(crossings_vec[1] - crossings_vec[0]) == 1)
    {
      seed.set_crossing(std::min(crossings_vec[0],crossings_vec[1]));
    }
    else
    {
      if(Verbosity()>1)
      {
        std::cout << "Warning: seed has multiple crossings within INTT clusters, setting to SHRTMAX" << std::endl;
      }
      seed.set_crossing(SHRT_MAX);
    }
  }
  else
  {
    if(crossings.size()>2 && Verbosity()>1)
    {
      std::cout << "Warning: seed has multiple crossings within INTT clusters, setting to SHRTMAX" << std::endl;
    }
    seed.set_crossing(SHRT_MAX);
  }
  if(Verbosity()>3)
  {
    std::cout << "after fit:" << std::endl;
    seed.identify();
  }
}

void PHCASiliconSeeding::publishSeeds(const std::vector<TrackSeed_v2>& seeds) const
{
  for (const auto& seed : seeds)
  {
    auto pseed = std::make_unique<TrackSeed_v2>(seed);
    if (Verbosity() > 4)
    {
      pseed->identify();
    }
    _track_map->insert(pseed.get());
  }
}

int PHCASiliconSeeding::Setup(PHCompositeNode* topNode)  // This is called by ::InitRun
{
  if(Verbosity()>0)
  {
    std::cout << "Called Setup" << std::endl;
  }
  if (Verbosity() > 0)
  {
    std::cout << "topNode:" << topNode << std::endl;
  }
  PHTrackSeeding::set_track_map_name(trackmapname);
  PHTrackSeeding::Setup(topNode);

  // geometry initialization
  int ret = InitializeGeometry(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }

  for (size_t i = _start_layer; i <= _end_layer; ++i)
  {
    if(i>=3 && i<=6)
    {
      // INTT z resolution constrains search windows to be at least 1 cm
      // (since INTT hits can be up to 2 cm in z)
      dZ_per_layer[i] = std::max<float>(_neighbor_phi_width,1.);
      max_dcaz_perlayer[i] = std::max<float>(_propagate_max_dcaz,1.);
      max_dcaxy_perlayer[i] = _propagate_max_dcaxy;
      dphi_per_layer[i] = _neighbor_phi_width;
    }
    else
    {
      dZ_per_layer[i] = _neighbor_z_width;
      max_dcaz_perlayer[i] = _propagate_max_dcaz;
      max_dcaxy_perlayer[i] = _propagate_max_dcaxy;
      dphi_per_layer[i] = _neighbor_phi_width;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHCASiliconSeeding::End()
{
  if (Verbosity() > 0)
  {
    std::cout << "Called End " << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
