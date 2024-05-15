/*!
 *  \file PHCASeeding.cc
 *  \brief Track seeding using ALICE-style "cellular automaton" (CA) algorithm
 *  \detail
 *  \author Michael Peters & Christof Roland
 */

#include "PHCASeeding.h"
#include "ALICEKF.h"
#include "GPUTPCTrackLinearisation.h"
#include "GPUTPCTrackParam.h"

// sPHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHTimer.h>  // for PHTimer
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

// tpc distortion correction
#include <tpc/TpcDistortionCorrectionContainer.h>

#include <phfield/PHFieldConfigv2.h>

#include <ffamodules/CDBInterface.h>

// trackbase_historic includes
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrCluster.h>  // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrClusterIterationMapv1.h>
#include <trackbase/TrkrDefs.h>  // for getLayer, clu...
#include <trackbase_historic/TrackSeedContainer.h>
#include <trackbase_historic/TrackSeed_v2.h>

// ROOT includes for debugging
#include <TFile.h>
#include <TNtuple.h>

// BOOST for combi seeding
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/policies/compare.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <algorithm>
#include <algorithm>  // for find
#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <numeric>
#include <unordered_set>
#include <utility>  // for pair, make_pair
#include <vector>


//#define _DEBUG_

#if defined(_DEBUG_)
#define LogDebug(exp) \
  if (Verbosity() > 0) std::cout << "DEBUG: " << __FILE__ << ": " << __LINE__ << ": " << exp
#else
#define LogDebug(exp) (void) 0
#endif

#define LogError(exp) \
  if (Verbosity() > 0) std::cout << "ERROR: " << __FILE__ << ": " << __LINE__ << ": " << exp
#define LogWarning(exp) \
  if (Verbosity() > 0) std::cout << "WARNING: " << __FILE__ << ": " << __LINE__ << ": " << exp

/* #define FIXME(str) std::cout <<" FIXME: " << str << std::endl; */

//#define _DEBUG_

// end


// apparently there is no builtin STL hash function for a std::array
// so to use std::unordered_set (essentially a hash table), we have to make our own hasher

namespace std
{
  template <typename T, size_t N>
  struct hash<std::array<T, N>>
  {
    using argument_type = std::array<T, N>;
    using result_type = size_t;

    result_type operator()(const argument_type& a) const
    {
      hash<T> hasher;
      result_type h = 0;
      for (result_type i = 0; i < N; ++i)
      {
        h = h * 31 + hasher(a[i]);
      }
      return h;
    }
  };
  template <typename A, typename B>
  struct hash<pair<A, B>>
  {
    using argument_type = pair<A, B>;
    using result_type = size_t;

    result_type operator()(const argument_type& a) const
    {
      hash<A> hashA;
      hash<B> hashB;
      return (hashA(a.first) * 31 + hashB(a.second));
    }
  };
}  // namespace std

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

  /// pseudo rapidity of Acts::Vector3
  inline double get_eta(const Acts::Vector3& position)
  {
    const double norm = std::sqrt(square(position.x()) + square(position.y()) + square(position.z()));
    return std::log((norm + position.z()) / (norm - position.z())) / 2;
  }

  ///@name utility convertion functions
  //@{

  PHCASeeding::coordKey fromPointKey(const PHCASeeding::pointKey& p)
  {
    return std::make_pair(std::array<float, 2>({p.first.get<0>(), p.first.get<1>()}), p.second);
  }

  /* std::vector<PHCASeeding::coordKey> fromPointKey(const std::vector<PHCASeeding::pointKey>& p) */
  /* { */
  /*   std::vector<PHCASeeding::coordKey> output; */
  /*   output.resize(p.size()); */
  /*   std::transform(p.begin(), p.end(), std::back_inserter(output), [](const PHCASeeding::pointKey& point) */
  /*                  { return fromPointKey(point); }); */
  /*   return output; */
  /* } */
  //@}

  inline double breaking_angle(double x1, double y1, double z1, double x2, double y2, double z2)
  {
    double l1 = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
    double l2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2);
    double sx = (x1 / l1 + x2 / l2);
    double sy = (y1 / l1 + y2 / l2);
    double sz = (z1 / l1 + z2 / l2);
    double dx = (x1 / l1 - x2 / l2);
    double dy = (y1 / l1 - y2 / l2);
    double dz = (z1 / l1 - z2 / l2);
    return 2 * atan2(sqrt(dx * dx + dy * dy + dz * dz), sqrt(sx * sx + sy * sy + sz * sz));
  }

}  // namespace

// using namespace ROOT::Minuit2;
namespace bgi = boost::geometry::index;

PHCASeeding::PHCASeeding(
    const std::string& name,
    unsigned int start_layer,
    unsigned int end_layer,
    unsigned int min_nhits_per_cluster,
    unsigned int min_clusters_per_track,
    unsigned int nlayers_maps,
    unsigned int nlayers_intt,
    unsigned int nlayers_tpc,
    float neighbor_phi_width,
    float neighbor_eta_width,
    float maxSinPhi,
    float cosTheta_limit)
  : PHTrackSeeding(name)
  , _nlayers_maps(nlayers_maps)
  , _nlayers_intt(nlayers_intt)
  , _nlayers_tpc(nlayers_tpc)
  , _start_layer(start_layer)
  , _end_layer(end_layer)
  , _min_nhits_per_cluster(min_nhits_per_cluster)
  , _min_clusters_per_track(min_clusters_per_track)
  , _neighbor_phi_width(neighbor_phi_width)
  , _neighbor_eta_width(neighbor_eta_width)
  , _max_sin_phi(maxSinPhi)
  , _cosTheta_limit(cosTheta_limit)
{
}

int PHCASeeding::InitializeGeometry(PHCompositeNode* topNode)
{
  tGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!tGeometry)
  {
    std::cout << PHWHERE << "No acts tracking geometry, can't proceed" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

Acts::Vector3 PHCASeeding::getGlobalPosition(TrkrDefs::cluskey key, TrkrCluster* cluster) const
{
  // get global position from Acts transform
  auto globalpos = tGeometry->getGlobalPosition(key, cluster);

  // check if TPC distortion correction are in place and apply
  if (m_dcc && !(_pp_mode))
  {
    globalpos = m_distortionCorrection.get_corrected_position(globalpos, m_dcc);
  }

  return globalpos;
}

void PHCASeeding::QueryTree(const bgi::rtree<PHCASeeding::pointKey, bgi::quadratic<16>>& rtree, double phimin, double etamin, double phimax, double etamax, std::vector<pointKey>& returned_values) const
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
    rtree.query(bgi::intersects(box(point(phimin, etamin), point(2 * M_PI, etamax))), std::back_inserter(returned_values));
    rtree.query(bgi::intersects(box(point(0., etamin), point(phimax, etamax))), std::back_inserter(returned_values));
  }
  else
  {
    rtree.query(bgi::intersects(box(point(phimin, etamin), point(phimax, etamax))), std::back_inserter(returned_values));
  }
}

std::pair<PHCASeeding::PositionMap, PHCASeeding::keysPerLayer> PHCASeeding::FillGlobalPositions() {
  keysPerLayer ckeys;
  PositionMap cachedPositions;
  cachedPositions.reserve(_cluster_map->size());  //avoid resizing mid-execution

  for (const auto& hitsetkey : _cluster_map->getHitSetKeys(TrkrDefs::TrkrId::tpcId))
  {
    auto range = _cluster_map->getClusters(hitsetkey);
    for (auto clusIter = range.first; clusIter != range.second; ++clusIter)
    {
      TrkrDefs::cluskey ckey = clusIter->first;
      TrkrCluster* cluster = clusIter->second;
      unsigned int layer = TrkrDefs::getLayer(ckey);
      if (layer < _start_layer || layer >= _end_layer)
      {
        if (Verbosity() > 0)
        {
          std::cout << "layer: " << layer << std::endl;
        }
        continue;
      }
      if (_iteration_map != nullptr && _n_iteration > 0)
      {
        if (_iteration_map->getIteration(ckey) > 0)
        {
          continue;  // skip hits used in a previous iteration
        }
      }

      // get global position, convert to Acts::Vector3 and store in map
      const Acts::Vector3 globalpos_d = getGlobalPosition(ckey, cluster);
#if defined(_CLUSTER_LOG_TUPOUT_)
      _tupclus_all->Fill(_nevent, globalpos_d.x(), globalpos_d.y(), globalpos_d.z());
#endif

      const Acts::Vector3 globalpos = {globalpos_d.x(), globalpos_d.y(), globalpos_d.z()};
      cachedPositions.insert(std::make_pair(ckey, globalpos));
      ckeys[layer-_start_layer].push_back(ckey);
    }
  }
  return std::make_pair(cachedPositions, ckeys);
}

std::vector<PHCASeeding::coordKey> PHCASeeding::FillTree(bgi::rtree<PHCASeeding::pointKey,bgi::quadratic<16>>& _rtree, const PHCASeeding::keyList& ckeys, const PHCASeeding::PositionMap& globalPositions, const int layer) 
{
  // Fill _rtree with the clusters in ckeys; remove duplicates, and return a vector of the coordKeys
  // Note that layer is only used for a cout statement
  int n_dupli = 0;
  std::vector<coordKey> coords;
  _rtree.clear();
  /* _rtree.reserve(ckeys.size()); */
  for (const auto& ckey : ckeys) {
    const auto& globalpos_d = globalPositions.at(ckey);
    const double clus_phi = get_phi(globalpos_d);
    const double clus_eta = get_eta(globalpos_d);
    if (Verbosity() > 0)
    {
      /* int layer = TrkrDefs::getLayer(ckey); */
      std::cout << "Found cluster " << ckey << " in layer " << layer << std::endl;
    }
    std::vector<pointKey> testduplicate;
    QueryTree(_rtree, clus_phi - 0.00001, clus_eta - 0.00001, clus_phi + 0.00001, clus_eta + 0.00001, testduplicate);
    if (!testduplicate.empty())
    {
      ++n_dupli;
      continue;
    }
    coords.push_back({{static_cast<float>(clus_phi), static_cast<float>(clus_eta)},ckey});
    t_fill->restart();
    _rtree.insert(std::make_pair(point(clus_phi, globalpos_d.z()), ckey));
    t_fill->stop();
  }
  if (Verbosity() > 1)
  {
      std::cout << "nhits in layer("<<layer<<"): " << coords.size() << std::endl;
  }
  if (Verbosity() > 0)
  {
    std::cout << "fill time: " << t_fill->get_accumulated_time() / 1000. << " sec" << std::endl;
  }
  if (Verbosity() > 0)
  {
    std::cout << "number of duplicates : " << n_dupli << std::endl;
  }
  return coords;
}

int PHCASeeding::Process(PHCompositeNode* /*topNode*/)
{
#if defined(_CLUSTER_LOG_TUPOUT_)
  _nevent += 1;
#endif
  if (Verbosity() > 1)
  {
    std::cout << " Process...  " << std::endl;
  }
  //  TFile fpara("CA_para.root", "RECREATE");
  if (_n_iteration > 0)
  {
    if (!_iteration_map)
    {
      std::cerr << PHWHERE << "Cluster Iteration Map missing, aborting." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  t_seed->restart();
  t_makebilinks->restart();

  PositionMap globalPositions;
  keysPerLayer ckeys;
  std::tie(globalPositions, ckeys) = FillGlobalPositions();

  t_seed->stop();
  if (Verbosity() > 0)
  {
    std::cout << "Initial RTree fill time: " << t_seed->get_accumulated_time() / 1000 << " s" << std::endl;
  }
  t_seed->restart();
  int numberofseeds = 0;
  numberofseeds += FindSeedsWithMerger(globalPositions,ckeys);
  t_seed->stop();
  if (Verbosity() > 0)
  {
    std::cout << "number of seeds " << numberofseeds << std::endl;
  }
  if (Verbosity() > 0)
  {
    std::cout << "Kalman filtering time: " << t_seed->get_accumulated_time() / 1000 << " s" << std::endl;
  }
  //  fpara.cd();
  //  fpara.Close();
  //  if(Verbosity()>0) std::cout << "fpara OK\n";
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHCASeeding::FindSeedsWithMerger(const PHCASeeding::PositionMap& globalPositions, const PHCASeeding::keysPerLayer& ckeys)
{
  t_seed->restart();

  keyLinksPerLayer bilinks = CreateBiLinks(globalPositions, ckeys);
  if (Verbosity() > 0) {
    t_makebilinks->stop();
    std::cout << "Time to make bilinks: " << t_makebilinks->elapsed() / 1000 << " s" << std::endl;
  }
  // std::pair<ckeylinks_per_layer, ckeylinks_per_layer> = CreateLinks(fromPointKey(allClusters), globalPositions); // put the bilinks finder into the end of create links --> don't have to save all the rows
  // std::vector<std::vector<keylink>> biLinks = FindBiLinks(links.first, links.second);
  

  if (Verbosity() > 0) { t_makeseeds->restart(); }
  std::vector<keyList> trackSeedKeyLists = FollowBiLinks(bilinks, globalPositions);
  if (Verbosity() > 0 )
  {
    t_makeseeds->stop();
    std::cout << "Time to make seeds: " << t_makeseeds->elapsed() / 1000 << " s" << std::endl;
  }
  std::vector<TrackSeed_v2> seeds = RemoveBadClusters(trackSeedKeyLists, globalPositions);

  publishSeeds(seeds);
  return seeds.size();
}

PHCASeeding::keyLinksPerLayer PHCASeeding::CreateBiLinks(const PHCASeeding::PositionMap& globalPositions, const PHCASeeding::keysPerLayer& ckeys)
{
  double cluster_find_time = 0;
  double rtree_query_time = 0;
  double transform_time = 0;
  double compute_best_angle_time = 0;
  double set_insert_time = 0;
  keyLinksPerLayer bilinks;

  std::array<std::vector<coordKey>,3> coord_arr;
  for (int layer=0;layer<3;++layer) {
    coord_arr[layer] = FillTree(_rtrees[layer], ckeys[layer], globalPositions, layer);
  }
  std::unordered_set<keyLink> previous_layer_above {};
  for (int layer_index=1; layer_index<(_NLAYERS_TPC-1);++layer_index) {
    // these lines of code will rotates through all three _rtree's in the array,
    // where the old upper becomes the new middle, the old middle the new lower,
    // and the old lower drops out and that _rtree is filled with the new upper
    int index_below = (layer_index-1)%3;
    int index       = (layer_index  )%3;
    int index_above = (layer_index+1)%3;

    auto& _rtree_below = _rtrees[index_below];
    auto& coord = coord_arr[index];
    auto& _rtree_above = _rtrees[index_above];
    // find the above and below links. Match the below links with the 
    // previous_layer_above
    std::vector<keyLink> belowLinks;
    std::vector<keyLink> aboveLinks;
    for (const auto& StartCluster : coord)
  {
    double StartPhi = StartCluster.first[0];
    const auto& globalpos = globalPositions.at(StartCluster.second);
    double StartX = globalpos(0);
    double StartY = globalpos(1);
    double StartZ = globalpos(2);
    t_seed->stop();
    cluster_find_time += t_seed->elapsed();
    t_seed->restart();
    LogDebug(" starting cluster:" << std::endl);
    LogDebug(" z: " << StartZ << std::endl);
    LogDebug(" phi: " << StartPhi << std::endl);
    // LogDebug(" layer: " << StartLayer << std::endl);

    std::vector<pointKey> ClustersAbove;
    std::vector<pointKey> ClustersBelow;

    QueryTree(_rtree_below,
              StartPhi - _neighbor_phi_width,
              StartZ - _neighbor_eta_width,
              // (double) StartLayer - 1.5,
              StartPhi + _neighbor_phi_width,
              StartZ + _neighbor_eta_width,
              // (double) StartLayer - 0.5,
              ClustersBelow);

    QueryTree(_rtree_above,
              StartPhi - _neighbor_phi_width,
              StartZ - _neighbor_eta_width,
              // (double) StartLayer + 0.5,
              StartPhi + _neighbor_phi_width,
              StartZ + _neighbor_eta_width,
              // (double) StartLayer + 1.5,
              ClustersAbove);

    t_seed->stop();
    rtree_query_time += t_seed->elapsed();
    t_seed->restart();
    LogDebug(" entries in below layer: " << ClustersBelow.size() << std::endl);
    LogDebug(" entries in above layer: " << ClustersAbove.size() << std::endl);
    std::vector<std::array<double, 3>> delta_below;
    std::vector<std::array<double, 3>> delta_above;
    delta_below.clear();
    delta_above.clear();
    delta_below.resize(ClustersBelow.size());
    delta_above.resize(ClustersAbove.size());
    // calculate (delta_eta, delta_phi) vector for each neighboring cluster

    std::transform(ClustersBelow.begin(), ClustersBelow.end(), delta_below.begin(),
                   [&](pointKey BelowCandidate)
                   {
	const auto& belowpos = globalPositions.at(BelowCandidate.second);
        return std::array<double,3>{belowpos(0)-StartX,
          belowpos(1)-StartY,
          belowpos(2)-StartZ}; });

    std::transform(ClustersAbove.begin(), ClustersAbove.end(), delta_above.begin(),
                   [&](pointKey AboveCandidate)
                   {
	const auto& abovepos = globalPositions.at(AboveCandidate.second);
        return std::array<double,3>{abovepos(0)-StartX,
          abovepos(1)-StartY,
          abovepos(2)-StartZ}; });
    t_seed->stop();
    transform_time += t_seed->elapsed();
    t_seed->restart();

    // find the three clusters closest to a straight line
    // (by maximizing the cos of the angle between the (delta_eta,delta_phi) vectors)
    double maxCosPlaneAngle = -0.95;
    // double minSumLengths = 1e9;
    std::vector<coordKey> bestBelowClusters;
    std::vector<coordKey> bestAboveClusters;
    for (size_t iAbove = 0; iAbove < delta_above.size(); ++iAbove)
    {
      for (size_t iBelow = 0; iBelow < delta_below.size(); ++iBelow)
      {
        //        double dotProduct = delta_below[iBelow][0]*delta_above[iAbove][0]+delta_below[iBelow][1]*delta_above[iAbove][1]+delta_below[iBelow][2]*delta_above[iAbove][2];
        //        double belowLength = sqrt(delta_below[iBelow][0]*delta_below[iBelow][0]+delta_below[iBelow][1]*delta_below[iBelow][1]+delta_below[iBelow][2]*delta_below[iBelow][2]);
        //        double aboveLength = sqrt(delta_above[iAbove][0]*delta_above[iAbove][0]+delta_above[iAbove][1]*delta_above[iAbove][1]+delta_above[iAbove][2]*delta_above[iAbove][2]);
        double angle = breaking_angle(
            delta_below[iBelow][0],
            delta_below[iBelow][1],
            delta_below[iBelow][2],
            delta_above[iAbove][0],
            delta_above[iAbove][1],
            delta_above[iAbove][2]);
        if (cos(angle) < maxCosPlaneAngle)
        {
          // maxCosPlaneAngle = cos(angle);
          // minSumLengths = belowLength+aboveLength;
          // bestBelowCluster = fromPointKey(ClustersBelow[iBelow]);
          // bestAboveCluster = fromPointKey(ClustersAbove[iAbove]);
          bestBelowClusters.push_back(fromPointKey(ClustersBelow[iBelow]));
          bestAboveClusters.push_back(fromPointKey(ClustersAbove[iAbove]));
        }
      }
    }
    /*
    if(mode == skip_layers::on)
    {
      if(maxCosPlaneAngle > _cosTheta_limit)
      {
        // if no triplet is sufficiently linear, then it's likely that there's a missing cluster
        // repeat search but skip one layer below
        std::vector<pointKey> clustersTwoLayersBelow;
        QueryTree(_rtree,
                StartPhi-_neighbor_phi_width,
                StartEta-_neighbor_eta_width,
                (double) StartLayer - 2.5,
                StartPhi+_neighbor_phi_width,
                StartEta+_neighbor_eta_width,
                (double) StartLayer - 1.5,
                clustersTwoLayersBelow);
        std::vector<std::array<double,3>> delta_2below;
        delta_2below.clear();
        delta_2below.resize(clustersTwoLayersBelow.size());
        std::transform(clustersTwoLayersBelow.begin(),clustersTwoLayersBelow.end(),delta_2below.begin(),
          [&](pointKey BelowCandidate){
            const auto& belowpos = globalPositions.at(BelowCandidate.second);
            return std::array<double,3>{(belowpos(0))-StartX,
              (belowpos(1))-StartY,
              (belowpos(2))-StartZ};});
        for(size_t iAbove = 0; iAbove<delta_above.size(); ++iAbove)
        {
          for(size_t iBelow = 0; iBelow<delta_2below.size(); ++iBelow)
          {
            double dotProduct = delta_2below[iBelow][0]*delta_above[iAbove][0]+delta_2below[iBelow][1]*delta_above[iAbove][1]+delta_2below[iBelow][2]*delta_above[iAbove][2];
            double belowSqLength = sqrt(delta_2below[iBelow][0]*delta_2below[iBelow][0]+delta_2below[iBelow][1]*delta_2below[iBelow][1]+delta_2below[iBelow][2]*delta_2below[iBelow][2]);
            double aboveSqLength = sqrt(delta_above[iAbove][0]*delta_above[iAbove][0]+delta_above[iAbove][1]*delta_above[iAbove][1]+delta_above[iAbove][2]*delta_above[iAbove][2]);
            double cosPlaneAngle = dotProduct / (belowSqLength*aboveSqLength);
            if(cosPlaneAngle < maxCosPlaneAngle)
            {
              maxCosPlaneAngle = cosPlaneAngle;
              bestBelowCluster = fromPointKey(clustersTwoLayersBelow[iBelow]);
              bestAboveCluster = fromPointKey(ClustersAbove[iAbove]);
            }
          }
        }
        // if no triplet is STILL sufficiently linear, then do the same thing, but skip one layer above
        if(maxCosPlaneAngle > _cosTheta_limit)
        {
          std::vector<pointKey> clustersTwoLayersAbove;
          QueryTree(_rtree,
                  StartPhi-_neighbor_phi_width,
                  StartEta-_neighbor_eta_width,
                  (double) StartLayer + 1.5,
                  StartPhi+_neighbor_phi_width,
                  StartEta+_neighbor_eta_width,
                  (double) StartLayer + 2.5,
                  clustersTwoLayersAbove);
          std::vector<std::array<double,3>> delta_2above;
          delta_2above.clear();
          delta_2above.resize(clustersTwoLayersAbove.size());
          std::transform(clustersTwoLayersAbove.begin(),clustersTwoLayersAbove.end(),delta_2above.begin(),
            [&](pointKey AboveCandidate){
              const auto& abovepos = globalPositions.at(AboveCandidate.second);
              return std::array<double,3>{(abovepos(0))-StartX,
                (abovepos(1))-StartY,
                (abovepos(2))-StartZ};});
          for(size_t iAbove = 0; iAbove<delta_2above.size(); ++iAbove)
          {
            for(size_t iBelow = 0; iBelow<delta_below.size(); ++iBelow)
            {
              double dotProduct = delta_below[iBelow][0]*delta_2above[iAbove][0]+delta_below[iBelow][1]*delta_2above[iAbove][1]+delta_below[iBelow][2]*delta_2above[iAbove][2];
              double belowSqLength = sqrt(delta_below[iBelow][0]*delta_below[iBelow][0]+delta_below[iBelow][1]*delta_below[iBelow][1]+delta_below[iBelow][2]*delta_below[iBelow][2]);
              double aboveSqLength = sqrt(delta_2above[iAbove][0]*delta_2above[iAbove][0]+delta_2above[iAbove][1]*delta_2above[iAbove][1]+delta_2above[iAbove][2]*delta_2above[iAbove][2]);
              double cosPlaneAngle = dotProduct / (belowSqLength*aboveSqLength);
              if(cosPlaneAngle < maxCosPlaneAngle)
              {
                maxCosPlaneAngle = cosPlaneAngle;
                bestBelowCluster = fromPointKey(ClustersBelow[iBelow]);
                bestAboveCluster = fromPointKey(clustersTwoLayersAbove[iAbove]);
              }
            }
          }
        }
      }
    }
    */
    t_seed->stop();
    compute_best_angle_time += t_seed->elapsed();
    t_seed->restart();
    // int layer_index = StartLayer - (_nlayers_intt + _nlayers_maps);
    // if(bestBelowCluster.second != 0) belowLinks[layer_index].insert(keylink{{*StartCluster,bestBelowCluster}});
    // if(bestAboveCluster.second != 0) aboveLinks[layer_index].insert(keylink{{*StartCluster,bestAboveCluster}});
    for (auto cluster : bestBelowClusters)
    {
      belowLinks.emplace_back(std::make_pair(StartCluster.second, cluster.second));
#if defined(_CLUSTER_LOG_TUPOUT_)
      auto& v = globalPositions.at(cluster.second);
      _tupclus_links->Fill(_nevent, TrkrDefs::getLayer(cluster.second), 0, v.x(), v.y(), v.z());
#endif
    }
    for (auto cluster : bestAboveClusters)
    {
      aboveLinks.emplace_back(std::make_pair(StartCluster.second, cluster.second));
#if defined(_CLUSTER_LOG_TUPOUT_)
      auto& v = globalPositions.at(cluster.second);
      _tupclus_links->Fill(_nevent, TrkrDefs::getLayer(cluster.second), 1, v.x(), v.y(), v.z());
#endif
    }
    t_seed->stop();
    set_insert_time += t_seed->elapsed();
    t_seed->restart();
    LogDebug(" max collinearity: " << maxCosPlaneAngle << std::endl);
  } // end loop filling belowLinks and aboveLinks
    // now find and fill bilinks for this layer
    if (!previous_layer_above.empty()) {
      for (const auto& b_link : belowLinks)
      {
        keyLink revkey = {b_link.second, b_link.first};
        if (previous_layer_above.find(revkey) != previous_layer_above.end())
        {
          bilinks[layer_index].push_back(b_link);
#if defined(_CLUSTER_LOG_TUPOUT_)
        auto& v0 = globalPositions.at(b_link.first);
        _tupclus_bilinks->Fill(_nevent, TrkrDefs::getLayer(b_link.first),0,v0.x(),v0.y(),v0.z());
#endif
        }
      }
    }

    // sorted bilinks by their second entry (the inner layer pointed to)
    std::sort(bilinks[layer_index].begin(), bilinks[layer_index].end(), 
        [](const keyLink& a, const keyLink& b) { return a.second<b.second; });

    // rotate data for the inspection of the next three rows
    if (layer_index < (_NLAYERS_TPC-2)) { // only do when there is at least one loop left
      previous_layer_above = std::unordered_set<keyLink>(aboveLinks.begin(), aboveLinks.end());
      // the _rtree and coord_arr in the old index_below are to be filled with the next
      // row above's data for the next layer's loop (i.e. at layer_index+2)
      coord_arr[index_below] = FillTree(_rtrees[index_below], ckeys[layer_index+2], globalPositions, layer_index+2);
    }
  }

  t_seed->stop();
  if (Verbosity() > 0)
  {
    std::cout << "triplet forming time: " << t_seed->get_accumulated_time() / 1000 << " s" << std::endl;
    std::cout << "starting cluster setup: " << cluster_find_time / 1000 << " s" << std::endl;
    std::cout << "RTree query: " << rtree_query_time / 1000 << " s" << std::endl;
    std::cout << "Transform: " << transform_time / 1000 << " s" << std::endl;
    std::cout << "Compute best triplet: " << compute_best_angle_time / 1000 << " s" << std::endl;
    std::cout << "Set insert: " << set_insert_time / 1000 << " s" << std::endl;
  }
  t_seed->restart();

  return bilinks; //std::make_pair(belowLinks, aboveLinks);
}

/* ckeylinks_per_layer PHCASeeding::FindBiLinks(const ckeylinks_per_layer& belowLinks, const ckeylinks_per_layer& aboveLinks) const */
/* { */
/*   // remove all triplets for which there isn't a mutual association between two clusters */
/*   ckeylinks_per_layer bidirectionalLinks; */
/*   bidirectionalLinks.resize(_nlayers_tpc); */


/*   for (int layer = _nlayers_tpc - 1; layer > 0; --layer) */
/*   { */
/*     for (const auto& belowLink : belowLinks[layer]) //= belowLinks[layer].begin(); belowLink != belowLinks[layer].end(); ++belowLink) */
/*     { */
/*       if ((*belowLink)[1].second == 0) */
/*       { */
/*         continue; */
/*       } */
/*       unsigned int end_layer_index = TrkrDefs::getLayer((*belowLink)[1].second) - (_nlayers_intt + _nlayers_maps); */
/*       keylink reversed = {(*belowLink)[1], (*belowLink)[0]}; */
/*       auto sameAboveLinkExists = aboveLinks[end_layer_index].find(reversed); */
/*       if (sameAboveLinkExists != aboveLinks[end_layer_index].end()) */
/*       { */
/*         bidirectionalLinks[layer].push_back((*belowLink)); */
/*       } */
/*     } */
/*   } */
/*   t_seed->stop(); */
/*   if (Verbosity() > 0) */
/*   { */
/*     std::cout << "bidirectional link forming time: " << t_seed->get_accumulated_time() / 1000 << " s" << std::endl; */
/*   } */
/*   t_seed->restart(); */

/*   return bidirectionalLinks; */
/* } */

double PHCASeeding::getMengerCurvature(TrkrDefs::cluskey a, TrkrDefs::cluskey b, TrkrDefs::cluskey c, const PHCASeeding::PositionMap& globalPositions) const
{
  // Menger curvature = 1/R for circumcircle of triangle formed by most recent three clusters
  // We use here 1/R = 2*sin(breaking angle)/(hypotenuse of triangle)
  auto& a_pos = globalPositions.at(a);
  auto& b_pos = globalPositions.at(b);
  auto& c_pos = globalPositions.at(c);
  double hypot_length = sqrt(square<double>(c_pos.x() - a_pos.x()) + square<double>(c_pos.y() - a_pos.y()) + square<double>(c_pos.z() - a_pos.z()));
  double break_angle = breaking_angle(
      a_pos.x() - b_pos.x(),
      a_pos.y() - b_pos.y(),
      a_pos.z() - b_pos.z(),
      c_pos.x() - b_pos.x(),
      c_pos.y() - b_pos.y(),
      c_pos.z() - b_pos.z());
  return 2 * sin(break_angle) / hypot_length;
}

std::vector<PHCASeeding::keyList> PHCASeeding::FollowBiLinks(const PHCASeeding::keyLinksPerLayer& bidirectionalLinks, const PHCASeeding::PositionMap& globalPositions) const
{

#if defined(_CLUSTER_LOG_TUPOUT_)
  for (auto& list: bidirectionalLinks) {
    for (auto& linkpair : list) {
        auto& v0 = globalPositions.at(linkpair.first);
        _tupclus_bilinks->Fill(_nevent, TrkrDefs::getLayer(linkpair.first),0,v0.x(),v0.y(),v0.z());
        auto& v1 = globalPositions.at(linkpair.second);
        _tupclus_bilinks->Fill(_nevent, TrkrDefs::getLayer(linkpair.second),0,v1.x(),v1.y(),v1.z());
    }
  }
#endif
  // follow bidirectional links to form lists of cluster keys
  // (to be fitted for track seed parameters)
  std::vector<keyList> trackSeedPairs;
  // get starting cluster keys, create a keyList for each
  // (only check last element of each pair because we start from the outer layers and go inward)
  for (unsigned int layer = 0; layer < _nlayers_tpc - 1; ++layer)
  {
    
    for (const auto& startCand : bidirectionalLinks[layer]) {
      //see if any link in an above layer is pointing to the starCand
      const auto& l_above = bidirectionalLinks[layer+1];
      if (!std::binary_search(l_above.begin(), l_above.end(), startCand,
            [](const keyLink& a, const keyLink& b) { return a.second < b.first; })) 
      {
        trackSeedPairs.push_back({startCand.first, startCand.second});
      }
    }
  }



/*     for (auto startCand = bidirectionalLinks[layer].begin(); startCand != bidirectionalLinks[layer].end(); ++startCand) */
/*     { */
/*       bool has_above_link = false; */
/*       unsigned int imax = 1; */
/*       if (layer == _nlayers_tpc - 2) */
/*       { */
/*         imax = 1; */
/*       } */
/*       for (unsigned int i = 1; i <= imax; i++) */
/*       { */
/*         has_above_link = has_above_link || std::any_of(bidirectionalLinks[layer + i].begin(), bidirectionalLinks[layer + i].end(), [&](keyLink k) */
/*           { return startCand->first == k.second; }); */
/*       } */
/*       //      for(std::vector<keylink>::iterator testlink = bidirectionalLinks[layer+1].begin(); !has_above_link && (testlink != bidirectionalLinks[layer+1].end()); ++testlink) */
/*       //      { */
/*       //        if((*startCand) == (*testlink)) continue; */
/*       //        if((*startCand)[0] == (*testlink)[1]) has_above_link = true; */
/*       //      } */
/*       if (!has_above_link) */
/*       { */
/*         trackSeedPairs.push_back({startCand->first, startCand->second}); */
/*       } */
/*     } */
/*   } */

  // form all possible starting 3-cluster tracks (we need that to calculate curvature)
  std::vector<keyList> trackSeedKeyLists;
  for (auto& trackKeyChain : trackSeedPairs)
  {
    TrkrDefs::cluskey trackHead = trackKeyChain.back();
    unsigned int trackHead_layer = TrkrDefs::getLayer(trackHead) - _start_layer;
    for (auto& testlink : bidirectionalLinks[trackHead_layer])
    {
      if (testlink.first == trackHead)
      {
        keyList trackSeedTriplet;
        trackSeedTriplet.push_back(trackKeyChain[0]);
        trackSeedTriplet.push_back(trackKeyChain[1]);
        trackSeedTriplet.push_back(testlink.second);
        trackSeedKeyLists.push_back(trackSeedTriplet);
#if defined(_CLUSTER_LOG_TUPOUT_)
      auto& c0 = trackKeyChain[0];
      auto& v0 = globalPositions.at(c0);
      _tupclus_seeds->Fill(_nevent, TrkrDefs::getLayer(c0), 0, v0.x(), v0.y(), v0.z());
      auto& c1 = trackKeyChain[1];
      auto& v1 = globalPositions.at(c1);
      _tupclus_seeds->Fill(_nevent, TrkrDefs::getLayer(c1), 1, v1.x(), v1.y(), v1.z());
      auto& c2 = testlink.second;
      auto& v2 = globalPositions.at(c2);
      _tupclus_seeds->Fill(_nevent, TrkrDefs::getLayer(c2), 2, v2.x(), v2.y(), v2.z());
#endif
      }
    }
  }

  t_seed->stop();
  if (Verbosity() > 0)
  {
    std::cout << "starting cluster finding time: " << t_seed->get_accumulated_time() / 1000 << " s" << std::endl;
  }
  t_seed->restart();
  // assemble track cluster chains from starting cluster keys (ordered from outside in)

  // std::cout << "STARTING SEED ASSEMBLY" << std::endl;

  //  std::vector<keyList> trackSeedKeyLists;
  std::vector<keyList> tempSeedKeyLists = trackSeedKeyLists;
  trackSeedKeyLists.clear();

  while (tempSeedKeyLists.size() > 0)
  {
    if (Verbosity() > 0)
    {
      std::cout << "temp size: " << tempSeedKeyLists.size() << std::endl;
    }
    if (Verbosity() > 0)
    {
      std::cout << "final size: " << trackSeedKeyLists.size() << std::endl;
    }
    std::vector<keyList> newtempSeedKeyLists;
    for (auto& seed : tempSeedKeyLists)
    {
      TrkrDefs::cluskey trackHead = seed.back();
      unsigned int trackHead_layer = TrkrDefs::getLayer(trackHead) - (_nlayers_intt + _nlayers_maps);
      // bool no_next_link = true;
      for (auto& link : bidirectionalLinks[trackHead_layer])
      {
        if (link.first != trackHead)
        {
          continue;
        }
        auto& head_pos = globalPositions.at(trackHead);
        auto& prev_pos = globalPositions.at(seed.rbegin()[1]);
        float x1 = head_pos.x();
        float y1 = head_pos.y();
        float z1 = head_pos.z();
        float x2 = prev_pos.x();
        float y2 = prev_pos.y();
        float z2 = prev_pos.z();
        float dr_12 = sqrt(x1 * x1 + y1 * y1) - sqrt(x2 * x2 + y2 * y2);
        TrkrDefs::cluskey testCluster = link.second;
        auto& test_pos = globalPositions.at(testCluster);
        float xt = test_pos.x();
        float yt = test_pos.y();
        float zt = test_pos.z();
        float new_dr = sqrt(xt * xt + yt * yt) - sqrt(x1 * x1 + y1 * y1);
        if (fabs((z1 - z2) / dr_12 - (zt - z1) / new_dr) > _clusadd_delta_dzdr_window)
        {
          continue;
        }
        auto& third_pos = globalPositions.at(seed.rbegin()[2]);
        float x3 = third_pos.x();
        float y3 = third_pos.y();
        float dr_23 = sqrt(x2 * x2 + y2 * y2) - sqrt(x3 * x3 + y3 * y3);
        float phi1 = atan2(y1, x1);
        float phi2 = atan2(y2, x2);
        float phi3 = atan2(y3, x3);
        float dphi12 = std::fmod(phi1 - phi2, M_PI);
        float dphi23 = std::fmod(phi2 - phi3, M_PI);
        float d2phidr2 = dphi12 / (dr_12 * dr_12) - dphi23 / (dr_23 * dr_23);
        float new_dphi = std::fmod(atan2(yt, xt) - atan2(y1, x1), M_PI);
        float new_d2phidr2 = new_dphi / (new_dr * new_dr) - dphi12 / (dr_12 * dr_12);
        if (seed.size() < 6 && fabs(d2phidr2 - new_d2phidr2) < _clusadd_delta_dphidr2_window)
        {
          //    no_next_link = false;
          keyList newseed = seed;
          newseed.push_back(link.second);
          newtempSeedKeyLists.push_back(newseed);
        }
      }
      if (seed.size() > 5)
      {
        trackSeedKeyLists.push_back(seed);
#if defined(_CLUSTER_LOG_TUPOUT_)
        for (int i=0;i<6;++i) {
          auto& v = globalPositions.at(seed[i]);
          _tupclus_grown_seeds->Fill(_nevent, TrkrDefs::getLayer(seed[i]), i, v.x(), v.y(), v.z());
        }
#endif
      }
    }
    if (Verbosity() > 0)
    {
      std::cout << "new temp size: " << newtempSeedKeyLists.size() << std::endl;
    }
    tempSeedKeyLists = newtempSeedKeyLists;
  }

  //  trackSeedKeyLists = tempSeedKeyLists;
  /*
    for(auto trackKeyChain = trackSeedKeyLists.begin(); trackKeyChain != trackSeedKeyLists.end(); ++trackKeyChain)
    {
      bool reached_end = false;
      while(!reached_end)
      {
        TrkrDefs::cluskey trackHead = trackKeyChain->back();
        TrkrDefs::cluskey secondToLast = trackKeyChain->rbegin()[1];
        TrkrDefs::cluskey thirdToLast = trackKeyChain->rbegin()[2];
        auto& head_pos = globalPositions.at(trackHead);
        auto& sec_pos = globalPositions.at(secondToLast);
        auto& third_pos = globalPositions.at(thirdToLast);
        double dz_avg = ((head_pos.z()-sec_pos.z())+(sec_pos.z()-third_pos.z()))/2.;
        double dx1 = head_pos.x()-sec_pos.x();
        double dy1 = head_pos.y()-sec_pos.y();
        double dx2 = sec_pos.x()-third_pos.x();
        double dy2 = sec_pos.y()-third_pos.y();
        double ddx = dx1-dx2;
        double ddy = dy1-dy2;
        double new_dx = dx1+ddx;
        double new_dy = dy1+ddy;
        double new_x = head_pos.x()+new_dx;
        double new_y = head_pos.y()+new_dy;
        double new_z = head_pos.z()+dz_avg;
        std::cout << "(x,y,z) = (" << head_pos.x() << ", " << head_pos.y() << ", " << head_pos.z() << ")" << std::endl;
        unsigned int trackHead_layer = TrkrDefs::getLayer(trackHead) - (_nlayers_intt + _nlayers_maps);
        std::cout << "layer " << trackHead_layer << std::endl;
        std::cout << "projected: (" << new_x << ", " << new_y << ", " << new_z << ")" << std::endl;
        TrkrDefs::cluskey nextCluster;
        double bestDist = 1e9;
        bool no_next_link = true;
        for(auto testlink = bidirectionalLinks[trackHead_layer].begin(); testlink != bidirectionalLinks[trackHead_layer].end(); ++testlink)
        {
          if((*testlink)[0].second==trackHead)
          {
            TrkrDefs::cluskey testCluster = (*testlink)[1].second;
            auto& test_pos = globalPositions.at(testCluster);
            std::cout << "test cluster: (" << test_pos.x() << ", " << test_pos.y() << ", " << test_pos.z() << ")" << std::endl;
            double distToNew = sqrt(square<double>(test_pos.x()-new_x)+square<double>(test_pos.y()-new_y)+square<double>(test_pos.z()-new_z));
            if(distToNew<bestDist)
            {
              std::cout << "current best" << std::endl;
              nextCluster = testCluster;
              bestDist = distToNew;
            }
            no_next_link = false;
          }
        }
        if(!no_next_link) trackKeyChain->push_back(nextCluster);
        if(no_next_link) reached_end = true;
      }
    }
  */
  t_seed->stop();
  if (Verbosity() > 0)
  {
    std::cout << "keychain assembly time: " << t_seed->get_accumulated_time() / 1000 << " s" << std::endl;
  }
  t_seed->restart();
  LogDebug(" track key chains assembled: " << trackSeedKeyLists.size() << std::endl);
  LogDebug(" track key chain lengths: " << std::endl);
  for (auto trackKeyChain = trackSeedKeyLists.begin(); trackKeyChain != trackSeedKeyLists.end(); ++trackKeyChain)
  {
    LogDebug(" " << trackKeyChain->size() << std::endl);
  }
  int jumpcount = 0;
  LogDebug(" track key associations:" << std::endl);
  for (auto& trackSeedKeyList : trackSeedKeyLists)
  {
    LogDebug(" seed " << i << ":" << std::endl);

    double lasteta = -100;
    double lastphi = -100;
    for (unsigned long& j : trackSeedKeyList)
    {
      const auto& globalpos = globalPositions.at(j); 
      const double clus_phi = get_phi(globalpos);
      const double clus_eta = get_eta(globalpos);
      const double etajump = clus_eta - lasteta;
      const double phijump = clus_phi - lastphi;
#if defined(_DEBUG_)
      unsigned int lay = TrkrDefs::getLayer(trackSeedKeyLists[i][j].second);
#endif
      if ((fabs(etajump) > 0.1 && lasteta != -100) || (fabs(phijump) > 1 && lastphi != -100))
      {
        LogDebug(" Eta or Phi jump too large! " << std::endl);
        ++jumpcount;
      }
      LogDebug(" (eta,phi,layer) = (" << clus_eta << "," << clus_phi << "," << lay << ") "
                                      << " (x,y,z) = (" << globalpos(0) << "," << globalpos(1) << "," << globalpos(2) << ")" << std::endl);

      if (Verbosity() > 0)
      {
        unsigned int lay = TrkrDefs::getLayer(j);
        std::cout << "  eta, phi, layer = (" << clus_eta << "," << clus_phi << "," << lay << ") "
                  << " (x,y,z) = (" << globalpos(0) << "," << globalpos(1) << "," << globalpos(2) << ")" << std::endl;
      }
      lasteta = clus_eta;
      lastphi = clus_phi;
    }
  }
  LogDebug(" Total large jumps: " << jumpcount << std::endl);
  t_seed->stop();
  if (Verbosity() > 0)
  {
    std::cout << "eta-phi sanity check time: " << t_seed->get_accumulated_time() / 1000 << " s" << std::endl;
  }
  t_seed->restart();
  return trackSeedKeyLists;
}

std::vector<TrackSeed_v2> PHCASeeding::RemoveBadClusters(const std::vector<PHCASeeding::keyList>& chains, const PHCASeeding::PositionMap& globalPositions) const
{
  if (Verbosity() > 0)
  {
    std::cout << "removing bad clusters" << std::endl;
  }
  std::vector<TrackSeed_v2> clean_chains;

  for (const auto& chain : chains)
  {
    if (chain.size() < 3)
    {
      continue;
    }
    if (Verbosity() > 0)
    {
      std::cout << "chain size: " << chain.size() << std::endl;
    }

    TrackFitUtils::position_vector_t xy_pts;
    for (const auto& ckey : chain)
    {
      const auto& global = globalPositions.at(ckey);
      xy_pts.emplace_back(global.x(), global.y());
    }

    // fit a circle through x,y coordinates
    const auto [R, X0, Y0] = TrackFitUtils::circle_fit_by_taubin(xy_pts);

    // skip chain entirely if fit fails
    if (std::isnan(R))
    {
      continue;
    }

    // calculate residuals
    const std::vector<double> xy_resid = TrackFitUtils::getCircleClusterResiduals(xy_pts, R, X0, Y0);

    // assign clusters to seed
    TrackSeed_v2 trackseed;
    for (unsigned long i : chain)
    {
      // if(xy_resid[i]>_xy_outlier_threshold) continue;
      trackseed.insert_cluster_key(i);
    }

    clean_chains.push_back(trackseed);
    if (Verbosity() > 0)
    {
      std::cout << "pushed clean chain with " << trackseed.size_cluster_keys() << " clusters" << std::endl;
    }
  }

  return clean_chains;
}

void PHCASeeding::publishSeeds(const std::vector<TrackSeed_v2>& seeds)
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

int PHCASeeding::Setup(PHCompositeNode* topNode) // This is called by ::InitRun
{
  //  if(Verbosity()>0)
  std::cout << "Called Setup" << std::endl;
  if (Verbosity() > 0)
  {
    std::cout << "topNode:" << topNode << std::endl;
  }
  PHTrackSeeding::Setup(topNode);

  // geometry initialization
  int ret = InitializeGeometry(topNode);
  if (ret != Fun4AllReturnCodes::EVENT_OK)
  {
    return ret;
  }

  // tpc distortion correction
  m_dcc = findNode::getClass<TpcDistortionCorrectionContainer>(topNode, "TpcDistortionCorrectionContainerStatic");
  if (m_dcc)
  {
    std::cout << "PHCASeeding::Setup - found static TPC distortion correction container" << std::endl;
  }

  t_fill = std::make_unique<PHTimer>("t_fill");
  t_fill->stop();

  t_seed = std::make_unique<PHTimer>("t_seed");
  t_seed->stop();

  t_makebilinks = std::make_unique<PHTimer>("t_makebilinks");
  t_makebilinks->stop();

  t_makeseeds = std::make_unique<PHTimer>("t_makeseeds");
  t_makeseeds->stop();

  //  fcfg.set_rescale(1);
  std::unique_ptr<PHField> field_map;
  if (_use_const_field)
  {
    PHFieldConfigv2 fcfg(0, 0, _const_field);
    field_map = std::unique_ptr<PHField>(PHFieldUtility::BuildFieldMap(&fcfg));
  }
  else
  {
    PHFieldConfigv1 fcfg;
    fcfg.set_field_config(PHFieldConfig::FieldConfigTypes::Field3DCartesian);
    if (std::filesystem::path(m_magField).extension() != ".root")
    {
      m_magField = CDBInterface::instance()->getUrl(m_magField);
    }
    fcfg.set_filename(m_magField);
    field_map = std::unique_ptr<PHField>(PHFieldUtility::BuildFieldMap(&fcfg));
  }

  fitter = std::make_unique<ALICEKF>(topNode, _cluster_map, field_map.get(), _fieldDir, _min_clusters_per_track, _max_sin_phi, Verbosity());
  fitter->useConstBField(_use_const_field);
  fitter->setConstBField(_const_field);
  fitter->useFixedClusterError(_use_fixed_clus_err);
  fitter->setFixedClusterError(0, _fixed_clus_err.at(0));
  fitter->setFixedClusterError(1, _fixed_clus_err.at(1));
  fitter->setFixedClusterError(2, _fixed_clus_err.at(2));
  
#if defined(_CLUSTER_LOG_TUPOUT_)
  std::cout << " Writing _CLUSTER_LOG_TUPOUT.root file " << std::endl;
  _f_clustering_process = new TFile("_CLUSTER_LOG_TUPOUT.root", "recreate");
  _tupclus_all         = new TNtuple("all",         "all clusters","event:x:y:z");
  _tupclus_links       = new TNtuple("links",       "links","event:layer:updown01:x:y:z");
  _tupclus_bilinks     = new TNtuple("bilinks",     "bilinks","event:layer:topbot01:x:y:z");
  _tupclus_seeds       = new TNtuple("seeds",       "3 bilink seeds cores","event:layer:seed012:x:y:z");
  _tupclus_grown_seeds = new TNtuple("grown_seeds", "grown seeds", "event:layer:seednum05:x:y:z");
#endif


  return Fun4AllReturnCodes::EVENT_OK;
}

int PHCASeeding::End()
{
  if (Verbosity() > 0)
  {
    std::cout << "Called End " << std::endl;
  }

#if defined(_CLUSTER_LOG_TUPOUT_)
  _f_clustering_process->cd();
  _tupclus_all         ->Write();
  _tupclus_links       ->Write();
  _tupclus_bilinks     ->Write();
  _tupclus_seeds       ->Write();
  _tupclus_grown_seeds ->Write();
  _f_clustering_process->Close();
#endif

  return Fun4AllReturnCodes::EVENT_OK;
}
