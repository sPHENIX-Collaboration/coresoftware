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

// trackbase_historic includes
#include <trackbase/ActsSurfaceMaps.h>
#include <trackbase/ActsTrackingGeometry.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack_v2.h>

#include <trackbase/TrkrCluster.h>  // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrDefs.h>  // for getLayer, clu...
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrClusterIterationMapv1.h>

//ROOT includes for debugging
#include <TFile.h>
#include <TNtuple.h>

//BOOST for combi seeding
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/policies/compare.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <utility>  // for pair, make_pair
#include <vector>
#include <algorithm> // for find
#include <unordered_set>
#include <memory>

//#define _DEBUG_

#if defined(_DEBUG_)
#define LogDebug(exp) if(Verbosity()>0) std::cout << "DEBUG: " << __FILE__ << ": " << __LINE__ << ": " << exp
#else
#define LogDebug(exp) (void)0
#endif

#define LogError(exp) if(Verbosity()>0) std::cout << "ERROR: " << __FILE__ << ": " << __LINE__ << ": " << exp
#define LogWarning(exp) if(Verbosity()>0) std::cout << "WARNING: " << __FILE__ << ": " << __LINE__ << ": " << exp

//#define _DEBUG_

//end

typedef bg::model::point<float, 3, bg::cs::cartesian> point;
typedef bg::model::box<point> box;
typedef std::pair<point, TrkrDefs::cluskey> pointKey;
typedef std::pair<std::array<float,3>,TrkrDefs::cluskey> coordKey;
typedef std::array<coordKey,2> keylink;
typedef std::vector<TrkrDefs::cluskey> keylist;

// apparently there is no builtin STL hash function for a std::array
// so to use std::unordered_set (essentially a hash table), we have to make our own hasher

namespace std
{
  template<typename T,size_t N>
  struct hash<std::array<T,N>>
  {
    typedef std::array<T,N> argument_type;
    typedef size_t result_type;

    result_type operator()(const argument_type& a) const
    {
      hash<T> hasher;
      result_type h = 0;
      for(result_type i = 0; i < N; ++i)
      {
        h = h * 31 + hasher(a[i]);
      }
      return h;
    }
  };
  template<typename A,typename B>
  struct hash<pair<A,B>>
  {
    typedef pair<A,B> argument_type;
    typedef size_t result_type;
    
    result_type operator()(const argument_type& a) const
    {
      hash<A> hashA;
      hash<B> hashB;
      return (hashA(a.first)*31+hashB(a.second));
    }
  }; 
}

// anonymous namespace for local functions
namespace
{
  // square
  template<class T> inline constexpr T square( const T& x ) { return x*x; }
  
  /// phi angle of Acts::Vector3D
  inline double get_phi( const Acts::Vector3F& position )
  {
    double phi = std::atan2( position.y(), position.x() );
    if( phi < 0 ) phi += 2.*M_PI;
    return phi;
  } 
  
  /// pseudo rapidity of Acts::Vector3D
  inline double get_eta( const Acts::Vector3F& position )
  {
    const double norm = std::sqrt( square(position.x()) + square(position.y()) + square(position.z()) );
    return std::log((norm+position.z())/(norm-position.z()))/2;
  }
  
  ///@name utility convertion functions
  //@{

  coordKey fromPointKey(const pointKey& p)
  { return std::make_pair(std::array<float,3>({p.first.get<0>(),p.first.get<1>(),p.first.get<2>()}),p.second); }

  std::vector<coordKey> fromPointKey(const std::vector<pointKey>& p)
  {
    std::vector<coordKey> output;
    output.resize(p.size());
    std::transform( p.begin(), p.end(), std::back_inserter( output ), []( const pointKey& point )
      { return fromPointKey(point); } );
    return output;
  }
  //@}
  
  double breaking_angle(double x1, double y1, double z1, double x2, double y2, double z2)
  {
    double l1 = sqrt(x1*x1+y1*y1+z1*z1);
    double l2 = sqrt(x2*x2+y2*y2+z2*z2);
    double sx = (x1/l1+x2/l2);
    double sy = (y1/l1+y2/l2);
    double sz = (z1/l1+z2/l2);
    double dx = (x1/l1-x2/l2);
    double dy = (y1/l1-y2/l2);
    double dz = (z1/l1-z2/l2);
    return 2*atan2(sqrt(dx*dx+dy*dy+dz*dz),sqrt(sx*sx+sy*sy+sz*sz));
  }

}

//using namespace ROOT::Minuit2;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

PHCASeeding::PHCASeeding(
    const std::string &name,
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

int PHCASeeding::InitializeGeometry(PHCompositeNode *topNode)
{
  tGeometry = findNode::getClass<ActsTrackingGeometry>(topNode,"ActsTrackingGeometry");
  if(!tGeometry)
    {
      std::cout << PHWHERE << "No acts tracking geometry, can't proceed" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  
  surfMaps = findNode::getClass<ActsSurfaceMaps>(topNode,"ActsSurfaceMaps");
  if(!surfMaps)
    {
      std::cout << PHWHERE << "No acts surface maps, can't proceed" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    
  return Fun4AllReturnCodes::EVENT_OK;
}

Acts::Vector3D PHCASeeding::getGlobalPosition( TrkrCluster* cluster ) const
{
  // get global position from Acts transform
  auto globalpos = m_transform.getGlobalPosition(cluster,
    surfMaps,
    tGeometry);

  // check if TPC distortion correction are in place and apply
  if( m_dcc ) { globalpos = m_distortionCorrection.get_corrected_position( globalpos, m_dcc ); }

  return globalpos;
}

void PHCASeeding::QueryTree(const bgi::rtree<pointKey, bgi::quadratic<16>> &rtree, double phimin, double etamin, double lmin, double phimax, double etamax, double lmax, std::vector<pointKey> &returned_values) const
{
  double phimin_2pi = phimin;
  double phimax_2pi = phimax;
  if (phimin < 0) phimin_2pi = 2*M_PI+phimin;
  if (phimax > 2*M_PI) phimax_2pi = phimax-2*M_PI;
  rtree.query(bgi::intersects(box(point(phimin_2pi, etamin, lmin), point(phimax_2pi, etamax, lmax))), std::back_inserter(returned_values));
}

PositionMap PHCASeeding::FillTree()
{ 
  t_fill->stop();
  int n_dupli = 0;
  int nlayer[60];

  PositionMap cachedPositions;

  for (int j = 0; j < 60; ++j) nlayer[j] = 0;
  auto hitsetrange = _hitsets->getHitSets(TrkrDefs::TrkrId::tpcId);
  for (auto hitsetitr = hitsetrange.first;
       hitsetitr != hitsetrange.second;
       ++hitsetitr){
    auto range = _cluster_map->getClusters(hitsetitr->first);
    for( auto clusIter = range.first; clusIter != range.second; ++clusIter ){

      TrkrCluster *cluster = clusIter->second;
      TrkrDefs::cluskey ckey = clusIter->first;
      unsigned int layer = TrkrDefs::getLayer(ckey);
      if (layer < _start_layer || layer >= _end_layer){
	if(Verbosity()>0) std::cout << "layer: " << layer << std::endl;
	continue;
      }
      if(_iteration_map != NULL && _n_iteration >0){
	if( _iteration_map->getIteration(ckey) > 0) 
	  continue; // skip hits used in a previous iteration
      }

      // get global position, convert to Acts::Vector3F and store in map
      const Acts::Vector3D globalpos_d = getGlobalPosition(cluster);
      const Acts::Vector3F globalpos = { (float) globalpos_d.x(), (float) globalpos_d.y(), (float) globalpos_d.z()};
      cachedPositions.insert(std::make_pair(ckey, globalpos));

      const double clus_phi = get_phi( globalpos );      
      const double clus_eta = get_eta( globalpos );
      const double clus_l = layer;  

      if(Verbosity() > 0) 
	std::cout << "Found cluster " << ckey << " in layer " << layer << std::endl;
      
      std::vector<pointKey> testduplicate;
      QueryTree(_rtree, clus_phi - 0.00001, clus_eta - 0.00001, layer - 0.5, clus_phi + 0.00001, clus_eta + 0.00001, layer + 0.5, testduplicate);
      if (!testduplicate.empty())
	{
	  ++n_dupli;
	  continue;
	}
      ++nlayer[layer];
      t_fill->restart();
      _rtree.insert(std::make_pair(point(clus_phi, clus_eta, clus_l), ckey));
      t_fill->stop();
    }
  }
  if(Verbosity()>1) for (int j = 0; j < 60; ++j) std::cout << "nhits in layer " << j << ":  " << nlayer[j] << std::endl;
  if(Verbosity()>0) std::cout << "fill time: " << t_fill->get_accumulated_time() / 1000. << " sec" << std::endl;
  if(Verbosity()>0) std::cout << "number of duplicates : " << n_dupli << std::endl;
  return cachedPositions;
}

int PHCASeeding::Process(PHCompositeNode */*topNode*/)
{
//  TFile fpara("CA_para.root", "RECREATE");
  if(_n_iteration>0){
    if (!_iteration_map){
      std::cerr << PHWHERE << "Cluster Iteration Map missing, aborting." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  t_seed->restart();

  _rtree.clear();
  PositionMap globalClusPositions = FillTree();
  t_seed->stop();
  if(Verbosity()>0) std::cout << "Initial RTree fill time: " << t_seed->get_accumulated_time() / 1000 << " s" << std::endl;
  t_seed->restart();
  int numberofseeds = 0;
  numberofseeds += FindSeedsWithMerger(globalClusPositions);
  t_seed->stop();
  if(Verbosity()>0) std::cout << "number of seeds " << numberofseeds << std::endl;
  if(Verbosity()>0) std::cout << "Kalman filtering time: " << t_seed->get_accumulated_time() / 1000 << " s" << std::endl;
//  fpara.cd();
//  fpara.Close();
//  if(Verbosity()>0) std::cout << "fpara OK\n";
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHCASeeding::FindSeedsWithMerger(const PositionMap& globalPositions)
{
  std::vector<pointKey> allClusters;
  std::vector<std::unordered_set<keylink>> belowLinks;
  std::vector<std::unordered_set<keylink>> aboveLinks;
  belowLinks.resize(_nlayers_tpc);
  aboveLinks.resize(_nlayers_tpc);
  QueryTree(_rtree,
            0, // phi
            -3, // eta
            _start_layer-0.5, // layer 
            2*M_PI, // phi
            3, // eta
            _end_layer+0.5, // layer
            allClusters);
  t_seed->stop();
  if(Verbosity()>0) std::cout << "allClusters search time: " << t_seed->get_accumulated_time() / 1000 << " s" << std::endl;
  LogDebug(" number of clusters: " << allClusters.size() << std::endl);
  t_seed->restart();

  std::pair<std::vector<std::unordered_set<keylink>>,std::vector<std::unordered_set<keylink>>> links = CreateLinks(fromPointKey(allClusters), globalPositions);
  std::vector<std::vector<keylink>> biLinks = FindBiLinks(links.first,links.second);
  std::vector<keylist> trackSeedKeyLists = FollowBiLinks(biLinks,globalPositions);
  std::vector<keylist> cleanSeedKeyLists = RemoveBadClusters(trackSeedKeyLists, globalPositions);
  std::vector<SvtxTrack_v2> seeds = fitter->ALICEKalmanFilter(cleanSeedKeyLists,true, globalPositions);
  publishSeeds(seeds);
  return seeds.size();
}

std::pair<std::vector<std::unordered_set<keylink>>,std::vector<std::unordered_set<keylink>>> PHCASeeding::CreateLinks(const std::vector<coordKey>& clusters, const PositionMap& globalPositions, int mode) const
{
  size_t nclusters = 0;

  double cluster_find_time = 0;
  double rtree_query_time = 0;
  double transform_time = 0;
  double compute_best_angle_time = 0;
  double set_insert_time = 0;

  std::vector<std::unordered_set<keylink>> belowLinks;
  std::vector<std::unordered_set<keylink>> aboveLinks;
  belowLinks.resize(_nlayers_tpc);
  aboveLinks.resize(_nlayers_tpc);

  for (auto StartCluster = clusters.begin(); StartCluster != clusters.end(); ++StartCluster)
  {
    nclusters++;
    // get clusters near this one in adjacent layers
    double StartPhi = StartCluster->first[0];
    double StartEta = StartCluster->first[1];
    unsigned int StartLayer = StartCluster->first[2];
    if(StartLayer < _start_layer) continue;
    if(StartLayer > _end_layer) continue;
    const auto& globalpos = globalPositions.at(StartCluster->second);
    double StartX = globalpos(0);
    double StartY = globalpos(1);
    double StartZ = globalpos(2);
    t_seed->stop();
    cluster_find_time += t_seed->elapsed();
    t_seed->restart();
    LogDebug(" starting cluster:" << std::endl);
    LogDebug(" eta: " << StartEta << std::endl);
    LogDebug(" phi: " << StartPhi << std::endl);
    LogDebug(" layer: " << StartLayer << std::endl);

    std::vector<pointKey> ClustersAbove;
    std::vector<pointKey> ClustersBelow;
    QueryTree(_rtree,
              StartPhi-_neighbor_phi_width,
              StartEta-_neighbor_eta_width,
              (double) StartLayer - 1.5,
              StartPhi+_neighbor_phi_width,
              StartEta+_neighbor_eta_width,
              (double) StartLayer - 0.5,
              ClustersBelow);
    QueryTree(_rtree,
              StartPhi-_neighbor_phi_width,
              StartEta-_neighbor_eta_width,
              (double) StartLayer + 0.5,
              StartPhi+_neighbor_phi_width,
              StartEta+_neighbor_eta_width,
              (double) StartLayer + 1.5,
              ClustersAbove);
    t_seed->stop();
    rtree_query_time += t_seed->elapsed();
    t_seed->restart();
    LogDebug(" entries in below layer: " << ClustersBelow.size() << std::endl);
    LogDebug(" entries in above layer: " << ClustersAbove.size() << std::endl);
    std::vector<std::array<double,3>> delta_below;
    std::vector<std::array<double,3>> delta_above;
    delta_below.clear();
    delta_above.clear();
    delta_below.resize(ClustersBelow.size());
    delta_above.resize(ClustersAbove.size());
    // calculate (delta_eta, delta_phi) vector for each neighboring cluster

    std::transform(ClustersBelow.begin(),ClustersBelow.end(),delta_below.begin(),
	      [&](pointKey BelowCandidate){
	const auto& belowpos = globalPositions.at(BelowCandidate.second);
        return std::array<double,3>{belowpos(0)-StartX,
          belowpos(1)-StartY,
          belowpos(2)-StartZ};});

    std::transform(ClustersAbove.begin(),ClustersAbove.end(),delta_above.begin(),
      [&](pointKey AboveCandidate){
	const auto& abovepos = globalPositions.at(AboveCandidate.second);
        return std::array<double,3>{abovepos(0)-StartX,
          abovepos(1)-StartY,
          abovepos(2)-StartZ};});
    t_seed->stop();
    transform_time += t_seed->elapsed();
    t_seed->restart();

    // find the three clusters closest to a straight line
    // (by maximizing the cos of the angle between the (delta_eta,delta_phi) vectors)
    double maxCosPlaneAngle = -0.9;
    //double minSumLengths = 1e9;
    coordKey bestBelowCluster = std::make_pair(std::array<float,3>({0.,0.,-1e9}),0);
    coordKey bestAboveCluster = std::make_pair(std::array<float,3>({0.,0.,-1e9}),0);
    for(size_t iAbove = 0; iAbove<delta_above.size(); ++iAbove)
    {
      for(size_t iBelow = 0; iBelow<delta_below.size(); ++iBelow)
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
        if(cos(angle) < maxCosPlaneAngle)
        {
          maxCosPlaneAngle = cos(angle);
          //minSumLengths = belowLength+aboveLength;
          bestBelowCluster = fromPointKey(ClustersBelow[iBelow]);
          bestAboveCluster = fromPointKey(ClustersAbove[iAbove]);
        }
      }
    }
    
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
    t_seed->stop();
    compute_best_angle_time += t_seed->elapsed();
    t_seed->restart();
    int layer_index = StartLayer - (_nlayers_intt + _nlayers_maps);
    if(bestBelowCluster.second != 0) belowLinks[layer_index].insert(keylink{{*StartCluster,bestBelowCluster}});
    if(bestAboveCluster.second != 0) aboveLinks[layer_index].insert(keylink{{*StartCluster,bestAboveCluster}});
    t_seed->stop();
    set_insert_time += t_seed->elapsed();
    t_seed->restart();
    LogDebug(" max collinearity: " << maxCosPlaneAngle << std::endl);
  }
  t_seed->stop();
  if(Verbosity()>0)
  {
    std::cout << "triplet forming time: " << t_seed->get_accumulated_time() / 1000 << " s" << std::endl;
    std::cout << "starting cluster setup: " << cluster_find_time / 1000 << " s" << std::endl;
    std::cout << "RTree query: " << rtree_query_time /1000 << " s" << std::endl;
    std::cout << "Transform: " << transform_time /1000 << " s" << std::endl;
    std::cout << "Compute best triplet: " << compute_best_angle_time /1000 << " s" << std::endl;
    std::cout << "Set insert: " << set_insert_time /1000 << " s" << std::endl;
  }
  t_seed->restart();

  return std::make_pair(belowLinks,aboveLinks);
}

std::vector<std::vector<keylink>> PHCASeeding::FindBiLinks(const std::vector<std::unordered_set<keylink>>& belowLinks, const std::vector<std::unordered_set<keylink>>& aboveLinks) const
{
  // remove all triplets for which there isn't a mutual association between two clusters
  std::vector<std::vector<keylink>> bidirectionalLinks;
  bidirectionalLinks.resize(_nlayers_tpc);
  for(int layer = _nlayers_tpc-1; layer > 0; --layer)
  {
    for(auto belowLink = belowLinks[layer].begin(); belowLink != belowLinks[layer].end(); ++belowLink)
    {
      if((*belowLink)[1].second==0) continue;
      unsigned int end_layer_index = TrkrDefs::getLayer((*belowLink)[1].second) - (_nlayers_intt + _nlayers_maps);
      keylink reversed = {(*belowLink)[1],(*belowLink)[0]};
      auto sameAboveLinkExists = aboveLinks[end_layer_index].find(reversed);
      if(sameAboveLinkExists != aboveLinks[end_layer_index].end())
      {
        bidirectionalLinks[layer].push_back((*belowLink));
      }
    }
  }
  t_seed->stop();
  if(Verbosity()>0) std::cout << "bidirectional link forming time: " << t_seed->get_accumulated_time() / 1000 << " s" << std::endl;
  t_seed->restart();

  return bidirectionalLinks;
}

std::vector<keylist> PHCASeeding::FollowBiLinks(const std::vector<std::vector<keylink>>& bidirectionalLinks, const PositionMap& globalPositions) const
{
  // follow bidirectional links to form lists of cluster keys
  // (to be fitted for track seed parameters)
  std::vector<keylist> trackSeedKeyLists;
  // get starting cluster keys, create a keylist for each
  // (only check last element of each pair because we start from the outer layers and go inward)
  for(unsigned int layer = 0; layer < _nlayers_tpc-1; ++layer)
  {
    for(auto startCand = bidirectionalLinks[layer].begin(); startCand != bidirectionalLinks[layer].end(); ++startCand)
    {
      bool has_above_link = false;
      unsigned int imax = 1;
      if(layer==_nlayers_tpc-2) imax = 1;
      for(unsigned int i=1;i<=imax;i++)
      {
        has_above_link = has_above_link || std::any_of(bidirectionalLinks[layer+i].begin(),bidirectionalLinks[layer+i].end(),[&](keylink k){return (*startCand)[0]==k[1];});
      }
//      for(std::vector<keylink>::iterator testlink = bidirectionalLinks[layer+1].begin(); !has_above_link && (testlink != bidirectionalLinks[layer+1].end()); ++testlink)
//      {
//        if((*startCand) == (*testlink)) continue;
//        if((*startCand)[0] == (*testlink)[1]) has_above_link = true;
//      } 
      if(!has_above_link)
      {
        trackSeedKeyLists.push_back({(*startCand)[0].second,(*startCand)[1].second});
      }
    }
  }
  t_seed->stop();
  if(Verbosity()>0) std::cout << "starting cluster finding time: " << t_seed->get_accumulated_time() / 1000 << " s" << std::endl;
  t_seed->restart();
  // assemble track cluster chains from starting cluster keys (ordered from outside in)
  for(auto trackKeyChain = trackSeedKeyLists.begin(); trackKeyChain != trackSeedKeyLists.end(); ++trackKeyChain)
  {
    bool reached_end = false;
    while(!reached_end)
    {
      TrkrDefs::cluskey trackHead = trackKeyChain->back();
      unsigned int trackHead_layer = TrkrDefs::getLayer(trackHead) - (_nlayers_intt + _nlayers_maps);
      bool no_next_link = true;
      for(auto testlink = bidirectionalLinks[trackHead_layer].begin(); testlink != bidirectionalLinks[trackHead_layer].end(); ++testlink)
      {
        if((*testlink)[0].second==trackHead)
        {
          trackKeyChain->push_back((*testlink)[1].second);
          no_next_link = false;
        }
      }
      if(no_next_link) reached_end = true;
    }
  }
  t_seed->stop();
  if(Verbosity()>0) std::cout << "keychain assembly time: " << t_seed->get_accumulated_time() / 1000 << " s" << std::endl;
  t_seed->restart();
  LogDebug(" track key chains assembled: " << trackSeedKeyLists.size() << std::endl);
  LogDebug(" track key chain lengths: " << std::endl);
  for(auto trackKeyChain = trackSeedKeyLists.begin(); trackKeyChain != trackSeedKeyLists.end(); ++trackKeyChain)
  {
    LogDebug(" " << trackKeyChain->size() << std::endl);
  }
  int jumpcount = 0;
  LogDebug(" track key associations:" << std::endl);
  for(size_t i=0;i<trackSeedKeyLists.size();++i)
  {
    LogDebug(" seed " << i << ":" << std::endl);

    double lasteta = -100;
    double lastphi = -100;
    for(size_t j=0;j<trackSeedKeyLists[i].size();++j)
    {
      const auto& globalpos = globalPositions.at(trackSeedKeyLists[i][j]);           
      const double clus_phi = get_phi( globalpos );
      const double clus_eta = get_eta( globalpos );
      const double etajump = clus_eta-lasteta;
      const double phijump = clus_phi-lastphi;
      #if defined(_DEBUG_) 
      unsigned int lay = TrkrDefs::getLayer(trackSeedKeyLists[i][j].second);
      #endif
      if((fabs(etajump)>0.1 && lasteta!=-100) || (fabs(phijump)>1 && lastphi!=-100))
	{
           LogDebug(" Eta or Phi jump too large! " << std::endl);
           ++jumpcount;
        }
      LogDebug(" (eta,phi,layer) = (" << clus_eta << "," << clus_phi << "," << lay << ") " <<
        " (x,y,z) = (" << globalpos(0) << "," << globalpos(1) << "," << globalpos(2) << ")" << std::endl);
      
      if(Verbosity() > 0)
	{
            unsigned int lay = TrkrDefs::getLayer(trackSeedKeyLists[i][j]);
            std::cout << "  eta, phi, layer = (" << clus_eta << "," << clus_phi << "," << lay << ") " <<
             " (x,y,z) = (" << globalpos(0) << "," << globalpos(1) << "," << globalpos(2) << ")" << std::endl;
        }
      lasteta = clus_eta;
      lastphi = clus_phi;
    }
  }
  LogDebug(" Total large jumps: " << jumpcount << std::endl);
  t_seed->stop();
  if(Verbosity()>0) std::cout << "eta-phi sanity check time: " << t_seed->get_accumulated_time() / 1000 << " s" << std::endl;
  t_seed->restart();
  return trackSeedKeyLists;
}

std::vector<keylist> PHCASeeding::RemoveBadClusters(const std::vector<keylist>& chains, const PositionMap& globalPositions) const
{
  if(Verbosity()>0) std::cout << "removing bad clusters" << std::endl;
  std::vector<keylist> clean_chains;

  for(const auto& chain : chains)
  {
    if(chain.size()<3) continue;
    keylist clean_chain;

    std::vector<std::pair<double,double>> xy_pts;
//     std::vector<std::pair<double,double>> rz_pts;

    for(const TrkrDefs::cluskey& ckey : chain)
    {
      const auto &global = globalPositions.at(ckey);
      double x = global(0);
      double y = global(1);
      xy_pts.push_back(std::make_pair(x,y));
//       double z = global(2);
//       rz_pts.push_back(std::make_pair(std::sqrt(square(x)+square(y)),z));
    }
    if(Verbosity()>0) std::cout << "chain size: " << chain.size() << std::endl;

//     double A = 0;
//     double B = 0;
//     fitter->line_fit(rz_pts,A,B);
//     const std::vector<double> rz_resid = fitter->GetLineClusterResiduals(rz_pts,A,B);
    double R = 0;
    double X0 = 0;
    double Y0 = 0;
    fitter->CircleFitByTaubin(xy_pts,R,X0,Y0);

    // skip chain entirely if fit fails
    /*
     * note: this is consistent with what the code was doing before
     * but in principle we could also keep the seed unchanged instead
     * this is to be studied independently
     */
    if( std::isnan( R ) ) continue;

    // calculate residuals
    const std::vector<double> xy_resid = fitter->GetCircleClusterResiduals(xy_pts,R,X0,Y0);
    for(size_t i=0;i<chain.size();i++)
    {
      if(xy_resid[i]>_xy_outlier_threshold) continue;
      clean_chain.push_back(chain[i]);
    }

    clean_chains.push_back(clean_chain);
    if(Verbosity()>0) std::cout << "pushed clean chain with " << clean_chain.size() << " clusters" << std::endl;
  }
  return clean_chains;
}


void PHCASeeding::publishSeeds(const std::vector<SvtxTrack_v2>& seeds)
{
  for( const auto&  seed:seeds )
  { _track_map->insert(&seed);}
}

int PHCASeeding::Setup(PHCompositeNode *topNode)
{
  if(Verbosity()>0) std::cout << "Called Setup" << std::endl;
  if(Verbosity()>0) std::cout << "topNode:" << topNode << std::endl;
  PHTrackSeeding::Setup(topNode);
  
  // geometry initialization
  int ret = InitializeGeometry(topNode);
  if(ret != Fun4AllReturnCodes::EVENT_OK)
    { return ret; }
    
  // tpc distortion correction
  m_dcc = findNode::getClass<TpcDistortionCorrectionContainer>(topNode,"TpcDistortionCorrectionContainer");
  if( m_dcc )
  { std::cout << "PHCASeeding::Setup - found TPC distortion correction container" << std::endl; }
  
  t_fill = std::make_unique<PHTimer>("t_fill");
  t_seed = std::make_unique<PHTimer>("t_seed");
  t_fill->stop();
  t_seed->stop();
  fitter = std::make_unique<ALICEKF>(topNode,_cluster_map,_fieldDir,_min_clusters_per_track,_max_sin_phi,Verbosity());
  fitter->useConstBField(_use_const_field);
  fitter->useFixedClusterError(_use_fixed_clus_err);
  fitter->setFixedClusterError(0,_fixed_clus_err.at(0));
  fitter->setFixedClusterError(1,_fixed_clus_err.at(1));
  fitter->setFixedClusterError(2,_fixed_clus_err.at(2));
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHCASeeding::End()
{
  if(Verbosity()>0) std::cout << "Called End " << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
