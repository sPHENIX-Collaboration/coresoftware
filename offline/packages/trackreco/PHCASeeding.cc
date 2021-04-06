/*!
 *  \file PHCASeeding.cc
 *  \brief Track seeding using ALICE-style "cellular automaton" (CA) algorithm
 *  \detail 
 *  \author Michael Peters & Christof Roland
 */

#include "PHCASeeding.h"
#include "GPUTPCTrackLinearisation.h"
#include "GPUTPCTrackParam.h"

// trackbase_historic includes
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrack_v2.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <trackbase/TrkrCluster.h>  // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>  // for getLayer, clu...
#include <trackbase/TrkrClusterHitAssocv2.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

// sPHENIX Geant4 includes
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>

// sPHENIX includes
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHTimer.h>  // for PHTimer
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

//ROOT includes for debugging
#include <TFile.h>
#include <TNtuple.h>
#include <TVector3.h>  // for TVector3

//BOOST for combi seeding
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/policies/compare.hpp>

#if __cplusplus < 201402L
#include <boost/make_unique.hpp>
#endif

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

// forward declarations
class PHCompositeNode;



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
  struct hash<array<T,N>>
  {
    typedef array<T,N> argument_type;
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

using namespace std;
//using namespace ROOT::Minuit2;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

PHCASeeding::PHCASeeding(
    const string &name,
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
    float Bz,
    float cosTheta_limit)
  : PHTrackSeeding(name)
  , _g4tracks(nullptr)
  , _g4vertexes(nullptr)
  , _svtxhitsmap(nullptr)
  , _hit_used_map(nullptr)
  , _hit_used_map_size(0)
  , _vertex(nullptr)
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
  , _Bz(Bz)
  , _cosTheta_limit(cosTheta_limit)
{
}

int PHCASeeding::InitializeGeometry(PHCompositeNode *topNode)
{
  PHG4CylinderCellGeomContainer *cellgeos = findNode::getClass<
      PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  PHG4CylinderGeomContainer *laddergeos = findNode::getClass<
      PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_INTT");
  PHG4CylinderGeomContainer *mapsladdergeos = findNode::getClass<
      PHG4CylinderGeomContainer>(topNode, "CYLINDERGEOM_MVTX");

  //_nlayers_seeding = _seeding_layer.size();
  //_radii.assign(_nlayers_seeding, 0.0);
  map<float, int> radius_layer_map;

  _radii_all.assign(60, 0.0);
  _layer_ilayer_map.clear();
  _layer_ilayer_map_all.clear();
  if (cellgeos)
  {
    PHG4CylinderCellGeomContainer::ConstRange layerrange =
        cellgeos->get_begin_end();
    for (PHG4CylinderCellGeomContainer::ConstIterator layeriter =
             layerrange.first;
         layeriter != layerrange.second; ++layeriter)
    {
      radius_layer_map.insert(
          make_pair(layeriter->second->get_radius(),
                    layeriter->second->get_layer()));
    }
  }

  if (laddergeos)
  {
    PHG4CylinderGeomContainer::ConstRange layerrange =
        laddergeos->get_begin_end();
    for (PHG4CylinderGeomContainer::ConstIterator layeriter =
             layerrange.first;
         layeriter != layerrange.second; ++layeriter)
    {
      radius_layer_map.insert(
          make_pair(layeriter->second->get_radius(),
                    layeriter->second->get_layer()));
    }
  }

  if (mapsladdergeos)
  {
    PHG4CylinderGeomContainer::ConstRange layerrange =
        mapsladdergeos->get_begin_end();
    for (PHG4CylinderGeomContainer::ConstIterator layeriter =
             layerrange.first;
         layeriter != layerrange.second; ++layeriter)
    {
      radius_layer_map.insert(
          make_pair(layeriter->second->get_radius(),
                    layeriter->second->get_layer()));
    }
  }
  for (map<float, int>::iterator iter = radius_layer_map.begin();
       iter != radius_layer_map.end(); ++iter)
  {
    _layer_ilayer_map_all.insert(make_pair(iter->second, _layer_ilayer_map_all.size()));

    /*if (std::find(_seeding_layer.begin(), _seeding_layer.end(),
                  iter->second) != _seeding_layer.end())
    {
      _layer_ilayer_map.insert(make_pair(iter->second, ilayer));
      ++ilayer;
      }*/
  }
  if (cellgeos)
  {
    PHG4CylinderCellGeomContainer::ConstRange begin_end =
        cellgeos->get_begin_end();
    PHG4CylinderCellGeomContainer::ConstIterator miter = begin_end.first;
    for (; miter != begin_end.second; ++miter)
    {
      PHG4CylinderCellGeom *geo = miter->second;
      _radii_all[_layer_ilayer_map_all[geo->get_layer()]] =
          geo->get_radius() + 0.5 * geo->get_thickness();

      /*if (_layer_ilayer_map.find(geo->get_layer()) != _layer_ilayer_map.end())
      {
        _radii[_layer_ilayer_map[geo->get_layer()]] =
            geo->get_radius();
	    }*/
    }
  }

  if (laddergeos)
  {
    PHG4CylinderGeomContainer::ConstRange begin_end =
        laddergeos->get_begin_end();
    PHG4CylinderGeomContainer::ConstIterator miter = begin_end.first;
    for (; miter != begin_end.second; ++miter)
    {
      PHG4CylinderGeom *geo = miter->second;
      _radii_all[_layer_ilayer_map_all[geo->get_layer()]] =
          geo->get_radius() + 0.5 * geo->get_thickness();

      /*if (_layer_ilayer_map.find(geo->get_layer()) != _layer_ilayer_map.end())
      {
        _radii[_layer_ilayer_map[geo->get_layer()]] = geo->get_radius();
	}*/
    }
  }

  if (mapsladdergeos)
  {
    PHG4CylinderGeomContainer::ConstRange begin_end =
        mapsladdergeos->get_begin_end();
    PHG4CylinderGeomContainer::ConstIterator miter = begin_end.first;
    for (; miter != begin_end.second; ++miter)
    {
      PHG4CylinderGeom *geo = miter->second;

      //if(geo->get_layer() > (int) _radii.size() ) continue;

      //			if (Verbosity() >= 2)
      //				geo->identify();

      //TODO
      _radii_all[_layer_ilayer_map_all[geo->get_layer()]] =
          geo->get_radius();

      /*if (_layer_ilayer_map.find(geo->get_layer()) != _layer_ilayer_map.end())
      {
        _radii[_layer_ilayer_map[geo->get_layer()]] = geo->get_radius();
	}*/
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

double PHCASeeding::phiadd(double phi1, double phi2)
{
  double s = phi1 + phi2;
  if (s > 2 * M_PI)
    return s - 2 * M_PI;
  else if (s < 0)
    return s + 2 * M_PI;
  else
    return s;
}

double PHCASeeding::phidiff(double phi1, double phi2)
{
  double d = phi1 - phi2;
  if (d > M_PI)
    return d - 2 * M_PI;
  else if (d < -M_PI)
    return d + 2 * M_PI;
  else
    return d;
}

void PHCASeeding::QueryTree(const bgi::rtree<pointKey, bgi::quadratic<16>> &rtree, double phimin, double etamin, double lmin, double phimax, double etamax, double lmax, std::vector<pointKey> &returned_values)
{
  double phimin_2pi = phimin;
  double phimax_2pi = phimax;
  if (phimin < 0) phimin_2pi = 2*M_PI+phimin;
  if (phimax > 2*M_PI) phimax_2pi = phimax-2*M_PI;
  rtree.query(bgi::intersects(box(point(phimin_2pi, etamin, lmin), point(phimax_2pi, etamax, lmax))), std::back_inserter(returned_values));
}

void PHCASeeding::FillTree()
{ 
  t_fill->stop();
  int n_dupli = 0;
  int nlayer[60];
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
	cout << "layer: " << layer << endl;
	continue;
      }/*
      if(!_use_truth_clusters)
	{
	  std::pair<std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator, std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator>   
	  hitrange = _cluster_hit_map->getHits(ckey);
	  unsigned int nhits = std::distance(hitrange.first,hitrange.second);
	  if(nhits<_min_nhits_per_cluster) continue;
	}
      */
      TVector3 vec(cluster->getPosition(0)-_vertex->get_x(), cluster->getPosition(1)-_vertex->get_y(), cluster->getPosition(2)-_vertex->get_z());
      double clus_phi = vec.Phi();
      if(clus_phi<0) clus_phi = 2*M_PI + clus_phi;
      //clus_phi -= 2 * M_PI * floor(clus_phi / (2 * M_PI));
      double clus_eta = vec.Eta();
      double clus_l = layer;  // _radii_all[layer];
      if(Verbosity() > 0) 
	std::cout << "Found cluster " << ckey << " in layer " << layer << std::endl;
      
      vector<pointKey> testduplicate;
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
  if(Verbosity()>1)for (int j = 0; j < 60; ++j) cout << "nhits in layer " << j << ":  " << nlayer[j] << endl;
  if(Verbosity()>0) std::cout << "fill time: " << t_fill->get_accumulated_time() / 1000. << " sec" << std::endl;
  if(Verbosity()>0) std::cout << "number of duplicates : " << n_dupli << std::endl;
}

pointKey PHCASeeding::makepointKey(TrkrDefs::cluskey k)
{
  TrkrCluster* cluster = _cluster_map->findCluster(k);
  TVector3 vec(cluster->getPosition(0)-_vertex->get_x(), cluster->getPosition(1)-_vertex->get_y(), cluster->getPosition(2)-_vertex->get_z());
  double clus_phi = vec.Phi();
  if(clus_phi<0) clus_phi = 2*M_PI + clus_phi;
  if(clus_phi>2*M_PI) clus_phi = clus_phi - 2*M_PI;
  double clus_eta = vec.Eta();
  double clus_l = TrkrDefs::getLayer(k);
  return std::make_pair(point(clus_phi,clus_eta,clus_l),k);
}

void PHCASeeding::FillTree(vector<pointKey> clusters)
{
  // WARNING: THIS VERSION DOES NOT TEST FOR DUPLICATES!!!
  // It's much faster that way (no more querying the RTree as it's being constructed) 
  // and we can get away with it because
  // the incoming vector is guaranteed by construction to be duplicate-free.
  t_fill->stop();
  int nlayer[60];
  for (int j = 0; j < 60; ++j) nlayer[j] = 0;

  for (vector<pointKey>::iterator iter = clusters.begin(); iter != clusters.end(); ++iter)
  {
    unsigned int layer = TrkrDefs::getLayer(iter->second);
    if(layer < _start_layer || layer >= _end_layer) continue;

    std::pair<std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator, std::multimap<TrkrDefs::cluskey, TrkrDefs::hitkey>::const_iterator> 
      hitrange = _cluster_hit_map->getHits(iter->second);

    unsigned int nhits = std::distance(hitrange.first,hitrange.second);
    if(nhits<_min_nhits_per_cluster){
      cout << "min hits fail" << endl;
      continue;
    }
    ++nlayer[layer];
    t_fill->restart();
    _rtree.insert(*iter);
    t_fill->stop();
  }

  if(Verbosity()>0) std::cout << "fill time: " << t_fill->get_accumulated_time() / 1000. << " sec" << std::endl;
}

pointKey PHCASeeding::toPointKey(coordKey v)
{
  return make_pair(point(v.first.at(0),v.first.at(1),v.first.at(2)),v.second);
}

vector<pointKey> PHCASeeding::toPointKey(vector<coordKey> v)
{
  vector<pointKey> output;
  output.resize(v.size());
  for(vector<coordKey>::iterator ck = v.begin();ck != v.end();++ck)
  {
    output.push_back(make_pair(point(ck->first.at(0),ck->first.at(1),ck->first.at(2)),ck->second));
  }
  return output;
}

coordKey PHCASeeding::fromPointKey(pointKey p)
{
  return make_pair(array<float,3>({p.first.get<0>(),p.first.get<1>(),p.first.get<2>()}),p.second);
}

vector<coordKey> PHCASeeding::fromPointKey(vector<pointKey> p)
{
  vector<coordKey> output;
  output.resize(p.size());
  for(vector<pointKey>::iterator pk = p.begin();pk != p.end();++pk)
  {
    output.push_back(make_pair(array<float,3>({pk->first.get<0>(),pk->first.get<1>(),pk->first.get<2>()}),pk->second));
  }
  return output;
}

int PHCASeeding::Process(PHCompositeNode *topNode)
{
//  TFile fpara("CA_para.root", "RECREATE");
  _vertex = _vertex_map->get(0);

  t_seed->restart();

  _rtree.clear();
  FillTree();
  t_seed->stop();
  if(Verbosity()>0) cout << "Initial RTree fill time: " << t_seed->get_accumulated_time() / 1000 << " s" << endl;
  t_seed->restart();
  int numberofseeds = 0;
  numberofseeds += FindSeedsWithMerger();
  t_seed->stop();
  //  if(Verbosity()>0)   
if(Verbosity()>1)  cout << "number of seeds " << numberofseeds << endl;
  if(Verbosity()>0) cout << "Kalman filtering time: " << t_seed->get_accumulated_time() / 1000 << " s" << endl;
//  fpara.cd();
//  fpara.Close();
//  if(Verbosity()>0) cout << "fpara OK\n";
  return Fun4AllReturnCodes::EVENT_OK;
}

vector<coordKey> PHCASeeding::FindLinkedClusters()
{
  vector<pointKey> allClusters;
  vector<unordered_set<keylink>> belowLinks;
  vector<unordered_set<keylink>> aboveLinks;
  belowLinks.resize(_nlayers_tpc);
  aboveLinks.resize(_nlayers_tpc);
  // get vector<pointKey> for all clusters in outer third of TPC
  QueryTree(_rtree,
            0, // phi
            -3, // eta
            _nlayers_maps+_nlayers_intt-0.5, // layer 
            2*M_PI, // phi
            3, // eta
            _nlayers_maps+_nlayers_intt+_nlayers_tpc+0.5, // layer
            allClusters);
  t_seed->stop();
  if(Verbosity()>0) cout << "allClusters search time: " << t_seed->get_accumulated_time() / 1000 << " s" << endl;
  LogDebug(" number of total clusters: " << allClusters.size() << endl);
  t_seed->restart();

  pair<vector<unordered_set<keylink>>,vector<unordered_set<keylink>>> links = CreateLinks(fromPointKey(allClusters));
  if(Verbosity()>0) cout << "created links\n";
  vector<vector<keylink>> bidirectionalLinks = FindBiLinks(links.first,links.second);
  if(Verbosity()>0) cout << "found bilinks\n";
  // extract involved clusters (and locations) from bi-links
  // std::set::insert automatically skips duplicates
  vector<coordKey> clusterCands;
  for(vector<keylink> link_list : bidirectionalLinks)
  {
    for(keylink link : link_list)
    {
      if(!any_of(clusterCands.begin(),clusterCands.end(),[&](coordKey k){return link[0]==k;})) clusterCands.push_back(link[0]);
      if(!any_of(clusterCands.begin(),clusterCands.end(),[&](coordKey k){return link[1]==k;})) clusterCands.push_back(link[1]);
    }
  }
  return clusterCands;
}

int PHCASeeding::FindSeedsWithMerger()
{
  vector<pointKey> allClusters;
  vector<unordered_set<keylink>> belowLinks;
  vector<unordered_set<keylink>> aboveLinks;
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
  if(Verbosity()>0) cout << "allClusters search time: " << t_seed->get_accumulated_time() / 1000 << " s" << endl;
  LogDebug(" number of clusters: " << allClusters.size() << endl);
  t_seed->restart();

  pair<vector<unordered_set<keylink>>,vector<unordered_set<keylink>>> links = CreateLinks(fromPointKey(allClusters));
if(Verbosity()>1)  cout << "links size: " << links.first.size() << "| " << links.second.size() << endl;
  vector<vector<keylink>> biLinks = FindBiLinks(links.first,links.second);
if(Verbosity()>1)  cout << "bilinks size: " << biLinks.size() << endl;
  vector<keylist> trackSeedKeyLists = FollowBiLinks(biLinks);
if(Verbosity()>1)  cout << "keylistvector size: " << trackSeedKeyLists.size() << endl;
//  if(Verbosity()>0)  std::cout << "seeds before merge: " << trackSeedKeyLists.size() << "\n";
//  vector<keylist> mergedSeedKeyLists = MergeSeeds(trackSeedKeyLists);
//  if(Verbosity()>0) std::cout << "seeds after merge round 1: " << mergedSeedKeyLists.size() << "\n";
//  mergedSeedKeyLists = MergeSeeds(mergedSeedKeyLists);
//  if(Verbosity()>0) std::cout << "seeds after merge round 2: " << mergedSeedKeyLists.size() << "\n";
  vector<SvtxTrack_v2> seeds = fitter->ALICEKalmanFilter(trackSeedKeyLists,true);
  publishSeeds(seeds);
  return seeds.size();
}

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

pair<vector<unordered_set<keylink>>,vector<unordered_set<keylink>>> PHCASeeding::CreateLinks(vector<coordKey> clusters, int mode)
{
  size_t nclusters = 0;

  double cluster_find_time = 0;
  double rtree_query_time = 0;
  double transform_time = 0;
  double compute_best_angle_time = 0;
  double set_insert_time = 0;

  vector<unordered_set<keylink>> belowLinks;
  vector<unordered_set<keylink>> aboveLinks;
  belowLinks.resize(_nlayers_tpc);
  aboveLinks.resize(_nlayers_tpc);

  for (vector<coordKey>::iterator StartCluster = clusters.begin(); StartCluster != clusters.end(); ++StartCluster)
  {
    nclusters++;
    // get clusters near this one in adjacent layers
    double StartPhi = StartCluster->first[0];
    double StartEta = StartCluster->first[1];
    unsigned int StartLayer = StartCluster->first[2];
    if(StartLayer < _start_layer) continue;
    if(StartLayer > _end_layer) continue;
    TrkrCluster* StartCl = _cluster_map->findCluster(StartCluster->second);
    double StartX = StartCl->getPosition(0)-_vertex->get_x();
    double StartY = StartCl->getPosition(1)-_vertex->get_y();
    double StartZ = StartCl->getPosition(2)-_vertex->get_z();
    t_seed->stop();
    cluster_find_time += t_seed->elapsed();
    t_seed->restart();
    LogDebug(" starting cluster:" << endl);
    LogDebug(" eta: " << StartEta << endl);
    LogDebug(" phi: " << StartPhi << endl);
    LogDebug(" layer: " << StartLayer << endl);

    vector<pointKey> ClustersAbove;
    vector<pointKey> ClustersBelow;
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
    LogDebug(" entries in below layer: " << ClustersBelow.size() << endl);
    LogDebug(" entries in above layer: " << ClustersAbove.size() << endl);
    vector<array<double,3>> delta_below;
    vector<array<double,3>> delta_above;
    delta_below.clear();
    delta_above.clear();
    delta_below.resize(ClustersBelow.size());
    delta_above.resize(ClustersAbove.size());
    // calculate (delta_eta, delta_phi) vector for each neighboring cluster

    transform(ClustersBelow.begin(),ClustersBelow.end(),delta_below.begin(),
      [&](pointKey BelowCandidate){
        TrkrCluster* BelowCl = _cluster_map->findCluster(BelowCandidate.second);
        return array<double,3>{(BelowCl->getPosition(0)-_vertex->get_x())-StartX,
          (BelowCl->getPosition(1)-_vertex->get_y())-StartY,
          (BelowCl->getPosition(2)-_vertex->get_z())-StartZ};});

    transform(ClustersAbove.begin(),ClustersAbove.end(),delta_above.begin(),
      [&](pointKey AboveCandidate){
        TrkrCluster* AboveCl = _cluster_map->findCluster(AboveCandidate.second);
        return array<double,3>{(AboveCl->getPosition(0)-_vertex->get_x())-StartX,
          (AboveCl->getPosition(1)-_vertex->get_y())-StartY,
          (AboveCl->getPosition(2)-_vertex->get_z())-StartZ};});
    t_seed->stop();
    transform_time += t_seed->elapsed();
    t_seed->restart();

    // find the three clusters closest to a straight line
    // (by maximizing the cos of the angle between the (delta_eta,delta_phi) vectors)
    double maxCosPlaneAngle = -0.9;
    //double minSumLengths = 1e9;
    coordKey bestBelowCluster = make_pair(array<float,3>({0.,0.,-1e9}),0);
    coordKey bestAboveCluster = make_pair(array<float,3>({0.,0.,-1e9}),0);
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
        vector<pointKey> clustersTwoLayersBelow;
        QueryTree(_rtree,
                StartPhi-_neighbor_phi_width,
                StartEta-_neighbor_eta_width,
                (double) StartLayer - 2.5,
                StartPhi+_neighbor_phi_width,
                StartEta+_neighbor_eta_width,
                (double) StartLayer - 1.5,
                clustersTwoLayersBelow);
        vector<array<double,3>> delta_2below;
        delta_2below.clear();
        delta_2below.resize(clustersTwoLayersBelow.size());
        transform(clustersTwoLayersBelow.begin(),clustersTwoLayersBelow.end(),delta_2below.begin(),
          [&](pointKey BelowCandidate){
            TrkrCluster* BelowCl = _cluster_map->findCluster(BelowCandidate.second);
            return array<double,3>{(BelowCl->getPosition(0)-_vertex->get_x())-StartX,
              (BelowCl->getPosition(1)-_vertex->get_y())-StartY,
              (BelowCl->getPosition(2)-_vertex->get_z())-StartZ};});
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
          vector<pointKey> clustersTwoLayersAbove;
          QueryTree(_rtree,
                  StartPhi-_neighbor_phi_width,
                  StartEta-_neighbor_eta_width,
                  (double) StartLayer + 1.5,
                  StartPhi+_neighbor_phi_width,
                  StartEta+_neighbor_eta_width,
                  (double) StartLayer + 2.5,
                  clustersTwoLayersAbove);
          vector<array<double,3>> delta_2above;
          delta_2above.clear();
          delta_2above.resize(clustersTwoLayersAbove.size());
          transform(clustersTwoLayersAbove.begin(),clustersTwoLayersAbove.end(),delta_2above.begin(),
            [&](pointKey AboveCandidate){
              TrkrCluster* AboveCl = _cluster_map->findCluster(AboveCandidate.second);
              return array<double,3>{(AboveCl->getPosition(0)-_vertex->get_x())-StartX,
                (AboveCl->getPosition(1)-_vertex->get_y())-StartY,
                (AboveCl->getPosition(2)-_vertex->get_z())-StartZ};});
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
    LogDebug(" max collinearity: " << maxCosPlaneAngle << endl);
  }
  t_seed->stop();
  if(Verbosity()>0)
  {
    cout << "triplet forming time: " << t_seed->get_accumulated_time() / 1000 << " s" << endl;
    cout << "starting cluster setup: " << cluster_find_time / 1000 << " s" << endl;
    cout << "RTree query: " << rtree_query_time /1000 << " s" << endl;
    cout << "Transform: " << transform_time /1000 << " s" << endl;
    cout << "Compute best triplet: " << compute_best_angle_time /1000 << " s" << endl;
    cout << "Set insert: " << set_insert_time /1000 << " s" << endl;
  }
  t_seed->restart();

  return make_pair(belowLinks,aboveLinks);
}

vector<vector<keylink>> PHCASeeding::FindBiLinks(vector<unordered_set<keylink>> belowLinks,vector<unordered_set<keylink>> aboveLinks)
{
  // remove all triplets for which there isn't a mutual association between two clusters
  vector<vector<keylink>> bidirectionalLinks;
  bidirectionalLinks.resize(_nlayers_tpc);
  for(int layer = _nlayers_tpc-1; layer > 0; --layer)
  {
    for(unordered_set<keylink>::iterator belowLink = belowLinks[layer].begin(); belowLink != belowLinks[layer].end(); ++belowLink)
    {
      if((*belowLink)[1].second==0) continue;
      unsigned int end_layer_index = TrkrDefs::getLayer((*belowLink)[1].second) - (_nlayers_intt + _nlayers_maps);
      keylink reversed = {(*belowLink)[1],(*belowLink)[0]};
      unordered_set<keylink>::iterator sameAboveLinkExists = aboveLinks[end_layer_index].find(reversed);
      if(sameAboveLinkExists != aboveLinks[end_layer_index].end())
      {
        bidirectionalLinks[layer].push_back((*belowLink));
      }
    }
  }
  t_seed->stop();
  if(Verbosity()>0) cout << "bidirectional link forming time: " << t_seed->get_accumulated_time() / 1000 << " s" << endl;
  t_seed->restart();

  return bidirectionalLinks;
}

vector<keylist> PHCASeeding::FollowBiLinks(vector<vector<keylink>> bidirectionalLinks)
{
  // follow bidirectional links to form lists of cluster keys
  // (to be fitted for track seed parameters)
  vector<keylist> trackSeedKeyLists;
  // get starting cluster keys, create a keylist for each
  // (only check last element of each pair because we start from the outer layers and go inward)
  for(unsigned int layer = 0; layer < _nlayers_tpc-1; ++layer)
  {
    for(vector<keylink>::iterator startCand = bidirectionalLinks[layer].begin(); startCand != bidirectionalLinks[layer].end(); ++startCand)
    {
      bool has_above_link = false;
      unsigned int imax = 1;
      if(layer==_nlayers_tpc-2) imax = 1;
      for(unsigned int i=1;i<=imax;i++)
      {
        has_above_link = has_above_link || any_of(bidirectionalLinks[layer+i].begin(),bidirectionalLinks[layer+i].end(),[&](keylink k){return (*startCand)[0]==k[1];});
      }
//      for(vector<keylink>::iterator testlink = bidirectionalLinks[layer+1].begin(); !has_above_link && (testlink != bidirectionalLinks[layer+1].end()); ++testlink)
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
  if(Verbosity()>0) cout << "starting cluster finding time: " << t_seed->get_accumulated_time() / 1000 << " s" << endl;
  t_seed->restart();
  // assemble track cluster chains from starting cluster keys (ordered from outside in)
  for(vector<keylist>::iterator trackKeyChain = trackSeedKeyLists.begin(); trackKeyChain != trackSeedKeyLists.end(); ++trackKeyChain)
  {
    bool reached_end = false;
    while(!reached_end)
    {
      TrkrDefs::cluskey trackHead = trackKeyChain->back();
      unsigned int trackHead_layer = TrkrDefs::getLayer(trackHead) - (_nlayers_intt + _nlayers_maps);
      bool no_next_link = true;
      for(vector<keylink>::iterator testlink = bidirectionalLinks[trackHead_layer].begin(); testlink != bidirectionalLinks[trackHead_layer].end(); ++testlink)
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
  if(Verbosity()>0) cout << "keychain assembly time: " << t_seed->get_accumulated_time() / 1000 << " s" << endl;
  t_seed->restart();
  LogDebug(" track key chains assembled: " << trackSeedKeyLists.size() << endl);
  LogDebug(" track key chain lengths: " << endl);
  for(vector<keylist>::iterator trackKeyChain = trackSeedKeyLists.begin(); trackKeyChain != trackSeedKeyLists.end(); ++trackKeyChain)
  {
    LogDebug(" " << trackKeyChain->size() << endl);
  }
  int jumpcount = 0;
  LogDebug(" track key associations:" << endl);
  for(size_t i=0;i<trackSeedKeyLists.size();++i)
  {
    LogDebug(" seed " << i << ":" << endl);

    double lasteta = -100;
    double lastphi = -100;
    for(size_t j=0;j<trackSeedKeyLists[i].size();++j)
    {
      TrkrCluster* cl = _cluster_map->findCluster(trackSeedKeyLists[i][j]);
      TVector3 vec(cl->getPosition(0)-_vertex->get_x(), cl->getPosition(1)-_vertex->get_y(), cl->getPosition(2)-_vertex->get_z());

      double clus_phi = vec.Phi();
      if(clus_phi<0) clus_phi = 2*M_PI + clus_phi;
      //clus_phi -= 2 * M_PI * floor(clus_phi / (2 * M_PI));
      double clus_eta = vec.Eta();
      double etajump = clus_eta-lasteta;
      double phijump = clus_phi-lastphi;
      #if defined(_DEBUG_) 
      unsigned int lay = TrkrDefs::getLayer(trackSeedKeyLists[i][j].second);
      #endif
      if((fabs(etajump)>0.1 && lasteta!=-100) || (fabs(phijump)>1 && lastphi!=-100))
	{
           LogDebug(" Eta or Phi jump too large! " << endl);
           ++jumpcount;
        }
      LogDebug(" (eta,phi,layer) = (" << clus_eta << "," << clus_phi << "," << lay << ") " <<
        " (x,y,z) = (" << cl->getPosition(0)-_vertex->get_x() << "," << cl->getPosition(1)-_vertex->get_y() << "," << cl->getPosition(2)-_vertex->get_z() << ")" << endl);
      
      if(Verbosity() > 0)
	{
            unsigned int lay = TrkrDefs::getLayer(trackSeedKeyLists[i][j]);
            std::cout << "  eta, phi, layer = (" << clus_eta << "," << clus_phi << "," << lay << ") " <<
             " (x,y,z) = (" << cl->getPosition(0)-_vertex->get_x() << "," << cl->getPosition(1)-_vertex->get_y() << "," << cl->getPosition(2)-_vertex->get_z() << ")" << std::endl;
        }
      lasteta = clus_eta;
      lastphi = clus_phi;
    }
  }
  LogDebug(" Total large jumps: " << jumpcount << endl);
  t_seed->stop();
  if(Verbosity()>0) cout << "eta-phi sanity check time: " << t_seed->get_accumulated_time() / 1000 << " s" << endl;
  t_seed->restart();
  return trackSeedKeyLists;
}

vector<keylist> PHCASeeding::MergeSeeds(vector<keylist> seeds)
{
  if(Verbosity()>0) std::cout << "entered merge\n";
  //initialize vector of flags specifying whether seed is used
  vector<bool> isUsed(seeds.size());
  std::fill(isUsed.begin(),isUsed.end(),false);
  if(Verbosity()>0) std::cout << "filled used vector\n";
  //get all seed ends
  vector<pointKey> seedEnds;
  vector<pointKey> frontEnds;
  vector<pointKey> backEnds;
  for(auto seed : seeds)
  {
    pointKey frontEnd = makepointKey(seed.front());
    pointKey backEnd = makepointKey(seed.back());
    seedEnds.push_back(frontEnd);
    seedEnds.push_back(backEnd);
    frontEnds.push_back(frontEnd);
    backEnds.push_back(backEnd);
  }
  if(Verbosity()>0) std::cout << "gotten seed ends\n";
  //make RTree with seed ends
  _rtree.clear();
  FillTree(seedEnds);
  if(Verbosity()>0) std::cout << "filled rtree\n";
  // find seeds that have similar eta, phi, with one layer between them
  vector<keylist> merged;
  for(size_t i=0;i<backEnds.size();i++)
  {
    pointKey backend = backEnds[i];
    vector<pointKey> brokenCandidates;
    double phi = backend.first.get<0>();
    double eta = backend.first.get<1>();
    double layer = backend.first.get<2>();
    QueryTree(_rtree,
      phi-_neighbor_phi_width,
      eta-_neighbor_eta_width,
      layer-3.5,
      phi+_neighbor_phi_width,
      eta+_neighbor_eta_width,
      layer-0.5,
      brokenCandidates);
    for(auto cand : brokenCandidates)
    {
      // find corresponding seed for starting cluster
      auto iter = std::find_if(frontEnds.begin(),frontEnds.end(),
        [&](const pointKey &k){return cand.second == k.second;});
      if(iter == frontEnds.end()) continue;
      int seed_index = std::distance(frontEnds.begin(),iter);
      // create merged seed
      keylist mergedseed;
      mergedseed.insert(mergedseed.end(),seeds[i].begin(),seeds[i].end());
      mergedseed.insert(mergedseed.end(),seeds[seed_index].begin(),seeds[seed_index].end());
      merged.push_back(mergedseed);
      // mark seeds as used
      isUsed[i] = true;
      isUsed[seed_index] = true;
    }
  }
  // include unused seeds as-is
  for(size_t i=0;i<seeds.size();i++)
  {
    if(!isUsed[i]) merged.push_back(seeds[i]);
  }
  return merged;
}

void PHCASeeding::publishSeeds(vector<SvtxTrack_v2> seeds)
{
  if(Verbosity()>1) cout << "publishing " << seeds.size() << " seeds" << endl;
  for(size_t i=0;i<seeds.size();i++)
  {
    _track_map->insert(&(seeds[i]));
  }
}

int PHCASeeding::Setup(PHCompositeNode *topNode)
{
  if(Verbosity()>0) cout << "Called Setup" << endl;
  if(Verbosity()>0) cout << "topNode:" << topNode << endl;
  PHTrackSeeding::Setup(topNode);
  InitializeGeometry(topNode);
 #if __cplusplus < 201402L
  t_fill = boost::make_unique<PHTimer>("t_fill");
  t_seed = boost::make_unique<PHTimer>("t_seed");
#else
  t_fill = std::make_unique<PHTimer>("t_fill");
  t_seed = std::make_unique<PHTimer>("t_seed");
#endif
  t_fill->stop();
  t_seed->stop();
  fitter = std::make_shared<ALICEKF>(topNode,_cluster_map,_fieldDir,_min_clusters_per_track,_max_sin_phi,Verbosity());
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHCASeeding::End()
{
  if(Verbosity()>0) cout << "Called End " << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
