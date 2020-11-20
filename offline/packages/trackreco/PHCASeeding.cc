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
#include <trackbase_historic/SvtxTrack_v1.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>

#include <trackbase/TrkrCluster.h>  // for TrkrCluster
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>  // for getLayer, clu...
#include <trackbase/TrkrClusterHitAssoc.h>

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
    float cluster_z_error,
    float cluster_alice_y_error,
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
  , _cluster_z_error(cluster_z_error)
  , _cluster_alice_y_error(cluster_alice_y_error)
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

  TrkrClusterContainer::ConstRange clusrange = _cluster_map->getClusters();

  for (TrkrClusterContainer::ConstIterator iter = clusrange.first; iter != clusrange.second; ++iter)
  {
    TrkrCluster *cluster = iter->second;
    TrkrDefs::cluskey ckey = iter->first;
    unsigned int layer = TrkrDefs::getLayer(ckey);
    if (layer < _start_layer || layer >= _end_layer) continue;

    if(!_use_truth_clusters)
      {
	TrkrClusterHitAssoc::ConstRange hitrange = _cluster_hit_map->getHits(ckey);
	unsigned int nhits = std::distance(hitrange.first,hitrange.second);
	if(nhits<_min_nhits_per_cluster) continue;
      }

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
    TrkrClusterHitAssoc::ConstRange hitrange = _cluster_hit_map->getHits(iter->second);
    unsigned int nhits = std::distance(hitrange.first,hitrange.second);
    if(nhits<_min_nhits_per_cluster) continue;
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
  if(Verbosity()>0)   cout << "number of seeds " << numberofseeds << endl;
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
  vector<vector<keylink>> biLinks = FindBiLinks(links.first,links.second);
  vector<keylist> trackSeedKeyLists = FollowBiLinks(biLinks);
  if(Verbosity()>0)  std::cout << "seeds before merge: " << trackSeedKeyLists.size() << "\n";
  vector<keylist> mergedSeedKeyLists = MergeSeeds(trackSeedKeyLists);
  if(Verbosity()>0) std::cout << "seeds after merge round 1: " << mergedSeedKeyLists.size() << "\n";
  mergedSeedKeyLists = MergeSeeds(mergedSeedKeyLists);
  if(Verbosity()>0) std::cout << "seeds after merge round 2: " << mergedSeedKeyLists.size() << "\n";
  int nseeds = ALICEKalmanFilter(mergedSeedKeyLists);
  return nseeds;
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

bool checknan(float val, std::string name, int num)
{
  if(std::isnan(val))
  {
    std::cout << "WARNING: " << val << " is NaN for seed " << num << ". Aborting this seed.\n";
  }
  return std::isnan(val);
}

int PHCASeeding::ALICEKalmanFilter(vector<keylist> trackSeedKeyLists)
{
  int nseeds = 0;
  if(Verbosity()>0) std::cout << "min clusters per track: " << _min_clusters_per_track << "\n";
  for(vector<keylist>::iterator trackKeyChain = trackSeedKeyLists.begin(); trackKeyChain != trackSeedKeyLists.end(); ++trackKeyChain)
  {
    if(trackKeyChain->size() < _min_clusters_per_track) continue;
    // get starting cluster from key
    TrkrCluster* startCluster = _cluster_map->findCluster(trackKeyChain->at(0));
    // Transform sPHENIX coordinates into ALICE-compatible coordinates
    double x0 = startCluster->getPosition(0);
    double y0 = startCluster->getPosition(1);
    double z0 = startCluster->getPosition(2);
    LogDebug("Initial (x,y,z): (" << x0 << "," << y0 << "," << z0 << ")" << endl);
    // ALICE x coordinate = distance from beampipe
    float alice_x0 = sqrt(x0*x0+y0*y0);
    float alice_y0 = 0;
    float alice_z0 = z0;
    // Initialize track and linearisation
    GPUTPCTrackParam trackSeed;
    trackSeed.InitParam();
    trackSeed.SetX(alice_x0);
    trackSeed.SetY(alice_y0);
    trackSeed.SetZ(alice_z0);
    float x = x0;
    float y = y0;
    #if defined(_DEBUG_)
    float z = z0;
    float alice_x = sqrt(x0*x0+y0*y0);
    #endif
    float trackCartesian_x = 0.;
    float trackCartesian_y = 0.;
    // float trackCartesian_z = 0.;
    // Pre-set momentum-based parameters to improve numerical stability
    TrkrCluster* SecondCluster = _cluster_map->findCluster(trackKeyChain->at(1));
    float second_x = SecondCluster->getPosition(0);
    float second_y = SecondCluster->getPosition(1);
    float second_z = SecondCluster->getPosition(2);
    float second_alice_x = sqrt(second_x*second_x+second_y*second_y);
    float delta_alice_x = second_alice_x - alice_x0;
    float first_phi = atan2(y0,x0);
    float second_alice_y = (second_x/cos(first_phi)-second_y/sin(first_phi))/(sin(first_phi)/cos(first_phi)+cos(first_phi)/sin(first_phi));
    float init_SinPhi = second_alice_y / sqrt(delta_alice_x*delta_alice_x + second_alice_y*second_alice_y);
    float delta_z = second_z - z0;
    float init_DzDs = delta_z / sqrt(delta_alice_x*delta_alice_x + second_alice_y*second_alice_y);
    trackSeed.SetSinPhi(init_SinPhi);
    LogDebug("Set initial SinPhi to " << init_SinPhi << endl);
    trackSeed.SetDzDs(init_DzDs);
    LogDebug("Set initial DzDs to " << init_DzDs << endl);
    GPUTPCTrackLinearisation trackLine(trackSeed);

    LogDebug(endl << endl << "------------------------" << endl << "seed size: " << trackKeyChain->size() << endl << endl << endl);
    int cluster_ctr = 1;
    // starting at second cluster, perform track propagation
    for(keylist::iterator clusterkey = next(trackKeyChain->begin()); clusterkey != trackKeyChain->end(); ++clusterkey)
    {
      LogDebug("cluster " << cluster_ctr << " -> " << cluster_ctr + 1 << endl);
      LogDebug("this cluster (x,y,z) = (" << x << "," << y << "," << z << ")" << endl);
      // get cluster from key
      TrkrCluster* nextCluster = _cluster_map->findCluster(*clusterkey);
      // find ALICE x-coordinate
      float nextCluster_x = nextCluster->getPosition(0);
      float nextCluster_y = nextCluster->getPosition(1);
      float nextCluster_z = nextCluster->getPosition(2);
      float nextAlice_x = sqrt(nextCluster_x*nextCluster_x+nextCluster_y*nextCluster_y);
      // rotate track coordinates to match orientation of next cluster
      float newPhi = atan2(nextCluster_y,nextCluster_x);
      LogDebug("new phi = " << newPhi << endl);
      float oldPhi = atan2(y,x);
      LogDebug("old phi = " << oldPhi << endl);
      float alpha = newPhi - oldPhi;
      LogDebug("alpha = " << alpha << endl);
      if(!trackSeed.Rotate(-alpha,trackLine,_max_sin_phi))
      {
        LogWarning("Rotate failed! Aborting for this seed...\n");
        break;
      }
      LogDebug("track coordinates (ALICE) after rotation: (" << trackSeed.GetX() << "," << trackSeed.GetY() << "," << trackSeed.GetZ() << ")" << endl);
      LogDebug("Transporting from " << alice_x << " to " << nextAlice_x << "...");
      if(!trackSeed.TransportToX(nextAlice_x,trackLine,_Bz,_max_sin_phi))
      {
        LogWarning("Transport failed! Aborting for this seed...\n");
        break;
      }
      // convert ALICE coordinates to sPHENIX cartesian coordinates, for debugging
      float predicted_alice_x = trackSeed.GetX();
      LogDebug("new track ALICE x = " << trackSeed.GetX() << endl);
      float predicted_alice_y = trackSeed.GetY();
      LogDebug("new track ALICE y = " << trackSeed.GetY() << endl);
      // float predicted_z = trackSeed.GetZ();
      LogDebug("new track z = " << trackSeed.GetZ() << endl);
      float cos_phi = x/sqrt(x*x+y*y);
      LogDebug("cos_phi = " << cos_phi << endl);
      float sin_phi = y/sqrt(x*x+y*y);
      LogDebug("sin phi = " << sin_phi << endl);
      trackCartesian_x = predicted_alice_x*cos_phi+predicted_alice_y*sin_phi;
      trackCartesian_y = predicted_alice_x*sin_phi-predicted_alice_y*cos_phi;
      // trackCartesian_z = predicted_z;
      LogDebug("Track transported to (x,y,z) = (" << trackCartesian_x << "," << trackCartesian_y << "," << trackCartesian_z << ")" << endl);
      LogDebug("Next cluster is at (x,y,z) = (" << nextCluster_x << "," << nextCluster_y << "," << nextCluster_z << ")" << endl);
      LogDebug("track coordinates (ALICE) after rotation: (" << trackSeed.GetX() << "," << trackSeed.GetY() << "," << trackSeed.GetZ() << ")" << endl);
      //float nextCluster_alice_y = (nextCluster_x/cos(newPhi) - nextCluster_y/sin(newPhi))/(tan(newPhi)+1./tan(newPhi));
      float nextCluster_alice_y = 0.;
      LogDebug("next cluster ALICE y = " << nextCluster_alice_y << endl);
      float y2_error = _cluster_alice_y_error*_cluster_alice_y_error;
      float z2_error = _cluster_z_error*_cluster_z_error;
      LogDebug("track ALICE SinPhi = " << trackSeed.GetSinPhi() << endl);
      // Apply Kalman filter
      if(!trackSeed.Filter(nextCluster_alice_y,nextCluster_z,y2_error,z2_error,_max_sin_phi))
      {
	if (Verbosity() >= 1)
	  LogError("Kalman filter failed for seed " << nseeds << "! Aborting for this seed..." << endl);
        break;
      }
      
      x = nextCluster_x;
      y = nextCluster_y;
      #if defined(_DEBUG_)
      z = nextCluster_z;
      alice_x = nextAlice_x;
      #endif
      ++cluster_ctr;
    }
    
    if(!trackSeed.CheckNumericalQuality())
    {
      cout << "ERROR: Track seed failed numerical quality check before conversion to sPHENIX coordinates! Skipping this one.\n";
      continue;
    } 
    
    //    pt:z:dz:phi:dphi:c:dc
    // Fill NT with track parameters
    // float StartEta = -log(tan(atan(z0/sqrt(x0*x0+y0*y0))));
    float track_pt = fabs( 1./(trackSeed.GetQPt()));
    if(checknan(track_pt,"pT",nseeds)) continue;
    float track_pterr = sqrt(trackSeed.GetErr2QPt())/(trackSeed.GetQPt()*trackSeed.GetQPt());
    if(checknan(track_pterr,"pT err",nseeds)) continue;
    LogDebug("Track pterr = " << track_pterr << endl);
    float track_z = trackSeed.GetZ();
    if(checknan(track_z,"z",nseeds)) continue;
    float track_zerr = sqrt(trackSeed.GetErr2Z());
    if(checknan(track_zerr,"zerr",nseeds)) continue;
    float track_phi = atan2(trackCartesian_y,trackCartesian_x);
    if(checknan(track_phi,"phi",nseeds)) continue;
    float last_cluster_phierr = _cluster_map->findCluster(trackKeyChain->back())->getPhiError();
    // phi error assuming error in track radial coordinate is zero
    float track_phierr = sqrt(pow(last_cluster_phierr,2)+(pow(trackSeed.GetX(),2)*trackSeed.GetErr2Y()) / 
      pow(pow(trackSeed.GetX(),2)+pow(trackSeed.GetY(),2),2));
    if(checknan(track_phierr,"phierr",nseeds)) continue;
    LogDebug("Track phi = " << track_phi << endl);
    LogDebug("Track phierr = " << track_phierr << endl);
    float track_curvature = trackSeed.GetKappa(_Bz);
    if(checknan(track_curvature,"curvature",nseeds)) continue;
    float track_curverr = sqrt(trackSeed.GetErr2QPt())*_Bz;
    if(checknan(track_curverr,"curvature error",nseeds)) continue;
    SvtxTrack_v1 track;
    track.set_id(nseeds);
    for (unsigned int j = 0; j < trackKeyChain->size(); ++j)
    {
      track.insert_cluster_key(trackKeyChain->at(j));
    }
    track.set_ndf(trackSeed.GetNDF());
    int track_charge = 0;
    if(trackSeed.GetQPt()<0) track_charge = -1 * _fieldDir;
    else track_charge = 1 * _fieldDir;
    if(Verbosity() > 0) std::cout << " trackSeed.GetQPt " << trackSeed.GetQPt() << "  track charge " << track_charge << std::endl; 
    track.set_charge(track_charge);
    TrkrCluster *cl = _cluster_map->findCluster(trackKeyChain->at(0));
    track.set_x(cl->getX());  //track.set_x(cl->getX());
    track.set_y(cl->getY());  //track.set_y(cl->getY());
    track.set_z(cl->getZ());  //track.set_z(cl->getZ());
    float s = sin(track_phi);
    float c = cos(track_phi);
    float p = trackSeed.GetSinPhi();
    if(checknan(p,"ALICE sinPhi",nseeds)) continue;
    float d = trackSeed.GetDzDs();
    if(checknan(d,"ALICE dz/ds",nseeds)) continue;
    float pY = track_pt*p;
    float pX = sqrt(track_pt*track_pt-pY*pY);
    track.set_px(pX*c-pY*s);
    track.set_py(pX*s+pY*c);
    track.set_pz(track_pt * trackSeed.GetDzDs()); 
    const float* cov = trackSeed.GetCov();
    for(int i=0;i<15;i++)
    {
      if(checknan(cov[i],"covariance element "+std::to_string(i),nseeds)) continue;
    }
    // make this into an actual Eigen matrix
    Eigen::Matrix<float,5,5> ecov;
    ecov(0,0)=cov[0];
    ecov(0,1)=cov[1];
    ecov(0,2)=cov[2];
    ecov(0,3)=cov[3];
    ecov(0,4)=cov[4];
    ecov(1,1)=cov[5];
    ecov(1,2)=cov[6];
    ecov(1,3)=cov[7];
    ecov(1,4)=cov[8];
    ecov(2,2)=cov[9];
    ecov(2,3)=cov[10];
    ecov(2,4)=cov[11];
    ecov(3,3)=cov[12];
    ecov(3,4)=cov[13];
    ecov(4,4)=cov[14];
    // symmetrize
    ecov(1,0)=ecov(0,1);
    ecov(2,0)=ecov(0,2);
    ecov(3,0)=ecov(0,3);
    ecov(4,0)=ecov(0,4);
    ecov(2,1)=ecov(1,2);
    ecov(3,1)=ecov(1,3);
    ecov(4,1)=ecov(1,4);
    ecov(3,2)=ecov(2,3);
    ecov(4,2)=ecov(2,4);
    ecov(4,3)=ecov(3,4);
    // make rotation matrix based on the following:
    // x = X*cos(track_phi) - Y*sin(track_phi)
    // y = X*sin(track_phi) + Y*cos(track_phi)
    // z = Z
    // pY = pt*sinphi
    // pX = sqrt(pt**2 - pY**2)
    // px = pX*cos(track_phi) - pY*sin(track_phi)
    // py = pX*sin(track_phi) + pY*cos(track_phi)
    // pz = pt*(dz/ds)
    Eigen::Matrix<float,6,5> J;
    J(0,0) = -s; // dx/dY
    J(0,1) = 0.; // dx/dZ
    J(0,2) = 0.; // dx/d(sinphi)
    J(0,3) = 0.; // dx/d(dz/ds)
    J(0,4) = 0.; // dx/d(Q/pt)

    J(1,0) = c;  // dy/dY
    J(1,1) = 0.; // dy/dZ
    J(1,2) = 0.; // dy/d(sinphi)
    J(1,3) = 0.; // dy/d(dz/ds)
    J(1,4) = 0.; // dy/d(Q/pt)

    J(2,0) = 0.; // dz/dY
    J(2,1) = 1.; // dz/dZ
    J(2,2) = 0.; // dz/d(sinphi)
    J(2,3) = 0.; // dz/d(dz/ds)
    J(2,4) = 0.; // dz/d(Q/pt)

    J(3,0) = 0.; // dpx/dY
    J(3,1) = 0.; // dpx/dZ
    J(3,2) = -track_pt*(p*c/sqrt(1-p*p)+s); // dpx/d(sinphi)
    J(3,3) = 0.; // dpx/d(dz/ds)
    J(3,4) = track_pt*track_pt*track_charge*(p*s-c*sqrt(1-p*p)); // dpx/d(Q/pt)

    J(4,0) = 0.; // dpy/dY
    J(4,1) = 0.; // dpy/dZ
    J(4,2) = track_pt*(c-p*s/sqrt(1-p*p)); // dpy/d(sinphi)
    J(4,3) = 0.; // dpy/d(dz/ds)
    J(4,4) = -track_pt*track_pt*track_charge*(p*c+s*sqrt(1-p*p)); // dpy/d(Q/pt)

    J(5,0) = 0.; // dpz/dY
    J(5,1) = 0.; // dpz/dZ
    J(5,2) = 0.; // dpz/d(sinphi)
    J(5,3) = track_pt; // dpz/d(dz/ds)
    J(5,4) = -track_pt*track_pt*track_charge*d; // dpz/d(Q/pt)

    for(int i=0;i<6;i++)
    {
      for(int j=0;j<5;j++)
      {
        if(checknan(J(i,j),"covariance rotator element ("+std::to_string(i)+","+std::to_string(j)+")",nseeds)) continue;
      }
    }

    // the heavy lifting happens here
    Eigen::Matrix<float,6,6> scov = J*ecov*J.transpose();
    
    // fill SvtxTrack covariance matrix with results
    for(int i=0;i<6;i++)
    {
      for(int j=0;j<6;j++)
      {
        track.set_error(i, j, scov(i,j));
      }
    }
/*
    // Proceed with the absolutely hellish coordinate transformation of the covariance matrix.
    // Derived from:
    // 1) Taking the Jacobian of the conversion from (Y,Z,SinPhi,DzDs,Q/Pt) to (x,y,z,px,py,pz)
    // 2) Computing (Jacobian)*(ALICE covariance matrix)*(transpose of Jacobian)
    track.set_error(0, 0, cov[0]*s*s);
    track.set_error(0, 1, -cov[0]*c*s);
    track.set_error(0, 2, -cov[1]*s);
    track.set_error(0, 3, cov[2]*s*s/q-cov[4]*s*(-c/(q*q)+p*s/(q*q)));
    track.set_error(0, 4, -cov[2]*c*s/q-cov[4]*s*(-c*p/(q*q)-s/(q*q)));
    track.set_error(0, 5, cov[4]*d*s/(q*q)-cov[3]*s/q);
    track.set_error(1, 1, cov[0]*c*c);
    track.set_error(1, 2, cov[1]*c);
    track.set_error(1, 3, -cov[2]*c*s/q+cov[4]*c*(-c/(q*q)+p*s/(q*q)));
    track.set_error(1, 4, cov[2]*c*c/q+cov[4]*c*(-c*p/(q*q)-s/(q*q)));
    track.set_error(1, 5, cov[4]*d*c/(q*q)+cov[3]*c/q);
    track.set_error(2, 2, cov[5]);
    track.set_error(2, 3, -cov[6]*s/q+cov[8]*(-c/(q*q)+p*s/(q*q)));
    track.set_error(2, 4, cov[6]*c/q+cov[8]*(-c*p/(q*q)-s/(q*q)));
    track.set_error(2, 5, -cov[8]*d/(q*q)+cov[7]/q);
    track.set_error(3, 3, cov[9]*s*s/(q*q)-cov[11]*(-c/(q*q*q)+p*s/(q*q*q)) + (-c/(q*q)+p*s/(q*q))*(-cov[11]*s/q+cov[14]*(-c/(q*q)+p*s/(q*q))));
    track.set_error(3, 4, -cov[9]*c*s/(q*q)+cov[11]*(-c/(q*q*q)+p*s/(q*q*q)) + (-c*p/(q*q)-s/(q*q))*(-cov[11]*s/q+cov[14]*(-c/(q*q)+p*s/(q*q))));
    track.set_error(3, 5, -cov[10]*s/(q*q)+cov[13]/q*(-c/(q*q)+p*s/(q*q))-d/(q*q)*(-cov[11]*s/q+cov[14]*(-c/(q*q)+p*s/(q*q))));
    track.set_error(4, 4, c/q*(c/q*cov[9]+cov[11]*(-c*p/(q*q)-s/(q*q)))+(-c*p/(q*q)-s/(q*q))*(c/q*cov[11]+cov[14]*(-c*p/(q*q)-s/(q*q))));
    track.set_error(4, 5, cov[10]*c/(q*q)+cov[13]/q*(-c*p/(q*q)-s/(q*q))-d/(q*q)*(c/q*cov[11]+cov[14]*(-c*p/(q*q)-s/(q*q))));
    track.set_error(5, 5, -d/(q*q)*(-d*cov[14]/(q*q)+cov[13]/q)-d*cov[13]/(q*q*q)+cov[12]/(q*q));
    // symmetrize covariance
    track.set_error(1, 0, track.get_error(0, 1));
    track.set_error(2, 0, track.get_error(0, 2));
    track.set_error(3, 0, track.get_error(0, 3));
    track.set_error(4, 0, track.get_error(0, 4));
    track.set_error(5, 0, track.get_error(0, 5));
    track.set_error(2, 1, track.get_error(1, 2));
    track.set_error(3, 1, track.get_error(1, 3));
    track.set_error(4, 1, track.get_error(1, 4));
    track.set_error(5, 1, track.get_error(1, 5));
    track.set_error(3, 2, track.get_error(2, 3));
    track.set_error(4, 2, track.get_error(2, 4));
    track.set_error(5, 2, track.get_error(2, 5));
    track.set_error(4, 3, track.get_error(3, 4));
    track.set_error(5, 3, track.get_error(3, 5));
    track.set_error(5, 4, track.get_error(4, 5));
*/
    if(!covIsPosDef(track))
    {
      repairCovariance(track);
    }

    if(Verbosity() > 0)  track.identify();
    _track_map->insert(&track);

    ++nseeds;
  }
  return nseeds;
}

Eigen::Matrix<float,6,6> PHCASeeding::getEigenCov(SvtxTrack_v1 &track)
{
  Eigen::Matrix<float,6,6> cov;
  for(int i=0;i<6;i++)
  {
    for(int j=0;j<6;j++)
    {
      cov(i,j) = track.get_error(i,j);
    }
  }
  return cov;
}

bool PHCASeeding::covIsPosDef(SvtxTrack_v1 &track)
{
  // put covariance matrix into Eigen container
  Eigen::Matrix<float,6,6> cov = getEigenCov(track);
  // attempt Cholesky decomposition
  Eigen::LLT<Eigen::Matrix<float,6,6>> chDec(cov);
  // if Cholesky decomposition does not exist, matrix is not positive definite
  return (chDec.info() != Eigen::NumericalIssue);
}

void PHCASeeding::repairCovariance(SvtxTrack_v1 &track)
{
  // find closest positive definite matrix
  Eigen::Matrix<float,6,6> cov = getEigenCov(track);
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<float,6,6>> solver(cov);
  Eigen::Matrix<float,6,1> D = solver.eigenvalues();
  Eigen::Matrix<float,6,6> Q = solver.eigenvectors();
  Eigen::Matrix<float,6,1> Dp = D.cwiseMax(1e-6);
  Eigen::Matrix<float,6,6> Z = Q*Dp.asDiagonal()*Q.transpose();
  // updates covariance matrix
  for(int i=0;i<6;i++)
  {
    for(int j=0;j<6;j++)
    {
      track.set_error(i,j,Z(i,j));
    }
  }
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
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHCASeeding::End()
{
  if(Verbosity()>0) cout << "Called End " << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
