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

// forward declarations
class PHCompositeNode;



//#define _DEBUG_

#if defined(_DEBUG_)
#define LogDebug(exp) std::cout << "DEBUG: " << __FILE__ << ": " << __LINE__ << ": " << exp
#else
#define LogDebug(exp) (void)0
#endif

#define LogError(exp) std::cout << "ERROR: " << __FILE__ << ": " << __LINE__ << ": " << exp
#define LogWarning(exp) std::cout << "WARNING: " << __FILE__ << ": " << __LINE__ << ": " << exp

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
/*
  template<typename T, size_t N, typename L, size_t M>
  struct hash<array<pair<array<T, N>,L>,M> >
  {
    typedef array<pair<array<T, N>,L>,M> argument_type;
    typedef size_t result_type;

    result_type operator()(const argument_type& a) const
    {
      hash<T> hasher;
      hash<L> other_hasher;
      hash<pair<array<T,N>,L> 
      result_type h = 0;
      for (result_type i = 0; i < N; ++i)
      {
        h = h * 31 + hasher(a[i]);
      }
      for (result_type i = 0; i < M; ++i)
      {
        h = h * 31 + other_hasher
      return h;
    }
  };
*/
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

/*
namespace std
{
  template<typename T>
  struct hash<pair<array<float,3>, T>>
  {
    typedef pair<array<float,3>,T> argument_type;
    typedef size_t result_type;
    result_type operator()(const argument_type& a) const
    {
      hash<T> hasher;
      hash<double> d_hasher;
      return d_hasher(a.first.get<0>())+d_hasher(a.first.get<1>())+d_hasher(a.first.get<2>())+hasher(a.second);
    }
  };
}
*/

bool PHCASeeding::samepoint(point q1,point q2)
{
  return (q1.get<0>() == q2.get<0>()) && (q1.get<1>() == q2.get<1>()) && (q1.get<2>() == q2.get<2>());
}

using namespace std;
//using namespace ROOT::Minuit2;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

PHCASeeding::PHCASeeding(
    const string &name,
    unsigned int nlayers_maps,
    unsigned int nlayers_intt,
    unsigned int nlayers_tpc,
    unsigned int start_layer,
    float cluster_z_error,
    float cluster_alice_y_error,
    float neighbor_phi_width,
    float neighbor_eta_width,
    float maxSinPhi,
    float Bz)
  : PHTrackSeeding(name)
  , _g4tracks(nullptr)
  , _g4vertexes(nullptr)
  , _svtxhitsmap(nullptr)
  , _hit_used_map(nullptr)
  , _hit_used_map_size(0)
  , phisr(0.005)
  , etasr(0.0035)
  , phist(0.001)
  , etast(0.003)
  , phixt(0.008)
  , etaxt(0.005)
  , _vertex(nullptr)
  , _nlayers_maps(nlayers_maps)
  , _nlayers_intt(nlayers_intt)
  , _nlayers_tpc(nlayers_tpc)
  , _start_layer(start_layer)
  , _cluster_z_error(cluster_z_error)
  , _cluster_alice_y_error(cluster_alice_y_error)
  , _neighbor_phi_width(neighbor_phi_width)
  , _neighbor_eta_width(neighbor_eta_width)
  , _max_sin_phi(maxSinPhi)
  , _Bz(Bz)
  , _cosTheta_limit(-0.98)
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

bool PHCASeeding::samepointKey(const pointKey &p1, const pointKey &p2)
{
  bool equalCluskey = (p1.second == p2.second);
  point q1 = p1.first;
  point q2 = p2.first;
  bool equalPoint = (q1.get<0>() == q2.get<0>()) && (q1.get<1>() == q2.get<1>()) && (q1.get<2>() == q2.get<2>());
  return equalCluskey && equalPoint;
}

void PHCASeeding::QueryTree(const bgi::rtree<pointKey, bgi::quadratic<16>> &rtree, double phimin, double etamin, double lmin, double phimax, double etamax, double lmax, std::vector<pointKey> &returned_values)
{
  rtree.query(bgi::intersects(box(point(phimin, etamin, lmin), point(phimax, etamax, lmax))), std::back_inserter(returned_values));
  if (phimin < 0) rtree.query(bgi::intersects(box(point(2 * M_PI + phimin, etamin, lmin), point(2 * M_PI, etamax, lmax))), std::back_inserter(returned_values));
  if (phimax > 2 * M_PI) rtree.query(bgi::intersects(box(point(0, etamin, lmin), point(phimax - 2 * M_PI, etamax, lmax))), std::back_inserter(returned_values));
}

void PHCASeeding::FillTree()
{
  //  bgi::rtree<pointKey, bgi::quadratic<16> > rtree;
  PHTimer *t_fill = new PHTimer("t_fill");
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
    if (layer < (_nlayers_maps + _nlayers_intt + _nlayers_tpc*2./3.)) continue;

    TVector3 vec(cluster->getPosition(0) - _vertex->get_x(), cluster->getPosition(1) - _vertex->get_y(), cluster->getPosition(2) - _vertex->get_z());

    double clus_phi = vec.Phi();
    clus_phi -= 2 * M_PI * floor(clus_phi / (2 * M_PI));
    double clus_eta = vec.Eta();
    double clus_l = layer;  // _radii_all[layer];

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

  std::cout << "fill time: " << t_fill->get_accumulated_time() / 1000. << " sec" << std::endl;
  std::cout << "number of duplicates : " << n_dupli << std::endl;
}

void PHCASeeding::FillTree(vector<pointKey> clusters)
{
  // WARNING: THIS VERSION DOES NOT TEST FOR DUPLICATES!!!
  // It's much faster that way (no more querying the RTree as it's being constructed) 
  // and we can get away with it because
  // the incoming vector is guaranteed by construction to be duplicate-free.
  PHTimer *t_fill = new PHTimer("t_fill");
  t_fill->stop();
  int nlayer[60];
  for (int j = 0; j < 60; ++j) nlayer[j] = 0;

  for (vector<pointKey>::iterator iter = clusters.begin(); iter != clusters.end(); ++iter)
  {
    unsigned int layer = TrkrDefs::getLayer(iter->second);
    if (layer < (_nlayers_maps + _nlayers_intt + _nlayers_tpc*2./3.)) continue;
    ++nlayer[layer];
    t_fill->restart();
    _rtree.insert(*iter);
    t_fill->stop();
  }

  std::cout << "fill time: " << t_fill->get_accumulated_time() / 1000. << " sec" << std::endl;
}

pointKey PHCASeeding::toPointKey(coordKey v)
{
  return make_pair(point(v.first.at(0),v.first.at(1),v.first.at(2)),v.second);
}

vector<pointKey> PHCASeeding::toPointKey(vector<coordKey> v)
{
  vector<pointKey> output;
  output.reserve(v.size());
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
  output.reserve(p.size());
  for(vector<pointKey>::iterator pk = p.begin();pk != p.end();++pk)
  {
    output.push_back(make_pair(array<float,3>({pk->first.get<0>(),pk->first.get<1>(),pk->first.get<2>()}),pk->second));
  }
  return output;
}

int PHCASeeding::Process(PHCompositeNode *topNode)
{
  TFile fpara("CA_para.root", "RECREATE");
  TNtuple *NT = new TNtuple("NT", "NT", "pt:dpt:z:dz:phi:dphi:c:dc:nhit");
  _vertex = _vertex_map->get(0);

  PHTimer *t_seed = new PHTimer("t_seed");
  t_seed->stop();
  t_seed->restart();

  _rtree.clear();
  FillTree();
  t_seed->stop();
  cout << "Initial RTree fill time: " << t_seed->get_accumulated_time() / 1000 << " s" << endl;
  t_seed->restart();
  int numberofseeds = 0;
  //vector<coordKey> clusterCands = FindLinkedClusters(NT,t_seed);
  t_seed->stop();
  cout << "Linked cluster find time: " << t_seed->get_accumulated_time() / 1000 << " s" << endl;
  t_seed->restart();
  //_rtree.clear();
  //FillTree(toPointKey(clusterCands));
  t_seed->stop();
  cout << "Reduced RTree fill time: " << t_seed->get_accumulated_time() / 1000 << " s" << endl;
  t_seed->restart();
  numberofseeds += FindSeedsLayerSkip(NT,t_seed);
  t_seed->stop();
  cout << "number of seeds " << numberofseeds << endl;
  cout << "Kalman filtering time: " << t_seed->get_accumulated_time() / 1000 << " s" << endl;
  fpara.cd();
  NT->Write();
  fpara.Close();
  
  return Fun4AllReturnCodes::EVENT_OK;
}

vector<coordKey> PHCASeeding::FindLinkedClusters(TNtuple* NT, PHTimer* t_seed)
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
            _nlayers_maps+_nlayers_intt+2./3.*_nlayers_tpc-0.5, // layer 
            2*M_PI, // phi
            3, // eta
            _nlayers_maps+_nlayers_intt+_nlayers_tpc+0.5, // layer
            allClusters);
  t_seed->stop();
  cout << "allClusters search time: " << t_seed->get_accumulated_time() / 1000 << " s" << endl;
  LogDebug(" number of total clusters: " << allClusters.size() << endl);
  t_seed->restart();

  pair<vector<unordered_set<keylink>>,vector<unordered_set<keylink>>> links = CreateLinks(fromPointKey(allClusters),t_seed);
  cout << "created links\n";
  vector<vector<keylink>> bidirectionalLinks = FindBiLinks(links.first,links.second,t_seed);
  cout << "found bilinks\n";
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

int PHCASeeding::FindSeedsLayerSkip(TNtuple* NT,PHTimer* t_seed)
{
  vector<pointKey> reducedClusters;
  vector<unordered_set<keylink>> belowLinks;
  vector<unordered_set<keylink>> aboveLinks;
  belowLinks.resize(_nlayers_tpc);
  aboveLinks.resize(_nlayers_tpc);
  // get vector<pointKey> for all clusters in outer third of TPC
  QueryTree(_rtree,
            0, // phi
            -3, // eta
            _nlayers_maps+_nlayers_intt+2./3.*_nlayers_tpc-0.5, // layer 
            2*M_PI, // phi
            3, // eta
            _nlayers_maps+_nlayers_intt+_nlayers_tpc+0.5, // layer
            reducedClusters);
  t_seed->stop();
  cout << "allClusters search time: " << t_seed->get_accumulated_time() / 1000 << " s" << endl;
  LogDebug(" number of reduced clusters: " << reducedClusters.size() << endl);
  t_seed->restart();

  pair<vector<unordered_set<keylink>>,vector<unordered_set<keylink>>> links = CreateLinks(fromPointKey(reducedClusters),t_seed);
  vector<vector<keylink>> biLinks = FindBiLinks(links.first,links.second,t_seed);
  vector<keylist> trackSeedKeyLists = FollowBiLinks(biLinks,t_seed);
  int nseeds = ALICEKalmanFilter(trackSeedKeyLists,NT,t_seed);
  return nseeds;
}

pair<vector<unordered_set<keylink>>,vector<unordered_set<keylink>>> PHCASeeding::CreateLinks(vector<coordKey> clusters, PHTimer* t_seed, int mode)
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
    t_seed->stop();
    t_seed->restart();
    // get clusters near this one in adjacent layers
    double StartPhi = StartCluster->first[0];
    double StartEta = StartCluster->first[1];
    unsigned int StartLayer = StartCluster->first[2];
    //if(StartLayer > _start_layer) continue;
    TrkrCluster* StartCl = _cluster_map->findCluster(StartCluster->second);
    double StartX = StartCl->getPosition(0);
    double StartY = StartCl->getPosition(1);
    double StartZ = StartCl->getPosition(2);
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
              (double) StartLayer - 3.5,
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
              (double) StartLayer + 3.5,
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
        return array<double,3>{BelowCl->getPosition(0)-StartX,
          BelowCl->getPosition(1)-StartY,
          BelowCl->getPosition(2)-StartZ};});

    transform(ClustersAbove.begin(),ClustersAbove.end(),delta_above.begin(),
      [&](pointKey AboveCandidate){
        TrkrCluster* AboveCl = _cluster_map->findCluster(AboveCandidate.second);
        return array<double,3>{AboveCl->getPosition(0)-StartX,
          AboveCl->getPosition(1)-StartY,
          AboveCl->getPosition(2)-StartZ};});
    t_seed->stop();
    transform_time += t_seed->elapsed();
    t_seed->restart();

    // find the three clusters closest to a straight line
    // (by maximizing the cos of the angle between the (delta_eta,delta_phi) vectors)
    double maxCosPlaneAngle = 1;
    double minSumLengths = 1e9;
    coordKey bestBelowCluster = make_pair(array<float,3>({0.,0.,-1e9}),0);
    coordKey bestAboveCluster = make_pair(array<float,3>({0.,0.,-1e9}),0);
    for(size_t iAbove = 0; iAbove<delta_above.size(); ++iAbove)
    {
      for(size_t iBelow = 0; iBelow<delta_below.size(); ++iBelow)
      {
        double dotProduct = delta_below[iBelow][0]*delta_above[iAbove][0]+delta_below[iBelow][1]*delta_above[iAbove][1]+delta_below[iBelow][2]*delta_above[iAbove][2];
        double belowLength = sqrt(delta_below[iBelow][0]*delta_below[iBelow][0]+delta_below[iBelow][1]*delta_below[iBelow][1]+delta_below[iBelow][2]*delta_below[iBelow][2]);
        double aboveLength = sqrt(delta_above[iAbove][0]*delta_above[iAbove][0]+delta_above[iAbove][1]*delta_above[iAbove][1]+delta_above[iAbove][2]*delta_above[iAbove][2]);
        double cosPlaneAngle = dotProduct / (belowLength*aboveLength);
        if(cosPlaneAngle < _cosTheta_limit && (belowLength+aboveLength) < minSumLengths)
        {
          //maxCosPlaneAngle = cosPlaneAngle;
          minSumLengths = belowLength+aboveLength;
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
            return array<double,3>{BelowCl->getPosition(0)-StartX,
              BelowCl->getPosition(1)-StartY,
              BelowCl->getPosition(2)-StartZ};});
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
              return array<double,3>{AboveCl->getPosition(0)-StartX,
                AboveCl->getPosition(1)-StartY,
                AboveCl->getPosition(2)-StartZ};});
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
  cout << "triplet forming time: " << t_seed->get_accumulated_time() / 1000 << " s" << endl;
  cout << "starting cluster setup: " << cluster_find_time / 1000 << " s" << endl;
  cout << "RTree query: " << rtree_query_time /1000 << " s" << endl;
  cout << "Transform: " << transform_time /1000 << " s" << endl;
  cout << "Compute best triplet: " << compute_best_angle_time /1000 << " s" << endl;
  cout << "Set insert: " << set_insert_time /1000 << " s" << endl;
  t_seed->restart();

  return make_pair(belowLinks,aboveLinks);
}

vector<vector<keylink>> PHCASeeding::FindBiLinks(vector<unordered_set<keylink>> belowLinks,vector<unordered_set<keylink>> aboveLinks,PHTimer* t_seed)
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
  cout << "bidirectional link forming time: " << t_seed->get_accumulated_time() / 1000 << " s" << endl;
  t_seed->restart();

  return bidirectionalLinks;
}

vector<keylist> PHCASeeding::FollowBiLinks(vector<vector<keylink>> bidirectionalLinks,PHTimer* t_seed)
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
      unsigned int imax = 2;
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
  cout << "starting cluster finding time: " << t_seed->get_accumulated_time() / 1000 << " s" << endl;
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
  cout << "keychain assembly time: " << t_seed->get_accumulated_time() / 1000 << " s" << endl;
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
      TVector3 vec(cl->getPosition(0) - _vertex->get_x(), cl->getPosition(1) - _vertex->get_y(), cl->getPosition(2) - _vertex->get_z());

      double clus_phi = vec.Phi();
      clus_phi -= 2 * M_PI * floor(clus_phi / (2 * M_PI));
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
        " (x,y,z) = (" << cl->getPosition(0) << "," << cl->getPosition(1) << "," << cl->getPosition(2) << ")" << endl);
      lasteta = clus_eta;
      lastphi = clus_phi;
    }
  }
  LogDebug(" Total large jumps: " << jumpcount << endl);
  t_seed->stop();
  cout << "eta-phi sanity check time: " << t_seed->get_accumulated_time() / 1000 << " s" << endl;
  t_seed->restart();
  return trackSeedKeyLists;
}

int PHCASeeding::ALICEKalmanFilter(vector<keylist> trackSeedKeyLists,TNtuple* NT,PHTimer* t_seed)
{
  int nseeds = 0;
  for(vector<keylist>::iterator trackKeyChain = trackSeedKeyLists.begin(); trackKeyChain != trackSeedKeyLists.end(); ++trackKeyChain)
  {
    if(trackKeyChain->size() < 5) continue;
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
      float oldPhi = atan(y/x);
      LogDebug("old phi = " << oldPhi << endl);
      float alpha = newPhi - oldPhi;
      LogDebug("alpha = " << alpha << endl);
      if(!trackSeed.Rotate(alpha,trackLine,_max_sin_phi))
      {
        LogError("Rotate failed! Aborting for this seed...");
        break;
      }
      LogDebug("track coordinates (ALICE) after rotation: (" << trackSeed.GetX() << "," << trackSeed.GetY() << "," << trackSeed.GetZ() << ")" << endl);
      LogDebug("Transporting from " << alice_x << " to " << nextAlice_x << "...");
      if(!trackSeed.TransportToX(nextAlice_x,trackLine,_Bz,_max_sin_phi))
      {
        LogError("Transport failed! Aborting for this seed...");
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
    /*
    if(!trackSeed.CheckNumericalQuality())
    {
      cout << "ERROR: Track seed failed numerical quality check before conversion to sPHENIX coordinates! Skipping this one.\n";
      continue;
    } 
    */
    //    pt:z:dz:phi:dphi:c:dc
    // Fill NT with track parameters
    // float StartEta = -log(tan(atan(z0/sqrt(x0*x0+y0*y0))));
    float track_pt = fabs( 1./(trackSeed.GetQPt()));
    LogDebug("Track pt = " << track_pt << endl);
    float track_pterr = sqrt(trackSeed.GetErr2QPt())/(trackSeed.GetQPt()*trackSeed.GetQPt());
    LogDebug("Track pterr = " << track_pterr << endl);
    float track_z = trackSeed.GetZ();
    float track_zerr = sqrt(trackSeed.GetErr2Z());
    float track_phi = atan(trackCartesian_y/trackCartesian_x);
    float last_cluster_phierr = _cluster_map->findCluster(trackKeyChain->back())->getPhiError();
    // phi error assuming error in track radial coordinate is zero
    float track_phierr = sqrt(pow(last_cluster_phierr,2)+(pow(trackSeed.GetX(),2)*trackSeed.GetErr2Y()) / 
      pow(pow(trackSeed.GetX(),2)+pow(trackSeed.GetY(),2),2));
    LogDebug("Track phi = " << track_phi << endl);
    LogDebug("Track phierr = " << track_phierr << endl);
    float track_curvature = trackSeed.GetKappa(_Bz);
    float track_curverr = sqrt(trackSeed.GetErr2QPt())*_Bz;
    NT->Fill(track_pt, track_pterr, track_z, track_zerr, track_phi, track_phierr, track_curvature, track_curverr, trackKeyChain->size());
    SvtxTrack_v1 track;
    track.set_id(nseeds);
    for (unsigned int j = 0; j < trackKeyChain->size(); ++j)
    {
      track.insert_cluster_key(trackKeyChain->at(j));
    }
    track.set_ndf(trackSeed.GetNDF());
    if(trackSeed.GetQPt()<0) track.set_charge(-1);
    else track.set_charge(1);
    TrkrCluster *cl = _cluster_map->findCluster(trackKeyChain->at(0));
    track.set_x(cl->getX());  //track.set_x(cl->getX());
    track.set_y(cl->getY());  //track.set_y(cl->getY());
    track.set_z(cl->getZ());  //track.set_z(cl->getZ());
    track.set_px(track_pt * cos(track_phi) - track_pt*trackSeed.GetSinPhi()*sin(track_phi));
    track.set_py(track_pt * sin(track_phi) + track_pt*trackSeed.GetSinPhi()*cos(track_phi));
    track.set_pz(track_pt * trackSeed.GetDzDs());
    float s = sin(track_phi);
    float c = cos(track_phi);
    float q = trackSeed.GetQPt();
    float p = trackSeed.GetSinPhi();
    float d = trackSeed.GetDzDs();
    const float* cov = trackSeed.GetCov();
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
    track.set_error(2, 1, track.get_error(1, 2));
    track.set_error(3, 1, track.get_error(1, 3));
    track.set_error(4, 1, track.get_error(1, 4));
    track.set_error(3, 2, track.get_error(2, 3));
    track.set_error(4, 2, track.get_error(2, 4));
    track.set_error(4, 3, track.get_error(3, 4));
    if(!covIsPosDef(track))
    {
      repairCovariance(track);
    }
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

int PHCASeeding::Setup(PHCompositeNode *topNode)
{
  cout << "Called Setup" << endl;
  cout << "topNode:" << topNode << endl;
  PHTrackSeeding::Setup(topNode);
  InitializeGeometry(topNode);
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHCASeeding::End()
{
  cout << "Called End " << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}
