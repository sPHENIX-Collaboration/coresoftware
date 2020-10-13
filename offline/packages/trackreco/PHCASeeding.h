#ifndef TRACKRECO_PHCASEEDING_H
#define TRACKRECO_PHCASEEDING_H

/*!
 *  \file PHCASeeding.cc
 *  \brief Track seeding using ALICE-style "cellular automaton" (CA) algorithm
 *  \detail 
 *  \author Michael Peters & Christof Roland
 */

//begin

#include "PHTrackSeeding.h"      // for PHTrackSeeding
#include <phool/PHTimer.h>

#include <trackbase/TrkrDefs.h>  // for cluskey
#include <trackbase/TrkrCluster.h>
#include <trackbase_historic/SvtxTrack_v1.h>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <TNtuple.h>

#include <boost/geometry/geometries/box.hpp>    // for box
#include <boost/geometry/geometries/point.hpp>  // for point
#include <boost/geometry/index/rtree.hpp>       // for ca

#include <cmath>     // for M_PI
#include <map>       // for map
#include <cstdint>  // for uint64_t
#include <string>    // for string
#include <utility>   // for pair
#include <vector>    // for vector
#include <set>


 
class PHCompositeNode;  // lines 196-196
class SvtxClusterMap;   // lines 202-202
class SvtxHitMap;       // lines 211-211
class SvtxTrackMap;     // lines 204-204
class SvtxVertex;
class SvtxVertexMap;    // lines 206-206

enum skip_layers {on, off};

//#define _USE_ALAN_FULL_VERTEXING_
#define _USE_ALAN_TRACK_REFITTING_

//#define _MEARGE_SEED_CLUSTER_
//#define _USE_ZERO_SEED_

//#define _USE_CONSTANT_SEARCH_WIN_

//#define _DO_FULL_FITTING_

//end

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
typedef bg::model::point<float, 3, bg::cs::cartesian> point;
typedef bg::model::box<point> box;
typedef std::pair<point, TrkrDefs::cluskey> pointKey;
typedef std::pair<std::array<float,3>, TrkrDefs::cluskey> coordKey;
typedef std::array<coordKey,2> keylink;
typedef std::vector<TrkrDefs::cluskey> keylist;


class PHCASeeding : public PHTrackSeeding
{
 public:
  PHCASeeding(
      const std::string &name = "PHCASeeding",
      unsigned int start_layer = 7,
      unsigned int end_layer = 55,
      unsigned int min_nhits_per_cluster = 2,
      unsigned int min_clusters_per_track = 20,
      const unsigned int nlayers_maps = 3,
      const unsigned int nlayers_intt = 4,
      const unsigned int nlayers_tpc = 48,
      float cluster_z_error = 0.015,
      float cluster_alice_y_error = 0.015,
      float neighbor_phi_width = .02,
      float neighbor_eta_width = .01,
      float maxSinPhi = 0.999,
      float Bz = 14*0.000299792458f,
      float cosTheta_limit = -0.8);

  virtual ~PHCASeeding()
  {
  }
  void SetLayerRange(unsigned int layer_low, unsigned int layer_up) {_start_layer = layer_low; _end_layer = layer_up;}
  void SetSearchWindow(float eta_width, float phi_width) {_neighbor_eta_width = eta_width; _neighbor_phi_width = phi_width;}
  void SetMinHitsPerCluster(int minHits) {_min_nhits_per_cluster = minHits;}
  void SetMinClustersPerTrack(int minClus) {_min_clusters_per_track = minClus;}

  void set_field_dir(const double rescale)
  {
    _fieldDir = -1;
    if(rescale > 0)
      _fieldDir = 1;     
  }

 protected:
  virtual int Setup(PHCompositeNode *topNode);
  virtual int Process(PHCompositeNode *topNode);
  int InitializeGeometry(PHCompositeNode *topNode);
  int FindSeedsLayerSkip(double cosTheta_limit,TNtuple* NT,PHTimer* t);
  virtual int End();

 private:
  /// fetch node pointers

  // node pointers
  SvtxTrackMap *_g4tracks;
  SvtxVertexMap *_g4vertexes;
  //nodes to get norm vector
  SvtxHitMap *_svtxhitsmap;
  int *_hit_used_map;
  int _hit_used_map_size;

  std::vector<float> _radii_all;

  double phiadd(double phi1, double phi2);
  double phidiff(double phi1, double phi2);
  void FillTree();
  void FillTree(std::vector<pointKey> clusters);
  std::vector<coordKey> FindLinkedClusters(TNtuple* NT, PHTimer* t_seed);
  int FindSeedsWithMerger(TNtuple* NT, PHTimer* t_seed);
  std::pair<std::vector<std::unordered_set<keylink>>,std::vector<std::unordered_set<keylink>>> CreateLinks(std::vector<coordKey> clusters, PHTimer* t_seed, int mode = skip_layers::off);
  std::vector<std::vector<keylink>> FindBiLinks(std::vector<std::unordered_set<keylink>> belowLinks, std::vector<std::unordered_set<keylink>> aboveLinks, PHTimer* t_seed);
  std::vector<keylist> FollowBiLinks(std::vector<std::vector<keylink>> bidirectionalLinks, PHTimer* t_seed);
  int ALICEKalmanFilter(std::vector<keylist> trackSeedKeyLists, TNtuple* NT, PHTimer* t_seed);
  void QueryTree(const bgi::rtree<pointKey, bgi::quadratic<16>> &rtree, double phimin, double etamin, double lmin, double phimax, double etamax, double lmax, std::vector<pointKey> &returned_values);
  pointKey toPointKey(coordKey v);
  std::vector<pointKey> toPointKey(std::vector<coordKey> v);
  coordKey fromPointKey(pointKey p);
  std::vector<coordKey> fromPointKey(std::vector<pointKey> p);
  Eigen::Matrix<float,6,6> getEigenCov(SvtxTrack_v1 &track);
  bool covIsPosDef(SvtxTrack_v1 &track);
  void repairCovariance(SvtxTrack_v1 &track);
  std::vector<keylist> MergeSeeds(std::vector<keylist> seeds, PHTimer* t_seed);
  pointKey makepointKey(TrkrDefs::cluskey k);


 private:
  std::map<int, unsigned int> _layer_ilayer_map_all;
  std::map<int, unsigned int> _layer_ilayer_map;

  //int _nlayers_all;
  //unsigned int _nlayers_seeding;
  //std::vector<int> _seeding_layer;
  SvtxVertex *_vertex;

  const unsigned int _nlayers_maps;
  const unsigned int _nlayers_intt;
  const unsigned int _nlayers_tpc;
  unsigned int _start_layer;
  unsigned int _end_layer;
  unsigned int _min_nhits_per_cluster;
  unsigned int _min_clusters_per_track;
  float _cluster_z_error;
  float _cluster_alice_y_error;
  float _neighbor_phi_width;
  float _neighbor_eta_width;
  float _max_sin_phi;
  float _Bz;
  float _cosTheta_limit;
  //std::vector<float> _radii_all;
  double _fieldDir = 1;

  bgi::rtree<pointKey, bgi::quadratic<16>> _rtree;
};

#endif
