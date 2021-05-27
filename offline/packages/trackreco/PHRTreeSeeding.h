#ifndef TRACKRECO_PHRTREESEEDING_H
#define TRACKRECO_PHRTREESEEDING_H

/*!
 *  \file RTreeSeeding.h
 *  \brief Progressive pattern recgnition based on GenFit Kalman filter
 *  \detail using Alan Dion's HelixHough for seeding, GenFit Kalman filter to do track propagation
 *  \author Christof Roland & Haiwang Yu
 */

//begin

#include "PHTrackSeeding.h"      // for PHTrackSeeding

#include <trackbase/TrkrDefs.h>  // for cluskey

#include <boost/geometry/geometries/box.hpp>    // for box
#include <boost/geometry/geometries/point.hpp>  // for point
#include <boost/geometry/index/rtree.hpp>       // for rtree

#include <cstdint>  // for uint64_t
#include <map>       // for map
#include <string>    // for string
#include <utility>   // for pair
#include <vector>    // for vector

class PHCompositeNode;  // lines 196-196
class SvtxClusterMap;   // lines 202-202
class SvtxHitMap;       // lines 211-211
class SvtxTrackMap;     // lines 204-204
class SvtxVertexMap;    // lines 206-206

#define LogDebug(exp) std::cout << "DEBUG: " << __FILE__ << ": " << __LINE__ << ": " << exp
#define LogError(exp) std::cout << "ERROR: " << __FILE__ << ": " << __LINE__ << ": " << exp
#define LogWarning(exp) std::cout << "WARNING: " << __FILE__ << ": " << __LINE__ << ": " << exp

//#define _DEBUG_

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

typedef uint64_t cluskey;

class PHRTreeSeeding : public PHTrackSeeding
{
 public:
  PHRTreeSeeding(
      const std::string &name = "PHRTreeSeeding",
      unsigned int nlayers_maps = 3,
      unsigned int nlayers_intt = 4,
      unsigned int nlayers_tpc = 48,
      unsigned int start_layer = 53);

  double chisq(const double *xx);
  //vector<TrkrCluster*> clusterpoints;
  void set_phi_scale(float scale) { _phi_scale = scale; }
  void set_z_scale(float scale) { _z_scale = scale; }

  ~PHRTreeSeeding() override
  {
  }

 protected:
  int Setup(PHCompositeNode *topNode) override;
  int GetNodes(PHCompositeNode *topNode);
  int Process();
  int Process(PHCompositeNode *topNode) override;
  int InitializeGeometry(PHCompositeNode *topNode);

  int End() override;

 private:
  /// fetch node pointers

  // node pointers
  SvtxClusterMap *_g4clusters;
  SvtxTrackMap *_g4tracks;
  SvtxVertexMap *_g4vertexes;
  TrkrClusterContainer *_cluster_map;
  //nodes to get norm vector
  SvtxHitMap *_svtxhitsmap;
  int *_hit_used_map;
  int _hit_used_map_size;

  //seed searching parameters
  double phisr, etasr, phist, etast, phixt, etaxt;

  std::vector<float> _radii_all;

  double phiadd(double phi1, double phi2);
  double phidiff(double phi1, double phi2);
  double costfunction(const double *xx);
  void FillTree();

  double pointKeyToTuple(pointKey *pK);
  void QueryTree(const bgi::rtree<pointKey, bgi::quadratic<16>> &rtree, double phimin, double etamin, double lmin, double phimax, double etamax, double lmax, std::vector<pointKey> &returned_values);

 private:
  std::map<int, unsigned int> _layer_ilayer_map_all;
  std::map<int, unsigned int> _layer_ilayer_map;

  //int _nlayers_all;
  //unsigned int _nlayers_seeding;
  //std::vector<int> _seeding_layer;

  unsigned int _nlayers_maps;
  unsigned int _nlayers_intt;
  unsigned int _nlayers_tpc;
  unsigned int _start_layer;
  float _phi_scale;
  float _z_scale;
  //std::vector<float> _radii_all;

  bgi::rtree<pointKey, bgi::quadratic<16>> _rtree;
};

#endif
