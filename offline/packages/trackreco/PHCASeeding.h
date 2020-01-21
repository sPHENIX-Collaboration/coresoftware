#ifndef TRACKRECO_PHCASEEDING_H
#define TRACKRECO_PHCASEEDING_H

/*!
 *  \file CaSeeding.h
 *  \brief Progressive pattern recgnition based on GenFit Kalman filter
 *  \detail using Alan Dion's HelixHough for seeding, GenFit Kalman filter to do track propagation
 *  \author Christof Roland & Haiwang Yu
 */

//begin

#include <trackbase/TrkrDefs.h>  // for cluskey
#include "PHTrackSeeding.h"      // for PHTrackSeeding
#include <trackbase/TrkrCluster.h>

#if !defined(__CINT__) || defined(__CLING__)
#include <boost/geometry/geometries/box.hpp>    // for box
#include <boost/geometry/geometries/point.hpp>  // for point
#include <boost/geometry/index/rtree.hpp>       // for ca
#endif

#include <stdint.h>  // for uint64_t
#include <map>       // for map
#include <string>    // for string
#include <utility>   // for pair
#include <vector>    // for vector

class PHCompositeNode;  // lines 196-196
class SvtxClusterMap;   // lines 202-202
class SvtxHitMap;       // lines 211-211
class SvtxTrackMap;     // lines 204-204
class SvtxVertex;
class SvtxVertexMap;    // lines 206-206

//#define _DEBUG_

#if defined(_DEBUG_)
#define LogDebug(exp) std::cout << "DEBUG: " << __FILE__ << ": " << __LINE__ << ": " << exp
#else
#define LogDebug(exp) (void)0
#endif

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

#if !defined(__CINT__) || defined(__CLING__)
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
typedef bg::model::point<float, 3, bg::cs::cartesian> point;
typedef bg::model::box<point> box;
typedef std::pair<point, TrkrDefs::cluskey> pointKey;
#endif

typedef uint64_t cluskey;

class PHCASeeding : public PHTrackSeeding
{
 public:
  PHCASeeding(
      const std::string &name = "PHCASeeding",
      unsigned int nlayers_maps = 3,
      unsigned int nlayers_intt = 4,
      unsigned int nlayers_tpc = 48,
      unsigned int start_layer = 53);

#if !defined(__CINT__) || defined(__CLING__)
  double chisq(const double *xx);
#endif
  void set_phi_scale(float scale) { _phi_scale = scale; }
  void set_z_scale(float scale) { _z_scale = scale; }

  virtual ~PHCASeeding()
  {
  }

 protected:
  int Setup(PHCompositeNode *topNode);
  int GetNodes(PHCompositeNode *topNode);
  int Process();
  int Process(PHCompositeNode *topNode);
  int InitializeGeometry(PHCompositeNode *topNode);

  int End();

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
  void FillTree();

#if !defined(__CINT__) || defined(__CLING__)
  double pointKeyToTuple(pointKey *pK);
  void QueryTree(const bgi::rtree<pointKey, bgi::quadratic<16>> &rtree, double phimin, double etamin, double lmin, double phimax, double etamax, double lmax, std::vector<pointKey> &returned_values);
#endif

 private:
  std::map<int, unsigned int> _layer_ilayer_map_all;
  std::map<int, unsigned int> _layer_ilayer_map;

  //int _nlayers_all;
  //unsigned int _nlayers_seeding;
  //std::vector<int> _seeding_layer;
  SvtxVertex *_vertex;

  unsigned int _nlayers_maps;
  unsigned int _nlayers_intt;
  unsigned int _nlayers_tpc;
  unsigned int _start_layer;
  float _phi_scale;
  float _z_scale;
  //std::vector<float> _radii_all;

#if !defined(__CINT__) || defined(__CLING__)
  bgi::rtree<pointKey, bgi::quadratic<16>> _rtree;
#endif  // __CINT__
};

#endif
