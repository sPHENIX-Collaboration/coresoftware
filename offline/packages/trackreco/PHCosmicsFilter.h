/*!
 *  \file RPHCosmicsFilter.h
 *  \RTree based hough tracking for cosmics
 *  \author Christof Roland 
 */


//begin

#include <fun4all/SubsysReco.h>
#include <trackbase/TrkrDefs.h>  // for cluskey
#include <trackbase/ActsGeometry.h>


//TrkrCluster includes
#include <trackbase/TrkrCluster.h>                      // for TrkrCluster
#include <trackbase/TrkrDefs.h>                         // for getLayer, clu...
#include <trackbase/TrkrClusterContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>                          // for SubsysReco

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>                                // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>                              // for PHObject
#include <phool/PHRandomSeed.h>
#include <phool/PHTimer.h>                               // for PHTimer
#include <phool/getClass.h>
#include <phool/phool.h>                                 // for PHWHERE

//ROOT includes for debugging
#include <TFile.h>
#include <TMatrixDSymfwd.h>                              // for TMatrixDSym
#include <TMatrixTSym.h>                                 // for TMatrixTSym
#include <TMatrixTUtils.h>                               // for TMatrixTRow
#include <TNtuple.h>
#include <TVector3.h>                                    // for TVector3
#include <TVectorDfwd.h>                                 // for TVectorD
#include <TVectorT.h>                                    // for TVectorT

// gsl
#include <gsl/gsl_rng.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>

//BOOST for combi seeding
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>

// standard includes
#include <algorithm>
#include <assert.h>                                      // for assert
#include <cfloat>
#include <climits>                                      // for UINT_MAX
#include <cmath>
#include <cstdlib>                                      // for NULL, exit
#include <fstream>
#include <iostream>
#include <iterator>                                      // for back_insert_...
#include <memory>
#include <tuple>

class PHField;
class TGeoManager;

using namespace Eigen;

#define LogDebug(exp) std::cout << "DEBUG: " << __FILE__ << ": " << __LINE__ << ": " << exp
#define LogError(exp) std::cout << "ERROR: " << __FILE__ << ": " << __LINE__ << ": " << exp
#define LogWarning(exp) std::cout << "WARNING: " << __FILE__ << ": " << __LINE__ << ": " << exp


using namespace std;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
//end



#ifndef PHCOSMICSFILTER_H
#define PHCOSMICSFILTER_H

#include "PHTrackSeeding.h"


#if !defined(__CINT__) || defined(__CLING__)
//BOOST for combi seeding
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>

#include <boost/geometry/index/rtree.hpp>
#endif




// standard includes
#include <cfloat>
#include <iostream>                    // for operator<<, endl, basic_ostream
#include <map>
#include <string>                      // for string
#include <vector>

// forward declarations
class PHCompositeNode;
class PHG4CellContainer;
class PHG4CylinderGeomContainer;
class PHG4HitContainer;
class PHTimer;
class sPHENIXSeedFinder;
class SvtxClusterMap;
class SvtxCluster;
class SvtxTrackMap;
class SvtxTrack;
class SvtxVertexMap;
class SvtxVertex;
class TNtuple;
class TFile;
class TRKR_CLUSTER;
class SvtxHitMap;
class TrackSeedContainer;


namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
typedef bg::model::point<float, 3, bg::cs::cartesian> point;
typedef bg::model::box<point> box;
typedef std::pair<point, TrkrDefs::cluskey> pointKey;

typedef uint64_t cluskey;




class PHCosmicsFilter : public SubsysReco
{
 public:
  PHCosmicsFilter(const std::string &name = "PHRTreeSeeding");

  double chisq(const double *xx);
  //vector<TrkrCluster*> clusterpoints;
  virtual ~PHCosmicsFilter()
  {
  }
  void set_write_debug_ntuple(bool b){_write_ntp = b;}
  void set_create_tracks(bool b){_create_tracks = b;}
  void set_max_distance_to_origin(float val){ _max_dist_to_origin = val;}
  void set_min_nclusters(int n){ _min_nclusters = n;}
  
 protected:
  int Setup(PHCompositeNode *topNode);
  int GetNodes(PHCompositeNode* topNode);
  int Init(PHCompositeNode *topNode) override;
  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;


 private:
  /// fetch node pointers
  // int GetNodes(PHCompositeNode *topNode);

   //static vector<tuple<double, double, double>> clusterpoints;
   /*static*/ //vector<TrkrCluster*> clusterpoints;

  // node pointers
  TrkrClusterContainer *_cluster_map = nullptr;
  //nodes to get norm vector

  double phiadd(double phi1, double phi2);
  double phidiff(double phi1, double phi2);
  double pointKeyToTuple(pointKey *pK);
  double costfunction(const double *xx);
  //double chisq(const double *xx);
  void get_stub(const bgi::rtree<pointKey, bgi::quadratic<16>> &rtree, float pointx, float pointy, int &count, double &slope, double &intercept);
  ActsGeometry *tGeometry{nullptr};
#ifndef __CINT__
 private:
  int createNodes(PHCompositeNode *topNode);
  std::string m_trackMapName = "TpcTrackSeedContainer";
  TrackSeedContainer *m_seedContainer = nullptr;

  //int _nlayers_all;
  //unsigned int _nlayers_seeding;
  //std::vector<int> _seeding_layer;

  unsigned int _nevent = 0;
  bool _write_ntp = false;
  bool _create_tracks = false;
  float _max_dist_to_origin = 0;
  unsigned int _min_nclusters = 50;
  TNtuple *_ntp_cos = nullptr;
  TNtuple *_ntp_stub = nullptr;
  TNtuple *_ntp_max = nullptr;
  TFile *_tfile = nullptr;
  //std::vector<float> _radii_all;


#endif  // __CINT__
};

#endif
