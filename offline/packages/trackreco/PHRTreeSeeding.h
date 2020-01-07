
/*!
 *  \file RTreeSeeding.h
 *  \brief Progressive pattern recgnition based on GenFit Kalman filter
 *  \detail using Alan Dion's HelixHough for seeding, GenFit Kalman filter to do track propagation
 *  \author Christof Roland & Haiwang Yu
 */


//begin

//#include "PHG4KalmanPatRec.h"

// trackbase_historic includes
#include <trackbase_historic/SvtxCluster.h>
#include <trackbase_historic/SvtxClusterMap.h>
#include <trackbase_historic/SvtxHitMap.h>
#include <trackbase_historic/SvtxHit.h>                  // for SvtxHit
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v1.h>
#include <trackbase_historic/SvtxTrackState.h>
#include <trackbase_historic/SvtxTrack_v1.h>
#include <trackbase_historic/SvtxVertex.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase_historic/SvtxVertexMap_v1.h>
#include <trackbase_historic/SvtxVertex_v1.h>

// sPHENIX Geant4 includes
#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>

#include <intt/CylinderGeomIntt.h>

#include <mvtx/CylinderGeom_Mvtx.h>

#include <g4bbc/BbcVertex.h>
#include <g4bbc/BbcVertexMap.h>

// sPHENIX includes

#include <phfield/PHFieldUtility.h>

#include <phgeom/PHGeomUtility.h>

//FIXME remove includes below after having real vertxing
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>



//TrkrCluster includes
#include <trackbase/TrkrCluster.h>                      // for TrkrCluster
#include <trackbase/TrkrDefs.h>                         // for getLayer, clu...
#include <trackbase/TrkrClusterContainer.h>


// Helix Hough includes
#include <HelixHough/HelixKalmanState.h>                 // for HelixKalmanS...
#include <HelixHough/HelixRange.h>
#include <HelixHough/SimpleHit3D.h>
#include <HelixHough/SimpleTrack3D.h>
#include <HelixHough/sPHENIXSeedFinder.h>                // for sPHENIXSeedF...
#include <HelixHough/VertexFinder.h>

#include <phgenfit/Fitter.h>
#include <phgenfit/Measurement.h>                        // for Measurement
#include <phgenfit/PlanarMeasurement.h>
#include <phgenfit/SpacepointMeasurement.h>
#include <phgenfit/Track.h>

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

// GenFit
#include <GenFit/EventDisplay.h>                         // for EventDisplay
#include <GenFit/MeasuredStateOnPlane.h>
#include <GenFit/RKTrackRep.h>
#include <GenFit/Track.h>

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
namespace genfit { class AbsTrackRep; }

using namespace Eigen;

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

using namespace std;
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
//end



#ifndef TRACKRECO_PHRTREESEEDING_H
#define TRACKRECO_PHRTREESEEDING_H

#include "PHTrackSeeding.h"

// Helix Hough includes
#if !defined(__CINT__) || defined(__CLING__)
#include <HelixHough/SimpleHit3D.h>
#include <HelixHough/SimpleTrack3D.h>
#include <HelixHough/VertexFinder.h>
#include <HelixHough/sPHENIXSeedFinder.h>
#include <Eigen/Core>                  // for Matrix
#endif


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

  virtual ~PHRTreeSeeding()
  {
  }

 protected:
  int Setup(PHCompositeNode *topNode);
  int GetNodes(PHCompositeNode* topNode);
  int Process();
  int Process(PHCompositeNode *topNode);
  int InitializeGeometry(PHCompositeNode *topNode);

  int End();
  
 private:
  /// fetch node pointers
  // int GetNodes(PHCompositeNode *topNode);

   //static vector<tuple<double, double, double>> clusterpoints;
   /*static*/ //vector<TrkrCluster*> clusterpoints;

  // node pointers
  SvtxClusterMap* _g4clusters;
  SvtxTrackMap* _g4tracks;
  SvtxVertexMap* _g4vertexes;
  TrkrClusterContainer *_cluster_map;
  //nodes to get norm vector
  SvtxHitMap* _svtxhitsmap;
  int* _hit_used_map;
  int _hit_used_map_size;


  //seed searching parameters
  double phisr,etasr,phist,etast,phixt,etaxt;

  std::vector<float> _radii_all;

  double phiadd(double phi1, double phi2);
  double phidiff(double phi1, double phi2);
  double pointKeyToTuple(pointKey *pK);
  double costfunction(const double *xx);
  //double chisq(const double *xx);
  void FillTree();
  void QueryTree(const bgi::rtree<pointKey, bgi::quadratic<16>> &rtree, double phimin, double etamin, double lmin, double phimax, double etamax, double lmax,std::vector<pointKey> &returned_values);


#ifndef __CINT__
 private:

  std::map<int, unsigned int> _layer_ilayer_map_all;
  std::map<int, unsigned int> _layer_ilayer_map;
  bgi::rtree<pointKey, bgi::quadratic<16> > _rtree;


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


#endif  // __CINT__
};

#endif
