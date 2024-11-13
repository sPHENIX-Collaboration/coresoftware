#ifndef TRACKRECO_PHHYBRIDSEEDING_H
#define TRACKRECO_PHHYBRIDSEEDING_H

/*!
 *  \file PHHybridSeeding.cc
 *  \brief Track Seeding using STAR "CA" algorithm and ALICE simplified Kalman filter
 *  \detail 
 *  \author Michael Peters
 */
//begin

#include "PHTrackSeeding.h"      // for PHTrackSeeding
#include "ALICEKF.h"

#include <trackbase/TrkrDefs.h>  // for cluskey
#include <trackbase/TrkrCluster.h>
#include <trackbase_historic/SvtxTrack_v3.h>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <cmath>     // for M_PI
#include <map>       // for map
#include <cstdint>  // for uint64_t
#include <string>    // for string
#include <utility>   // for pair
#include <vector>    // for vector
#include <memory>

 
class PHCompositeNode;  // lines 196-196
class SvtxClusterMap;   // lines 202-202
class SvtxHitMap;       // lines 211-211
class SvtxTrackMap;     // lines 204-204
class SvtxVertex;
class SvtxVertexMap;    // lines 206-206

//end

class PHHybridSeeding : public PHTrackSeeding
{
 public:
  PHHybridSeeding(
      const std::string &name = "PHHybridSeeding",
      double max_sin_phi = 1000.,
      double fieldDir = 1,
      double search_radius1 = 3.,
      double search_angle1 = M_PI/8.,
      size_t min_track_size1 = 10,
      double search_radius2 = 6.,
      double search_angle2 = M_PI/8.,
      size_t min_track_size2 = 5,
      size_t nthreads = 1
      );

  ~PHHybridSeeding() override
  {
  }
  
  void set_field_dir(const double rescale)
  {

    if(rescale > 0)
      _fieldDir = -1;
    else
    {
      _fieldDir = 1;
      //_Bz = -1*_Bz;
    }
  
  }
  void setSearchRadius(double r1, double r2) {_search_radius1 = r1; _search_radius2 = r2;}
  void setSearchAngle(double a1, double a2) {_search_angle1 = a1; _search_angle2 = a2;}
  void setMinTrackSize(size_t n1, size_t n2) {_min_track_size1 = n1; _min_track_size2 = n2;}
  void setNThreads(size_t n) {_nthreads = n;} 
  void setMinFitTrackSize(size_t s) {_min_fit_track_size = s;}

 protected:
  int Setup(PHCompositeNode *topNode) override;
  int Process(PHCompositeNode *topNode) override;
  int InitializeGeometry(PHCompositeNode *topNode);
  int End() override;

 private:
  /// fetch node pointers

  // node pointers
  //nodes to get norm vector

  std::vector<double> _radii_all;

  std::vector<double> _vertex_ids;
  std::vector<double> _vertex_x;
  std::vector<double> _vertex_y;
  std::vector<double> _vertex_z;
  std::vector<double> _vertex_xerr;
  std::vector<double> _vertex_yerr;
  std::vector<double> _vertex_zerr;

  std::map<int, unsigned int> _layer_ilayer_map_all;
  std::map<int, unsigned int> _layer_ilayer_map;

  void publishSeeds(std::vector<SvtxTrack_v3> seeds);

  double _max_sin_phi;
  double _fieldDir;
  double _search_radius1;
  double _search_angle1;
  size_t _min_track_size1;
  double _search_radius2;
  double _search_angle2;
  size_t _min_track_size2;
  size_t _nthreads;
  size_t _min_fit_track_size = 5;
  std::shared_ptr<ALICEKF> fitter;
};

#endif
