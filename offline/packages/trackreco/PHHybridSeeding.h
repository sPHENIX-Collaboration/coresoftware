#ifndef TRACKRECO_PHHYBRIDSEEDING_H
#define TRACKRECO_PHHYBRIDSEEDING_H

/*!
 *  \file PHHybridSeeding.cc
 *  \brief Track Seeding using STAR "CA" algorithm and ALICE simplified Kalman filter
 *  \detail 
 *  \author Michael Peters & Christof Roland
 */
//begin

#include "PHTrackSeeding.h"      // for PHTrackSeeding

#include <trackbase/TrkrDefs.h>  // for cluskey
#include <trackbase/TrkrCluster.h>
#include <trackbase_historic/SvtxTrack_v1.h>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <cmath>     // for M_PI
#include <map>       // for map
#include <cstdint>  // for uint64_t
#include <string>    // for string
#include <utility>   // for pair
#include <vector>    // for vector


 
class PHCompositeNode;  // lines 196-196
class SvtxClusterMap;   // lines 202-202
class SvtxHitMap;       // lines 211-211
class SvtxTrackMap;     // lines 204-204
class SvtxVertex;
class SvtxVertexMap;    // lines 206-206

//end

typedef std::vector<TrkrDefs::cluskey> keylist;

class PHHybridSeeding : public PHTrackSeeding
{
 public:
  PHHybridSeeding(
      const std::string &name = "PHHybridSeeding",
      unsigned int min_clusters_per_track = 20,
      float cluster_z_error = 0.015,
      float cluster_alice_y_error = 0.015,
      float maxSinPhi = 0.999,
      float Bz = 14*0.000299792458f,
      float search_radius1 = 10.,
      float search_angle1 = M_PI/8.,
      size_t min_track_size1 = 10,
      float search_radius2 = 12.,
      float search_angle2 = M_PI/8.,
      size_t min_track_size2 = 6,
      size_t nthreads = 1
      );

  virtual ~PHHybridSeeding()
  {
  }
  
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

  std::vector<SvtxTrack_v1> ALICEKalmanFilter(std::vector<keylist> trackSeedKeyLists, bool use_nhits_limit = true);
  Eigen::Matrix<float,6,6> getEigenCov(SvtxTrack_v1 &track);
  bool covIsPosDef(SvtxTrack_v1 &track);
  void repairCovariance(SvtxTrack_v1 &track);
  void publishSeeds(std::vector<SvtxTrack_v1> seeds);
  bool checknan(float val, std::string name, int num);
  void setSearchRadius(float r1, float r2) {_search_radius1 = r1; _search_radius2 = r2;}
  void setSearchAngle(float a1, float a2) {_search_angle1 = a1; _search_angle2 = a2;}
  void setMinTrackSize(size_t n1, size_t n2) {_min_track_size1 = n1; _min_track_size2 = n2;}
  void setNThreads(size_t n) {_nthreads = n;}

 private:
  std::map<int, unsigned int> _layer_ilayer_map_all;
  std::map<int, unsigned int> _layer_ilayer_map;
  SvtxVertex *_vertex;

  unsigned int _min_clusters_per_track;
  float _cluster_z_error;
  float _cluster_alice_y_error;
  float _max_sin_phi;
  float _Bz;
  float _cosTheta_limit;
  float _search_radius1;
  float _search_angle1;
  size_t _min_track_size1;
  float _search_radius2;
  float _search_angle2;
  size_t _min_track_size2;
  size_t _nthreads;
  double _fieldDir = 1;
};

#endif
