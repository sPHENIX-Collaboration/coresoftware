/*!
 *  \file PrelimDistortionCorrection.h
 *  \brief		Makes preliminary distortion corrections when crossing is unknown
 *  \author Tony Frawley
 */

#ifndef TRACKRECO_PRELIMDISTORTIONCORRECTION_H
#define TRACKRECO_PRELIMDISTORTIONCORRECTION_H

#include "ALICEKF.h"
#include "nanoflann.hpp"

// PHENIX includes
#include <fun4all/SubsysReco.h>
#include <tpc/TpcDistortionCorrection.h>
#include <trackbase/TrkrDefs.h>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>

#include <Eigen/Core>

// STL includes
#include <memory>
#include <string>
#include <vector>

class ActsGeometry;
class PHCompositeNode;
class PHField;
class TpcDistortionCorrectionContainer;
class TrkrClusterContainer;
class TrkrClusterIterationMapv1;
class SvtxTrackMap;
class TrackSeedContainer;
class TrackSeed;

using PositionMap = std::map<TrkrDefs::cluskey, Acts::Vector3>;

class PrelimDistortionCorrection : public SubsysReco
{
 public:
  PrelimDistortionCorrection(const std::string &name = "PrelimDistortionCorrection");
  ~PrelimDistortionCorrection() override = default;

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void set_field_dir(const double rescale)
  {
    _fieldDir = 1;
    if(rescale > 0)
      { _fieldDir = -1; }
  }
  void set_max_window(double s){_max_dist = s;}
  void useConstBField(bool opt){_use_const_field = opt;}
  void setConstBField(float b) { _const_field = b; }
  void useFixedClusterError(bool opt){_use_fixed_clus_err = opt;}
  void setFixedClusterError(int i, double val){_fixed_clus_err.at(i) = val;}
  void use_truth_clusters(bool truth)
  { _use_truth_clusters = truth; }
  void SetIteration(int iter){_n_iteration = iter;}
  void set_pp_mode(bool mode) {_pp_mode = mode;}
 
 private:

  /// tpc distortion correction utility class
  TpcDistortionCorrection m_distortionCorrection;

  bool _use_truth_clusters = false;
  
  /// fetch node pointers
  int get_nodes(PHCompositeNode *topNode);
  std::vector<double> radii;
  std::vector<double> _vertex_x;
  std::vector<double> _vertex_y;
  std::vector<double> _vertex_z;
  std::vector<double> _vertex_xerr;
  std::vector<double> _vertex_yerr;
  std::vector<double> _vertex_zerr;
  std::vector<double> _vertex_ids;
  //double _Bz = 1.4*_Bzconst;
  double _max_dist = .05;
  size_t _min_clusters_per_track = 3;
  double _fieldDir = -1;
  double _max_sin_phi = 1.;
  bool _pp_mode = false;

  TrkrClusterContainer *_cluster_map = nullptr;

  TrackSeedContainer *_track_map = nullptr;

  std::unique_ptr<PHField> _field_map = nullptr;
  
  /// acts geometry
  ActsGeometry *_tgeometry = nullptr;

  /// distortion correction container
  TpcDistortionCorrectionContainer* m_dcc = nullptr;

  /// get global position for a given cluster
  /**
   * uses ActsTransformation to convert cluster local position into global coordinates
   * incorporates TPC distortion correction, if present
   */
  Acts::Vector3 getGlobalPosition(TrkrDefs::cluskey, TrkrCluster*) const;

  PositionMap PrepareKDTrees();

  template <typename T>
  struct KDPointCloud
  {
    KDPointCloud<T>(){}
    std::vector<std::vector<T>> pts;
    inline size_t kdtree_get_point_count() const
    {
      return pts.size();
    }
    inline T kdtree_distance(const T* p1, const size_t idx_p2, size_t /*size*/) const
    {
      const T d0 = p1[0] - pts[idx_p2][0];
      const T d1 = p1[1] - pts[idx_p2][1];
      const T d2 = p1[2] - pts[idx_p2][2];
      return d0 * d0 + d1 * d1 + d2 * d2;
    }
    inline T kdtree_get_pt(const size_t idx, int dim) const
    {
      if (dim == 0)
        return pts[idx][0];
      else if (dim == 1)
        return pts[idx][1];
      else
        return pts[idx][2];
    }
    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /*bb*/) const
    {
      return false;
    }
  };
  std::vector<std::shared_ptr<KDPointCloud<double>>> _ptclouds;
  std::vector<std::shared_ptr<nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, KDPointCloud<double>>, KDPointCloud<double>,3>>> _kdtrees;
  std::unique_ptr<ALICEKF> fitter;
  double get_Bz(double x, double y, double z) const;
  void publishSeeds(std::vector<TrackSeed_v1>& seeds, PositionMap &positions);
  void publishSeeds(const std::vector<TrackSeed_v1>&);
//   void MoveToVertex();

  bool _use_const_field = false;
  float _const_field = 1.4;
  bool _use_fixed_clus_err = false;
  std::array<double,3> _fixed_clus_err = {.1,.1,.1};
  TrkrClusterIterationMapv1* _iteration_map = nullptr;
  int _n_iteration = 0;

};

#endif
