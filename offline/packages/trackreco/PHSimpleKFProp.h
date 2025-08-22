/*!
 *  \file PHSimpleKFProp.h
 *  \brief		kalman filter based propagator
 *  \author Michael Peters & Christof Roland
 */

#ifndef TRACKRECO_PHSIMPLEKFPROP_H
#define TRACKRECO_PHSIMPLEKFPROP_H

#include "ALICEKF.h"
#include "nanoflann.hpp"

// PHENIX includes
#include <tpc/TpcGlobalPositionWrapper.h>

#include <trackbase/TrkrDefs.h>

#include <fun4all/SubsysReco.h>

#include <phool/PHTimer.h>

#include <Acts/MagneticField/MagneticFieldProvider.hpp>

#include <Eigen/Core>

// STL includes
#include <limits>
#include <memory>
#include <string>
#include <vector>

class ActsGeometry;
class PHCompositeNode;
class TrkrClusterContainer;
class TrkrClusterIterationMapv1;
class SvtxTrackMap;
class TrackSeedContainer;
class TrackSeed;

using PositionMap = std::map<TrkrDefs::cluskey, Acts::Vector3>;

class PHSimpleKFProp : public SubsysReco
{
 public:
  PHSimpleKFProp(const std::string& name = "PHSimpleKFProp");
  ~PHSimpleKFProp() override = default;

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;

  // noop
  void set_field_dir(const double)
  {}

  // noop
  void magFieldFile(const std::string&)
  {}

  // noop
  void useConstBField(bool)
  {}

  // noop
  void setConstBField(float)
  {}

  void ghostRejection(bool set_value = true) { m_ghostrejection = set_value; }
  void set_max_window(double s) { _max_dist = s; }
  void useFixedClusterError(bool opt) { _use_fixed_clus_err = opt; }
  void setFixedClusterError(int i, double val) { _fixed_clus_err.at(i) = val; }
  void use_truth_clusters(bool truth)
  {
    _use_truth_clusters = truth;
  }
  void SetIteration(int iter) { _n_iteration = iter; }
  void set_pp_mode(bool mode) { _pp_mode = mode; }
  void set_max_seeds(unsigned int ui) { _max_seeds = ui; }
  enum class PropagationDirection
  {
    Outward,
    Inward
  };

  void setNeonFraction(double frac) { Ne_frac = frac; };
  void setArgonFraction(double frac) { Ar_frac = frac; };
  void setCF4Fraction(double frac) { CF4_frac = frac; };
  void setNitrogenFraction(double frac) { N2_frac = frac; };
  void setIsobutaneFraction(double frac) { isobutane_frac = frac; };

  void set_ghost_phi_cut(double d) { _ghost_phi_cut = d; }
  void set_ghost_eta_cut(double d) { _ghost_eta_cut = d; }
  void set_ghost_x_cut(double d) { _ghost_x_cut = d; }
  void set_ghost_y_cut(double d) { _ghost_y_cut = d; }
  void set_ghost_z_cut(double d) { _ghost_z_cut = d; }

  // number of threads
  void set_num_threads(int value) { m_num_threads = value; }

 private:
  bool _use_truth_clusters = false;
  bool m_ghostrejection = true;
  /// fetch node pointers
  int get_nodes(PHCompositeNode* topNode);
  std::vector<double> radii;
  std::vector<double> _vertex_x;
  std::vector<double> _vertex_y;
  std::vector<double> _vertex_z;
  std::vector<double> _vertex_xerr;
  std::vector<double> _vertex_yerr;
  std::vector<double> _vertex_zerr;
  std::vector<double> _vertex_ids;
  double _Bzconst = 10. * 0.000299792458f;
  // double _Bz = 1.4*_Bzconst;
  double _max_dist = .05;
  size_t _min_clusters_per_track = 3;

  //! internal field conversion to get the right particle sign

  double _max_sin_phi = 1.;
  bool _pp_mode = false;
  unsigned int _max_seeds = 0;

  TrkrClusterContainer* _cluster_map = nullptr;

  TrackSeedContainer* _track_map = nullptr;

  /// acts geometry
  ActsGeometry* m_tgeometry = nullptr;

  /// global position wrapper
  TpcGlobalPositionWrapper m_globalPositionWrapper;

  /// get global position for a given cluster
  /**
   * uses ActsTransformation to convert cluster local position into global coordinates
   * incorporates TPC distortion correction, if present
   */
  Acts::Vector3 getGlobalPosition(TrkrDefs::cluskey, TrkrCluster*) const;

  PositionMap PrepareKDTrees();

  bool TransportAndRotate(
    double old_layer,
    double new_layer,
    double& phi,
    GPUTPCTrackParam& kftrack,
    GPUTPCTrackParam::GPUTPCTrackFitParam& fp) const;

  bool PropagateStep(
    unsigned int& current_layer,
    double& current_phi,
    PropagationDirection& direction,
    std::vector<TrkrDefs::cluskey>& propagated_track,
    const std::vector<TrkrDefs::cluskey>& ckeys,
    GPUTPCTrackParam& kftrack,
    GPUTPCTrackParam::GPUTPCTrackFitParam& fp,
    const PositionMap& globalPositions) const;

  // TrackSeed objects store clusters in order of increasing cluster key (std::set<TrkrDefs::cluskey>),
  // which means we have to have a way to directly pass a list of clusters in order to extend looping tracks
  std::vector<TrkrDefs::cluskey> PropagateTrack(TrackSeed* track, PropagationDirection direction, GPUTPCTrackParam& aliceSeed, const PositionMap& globalPositions) const;
  std::vector<TrkrDefs::cluskey> PropagateTrack(TrackSeed* track, std::vector<TrkrDefs::cluskey>& ckeys, PropagationDirection direction, GPUTPCTrackParam& aliceSeed, const PositionMap& globalPositions) const;
  std::vector<std::vector<TrkrDefs::cluskey>> RemoveBadClusters(const std::vector<std::vector<TrkrDefs::cluskey>>& seeds, const PositionMap& globalPositions) const;

  template <typename T>
  struct KDPointCloud
  {
    KDPointCloud() = default;

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
    { return false; }

  };

  std::vector<std::shared_ptr<KDPointCloud<double>>> _ptclouds;

  std::vector<std::shared_ptr<nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, KDPointCloud<double>>, KDPointCloud<double>, 3>>> _kdtrees;

  std::unique_ptr<ALICEKF> fitter;

  void rejectAndPublishSeeds(std::vector<TrackSeed_v2>& seeds, const PositionMap& positions, std::vector<float>& trackChi2);

  void publishSeeds(const std::vector<TrackSeed_v2>&);

  int _max_propagation_steps = 200;

  //! fixed cluster errors
  bool _use_fixed_clus_err = false;

  //! fixed cluster errors
  std::array<double, 3> _fixed_clus_err = {.2, .2, .5};

  //! keep track of used clusters in previous iteration for iterative tracking
  TrkrClusterIterationMapv1* _iteration_map = nullptr;

  //! iteration number
  int _n_iteration = 0;

  //!@name TPC gas. Needed for Kalman Filter
  //@{
  double Ne_frac = 0.00;
  double Ar_frac = 0.75;
  double CF4_frac = 0.20;
  double N2_frac = 0.00;
  double isobutane_frac = 0.05;
  //@}

  //!@name ghost rejection cuts
  //@{
  double _ghost_phi_cut = std::numeric_limits<double>::max();
  double _ghost_eta_cut = std::numeric_limits<double>::max();
  double _ghost_x_cut = std::numeric_limits<double>::max();
  double _ghost_y_cut = std::numeric_limits<double>::max();
  double _ghost_z_cut = std::numeric_limits<double>::max();
  //@}

  //! number of threads
  /**
   * default is 0. This corresponds to allocating as many threads as available on the host
   * on CONDOR, by default, this results in only one thread allocated
   * being allocated unless several cores are allocated to the job
   * via request_cpus directive in the jdf
   */
  int m_num_threads = 0;

};

#endif
