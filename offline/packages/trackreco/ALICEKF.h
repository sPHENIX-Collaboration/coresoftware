#ifndef ALICEKF_H
#define ALICEKF_H

#include <phfield/PHField.h>
#include <trackbase/ClusterErrorPara.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase_historic/TrackSeed_v2.h>
#include "GPUTPCTrackParam.h"

#include <Acts/Definitions/Algebra.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <string>
#include <utility>
#include <vector>

using PositionMap = std::map<TrkrDefs::cluskey, Acts::Vector3>;
using TrackSeedAliceSeedMap = std::pair<std::vector<TrackSeed_v2>, std::vector<GPUTPCTrackParam>>;

class ALICEKF
{
 public:
  ALICEKF( TrkrClusterContainer* cmap, PHField* B, unsigned int min_clusters, double max_sin_phi, int verbosity);

  ~ALICEKF() = default;

  explicit ALICEKF(const ALICEKF&) = delete;
  ALICEKF& operator=(const ALICEKF&) = delete;

  void setNeonFraction(double frac) { Ne_frac = frac; };
  void setArgonFraction(double frac) { Ar_frac = frac; };
  void setCF4Fraction(double frac) { CF4_frac = frac; };
  void setNitrogenFraction(double frac) { N2_frac = frac; };
  void setIsobutaneFraction(double frac) { isobutane_frac = frac; };

  bool TransportAndRotate(double old_radius, double new_radius, double& phi, GPUTPCTrackParam& kftrack, GPUTPCTrackParam::GPUTPCTrackFitParam& fp) const;
  bool FilterStep(TrkrDefs::cluskey ckey, std::vector<TrkrDefs::cluskey>& keys, double& current_phi, GPUTPCTrackParam& kftrack, GPUTPCTrackParam::GPUTPCTrackFitParam& fp, const PositionMap& globalPositions) const;

  TrackSeedAliceSeedMap ALICEKalmanFilter(const std::vector<std::vector<TrkrDefs::cluskey>>& chains, bool use_nhits_limit, const PositionMap& globalPositions, std::vector<float>& trackChi2) const;
  bool covIsPosDef(Eigen::Matrix<double, 6, 6>& cov) const;
  void repairCovariance(Eigen::Matrix<double, 6, 6>& cov) const;
  bool checknan(double val, const std::string& msg, int num) const;
  double get_Bz(double x, double y, double z) const;

  //! used for fast momentum calculation
  void setConstBField(float b) { _const_field = b; }

void useFixedClusterError(bool opt) { _use_fixed_clus_error = opt; }
  void setFixedClusterError(int i, double val) { _fixed_clus_error.at(i) = val; }
  double getClusterError(TrkrCluster* c, TrkrDefs::cluskey key, Acts::Vector3 global, int i, int j) const;
  std::vector<double> GetCircleClusterResiduals(const std::vector<std::pair<double, double>>& pts, double R, double X0, double Y0) const;
  std::vector<double> GetLineClusterResiduals(const std::vector<std::pair<double, double>>& pts, double A, double B) const;
  double get_Bzconst() const { return _Bzconst; }

 private:
  int Verbosity() const
  { return _v; }

  PHField* _B = nullptr;
  size_t _min_clusters_per_track = 20;
  TrkrClusterContainer* _cluster_map = nullptr;

  int _v = 0;
  static constexpr double _Bzconst = 10. * 0.000299792458f;
  double _max_sin_phi = 1.;

  //! used for fast momentum calculation
  float _const_field = 1.4;

  // cluster error parametrization
  std::unique_ptr<ClusterErrorPara> _ClusErrPara;

  // cluster errors
  bool _use_fixed_clus_error = true;
  std::array<double, 3> _fixed_clus_error = {.2, .2, .5};

  double Ne_frac = 0.00;
  double Ar_frac = 0.75;
  double CF4_frac = 0.20;
  double N2_frac = 0.00;
  double isobutane_frac = 0.05;
};

#endif
