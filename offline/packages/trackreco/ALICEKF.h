#ifndef ALICEKF_H
#define ALICEKF_H

#include <trackbase_historic/SvtxTrack_v3.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrCluster.h>
#include <phfield/PHField.h>
#include <phfield/PHFieldUtility.h>

#include <Acts/Definitions/Algebra.hpp>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <string>
#include <utility>

using PositionMap = std::map<TrkrDefs::cluskey, Acts::Vector3>;

class ALICEKF
{
  public:
  ALICEKF(PHCompositeNode *topNode, 
          TrkrClusterContainer* cmap, 
          double fieldDir,
          unsigned int min_clusters,
          double max_sin_phi,
          int verbosity)
  {
    _B = PHFieldUtility::GetFieldMapNode(nullptr,topNode);
    _cluster_map = cmap;
    _fieldDir = fieldDir;
    _max_sin_phi = max_sin_phi;
    _v = verbosity;
    _min_clusters_per_track = min_clusters;
  }
  std::vector<SvtxTrack_v3> ALICEKalmanFilter(const std::vector<std::vector<TrkrDefs::cluskey>>& chains, bool use_nhits_limit, const PositionMap& globalPositions) const;
  Eigen::Matrix<double,6,6> getEigenCov(const SvtxTrack_v3 &track) const;
  bool covIsPosDef(const SvtxTrack_v3 &track) const;
  void repairCovariance(SvtxTrack_v3 &track) const;
  bool checknan(double val, const std::string &msg, int num) const;
  double get_Bz(double x, double y, double z) const;
  void CircleFitByTaubin(const std::vector<std::pair<double,double>>& pts, double &R, double &X0, double &Y0) const;
  void useConstBField(bool opt) {_use_const_field = opt;}
  void useFixedClusterError(bool opt) {_use_fixed_clus_error = opt;}
  void setFixedClusterError(int i,double val) {_fixed_clus_error.at(i)=val;}
  double getClusterError(TrkrCluster* c, Acts::Vector3 global, int i, int j) const;
  void line_fit(const std::vector<std::pair<double,double>>& pts, double& a, double& b) const;
  std::vector<double> GetCircleClusterResiduals(const std::vector<std::pair<double,double>>& pts, double R, double X0, double Y0) const;
  std::vector<double> GetLineClusterResiduals(const std::vector<std::pair<double,double>>& pts, double A, double B) const; 
  private:
  PHField* _B = nullptr;
  size_t _min_clusters_per_track = 20;
  TrkrClusterContainer* _cluster_map = nullptr;
  int Verbosity() const
  { return _v; }
  
  int _v = 0;
  double _Bzconst = 10*0.000299792458f;
  double _fieldDir = -1;
  double _max_sin_phi = 1.;
  bool _use_const_field = false;
  bool _use_fixed_clus_error = true;
  std::array<double,3> _fixed_clus_error = {.1,.1,.1};
};

#endif
