#ifndef ALICEKF_H
#define ALICEKF_H

#include <trackbase_historic/SvtxTrack_v1.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrClusterContainer.h>
#include <phfield/PHField.h>
#include <phfield/PHFieldUtility.h>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <vector>
#include <string>
#include <utility>

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
  ~ALICEKF(){}
  std::vector<SvtxTrack_v1> ALICEKalmanFilter(std::vector<std::vector<TrkrDefs::cluskey>> chains, bool use_nhits_limit);
  Eigen::Matrix<double,6,6> getEigenCov(SvtxTrack_v1 &track);
  bool covIsPosDef(SvtxTrack_v1 &track);
  void repairCovariance(SvtxTrack_v1 &track);
  bool checknan(double val, std::string msg, int num);
  double get_Bz(double x, double y, double z);
  void CircleFitByTaubin(std::vector<std::pair<double,double>> pts, double &R, double &X0, double &Y0);
  private:
  PHField* _B;
  size_t _min_clusters_per_track = 20;
  TrkrClusterContainer* _cluster_map = nullptr;
  int Verbosity(){ return _v; }
  int _v = 0;
  double _Bzconst = 10*0.000299792458f;
  double _fieldDir = -1;
  double _max_sin_phi = 1.;
  bool _use_const_field = false;
};

#endif
