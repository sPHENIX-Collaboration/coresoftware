// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHTPCTRACKSEEDVERTEXASSOC_H
#define PHTPCTRACKSEEDVERTEXASSOC_H

#include <trackreco/PHTrackPropagating.h>

#include <string>
#include <vector>

class PHCompositeNode;
class SvtxTrackMap;
class SvtxTrack;
class TrkrCluster;
class TF1;


class PHTpcTrackSeedVertexAssoc : public PHTrackPropagating
{
 public:

  PHTpcTrackSeedVertexAssoc(const std::string &name = "PHTpcTrackSeedVertexAssoc");

  virtual ~PHTpcTrackSeedVertexAssoc();

  void set_field_dir(const double rescale)
  {
    _fieldDir = -1;
    if(rescale > 0)
      _fieldDir = 1;     
  }
  void set_field(const std::string &field) { _field = field;}

  void reject_xy_outliers(const bool reject){_reject_xy_outliers = reject;}  
  void reject_z_outliers(const bool reject){_reject_z_outliers = reject;}

  void set_xy_residual_cut(const double cut){_xy_residual_cut = cut;}
  void set_z_residual_cut(const double cut){_z_residual_cut = cut;}
  void set_refit(const bool refit) {_refit = refit;}

 protected:
  int Setup(PHCompositeNode* topNode) override;

  int Process() override;

  int End() override;
  
 private:

  int GetNodes(PHCompositeNode* topNode);

  void  line_fit_clusters(std::vector<TrkrCluster*> clusters, double &a, double &b);
  void  line_fit(std::vector<std::pair<double,double>> points, double &a, double &b);
  void CircleFitByTaubin (std::vector<std::pair<double,double>> points, double &R, double &X0, double &Y0);
  std::vector<double> GetCircleClusterResiduals(std::vector<std::pair<double,double>> points, double R, double X0, double Y0);
  std::vector<double> GetLineClusterResiduals(std::vector<std::pair<double,double>> points, double A, double B);
  std::vector<TrkrCluster*> getTrackClusters(SvtxTrack *_tracklet_tpc);
						    
  std::string _track_map_name_silicon;

  
  SvtxTrack *_tracklet_tpc{nullptr};

  unsigned int _min_tpc_layer = 7;
  unsigned int _max_tpc_layer = 54;

  double _z_proj= 0;

  bool _reject_xy_outliers = false;
  bool _reject_z_outliers = false;
  bool _refit  = false;

  double _xy_residual_cut = 0.06;
  double _z_residual_cut = 0.15;


  std::string _field;
  int _fieldDir = -1;

};

#endif // PHTRUTHSILICONASSOCIATION_H
