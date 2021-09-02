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

  ~PHTpcTrackSeedVertexAssoc() override;

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
  void findRoot(const double R, const double X0, const double Y0, double& x, double& y);
						    
  std::string _track_map_name_silicon;

  
  SvtxTrack *_tracklet_tpc{nullptr};

  unsigned int _min_tpc_layer = 7;
  unsigned int _max_tpc_layer = 54;

  double _z_proj= 0;

};

#endif // PHTRACKSEEDVERTEXASSOCIATION_H
