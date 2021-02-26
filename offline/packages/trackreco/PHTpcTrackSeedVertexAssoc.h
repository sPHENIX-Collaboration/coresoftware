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

 protected:
  int Setup(PHCompositeNode* topNode) override;

  int Process() override;

  int End() override;
  
 private:

  int GetNodes(PHCompositeNode* topNode);

  void  line_fit_clusters(std::vector<TrkrCluster*> clusters, double &a, double &b);
  void  line_fit(std::vector<std::pair<double,double>> points, double &a, double &b);
  void CircleFitByTaubin (std::vector<std::pair<double,double>> points, double &R, double &X0, double &Y0);
						    
  std::string _track_map_name_silicon;

  
  SvtxTrackMap *_track_map_silicon{nullptr};
  SvtxTrack *_tracklet_tpc{nullptr};

  unsigned int _min_tpc_layer = 7;
  unsigned int _max_tpc_layer = 54;

  double _z_proj= 0;

  std::string _field;
  int _fieldDir = -1;

};

#endif // PHTRUTHSILICONASSOCIATION_H
