#ifndef RAWCLUSTERBUILDERFWD_H__
#define RAWCLUSTERBUILDERFWD_H__

#include <fun4all/SubsysReco.h>
#include <string>

class PHCompositeNode;
class RawCluster;
class RawClusterContainer;
class RawTowerContainer;
class RawTowerGeomContainer;

class RawClusterBuilderFwd : public SubsysReco {

 public:
  RawClusterBuilderFwd(const std::string& name = "RawClusterBuilder"); 
  virtual ~RawClusterBuilderFwd() {}

  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
  void Detector(const std::string &d) {detector = d;}

  void set_threshold_energy(const float e) {_min_tower_e = e;}
  void checkenergy(const int i = 1) {chkenergyconservation = i;}

 private:
  void CreateNodes(PHCompositeNode *topNode);
  bool CorrectPhi(RawCluster* cluster, RawTowerContainer* towers, RawTowerGeomContainer *towergemom);

  RawClusterContainer* _clusters;

  float _min_tower_e;
  int chkenergyconservation;

  std::string detector;
  std::string ClusterNodeName;

};

#endif /* RAWCLUSTERBUILDERFWD_H__ */
