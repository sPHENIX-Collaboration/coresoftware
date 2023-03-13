#ifndef CALORECO_RAWCLUSTERBUILDERFWD_H
#define CALORECO_RAWCLUSTERBUILDERFWD_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class RawCluster;
class RawClusterContainer;
class RawTowerContainer;
class RawTowerGeomContainer;

class RawClusterBuilderFwd : public SubsysReco
{
 public:
  explicit RawClusterBuilderFwd(const std::string &name = "RawClusterBuilder");
  ~RawClusterBuilderFwd() override {}

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;
  void Detector(const std::string &d) { detector = d; }

  void set_threshold_energy(const float e) { _min_tower_e = e; }
  void checkenergy(const int i = 1) { chkenergyconservation = i; }

 private:
  void CreateNodes(PHCompositeNode *topNode);
  bool CorrectPhi(RawCluster *cluster, RawTowerContainer *towers, RawTowerGeomContainer *towergemom);

  RawClusterContainer *_clusters;

  float _min_tower_e;
  int chkenergyconservation;

  std::string detector;
  std::string ClusterNodeName;
};

#endif
