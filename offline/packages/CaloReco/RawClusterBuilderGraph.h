#ifndef CALORECO_RAWCLUSTERBUILDER_H
#define CALORECO_RAWCLUSTERBUILDER_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class RawClusterContainer;

class RawClusterBuilderGraph : public SubsysReco
{
 public:
  explicit RawClusterBuilderGraph(const std::string &name = "RawClusterBuilderGraph");
  ~RawClusterBuilderGraph() override {}

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;
  void Detector(const std::string &d) { detector = d; }

  void set_threshold_energy(const float e) { _min_tower_e = e; }
  void checkenergy(const int i = 1) { chkenergyconservation = i; }

 private:
  void CreateNodes(PHCompositeNode *topNode);

  RawClusterContainer *_clusters;

  float _min_tower_e;
  int chkenergyconservation;

  std::string detector;
  std::string ClusterNodeName;
};

#endif
