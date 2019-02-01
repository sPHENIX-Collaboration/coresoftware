#ifndef CALORECO_RAWCLUSTERBUILDERTEMPLATE_H
#define CALORECO_RAWCLUSTERBUILDERTEMPLATE_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class RawCluster;
class RawClusterContainer;
class RawTowerContainer;
class RawTowerGeomContainer;
class BEmcRecCEMC;

class RawClusterBuilderTemplate : public SubsysReco {

 public:
  RawClusterBuilderTemplate(const std::string& name = "RawClusterBuilderGraph");
  virtual ~RawClusterBuilderTemplate();

  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
  void Detector(const std::string &d) {detector = d;}

  void set_threshold_energy(const float e) {_min_tower_e = e;}
  void setEnergyNorm(float norm) {fEnergyNorm = norm;}
  void checkenergy(const int i = 1) {chkenergyconservation = i;}

 private:
  void CreateNodes(PHCompositeNode *topNode);
  bool CorrectPhi(RawCluster* cluster, RawTowerContainer* towers, RawTowerGeomContainer *towergemom);
  bool Cell2Abs(RawTowerGeomContainer *towergeom, float phiC, float etaC, float& phi, float& eta);

  RawClusterContainer* _clusters;

  BEmcRecCEMC* bemc;
  float fEnergyNorm;

  float _min_tower_e;
  int chkenergyconservation;

  std::string detector;
  std::string ClusterNodeName;

};

#endif /* RawClusterBuilderTemplate_H__ */
