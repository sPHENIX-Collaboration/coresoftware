#ifndef RAWCLUSTERBUILDERV1_H__
#define RAWCLUSTERBUILDERV1_H__

#include <fun4all/SubsysReco.h>
#include <string>

class PHCompositeNode;
class RawCluster;
class RawClusterContainer;
class RawTowerContainer;
class RawTowerGeomContainer;
class BEmcRec;

class RawClusterBuilderv1 : public SubsysReco {

 public:
  RawClusterBuilderv1(const std::string& name = "RawClusterBuilder"); 
  virtual ~RawClusterBuilderv1();

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

  BEmcRec* bemc;
  float fEnergyNorm;

  float _min_tower_e;
  int chkenergyconservation;

  std::string detector;
  std::string ClusterNodeName;

};

#endif /* RAWCLUSTERBUILDERV1_H__ */
