#ifndef CALORECO_RAWCLUSTERBUILDERTEMPLATE_H
#define CALORECO_RAWCLUSTERBUILDERTEMPLATE_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class RawClusterContainer;
class RawTowerGeomContainer;
class BEmcRec;

class RawClusterBuilderTemplate : public SubsysReco
{
 public:
  RawClusterBuilderTemplate(const std::string& name = "RawClusterBuilderTemplate");
  ~RawClusterBuilderTemplate() override;

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  void Detector(const std::string& d);

  void SetCylindricalGeometry();
  void SetPlanarGeometry();
  void PrintGeometry() { bPrintGeom = true; }  // Prints it at InitRun time
  void PrintCylGeom(RawTowerGeomContainer* towergeom, const std::string& fname);
  void SetProfileProb(bool pprob) { bProfProb = pprob; }

  void set_threshold_energy(const float e) { _min_tower_e = e; }
  void setEnergyNorm(float norm) { fEnergyNorm = norm; }
  void checkenergy(const int i = 1) { chkenergyconservation = i; }
  void LoadProfile(const std::string& fname);

 private:
  void CreateNodes(PHCompositeNode* topNode);
  bool Cell2Abs(RawTowerGeomContainer* towergeom, float phiC, float etaC, float& phi, float& eta);

  RawClusterContainer* _clusters;
  //  BEmcProfile *_emcprof;

  BEmcRec* bemc;
  float fEnergyNorm;

  float _min_tower_e;
  int chkenergyconservation;

  std::string detector;
  std::string ClusterNodeName;

  int BINX0;
  int NBINX;
  int BINY0;
  int NBINY;

  bool bPrintGeom;
  bool bProfProb;
};

#endif /* RawClusterBuilderTemplate_H__ */
