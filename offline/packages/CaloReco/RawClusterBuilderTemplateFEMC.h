#ifndef CALORECO_RAWCLUSTERBUILDERTEMPLATEFEMC_H
#define CALORECO_RAWCLUSTERBUILDERTEMPLATEFEMC_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class RawClusterContainer;
class RawTowerGeomContainer;
class BEmcRecFEMC;
class BEmcProfile;

class RawClusterBuilderTemplateFEMC : public SubsysReco
{
 public:
  RawClusterBuilderTemplateFEMC(const std::string &name = "RawClusterBuilderGraph");
  virtual ~RawClusterBuilderTemplateFEMC();

  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
  void Detector(const std::string &d) { detector = d; }

  void set_threshold_energy(const float e) { _min_tower_e = e; }
  void setEnergyNorm(float norm) { fEnergyNorm = norm; }
  void checkenergy(const int i = 1) { chkenergyconservation = i; }
  void LoadProfile(const char *fname);

 private:
  void CreateNodes(PHCompositeNode *topNode);
  bool Cell2Abs(RawTowerGeomContainer *towergeom, float phiC, float etaC, float &phi, float &eta);

  RawClusterContainer *_clusters;
  BEmcProfile *_emcprof;

  BEmcRecFEMC *bemc;
  float fEnergyNorm;

  float _min_tower_e;
  int chkenergyconservation;

  std::string detector;
  std::string ClusterNodeName;

  int BINX0;
  int NBINX;
  int BINY0;
  int NBINY;
  float Zcenter;
};

#endif
