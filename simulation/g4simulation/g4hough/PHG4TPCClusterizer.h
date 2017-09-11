#ifndef __PHG4TPCCLUSTERIZER_H__
#define __PHG4TPCCLUSTERIZER_H__

#include <fun4all/SubsysReco.h>
#include <vector>
#include <limits.h>

class PHG4CylinderCellGeom;
class TH1F;
class TProfile2D;
class TStopwatch;

class PHG4TPCClusterizer : public SubsysReco {
 public:
  PHG4TPCClusterizer(const char *name = "PHG4SvtxClusterizer");
  ~PHG4TPCClusterizer();

  int Init(PHCompositeNode *topNode) { return 0; }
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode) { return 0; }

  void setEnergyCut(float val) { fEnergyCut = val; }
  void setFitWindowSigmas(float rp, float rz) { fDCT = rp; fDCL = rz; }
  void setFitWindowMax(int rp, int rz) { fFitRangeMP = rp; fFitRangeMZ = rz; }
  void setRangeLayers(unsigned int minLayer, unsigned int maxLayer) {fMinLayer=minLayer; fMaxLayer=maxLayer;}
  void setFitEnergyThreshold(float val) { fFitEnergyThreshold = val; }
  void setShapingRMSLead(float val) {fShapingLead = val;}
  void setShapingRMSTail(float val) {fShapingTail = val;}

 private:
  void reset();
  int wrap_phibin(int bin);
  bool is_local_maximum(int phi, int z);
  void fit(int pbin, int zbin, int& nhits_tot);
  float fit_p_mean() {return fFitSumP/fFitW+fFitP0;}
  float fit_z_mean() {return fFitSumZ/fFitW+fFitZ0;}

  float fit_p_cov() {return fFitSumP2/fFitW-fFitSumP/fFitW*fFitSumP/fFitW;}
  float fit_z_cov() {return fFitSumZ2/fFitW-fFitSumZ/fFitW*fFitSumZ/fFitW;}
  float fit_pz_cov() {return fFitSumPZ/fFitW-fFitSumP/fFitW*fFitSumZ/fFitW;}

  std::vector<int> fNHitsPerZ;
  std::vector<float> fAmps;
  std::vector<int> fCellIDs;
  int fNPhiBins;
  int fNZBins;
  PHG4CylinderCellGeom *fGeoLayer;

  float fFitW;
  float fFitSumP;
  float fFitSumZ;
  float fFitSumP2;
  float fFitSumZ2;
  float fFitSumPZ;
  float fFitP0;
  float fFitZ0;
  int fFitRangeP;
  int fFitRangeZ;
  int fFitRangeMP;
  int fFitRangeMZ;
  float fFitEnergyThreshold;
  float fFitSizeP;
  int fFitSizeZ;
  float fShapingLead;
  float fShapingTail;
  unsigned int fMinLayer;
  unsigned int fMaxLayer;
  float fEnergyCut;

  float fDCT;
  float fDCL;

  float _inv_sqrt12;
  float _twopi;

  TH1F *fHClusterEnergy;
  TProfile2D *fHClusterSizePP;
  TProfile2D *fHClusterSizeZZ;
  TProfile2D *fHClusterErrorPP;
  TProfile2D *fHClusterErrorZZ;
  TProfile2D *fHClusterDensity;
  TProfile2D *fHClusterSizePP2;
  TProfile2D *fHClusterSizeZZ2;
  TProfile2D *fHClusterErrorPP2;
  TProfile2D *fHClusterErrorZZ2;
  TProfile2D *fHClusterDensity2;
  TProfile2D *fHClusterWindowP;
  TProfile2D *fHClusterWindowZ;
  TStopwatch *fSW;
  TH1F *fHTime;

};

#endif
