#ifndef __PHG4TPCCLUSTERIZER_H__
#define __PHG4TPCCLUSTERIZER_H__

#include <g4detectors/PHG4CellDefs.h>
#include <fun4all/SubsysReco.h>
#include <RVersion.h>
#include <vector>
#include <limits.h>

class PHG4CylinderCellGeom;
class TH1F;
class TProfile2D;
class TStopwatch;

class PHG4TPCClusterizer : public SubsysReco {
 public:
  PHG4TPCClusterizer(const char *name = "PHG4TPCClusterizer");
  ~PHG4TPCClusterizer();

  int Init(PHCompositeNode *topNode) { return 0; }
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode) { return 0; }

  void setEnergyCut(float val) { fEnergyCut = val; }
  void setFitWindowSigmas(float rp, float rz) { fDCT = rp; fDCL = rz; }
  void setFitWindowMax(int rp, int rz) { fFitRangeMP = rp; fFitRangeMZ = rz; }
  void setRangeLayers(unsigned int minLayer, unsigned int maxLayer) {fMinLayer=minLayer; fMaxLayer=maxLayer;}
  void setClusterCut(float val) { fClusterCut = val; }
  void setPedestal(float val) { fPedestal = val; }
  void setShapingRMSLead(float val) {fShapingLead = val;}
  void setShapingRMSTail(float val) {fShapingTail = val;}
  void setClusterWindow(float val)  {fClusterWindow = val;}
  void setClusterZSplit(bool val)   {fClusterZSplit = val;}
  void setDeconvolutionMode(bool val)   {fDeconMode = val;}

 private:
  void reset();
  void prepare_layer(float radius);
  void deconvolution();
  int  wrap_phibin(int bin);
  bool is_local_maximum(int phi, int z);
  void find_z_range(int zbin, int phibin, int zmax, float peak, int& zup, int& zdown);
  void find_phi_range(int zbin, int phibin, int phimax, float peak, int& phiup, int& phidown);

  void fit(int pbin, int zbin, int& nhits_tot);
  // FitSumP = weighted sum of dphi values in cluster (ee*dphi)
  // FitSumZ = weighted sum of dz values in cluster
  // fFitW is the sum of weights in the cluster ( sigma(ee) )
  float fit_p_mean() {return fFitSumP/fFitW+fFitP0;}
  float fit_z_mean() {return fFitSumZ/fFitW+fFitZ0;}

  // FitSumP2 = weighted sum of dphi*dphi values in cluster (ee*dphi^2)
  // FitSumZ2 = weighted sum of dz*dz values in cluster
  // So fit_p_cov = sigma(dphi^2*ee)/sigma(ee) - ( sigma(dphi*ee)^2 / sigma(ee)^2 ) = weighted mean of dphi^2 - (weighted mean of dphi)^2 
  float fit_p_cov() {return fFitSumP2/fFitW-fFitSumP/fFitW*fFitSumP/fFitW;}
  float fit_z_cov() {return fFitSumZ2/fFitW-fFitSumZ/fFitW*fFitSumZ/fFitW;}
  float fit_pz_cov() {return fFitSumPZ/fFitW-fFitSumP/fFitW*fFitSumZ/fFitW;}

  std::vector<int> fNHitsPerZ;
  std::vector<float> fAmps;
  std::vector<unsigned int> fCellIDs;
  std::vector<int> fCellz;
  std::vector<int> fCellphi;
  
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
  float fClusterCut;
  float fPedestal;;
  float fFitSizeP;
  int fFitSizeZ;
  float fShapingLead;
  float fShapingTail;
  unsigned int fMinLayer;
  unsigned int fMaxLayer;
  float fEnergyCut;
  float fClusterWindow;
  bool  fClusterZSplit;
  bool  fDeconMode;
  float fDCT;
  float fDCL;
  float _inv_sqrt12;
  float _twopi;
  float zz_shaping_correction;

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
#ifndef __CINT__
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 10, 4)
  double **fSource;
  double **fResponse;
#else
  float **fSource;
  float **fResponse;
#endif
#endif
};

#endif
