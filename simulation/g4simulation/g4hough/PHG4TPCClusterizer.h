#ifndef __PHG4TPCCLUSTERIZER__
#define __PHG4TPCCLUSTERIZER__

#include <fun4all/SubsysReco.h>
#include <vector>
#include <limits.h>

class PHG4CylinderCellGeom;

class PHG4TPCClusterizer : public SubsysReco {
 public:
  PHG4TPCClusterizer(const char *name = "PHG4SvtxClusterizer");
  ~PHG4TPCClusterizer();

  int Init(PHCompositeNode *topNode) { return 0; }
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode) { return 0; }

  void setEnergyCut(float val) { fEnergyCut = val; }
  void setFitWindow(int rp, int rz) { fFitRangeP = rp; fFitRangeZ = rz; }
  void setRangeLayers(unsigned int minLayer, unsigned int maxLayer) {fMinLayer=minLayer; fMaxLayer=maxLayer;}
  void setFitEnergyThreshold(float val) { fFitEnergyThreshold = val; }

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
  float fFitEnergyThreshold;
  float fFitSizeP;
  int fFitSizeZ;

  unsigned int fMinLayer;
  unsigned int fMaxLayer;
  float fEnergyCut;

  float _inv_sqrt12;
};

#endif
