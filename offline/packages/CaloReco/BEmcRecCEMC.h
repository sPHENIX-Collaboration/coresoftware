#ifndef CALORECO_BEMCRECCEMC_H
#define CALORECO_BEMCRECCEMC_H

#include "BEmcRec.h"

class BEmcProfile;

class BEmcRecCEMC : public BEmcRec
{
 public:
  BEmcRecCEMC();
  virtual ~BEmcRecCEMC();
  void CorrectEnergy(float energy, float x, float y, float *ecorr);
  void CorrectECore(float ecore, float x, float y, float *ecorecorr);
  void CorrectPosition(float energy, float x, float y, float& xcorr, float& ycorr);
  void CorrectShowerDepth(float energy, float x, float y, float z, float& xc, float& yc, float& zc );

  void LoadProfile(const std::string &fname) override;
  float GetProb(std::vector<EmcModule> HitList, float e, float xg, float yg, float zg, float &chi2, int &ndf);

 private:
  BEmcProfile *_emcprof;
};

#endif
