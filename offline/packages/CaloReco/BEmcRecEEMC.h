#ifndef CALORECO_BEMCRECEEMC_H
#define CALORECO_BEMCRECEEMC_H

#include "BEmcRec.h"

class BEmcProfile;

class BEmcRecEEMC : public BEmcRec
{
 public:
  BEmcRecEEMC();
  ~BEmcRecEEMC();
  void CorrectEnergy(float energy, float x, float y, float *ecorr);
  void CorrectECore(float ecore, float x, float y, float *ecorecorr);
  void CorrectPosition(float energy, float x, float y, float& xcorr, float& ycorr);
  void CorrectShowerDepth(float energy, float x, float y, float z, float& xc, float& yc, float& zc );
  static float GetImpactAngle(float e, float x, float y);

  void LoadProfile(const char *fname);
  float GetProb(std::vector<EmcModule> HitList, float e, float xg, float yg, float zg, float &chi2, int &ndf);

 private:
  BEmcProfile *_emcprof;
};

#endif
