#ifndef CALORECO_BEMCRECFEMC_H
#define CALORECO_BEMCRECFEMC_H

#include "BEmcRec.h"

class BEmcRecFEMC : public BEmcRec
{
 public:
  BEmcRecFEMC(){};
  ~BEmcRecFEMC() {}
  void CorrectEnergy(float energy, float x, float y, float *ecorr);
  static float GetImpactAngle(float e, float x, float y);
  void CorrectPosition(float energy, float x, float y, float *xcorr, float *ycorr);
  void CorrectECore(float ecore, float x, float y, float *ecorecorr);
  void Tower2Global(float E, float xC, float yC, float &xA, float &yA, float &zA);
};

#endif
