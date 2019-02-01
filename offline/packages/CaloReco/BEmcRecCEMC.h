#ifndef CALORECO_BEMCRECCEMC_H
#define CALORECO_BEMCRECCEMC_H

#include "BEmcRec.h"

class BEmcRecCEMC : public BEmcRec
{
 public:
  BEmcRecCEMC(){};
  ~BEmcRecCEMC() {}
  void CorrectEnergy(float energy, float x, float y, float *ecorr);
  void CorrectPosition(float energy, float x, float y, float *xcorr, float *ycorr);
  void CorrectECore(float ecore, float x, float y, float *ecorecorr);
  void Tower2Global(float E, float xC, float yC, float &xA, float &yA, float &zA);
};

#endif
