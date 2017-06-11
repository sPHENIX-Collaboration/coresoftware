#ifndef __QPILEUPTOY_H__
#define __QPILEUPTOY_H__

//===========================================================
/// \file QpileUpToy.h
/// \brief Initial charge density from phenomenological model
/// \author Carlos Perez Lara
//===========================================================

#include "QPileUp.h"

class QPileUpToy:public QPileUp {
 public:
  QPileUpToy(float gasf/*[Vs]*/, float mult, float rate/*[Hz]*/, float eps, float rad=2.0);
  virtual void Make();

 protected:
  float fGasFactor;
  float fMultiplicity;
  float fDAQRate;
  float fEPS;
  float fRad;
};

#endif /* __QPILEUPTOY_H__ */
