#ifndef __QPILEUP_H__
#define __QPILEUP_H__

//===========================================================
/// \file QPileUp.h
/// \brief Base class for initial charge density
/// \author Carlos Perez Lara
//===========================================================

#include <string>
#include "TH3F.h"

class QPileUp {
 public:
  QPileUp();
  ~QPileUp();
  void SetDebugLevel(int n) {fDebug=n;}
  virtual void Make();
  void TPCDimensions(float irad, float orad, float hzet) {fInnerRadius=irad; fOutterRadius=orad; fHalfLength=hzet;}
  void TPCGridSize(int nr, int np, int nz) {fNRadialSteps=nr; fNAzimuthalSteps=np; fNLongitudinalSteps=nz;}
  void OutputFileName(std::string a) {fFileNameRoot=a;}

 protected:
  void InitMaps();
  void SaveMaps();
  int fDebug;

  TH3F *fRho;

  float fInnerRadius;
  float fOutterRadius;
  float fHalfLength;
  int fNRadialSteps;
  int fNAzimuthalSteps;
  int fNLongitudinalSteps;
  std::string fFileNameRoot;
};

#endif /* __QPILEUP_H__ */
