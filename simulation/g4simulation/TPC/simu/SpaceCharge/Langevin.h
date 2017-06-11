#ifndef __LANGEVIN_H__
#define __LANGEVIN_H__

//===========================================================
/// \file Langevin.h
/// \brief Langevin solution for charge transport
/// \author Carlos Perez Lara
//===========================================================

#include "TH3F.h"

class Langevin {
 public:
  Langevin();
  virtual ~Langevin();
  void ReadFile();
  void SetDebugLevel(int n) {fDebug=n;}
  void Make();
  void SetMirrorZ() {fMirrorZ=true;}
  void TPCDimensions(float irad, float orad, float hzet) {fInnerRadius=irad; fOutterRadius=orad; fHalfLength=hzet;}
  void TPCGridSize(int nr, int np, int nz) {fNRadialSteps=nr; fNAzimuthalSteps=np; fNLongitudinalSteps=nz;}
  void OutputFileName(std::string a) {fFileNameRoot=a;}

 protected:
  void InitMaps();
  void SaveMaps();
  int fDebug;

  TH3F *fEr;
  TH3F *fEp;
  TH3F *fEz;
  TH3F *fDeltaR;
  TH3F *fRDeltaPHI;

  bool fMirrorZ;
  float fInnerRadius;
  float fOutterRadius;
  float fHalfLength;
  int fNRadialSteps;
  int fNAzimuthalSteps;
  int fNLongitudinalSteps;
  std::string fFileNameRoot;
};

#endif /* __Langevin_H__ */
