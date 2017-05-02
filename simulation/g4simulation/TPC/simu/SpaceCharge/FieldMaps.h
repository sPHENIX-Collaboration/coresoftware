#ifndef __FIELDMAPS_H__
#define __FIELDMAPS_H__

//===========================================================
/// \file FieldMaps.h
/// \brief Base class for E field distortions
/// \author Carlos Perez Lara
//===========================================================

#include "TH3F.h"

class FieldMaps {
 public:
  FieldMaps();
  virtual ~FieldMaps();
  void ReadFile();
  void SetDebugLevel(int n) {fDebug=n;}
  void TPCDimensions(float irad, float orad, float hzet) {fInnerRadius=irad; fOutterRadius=orad; fHalfLength=hzet;}
  void TPCGridSize(int nr, int np, int nz) {fNRadialSteps=nr; fNAzimuthalSteps=np; fNLongitudinalSteps=nz;}
  void OutputFileName(std::string a) {fFileNameRoot=a;}
  void LSFileName(std::string a) {fLSNameRoot=a;}
  void MirrorZ() {fMirrorZ=true;}
  void Make(int n=-1);
  virtual void ComputeE()=0;

 protected:
  void InitMaps();
  void SaveMaps();
  float ReadCharge(float rprime, float phiprime, float zprime, float dr, float dphi,float dz);

  int fDebug;
  int fRadialBin;

  TH3F *fEr;
  TH3F *fEp;
  TH3F *fEz;
  TH3F *fRho;

  float fInnerRadius;
  float fOutterRadius;
  float fHalfLength;
  int fNRadialSteps;
  int fNAzimuthalSteps;
  int fNLongitudinalSteps;
  std::string fFileNameRoot;
  std::string fLSNameRoot;
  bool fMirrorZ;
};

#endif /* __FIELDMAPS_H__ */
