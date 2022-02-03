//Shifter code
#ifndef SHIFTER_H
#define SHIFTER_H

#include <TString.h>

class TFile;
class TVector3;
class TH3F;

class Shifter {
public:
  Shifter(TString truthfilename) :Shifter(truthfilename, "" ){return;} ;
  Shifter(TString truthfilename, TString correctionfilename);
  
  TVector3 Shift(TVector3& position);
  TVector3 ShiftForward(TVector3 position); //only shift with forward histogram
  TVector3 ShiftBack(TVector3 position); //
  TFile *forward=0, *back=0, *average=0;
  bool hasTruth, hasCorrection;
  TH3F *hX = 0;
  TH3F *hY = 0;
  TH3F *hZ = 0;
  TH3F *hR = 0;
  TH3F *hPhi = 0;
  TH3F *hXave = 0;
  TH3F *hYave = 0;
  TH3F *hZave = 0;
  TH3F *hRave = 0;
  TH3F *hPhiave = 0;
  TH3F *hXBack = 0;
  TH3F *hYBack = 0;
  TH3F *hZBack = 0;
};

#endif
