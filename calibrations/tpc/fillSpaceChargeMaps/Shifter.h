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
  
  TVector3 Shift(TVector3 position);
  TVector3 ShiftForward(TVector3 position); //only shift with forward histogram
  TVector3 ShiftBack(TVector3 position); //
  TFile *forward, *back, *average;
  bool hasTruth, hasCorrection;
  TH3F *hX, *hY, *hZ, *hR, *hPhi, *hXave, *hYave, *hZave, *hRave, *hPhiave, *hXBack, *hYBack, *hZBack;  
};

#endif
