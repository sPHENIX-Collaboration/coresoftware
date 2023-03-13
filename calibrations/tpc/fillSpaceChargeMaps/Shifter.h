//Shifter code
#ifndef FILLSPACECHARGEMAPS_SHIFTER_H
#define FILLSPACECHARGEMAPS_SHIFTER_H

#include <string>

class TFile;
class TVector3;
class TH3;

class Shifter
{
 public:
  explicit Shifter(const std::string &truthfilename, const std::string &correctionfilename = "");

  TVector3 Shift(const TVector3 &position);
  TVector3 ShiftForward(const TVector3 &position);  //only shift with forward histogram
  TVector3 ShiftBack(const TVector3 &position);     //
  TFile *forward = nullptr;
  TFile *back = nullptr;
  TFile *average = nullptr;
  bool hasTruth = false;
  bool hasCorrection = false;
  TH3 *hX = nullptr;
  TH3 *hY = nullptr;
  TH3 *hZ = nullptr;
  TH3 *hR = nullptr;
  TH3 *hPhi = nullptr;
  TH3 *hXave = nullptr;
  TH3 *hYave = nullptr;
  TH3 *hZave = nullptr;
  TH3 *hRave = nullptr;
  TH3 *hPhiave = nullptr;
  TH3 *hXBack = nullptr;
  TH3 *hYBack = nullptr;
  TH3 *hZBack = nullptr;
};

#endif
