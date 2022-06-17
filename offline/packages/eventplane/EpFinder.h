// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EVENTPLANE_EPFINDER_H
#define EVENTPLANE_EPFINDER_H

#include <TString.h>

#include <array>  // for array
#include <vector>

class EpInfo;

/*
\author Mike Lisa
\date 23 June 2018
\adjusted for sPHENIX by J. Lajoie
\24 August 2019
\Further adjustments for a number of sPHENIX
\detectors by Ejiro
\04 June 2022
*/

typedef struct
{
  float nMip;
  double phi;

} EpHit;

class EpFinder
{
 public:
  EpFinder(int nEventTypeBins = 10, unsigned int order = 3);
  ~EpFinder(){};

  void SetnMipThreshold(const double thresh) { mThresh = thresh; };
  void SetMaxTileWeight(const double MAX) { mMax = MAX; };
  void Results(const std::vector<EpHit> &EpdHits, int EventTypeID, EpInfo *epinfo);
  TString Report();
  void ResetEvent();  // clear current event for safety

 private:
  double GetPsiInRange(const double Qx, const double Qy, const unsigned int order) const;
  double mThresh = 0.;
  double mMax = 100.;
  int mNumberOfEventTypeBins = 0;
  unsigned int m_MaxOrder = 0;

  std::vector<std::array<double, 1>> TotalWeight4Side;
  std::vector<std::vector<double>> QrawOneSide;
  std::vector<double> WheelSumWeightsRaw;
  std::vector<double> PsiRaw;
};

#endif  // EVENTPLANE_EPFINDER_H
