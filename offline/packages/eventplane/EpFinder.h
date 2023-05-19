// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EVENTPLANE_EPFINDER_H
#define EVENTPLANE_EPFINDER_H

#include <TString.h>
#include <string>

#include <array>  // for array
#include <vector>
#include "TH2D.h"

class EpInfo;


typedef struct
{
  float nMip;
  double phi;
  int ix;
  int iy;
  int wheel;
  std::vector<std::pair<int,int>> *sameRing;

} EpHit;

class EpFinder
{
 public:
  EpFinder(int nEventTypeBins = 10, unsigned int order = 3, int ns = 0);
  ~EpFinder(){};

  void Finish(int thishist = 0);
  void Results(const std::vector<EpHit> &Hits, int EventTypeID, EpInfo *epinfo);
  void ResultsEPD(const std::vector<EpHit> &EpdHits, int EventTypeID, EpInfo *epinfo);
  TString Report();
  void ResetEvent();  // clear current event for safety

 private:
  double GetPsiInRange(const double Qx, const double Qy, const unsigned int order) const;
  int mNumberOfEventTypeBins = 0;
  unsigned int m_MaxOrder = 0;

  std::vector<std::array<double, 2>> TotalWeight4Side;
  std::vector<std::vector<double>> QrawOneSide;
  std::vector<std::vector<double>> QphiWeightedOneSide;
  std::vector<double> WheelSumWeightsRaw;
  std::vector<double> WheelSumWeightsPhiWeighted;
  std::vector<double> PsiRaw;
  std::vector<double> PsiPhiWeighted;

  TH2D* mPhiWeightInput[2];
  TH2D* mPhiWeightOutput[2];
  TH2D* mPhiAveraged[2];
};

#endif  // EVENTPLANE_EPFINDER_H
