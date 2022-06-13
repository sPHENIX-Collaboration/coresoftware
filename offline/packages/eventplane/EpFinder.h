// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EPFINDER_H
#define EPFINDER_H

#include <fun4all/SubsysReco.h>

#include <string>
#include "TH2D.h"
#include "TVector3.h"
#include "EpInfo.h"
#include <vector>
#include <utility>

class TVector3;
class TProfile;
class TProfile2D;
class PHCompositeNode;

/*
\author Mike Lisa
\date 23 June 2018
\adjusted for sPHENIX by J. Lajoie
\24 August 2019
\Further adjustments for a number of sPHENIX
\detectors by Ejiro
\04 June 2022
*/

#define _EpTermsMax 6

typedef struct{
    
 float nMip;
 double phi;
 int ix;
 int iy;
 std::vector<std::pair<int,int>> *samePhi;
  
} EpHit;

class EpFinder {
  public:
    
    
  EpFinder(int nEventTypeBins=10);
  ~EpFinder() {};

  void SetnMipThreshold(double thresh){mThresh=thresh;};
  void SetMaxTileWeight(double MAX){mMax=MAX;};
  void Finish();
  EpInfo Results(std::vector<EpHit> *EpdHits, int EventTypeID);
  TString Report();

 private:
    
  bool OrderOutsideRange(int order);
  double GetPsiInRange(double Qx, double Qy, int order);
  double mThresh;
  double mMax;
  int mNumberOfEventTypeBins;  
};

#endif // EPFINDER_H
