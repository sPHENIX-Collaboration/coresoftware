#ifndef EpInfo_H
#define EpInfo_H

#define _EpOrderMax 3

#include "TVector2.h"

#include <phool/PHObject.h>

class EpInfo : public PHObject {

 friend class EpFinder;

 public:
  EpInfo();
  ~EpInfo(){/* no op */};

  void Reset(); 


  TVector2 RawQ(int order);
 
  TVector2 PhiWeightedQ(int order);

  double RawPsi(int order);

  double PhiWeightedPsi(int order);

  double PhiWeightedAndShiftedPsi(int order);

  double SWRaw(int order);


  double PsiRaw[_EpOrderMax];
 
 private:

  bool ArgumentOutOfBounds(int order);

  double Range(double psi, int order);           /// puts angle psi into range (0,2pi/n)

  double QrawOneSide[_EpOrderMax][2];           /// indices: [order][x,y]
  double WheelSumWeightsRaw[_EpOrderMax];       /// indices: [order]

};

#endif
