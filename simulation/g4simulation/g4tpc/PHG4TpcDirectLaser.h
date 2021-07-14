#ifndef PHG4TPCDIRECTLASER_H
#define PHG4TPCDIRECTLASER_H

#include <iostream>
#include <cmath>
#include <vector>
#include "TMath.h"
#include "TVector3.h"

//from phg4tpcsteppingaction.cc
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Hitv1.h>
//R__LOAD_LIBRARY(libphg4hit.so)


// all distances in mm, all angles in rad
// class that generates stripes and dummy hit coordinates
// stripes have width of one mm, length of one pad width, and are centered in middle of sector gaps

using namespace std;

class PHG4TpcDirectLaser {
public:
  PHG4TpcDirectLaser(); //default constructor

  double begin_CM, end_CM; // inner and outer radii of field cages/TPC
  double halfwidth_CM; //half the thickness of the CM;
  double ifc,ofc;
  
  vector<PHG4Hitv1*> PHG4Hits;

  void SetPhiStepping(int n, float min,float max);
  void SetThetaStepping(int n, float min,float max);
  int GetNpatternSteps(){return nPhiSteps*nThetaSteps;};
  void AimToThetaPhi(float theta, float phi);
  void AimToPatternStep(int n);
  void AimToNextPatternStep(){if (nTotalSteps>1)AimToPatternStep(currentPatternStep+1);};
  
private:
  static const int nLasers = 4; //per side
  const double mm = 0.10;
  const double cm = 1.0;
  const float maxHitLength=1.0;//1cm.

  int nPhiSteps=1;
  int nThetaSteps=1;
  int nTotalSteps=1;
  int currentPatternStep=0;
  float minPhi=0;
  float maxPhi=0;
  float minTheta=0;
  float maxTheta=0;

  TVector3 GetCmStrike(TVector3 start, TVector3 direction);
  TVector3 GetFieldcageStrike(TVector3 start, TVector3 direction);
  TVector3 GetCylinderStrike(TVector3 s, TVector3 v, float radius);
  
  int nElectrons;
 
  void AppendLaserTrack(float theta, float phi, int laser);
  void ClearHits();
};


#endif
