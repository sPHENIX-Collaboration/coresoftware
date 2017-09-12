// Stores  one hit in TPC fiducial volume
// Author: Carlos Perez
#include <cmath>
#include "TPCConstants.h"
#include "TPCHit.h"

using namespace TPCConstants;

//=====
TPCHit::TPCHit() : fHit(NULL) {
}
//=====
TPCHit::~TPCHit() {}
//=====
float TPCHit::GetR() {
  float x = GetX();
  float y = GetY();
  return std::sqrt( x*x + y*y );
}
//=====
float TPCHit::GetPhi() {
  float x = GetX();
  float y = GetY();
  return kPi + std::atan2( -y, -x );
}
