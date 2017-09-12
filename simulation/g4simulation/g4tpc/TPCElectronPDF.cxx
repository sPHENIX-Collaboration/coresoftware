// Stores the electron probability distribution function
// Author: Carlos Perez
#include <TPCConstants.h>
#include <TPCHit.h>
#include "TPCElectronPDF.h"
#include <cmath>

using namespace TPCConstants;

//=====
TPCElectronPDF::TPCElectronPDF()
{
}
//=====
TPCElectronPDF::~TPCElectronPDF()
{
}
//=====
void TPCElectronPDF::AddElectron(float w)
{
  fW.push_back(w);
  fDX.push_back(0.0);
  fDY.push_back(0.0);
  fDZ.push_back(0.0);
  fT.push_back(0.0);
  fMST.push_back(0.0);
  fMSL0.push_back(0.0);
  fMSL1.push_back(0.0);
}
//=====
void TPCElectronPDF::Reset()
{
  fHit = NULL;
  fW.clear();
  fDX.clear();
  fDY.clear();
  fDZ.clear();
  fT.clear();
  fMST.clear();
  fMSL0.clear();
  fMSL1.clear();
}
//=====
void TPCElectronPDF::CopyFrom(TPCHit *th)
{
  Reset();
  TPCHit::CopyFrom(th);
}
//=====
float TPCElectronPDF::GetR(unsigned int i) {
  if(i>=GetN()) return 0;
  float x = GetX() + GetDx(i);
  float y = GetY() + GetDy(i);
  return std::sqrt( x*x + y*y );
}
//=====
float TPCElectronPDF::GetPhi(unsigned int i) {
  if(i>=GetN()) return 0;
  float x = GetX() + GetDx(i);
  float y = GetY() + GetDy(i);
  return kPi+std::atan2(-y,-x);
}
