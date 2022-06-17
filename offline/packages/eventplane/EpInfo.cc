#include "EpInfo.h"

#include <TVector2.h>  // for TVector2

#include <cmath>  // for fabs, M_PI
#include <iostream>

EpInfo::EpInfo()
{
  QrawOneSide.resize(_EpOrderMax);
  for (auto &vec: QrawOneSide)
  {
    vec.resize(2);
  }
  WheelSumWeightsRaw.resize(_EpOrderMax);
  Reset();
}

void EpInfo::Reset()
{
  for (auto &vec: QrawOneSide)
  {
    std::fill(vec.begin(),vec.end(),NAN);
  }
/*
  for (unsigned int iorder = 0; iorder < MaxOrder(); iorder++)
  {
    for (int xy = 0; xy < 2; xy++)
    {
      QrawOneSide[iorder][xy] = 0.0;
    }
  }
*/
  for (unsigned int iorder = 0; iorder < MaxOrder(); iorder++)
  {
    PsiRaw[iorder] = -999.0;
    WheelSumWeightsRaw[iorder] = -999.0;
  }
}

void EpInfo::InitializeToZero()
{
  std::fill(std::begin(WheelSumWeightsRaw), std::end(WheelSumWeightsRaw), 0.);
}

// ===================== Access to Q-vectors ==========================

//------------------------------ Raw Q --------------------------------
std::pair<double, double> EpInfo::RawQ(unsigned int order)
{
  if (ArgumentOutOfBounds(order))
  {
    return std::make_pair(NAN,NAN);
  }
  return std::make_pair(QrawOneSide[order - 1][0], QrawOneSide[order - 1][1]);
}

// --------------------- Wheel sum-of-weights, raw ----------------------
double EpInfo::SWRaw(unsigned int order)
{
  if (ArgumentOutOfBounds(order)) return -999;
  return WheelSumWeightsRaw[order - 1];
}

// ===================== Access to Event-plane angles ====================

//------------------------- raw EP angles --------------------------------
double EpInfo::RawPsi(unsigned int order)
{
  if (ArgumentOutOfBounds(order)) return -999;
  return Range(PsiRaw[order - 1], order);
}
//-----------------------------------------------------------------------

//----- Simple method to put angles in a convenient range: (0,2pi/n) ----
double EpInfo::Range(double psi, unsigned int order)
{
  if (ArgumentOutOfBounds(order)) return -999;
  double wrap = (2.0 * M_PI) / (double) order;
  if (psi < 0.0)
    psi += (1.0 + (int) (fabs(psi) / wrap)) * wrap;
  else
  {
    if (psi > wrap) psi -= ((int) (psi / wrap)) * wrap;
  }
  return psi;
}

//--------- protection against argument out-of-bounds -------
bool EpInfo::ArgumentOutOfBounds(unsigned int order)
{
  if ((order < 1) || ((unsigned int)order > MaxOrder()))
  {
    std::cout << "\n *** Invalid order requested " << order << "***\n";
    std::cout << "  order must be between 1 (for first-order EP) and " << MaxOrder()
              << ".    To change the upuper limit, edit StEpdUtil/EpInfo.h\n";
    std::cout << "  I will now return you an invalid result.  Have a nice day\n";
    return true;
  }
  return false;
}
