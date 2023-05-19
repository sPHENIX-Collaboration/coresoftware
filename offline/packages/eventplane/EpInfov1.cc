#include "EpInfov1.h"

#include <algorithm>  // for fill
#include <cmath>      // for fabs, M_PI
#include <iostream>
#include <memory>       // for allocator_traits<>::value_type
#include <type_traits>  // for __decay_and_strip<>::__type

void EpInfov1::Reset()
{
  for (auto &vec : QrawOneSide)
  {
    std::fill(vec.begin(), vec.end(), NAN);
  }
    
 for (auto &vec : QphiWeightedOneSide)
 {
    std::fill(vec.begin(), vec.end(), NAN);
 }
  std::fill(PsiRaw.begin(), PsiRaw.end(), NAN);
  std::fill(PsiPhiWeighted.begin(), PsiPhiWeighted.end(), NAN);
  std::fill(WheelSumWeightsRaw.begin(), WheelSumWeightsRaw.end(), NAN);
  std::fill(WheelSumWeightsPhiWeighted.begin(), WheelSumWeightsPhiWeighted.end(), NAN);

}

// ===================== Access to Q-vectors ==========================

//------------------------------ Raw Q --------------------------------
std::pair<double, double> EpInfov1::RawQ(unsigned int order) const
{
  if (ArgumentOutOfBounds(order))
  {
    return std::make_pair(NAN, NAN);
  }
  return std::make_pair(QrawOneSide.at(order - 1).at(0), QrawOneSide.at(order - 1).at(1));
}

//------------------------ phi-weighted Q -------------------------------

std::pair<double, double> EpInfov1::PhiWeightedQ(unsigned int order) const
{
  if (ArgumentOutOfBounds(order))
  {
    return std::make_pair(NAN, NAN);
  }
  return std::make_pair(QphiWeightedOneSide.at(order - 1).at(0), QphiWeightedOneSide.at(order - 1).at(1));
}

// --------------------- Wheel sum-of-weights, raw ----------------------
double EpInfov1::SWRaw(unsigned int order) const
{
  if (ArgumentOutOfBounds(order)) return NAN;
  return WheelSumWeightsRaw.at(order - 1);
}

// --------------------- Wheel sum-of-weights, phi-weighted ---------------
double EpInfov1::SWPhiWeighted(unsigned int order) const
{
  if (ArgumentOutOfBounds(order)) return NAN;
  return WheelSumWeightsPhiWeighted.at(order - 1);
}
// ===================== Access to Event-plane angles ====================

//------------------------- raw EP angles --------------------------------
double EpInfov1::RawPsi(unsigned int order) const
{
  if (ArgumentOutOfBounds(order)) return NAN;
  return Range(PsiRaw.at(order - 1), order);
}

//-------------------- phi-weighted EP angles ---------------------------
double EpInfov1::PhiWeightedPsi(unsigned int order) const
{
  if (ArgumentOutOfBounds(order)) return NAN;
  return Range(PsiPhiWeighted.at(order - 1), order);
}

//-----------------------------------------------------------------------

//----- Simple method to put angles in a convenient range: (0,2pi/n) ----
double EpInfov1::Range(double psi, unsigned int order) const
{
  if (ArgumentOutOfBounds(order)) return NAN;
  double wrap = (2.0 * M_PI) / (double) order;
  if (psi < 0.0)
  {
    psi += (1.0 + (int) (fabs(psi) / wrap)) * wrap;
  }
  else
  {
    if (psi > wrap) psi -= ((int) (psi / wrap)) * wrap;
  }
  return psi;
}

//--------- protection against argument out-of-bounds -------
bool EpInfov1::ArgumentOutOfBounds(unsigned int order) const
{
  if ((order < 1) || ((unsigned int) order > MaxOrder()))
  {
    std::cout << "\n *** Invalid order requested " << order << "***\n";
    std::cout << "  order must be between 1 (for first-order EP) and " << MaxOrder()
              << ".    To change the upuper limit, edit StEpdUtil/EpInfo.h\n";
    std::cout << "  I will now return you an invalid result.  Have a nice day\n";
    return true;
  }
  return false;
}

void EpInfov1::CopyPsiRaw(const std::vector<double> &vec)
{
  if (PsiRaw.size() != vec.size())
  {
    PsiRaw.resize(vec.size());
  }
  PsiRaw = vec;
}

void EpInfov1::CopyPsiPhiWeighted(const std::vector<double> &vec)
{
  if (PsiPhiWeighted.size() != vec.size())
  {
      PsiPhiWeighted.resize(vec.size());
  }
    PsiPhiWeighted = vec;
}

void EpInfov1::CopyWheelSumWeightsRaw(const std::vector<double> &vec)
{
  if (WheelSumWeightsRaw.size() != vec.size())
  {
    WheelSumWeightsRaw.resize(vec.size());
  }
  WheelSumWeightsRaw = vec;
}

void EpInfov1::CopyWheelSumWeightsPhiWeighted(const std::vector<double> &vec)
{
  if (WheelSumWeightsPhiWeighted.size() != vec.size())
  {
      WheelSumWeightsPhiWeighted.resize(vec.size());
  }
    WheelSumWeightsPhiWeighted = vec;
}

void EpInfov1::CopyQrawOneSide(const std::vector<std::vector<double>> &vecvec)
{
  if (QrawOneSide.size() != vecvec.size())
  {
    QrawOneSide.resize(vecvec.size());
    int i = 0;
    for (auto &vec : vecvec)
    {
      QrawOneSide.at(i).resize(vec.size());
      i++;
    }
  }
  QrawOneSide = vecvec;
}

void EpInfov1::CopyQphiWeightedOneSide(const std::vector<std::vector<double>> &vecvec)
{
  if (QphiWeightedOneSide.size() != vecvec.size())
  {
      QphiWeightedOneSide.resize(vecvec.size());
    int i = 0;
    for (auto &vec : vecvec)
    {
        QphiWeightedOneSide.at(i).resize(vec.size());
      i++;
    }
  }
    QphiWeightedOneSide = vecvec;
}
