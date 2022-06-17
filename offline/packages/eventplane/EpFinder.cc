#include "EpFinder.h"

#include "EpInfo.h"  // for EpInfo

#include <algorithm>  // for fill
#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

using namespace std;

EpFinder::EpFinder(int nEventTypeBins, unsigned int maxorder)
  : mNumberOfEventTypeBins(nEventTypeBins)
  , m_MaxOrder(maxorder)
{
  TotalWeight4Side.resize(m_MaxOrder);
  WheelSumWeightsRaw.resize(m_MaxOrder);
  PsiRaw.resize(m_MaxOrder);
  QrawOneSide.resize(m_MaxOrder);
  for (auto &vec : QrawOneSide)
  {
    vec.resize(2);
  }
}

//==================================================================================================================
void EpFinder::Results(const std::vector<EpHit> &EpdHits, int EventTypeId, EpInfo *epinfo)
{
  if ((EventTypeId < 0) || (EventTypeId >= mNumberOfEventTypeBins))
  {
    cout << "You are asking for an undefined EventType - fail!\n";
    assert(0);
  }

  //--------------------------------
  // begin loop over hits
  //--------------------------------

  for (unsigned int hit = 0; hit < EpdHits.size(); hit++)
  {
    float nMip = EpdHits.at(hit).nMip;
    if (nMip < mThresh) continue;
    double TileWeight = (nMip < mMax) ? nMip : mMax;
    double phi = EpdHits.at(hit).phi;

    //--------------------------------
    // now calculate Q-vectors
    //--------------------------------

    for (unsigned int order = 1; order < m_MaxOrder + 1; order++)
    {
      double etaWeight = 1.0;  // not implemented - JGL 8/27/2019

      //yes the fabs() makes sense.  The sign in the eta weight is equivalent to a trigonometric phase. --JGL
      TotalWeight4Side[order - 1][0] += fabs(etaWeight) * TileWeight;

      double Cosine = cos(phi * (double) order);
      double Sine = sin(phi * (double) order);
      QrawOneSide[order - 1][0] += etaWeight * TileWeight * Cosine;
      QrawOneSide[order - 1][1] += etaWeight * TileWeight * Sine;
    }
  }  // loop over hits

  // Weights used, so you can "un-normalize" the ring-by-ring Q-vectors.
  for (unsigned int order = 1; order < m_MaxOrder + 1; order++)
  {
    WheelSumWeightsRaw[order - 1] = TotalWeight4Side[order - 1][0];
  }

  // at this point, we are finished with the detector hits, and deal only with the Q-vectors,

  // first, normalize the Q's...
  for (unsigned int order = 1; order < m_MaxOrder + 1; order++)
  {
    if (TotalWeight4Side[order - 1][0] > 0.0001)
    {
      QrawOneSide[order - 1][0] /= TotalWeight4Side[order - 1][0];
      QrawOneSide[order - 1][1] /= TotalWeight4Side[order - 1][0];
    }
  }

  // at this point, we are finished with the Q-vectors and just use them to get angles Psi

  //---------------------------------
  // Calculate unshifted EP angles
  //---------------------------------
  for (unsigned int order = 1; order < m_MaxOrder + 1; order++)
  {
    PsiRaw[order - 1] = GetPsiInRange(QrawOneSide[order - 1][0], QrawOneSide[order - 1][1], order);
  }  // loop over order
     // copy results to i/o object
  if (epinfo)
  {
    epinfo->CopyQrawOneSide(QrawOneSide);
    epinfo->CopyWheelSumWeightsRaw(WheelSumWeightsRaw);
    epinfo->CopyPsiRaw(PsiRaw);
  }
  return;
}

double EpFinder::GetPsiInRange(double Qx, double Qy, unsigned int order) const
{
  double temp;
  if ((Qx == 0.0) || (Qy == 0.0))
    temp = -999;
  else
  {
    temp = atan2(Qy, Qx) / ((double) order);
    double AngleWrapAround = (2.0 * M_PI) / (double) order;
    if (temp < 0.0)
      temp += AngleWrapAround;
    else if (temp > AngleWrapAround)
      temp -= AngleWrapAround;
  }
  return temp;
}

TString EpFinder::Report()
{
  TString rep = Form("This is the EpFinder Report:\n");
  rep += Form("Number of EventType bins = %d\n", mNumberOfEventTypeBins);
  rep += Form("Threshold (in MipMPV units) = %f  and MAX weight = %f\n", mThresh, mMax);
  return rep;
}

void EpFinder::ResetEvent()
{
  for (auto &vec : TotalWeight4Side)
  {
    vec.fill(0);
  }
  for (auto &vec : QrawOneSide)
  {
    std::fill(vec.begin(), vec.end(), 0);
  }
}
