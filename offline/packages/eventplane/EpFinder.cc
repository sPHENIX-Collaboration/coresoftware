#include "EpFinder.h"
#include "TVector3.h"
#include <cassert>
#include <iostream>
#include <vector>

#include <fun4all/PHTFileServer.h>

using namespace std;


EpFinder::EpFinder(int nEventTypeBins) : mThresh(0.2), mMax(15.0)
{
   mNumberOfEventTypeBins = nEventTypeBins;
}

void EpFinder::Finish()
{
  cout << "EpFinder is finished!\n\n";
}

//==================================================================================================================
EpInfo EpFinder::Results(std::vector<EpHit> *EpdHits, int EventTypeId)
 {

  if ((EventTypeId<0)||(EventTypeId>=mNumberOfEventTypeBins)){
    cout << "You are asking for an undefined EventType - fail!\n";
    assert(0);
 }

  EpInfo result;

    
  // For normalizing Q-vectors
  //phi weighting not implemented --Ejiro

  double TotalWeight4Side[_EpOrderMax][1];
    
  for (int phiWeightedOrNo=0; phiWeightedOrNo<1; phiWeightedOrNo++)
  {
    for (int order=1; order<_EpOrderMax+1; order++)
    {
      TotalWeight4Side[order-1][phiWeightedOrNo] = 0;
    }
  }

    //--------------------------------
    // begin loop over hits
    //--------------------------------
    
  for (unsigned int hit=0; hit<EpdHits->size(); hit++){
        
    float nMip = EpdHits->at(hit).nMip;
    if (nMip<mThresh) continue;
    double TileWeight = (nMip<mMax)?nMip:mMax;
    double phi = EpdHits->at(hit).phi;

    //--------------------------------
    // now calculate Q-vectors
    //--------------------------------


    for (int order=1; order<_EpOrderMax+1; order++)
    {
        
      double etaWeight = 1.0; // not implemented - JGL 8/27/2019
      
      //yes the fabs() makes sense.  The sign in the eta weight is equivalent to a trigonometric phase. --JGL
      TotalWeight4Side[order-1][0] += fabs(etaWeight) * TileWeight;
        

      double Cosine = cos(phi*(double)order);
      double Sine   = sin(phi*(double)order);

      result.QrawOneSide[order-1][0]      += etaWeight * TileWeight * Cosine;
      result.QrawOneSide[order-1][1]      += etaWeight * TileWeight * Sine;

    }
  }  // loop over hits


    
  // Weights used, so you can "un-normalize" the ring-by-ring Q-vectors.
  for (int order=1; order<_EpOrderMax+1; order++){
    result.WheelSumWeightsRaw[order-1]         = TotalWeight4Side[order-1][0];
  }

  // at this point, we are finished with the detector hits, and deal only with the Q-vectors,

  // first, normalize the Q's...
  for (int order=1; order<_EpOrderMax+1; order++){
    if (TotalWeight4Side[order-1][0]>0.0001){
      result.QrawOneSide[order-1][0] /= TotalWeight4Side[order-1][0];
      result.QrawOneSide[order-1][1] /= TotalWeight4Side[order-1][0];
    }
  }

  // at this point, we are finished with the Q-vectors and just use them to get angles Psi

  //---------------------------------
  // Calculate unshifted EP angles
  //---------------------------------
  for (int order=1; order<_EpOrderMax+1; order++)
  {
    result.PsiRaw[order-1]  = GetPsiInRange(result.QrawOneSide[order-1][0],result.QrawOneSide[order-1][1],order);
  } // loop over order

  return result;
}

double EpFinder::GetPsiInRange(double Qx, double Qy, int order){
  double temp;
  if ((Qx==0.0)||(Qy==0.0)) temp=-999;
  else{
    temp = atan2(Qy,Qx)/((double)order);
    double AngleWrapAround = (2.0*M_PI)/(double)order;
    if (temp<0.0) temp+= AngleWrapAround;
    else if (temp>AngleWrapAround) temp -= AngleWrapAround;
  }
  return temp;
}

bool EpFinder::OrderOutsideRange(int order){
  if (order < 1) {
    cout << "\n*** Invalid order specified ***\n";
    cout << "order must be 1 (for first-order plane) or higher\n";
    return true;
  }
  if (order > _EpOrderMax) {
    cout << "\n*** Invalid order specified ***\n";
    cout << "Maximum order=" << _EpOrderMax << ". To change, edit StEpdUtil/StEpdEpInfo.h\n";
    return true;
  }
  return false;
}

TString EpFinder::Report(){
  TString rep = Form("This is the EpFinder Report:\n");
  rep += Form("Number of EventType bins = %d\n",mNumberOfEventTypeBins);
  rep += Form("Threshold (in MipMPV units) = %f  and MAX weight = %f\n",mThresh,mMax);
  return rep;
}


