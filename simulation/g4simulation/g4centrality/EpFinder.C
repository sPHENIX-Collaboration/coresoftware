#include "EpFinder.h"
#include "TVector3.h"
#include "TH2D.h"
#include "TFile.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TClonesArray.h"
#include "TMath.h"
#include <cassert>
#include <iostream>
#include <vector>

#include <fun4all/PHTFileServer.h>

using namespace std;


EpFinder::EpFinder(int nEventTypeBins, char const* OutFileName, char const* CorrectionFile, int pbinsx, int pbinsy) : mThresh(0.0), mMax(100.0)
{

/*
  cout << "\n**********\n*  Welcome to the Event Plane finder.\n"
       << "*  Please note that when you specify 'order' as an argument to a method,\n"
       << "*  then 1=first-order plane, 2=second-order plane, etc.\n"
       << "*  This code is currently configured to go up to order=" << _EpOrderMax << "\n"; 
  cout << "**********\n"; 
*/


  mNumberOfEventTypeBins = nEventTypeBins;

  // ----------------------------------- Stuff read from the "Correction File" -----------------------------------------                                                    
  // "Shift correction" histograms that we INPUT and apply here                                                                                                             
   gErrorIgnoreLevel = 5000;
   TFile *mCorrectionInputFile = new TFile(CorrectionFile,"READ");

  if (mCorrectionInputFile->IsZombie()) {
//   std::cout << "EPFinder: Error opening file with Correction Histograms" << std::endl;
//   std::cout << "EPFinder: no correction file provided, no weighting will be used" << std::endl;
    for (int order=1; order<_EpOrderMax+1; order++){
      mEpShiftInput_sin[order-1] = NULL;
      mEpShiftInput_cos[order-1] = NULL;
      mAveCosDeltaPsi[order-1] = NULL;
    }
    mPhiWeightInput = NULL;
  }
  else{
    for (int order=1; order<_EpOrderMax+1; order++){
      mEpShiftInput_sin[order-1] = (TProfile2D*)mCorrectionInputFile->Get(Form("EpShiftPsi%d_sin",order));
      mEpShiftInput_cos[order-1] = (TProfile2D*)mCorrectionInputFile->Get(Form("EpShiftPsi%d_cos",order));
    }
    mPhiWeightInput = (TH2D*)mCorrectionInputFile->Get(Form("PhiWeight"));
  }

  mCorrectionInputFile->Close();
  delete mCorrectionInputFile;

  // ----------------------------------- Stuff written to the "Correction File" -----------------------------------------                                                   
  // "Shift correction" histograms that we produce and OUTPUT                                                                                                               
  OutFileNameString = OutFileName;  
  //PHTFileServer::get().open(OutFileName, "RECREATE");

  for (int order=1; order<_EpOrderMax+1; order++){
    mEpShiftOutput_sin[order-1] = new TProfile2D(Form("EpShiftPsi%d_sin",order),Form("EpShiftPsi%d_sin",order),
						  _EpTermsMax,0.5,1.0*_EpTermsMax+.5,nEventTypeBins,-0.5,(double)nEventTypeBins-0.5,-1.0,1.0);
    mEpShiftOutput_cos[order-1] = new TProfile2D(Form("EpShiftPsi%d_cos",order),Form("EpShiftPsi%d_cos",order),
						  _EpTermsMax,0.5,1.0*_EpTermsMax+.5,nEventTypeBins,-0.5,(double)nEventTypeBins-0.5,-1.0,1.0);
  }
  // Phi weighting histograms based on requested binning
  if(pbinsx<=0) pbinsx = 1; 
  if(pbinsy<=0) pbinsy = 1; 
  
  // binning tuned for FEMC, to be made generic in future JGL 8/28/2019
  // bins are ix,iy 
//  mPhiWeightOutput   = new TH2D(Form("PhiWeight"),Form("Phi Weight"),pbinsx,-0.5,(pbinsx-0.5),pbinsy,-0.5,(pbinsy-0.5));
    
  // just for normalization. discard after use
//  mPhiAveraged       = new TH2D(Form("PhiAveraged"),Form("Phi Average"),pbinsx,-0.5,(pbinsx-0.5),pbinsy,-0.5,(pbinsy-0.5));
 
}

void EpFinder::Finish(){

//  mPhiWeightOutput->Divide(mPhiAveraged);
//  delete mPhiAveraged;

  //PHTFileServer::get().cd(OutFileNameString.data());
  //PHTFileServer::get().write(OutFileNameString.data());

  //cout << "EpFinder is finished!\n\n";
}

//==================================================================================================================
EpInfo EpFinder::Results(std::vector<EpHit> *EpdHits, int EventTypeId){

  if ((EventTypeId<0)||(EventTypeId>=mNumberOfEventTypeBins)){
    cout << "You are asking for an undefined EventType - fail!\n";
    assert(0);
  }

  EpInfo result;

  double pi = M_PI;

  // This below is for normalizing Q-vectors
  double TotalWeight4Side[_EpOrderMax][2];       // for normalizing Q-vector: order, (nonPhiWeighted or PhiWeighted)  ** depends on Order because eta-weight depends on order
  for (int phiWeightedOrNo=0; phiWeightedOrNo<2; phiWeightedOrNo++){
    for (int order=1; order<_EpOrderMax+1; order++){
      TotalWeight4Side[order-1][phiWeightedOrNo] = 0;
    }
  }

  //--------------------------------- begin loop over hits ---------------------------------------
  for (unsigned int hit=0; hit<EpdHits->size(); hit++){
        
    float nMip = EpdHits->at(hit).nMip;
    if (nMip<mThresh) continue;
    double TileWeight = (nMip<mMax)?nMip:mMax;
    int idx_x =  EpdHits->at(hit).ix;
    int idx_y =  EpdHits->at(hit).iy;
    double phi = EpdHits->at(hit).phi;

    //---------------------------------
    // fill Phi Weight histograms to be used in next iteration (if desired)
    // Obviously, do this BEFORE phi weighting!
    //---------------------------------

  /*
    if(EpdHits->at(hit).samePhi){
//      mPhiWeightOutput->Fill(idx_x,idx_y,TileWeight);
      for(unsigned int i = 0; i<EpdHits->at(hit).samePhi->size(); i++){
	float x = EpdHits->at(hit).samePhi->at(i).first; 
	float y = EpdHits->at(hit).samePhi->at(i).second; 
//	mPhiAveraged->Fill(x,y,TileWeight/EpdHits->at(hit).samePhi->size());
      }
    }
    */
    
    //--------------------------------
    // now calculate Q-vectors
    //--------------------------------

    double PhiWeightedTileWeight = TileWeight;
    if (mPhiWeightInput) PhiWeightedTileWeight /= mPhiWeightInput->GetBinContent(idx_x+1,idx_y+1); 

    for (int order=1; order<_EpOrderMax+1; order++){
      double etaWeight = 1.0; // not implemented - JGL 8/27/2019
      TotalWeight4Side[order-1][0] += fabs(etaWeight) * TileWeight;             // yes the fabs() makes sense.  The sign in the eta weight is equivalent to a trigonometric phase.
      TotalWeight4Side[order-1][1] += fabs(etaWeight) * PhiWeightedTileWeight;  // yes the fabs() makes sense.  The sign in the eta weight is equivalent to a trigonometric phase.

      double Cosine = cos(phi*(double)order);
      double Sine   = sin(phi*(double)order);

      result.QrawOneSide[order-1][0]      += etaWeight * TileWeight * Cosine;
      result.QrawOneSide[order-1][1]      += etaWeight * TileWeight * Sine;

      result.QphiWeightedOneSide[order-1][0]      += etaWeight * PhiWeightedTileWeight * Cosine;
      result.QphiWeightedOneSide[order-1][1]      += etaWeight * PhiWeightedTileWeight * Sine;

    }
  }  // loop over hits

  //--------------------------------- end loop over hits ---------------------------------------

  // Weights used, so you can "un-normalize" the ring-by-ring Q-vectors.  
  for (int order=1; order<_EpOrderMax+1; order++){
    result.WheelSumWeightsRaw[order-1]         = TotalWeight4Side[order-1][0];
    result.WheelSumWeightsPhiWeighted[order-1] = TotalWeight4Side[order-1][1];
  }

  // at this point, we are finished with the detector hits, and deal only with the Q-vectors,

  // first, normalize the Q's...
  for (int order=1; order<_EpOrderMax+1; order++){
    if (TotalWeight4Side[order-1][0]>0.0001){
      result.QrawOneSide[order-1][0] /= TotalWeight4Side[order-1][0];
      result.QrawOneSide[order-1][1] /= TotalWeight4Side[order-1][0];
    }
    if (TotalWeight4Side[order-1][1]>0.0001){
      result.QphiWeightedOneSide[order-1][0] /= TotalWeight4Side[order-1][1];
      result.QphiWeightedOneSide[order-1][1] /= TotalWeight4Side[order-1][1];
    }
  }

  // at this point, we are finished with the Q-vectors and just use them to get angles Psi

  //---------------------------------
  // Calculate unshifted EP angles
  //---------------------------------
  for (int order=1; order<_EpOrderMax+1; order++){
    result.PsiRaw[order-1]                       = GetPsiInRange(result.QrawOneSide[order-1][0],result.QrawOneSide[order-1][1],order);
    result.PsiPhiWeighted[order-1]               = GetPsiInRange(result.QphiWeightedOneSide[order-1][0],result.QphiWeightedOneSide[order-1][1],order);
  } // loop over order

  //---------------------------------
  // Now shift
  //---------------------------------
  for (int order=1; order<_EpOrderMax+1; order++){
    result.PsiPhiWeightedAndShifted[order-1] = result.PsiPhiWeighted[order-1];
    if (mEpShiftInput_sin[order-1]) {
      for (int i=1; i<=_EpTermsMax; i++){
	double tmp = (double)(order*i);
	double sinAve = mEpShiftInput_sin[order-1]->GetBinContent(i,EventTypeId+1);    /// note the "+1" since EventTypeId begins at zero
	double cosAve = mEpShiftInput_cos[order-1]->GetBinContent(i,EventTypeId+1);    /// note the "+1" since EventTypeId begins at zero
	result.PsiPhiWeightedAndShifted[order-1] +=
	  2.0*(cosAve*sin(tmp*result.PsiPhiWeighted[order-1]) - sinAve*cos(tmp*result.PsiPhiWeighted[order-1]))/tmp;
      }
      double AngleWrapAround = 2.0*pi/(double)order;
      if (result.PsiPhiWeightedAndShifted[order-1]<0) result.PsiPhiWeightedAndShifted[order-1] += AngleWrapAround;
      else if (result.PsiPhiWeightedAndShifted[order-1]>AngleWrapAround) result.PsiPhiWeightedAndShifted[order-1] -= AngleWrapAround;
    }
  }

  //---------------------------------
  // Now calculate shift histograms for a FUTURE run (if you want it)
  //---------------------------------
  for (int i=1; i<=_EpTermsMax; i++){
    for (int order=1; order<_EpOrderMax+1; order++){
      double tmp = (double)(order*i);
      mEpShiftOutput_sin[order-1]->Fill(i,EventTypeId,sin(tmp*result.PsiPhiWeighted[order-1]));
      mEpShiftOutput_cos[order-1]->Fill(i,EventTypeId,cos(tmp*result.PsiPhiWeighted[order-1]));
    }
  }

  return result;
}

double EpFinder::GetPsiInRange(double Qx, double Qy, int order){
  double temp;
  if ((Qx==0.0)||(Qy==0.0)) temp=-999;
  else{
    temp = atan2(Qy,Qx)/((double)order);
    double AngleWrapAround = 2.0*M_PI/(double)order;
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

