#include "EpFinder.h"

#include "EpInfo.h"  // for EpInfo

#include "TH2D.h"
#include "TFile.h"
#include <algorithm>  // for fill
#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>
#include <fstream>

#include <cdbobjects/CDBHistos.h>

EpFinder::EpFinder(int nEventTypeBins, unsigned int maxorder, int ns)
  : mNumberOfEventTypeBins(nEventTypeBins)
  , m_MaxOrder(maxorder)
{
    TotalWeight4Side.resize(m_MaxOrder);
    WheelSumWeightsRaw.resize(m_MaxOrder);
    WheelSumWeightsPhiWeighted.resize(m_MaxOrder);
    PsiRaw.resize(m_MaxOrder);
    PsiPhiWeighted.resize(m_MaxOrder);
    QrawOneSide.resize(m_MaxOrder);
    QphiWeightedOneSide.resize(m_MaxOrder);
    
    for (auto &vec : QrawOneSide)
    {
        vec.resize(2);
    }
    
    for (auto &vec : QphiWeightedOneSide)
    {
        vec.resize(2);
    }
    
    //---------------------------------------------------
    int pbinsx = 16; int pbinsy = 24; 
    mPhiWeightOutput[ns]   = new TH2D(Form("PhiWeightNS%d",ns),Form("sEPD Tile Weight divided by Averaged NS=%d",ns),pbinsx,-0.5,(pbinsx-0.5),pbinsy,-0.5,(pbinsy-0.5));
    mPhiAveraged[ns]       = new TH2D(Form("PhiAveraged%d",ns),Form("sEPD PhiAverage NS=%d",ns),pbinsx,-0.5,(pbinsx-0.5),pbinsy,-0.5,(pbinsy-0.5));
    mPhiWeightInput[ns]    = new TH2D(Form("PhiAveragedIn%d",ns),Form("sEPD PhiAverage In NS=%d",ns),pbinsx,-0.5,(pbinsx-0.5),pbinsy,-0.5,(pbinsy-0.5));

    CDBHistos *cdbhistosIn = new CDBHistos(Form("sEPD_PhiWeights_%d_Input.root", ns));
    cdbhistosIn->LoadCalibrations();
    if(!(cdbhistosIn->getHisto(Form("PhiWeightNS%d", ns))))
    {
        mPhiWeightInput[ns] = NULL;
    }
    else
    {
        //cloning mPhiWeightOutput[ns] since can't simply do mPhiWeightInput[ns] = cdbhistosIn->getHisto(Form("PhiWeightNS%d", ns)); CDBhistos TH1 restriction =(
        for(int etabin = 1; etabin <= pbinsx; etabin ++)
        {
            for(int phibin = 1; phibin <= pbinsy; phibin ++)
            {
                double histcontent = cdbhistosIn->getHisto(Form("PhiWeightNS%d", ns))->GetBinContent(etabin,phibin);
                mPhiWeightInput[ns]->SetBinContent(etabin, phibin, histcontent);
            }
        }
    }
   
   cdbhistosIn->Print(); 
   
}
void EpFinder::Finish(int thishist)
{

    mPhiWeightOutput[thishist]->Divide(mPhiAveraged[thishist]);
    delete mPhiAveraged[thishist];
    
    CDBHistos *cdbhistosOut = new CDBHistos(Form("sEPD_PhiWeights_%d_Output.root", thishist));
    cdbhistosOut->registerHisto(mPhiWeightOutput[thishist]);
    cdbhistosOut->WriteCDBHistos();
     
    std::cout << "EpFinder is finished!"<<std::endl;
}

//==================================================================================================================
void EpFinder::ResultsEPD(const std::vector<EpHit> &EpdHits, int EventTypeId, EpInfo *epinfo)
{
  if ((EventTypeId < 0) || (EventTypeId >= mNumberOfEventTypeBins))
  {
    std::cout << "You are asking for an undefined EventType - fail!" << std::endl;
    assert(0);
  }

  //--------------------------------
  // begin loop over hits
  //--------------------------------

  for (unsigned int hit = 0; hit < EpdHits.size(); hit++)
  {
    float nMip = EpdHits.at(hit).nMip;
    double TileWeight = nMip;
    double phi = EpdHits.at(hit).phi;
    int idx_x =  EpdHits.at(hit).ix;
    int idx_y =  EpdHits.at(hit).iy;
    int _ns =  EpdHits.at(hit).wheel;

      
    if(EpdHits.at(hit).sameRing){
    mPhiWeightOutput[_ns]->Fill(idx_x,idx_y,TileWeight);
    for(unsigned int i = 0; i<EpdHits.at(hit).sameRing->size(); i++){
    float x = EpdHits.at(hit).sameRing->at(i).first;
    float y = EpdHits.at(hit).sameRing->at(i).second;
    mPhiAveraged[_ns]->Fill(x,y,TileWeight/EpdHits.at(hit).sameRing->size());
     }
    }
      

    //--------------------------------
    // now calculate Q-vectors
    //--------------------------------
    
    double PhiWeightedTileWeight = TileWeight;
      if (mPhiWeightInput[_ns]) PhiWeightedTileWeight /= mPhiWeightInput[_ns]->GetBinContent(idx_x+1,idx_y+1);
    
    for (unsigned int order = 1; order < m_MaxOrder + 1; order++)
    {
     
      double etaWeight = 1.0;

      TotalWeight4Side[order - 1][0] += fabs(etaWeight) * TileWeight;
      TotalWeight4Side[order-1][1] += fabs(etaWeight) * PhiWeightedTileWeight;

      double Cosine = cos(phi * (double) order);
      double Sine = sin(phi * (double) order);
      QrawOneSide[order - 1][0] += etaWeight * TileWeight * Cosine;
      QrawOneSide[order - 1][1] += etaWeight * TileWeight * Sine;
        
      QphiWeightedOneSide[order - 1][0] += etaWeight * PhiWeightedTileWeight * Cosine;
      QphiWeightedOneSide[order - 1][1] += etaWeight * PhiWeightedTileWeight * Sine;
        
    }
  }  // loop over hits

  // Weights used, so you can "un-normalize" the ring-by-ring Q-vectors.
  for (unsigned int order = 1; order < m_MaxOrder + 1; order++)
  {
    WheelSumWeightsRaw[order - 1] = TotalWeight4Side[order - 1][0];
    WheelSumWeightsPhiWeighted[order - 1] = TotalWeight4Side[order - 1][1];

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

    if (TotalWeight4Side[order-1][1]>0.0001){
        QphiWeightedOneSide[order-1][0] /= TotalWeight4Side[order-1][1];
        QphiWeightedOneSide[order-1][1] /= TotalWeight4Side[order-1][1];
    }

  }

  // at this point, we are finished with the Q-vectors and just use them to get angles Psi

  //---------------------------------
  // Calculate unshifted EP angles
  //---------------------------------
  for (unsigned int order = 1; order < m_MaxOrder + 1; order++)
  {
    PsiRaw[order - 1] = GetPsiInRange(QrawOneSide[order - 1][0], QrawOneSide[order - 1][1], order);
    PsiPhiWeighted[order-1] = GetPsiInRange(QphiWeightedOneSide[order-1][0],QphiWeightedOneSide[order-1][1],order);
  } // loop over order


  // copy results to i/o object
  if (epinfo)
  {
    epinfo->CopyQrawOneSide(QrawOneSide);
    epinfo->CopyQphiWeightedOneSide(QphiWeightedOneSide);
    epinfo->CopyWheelSumWeightsRaw(WheelSumWeightsRaw);
    epinfo->CopyWheelSumWeightsPhiWeighted(WheelSumWeightsPhiWeighted);
    epinfo->CopyPsiRaw(PsiRaw);
    epinfo->CopyPsiPhiWeighted(PsiPhiWeighted);
  }

  return;
}

double EpFinder::GetPsiInRange(double Qx, double Qy, unsigned int order) const
{
  double temp;
  if ((Qx == 0.0) || (Qy == 0.0))
    temp = NAN;
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
    
 for (auto &vec : QphiWeightedOneSide)
 {
    std::fill(vec.begin(), vec.end(), 0);
 }
    
}
