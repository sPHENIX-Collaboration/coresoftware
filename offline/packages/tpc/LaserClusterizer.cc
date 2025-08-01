#include "LaserClusterizer.h"

#include "LaserEventInfo.h"

#include <trackbase/LaserCluster.h>
#include <trackbase/LaserClusterContainer.h>
#include <trackbase/LaserClusterContainerv1.h>
#include <trackbase/LaserClusterv2.h>
#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>  // for hitkey, getLayer
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <ffaobjects/EventHeader.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>  // for PHIODataNode
#include <phool/PHNode.h>        // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/PHTimer.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <Acts/Definitions/Units.hpp>
#include <Acts/Surfaces/Surface.hpp>

#include <TF1.h>
#include <TFile.h>
#include <TH3.h>
//#include <TF3.h>
//#include <ROOT/Fit/Fitter.h>
#include <Fit/Fitter.h>
//#include <Math/Functor.h>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <array>
#include <cmath>  // for sqrt, cos, sin
#include <iostream>
#include <limits>
#include <map>  // for _Rb_tree_cons...
#include <string>
#include <utility>  // for pair
#include <vector>
#include <set>
#include <queue>
#include <tuple>

#include <pthread.h>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

using point = bg::model::point<float, 3, bg::cs::cartesian>;
using box = bg::model::box<point>;
using specHitKey = std::pair<TrkrDefs::hitkey, TrkrDefs::hitsetkey>;
using adcKey = std::pair<unsigned int, specHitKey>;
using pointKeyLaser = std::pair<point, specHitKey>;
using hitData = std::pair<point, adcKey>;


int layerMins[3] = {7,23,39};
int layerMaxes[3] = {22, 38, 54};



namespace
{
  struct thread_data
  {
    PHG4TpcCylinderGeomContainer *geom_container = nullptr;
    ActsGeometry *tGeometry = nullptr;
    std::vector<TrkrHitSet *> hitsets;
    std::vector<unsigned int> layers;
    bool side = false;
    unsigned int sector = 0;
    unsigned int module = 0;
    std::vector<LaserCluster *> cluster_vector;
    std::vector<TrkrDefs::cluskey> cluster_key_vector;
    double adc_threshold = 74.4;
    int peakTimeBin = 325;
    int layerMin = 1;
    int layerMax = 1;
    double tdriftmax = 0;
    int eventNum = 0;
    int Verbosity = 0;
    TH3D *hitHist = nullptr;
    bool doFitting = false;
  };

  pthread_mutex_t mythreadlock;

  const std::vector<point> neighborOffsets = {
    point(1, 0, 0), point(-1, 0, 0),
    point(0, 1, 0), point(0, -1, 0),
    point(0, 0, 1), point(0, 0, -1),
    point(0, 0, 2), point(0, 0, -2)
  };


  double layerFunction(double *x, double *par)
  {
    double A = par[0];
    double mu = par[1];


    double binCenter = round(x[0]);
    double overlapLow = std::max(binCenter - 0.5, mu - 0.5);
    double overlapHigh = std::min(binCenter + 0.5, mu + 0.5);
    double overlap = overlapHigh - overlapLow;
    if(overlap <= 0.0)
    {
      return 0.0;
    }
    return A*overlap;
    /*
    if(fabs(x[0] - mu) < 1)
    {
      double frac = 1-fabs(mu - round(x[0]));
      return A*(frac);
    }
    return 0.0;
    */


  }

  double phiFunction(double *x, double *par)
  {
    if(par[2] < 0.0)
    {
      return 0.0;
    }
    return par[0] * TMath::Gaus(x[0],par[1],par[2],false);
  }

  double timeFunction(double *x, double *par)
  {
    if(par[2] < 0.0)
    {
      return 0.0;
    }
    double g = TMath::Gaus(x[0],par[1],par[2],true);
    double cdf = 1 + TMath::Erfc(par[3]*(x[0]-par[1])/(sqrt(2.0)*par[2]));
    return par[0]*g*cdf;
  }

  void findConnectedRegions3(std::vector<hitData> &clusHits, std::pair<TrkrDefs::hitkey, TrkrDefs::hitsetkey> &maxKey)
  {
    std::vector<std::vector<hitData>> regions;

    std::vector<hitData> unvisited;
    for(auto &clusHit : clusHits)
    {
      unvisited.push_back(clusHit);
    }

    while(!unvisited.empty())
    {
      std::vector<hitData> region;
      std::queue<hitData> q;

      unsigned int mIndex = 0;
      int i=0;
      for(auto hit : unvisited)
      {
        if(hit.second.second.first == maxKey.first && hit.second.second.second == maxKey.second)
        {
          mIndex = i;
          break;
        }
        i++;
      }

      auto seed = unvisited[mIndex];
      unvisited.erase(unvisited.begin()+mIndex);
      q.push(seed);
      region.push_back(seed);

      while(!q.empty())
      {
        float ix = q.front().first.get<0>();
        float iy = q.front().first.get<1>();
        float iz = q.front().first.get<2>();
        q.pop();

        for (auto neigh : neighborOffsets)
        {
          float nx = ix + neigh.get<0>();
          float ny = iy + neigh.get<1>();
          float nz = iz + neigh.get<2>();

          for(unsigned int v=0; v<unvisited.size(); v++)
          {

            if(fabs(unvisited[v].first.get<0>() - nx) < 0.01 && fabs(unvisited[v].first.get<1>() - ny) < 0.01 && fabs(unvisited[v].first.get<2>() - nz) < 0.01)
            {
              auto newSeed = unvisited[v];
              unvisited.erase(unvisited.begin()+v);
              q.push(newSeed);
              region.push_back(newSeed);
              break;
            }
          }
        }
      }
      regions.push_back(region);

    }

    clusHits.clear();
    for(auto hit : regions[0])
    {
      clusHits.push_back(hit);
    }

  }

   
  void remove_hits(std::vector<hitData> &clusHits, bgi::rtree<hitData, bgi::quadratic<16>> &rtree, std::multimap<unsigned int, pointKeyLaser> &adcMap)
  {
    for (auto &clusHit : clusHits)
    {
      auto spechitkey = clusHit.second.second;

      rtree.remove(clusHit);

      for (auto iterAdc = adcMap.begin(); iterAdc != adcMap.end();)
      {
	      if(iterAdc->second.second == spechitkey)
	      {
	        iterAdc = adcMap.erase(iterAdc);
	        break;
	      }
	      else
	      {
	        ++iterAdc;
	      }
      }
    }

  }

  void calc_cluster_parameter(std::vector<hitData> &clusHits, thread_data &my_data, std::pair<TrkrDefs::hitkey, TrkrDefs::hitsetkey> maxADCKey)
  {



    findConnectedRegions3(clusHits, maxADCKey);


    double rSum = 0.0;
    double phiSum = 0.0;
    double tSum = 0.0;
    
    double layerSum = 0.0;
    double iphiSum = 0.0;
    double itSum = 0.0;
    
    double adcSum = 0.0;
    
    double maxAdc = 0.0;
    TrkrDefs::hitsetkey maxKey = 0;
    
    unsigned int nHits = clusHits.size();
    
    auto *clus = new LaserClusterv2;
    
    int meanSide = 0;
    
    std::vector<float> usedLayer;
    std::vector<float> usedIPhi;
    std::vector<float> usedIT;
    
    double meanLayer = 0.0;
    double meanIPhi = 0.0;
    double meanIT = 0.0;

    for (auto &clusHit : clusHits)
    {
      float coords[3] = {clusHit.first.get<0>(), clusHit.first.get<1>(), clusHit.first.get<2>()};
      std::pair<TrkrDefs::hitkey, TrkrDefs::hitsetkey> spechitkey = clusHit.second.second;
      unsigned int adc = clusHit.second.first;

      int side = TpcDefs::getSide(spechitkey.second);
          
      if (side)
      {
	      meanSide++;
      }
      else
      {
	      meanSide--;
      }
      
      PHG4TpcCylinderGeom *layergeom = my_data.geom_container->GetLayerCellGeom((int) coords[0]);
      
      double r = layergeom->get_radius();
      double phi = layergeom->get_phi(coords[1], side);
      double t = layergeom->get_zcenter(fabs(coords[2]));
      
      double hitzdriftlength = t * my_data.tGeometry->get_drift_velocity();
      double hitZ = my_data.tdriftmax * my_data.tGeometry->get_drift_velocity() - hitzdriftlength;
      
	  
	    bool foundLayer = false;
	    for (float i : usedLayer)
	    {
	      if (coords[0] == i)
	      {
	        foundLayer = true;
	        break;
	      }
	    }
	  
	    if (!foundLayer)
	    {
	      usedLayer.push_back(coords[0]);
	    }
	  
	    bool foundIPhi = false;
	    for (float i : usedIPhi)
	    {
	      if (coords[1] == i)
	      {
	        foundIPhi = true;
	        break;
	      }
	    }
	  
	    if (!foundIPhi)
	    {
	      usedIPhi.push_back(coords[1]);
	    }
	  
	    bool foundIT = false;
	    for (float i : usedIT)
	    {
	      if (coords[2] == i)
	      {
	        foundIT = true;
	        break;
	      }
	    }
	  
	    if (!foundIT)
	    {
	      usedIT.push_back(coords[2]);
	    }
	  
	    clus->addHit();
	    clus->setHitLayer(clus->getNhits() - 1, coords[0]);
	    clus->setHitIPhi(clus->getNhits() - 1, coords[1]);
	    clus->setHitIT(clus->getNhits() - 1, coords[2]);
	    clus->setHitX(clus->getNhits() - 1, r * cos(phi));
	    clus->setHitY(clus->getNhits() - 1, r * sin(phi));
	    clus->setHitZ(clus->getNhits() - 1, hitZ);
	    clus->setHitAdc(clus->getNhits() - 1, (float) adc);
	  
	    rSum += r * adc;
	    phiSum += phi * adc;
	    tSum += t * adc;
	  
	    layerSum += coords[0] * adc;
	    iphiSum += coords[1] * adc;
	    itSum += coords[2] * adc;
	  
	    meanLayer += coords[0];
	    meanIPhi += coords[1];
	    meanIT += coords[2];
	  
	    adcSum += adc;
	  
	    if (adc > maxAdc)
	    {
	      maxAdc = adc;
	      maxKey = spechitkey.second;
	    }
	  
	  }
    
    if (nHits == 0)
    {
      return;
    }

    
    double clusR = rSum / adcSum;
    double clusPhi = phiSum / adcSum;
    double clusT = tSum / adcSum;
    double zdriftlength = clusT * my_data.tGeometry->get_drift_velocity();
    
    double clusX = clusR * cos(clusPhi);
    double clusY = clusR * sin(clusPhi);
    double clusZ = my_data.tdriftmax * my_data.tGeometry->get_drift_velocity() - zdriftlength;
    if (meanSide < 0)
    {
      clusZ = -clusZ;
      for (int i = 0; i < (int) clus->getNhits(); i++)
      {
	      clus->setHitZ(i, -1 * clus->getHitZ(i));
      }
    }
    
    std::sort(usedLayer.begin(), usedLayer.end());
    std::sort(usedIPhi.begin(), usedIPhi.end());
    std::sort(usedIT.begin(), usedIT.end());

    meanLayer = meanLayer / nHits;
    meanIPhi = meanIPhi / nHits;
    meanIT = meanIT / nHits;
    
    double sigmaLayer = 0.0;
    double sigmaIPhi = 0.0;
    double sigmaIT = 0.0;
    
    double sigmaWeightedLayer = 0.0;
    double sigmaWeightedIPhi = 0.0;
    double sigmaWeightedIT = 0.0;
    
    pthread_mutex_lock(&mythreadlock);
    my_data.hitHist = new TH3D(Form("hitHist_event%d_side%d_sector%d_module%d_cluster%d",my_data.eventNum,(int)my_data.side,(int)my_data.sector,(int)my_data.module,(int)my_data.cluster_vector.size()),";layer;iphi;it",usedLayer.size()+2,usedLayer[0]-1.5,*usedLayer.rbegin()+1.5,usedIPhi.size()+2,usedIPhi[0]-1.5,*usedIPhi.rbegin()+1.5,usedIT.size()+2,usedIT[0]-1.5,*usedIT.rbegin()+1.5);


    //TH3D *hitHist = new TH3D(Form("hitHist_event%d_side%d_sector%d_module%d_cluster%d",my_data.eventNum,(int)my_data.side,(int)my_data.sector,(int)my_data.module,(int)my_data.cluster_vector.size()),";layer;iphi;it",usedLayer.size()+2,usedLayer[0]-1.5,*usedLayer.rbegin()+1.5,usedIPhi.size()+2,usedIPhi[0]-1.5,*usedIPhi.rbegin()+1.5,usedIT.size()+2,usedIT[0]-1.5,*usedIT.rbegin()+1.5);

    for (int i = 0; i < (int) clus->getNhits(); i++)
    {

      my_data.hitHist->Fill(clus->getHitLayer(i), clus->getHitIPhi(i), clus->getHitIT(i), clus->getHitAdc(i));

      sigmaLayer += pow(clus->getHitLayer(i) - meanLayer, 2);
      sigmaIPhi += pow(clus->getHitIPhi(i) - meanIPhi, 2);
      sigmaIT += pow(clus->getHitIT(i) - meanIT, 2);
      
      sigmaWeightedLayer += clus->getHitAdc(i) * pow(clus->getHitLayer(i) - (layerSum / adcSum), 2);
      sigmaWeightedIPhi += clus->getHitAdc(i) * pow(clus->getHitIPhi(i) - (iphiSum / adcSum), 2);
      sigmaWeightedIT += clus->getHitAdc(i) * pow(clus->getHitIT(i) - (itSum / adcSum), 2);
    }

    bool fitSuccess = false;
    ROOT::Fit::Fitter *fit3D = new ROOT::Fit::Fitter;

    if(my_data.doFitting)
    {

      double par_init[7] = {
        maxAdc,
        meanLayer,
        meanIPhi, 0.75,
        meanIT, 0.5, 1
      };

      double satThreshold = 900.0;
      double sigma_ADC = 20.0;

      auto nll = [&](const double* par)
      {
        double nll_val = 0.0;

        int nx = my_data.hitHist->GetNbinsX();
        int ny = my_data.hitHist->GetNbinsY();
        int nz = my_data.hitHist->GetNbinsZ();

        double parLayer[2] = {1.0,par[1]};
        double parPhi[4] = {1.0,par[2],par[3]};
        double parTime[4] = {1.0,par[4],par[5],par[6]};

        double xyz[3];

        for (int i = 1; i <= nx; ++i)
        {
          xyz[0] = my_data.hitHist->GetXaxis()->GetBinCenter(i);  
          for (int j = 1; j <= ny; ++j)
          {
            xyz[1] = my_data.hitHist->GetYaxis()->GetBinCenter(j);    
            for (int k = 1; k <= nz; ++k)
            {
              xyz[2] = my_data.hitHist->GetZaxis()->GetBinCenter(k);
              double observed = my_data.hitHist->GetBinContent(i, j, k);

              double expected = par[0]*layerFunction(&xyz[0], parLayer)*phiFunction(&xyz[1], parPhi)*timeFunction(&xyz[2], parTime);

              if(observed <= my_data.adc_threshold)
              {
                double arg = (expected - my_data.adc_threshold) / (sqrt(2.0) * sigma_ADC);
                double tail_prob = 0.5 * TMath::Erfc(arg);
                nll_val -= log(tail_prob + 1e-12);
              }
              else if(observed < satThreshold)
              {
                double resid = (observed - expected) / sigma_ADC;
                nll_val += 0.5 * (resid * resid + log(2 * TMath::Pi() * sigma_ADC * sigma_ADC));
              }
              else if(observed >= satThreshold)
              {
                double arg = (satThreshold - expected) / (sqrt(2.0) * sigma_ADC);
                double tail_prob = 0.5 * TMath::Erfc(arg);
                nll_val -= log(tail_prob + 1e-12);
              }
            }
          }
        }
        return nll_val;
      };

      fit3D->SetFCN(7, nll, par_init);

      fit3D->Config().ParSettings(0).SetName("amp");
      fit3D->Config().ParSettings(0).SetStepSize(10);
      fit3D->Config().ParSettings(0).SetLimits(0,5000);
      fit3D->Config().ParSettings(1).SetName("mu_layer");
      fit3D->Config().ParSettings(1).SetStepSize(0.1);
      fit3D->Config().ParSettings(1).SetLimits(usedLayer[0],*usedLayer.rbegin());
      fit3D->Config().ParSettings(2).SetName("mu_phi");
      fit3D->Config().ParSettings(2).SetStepSize(0.1);
      fit3D->Config().ParSettings(2).SetLimits(usedIPhi[0],*usedIPhi.rbegin());
      fit3D->Config().ParSettings(3).SetName("sig_phi");
      fit3D->Config().ParSettings(3).SetStepSize(0.1);
      fit3D->Config().ParSettings(3).SetLimits(0.01,2);
      fit3D->Config().ParSettings(4).SetName("mu_t");
      fit3D->Config().ParSettings(4).SetStepSize(0.1);
      fit3D->Config().ParSettings(4).SetLimits(usedIT[0],*usedIT.rbegin());
      fit3D->Config().ParSettings(5).SetName("sig_t");
      fit3D->Config().ParSettings(5).SetStepSize(0.1);
      fit3D->Config().ParSettings(5).SetLimits(0.01,10);
      fit3D->Config().ParSettings(6).SetName("lambda_t");
      fit3D->Config().ParSettings(6).SetStepSize(0.01);
      fit3D->Config().ParSettings(6).SetLimits(0,5);

 
      if(usedLayer.size() == 1)
      {
        fit3D->Config().ParSettings(1).Fix();
      }

      fitSuccess = fit3D->FitFCN();


      if (my_data.Verbosity > 2)
      {
        std::cout << "fit success: " << fitSuccess << std::endl;
      }
    }
    pthread_mutex_unlock(&mythreadlock);



    if(my_data.doFitting && fitSuccess)
    {

      const ROOT::Fit::FitResult& result = fit3D->Result();


      PHG4TpcCylinderGeom *layergeomLow = my_data.geom_container->GetLayerCellGeom((int) floor(result.Parameter(1)));
      PHG4TpcCylinderGeom *layergeomHigh = my_data.geom_container->GetLayerCellGeom((int) ceil(result.Parameter(1)));

      double RLow = layergeomLow->get_radius();
      double RHigh = layergeomHigh->get_radius();

      double phiHigh_RLow = -999.0;
      if(ceil(result.Parameter(2)) < layergeomLow->get_phibins())
      {
        phiHigh_RLow = layergeomLow->get_phi(ceil(result.Parameter(2)), (meanSide < 0 ? 0 : 1));
      }
      double phiHigh_RHigh = -999.0;
      if(ceil(result.Parameter(2)) < layergeomHigh->get_phibins())
      {
        phiHigh_RHigh = layergeomHigh->get_phi(ceil(result.Parameter(2)), (meanSide < 0 ? 0 : 1));
      }

      double phiLow_RLow = layergeomLow->get_phi(floor(result.Parameter(2)), (meanSide < 0 ? 0 : 1));
      double phiLow_RHigh = layergeomHigh->get_phi(floor(result.Parameter(2)), (meanSide < 0 ? 0 : 1));

      double meanR = (result.Parameter(1) - floor(result.Parameter(1))) * (RHigh - RLow) + RLow;

      double meanPhi_RLow = ((result.Parameter(2) - floor(result.Parameter(2)))) * (phiHigh_RLow - phiLow_RLow) + phiLow_RLow;
      double meanPhi_RHigh = ((result.Parameter(2) - floor(result.Parameter(2)))) * (phiHigh_RHigh - phiLow_RHigh) + phiLow_RHigh;

      double meanPhi = 0.5*(meanPhi_RLow + meanPhi_RHigh);
      if(phiHigh_RLow == -999.0 && phiHigh_RHigh != -999.0)
      {
        meanPhi = meanPhi_RHigh;
      }
      else if(phiHigh_RLow != -999.0 && phiHigh_RHigh == -999.0)
      {
        meanPhi = meanPhi_RLow;
      }
      

      if(phiHigh_RLow == -999.0 && phiHigh_RHigh == -999.0)
      {
        clus->setAdc(adcSum);
        clus->setX(clusX);
        clus->setY(clusY);
        clus->setZ(clusZ);
        clus->setFitMode(false);
        clus->setLayer(layerSum / adcSum);
        clus->setIPhi(iphiSum / adcSum);
        clus->setIT(itSum / adcSum);
        clus->setNLayers(usedLayer.size());
        clus->setNIPhi(usedIPhi.size());
        clus->setNIT(usedIT.size());
        clus->setSDLayer(sqrt(sigmaLayer / nHits));
        clus->setSDIPhi(sqrt(sigmaIPhi / nHits));
        clus->setSDIT(sqrt(sigmaIT / nHits));
        clus->setSDWeightedLayer(sqrt(sigmaWeightedLayer / adcSum));
        clus->setSDWeightedIPhi(sqrt(sigmaWeightedIPhi / adcSum));
        clus->setSDWeightedIT(sqrt(sigmaWeightedIT / adcSum));
      }
      else
      { 
        clus->setAdc(adcSum);
        clus->setX(meanR*cos(meanPhi));
        clus->setY(meanR*sin(meanPhi));
        clus->setZ(clusZ);
        clus->setFitMode(true);
        clus->setLayer(result.Parameter(1));
        clus->setIPhi(result.Parameter(2));
        clus->setIT(result.Parameter(4));
        clus->setNLayers(usedLayer.size());
        clus->setNIPhi(usedIPhi.size());
        clus->setNIT(usedIT.size());
        clus->setSDLayer(sqrt(sigmaLayer / nHits));
        clus->setSDIPhi(sqrt(sigmaIPhi / nHits));
        clus->setSDIT(sqrt(sigmaIT / nHits));
        clus->setSDWeightedLayer(sqrt(sigmaWeightedLayer / adcSum));
        clus->setSDWeightedIPhi(result.Parameter(3));
        clus->setSDWeightedIT(result.Parameter(5));
      }
    }
    else
    {
      clus->setAdc(adcSum);
      clus->setX(clusX);
      clus->setY(clusY);
      clus->setZ(clusZ);
      clus->setFitMode(false);
      clus->setLayer(layerSum / adcSum);
      clus->setIPhi(iphiSum / adcSum);
      clus->setIT(itSum / adcSum);
      clus->setNLayers(usedLayer.size());
      clus->setNIPhi(usedIPhi.size());
      clus->setNIT(usedIT.size());
      clus->setSDLayer(sqrt(sigmaLayer / nHits));
      clus->setSDIPhi(sqrt(sigmaIPhi / nHits));
      clus->setSDIT(sqrt(sigmaIT / nHits));
      clus->setSDWeightedLayer(sqrt(sigmaWeightedLayer / adcSum));
      clus->setSDWeightedIPhi(sqrt(sigmaWeightedIPhi / adcSum));
      clus->setSDWeightedIT(sqrt(sigmaWeightedIT / adcSum));
    }

    const auto ckey = TrkrDefs::genClusKey(maxKey, my_data.cluster_vector.size());
    my_data.cluster_vector.push_back(clus);
    my_data.cluster_key_vector.push_back(ckey);
    
    if(fit3D)
    {
      delete fit3D;
    }

    if(my_data.hitHist)
    {
      delete my_data.hitHist;
      my_data.hitHist = nullptr;
    }
    

  }
  

  void ProcessModuleData(thread_data *my_data)
  {
    
    if (my_data->Verbosity > 2)
    {
      pthread_mutex_lock(&mythreadlock);
      std::cout << "working on side: " << my_data->side << "   sector: " << my_data->sector << "   module: " << my_data->module << std::endl;
      pthread_mutex_unlock(&mythreadlock);
    }

    bgi::rtree<hitData, bgi::quadratic<16>> rtree;

    std::multimap<unsigned int, pointKeyLaser> adcMap;

    if (my_data->hitsets.size() == 0)
    {
      return;
    }
    for(int i=0; i<(int)my_data->hitsets.size(); i++)
    {
      auto *hitset = my_data->hitsets[i];
      unsigned int layer = my_data->layers[i];
      bool side = my_data->side;
      unsigned int sector = my_data->sector;

      TrkrDefs::hitsetkey hitsetKey = TpcDefs::genHitSetKey(layer, sector, (int)side);

      TrkrHitSet::ConstRange hitrangei = hitset->getHits();

      for (TrkrHitSet::ConstIterator hitr = hitrangei.first; hitr != hitrangei.second; ++hitr)
      {
	      float_t fadc = hitr->second->getAdc();
	      unsigned short adc = 0;
	      if (fadc > my_data->adc_threshold)
	      {
	        adc = (unsigned short) fadc;
	      }
	      else
	      {
	        continue;
	      }
	
	      int iphi = TpcDefs::getPad(hitr->first);
	      int it = TpcDefs::getTBin(hitr->first);
        
	      if(fabs(it - my_data->peakTimeBin) > 5)
	      {
	        continue;
	      }

	      point coords = point((int) layer, iphi, it);

	      std::vector<hitData> testduplicate;
	      rtree.query(bgi::intersects(box(point(layer - 0.001, iphi - 0.001, it - 0.001),
					point(layer + 0.001, iphi + 0.001, it + 0.001))),
		    std::back_inserter(testduplicate));
	      if (!testduplicate.empty())
	      {
	        testduplicate.clear();
	        continue;
	      }

	      TrkrDefs::hitkey hitKey = TpcDefs::genHitKey(iphi, it);
	
	      auto spechitkey = std::make_pair(hitKey, hitsetKey);
	      pointKeyLaser coordsKey = std::make_pair(coords, spechitkey);
	      adcMap.insert(std::make_pair(adc, coordsKey));
        auto adckey = std::make_pair(adc, spechitkey);
	      rtree.insert(std::make_pair(point(1.0*layer, 1.0*iphi, 1.0*it), adckey));
      }
    }
    //finished filling rtree

    while (adcMap.size() > 0)
    {
      auto iterKey = adcMap.rbegin();
      if(iterKey == adcMap.rend())
      {
	      break;
      }
      

      auto coords = iterKey->second.first;
      int layer = coords.get<0>();
      int iphi = coords.get<1>();
      int it = coords.get<2>();
   
      if (my_data->Verbosity > 2)
      {
        pthread_mutex_lock(&mythreadlock);
        std::cout << "working on cluster " << my_data->cluster_vector.size() << "   side: " << my_data->side << "   sector: " << my_data->sector << "   module: " << (layer<23 ? 1 : (layer<39 ? 2 : 3) ) << std::endl;
        pthread_mutex_unlock(&mythreadlock);

      }

      std::vector<hitData> clusHits;

      rtree.query(bgi::intersects(box(point(layer - my_data->layerMin, iphi - 6, it - 5), point(layer + my_data->layerMax, iphi + 6, it + 5))), std::back_inserter(clusHits));

      calc_cluster_parameter(clusHits, *my_data, iterKey->second.second);

      remove_hits(clusHits, rtree, adcMap);

    }
  }

  void *ProcessModule(void *threadarg)
  {
    auto my_data = static_cast<thread_data *>(threadarg);
    ProcessModuleData(my_data);
    pthread_exit(nullptr);
  }
} //namespace

LaserClusterizer::LaserClusterizer(const std::string &name)
  : SubsysReco(name)
{
}

int LaserClusterizer::InitRun(PHCompositeNode *topNode)
{
  TH1::AddDirectory(kFALSE);

  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // Create the Cluster node if required
  std::string laserClusterNodeName = "LASER_CLUSTER";
  if (m_lamination)
  {
    laserClusterNodeName = "LAMINATION_CLUSTER";
  }
  auto laserclusters = findNode::getClass<LaserClusterContainer>(dstNode, laserClusterNodeName);
  if (!laserclusters)
  {
    PHNodeIterator dstiter(dstNode);
    PHCompositeNode *DetNode =
        dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

    laserclusters = new LaserClusterContainerv1;
    PHIODataNode<PHObject> *LaserClusterContainerNode =
        new PHIODataNode<PHObject>(laserclusters, laserClusterNodeName, "PHObject");
    DetNode->addNode(LaserClusterContainerNode);
  }
  
  m_geom_container =
      findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!m_geom_container)
  {
    std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  // get the first layer to get the clock freq
  AdcClockPeriod = m_geom_container->GetFirstLayerCellGeom()->get_zstep();
  m_tdriftmax = AdcClockPeriod * NZBinsSide;

  return Fun4AllReturnCodes::EVENT_OK;
}

int LaserClusterizer::process_event(PHCompositeNode *topNode)
{
  eventHeader = findNode::getClass<EventHeader>(topNode, "EventHeader");
  if (!eventHeader)
  {
    std::cout << PHWHERE << " EventHeader Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_event = eventHeader->get_EvtSequence();

  if (Verbosity() > 1)
  {
    std::cout << "LaserClusterizer::process_event working on event " << m_event << std::endl;
  }

  m_laserEventInfo = findNode::getClass<LaserEventInfo>(topNode, "LaserEventInfo");
  if (!m_laserEventInfo)
  {
    std::cout << PHWHERE << "ERROR: Can't find node LaserEventInfo" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if(!m_laserEventInfo->isLaserEvent() && !m_laserEventInfo->isGl1LaserEvent() && !m_laserEventInfo->isGl1LaserPileupEvent())
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  if (Verbosity() > 1)
  {
    std::cout << "LaserClusterizer::process_event laser event found" << std::endl;
  }

  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  // get node containing the digitized hits
  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!m_hits)
  {
    std::cout << PHWHERE << "ERROR: Can't find node TRKR_HITSET" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  
  // get node for clusters
  std::string laserClusterNodeName = "LASER_CLUSTER";
  if (m_lamination)
  {
    laserClusterNodeName = "LAMINATION_CLUSTER";
  }
  m_clusterlist = findNode::getClass<LaserClusterContainer>(topNode, laserClusterNodeName);
  if (!m_clusterlist)
  {
    std::cout << PHWHERE << " ERROR: Can't find " << laserClusterNodeName << "." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_tGeometry = findNode::getClass<ActsGeometry>(topNode,
                                                 "ActsGeometry");
  if (!m_tGeometry)
  {
    std::cout << PHWHERE
              << "ActsGeometry not found on node tree. Exiting"
              << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  TrkrHitSetContainer::ConstRange hitsetrange = m_hits->getHitSets(TrkrDefs::TrkrId::tpcId);;
  
  struct thread_pair_t
  {
    pthread_t thread{};
    thread_data data;
  };

  std::vector<thread_pair_t> threads;
  threads.reserve(72);

  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  if (pthread_mutex_init(&mythreadlock, nullptr) != 0)
  {
    std::cout << std::endl << " mutex init failed" << std::endl;
    return 1;
  }
  
  for (unsigned int sec=0; sec<12; sec++)
  {
    for (int s=0; s<2; s++)
    {
      for (unsigned int mod=0; mod<3; mod++)
      {
	
        if(Verbosity() > 2)
        {
          std::cout << "making thread for side: " << s << "   sector: " << sec << "   module: " << mod << std::endl;
        }

	      thread_pair_t &thread_pair = threads.emplace_back();

	      std::vector<TrkrHitSet *> hitsets;
	      std::vector<unsigned int> layers;

	      std::vector<LaserCluster *> cluster_vector;
	      std::vector<TrkrDefs::cluskey> cluster_key_vector;

	      for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
	          hitsetitr != hitsetrange.second;
	          ++hitsetitr)
	      {
	        unsigned int layer = TrkrDefs::getLayer(hitsetitr->first);
	        int side = TpcDefs::getSide(hitsetitr->first);
	        unsigned int sector = TpcDefs::getSectorId(hitsetitr->first);
	        if (sector != sec || side != s)
	        {
	          continue;
	        }
	        if ((mod==0 && (layer<7 || layer>22)) || (mod==1 && (layer<=22 || layer>38) ) || (mod==2 && (layer<=38 || layer>54)))
	        {
	          continue;
	        }

	        TrkrHitSet *hitset = hitsetitr->second;
	  
	        hitsets.push_back(hitset);
	        layers.push_back(layer);
	  
      	}

	      thread_pair.data.geom_container = m_geom_container;
	      thread_pair.data.tGeometry = m_tGeometry;
	      thread_pair.data.hitsets = hitsets;
	      thread_pair.data.layers = layers;
	      thread_pair.data.side = (bool)s;
	      thread_pair.data.sector = sec;
        thread_pair.data.module = mod;
	      thread_pair.data.cluster_vector = cluster_vector;
	      thread_pair.data.cluster_key_vector = cluster_key_vector;
	      thread_pair.data.adc_threshold = m_adc_threshold;
	      thread_pair.data.peakTimeBin = m_laserEventInfo->getPeakSample(s);
	      thread_pair.data.layerMin = 3;
	      thread_pair.data.layerMax = 3;
	      thread_pair.data.tdriftmax = m_tdriftmax;
        thread_pair.data.eventNum = m_event;
        thread_pair.data.Verbosity = Verbosity();
        thread_pair.data.hitHist = nullptr;
        thread_pair.data.doFitting = m_do_fitting;

	      int rc;
	      rc = pthread_create(&thread_pair.thread, &attr, ProcessModule, (void *) &thread_pair.data);

	      if (rc)
	      {
	        std::cout << "Error:unable to create thread," << rc << std::endl;
	      }

	      if (m_do_sequential)
	      {
	        //wait for termination of thread
	        int rc2 = pthread_join(thread_pair.thread, nullptr);
	        if (rc2)
	        {
	          std::cout << "Error:unable to join," << rc2 << std::endl;
	        }

      	  //add clusters from thread to laserClusterContainer
	        const auto &data(thread_pair.data);
	        for(int index = 0; index < (int) data.cluster_vector.size(); ++index)
	        {
	          auto cluster = data.cluster_vector[index];
	          const auto ckey = data.cluster_key_vector[index];
	    
      	    m_clusterlist->addClusterSpecifyKey(ckey, cluster);
	        }
	      }
      }
    }
  }
  
  pthread_attr_destroy(&attr);

  if (!m_do_sequential)
  {
    for (const auto & thread_pair : threads)
    {
      int rc2 = pthread_join(thread_pair.thread, nullptr);
      if (rc2)
      {
	      std::cout << "Error:unable to join," << rc2 << std::endl;
      }
      
      //const auto &data(thread_pair.data);
      
      for(int index = 0; index < (int) thread_pair.data.cluster_vector.size(); ++index)
      {
	      auto cluster = thread_pair.data.cluster_vector[index];
	      const auto ckey = thread_pair.data.cluster_key_vector[index];
	
	      m_clusterlist->addClusterSpecifyKey(ckey, cluster);
      }
    }
  }

  threads.clear();

  if (Verbosity() > 1)
  {
    std::cout << "LaserClusterizer::process_event " << m_clusterlist->size() << " clusters found" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;

}
