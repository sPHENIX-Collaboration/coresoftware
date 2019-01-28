#include "TpcClusterizer.h"

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterv1.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include<TpcDefs.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderGeom.h>
#include <g4detectors/PHG4CylinderGeomContainer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHTypedNodeIterator.h>
#include <phool/getClass.h>

#include <TMath.h>

#include <TCanvas.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TProfile2D.h>
#include <TSpectrum2.h>
#include <TStopwatch.h>

#include <TMatrixF.h>

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <vector>

using namespace std;

TpcClusterizer::TpcClusterizer(const char* name)
  : zz_shaping_correction(0.0754)
  , pedestal(74.4)
{
}
//===================
TpcClusterizer::~TpcClusterizer()
{
}

bool TpcClusterizer::is_local_maximum(int phibin, int zbin, std::vector<std::vector<double>> &adcval)
{
  bool retval = true;
  double centval = adcval[phibin][zbin];
  //cout << "enter is_local_maximum for phibin " << phibin << " zbin " << zbin << " adcval " << centval <<  endl; 
  
  // search contiguous adc values for a larger signal 
  for(int iz = zbin - 4; iz <= zbin+4;iz++)
    for(int iphi = phibin -2; iphi <= phibin + 2; iphi++) 
      {
	//cout << " is_local_maximum: iphi " <<  iphi << " iz " << iz << " adcval " << adcval[iphi][iz] << endl;

	if(adcval[iphi][iz] > centval)
	  {
	    retval = false;
	  }
      }
  
	//if(retval)  cout << "**********************   success: returning " << retval << endl;
  
  return retval;
}
	
void TpcClusterizer::get_cluster(int phibin, int zbin, int &phiup, int &phidown, int &zup, int &zdown,  std::vector<std::vector<double>> &adcval)
{
  // search along phi at the peak in z

  for(int iphi = phibin+1; iphi<phibin+4; iphi++)
    {
      if(adcval[iphi][zbin] == 0)
	break;

      if( adcval[iphi][zbin] <= adcval[iphi-1][zbin] )
	phiup++;
      else
	break;
    }

  for(int iphi = phibin-1; iphi>phibin-4; iphi--)
    {
      if(adcval[iphi][zbin] == 0)
	break;

      if(adcval[iphi][zbin] <= adcval[iphi+1][zbin])
	phidown++;
      else
	break;
    }

  // search along z at the peak in phi

  for(int iz = zbin+1; iz<zbin+6; iz++)
    {
      if(adcval[phibin][iz] == 0)
	break;

      if(adcval[phibin][iz] <= adcval[phibin][iz-1])
	zup++;
      else
	break;
    }

  for(int iz =zbin-1; iz>zbin-6; iz--)
    {
      if(adcval[phibin][iz] == 0)
	break;

      if(adcval[phibin][iz] <= adcval[phibin][iz+1])
	zdown++;
      else
	break;
    }
}

int TpcClusterizer::InitRun(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcClusterizer::process_event(PHCompositeNode* topNode)
{
  //if (Verbosity() > 1000) 
    std::cout << "TpcClusterizer::Process_Event" << std::endl;

  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    cout << PHWHERE << "DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // get node containing the digitized hits
  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!m_hits)
  {
    cout << PHWHERE << "ERROR: Can't find node TRKR_HITSET" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // get node for clusters
  m_clusterlist = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!m_clusterlist)
  {
    cout << PHWHERE << " ERROR: Can't find TRKR_CLUSTER." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // get node for cluster hit associations
  m_clusterhitassoc = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  if (!m_clusterhitassoc)
  {
    cout << PHWHERE << " ERROR: Can't find TRKR_CLUSTERHITASSOC" << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PHG4CylinderCellGeomContainer* geom_container =
      findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!geom_container)
  {
    std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // The hits are stored in hitsets, where each hitset contains all hits in a given TPC readout (layer, sector, side), so clusters are confined to a hitset
  // The TPC clustering is more complicated than for the silicon, because we have to deal with overlapping clusters

  // loop over the TPC HitSet objects 
  TrkrHitSetContainer::ConstRange hitsetrange =  m_hits->getHitSets(TrkrDefs::TrkrId::tpcId);
  for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
       hitsetitr != hitsetrange.second;
       ++hitsetitr)
    {
      TrkrHitSet *hitset = hitsetitr->second;
      //if(Verbosity() > 1) 
	cout << "TpcClusterizer process hitsetkey " << hitsetitr->first << endl;
      if (Verbosity() > 2) hitset->identify();
      
      // we have a single hitset, get the info that identifies the module
      int layer = TrkrDefs::getLayer(hitsetitr->first);
      // int sector = TpcDefs::getSectorId(hitsetitr->first);
      int side = TpcDefs::getSide(hitsetitr->first);
      
      // we will need the geometry object for this layer to get the global position	
      PHG4CylinderCellGeom* layergeom = geom_container->GetLayerCellGeom(layer);
      int NPhiBins = layergeom->get_phibins();
      int NZBins = layergeom->get_zbins();
      if(side == 0)
	{
	  NZBinsMin = 0;
	  NZBinsMax = NZBins/2;
	}
      else
	{
	  NZBinsMin = NZBins/2 + 1;
	  NZBinsMax = NZBins;
	}

      // for convenience, create a 2D vector to store adc values in and initialize to zero
      std::vector<std::vector<double>> adcval (NPhiBins, std::vector<double>(NZBins, 0) );
       
      TrkrHitSet::ConstRange hitrangei = hitset->getHits();
      for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
	   hitr != hitrangei.second;
	   ++hitr)
	{
	  int phibin = TpcDefs::getPad(hitr->first);
	  int zbin = TpcDefs::getTBin(hitr->first);
	  if(hitr->second->getAdc() > 0)
	    adcval[phibin][zbin] = (double) hitr->second->getAdc() - pedestal;
	  //cout << " add hit with phibin " << phibin << " zbin " << zbin << " adcval " << adcval[phibin][zbin] << endl;
	}

      vector<int> phibinlo;
      vector<int> phibinhi;
      vector<int> zbinlo;
      vector<int> zbinhi;
      // we want to search the hit list for local maxima in phi-z space and cluster around them      
      for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
	   hitr != hitrangei.second;
	   ++hitr)
	{
	  int phibin = TpcDefs::getPad(hitr->first);
	  int zbin = TpcDefs::getTBin(hitr->first);
	  if(!is_local_maximum(phibin, zbin, adcval)) continue;
	  
	  int phiup= 0;
	  int phidown = 0;
	  int zup= 0;
	  int zdown = 0;
	  // cluster the hits around this local maximum
	  get_cluster(phibin, zbin, phiup, phidown, zup, zdown, adcval);

	  // Add this cluster to a vector of clusters for later analysis
	  phibinlo.push_back(phibin - phidown);
	  phibinhi.push_back(phibin + phiup);
	  zbinlo.push_back(zbin - zdown);
	  zbinhi.push_back(zbin + zup);
	  //cout << " cluster found with zbin " << zbin << " zup " << zup << " zdown " << zdown 
	  //   << " phibin " << phibin << " phiup " << phiup << " phidown " << phidown << endl; 
	}  // end loop over hits in this hiset

      // Now we analyze the clusters to get their parameters
      for(unsigned int iclus = 0; iclus < phibinlo.size(); iclus++)
	{  	  
	  // create the cluster entry directly in the node tree
	  TrkrDefs::cluskey ckey = TpcDefs::genClusKey(hitset->getHitSetKey(), iclus);
	  TrkrClusterv1 *clus = static_cast<TrkrClusterv1 *>((m_clusterlist->findOrAddCluster(ckey))->second);

	  // loop over the hits in this cluster
	  double xsum = 0.0;
	  double ysum = 0.0;
	  double zsum = 0.0;
	  int adc_sum = 0;
	  for(int iphi = phibinlo[iclus]; iphi < phibinhi[iclus]; iphi++)
	    {
	      for(int iz = zbinlo[iclus]; iz < zbinhi[iclus]; iz++)
		{
		  // find the center of the pixel in world coords
		  double phi_center = layergeom->get_phicenter(iphi);
		  double radius = layergeom->get_radius();
		  double x = radius * cos(phi_center);
		  double y = radius * sin(phi_center);
		  double z = layergeom->get_zcenter(iz);
	  
		  xsum += (double) adcval[iphi][iz] * x;
		  ysum += (double) adcval[iphi][iz] * y;
		  zsum += (double) adcval[iphi][iz] * z;
		  adc_sum += adcval[iphi][iz];
		}
	    }
	  
	  // This is the global position
	  double clusx = xsum / (double) adc_sum;
	  double clusy = ysum / (double) adc_sum;
	  double clusz = zsum / (double) adc_sum;
	  
	  // Fill in the cluster details
	  //================
	  clus->setAdc(adc_sum);
	  clus->setPosition(0, clusx);
	  clus->setPosition(1, clusy);
	  clus->setPosition(2, clusz);
	  clus->setGlobal();

	  double radius = layergeom->get_radius();
	  double phi_size = (double) (phibinlo[iclus] - phibinhi[iclus]) * layergeom->get_phistep();
	  double z_size = (double) (phibinlo[iclus] - phibinhi[iclus]) * layergeom->get_zstep();
	  double phi_err = layergeom->get_phistep() / sqrt(12.0);
	  double z_err = layergeom->get_zstep() / sqrt(12.0);
	  double phi = atan(clusy/clusx);

	  TMatrixF DIM(3, 3);
	  DIM[0][0] = 0.0;
	  DIM[0][1] = 0.0;
	  DIM[0][2] = 0.0;
	  DIM[1][0] = 0.0;
	  DIM[1][1] = phi_size / radius * phi_size / radius;  //cluster_v1 expects polar coordinates covariance
	  DIM[1][2] = 0.0;
	  DIM[2][0] = 0.0;
	  DIM[2][1] = 0.0;
	  DIM[2][2] = z_size * z_size;
	  
	  TMatrixF ERR(3, 3);
	  ERR[0][0] = 0.0;
	  ERR[0][1] = 0.0;
	  ERR[0][2] = 0.0;
	  ERR[1][0] = 0.0;
	  ERR[1][1] = phi_err * phi_err;  //cluster_v1 expects rad, arc, z as elementsof covariance
	  ERR[1][2] = 0.0;
	  ERR[2][0] = 0.0;
	  ERR[2][1] = 0.0;
	  ERR[2][2] = z_err * z_err;
	  
	  TMatrixF ROT(3, 3);
	  ROT[0][0] = cos(phi);
	  ROT[0][1] = -sin(phi);
	  ROT[0][2] = 0.0;
	  ROT[1][0] = sin(phi);
	  ROT[1][1] = cos(phi);
	  ROT[1][2] = 0.0;
	  ROT[2][0] = 0.0;
	  ROT[2][1] = 0.0;
	  ROT[2][2] = 1.0;
	  
	  TMatrixF ROT_T(3, 3);
	  ROT_T.Transpose(ROT);
	  
	  TMatrixF COVAR_DIM(3, 3);
	  COVAR_DIM = ROT * DIM * ROT_T;
	  
	  clus->setSize(0, 0, COVAR_DIM[0][0]);
	  clus->setSize(0, 1, COVAR_DIM[0][1]);
	  clus->setSize(0, 2, COVAR_DIM[0][2]);
	  clus->setSize(1, 0, COVAR_DIM[1][0]);
	  clus->setSize(1, 1, COVAR_DIM[1][1]);
	  clus->setSize(1, 2, COVAR_DIM[1][2]);
	  clus->setSize(2, 0, COVAR_DIM[2][0]);
	  clus->setSize(2, 1, COVAR_DIM[2][1]);
	  clus->setSize(2, 2, COVAR_DIM[2][2]);
	  //cout << " covar_dim[2][2] = " <<  COVAR_DIM[2][2] << endl;
	  
	  TMatrixF COVAR_ERR(3, 3);
	  COVAR_ERR = ROT * ERR * ROT_T;
	  
	  clus->setError(0, 0, COVAR_ERR[0][0]);
	  clus->setError(0, 1, COVAR_ERR[0][1]);
	  clus->setError(0, 2, COVAR_ERR[0][2]);
	  clus->setError(1, 0, COVAR_ERR[1][0]);
	  clus->setError(1, 1, COVAR_ERR[1][1]);
	  clus->setError(1, 2, COVAR_ERR[1][2]);
	  clus->setError(2, 0, COVAR_ERR[2][0]);
	  clus->setError(2, 1, COVAR_ERR[2][1]);
	  clus->setError(2, 2, COVAR_ERR[2][2]);
	  
	} // end loop over clusters for this hitset
    }  // end loop over hitsets

  cout << "Dump clusters after TpcClusterizer" << endl;
  m_clusterlist->identify();
  
  return Fun4AllReturnCodes::EVENT_OK;
}
