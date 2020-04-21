#include "TpcClusterizer.h"

#include "TpcDefs.h"

#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssoc.h>
#include <trackbase/TrkrClusterv1.h>
#include <trackbase/TrkrDefs.h>  // for hitkey, getLayer
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>                         // for SubsysReco

#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>                         // for PHIODataNode
#include <phool/PHNode.h>                               // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>                             // for PHObject
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <TMatrixFfwd.h>    // for TMatrixF
#include <TMatrixT.h>       // for TMatrixT, ope...
#include <TMatrixTUtils.h>  // for TMatrixTRow

#include <TNtuple.h>  
#include <TFile.h>  

#include <cmath>  // for sqrt, cos, sin
#include <iostream>
#include <map>  // for _Rb_tree_cons...
#include <string>
#include <utility>  // for pair
#include <vector>

using namespace std;

TpcClusterizer::TpcClusterizer(const string &name)
  : SubsysReco(name)
  , m_hits(nullptr)
  , m_clusterlist(nullptr)
  , m_clusterhitassoc(nullptr)
  , zz_shaping_correction(0.0754)
  , pedestal(74.4)
  , NPhiBinsMax(0)
  , NPhiBinsMin(0)
  , NZBinsMax(0)
  , NZBinsMin(0)
  , hit_nt(nullptr)
  , cluster_nt(nullptr)
{
}

//===================
bool TpcClusterizer::is_local_maximum(int phibin, int zbin, std::vector<std::vector<double>> &adcval)
{
  bool retval = true;
  double centval = adcval[phibin][zbin];
  
  // search contiguous adc values for a larger signal
  for (int iz = zbin - 1; iz <= zbin + 1; iz++)
    for (int iphi = phibin - 1; iphi <= phibin + 1; iphi++)
      {
	if (iz >= NZBinsMax) continue;
	if (iz < NZBinsMin) continue;
	
	if (iphi >= NPhiBinsMax) continue;
	if (iphi < NPhiBinsMin) continue;
	
	if (adcval[iphi][iz] > centval)
	  {
	    retval = false;
	  }
      }
  
  return retval;
}

void TpcClusterizer::get_cluster(int phibin, int zbin, int &phiup, int &phidown, int &zup, int &zdown, std::vector<std::vector<double>> &adcval)
{
  // search along phi at the peak in z

  for (int iphi = phibin + 1; iphi < phibin + 4; iphi++)
  {
    if (iphi >= NPhiBinsMax) continue;

    if (adcval[iphi][zbin] <= 0)
      break;

    if (adcval[iphi][zbin] <= adcval[iphi - 1][zbin])
      phiup++;
    else
      break;
  }

  for (int iphi = phibin - 1; iphi > phibin - 4; iphi--)
  {
    if (iphi < NPhiBinsMin) continue;

    if (adcval[iphi][zbin] <= 0)
      break;

    if (adcval[iphi][zbin] <= adcval[iphi + 1][zbin])
      phidown++;
    else
      break;
  }

  // search along z at the peak in phi

  for (int iz = zbin + 1; iz < zbin + 10; iz++)
  {
    if (iz >= NZBinsMax) continue;

    if (adcval[phibin][iz] <= 0)
      break;

    if (adcval[phibin][iz] <= adcval[phibin][iz - 1])
      zup++;
    else
      break;
  }

  for (int iz = zbin - 1; iz > zbin - 10; iz--)
  {
    if (iz < NZBinsMin) continue;

    if (adcval[phibin][iz] <= 0)
      break;

    if (adcval[phibin][iz] <= adcval[phibin][iz + 1])
      zdown++;
    else
      break;
  }
}

int TpcClusterizer::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    cout << PHWHERE << "DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // Create the Cluster node if required
  TrkrClusterContainer *trkrclusters = findNode::getClass<TrkrClusterContainer>(dstNode, "TRKR_CLUSTER");
  if (!trkrclusters)
  {
    PHNodeIterator dstiter(dstNode);
    PHCompositeNode *DetNode =
        dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

    trkrclusters = new TrkrClusterContainer();
    PHIODataNode<PHObject> *TrkrClusterContainerNode =
        new PHIODataNode<PHObject>(trkrclusters, "TRKR_CLUSTER", "PHObject");
    DetNode->addNode(TrkrClusterContainerNode);
  }

  TrkrClusterHitAssoc *clusterhitassoc = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  if (!clusterhitassoc)
  {
    PHNodeIterator dstiter(dstNode);
    PHCompositeNode *DetNode =
        dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

    clusterhitassoc = new TrkrClusterHitAssoc();
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(clusterhitassoc, "TRKR_CLUSTERHITASSOC", "PHObject");
    DetNode->addNode(newNode);
  }

  if(Verbosity() > 10)
    {
      hit_nt = new TNtuple("hit_nt", "TPC hits", "phibin:zbin:layer:adc");
      cluster_nt = new TNtuple("cluster_nt", "TPC hits", "phibin:zbin:layer:adc");
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcClusterizer::process_event(PHCompositeNode *topNode)
{
  int print_layer = 18;

  if (Verbosity() > 1000)
    std::cout << "TpcClusterizer::Process_Event" << std::endl;

  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
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

  PHG4CylinderCellGeomContainer *geom_container =
      findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!geom_container)
  {
    std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // The hits are stored in hitsets, where each hitset contains all hits in a given TPC readout (layer, sector, side), so clusters are confined to a hitset
  // The TPC clustering is more complicated than for the silicon, because we have to deal with overlapping clusters

  // loop over the TPC HitSet objects
  TrkrHitSetContainer::ConstRange hitsetrange = m_hits->getHitSets(TrkrDefs::TrkrId::tpcId);
  for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
       hitsetitr != hitsetrange.second;
       ++hitsetitr)
  {
    TrkrHitSet *hitset = hitsetitr->second;
    int layer = TrkrDefs::getLayer(hitsetitr->first);
    if (Verbosity() > 1)
      if (layer == print_layer)
      {
	cout << "TpcClusterizer process hitsetkey " << hitsetitr->first
	     << " layer " << (int) TrkrDefs::getLayer(hitsetitr->first)
	     << " side " << (int) TpcDefs::getSide(hitsetitr->first)
	     << " sector " << (int) TpcDefs::getSectorId(hitsetitr->first)
	     << endl;
	if (Verbosity() > 5) hitset->identify();
      }

    // we have a single hitset, get the info that identifies the module
    // int sector = TpcDefs::getSectorId(hitsetitr->first);
    int side = TpcDefs::getSide(hitsetitr->first);

    // we will need the geometry object for this layer to get the global position
    PHG4CylinderCellGeom *layergeom = geom_container->GetLayerCellGeom(layer);
    int NPhiBins = layergeom->get_phibins();
    NPhiBinsMin = 0;
    NPhiBinsMax = NPhiBins;

    int NZBins = layergeom->get_zbins();
    if (side == 0)
    {
      NZBinsMin = 0;
      NZBinsMax = NZBins / 2;
    }
    else
    {
      NZBinsMin = NZBins / 2 + 1;
      NZBinsMax = NZBins;
    }

    // for convenience, create a 2D vector to store adc values in and initialize to zero
    std::vector<std::vector<double>> adcval(NPhiBins, std::vector<double>(NZBins, 0));

    TrkrHitSet::ConstRange hitrangei = hitset->getHits();
    for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
         hitr != hitrangei.second;
         ++hitr)
    {
      int phibin = TpcDefs::getPad(hitr->first);
      int zbin = TpcDefs::getTBin(hitr->first);
      if (hitr->second->getAdc() > 0)
	{
	  adcval[phibin][zbin] = (double) hitr->second->getAdc() - pedestal;

	  if (Verbosity() > 2)
	    if (layer == print_layer)
	      cout << " add hit in layer " << layer << " with phibin " << phibin << " zbin " << zbin << " adcval " << adcval[phibin][zbin] << endl;
	  
	  if(Verbosity() > 10)
	    //if (layer == print_layer)
	    hit_nt->Fill(phibin, zbin, layer, adcval[phibin][zbin]);
	}
    }

    int phibin_last = -1;
    int zbin_last = -1;
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
      if (!is_local_maximum(phibin, zbin, adcval)) continue;

      // eliminate double counting when two contiguous adcvals are identical (both bins will register as local maximum)
      if(phibin == phibin_last -1 || phibin == phibin_last || phibin == phibin_last + 1)
	if(zbin == zbin_last -1 || zbin == zbin_last || zbin == zbin_last + 1)
	  continue;
      phibin_last = phibin;
      zbin_last = zbin;

      int phiup = 0;
      int phidown = 0;
      int zup = 0;
      int zdown = 0;
      // cluster the hits around this local maximum
      get_cluster(phibin, zbin, phiup, phidown, zup, zdown, adcval);

      if (phiup == 0 && phidown == 0 && zup == 0 and zdown == 0) continue;  // ignore isolated noise hit

      // Add this cluster to a vector of clusters for later analysis
      phibinlo.push_back(phibin - phidown);
      phibinhi.push_back(phibin + phiup);
      zbinlo.push_back(zbin - zdown);
      zbinhi.push_back(zbin + zup);

      if(Verbosity() > 10) 
	//if (layer == print_layer)
	cluster_nt->Fill(phibin, zbin, layer, adcval[phibin][zbin]);

      if (Verbosity() > 2)
        if (layer == print_layer)
          cout << " cluster found in layer " << layer << " around hitkey " << hitr->first << " with zbin " << zbin << " zup " << zup << " zdown " << zdown
               << " phibin " << phibin << " phiup " << phiup << " phidown " << phidown << endl;
    }  // end loop over hits in this hitset

    // Now we analyze the clusters to get their parameters
    for (unsigned int iclus = 0; iclus < phibinlo.size(); iclus++)
    {
      //cout << "TpcClusterizer: process cluster iclus = " << iclus <<  " in layer " << layer << endl;
      // loop over the hits in this cluster
      double zsum = 0.0;
      double phi_sum = 0.0;
      double adc_sum = 0.0;
      double radius = layergeom->get_radius();  // returns center of layer
      if (Verbosity() > 2)
        if (layer == print_layer)
        {
          cout << "iclus " << iclus << " layer " << layer << " radius " << radius << endl;
          cout << "    z bin range " << zbinlo[iclus] << " to " << zbinhi[iclus] << " phibin range " << phibinlo[iclus] << " to " << phibinhi[iclus] << endl;
        }

      std::vector<TrkrDefs::hitkey> hitkeyvec;
      for (int iphi = phibinlo[iclus]; iphi <= phibinhi[iclus]; iphi++)
      {
        for (int iz = zbinlo[iclus]; iz <= zbinhi[iclus]; iz++)
        {
          double phi_center = layergeom->get_phicenter(iphi);
          double z = layergeom->get_zcenter(iz);

          zsum += z * adcval[iphi][iz];
          phi_sum += phi_center * adcval[iphi][iz];
          adc_sum += adcval[iphi][iz];

          // capture the hitkeys for all non-zero adc values
          if (adcval[iphi][iz] != 0)
          {
            TrkrDefs::hitkey hitkey = TpcDefs::genHitKey(iphi, iz);
            hitkeyvec.push_back(hitkey);
          }
        }
      }

      if (adc_sum < 10) continue;  // skip obvious noise "clusters"

      // This is the global position
      double clusphi = phi_sum / adc_sum;
      double clusz = zsum / adc_sum;

      // create the cluster entry directly in the node tree
      TrkrDefs::cluskey ckey = TpcDefs::genClusKey(hitset->getHitSetKey(), iclus);
      TrkrClusterv1 *clus = static_cast<TrkrClusterv1 *>((m_clusterlist->findOrAddCluster(ckey))->second);

      double phi_size = (double) (phibinhi[iclus] - phibinlo[iclus] + 1) * radius * layergeom->get_phistep();
      double z_size = (double) (zbinhi[iclus] - zbinlo[iclus] + 1) * layergeom->get_zstep();

      // Estimate the errors
      double dphi2_adc = 0.0;
      double dphi_adc = 0.0;
      double dz2_adc = 0.0;
      double dz_adc = 0.0;
      for (int iz = zbinlo[iclus]; iz <= zbinhi[iclus]; iz++)
      {
        for (int iphi = phibinlo[iclus]; iphi <= phibinhi[iclus]; iphi++)
        {
          double dphi = layergeom->get_phicenter(iphi) - clusphi;
          dphi2_adc += dphi * dphi * adcval[iphi][iz];
          dphi_adc += dphi * adcval[iphi][iz];

          double dz = layergeom->get_zcenter(iz) - clusz;
          dz2_adc += dz * dz * adcval[iphi][iz];
          dz_adc += dz * adcval[iphi][iz];
        }
      }
      double phi_cov = (dphi2_adc / adc_sum - dphi_adc * dphi_adc / (adc_sum * adc_sum));
      double z_cov = dz2_adc / adc_sum - dz_adc * dz_adc / (adc_sum * adc_sum);

      //cout << " layer " << layer << " z_cov " << z_cov << " dz2_adc " << dz2_adc << " adc_sum " <<  adc_sum << " dz_adc " << dz_adc << endl;

      // phi_cov = (weighted mean of dphi^2) - (weighted mean of dphi)^2,  which is essentially the weighted mean of dphi^2. The error is then:
      // e_phi = sigma_dphi/sqrt(N) = sqrt( sigma_dphi^2 / N )  -- where N is the number of samples of the distribution with standard deviation sigma_dphi
      //    - N is the number of electrons that drift to the readout plane
      // We have to convert (sum of adc units for all bins in the cluster) to number of ionization electrons N
      // Conversion gain is 20 mV/fC - relates total charge collected on pad to PEAK voltage out of ADC. The GEM gain is assumed to be 2000
      // To get equivalent charge per Z bin, so that summing ADC input voltage over all Z bins returns total input charge, divide voltages by 2.4 for 80 ns SAMPA
      // Equivalent charge per Z bin is then  (ADU x 2200 mV / 1024) / 2.4 x (1/20) fC/mV x (1/1.6e-04) electrons/fC x (1/2000) = ADU x 0.14

      double phi_err = radius * sqrt(phi_cov / (adc_sum * 0.14));
      if (phi_err == 0.0)  // a single phi bin will cause this
        phi_err = radius * layergeom->get_phistep() / sqrt(12.0);

      double z_err = sqrt(z_cov / (adc_sum * 0.14));
      if (z_err == 0.0)
        z_err = layergeom->get_zstep() / sqrt(12.0);

      // This corrects the bias introduced by the asymmetric SAMPA chip shaping - assumes 80 ns shaping time
      if (clusz < 0)
        clusz -= zz_shaping_correction;
      else
        clusz += zz_shaping_correction;

      // Fill in the cluster details
      //================
      clus->setAdc(adc_sum);
      clus->setPosition(0, radius * cos(clusphi));
      clus->setPosition(1, radius * sin(clusphi));
      clus->setPosition(2, clusz);
      clus->setGlobal();

      TMatrixF DIM(3, 3);
      DIM[0][0] = 0.0;
      DIM[0][1] = 0.0;
      DIM[0][2] = 0.0;
      DIM[1][0] = 0.0;
      //DIM[1][1] = pow(0.5 * phi_size / radius,2);  //cluster_v1 expects polar coordinates covariance
      DIM[1][1] = pow(0.5 * phi_size,2);  //cluster_v1 expects polar coordinates covariance
      //DIM[1][1] = pow(phi_size / radius,2);  //cluster_v1 expects polar coordinates covariance
      DIM[1][2] = 0.0;
      DIM[2][0] = 0.0;
      DIM[2][1] = 0.0;
      DIM[2][2] = pow(0.5 * z_size,2);
      //DIM[2][2] = pow(z_size,2);


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
      ROT[0][0] = cos(clusphi);
      ROT[0][1] = -sin(clusphi);
      ROT[0][2] = 0.0;
      ROT[1][0] = sin(clusphi);
      ROT[1][1] = cos(clusphi);
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

      /*
	  for(int i=0;i<3;i++)
	    for(int j=0;j<3;j++)
	      cout << "    i " << i << " j " << j << " clusphi " << clusphi << " clus error " << clus->getError(i,j) 
		   << " clus size " << clus->getSize(i,j) << " ROT " << ROT[i][j] << " ROT_T " << ROT_T[i][j] << " COVAR_ERR " << COVAR_ERR[i][j] << endl;
	  */

      // Add the hit associations to the TrkrClusterHitAssoc node
      // we need the cluster key and all associated hit keys (note: the cluster key includes the hitset key)
      for (unsigned int i = 0; i < hitkeyvec.size(); i++)
      {
        m_clusterhitassoc->addAssoc(ckey, hitkeyvec[i]);
      }

    }  // end loop over clusters for this hitset
  }    // end loop over hitsets

  if (Verbosity() > 100)
  {
    cout << "Dump clusters after TpcClusterizer" << endl;
    m_clusterlist->identify();
  }

  if (Verbosity() > 100)
  {
    cout << "Dump cluster hit associations after TpcClusterizer" << endl;
    m_clusterhitassoc->identify();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcClusterizer::End(PHCompositeNode *topNode)
{
  if (Verbosity() > 10)
  {
    TFile *outf = new TFile("cluster_nt_out.root", "recreate");
    outf->WriteTObject(hit_nt);
    outf->WriteTObject(cluster_nt);
    outf->Close();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
