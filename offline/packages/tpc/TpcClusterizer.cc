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

namespace 
{
  template<class T> inline constexpr T square( const T& x ) { return x*x; }
}

using namespace std;

typedef std::pair<int, int> iphiz;
typedef std::pair<double, iphiz> ihit;

TpcClusterizer::TpcClusterizer(const string &name)
  : SubsysReco(name)
  , m_hits(nullptr)
  , m_clusterlist(nullptr)
  , m_clusterhitassoc(nullptr)
  , zz_shaping_correction(0.0754)
  , pedestal(74.4)
  , SectorFiducialCut(0.5)
  , NSearch(2)
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
  for (int iz = zbin - NSearch; iz <= zbin + NSearch; iz++)
    for (int iphi = phibin - NSearch; iphi <= phibin + NSearch; iphi++)
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

bool TpcClusterizer::is_in_sector_boundary(int phibin, int sector, PHG4CylinderCellGeom *layergeom)
{
  bool reject_it = false;

  // sector boundaries occur every 1/12 of the full phi bin range  
  int PhiBins = layergeom->get_phibins();
  int PhiBinsSector = PhiBins/12;

  double radius = layergeom->get_radius();
  double PhiBinSize = 2.0* radius * M_PI / (double) PhiBins;

  // sector starts where?
  int sector_lo = sector * PhiBinsSector;
  int sector_hi = sector_lo + PhiBinsSector - 1;

  int sector_fiducial_bins = (int) (SectorFiducialCut / PhiBinSize);

  if(phibin < sector_lo + sector_fiducial_bins || phibin > sector_hi - sector_fiducial_bins)
    {
      reject_it = true;
      /*
      int layer = layergeom->get_layer();
      cout << " local maximum is in sector fiducial boundary: layer " << layer << " radius " << radius << " sector " << sector 
      << " PhiBins " << PhiBins << " sector_fiducial_bins " << sector_fiducial_bins
      << " PhiBinSize " << PhiBinSize << " phibin " << phibin << " sector_lo " << sector_lo << " sector_hi " << sector_hi << endl;  
      */
    }

  return reject_it;
}

void TpcClusterizer::remove_hit(double adc, int phibin, int zbin, std::multimap<double, ihit> &all_hit_map, std::vector<std::vector<double>> &adcval)
{
  typedef multimap<double, ihit>::iterator hit_iterator;
  pair<hit_iterator, hit_iterator> iterpair = all_hit_map.equal_range(adc);
  hit_iterator it = iterpair.first;
  for (; it != iterpair.second; ++it) {
    if (it->second.second.first == phibin && it->second.second.second == zbin) { 
      if(Verbosity()>10) cout << " found it, erasing it "<< endl;
      all_hit_map.erase(it);
      break;
    }
  }
  adcval[phibin][zbin] = 0;
}

void TpcClusterizer::remove_hits(std::vector<ihit> &ihit_list, std::multimap<double, ihit> &all_hit_map,std::vector<std::vector<double>> &adcval )
{
  for(auto iter = ihit_list.begin(); iter != ihit_list.end();++iter){
    double adc = iter->first; 
    int phibin = iter->second.first;
    int zbin   = iter->second.second;
    if(Verbosity()>10) cout << "remove adc: " << adc << " phi: " << phibin << " zbin " << zbin << endl; 
    remove_hit(adc,phibin,zbin,all_hit_map,adcval);

  }
}

void TpcClusterizer::calc_cluster_parameter(std::vector<ihit> &ihit_list,int iclus, PHG4CylinderCellGeom *layergeom, TrkrHitSet *hitset)
{

  ///HERE TODO
  //cout << "TpcClusterizer: process cluster iclus = " << iclus <<  " in layer " << layer << endl;
  // loop over the hits in this cluster
  double z_sum = 0.0;
  double phi_sum = 0.0;
  double adc_sum = 0.0;
  double z2_sum = 0.0;
  double phi2_sum = 0.0;

  double radius = layergeom->get_radius();  // returns center of layer
    
  int phibinhi = -1;
  int phibinlo = 666666;
  int zbinhi = -1;
  int zbinlo = 666666;
  int clus_size = ihit_list.size();

  if(clus_size == 1) return;

  std::vector<TrkrDefs::hitkey> hitkeyvec;
  for(auto iter = ihit_list.begin(); iter != ihit_list.end();++iter){
    double adc = iter->first; 

    if (adc <= 0) continue;

    int iphi = iter->second.first;
    int iz   = iter->second.second;
    if(iphi > phibinhi) phibinhi = iphi;
    if(iphi < phibinlo) phibinlo = iphi;
    if(iz > zbinhi) zbinhi = iz;
    if(iz < zbinlo) zbinlo = iz;

    //    if(Verbosity()>10) cout << "remove adc: " << adc << " phi: " << iphi << " zbin " << iz << endl; 

    // update phi sums
    double phi_center = layergeom->get_phicenter(iphi);
    phi_sum += phi_center * adc;
    phi2_sum += square(phi_center)*adc;

    // update z sums
    double z = layergeom->get_zcenter(iz);	  
    z_sum += z * adc;
    z2_sum += square(z)*adc;

    adc_sum += adc;

    // capture the hitkeys for all non-zero adc values
    TrkrDefs::hitkey hitkey = TpcDefs::genHitKey(iphi, iz);
    hitkeyvec.push_back(hitkey);
  }
  if (adc_sum < 10) return;  // skip obvious noise "clusters"
  
  // This is the global position
  double clusphi = phi_sum / adc_sum;
  double clusz = z_sum / adc_sum;
  
  const double phi_cov = phi2_sum/adc_sum - square(clusphi);
  const double z_cov = z2_sum/adc_sum - square(clusz);
  
  // create the cluster entry directly in the node tree
  TrkrDefs::cluskey ckey = TpcDefs::genClusKey(hitset->getHitSetKey(), iclus);
  TrkrClusterv1 *clus = static_cast<TrkrClusterv1 *>((m_clusterlist->findOrAddCluster(ckey))->second);
  
  //  int phi_nsize = phibinhi - phibinlo + 1;
  //  int z_nsize   = zbinhi   - zbinlo + 1;
  
  double phi_size = (double) (phibinhi - phibinlo + 1) * radius * layergeom->get_phistep();
  double z_size = (double) (zbinhi - zbinlo + 1) * layergeom->get_zstep();

  // Estimate the errors
  const double phi_err_square = (phibinhi == phibinlo) ?
    square(radius*layergeom->get_phistep())/12:
    square(radius)*phi_cov/(adc_sum*0.14);
  
  const double z_err_square = (zbinhi == zbinlo) ?
    square(layergeom->get_zstep())/12:
    z_cov/(adc_sum*0.14);

  //cout << "   layer " << layer << " z_cov " << z_cov << " dz2_adc " << dz2_adc << " adc_sum " <<  adc_sum << " dz_adc " << dz_adc << endl;
  
  // phi_cov = (weighted mean of dphi^2) - (weighted mean of dphi)^2,  which is essentially the weighted mean of dphi^2. The error is then:
  // e_phi = sigma_dphi/sqrt(N) = sqrt( sigma_dphi^2 / N )  -- where N is the number of samples of the distribution with standard deviation sigma_dphi
  //    - N is the number of electrons that drift to the readout plane
  // We have to convert (sum of adc units for all bins in the cluster) to number of ionization electrons N
  // Conversion gain is 20 mV/fC - relates total charge collected on pad to PEAK voltage out of ADC. The GEM gain is assumed to be 2000
  // To get equivalent charge per Z bin, so that summing ADC input voltage over all Z bins returns total input charge, divide voltages by 2.4 for 80 ns SAMPA
  // Equivalent charge per Z bin is then  (ADU x 2200 mV / 1024) / 2.4 x (1/20) fC/mV x (1/1.6e-04) electrons/fC x (1/2000) = ADU x 0.14

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
  DIM[1][1] = pow(0.5 * phi_size,2);  //cluster_v1 expects 1/2 of actual size
  DIM[1][2] = 0.0;
  DIM[2][0] = 0.0;
  DIM[2][1] = 0.0;
  DIM[2][2] = pow(0.5 * z_size,2);
  
  TMatrixF ERR(3, 3);
  ERR[0][0] = 0.0;
  ERR[0][1] = 0.0;
  ERR[0][2] = 0.0;
  ERR[1][0] = 0.0;
  ERR[1][1] = phi_err_square;  //cluster_v1 expects rad, arc, z as elementsof covariance
  ERR[1][2] = 0.0;
  ERR[2][0] = 0.0;
  ERR[2][1] = 0.0;
  ERR[2][2] = z_err_square;
  
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
  
  // Add the hit associations to the TrkrClusterHitAssoc node
  // we need the cluster key and all associated hit keys (note: the cluster key includes the hitset key)
  for (unsigned int i = 0; i < hitkeyvec.size(); i++)
    {
      m_clusterhitassoc->addAssoc(ckey, hitkeyvec[i]);
    }

}

void TpcClusterizer::print_cluster(std::vector<ihit> &ihit_list)
{
  for(auto iter = ihit_list.begin(); iter != ihit_list.end();++iter){
    double adc = iter->first; 
    int phibin = iter->second.first;
    int zbin   = iter->second.second;
    if(Verbosity()>10) cout << "  iz: " << zbin << " iphi: " << phibin << " adc: " << adc << endl;
  }
}

void TpcClusterizer::find_z_range(int phibin, int zbin, std::vector<std::vector<double>> &adcval, int& zdown, int& zup){

  int FitRangeZ = 5;
  zup = 0;
  zdown = 0;
  for(int iz=0; iz< FitRangeZ; iz++){
    int cz = zbin + iz;

    if(Verbosity()>10) cout << " cz " << cz << " adc: " << adcval[phibin][cz] << " phibin: " << phibin << " zbin " << cz << endl;
    if(cz <= 0 || cz >= NZBinsMax){
      // zup = iz;
      break; // truncate edge
    }
    
    //break when below minimum
    if(adcval[phibin][cz] <= 0) {
      //zup = iz;
      if(Verbosity() > 1000) cout << " failed threshold cut, set izup to " << zup << endl;
      break;
    }
    //check local minima and break at minimum.
    if(cz<NZBinsMax-4){//make sure we stay clear from the edge
      if(adcval[phibin][cz]+adcval[phibin][cz+1] < 
	 adcval[phibin][cz+2]+adcval[phibin][cz+3]){//rising again
	zup = iz+1;
	break;
      }
    }
    zup = iz;
  }
  for(int iz=0; iz< FitRangeZ; iz++){
    int cz = zbin - iz;
    if(Verbosity()>10) cout << " cz " << cz << " adc: " << adcval[phibin][cz] << " phibin: " << phibin << " zbin " << cz << endl;
    if(cz <= 0 || cz >= NZBinsMax){
      //      zdown = iz;
      break; // truncate edge
    }
    if(adcval[phibin][cz] <= 0) {
      // zdown = iz;
      if(Verbosity() > 1000) cout << " failed threshold cut, set izup to " << zup << endl;
      break;
    }

    if(cz>4){//make sure we stay clear from the edge
      if(adcval[phibin][cz]+adcval[phibin][cz-1] < 
	 adcval[phibin][cz-2]+adcval[phibin][cz-3]){//rising again
	zdown = iz+1;
	break;
      }
    }
    zdown = iz;
  }
}

void TpcClusterizer::find_phi_range(int phibin, int zbin, std::vector<std::vector<double>> &adcval, int& phidown, int& phiup){

  int FitRangePHI = 3;
  phidown = 0;
  phiup = 0;
  for(int iphi=0; iphi< FitRangePHI; iphi++){
    int cphi = phibin + iphi;
    if(Verbosity()>10) cout << " cphi " << cphi << " adc: " << adcval[cphi][zbin] << " phibin: " << phibin << " zbin " << zbin << endl;
    if(cphi < 0 || cphi >= NPhiBinsMax){
      // phiup = iphi;
      break; // truncate edge
    }
    // consider only the peak bin in phi when searching for Z limit     
    
    //break when below minimum
    if(adcval[cphi][zbin] <= 0) {
      // phiup = iphi;
      if(Verbosity() > 1000) cout << " failed threshold cut, set iphiup to " << phiup << endl;
      break;
    }
    //check local minima and break at minimum.
    if(cphi<NPhiBinsMax-4){//make sure we stay clear from the edge
      if(adcval[cphi][zbin]+adcval[cphi+1][zbin] < 
	 adcval[cphi+2][zbin]+adcval[cphi+3][zbin]){//rising again
	phiup = iphi+1;
	break;
      }
    }
    phiup = iphi;
  }
  if(Verbosity()>10) cout << " phiup " << phiup << endl; 
  for(int iphi=0; iphi< FitRangePHI; iphi++){
    int cphi = phibin - iphi;
    if(Verbosity()>10) cout << " cphi " << cphi << " adc: " << adcval[cphi][zbin] << " phibin: " << phibin << " zbin " << zbin << endl;
    if(cphi < 0 || cphi >= NPhiBinsMax){
      // phidown = iphi;
      break; // truncate edge
    }
    
    if(adcval[cphi][zbin] <= 0) {
      //phidown = iphi;
      if(Verbosity() > 1000) cout << " failed threshold cut, set iphiup to " << phiup << endl;
      break;
    }

    if(cphi>4){//make sure we stay clear from the edge
      if(adcval[cphi][zbin]+adcval[cphi-1][zbin] < 
	 adcval[cphi-2][zbin]+adcval[cphi-3][zbin]){//rising again
	phidown = iphi+1;
	break;
      }
    }
    phidown = iphi;
  }
}

void TpcClusterizer::get_cluster(int phibin, int zbin, std::vector<std::vector<double>> &adcval, std::vector<ihit> &ihit_list)
{
  // search along phi at the peak in z
 
  int zup =0;
  int zdown =0;
  find_z_range(phibin, zbin, adcval, zdown, zup);
  if(Verbosity()>10) cout << " zbin: " << zbin << " zdown: " << zdown << " zup: " << zup << " phi " << phibin << endl;
  //now we have the z extent of the cluster, go find the phi edges

  for(int iz=zbin - zdown ; iz<= zbin + zup; iz++){
    int phiup = 0;
    int phidown = 0;
    find_phi_range(phibin, iz, adcval, phidown, phiup);
    if(Verbosity()>10) cout << "phibin: " << phibin << " zbin: " << " phidown " << phidown << " phiup " << phiup  << endl;
    for (int iphi = phibin - phidown; iphi <= (phibin + phiup); iphi++){
      iphiz iCoord(make_pair(iphi,iz));
      ihit  thisHit(adcval[iphi][iz],iCoord);
      ihit_list.push_back(thisHit);
    }
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
    if (Verbosity() > 2)
      if (layer == print_layer)
      {
	cout << endl << "TpcClusterizer process hitsetkey " << hitsetitr->first
	     << " layer " << (int) TrkrDefs::getLayer(hitsetitr->first)
	     << " side " << (int) TpcDefs::getSide(hitsetitr->first)
	     << " sector " << (int) TpcDefs::getSectorId(hitsetitr->first)
	     << " nhits " << (int) hitset->size() 
	     << endl;
	//	if (Verbosity() > 5) hitset->identify();
      }

    // we have a single hitset, get the info that identifies the module
    //    int sector = TpcDefs::getSectorId(hitsetitr->first);
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

    std::multimap<double, ihit> all_hit_map;
    
    //put all hits in the all_hit_map (sorted by adc)
    //start with highest adc hit
    // -> cluster around it
    // -> vector of hits
    // -> calculate cluster parameters
    // -> add hits to truth association
    // remove hits from all_hit_map
    // repeat untill all_hit_map empty

    TrkrHitSet::ConstRange hitrangei = hitset->getHits();
    if(Verbosity()>10) cout << " New Hitset " << endl;
    for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
         hitr != hitrangei.second;
         ++hitr)
    {

      int phibin = TpcDefs::getPad(hitr->first);
      int zbin = TpcDefs::getTBin(hitr->first);
      double adc =  (double) hitr->second->getAdc() - pedestal;
      
      if(Verbosity()>10) cout << " iphi: " << phibin << " iz: " << zbin << " adc: "  << adc << endl;

      if (hitr->second->getAdc() > 0)
	{
	  if(adc>0){
	    iphiz iCoord(make_pair(phibin,zbin));
	    ihit  thisHit(adc,iCoord);
	    all_hit_map.insert(make_pair(adc, thisHit));
	  }
	  adcval[phibin][zbin] = (double) hitr->second->getAdc() - pedestal;

	  if (Verbosity() > 2)
	    if (layer == print_layer)
	      cout << " add hit in layer " << layer << " with phibin " << phibin << " zbin " << zbin << " adcval " << adcval[phibin][zbin] << endl;
	}
      
    }
    //put all hits in the all_hit_map (sorted by adc)

    int nclus = 0;
    while(all_hit_map.size()>0){
      auto iter = all_hit_map.rbegin();
      if(iter == all_hit_map.rend()) break;

      ihit hiHit = iter->second;
      //start with highest adc hit
      if(Verbosity()>10) cout << "  test entries: " << all_hit_map.size() << " adc: " << hiHit.first << " iphi: " << hiHit. second.first << " iz: " << hiHit.second.second << endl;
      //      double adc = hiHit.first;
      int iphi = hiHit.second.first;
      int iz = hiHit.second.second;

      //put all hits in the all_hit_map (sorted by adc)
      //start with highest adc hit
      // -> cluster around it and get vector of hits
      std::vector<ihit> ihit_list;
      get_cluster(iphi, iz, adcval, ihit_list);
      if(Verbosity()>10) 
	cout << " cluster size: " << ihit_list.size() << " #clusters: " << nclus<< endl;
      // -> calculate cluster parameters
      // -> add hits to truth association
      // remove hits from all_hit_map
      // repeat untill all_hit_map empty
      //      print_cluster(ihit_list);
      calc_cluster_parameter(ihit_list,nclus++, layergeom, hitset);
      remove_hits(ihit_list,all_hit_map,adcval);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcClusterizer::End(PHCompositeNode *topNode)
{
  /*
  if (Verbosity() > 10)
  {
    TFile *outf = new TFile("cluster_nt_out.root", "recreate");
    outf->WriteTObject(hit_nt);
    outf->WriteTObject(cluster_nt);
    outf->Close();
  }
  */
  return Fun4AllReturnCodes::EVENT_OK;
}
