#include "TpcClusterizer.h"

#include "TpcDefs.h"

#include <trackbase/TrkrClusterContainerv4.h>
#include <trackbase/TrkrClusterv3.h>
#include <trackbase/TrkrClusterv4.h>
#include <trackbase/TrkrClusterHitAssocv3.h>
#include <trackbase/TrkrDefs.h>  // for hitkey, getLayer
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>                         // for SubsysReco

#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>

#include <Acts/Definitions/Units.hpp>
#include <Acts/Surfaces/Surface.hpp>

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

#include <TFile.h>  

#include <cmath>  // for sqrt, cos, sin
#include <iostream>
#include <map>  // for _Rb_tree_cons...
#include <string>
#include <utility>  // for pair
#include <array>
#include <vector>
#include <limits>
// Terra incognita....
#include <pthread.h>

namespace 
{
  template<class T> inline constexpr T square( const T& x ) { return x*x; }
  
  using assoc = std::pair<TrkrDefs::cluskey, TrkrDefs::hitkey>;

  struct ihit
  {
    unsigned short iphi = 0;
    unsigned short iz = 0;
    unsigned short adc = 0;
    unsigned short edge = 0;
  };

  struct thread_data 
  {
    PHG4CylinderCellGeom *layergeom = nullptr;
    TrkrHitSet *hitset = nullptr;
    ActsGeometry *tGeometry = nullptr;
    unsigned int layer = 0;
    int side = 0;
    unsigned int sector = 0;
    float pedestal = 0;
    bool do_assoc = true;
    bool do_wedge_emulation = true;
    unsigned short phibins = 0;
    unsigned short phioffset = 0;
    unsigned short zbins = 0;
    unsigned short zoffset = 0;
    unsigned short maxHalfSizeZ = 0;
    unsigned short maxHalfSizePhi = 0;
    double m_drift_velocity_scale = 1.0;
    double par0_neg = 0;
    double par0_pos = 0;
    int cluster_version = 3;
    std::vector<assoc> association_vector;
    std::vector<TrkrCluster*> cluster_vector;
    int verbosity = 0;
  };
  
  pthread_mutex_t mythreadlock;
  
  void remove_hit(double adc, int phibin, int zbin, int edge, std::multimap<unsigned short, ihit> &all_hit_map, std::vector<std::vector<unsigned short>> &adcval)
  {
    typedef std::multimap<unsigned short, ihit>::iterator hit_iterator;
    std::pair<hit_iterator, hit_iterator> iterpair = all_hit_map.equal_range(adc);
    hit_iterator it = iterpair.first;
    for (; it != iterpair.second; ++it) {
      if (it->second.iphi == phibin && it->second.iz == zbin) { 
	all_hit_map.erase(it);
	break;
      }
    }
    if(edge)
      adcval[phibin][zbin] = USHRT_MAX;
    else
      adcval[phibin][zbin] = 0;
  }
  
  void remove_hits(std::vector<ihit> &ihit_list, std::multimap<unsigned short, ihit> &all_hit_map,std::vector<std::vector<unsigned short>> &adcval)
  {
    for(auto iter = ihit_list.begin(); iter != ihit_list.end();++iter){
      unsigned short adc    = iter->adc; 
      unsigned short phibin = iter->iphi;
      unsigned short zbin   = iter->iz;
      unsigned short edge   = iter->edge;
      remove_hit(adc,phibin,zbin,edge,all_hit_map,adcval);
    }
  }
  
  void find_z_range(int phibin, int zbin, const thread_data& my_data, const std::vector<std::vector<unsigned short>> &adcval, int& zdown, int& zup, int &touch, int &edge){
	
    const int FitRangeZ = (int) my_data.maxHalfSizeZ;
    const int NZBinsMax = (int) my_data.zbins;
    zup = 0;
    zdown = 0;
    for(int iz=0; iz< FitRangeZ; iz++){
      int cz = zbin + iz;
      
      if(cz <= 0 || cz >= NZBinsMax){
	// zup = iz;
	edge++;
	break; // truncate edge
      }
      
      if(adcval[phibin][cz] <= 0) {
	break;
      }
      if(adcval[phibin][cz] == USHRT_MAX) {
	touch++;
	break;
      }
      //check local minima and break at minimum.
      if(cz<NZBinsMax-4){//make sure we stay clear from the edge
	if(adcval[phibin][cz]+adcval[phibin][cz+1] < 
	   adcval[phibin][cz+2]+adcval[phibin][cz+3]){//rising again
	  zup = iz+1;
	  touch++;
	  break;
	      }
      }
      zup = iz;
    }
    for(int iz=0; iz< FitRangeZ; iz++){
      int cz = zbin - iz;
      if(cz <= 0 || cz >= NZBinsMax){
	//      zdown = iz;
	edge++;
	break; // truncate edge
      }
      if(adcval[phibin][cz] <= 0) {
	break;
      }
      if(adcval[phibin][cz] == USHRT_MAX) {
	touch++;
	break;
      }
      if(cz>4){//make sure we stay clear from the edge
	if(adcval[phibin][cz]+adcval[phibin][cz-1] < 
	   adcval[phibin][cz-2]+adcval[phibin][cz-3]){//rising again
	  zdown = iz+1;
	  touch++;
	  break;
	}
      }
      zdown = iz;
    }
    return;
  }
	
  void find_phi_range(int phibin, int zbin, const thread_data& my_data, const std::vector<std::vector<unsigned short>> &adcval, int& phidown, int& phiup, int &touch, int &edge)
  {
	
    int FitRangePHI = (int) my_data.maxHalfSizePhi;
    int NPhiBinsMax = (int) my_data.phibins;
    phidown = 0;
    phiup = 0;
    for(int iphi=0; iphi< FitRangePHI; iphi++){
      int cphi = phibin + iphi;
      if(cphi < 0 || cphi >= NPhiBinsMax){
	// phiup = iphi;
	edge++;
	break; // truncate edge
      }
      
      //break when below minimum
      if(adcval[cphi][zbin] <= 0) {
	// phiup = iphi;
	break;
      }
      if(adcval[cphi][zbin] == USHRT_MAX) {
	touch++;
	break;
      }
      //check local minima and break at minimum.
      if(cphi<NPhiBinsMax-4){//make sure we stay clear from the edge
	if(adcval[cphi][zbin]+adcval[cphi+1][zbin] < 
	   adcval[cphi+2][zbin]+adcval[cphi+3][zbin]){//rising again
	  phiup = iphi+1;
	  touch++;
	  break;
	}
      }
      phiup = iphi;
    }
    
    for(int iphi=0; iphi< FitRangePHI; iphi++){
      int cphi = phibin - iphi;
      if(cphi < 0 || cphi >= NPhiBinsMax){
	// phidown = iphi;
	edge++;
	break; // truncate edge
      }
      
      if(adcval[cphi][zbin] <= 0) {
	//phidown = iphi;
	break;
      }
      if(adcval[cphi][zbin] == USHRT_MAX) {
	touch++;
	break;
      }
      if(cphi>4){//make sure we stay clear from the edge
	if(adcval[cphi][zbin]+adcval[cphi-1][zbin] < 
	   adcval[cphi-2][zbin]+adcval[cphi-3][zbin]){//rising again
	  phidown = iphi+1;
	  touch++;
	  break;
	}
      }
      phidown = iphi;
    }
    return;
  }
	
  void get_cluster(int phibin, int zbin, const thread_data& my_data, const std::vector<std::vector<unsigned short>> &adcval, std::vector<ihit> &ihit_list, int &touch, int &edge)
	{
	  // search along phi at the peak in z
	 
	  int zup =0;
	  int zdown =0;
	  find_z_range(phibin, zbin, my_data, adcval, zdown, zup, touch, edge);
	  //now we have the z extent of the cluster, go find the phi edges
	
	  for(int iz=zbin - zdown ; iz<= zbin + zup; iz++){
	    int phiup = 0;
	    int phidown = 0;
	    find_phi_range(phibin, iz, my_data, adcval, phidown, phiup, touch, edge);
	    for (int iphi = phibin - phidown; iphi <= (phibin + phiup); iphi++){
	      if(adcval[iphi][iz]>0 && adcval[iphi][iz]!=USHRT_MAX){
		ihit hit;
		hit.iphi = iphi;
		hit.iz = iz;
		hit.adc = adcval[iphi][iz];
		if(touch>0){
		  if((iphi == (phibin - phidown))||
		     (iphi == (phibin + phiup))){
		    hit.edge = 1;
		  }
		}
		ihit_list.push_back(hit);
	      }
	    }
	  }
	  return;
	}
  
    void calc_cluster_parameter(const std::vector<ihit> &ihit_list, thread_data& my_data, int ntouch, int nedge )
    {
    
      // get z range from layer geometry
      /* these are used for rescaling the drift velocity */
      const double z_min = -105.5;
      const double z_max = 105.5;
      
      // loop over the hits in this cluster
      double z_sum = 0.0;
      double phi_sum = 0.0;
      double adc_sum = 0.0;
      double z2_sum = 0.0;
      double phi2_sum = 0.0;
      
      double radius = my_data.layergeom->get_radius();  // returns center of layer
      
      int phibinhi = -1;
      int phibinlo = 666666;
      int zbinhi = -1;
      int zbinlo = 666666;
      int clus_size = ihit_list.size();
      
      if(clus_size == 1) return;
      
      std::vector<TrkrDefs::hitkey> hitkeyvec;
      for(auto iter = ihit_list.begin(); iter != ihit_list.end();++iter){
	double adc = iter->adc; 
	
	if (adc <= 0) continue;
	
	int iphi = iter->iphi + my_data.phioffset;
	int iz   = iter->iz + my_data.zoffset;
	if(iphi > phibinhi) phibinhi = iphi;
	if(iphi < phibinlo) phibinlo = iphi;
	if(iz > zbinhi) zbinhi = iz;
	if(iz < zbinlo) zbinlo = iz;
	
	// update phi sums
	double phi_center = my_data.layergeom->get_phicenter(iphi);
	phi_sum += phi_center * adc;
	phi2_sum += square(phi_center)*adc;
	
	// update z sums
	double z = my_data.layergeom->get_zcenter(iz);
	
	// apply drift velocity scale
	/* this formula ensures that z remains unchanged when located on one of the readout plane, at either z_max or z_min */
	z =
	  z*my_data.m_drift_velocity_scale +
	  (z>0 ? z_max:z_min)*(1.-my_data.m_drift_velocity_scale);
	
	z_sum += z*adc;
	z2_sum += square(z)*adc;
	
	adc_sum += adc;
	
	// capture the hitkeys for all adc values above a certain threshold
	TrkrDefs::hitkey hitkey = TpcDefs::genHitKey(iphi, iz);
	// if(adc>5)
	hitkeyvec.push_back(hitkey);
      }
      if (adc_sum < 10){
	hitkeyvec.clear();
	return;  // skip obvious noise "clusters"
      }  
      // This is the global position
      double clusphi = phi_sum / adc_sum;
      double clusz = z_sum / adc_sum;
      float clusx = radius * cos(clusphi);
      float clusy = radius * sin(clusphi);
      
      const double phi_cov = phi2_sum/adc_sum - square(clusphi);
      const double z_cov = z2_sum/adc_sum - square(clusz);
       // Get the surface key to find the surface from the 

      TrkrDefs::hitsetkey tpcHitSetKey = TpcDefs::genHitSetKey( my_data.layer, my_data.sector, my_data.side );      
      Acts::Vector3 global(clusx, clusy, clusz);
      TrkrDefs::subsurfkey subsurfkey = 0;
    
      Surface surface = my_data.tGeometry->get_tpc_surface_from_coords(
         tpcHitSetKey,
	 global,
	 subsurfkey);
    
      if(!surface)
	{
	  /// If the surface can't be found, we can't track with it. So 
	  /// just return and don't add the cluster to the container
	  hitkeyvec.clear();
	  return;
	}
      
      // Estimate the errors
      const double phi_err_square = (phibinhi == phibinlo) ?
	square(radius*my_data.layergeom->get_phistep())/12:
	square(radius)*phi_cov/(adc_sum*0.14);
      
      const double z_err_square = (zbinhi == zbinlo) ?
	square(my_data.layergeom->get_zstep())/12:
	z_cov/(adc_sum*0.14);
      
      char zsize = zbinhi - zbinlo + 1;
      char phisize = phibinhi - phibinlo + 1;
      // phi_cov = (weighted mean of dphi^2) - (weighted mean of dphi)^2,  which is essentially the weighted mean of dphi^2. The error is then:
      // e_phi = sigma_dphi/sqrt(N) = sqrt( sigma_dphi^2 / N )  -- where N is the number of samples of the distribution with standard deviation sigma_dphi
      //    - N is the number of electrons that drift to the readout plane
      // We have to convert (sum of adc units for all bins in the cluster) to number of ionization electrons N
      // Conversion gain is 20 mV/fC - relates total charge collected on pad to PEAK voltage out of ADC. The GEM gain is assumed to be 2000
      // To get equivalent charge per Z bin, so that summing ADC input voltage over all Z bins returns total input charge, divide voltages by 2.4 for 80 ns SAMPA
      // Equivalent charge per Z bin is then  (ADU x 2200 mV / 1024) / 2.4 x (1/20) fC/mV x (1/1.6e-04) electrons/fC x (1/2000) = ADU x 0.14
      clusz -= (clusz<0) ? my_data.par0_neg:my_data.par0_pos;
      
      Acts::Vector3 center = surface->center(my_data.tGeometry->geometry().geoContext)/Acts::UnitConstants::cm;
  
      // no conversion needed, only used in acts
      Acts::Vector3 normal = surface->normal(my_data.tGeometry->geometry().geoContext);
      double clusRadius = sqrt(clusx * clusx + clusy * clusy);
      double rClusPhi = clusRadius * clusphi;
      double surfRadius = sqrt(center(0)*center(0) + center(1)*center(1));
      double surfPhiCenter = atan2(center[1], center[0]);
      double surfRphiCenter = surfPhiCenter * surfRadius;
      double surfZCenter = center[2];
     
      auto local = surface->globalToLocal(my_data.tGeometry->geometry().geoContext,
					  global * Acts::UnitConstants::cm,
					  normal);
      Acts::Vector2 localPos;

      // Prefer Acts transformation since we build the TPC surfaces manually
      if(local.ok())
	{
	  localPos = local.value() / Acts::UnitConstants::cm;
	}
      else
	{
	  /// otherwise take the manual calculation
	  localPos(0) = rClusPhi - surfRphiCenter;
	  localPos(1) = clusz - surfZCenter; 
	}
      
      
      // we need the cluster key and all associated hit keys (note: the cluster key includes the hitset key)
      
      if(my_data.cluster_version==3){
	
	// Fill in the cluster details
	//================
	auto clus = new TrkrClusterv3;
	//auto clus = std::make_unique<TrkrClusterv3>();
	clus->setAdc(adc_sum);      
	clus->setSubSurfKey(subsurfkey);      
	clus->setLocalX(localPos(0));
	clus->setLocalY(localPos(1));
	clus->setActsLocalError(0,0, phi_err_square);
	clus->setActsLocalError(1,0, 0);
	clus->setActsLocalError(0,1, 0);
	clus->setActsLocalError(1,1, z_err_square);
	my_data.cluster_vector.push_back(clus);
      }else if(my_data.cluster_version==4){
	auto clus = new TrkrClusterv4;
	//auto clus = std::make_unique<TrkrClusterv3>();
	clus->setAdc(adc_sum);  
	clus->setOverlap(ntouch);
	clus->setEdge(nedge);
	clus->setPhiSize(phisize);
	clus->setZSize(zsize);
	clus->setSubSurfKey(subsurfkey);      
	clus->setLocalX(localPos(0));
	clus->setLocalY(localPos(1));
	my_data.cluster_vector.push_back(clus);
	
      }
      
      //      if(my_data.do_assoc && my_data.clusterhitassoc){
      if(my_data.do_assoc)
	{
        // get cluster index in vector. It is used to store associations, and build relevant cluster keys when filling the containers
        uint32_t index = my_data.cluster_vector.size()-1;
        for (unsigned int i = 0; i < hitkeyvec.size(); i++){
          my_data.association_vector.emplace_back(index, hitkeyvec[i]);
        }
      }
      hitkeyvec.clear();
    }
  
  void *ProcessSector(void *threadarg) {

    auto my_data = static_cast<thread_data*>(threadarg);
    const auto& pedestal  = my_data->pedestal;
    const auto& phibins   = my_data->phibins;
    const auto& phioffset = my_data->phioffset;
    const auto& zbins     = my_data->zbins ;
    const auto& zoffset   = my_data->zoffset ;
    const auto& layer   = my_data->layer ;
	
    TrkrHitSet *hitset = my_data->hitset;
    TrkrHitSet::ConstRange hitrangei = hitset->getHits();
    
    // for convenience, create a 2D vector to store adc values in and initialize to zero
    std::vector<std::vector<unsigned short>> adcval(phibins, std::vector<unsigned short>(zbins, 0));
    std::multimap<unsigned short, ihit> all_hit_map;
    std::vector<ihit> hit_vect;

    int zbinmax = 498;
    int zbinmin = 0;
    if(my_data->do_wedge_emulation){
      if(layer>=7 && layer <22){
	int etacut = 249 - ((50+(layer-7))/105.5)*249;
	zbinmin = etacut;
	zbinmax -= etacut;
      }
      if(layer>=22 && layer <=48){
	int etacut = 249 - ((65+((40.5/26)*(layer-22)))/105.5)*249;
	zbinmin = etacut;
	zbinmax -= etacut;
      }
      // std::cout << " layer: " << layer << " zbin limit " << zbinmin << " | " << zbinmax <<std::endl;
    }
    for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
	 hitr != hitrangei.second;
	 ++hitr){
      if( TpcDefs::getPad(hitr->first) - phioffset < 0 ){
	//std::cout << "WARNING phibin out of range: " << TpcDefs::getPad(hitr->first) - phioffset << " | " << phibins << std::endl;
	continue;
      }
      if( TpcDefs::getTBin(hitr->first) - zoffset < 0 ){
	//std::cout << "WARNING zbin out of range: " << TpcDefs::getTBin(hitr->first) - zoffset  << " | " << zbins <<std::endl;
      }
      unsigned short phibin = TpcDefs::getPad(hitr->first) - phioffset;
      unsigned short zbin = TpcDefs::getTBin(hitr->first) - zoffset;
      unsigned short zbinorg = TpcDefs::getTBin(hitr->first);
      if(phibin>=phibins){
	//std::cout << "WARNING phibin out of range: " << phibin << " | " << phibins << std::endl;
	continue;
      }
      if(zbin>=zbins){
	//std::cout << "WARNING z bin out of range: " << zbin << " | " << zbins << std::endl;
	continue;
      }
      if(zbinorg>zbinmax||zbinorg<zbinmin)
	continue;
      float_t fadc = (hitr->second->getAdc()) - pedestal; // proper int rounding +0.5
      unsigned short adc = 0;
      if(fadc>0) adc =  (unsigned short) fadc;
      if(phibin >= phibins) continue;
      if(zbin   >= zbins) continue; // zbin is unsigned int, <0 cannot happen
      
      if(adc>0){
	if(adc>5){
	  ihit  thisHit;
	
	  thisHit.iphi = phibin;
	  thisHit.iz = zbin;
	  thisHit.adc = adc;
	  thisHit.edge = 0;
	  all_hit_map.insert(std::make_pair(adc, thisHit));
	}
	adcval[phibin][zbin] = (unsigned short) adc;
      }
    }
    
    while(all_hit_map.size()>0){
      
      auto iter = all_hit_map.rbegin();
      if(iter == all_hit_map.rend()){
	break;
      }
      ihit hiHit = iter->second;
      int iphi = hiHit.iphi;
      int iz = hiHit.iz;
      
      //put all hits in the all_hit_map (sorted by adc)
      //start with highest adc hit
      // -> cluster around it and get vector of hits
      std::vector<ihit> ihit_list;
      int ntouch = 0;
      int nedge  =0;
      get_cluster(iphi, iz, *my_data, adcval, ihit_list, ntouch, nedge );
      
      // -> calculate cluster parameters
      // -> add hits to truth association
      // remove hits from all_hit_map
      // repeat untill all_hit_map empty
      calc_cluster_parameter(ihit_list, *my_data, ntouch, nedge );
      remove_hits(ihit_list,all_hit_map, adcval);
      ihit_list.clear();
    }

    pthread_exit(nullptr);
  }
}

TpcClusterizer::TpcClusterizer(const std::string &name)
  : SubsysReco(name)
{}

bool TpcClusterizer::is_in_sector_boundary(int phibin, int sector, PHG4CylinderCellGeom *layergeom) const
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
      std::cout << " local maximum is in sector fiducial boundary: layer " << layer << " radius " << radius << " sector " << sector 
      << " PhiBins " << PhiBins << " sector_fiducial_bins " << sector_fiducial_bins
      << " PhiBinSize " << PhiBinSize << " phibin " << phibin << " sector_lo " << sector_lo << " sector_hi " << sector_hi << std::endl;  
      */
    }

  return reject_it;
}


int TpcClusterizer::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // Create the Cluster node if required
  auto trkrclusters = findNode::getClass<TrkrClusterContainer>(dstNode, "TRKR_CLUSTER");
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

    trkrclusters = new TrkrClusterContainerv4;
    PHIODataNode<PHObject> *TrkrClusterContainerNode =
        new PHIODataNode<PHObject>(trkrclusters, "TRKR_CLUSTER", "PHObject");
    DetNode->addNode(TrkrClusterContainerNode);
  }

  auto clusterhitassoc = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
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

    clusterhitassoc = new TrkrClusterHitAssocv3;
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(clusterhitassoc, "TRKR_CLUSTERHITASSOC", "PHObject");
    DetNode->addNode(newNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcClusterizer::process_event(PHCompositeNode *topNode)
{
  //  int print_layer = 18;

  if (Verbosity() > 1000)
    std::cout << "TpcClusterizer::Process_Event" << std::endl;

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
  m_clusterlist = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!m_clusterlist)
  {
    std::cout << PHWHERE << " ERROR: Can't find TRKR_CLUSTER." << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // get node for cluster hit associations
  m_clusterhitassoc = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  if (!m_clusterhitassoc)
  {
    std::cout << PHWHERE << " ERROR: Can't find TRKR_CLUSTERHITASSOC" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PHG4CylinderCellGeomContainer *geom_container =
      findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!geom_container)
  {
    std::cout << PHWHERE << "ERROR: Can't find node CYLINDERCELLGEOM_SVTX" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_tGeometry = findNode::getClass<ActsGeometry>(topNode,
						 "ActsGeometry");
  if(!m_tGeometry)
    {
      std::cout << PHWHERE
		<< "ActsGeometry not found on node tree. Exiting"
		<< std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

  // The hits are stored in hitsets, where each hitset contains all hits in a given TPC readout (layer, sector, side), so clusters are confined to a hitset
  // The TPC clustering is more complicated than for the silicon, because we have to deal with overlapping clusters

  // loop over the TPC HitSet objects
  TrkrHitSetContainer::ConstRange hitsetrange = m_hits->getHitSets(TrkrDefs::TrkrId::tpcId);
  const int num_hitsets = std::distance(hitsetrange.first,hitsetrange.second);

  // create structure to store given thread and associated data
  struct thread_pair_t
  {
    pthread_t thread;
    thread_data data;
  };
  
  // create vector of thread pairs and reserve the right size upfront to avoid reallocation
  std::vector<thread_pair_t> threads;
  threads.reserve( num_hitsets );

  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  
  if (pthread_mutex_init(&mythreadlock, nullptr) != 0)
    {
      printf("\n mutex init failed\n");
      return 1;
    }
  int count = 0;
  for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
       hitsetitr != hitsetrange.second;
       ++hitsetitr)
  {
    //if(count>0)continue;
    //    const auto hitsetid = hitsetitr->first;
    //    std::cout << " starting thread # " << count << std::endl;
    TrkrHitSet *hitset = hitsetitr->second;
    unsigned int layer = TrkrDefs::getLayer(hitsetitr->first);
    int side = TpcDefs::getSide(hitsetitr->first);
    unsigned int sector= TpcDefs::getSectorId(hitsetitr->first);
    PHG4CylinderCellGeom *layergeom = geom_container->GetLayerCellGeom(layer);
    
    // instanciate new thread pair, at the end of thread vector
    thread_pair_t& thread_pair = threads.emplace_back();
    
    thread_pair.data.layergeom = layergeom;
    thread_pair.data.hitset = hitset;
    thread_pair.data.layer = layer;
    thread_pair.data.pedestal = pedestal;
    thread_pair.data.sector = sector;
    thread_pair.data.side = side;
    thread_pair.data.do_assoc = do_hit_assoc;
    thread_pair.data.do_wedge_emulation = do_wedge_emulation;
    thread_pair.data.tGeometry = m_tGeometry;
    thread_pair.data.maxHalfSizeZ =  MaxClusterHalfSizeZ;
    thread_pair.data.maxHalfSizePhi = MaxClusterHalfSizePhi;
    thread_pair.data.m_drift_velocity_scale = m_drift_velocity_scale;
    thread_pair.data.par0_neg = par0_neg;
    thread_pair.data.par0_pos = par0_pos;
    thread_pair.data.cluster_version = cluster_version;
    thread_pair.data.verbosity = Verbosity();
    
    unsigned short NPhiBins = (unsigned short) layergeom->get_phibins();
    unsigned short NPhiBinsSector = NPhiBins/12;
    unsigned short NZBins = (unsigned short)layergeom->get_zbins();
    //unsigned short NZBinsSide = NZBins/2;
    unsigned short NZBinsSide = NZBins;
    unsigned short NZBinsMin = 0;
    unsigned short PhiOffset = NPhiBinsSector * sector;

    /*
    if (side == 0){
      NZBinsMin = 0;
    }
    else{
      NZBinsMin = NZBins / 2;
    }
    */
    unsigned short ZOffset = NZBinsMin;

    thread_pair.data.phibins   = NPhiBinsSector;
    thread_pair.data.phioffset = PhiOffset;
    thread_pair.data.zbins     = NZBinsSide;
    thread_pair.data.zoffset   = ZOffset ;

    int rc = pthread_create(&thread_pair.thread, &attr, ProcessSector, (void *)&thread_pair.data);
    if (rc) {
      std::cout << "Error:unable to create thread," << rc << std::endl;
    }
    count++;
  }
  
  pthread_attr_destroy(&attr);
  count =0;
  // wait for completion of all threads
  for( const auto& thread_pair:threads )
  { 
    int rc2 = pthread_join(thread_pair.thread, nullptr);
    if (rc2) 
    { std::cout << "Error:unable to join," << rc2 << std::endl; }

    // get the hitsetkey from thread data
    const auto& data( thread_pair.data );
    const auto hitsetkey = TpcDefs::genHitSetKey( data.layer, data.sector, data.side );      

    // copy clusters to map
    for( uint32_t index = 0; index < data.cluster_vector.size(); ++index )
    {
      // generate cluster key
      const auto ckey = TrkrDefs::genClusKey( hitsetkey, index );

      // get cluster
      auto cluster = data.cluster_vector[index];

      // insert in map
      m_clusterlist->addClusterSpecifyKey(ckey, cluster);
    }

    // copy hit associations to map
    for( const auto& [index,hkey]:thread_pair.data.association_vector)
    { 
      // generate cluster key
      const auto ckey = TrkrDefs::genClusKey( hitsetkey, index );

      // add to association table
      m_clusterhitassoc->addAssoc(ckey,hkey); 
    }

  }

  if (Verbosity() > 0)
    std::cout << "TPC Clusterizer found " << m_clusterlist->size() << " Clusters "  << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcClusterizer::End(PHCompositeNode */*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
