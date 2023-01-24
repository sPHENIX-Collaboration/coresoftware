#include "TpcClusterizer.h"

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrClusterContainerv4.h>
#include <trackbase/TrkrClusterv3.h>
#include <trackbase/TrkrClusterv4.h>
#include <trackbase/TrkrClusterHitAssocv3.h>
#include <trackbase/TrkrDefs.h>  // for hitkey, getLayer
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <trackbase/RawHit.h>
#include <trackbase/RawHitSet.h>
#include <trackbase/RawHitSetv1.h>
#include <trackbase/RawHitSetContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>                         // for SubsysReco

#include <g4detectors/PHG4TpcCylinderGeom.h>
#include <g4detectors/PHG4TpcCylinderGeomContainer.h>

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
    unsigned short it = 0;
    unsigned short adc = 0;
    unsigned short edge = 0;
  };

  struct thread_data 
  {
    PHG4TpcCylinderGeom *layergeom = nullptr;
    TrkrHitSet *hitset = nullptr;
    RawHitSetv1 *rawhitset = nullptr;
    ActsGeometry *tGeometry = nullptr;
    unsigned int layer = 0;
    int side = 0;
    unsigned int sector = 0;
    float radius = 0;
    float drift_velocity = 0;
    unsigned short pads_per_sector = 0;
    float phistep = 0;
    float pedestal = 0;
    float threshold = 0;
    bool do_assoc = true;
    bool do_wedge_emulation = true;
    bool do_singles = true;
    unsigned short phibins = 0;
    unsigned short phioffset = 0;
    unsigned short tbins = 0;
    unsigned short toffset = 0;
    unsigned short maxHalfSizeT = 0;
    unsigned short maxHalfSizePhi = 0;
    double m_tdriftmax = 0;
    double sampa_tbias = 0;
    int cluster_version = 4;
    std::vector<assoc> association_vector;
    std::vector<TrkrCluster*> cluster_vector;
    int verbosity = 0;
  };
  
  pthread_mutex_t mythreadlock;
  
  void remove_hit(double adc, int phibin, int tbin, int edge, std::multimap<unsigned short, ihit> &all_hit_map, std::vector<std::vector<unsigned short>> &adcval)
  {
    typedef std::multimap<unsigned short, ihit>::iterator hit_iterator;
    std::pair<hit_iterator, hit_iterator> iterpair = all_hit_map.equal_range(adc);
    hit_iterator it = iterpair.first;
    for (; it != iterpair.second; ++it) {
      if (it->second.iphi == phibin && it->second.it == tbin) { 
	all_hit_map.erase(it);
	break;
      }
    }
    if(edge)
      adcval[phibin][tbin] = USHRT_MAX;
    else
      adcval[phibin][tbin] = 0;
  }
  
  void remove_hits(std::vector<ihit> &ihit_list, std::multimap<unsigned short, ihit> &all_hit_map,std::vector<std::vector<unsigned short>> &adcval)
  {
    for(auto iter = ihit_list.begin(); iter != ihit_list.end();++iter){
      unsigned short adc    = iter->adc; 
      unsigned short phibin = iter->iphi;
      unsigned short tbin   = iter->it;
      unsigned short edge   = iter->edge;
      remove_hit(adc,phibin,tbin,edge,all_hit_map,adcval);
    }
  }
  
  void find_t_range(int phibin, int tbin, const thread_data& my_data, const std::vector<std::vector<unsigned short>> &adcval, int& tdown, int& tup, int &touch, int &edge){
	
    const int FitRangeT= (int) my_data.maxHalfSizeT;
    const int NTBinsMax = (int) my_data.tbins;
    tup = 0;
    tdown = 0;
    for(int it=0; it< FitRangeT; it++){
      int ct = tbin + it;
      
      if(ct <= 0 || ct >= NTBinsMax){
	// tup = it;
	edge++;
	break; // truncate edge
      }
      
      if(adcval[phibin][ct] <= 0) {
	break;
      }
      if(adcval[phibin][ct] == USHRT_MAX) {
	touch++;
	break;
      }
      //check local minima and break at minimum.
      if(ct<NTBinsMax-4){//make sure we stay clear from the edge
	if(adcval[phibin][ct]+adcval[phibin][ct+1] < 
	   adcval[phibin][ct+2]+adcval[phibin][ct+3]){//rising again
	  tup = it+1;
	  touch++;
	  break;
	      }
      }
      tup = it;
    }
    for(int it=0; it< FitRangeT; it++){
      int ct = tbin - it;
      if(ct <= 0 || ct >= NTBinsMax){
	//      tdown = it;
	edge++;
	break; // truncate edge
      }
      if(adcval[phibin][ct] <= 0) {
	break;
      }
      if(adcval[phibin][ct] == USHRT_MAX) {
	touch++;
	break;
      }
      if(ct>4){//make sure we stay clear from the edge
	if(adcval[phibin][ct]+adcval[phibin][ct-1] < 
	   adcval[phibin][ct-2]+adcval[phibin][ct-3]){//rising again
	  tdown = it+1;
	  touch++;
	  break;
	}
      }
      tdown = it;
    }
    return;
  }
	
  void find_phi_range(int phibin, int tbin, const thread_data& my_data, const std::vector<std::vector<unsigned short>> &adcval, int& phidown, int& phiup, int &touch, int &edge)
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
      if(adcval[cphi][tbin] <= 0) {
	// phiup = iphi;
	break;
      }
      if(adcval[cphi][tbin] == USHRT_MAX) {
	touch++;
	break;
      }
      //check local minima and break at minimum.
      if(cphi<NPhiBinsMax-4){//make sure we stay clear from the edge
	if(adcval[cphi][tbin]+adcval[cphi+1][tbin] < 
	   adcval[cphi+2][tbin]+adcval[cphi+3][tbin]){//rising again
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
      
      if(adcval[cphi][tbin] <= 0) {
	//phidown = iphi;
	break;
      }
      if(adcval[cphi][tbin] == USHRT_MAX) {
	touch++;
	break;
      }
      if(cphi>4){//make sure we stay clear from the edge
	if(adcval[cphi][tbin]+adcval[cphi-1][tbin] < 
	   adcval[cphi-2][tbin]+adcval[cphi-3][tbin]){//rising again
	  phidown = iphi+1;
	  touch++;
	  break;
	}
      }
      phidown = iphi;
    }
    return;
  }
	
  void get_cluster(int phibin, int tbin, const thread_data& my_data, const std::vector<std::vector<unsigned short>> &adcval, std::vector<ihit> &ihit_list, int &touch, int &edge)
	{
	  // search along phi at the peak in t
	 
	  int tup =0;
	  int tdown =0;
	  find_t_range(phibin, tbin, my_data, adcval, tdown, tup, touch, edge);
	  //now we have the t extent of the cluster, go find the phi edges
	
	  for(int it=tbin - tdown ; it<= tbin + tup; it++){
	    int phiup = 0;
	    int phidown = 0;
	    find_phi_range(phibin, it, my_data, adcval, phidown, phiup, touch, edge);
	    for (int iphi = phibin - phidown; iphi <= (phibin + phiup); iphi++){
	      if(adcval[iphi][it]>0 && adcval[iphi][it]!=USHRT_MAX){
		ihit hit;
		hit.iphi = iphi;
		hit.it = it;
		hit.adc = adcval[iphi][it];
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
      //
      // get z range from layer geometry
      /* these are used for rescaling the drift velocity */
      //const double z_min = -105.5;
      //const double z_max = 105.5;
      // std::cout << "calc clus" << std::endl;    
      // loop over the hits in this cluster
      double t_sum = 0.0;
      //double phi_sum = 0.0;
      double adc_sum = 0.0;
      double t2_sum = 0.0;
      // double phi2_sum = 0.0;

      double iphi_sum = 0.0;
      double iphi2_sum = 0.0;

      double radius = my_data.layergeom->get_radius();  // returns center of layer
      
      int phibinhi = -1;
      int phibinlo = 666666;
      int tbinhi = -1;
      int tbinlo = 666666;
      int clus_size = ihit_list.size();
      
      if(clus_size == 1) return;
      //      std::cout << "process list" << std::endl;    
      std::vector<TrkrDefs::hitkey> hitkeyvec;
      for(auto iter = ihit_list.begin(); iter != ihit_list.end();++iter){
	double adc = iter->adc; 
	
	if (adc <= 0) continue;
	
	int iphi = iter->iphi + my_data.phioffset;
	int it   = iter->it + my_data.toffset;
	if(iphi > phibinhi) phibinhi = iphi;
	if(iphi < phibinlo) phibinlo = iphi;
	if(it > tbinhi) tbinhi = it;
	if(it < tbinlo) tbinlo = it;
	
	// update phi sums
	//	double phi_center = my_data.layergeom->get_phicenter(iphi);
	
	//phi_sum += phi_center * adc;
	//phi2_sum += square(phi_center)*adc;
	//	std::cout << "phi_center: " << phi_center << " adc: " << adc <<std::endl;
	iphi_sum += iphi * adc;
	iphi2_sum += square(iphi)*adc;

	// update t sums
	double t = my_data.layergeom->get_zcenter(it);
	t_sum += t*adc;
	t2_sum += square(t)*adc;
	
	adc_sum += adc;
	
	// capture the hitkeys for all adc values above a certain threshold
	TrkrDefs::hitkey hitkey = TpcDefs::genHitKey(iphi, it);
	// if(adc>5)
	hitkeyvec.push_back(hitkey);
      }
      //      std::cout << "done process list" << std::endl;
      if (adc_sum < 10){
	hitkeyvec.clear();
	return;  // skip obvious noise "clusters"
      }  
      // This is the global position
      double clusiphi = iphi_sum / adc_sum;
      double clusphi = my_data.layergeom->get_phi(clusiphi);

      float clusx = radius * cos(clusphi);
      float clusy = radius * sin(clusphi);
      double clust = t_sum / adc_sum;
      // needed for surface identification
      double zdriftlength = clust * my_data.tGeometry->get_drift_velocity();
      // convert z drift length to z position in the TPC
      double clusz  =  my_data.m_tdriftmax * my_data.tGeometry->get_drift_velocity() - zdriftlength; 
      if(my_data.side == 0) 
	clusz = -clusz;

      const double phi_cov = (iphi2_sum/adc_sum - square(clusiphi))* pow(my_data.layergeom->get_phistep(),2);
      const double t_cov = t2_sum/adc_sum - square(clust);

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
      
      const double t_err_square = (tbinhi == tbinlo) ?
	square(my_data.layergeom->get_zstep())/12:
	t_cov/(adc_sum*0.14);
      
      char tsize = tbinhi - tbinlo + 1;
      char phisize = phibinhi - phibinlo + 1;
      // phi_cov = (weighted mean of dphi^2) - (weighted mean of dphi)^2,  which is essentially the weighted mean of dphi^2. The error is then:
      // e_phi = sigma_dphi/sqrt(N) = sqrt( sigma_dphi^2 / N )  -- where N is the number of samples of the distribution with standard deviation sigma_dphi
      //    - N is the number of electrons that drift to the readout plane
      // We have to convert (sum of adc units for all bins in the cluster) to number of ionization electrons N
      // Conversion gain is 20 mV/fC - relates total charge collected on pad to PEAK voltage out of ADC. The GEM gain is assumed to be 2000
      // To get equivalent charge per T bin, so that summing ADC input voltage over all T bins returns total input charge, divide voltages by 2.4 for 80 ns SAMPA
      // Equivalent charge per T bin is then  (ADU x 2200 mV / 1024) / 2.4 x (1/20) fC/mV x (1/1.6e-04) electrons/fC x (1/2000) = ADU x 0.14

      // SAMPA shaping bias correction
      clust = clust + my_data.sampa_tbias;

      /// convert to Acts units
      global *= Acts::UnitConstants::cm;
      //std::cout << "transform" << std::endl;
      Acts::Vector3 local = surface->transform(my_data.tGeometry->geometry().getGeoContext()).inverse() * global;
      local /= Acts::UnitConstants::cm;     
      //std::cout << "done transform" << std::endl;
      // we need the cluster key and all associated hit keys (note: the cluster key includes the hitset key)
      
      if(my_data.cluster_version==3){
	
	// Fill in the cluster details
	//================
	auto clus = new TrkrClusterv3;
	//auto clus = std::make_unique<TrkrClusterv3>();
	clus->setAdc(adc_sum);      
	clus->setSubSurfKey(subsurfkey);      
	clus->setLocalX(local(0));
	clus->setLocalY(clust);
	clus->setActsLocalError(0,0, phi_err_square);
	clus->setActsLocalError(1,0, 0);
	clus->setActsLocalError(0,1, 0);
	clus->setActsLocalError(1,1, t_err_square * pow(my_data.tGeometry->get_drift_velocity(),2));
	my_data.cluster_vector.push_back(clus);
      }else if(my_data.cluster_version==4){
	//	std::cout << "clus num" << my_data.cluster_vector.size() << " X " << local(0) << " Y " << clust << std::endl;
	if(sqrt(phi_err_square) > 0.01){
	auto clus = new TrkrClusterv4;
	//auto clus = std::make_unique<TrkrClusterv3>();
	clus->setAdc(adc_sum);  
	clus->setOverlap(ntouch);
	clus->setEdge(nedge);
	clus->setPhiSize(phisize);
	clus->setZSize(tsize);
	clus->setSubSurfKey(subsurfkey);      
	clus->setLocalX(local(0));
	clus->setLocalY(clust);
	my_data.cluster_vector.push_back(clus);
	}
      }
      //std::cout << "end clus out" << std::endl;
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
      //      std::cout << "done calc" << std::endl;
    }
  
  void ProcessSectorData(thread_data* my_data) {

    const auto& pedestal  = my_data->pedestal;
    const auto& phibins   = my_data->phibins;
    const auto& phioffset = my_data->phioffset;
    const auto& tbins     = my_data->tbins ;
    const auto& toffset   = my_data->toffset ;
    const auto& layer   = my_data->layer ;
    //    int nhits = 0;
    // for convenience, create a 2D vector to store adc values in and initialize to zero
    std::vector<std::vector<unsigned short>> adcval(phibins, std::vector<unsigned short>(tbins, 0));
    std::multimap<unsigned short, ihit> all_hit_map;
    std::vector<ihit> hit_vect;

    int tbinmax = 498;
    int tbinmin = 0;
    if(my_data->do_wedge_emulation){
      if(layer>=7 && layer <22){
	int etacut = 249 - ((50+(layer-7))/105.5)*249;
	tbinmin = etacut;
	tbinmax -= etacut;
      }
      if(layer>=22 && layer <=48){
	int etacut = 249 - ((65+((40.5/26)*(layer-22)))/105.5)*249;
	tbinmin = etacut;
	tbinmax -= etacut;
      }
    }

    if( my_data->hitset!=nullptr){
      TrkrHitSet *hitset = my_data->hitset;
      TrkrHitSet::ConstRange hitrangei = hitset->getHits();
      
      for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
	   hitr != hitrangei.second;
	   ++hitr){
	
	if( TpcDefs::getPad(hitr->first) - phioffset < 0 ){
	  //std::cout << "WARNING phibin out of range: " << TpcDefs::getPad(hitr->first) - phioffset << " | " << phibins << std::endl;
	  continue;
	}
	if( TpcDefs::getTBin(hitr->first) - toffset < 0 ){
	  //std::cout << "WARNING tbin out of range: " << TpcDefs::getTBin(hitr->first) - toffset  << " | " << tbins <<std::endl;
	}
	unsigned short phibin = TpcDefs::getPad(hitr->first) - phioffset;
	unsigned short tbin = TpcDefs::getTBin(hitr->first) - toffset;
	unsigned short tbinorg = TpcDefs::getTBin(hitr->first);
	if(phibin>=phibins){
	  //std::cout << "WARNING phibin out of range: " << phibin << " | " << phibins << std::endl;
	  continue;
	}
	if(tbin>=tbins){
	  //std::cout << "WARNING z bin out of range: " << tbin << " | " << tbins << std::endl;
	  continue;
	}
	if(tbinorg>tbinmax||tbinorg<tbinmin)
	  continue;
	float_t fadc = (hitr->second->getAdc()) - pedestal; // proper int rounding +0.5
	unsigned short adc = 0;
	if(fadc>0) adc =  (unsigned short) fadc;
	if(phibin >= phibins) continue;
	if(tbin   >= tbins) continue; // tbin is unsigned int, <0 cannot happen
	
	if(adc>0){
	  if(adc>(5+my_data->threshold)){
	    ihit  thisHit;
	    
	    thisHit.iphi = phibin;
	    thisHit.it = tbin;
	    thisHit.adc = adc;
	    thisHit.edge = 0;
	    all_hit_map.insert(std::make_pair(adc, thisHit));
	  }
	  if(adc>my_data->threshold){
	    adcval[phibin][tbin] = (unsigned short) adc;
	  }
	}
      }
    }else  if( my_data->rawhitset!=nullptr){
      RawHitSetv1 *hitset = my_data->rawhitset;
      /*std::cout << "Layer: " << my_data->layer 
		<< "Side: " << my_data->side
		<< "Sector: " << my_data->sector
		<< " nhits:  " << hitset.size()
		<< std::endl;
      */
      for(int nphi= 0; nphi < phibins;nphi++){
	//	nhits += hitset->m_tpchits[nphi].size();
	if(hitset->m_tpchits[nphi].size()==0) continue;

	int pindex = 0;
	for(unsigned int nt = 0;nt<hitset->m_tpchits[nphi].size();nt++){
	  unsigned short val = hitset->m_tpchits[nphi][nt];
	  
	  if(val==0)
	    pindex++;
	  else{
	    if(nt==0){
	      if(val>5){
		ihit  thisHit;
		thisHit.iphi = nphi;
		thisHit.it = pindex;
		thisHit.adc = val;
		thisHit.edge = 0;
		all_hit_map.insert(std::make_pair(val, thisHit));
	      }
	      adcval[nphi][pindex++]=val;
	    }else{
	      if((hitset->m_tpchits[nphi][nt-1]==0)&&(hitset->m_tpchits[nphi][nt+1]==0))//found zero count
		pindex+=val;
	      else{
		if(val>5){
		  ihit  thisHit;
		  thisHit.iphi = nphi;
		  thisHit.it = pindex;
		  thisHit.adc = val;
		  thisHit.edge = 0;
		  all_hit_map.insert(std::make_pair(val, thisHit));
		}
		adcval[nphi][pindex++]=val;
	      }
	    }
	  }
	}
      }
    }
   
    if(my_data->do_singles){
      for(auto ahit:all_hit_map){
	ihit hiHit = ahit.second;
	int iphi = hiHit.iphi;
	int it = hiHit.it;
	unsigned short edge = hiHit.edge;
	double adc = hiHit.adc;
	if(it>0&&it<tbins){
	  if(adcval[iphi][it-1]==0&&
	     adcval[iphi][it+1]==0){
	    remove_hit(adc, iphi, it, edge, all_hit_map, adcval);
	    
	  }
	}
      }
    }
    
    // std::cout << "done filling " << std::endl;
    while(all_hit_map.size()>0){
      //std::cout << "all hit map size: " << all_hit_map.size() << std::endl;
      auto iter = all_hit_map.rbegin();
      if(iter == all_hit_map.rend()){
	break;
      }
      ihit hiHit = iter->second;
      int iphi = hiHit.iphi;
      int it = hiHit.it;
      //put all hits in the all_hit_map (sorted by adc)
      //start with highest adc hit
      // -> cluster around it and get vector of hits
      std::vector<ihit> ihit_list;
      int ntouch = 0;
      int nedge  =0;
      get_cluster(iphi, it, *my_data, adcval, ihit_list, ntouch, nedge );
      
      // -> calculate cluster parameters
      // -> add hits to truth association
      // remove hits from all_hit_map
      // repeat untill all_hit_map empty
      calc_cluster_parameter(ihit_list, *my_data, ntouch, nedge );
      remove_hits(ihit_list,all_hit_map, adcval);
      ihit_list.clear();
    }
    /*    if( my_data->rawhitset!=nullptr){
      RawHitSetv1 *hitset = my_data->rawhitset;
      std::cout << "Layer: " << my_data->layer 
		<< " Side: " << my_data->side
		<< " Sector: " << my_data->sector
		<< " nhits:  " << hitset->size()
		<< " nhits coutn :  " << nhits
		<< " nclus: " << my_data->cluster_vector.size()
		<< std::endl;
    }
    */
    // pthread_exit(nullptr);
  }
  void *ProcessSector(void *threadarg) {

    auto my_data = static_cast<thread_data*>(threadarg);
    ProcessSectorData(my_data);
    pthread_exit(nullptr);
  }
}

TpcClusterizer::TpcClusterizer(const std::string &name)
  : SubsysReco(name)
{}

bool TpcClusterizer::is_in_sector_boundary(int phibin, int sector, PHG4TpcCylinderGeom *layergeom) const
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
  if(!do_read_raw){
    // get node containing the digitized hits
    m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
    if (!m_hits)
      {
	std::cout << PHWHERE << "ERROR: Can't find node TRKR_HITSET" << std::endl;
	return Fun4AllReturnCodes::ABORTRUN;
      }
  }else{
    // get node containing the digitized hits
    m_rawhits = findNode::getClass<RawHitSetContainer>(topNode, "TRKR_RAWHITSET");
    if (!m_rawhits)
      {
	std::cout << PHWHERE << "ERROR: Can't find node TRKR_HITSET" << std::endl;
	return Fun4AllReturnCodes::ABORTRUN;
      }
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

  PHG4TpcCylinderGeomContainer *geom_container =
      findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
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

  TrkrHitSetContainer::ConstRange hitsetrange;
  RawHitSetContainer::ConstRange rawhitsetrange;
  int num_hitsets = 0;

  if(!do_read_raw){
    hitsetrange = m_hits->getHitSets(TrkrDefs::TrkrId::tpcId);
    num_hitsets = std::distance(hitsetrange.first,hitsetrange.second);
  }else{
    rawhitsetrange = m_rawhits->getHitSets(TrkrDefs::TrkrId::tpcId);
    num_hitsets = std::distance(rawhitsetrange.first,rawhitsetrange.second);
  }

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

  if(!do_read_raw){
    for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
	 hitsetitr != hitsetrange.second;
	 ++hitsetitr)
      {
      	//if(count>0)continue;
	TrkrHitSet *hitset = hitsetitr->second;
	unsigned int layer = TrkrDefs::getLayer(hitsetitr->first);
	int side = TpcDefs::getSide(hitsetitr->first);
	unsigned int sector= TpcDefs::getSectorId(hitsetitr->first);
	PHG4TpcCylinderGeom *layergeom = geom_container->GetLayerCellGeom(layer);
	
	// instanciate new thread pair, at the end of thread vector
	thread_pair_t& thread_pair = threads.emplace_back();
	
	thread_pair.data.layergeom = layergeom;
	thread_pair.data.hitset = hitset;
	thread_pair.data.rawhitset = nullptr;
	thread_pair.data.layer = layer;
	thread_pair.data.pedestal = pedestal;
	thread_pair.data.threshold = threshold;
	thread_pair.data.sector = sector;
	thread_pair.data.side = side;
	thread_pair.data.do_assoc = do_hit_assoc;
	thread_pair.data.do_wedge_emulation = do_wedge_emulation;
	thread_pair.data.do_singles = do_singles;
	thread_pair.data.tGeometry = m_tGeometry;
	thread_pair.data.maxHalfSizeT =  MaxClusterHalfSizeT;
	thread_pair.data.maxHalfSizePhi = MaxClusterHalfSizePhi;
	thread_pair.data.sampa_tbias = m_sampa_tbias;
	thread_pair.data.cluster_version = cluster_version;
	thread_pair.data.verbosity = Verbosity();
	
	unsigned short NPhiBins = (unsigned short) layergeom->get_phibins();
	unsigned short NPhiBinsSector = NPhiBins/12;
	unsigned short NTBins = (unsigned short)layergeom->get_zbins();
	unsigned short NTBinsSide = NTBins;
	unsigned short NTBinsMin = 0;
	unsigned short PhiOffset = NPhiBinsSector * sector;
	unsigned short TOffset = NTBinsMin;
	
	m_tdriftmax = AdcClockPeriod * NTBins / 2.0;  
	thread_pair.data.m_tdriftmax = m_tdriftmax;
	
	thread_pair.data.phibins   = NPhiBinsSector;
	thread_pair.data.phioffset = PhiOffset;
	thread_pair.data.tbins     = NTBinsSide;
	thread_pair.data.toffset   = TOffset ;

	thread_pair.data.radius = layergeom->get_radius();
	thread_pair.data.drift_velocity = m_tGeometry->get_drift_velocity();
	thread_pair.data.pads_per_sector = 0;
	thread_pair.data.phistep = 0;
	int rc;
	rc = pthread_create(&thread_pair.thread, &attr, ProcessSector, (void *)&thread_pair.data);

	if (rc) {
	  std::cout << "Error:unable to create thread," << rc << std::endl;
	}
	if(do_sequential){
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
	count++;
      }
  }else{

    for (RawHitSetContainer::ConstIterator hitsetitr = rawhitsetrange.first;
	 hitsetitr != rawhitsetrange.second;
	 ++hitsetitr){ 
      //	if(count>0)continue;
	//    const auto hitsetid = hitsetitr->first;
      //	std::cout << " starting thread # " << count << std::endl;

	RawHitSet *hitset = hitsetitr->second;
      unsigned int layer = TrkrDefs::getLayer(hitsetitr->first);
      int side = TpcDefs::getSide(hitsetitr->first);
      unsigned int sector= TpcDefs::getSectorId(hitsetitr->first);
      PHG4TpcCylinderGeom *layergeom = geom_container->GetLayerCellGeom(layer);
      
      // instanciate new thread pair, at the end of thread vector
      thread_pair_t& thread_pair = threads.emplace_back();
      
      thread_pair.data.layergeom = layergeom;
      thread_pair.data.hitset = nullptr;
      thread_pair.data.rawhitset = dynamic_cast<RawHitSetv1 *>(hitset);
      thread_pair.data.layer = layer;
      thread_pair.data.pedestal = pedestal;
      thread_pair.data.sector = sector;
      thread_pair.data.side = side;
      thread_pair.data.do_assoc = do_hit_assoc;
      thread_pair.data.do_wedge_emulation = do_wedge_emulation;
      thread_pair.data.tGeometry = m_tGeometry;
      thread_pair.data.maxHalfSizeT =  MaxClusterHalfSizeT;
      thread_pair.data.maxHalfSizePhi = MaxClusterHalfSizePhi;
      thread_pair.data.sampa_tbias = m_sampa_tbias;
      thread_pair.data.cluster_version = cluster_version;
      thread_pair.data.verbosity = Verbosity();

      unsigned short NPhiBins = (unsigned short) layergeom->get_phibins();
      unsigned short NPhiBinsSector = NPhiBins/12;
      unsigned short NTBins = (unsigned short)layergeom->get_zbins();
      unsigned short NTBinsSide = NTBins;
      unsigned short NTBinsMin = 0;
      unsigned short PhiOffset = NPhiBinsSector * sector;
      unsigned short TOffset = NTBinsMin;

      m_tdriftmax = AdcClockPeriod * NTBins / 2.0;  
      thread_pair.data.m_tdriftmax = m_tdriftmax;

      thread_pair.data.phibins   = NPhiBinsSector;
      thread_pair.data.phioffset = PhiOffset;
      thread_pair.data.tbins     = NTBinsSide;
      thread_pair.data.toffset   = TOffset ;

      /*
      PHG4TpcCylinderGeom *testlayergeom = geom_container->GetLayerCellGeom(32);
      for( float iphi = 1408; iphi < 1408+ 128;iphi+=0.1){
	double clusiphi = iphi;
	double clusphi = testlayergeom->get_phi(clusiphi);
	double radius = layergeom->get_radius(); 
	float clusx = radius * cos(clusphi);
	float clusy = radius * sin(clusphi);
	float clusz  = -37.524;
	
	TrkrDefs::hitsetkey tpcHitSetKey = TpcDefs::genHitSetKey( 32,11, 0 );      
	Acts::Vector3 global(clusx, clusy, clusz);
	TrkrDefs::subsurfkey subsurfkey = 0;

	Surface surface = m_tGeometry->get_tpc_surface_from_coords(
								   tpcHitSetKey,
								   global,
								   subsurfkey);
	std::cout << " iphi: " << iphi << " clusphi: " << clusphi << " surfkey " << subsurfkey << std::endl;
	//	std::cout << "surfkey" << subsurfkey << std::endl;
      }
      continue;
      */
      int rc = 0;
      //      if(layer==32)
      rc = pthread_create(&thread_pair.thread, &attr, ProcessSector, (void *)&thread_pair.data);	  
      //      else
      //continue;
      
      if (rc) {
	std::cout << "Error:unable to create thread," << rc << std::endl;
      }
      
      if(do_sequential){
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
      count++;
    }
  }

  
  pthread_attr_destroy(&attr);
  count =0;
  // wait for completion of all threads
  if(!do_sequential){
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
	    //std::cout << "X: " << cluster->getLocalX() << "Y: " << cluster->getLocalY() << std::endl;
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
  }

  if (Verbosity() > 0)
    std::cout << "TPC Clusterizer found " << m_clusterlist->size() << " Clusters "  << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcClusterizer::End(PHCompositeNode */*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
