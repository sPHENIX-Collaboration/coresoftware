#include "TpcSimpleClusterizer.h"

#include <trackbase/TpcDefs.h>

#include <trackbase/TrkrClusterContainerv4.h>
#include <trackbase/TrkrClusterv3.h>
#include <trackbase/TrkrClusterHitAssocv3.h>
#include <trackbase/TrkrDefs.h>  // for hitkey, getLayer
#include <trackbase/TrkrHitv2.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

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
// Terra incognita....
#include <pthread.h>

namespace 
{
  template<class T> inline constexpr T square( const T& x ) { return x*x; }

  using iphiz = std::pair<unsigned short, unsigned short>;
  using ihit = std::pair<unsigned short, iphiz>;
  using assoc = std::pair<uint32_t, TrkrDefs::hitkey> ;
	
  struct thread_data 
  {
    PHG4TpcCylinderGeom *layergeom = nullptr;
    TrkrHitSet *hitset = nullptr;
    ActsGeometry *tGeometry = nullptr;
    unsigned int layer = 0;
    int side = 0;
    unsigned int sector = 0;
    float pedestal = 0;
    bool do_assoc = true;
    unsigned short phibins = 0;
    unsigned short phioffset = 0;
    unsigned short zbins = 0;
    unsigned short zoffset = 0;
    double par0_neg = 0;
    double par0_pos = 0;
    std::vector<assoc> association_vector;
    std::vector<TrkrCluster*> cluster_vector;
  };
  
	pthread_mutex_t mythreadlock;
	
	void remove_hit(double adc, int phibin, int zbin, std::multimap<unsigned short, ihit> &all_hit_map, std::vector<std::vector<unsigned short>> &adcval)
	{
	  typedef std::multimap<unsigned short, ihit>::iterator hit_iterator;
	  std::pair<hit_iterator, hit_iterator> iterpair = all_hit_map.equal_range(adc);
	  hit_iterator it = iterpair.first;
	  for (; it != iterpair.second; ++it) {
	    if (it->second.second.first == phibin && it->second.second.second == zbin) { 
	      all_hit_map.erase(it);
	      break;
	    }
	  }
	  adcval[phibin][zbin] = 0;
	}
	
	void remove_hits(std::vector<ihit> &ihit_list, std::multimap<unsigned short, ihit> &all_hit_map,std::vector<std::vector<unsigned short>> &adcval)
	{
	  for(auto iter = ihit_list.begin(); iter != ihit_list.end();++iter){
	    unsigned short adc    = iter->first; 
	    unsigned short phibin = iter->second.first;
	    unsigned short zbin   = iter->second.second;
	    remove_hit(adc,phibin,zbin,all_hit_map,adcval);
	  }
	}
	
	void get_cluster(int phibin, int zbin, const std::vector<std::vector<unsigned short>> &adcval, std::vector<ihit> &ihit_list)
	{
	  // on hit = one cluster
	  const int& iphi = phibin;
	  const int& iz = zbin;
	  iphiz iCoord(std::make_pair(iphi,iz));
	  ihit  thisHit(adcval[iphi][iz],iCoord);
	  ihit_list.push_back(thisHit);
	}
		
	void calc_cluster_parameter(std::vector<ihit> &ihit_list, thread_data& my_data)
	{
	
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
	
	  std::vector<TrkrDefs::hitkey> hitkeyvec;
	  for(auto iter = ihit_list.begin(); iter != ihit_list.end();++iter){
	    double adc = iter->first; 
	
	    if (adc <= 0) continue;
	
	    int iphi = iter->second.first + my_data.phioffset;
	    int iz   = iter->second.second + my_data.zoffset;
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
	    z_sum += z * adc;
	    z2_sum += square(z)*adc;
	
	    adc_sum += adc;
	
	    // capture the hitkeys for all adc values above a certain threshold
	    TrkrDefs::hitkey hitkey = TpcDefs::genHitKey(iphi, iz);
	    // if(adc>5)
	    hitkeyvec.push_back(hitkey);
	  }
	  // if (adc_sum < 10) return;  // skip obvious noise "clusters"
	  
	  // This is the global position
	  double clusphi = phi_sum / adc_sum;
	  double clusz = z_sum / adc_sum;
	  const double clusx = radius * std::cos(clusphi);
	  const double clusy = radius * std::sin(clusphi);
    
	  const double phi_cov = phi2_sum/adc_sum - square(clusphi);
	  const double z_cov = z2_sum/adc_sum - square(clusz);
	
	  // create the cluster entry directly in the node tree	
	  // Estimate the errors
	  const double phi_err_square = (phibinhi == phibinlo) ?
	    square(radius*my_data.layergeom->get_phistep())/12:
	    square(radius)*phi_cov/(adc_sum*0.14);
	  
	  const double z_err_square = (zbinhi == zbinlo) ?
	    square(my_data.layergeom->get_zstep())/12:
	    z_cov/(adc_sum*0.14);
	
	  // phi_cov = (weighted mean of dphi^2) - (weighted mean of dphi)^2,  which is essentially the weighted mean of dphi^2. The error is then:
	  // e_phi = sigma_dphi/sqrt(N) = sqrt( sigma_dphi^2 / N )  -- where N is the number of samples of the distribution with standard deviation sigma_dphi
	  //    - N is the number of electrons that drift to the readout plane
	  // We have to convert (sum of adc units for all bins in the cluster) to number of ionization electrons N
	  // Conversion gain is 20 mV/fC - relates total charge collected on pad to PEAK voltage out of ADC. The GEM gain is assumed to be 2000
	  // To get equivalent charge per Z bin, so that summing ADC input voltage over all Z bins returns total input charge, divide voltages by 2.4 for 80 ns SAMPA
	  // Equivalent charge per Z bin is then  (ADU x 2200 mV / 1024) / 2.4 x (1/20) fC/mV x (1/1.6e-04) electrons/fC x (1/2000) = ADU x 0.14
		
    // cluster z correction
    clusz -= (clusz<0) ? my_data.par0_neg:my_data.par0_pos;
	

    // create cluster and fill
    auto clus = new TrkrClusterv3;
    clus->setAdc(adc_sum);
	  
	  /// Get the surface key to find the surface from the map
	  TrkrDefs::hitsetkey tpcHitSetKey = TpcDefs::genHitSetKey(my_data.layer, my_data.sector, my_data.side);
	
	  Acts::Vector3 global(clusx, clusy, clusz);
	  
	  TrkrDefs::subsurfkey subsurfkey;
	  Surface surface = my_data.tGeometry->get_tpc_surface_from_coords(
            tpcHitSetKey,
	    global,
	    subsurfkey);
	
	  if(!surface)
	    {
	      /// If the surface can't be found, we can't track with it. So 
	      /// just return and don't add the cluster to the container
	      return;
	    }
	
	  clus->setSubSurfKey(subsurfkey);
	
	  Acts::Vector3 center = surface->center(my_data.tGeometry->geometry().getGeoContext()) 
	    / Acts::UnitConstants::cm;
	  
	  /// no conversion needed, only used in acts
	  const Acts::Vector3 normal = surface->normal(my_data.tGeometry->geometry().getGeoContext());
   const double clusRadius = std::sqrt(square(clusx) + square(clusy));
	  const double rClusPhi = clusRadius * clusphi;
   const double surfRadius = sqrt(center(0)*center(0) + center(1)*center(1));
   const double surfPhiCenter = atan2(center[1], center[0]);
	  const double surfRphiCenter = surfPhiCenter * surfRadius;
   const double surfZCenter = center[2];
   auto local = surface->globalToLocal(my_data.tGeometry->geometry().getGeoContext(),
				       global * Acts::UnitConstants::cm,
				       normal);
	  Acts::Vector2 localPos;
	  
	  /// Prefer Acts transformation since we build the TPC surfaces manually
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
	      
	  clus->setLocalX(localPos(0));
	  clus->setLocalY(localPos(1));
   clus->setActsLocalError(0,0, phi_err_square);
   clus->setActsLocalError(1,0, 0);
   clus->setActsLocalError(0,1, 0);
   clus->setActsLocalError(1,1, z_err_square);

    my_data.cluster_vector.push_back(clus);

	  // Add the hit associations to the TrkrClusterHitAssoc node
	  // we need the cluster key and all associated hit keys (note: the cluster key includes the hitset key)

    if( my_data.do_assoc ) 
    {
      // get cluster index in vector. It is used to store associations, and build relevant cluster keys when filling the containers
      uint32_t index = my_data.cluster_vector.size()-1;
      for (unsigned int i = 0; i < hitkeyvec.size(); i++){
        my_data.association_vector.push_back(std::make_pair(index, hitkeyvec[i]));
	    }
	  }
	}
	
	void *ProcessSector(void *threadarg) {
    
    auto my_data = (struct thread_data *) threadarg;

    const auto& pedestal = my_data->pedestal;
    const auto& phibins   = my_data->phibins;
    const auto& phioffset = my_data->phioffset;
    const auto& zbins     = my_data->zbins ;
    const auto& zoffset   = my_data->zoffset ;

	   TrkrHitSet *hitset = my_data->hitset;
	   TrkrHitSet::ConstRange hitrangei = hitset->getHits();
	
	   // for convenience, create a 2D vector to store adc values in and initialize to zero
	   std::vector<std::vector<unsigned short>> adcval(phibins, std::vector<unsigned short>(zbins, 0));
	   std::multimap<unsigned short, ihit> all_hit_map;
	   std::vector<ihit> hit_vect;
	
	   for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
		hitr != hitrangei.second;
		++hitr){
	     unsigned short phibin = TpcDefs::getPad(hitr->first) - phioffset;
	     unsigned short zbin = TpcDefs::getTBin(hitr->first) - zoffset;
	     
	     float_t fadc = (hitr->second->getAdc()) - pedestal; // proper int rounding +0.5
	     //std::cout << " layer: " << my_data->layer  << " phibin " << phibin << " zbin " << zbin << " fadc " << hitr->second->getAdc() << " pedestal " << pedestal << " fadc " << std::endl
	
	     unsigned short adc = 0;
	     if(fadc>0) 
	       adc =  (unsigned short) fadc;
	     
	//     if(phibin < 0) continue; // phibin is unsigned int, <0 cannot happen
	     if(phibin >= phibins) continue;
	//     if(zbin   < 0) continue;
	     if(zbin   >= zbins) continue; // zbin is unsigned int, <0 cannot happen
	
	     if(adc>0){
	       iphiz iCoord(std::make_pair(phibin,zbin));
	       ihit  thisHit(adc,iCoord);
	       if(adc>5){
		 all_hit_map.insert(std::make_pair(adc, thisHit));
	       }
	       //adcval[phibin][zbin] = (unsigned short) adc;
	       adcval[phibin][zbin] = (unsigned short) adc;
	     }
	   }
	
	   while(all_hit_map.size()>0){
	
	     auto iter = all_hit_map.rbegin();
	     if(iter == all_hit_map.rend()){
	       break;
	     }
	     ihit hiHit = iter->second;
	     int iphi = hiHit.second.first;
	     int iz = hiHit.second.second;
	     
	     //put all hits in the all_hit_map (sorted by adc)
	     //start with highest adc hit
	     // -> cluster around it and get vector of hits
	     std::vector<ihit> ihit_list;
	     get_cluster(iphi, iz, adcval, ihit_list);
	
	     // -> calculate cluster parameters
	     // -> add hits to truth association
	     // remove hits from all_hit_map
	     // repeat untill all_hit_map empty
	     calc_cluster_parameter(ihit_list, *my_data );
	     remove_hits(ihit_list,all_hit_map, adcval);
	   }
	   pthread_exit(nullptr);
	}
}

TpcSimpleClusterizer::TpcSimpleClusterizer(const std::string &name)
  : SubsysReco(name)
{}

bool TpcSimpleClusterizer::is_in_sector_boundary(int phibin, int sector, PHG4TpcCylinderGeom *layergeom) const
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


int TpcSimpleClusterizer::InitRun(PHCompositeNode *topNode)
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

int TpcSimpleClusterizer::process_event(PHCompositeNode *topNode)
{
  //  int print_layer = 18;

  if (Verbosity() > 1000)
    std::cout << "TpcSimpleClusterizer::Process_Event" << std::endl;

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
  
  for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
       hitsetitr != hitsetrange.second;
       ++hitsetitr)
  {
    TrkrHitSet *hitset = hitsetitr->second;
    unsigned int layer = TrkrDefs::getLayer(hitsetitr->first);
    int side = TpcDefs::getSide(hitsetitr->first);
    unsigned int sector= TpcDefs::getSectorId(hitsetitr->first);
    PHG4TpcCylinderGeom *layergeom = geom_container->GetLayerCellGeom(layer);
    
    // instanciate new thread pair, at the end of thread vector
    thread_pair_t& thread_pair = threads.emplace_back();
    
    thread_pair.data.layergeom = layergeom;
    thread_pair.data.hitset = hitset;
    thread_pair.data.layer = layer;
    thread_pair.data.pedestal = pedestal;
    thread_pair.data.sector = sector;
    thread_pair.data.side = side;
    thread_pair.data.do_assoc = do_hit_assoc;
    thread_pair.data.tGeometry = m_tGeometry;
    thread_pair.data.par0_neg = par0_neg;
    thread_pair.data.par0_pos = par0_pos;

    unsigned short NPhiBins = (unsigned short) layergeom->get_phibins();
    unsigned short NPhiBinsSector = NPhiBins/12;
    unsigned short NZBins = (unsigned short)layergeom->get_zbins();
    unsigned short NZBinsSide = NZBins/2;
    unsigned short NZBinsMin = 0;
    unsigned short PhiOffset = NPhiBinsSector * sector;

    if (side == 0){
      NZBinsMin = 0;
    }
    else{
      NZBinsMin = NZBins / 2;
    }

    unsigned short ZOffset = NZBinsMin;

    thread_pair.data.phibins   = NPhiBinsSector;
    thread_pair.data.phioffset = PhiOffset;
    thread_pair.data.zbins     = NZBinsSide;
    thread_pair.data.zoffset   = ZOffset ;

    int rc = pthread_create(&thread_pair.thread, &attr, ProcessSector, (void *)&thread_pair.data);
    if (rc) {
      std::cout << "Error:unable to create thread," << rc << std::endl;
    }
  }
  
  pthread_attr_destroy(&attr);

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

int TpcSimpleClusterizer::End(PHCompositeNode */*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
