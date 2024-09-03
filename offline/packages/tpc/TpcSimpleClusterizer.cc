#include "TpcSimpleClusterizer.h"

#include <trackbase/TpcDefs.h>

#include <trackbase/TrkrClusterContainerv4.h>
#include <trackbase/TrkrClusterv5.h>
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
  using assoc = std::pair<uint32_t, TrkrDefs::hitkey> ;

  struct ihit
  {
    unsigned short iphi = 0;
    unsigned short it = 0;
    unsigned short adc = 0;
  };

  struct thread_data
  {
    PHG4TpcCylinderGeom *layergeom = nullptr;
    TrkrHitSet *hitset = nullptr;
    ActsGeometry *tGeometry = nullptr;
    unsigned int layer = 0;
    int side = 0;
    unsigned int sector = 0;
    float drift_velocity = 0;

    float pedestal = 0;
    bool do_assoc = true;

    unsigned short phibins = 0;
    unsigned short phioffset = 0;

    unsigned short tbins = 0;
    unsigned short toffset = 0;
    double m_tdriftmax = 0;
    double sampa_tbias = 0;

    std::vector<assoc> association_vector;
    std::vector<TrkrCluster*> cluster_vector;
  };

	pthread_mutex_t mythreadlock;

  void remove_hit(double adc, int phibin, int tbin, std::multimap<unsigned short, ihit> &all_hit_map, std::vector<std::vector<unsigned short>> &adcval)
  {
    const auto iterpair = all_hit_map.equal_range(adc);
    for ( auto it = iterpair.first; it != iterpair.second; ++it)
    {
      if (it->second.iphi == phibin && it->second.it == tbin)
      {
        all_hit_map.erase(it);
        break;
      }
    }
    adcval[phibin][tbin] = 0;
  }

  void remove_hits( const std::vector<ihit> &ihit_list, std::multimap<unsigned short, ihit> &all_hit_map,std::vector<std::vector<unsigned short>> &adcval)
  {
    for( const auto& ihit:ihit_list )
    {
      const unsigned short& adc    = ihit.adc;
      const unsigned short& phibin = ihit.iphi;
      const unsigned short& tbin   = ihit.it;
      remove_hit(adc,phibin,tbin,all_hit_map,adcval);
    }
  }

	void get_cluster(int iphi, int it, const std::vector<std::vector<unsigned short>> &adcval, std::vector<ihit> &ihit_list)
	{
	  // on hit = one cluster
    ihit hit;
    hit.iphi = iphi;
    hit.it = it;
    hit.adc = adcval[iphi][it];
	  ihit_list.push_back(hit);
	}

	void calc_cluster_parameter(std::vector<ihit> &ihit_list, thread_data& my_data)
  {
    double t_sum = 0.0;
    double adc_sum = 0.0;
    double t2_sum = 0.0;

    double iphi_sum = 0.0;
    double iphi2_sum = 0.0;

    double radius = my_data.layergeom->get_radius();  // returns center of layer

    int phibinhi = -1;
    int phibinlo = 666666;
    int tbinhi = -1;
    int tbinlo = 666666;
    int max_adc  = 0;

    std::vector<TrkrDefs::hitkey> hitkeyvec;

    // keep track of the hit locations in a given cluster
    std::map<int,unsigned int> m_phi {};
    std::map<int,unsigned int> m_z   {};

    for(auto iter = ihit_list.begin(); iter != ihit_list.end();++iter)
    {

      double adc = iter->adc;

      if (adc <= 0) continue;
      if(adc > max_adc) max_adc = adc;
      int iphi = iter->iphi + my_data.phioffset;
      int it   = iter->it + my_data.toffset;
      if(iphi > phibinhi) phibinhi = iphi;
      if(iphi < phibinlo) phibinlo = iphi;
      if(it > tbinhi) tbinhi = it;
      if(it < tbinlo) tbinlo = it;

      // update phi sums
      iphi_sum += iphi * adc;
      iphi2_sum += square(iphi)*adc;

      // update t sums
      double t = my_data.layergeom->get_zcenter(it);
      t_sum += t*adc;
      t2_sum += square(t)*adc;

      adc_sum += adc;

      // capture the hitkeys for all adc values above a certain threshold
      TrkrDefs::hitkey hitkey = TpcDefs::genHitKey(iphi, it);
      hitkeyvec.push_back(hitkey);
    }

    // This is the global position
    double clusiphi = iphi_sum / adc_sum;
    double clusphi = my_data.layergeom->get_phi(clusiphi);

    float clusx = radius * cos(clusphi);
    float clusy = radius * sin(clusphi);
    double clust = t_sum / adc_sum;

    // convert z drift length to z position in the TPC
    // needed for surface identification
    double zdriftlength = clust * my_data.tGeometry->get_drift_velocity();
    double clusz  =  my_data.m_tdriftmax * my_data.tGeometry->get_drift_velocity() - zdriftlength;
    if(my_data.side == 0) clusz = -clusz;

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

    auto local = surface->transform(my_data.tGeometry->geometry().getGeoContext()).inverse() * global;
    local /= Acts::UnitConstants::cm;

    {
      auto clus = new TrkrClusterv5;
      clus->setAdc(adc_sum);
      clus->setMaxAdc(max_adc);
      clus->setPhiSize(phisize);
      clus->setZSize(tsize);
      clus->setSubSurfKey(subsurfkey);
      clus->setLocalX(local(0));
      clus->setLocalY(clust);
      clus->setPhiError(sqrt(phi_err_square));
      clus->setZError(sqrt(t_err_square * pow(my_data.tGeometry->get_drift_velocity(),2)));
      my_data.cluster_vector.push_back(clus);
    }

    if(my_data.do_assoc)
    {
      // get cluster index in vector. It is used to store associations, and build relevant cluster keys when filling the containers
      uint32_t index = my_data.cluster_vector.size()-1;
      for (unsigned int i = 0; i < hitkeyvec.size(); i++)
      {
        my_data.association_vector.emplace_back(index, hitkeyvec[i]);
      }

    }

    hitkeyvec.clear();

  }

	void *ProcessSector(void *threadarg)
  {

    auto my_data = (struct thread_data *) threadarg;

    const auto& pedestal  = my_data->pedestal;
    const auto& phibins   = my_data->phibins;
    const auto& phioffset = my_data->phioffset;
    const auto& tbins     = my_data->tbins ;
    const auto& toffset   = my_data->toffset ;

    std::vector<std::vector<unsigned short>> adcval(phibins, std::vector<unsigned short>(tbins, 0));
    std::multimap<unsigned short, ihit> all_hit_map;
    std::vector<ihit> hit_vect;

    int tbinmax = 498;
    int tbinmin = 0;

    if( my_data->hitset!=nullptr)
    {
      auto& hitset = my_data->hitset;
      const auto& hitrangei = hitset->getHits();

      for (TrkrHitSet::ConstIterator hitr = hitrangei.first; hitr != hitrangei.second; ++hitr)
      {

        if( TpcDefs::getPad(hitr->first) - phioffset < 0 ) continue;

        const unsigned short phibin = TpcDefs::getPad(hitr->first) - phioffset;
        if(phibin>=phibins) continue;

        const unsigned short tbin = TpcDefs::getTBin(hitr->first) - toffset;
        if(tbin>=tbins) continue;

        const unsigned short tbinorg = TpcDefs::getTBin(hitr->first);
        if(tbinorg>tbinmax||tbinorg<tbinmin) continue;

        float_t fadc = (hitr->second->getAdc()) - pedestal; // proper int rounding +0.5
        unsigned short adc = 0;
        if(fadc>0) adc =  (unsigned short) fadc;

        if(adc>0)
        {
          ihit  thisHit;
          thisHit.iphi = phibin;
          thisHit.it = tbin;
          thisHit.adc = adc;
          all_hit_map.insert(std::make_pair(adc, thisHit));

          adcval[phibin][tbin] = (unsigned short) adc;

        }
      }

    }

    while(!all_hit_map.empty())
    {

      auto iter = all_hit_map.rbegin();

      const ihit hiHit = iter->second;
      const int iphi = hiHit.iphi;
      const int it = hiHit.it;

      //put all hits in the all_hit_map (sorted by adc)
      //start with highest adc hit
      // -> cluster around it and get vector of hits
      std::vector<ihit> ihit_list;
      get_cluster(iphi, it, adcval, ihit_list );

      // -> calculate cluster parameters
      // -> add hits to truth association
      // remove hits from all_hit_map
      // repeat untill all_hit_map empty
      calc_cluster_parameter( ihit_list, *my_data );
      remove_hits(ihit_list,all_hit_map, adcval);
      ihit_list.clear();
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
    std::cout << std::endl << " mutex init failed" << std::endl;
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

    thread_pair.data.drift_velocity = m_tGeometry->get_drift_velocity();

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
