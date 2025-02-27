#include "LaserClusterizer.h"

#include "LaserEventInfo.h"

#include <trackbase/LaserCluster.h>
#include <trackbase/LaserClusterContainer.h>
#include <trackbase/LaserClusterContainerv1.h>
#include <trackbase/LaserClusterv1.h>
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

#include <pthread.h>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

using point = bg::model::point<float, 3, bg::cs::cartesian>;
using box = bg::model::box<point>;
using specHitKey = std::pair<TrkrDefs::hitkey, TrkrDefs::hitsetkey>;
using pointKeyLaser = std::pair<point, specHitKey>;

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
    std::vector<LaserCluster *> cluster_vector;
    std::vector<TrkrDefs::cluskey> cluster_key_vector;
    double adc_threshold = 74.4;
    int peakTimeBin = 325;
    int layerMin = 1;
    int layerMax = 1;
    double tdriftmax = 0;
  };

  pthread_mutex_t mythreadlock;
  
  void remove_hits(std::vector<pointKeyLaser> &clusHits, bgi::rtree<pointKeyLaser, bgi::quadratic<16>> &rtree, std::multimap<unsigned int, std::pair<std::pair<TrkrDefs::hitkey, TrkrDefs::hitsetkey>, std::array<int, 3>>> &adcMap)
  {
    for (auto &clusHit : clusHits)
    {
      auto spechitkey = clusHit.second;

      rtree.remove(clusHit);

      for (auto iterAdc = adcMap.begin(); iterAdc != adcMap.end();)
      {
	if(iterAdc->second.first == spechitkey)
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

  void calc_cluster_parameter(std::vector<LaserCluster *> &laserClusters, std::vector<TrkrDefs::cluskey> &laserClusterKeys, std::vector<pointKeyLaser> &clusHits, std::multimap<unsigned int, std::pair<std::pair<TrkrDefs::hitkey, TrkrDefs::hitsetkey>, std::array<int, 3>>> &adcMap, thread_data &my_data)
  {
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
    
    auto *clus = new LaserClusterv1;
    
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
      std::pair<TrkrDefs::hitkey, TrkrDefs::hitsetkey> spechitkey = clusHit.second;
      
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
      double phi = layergeom->get_phi(coords[1]);
      double t = layergeom->get_zcenter(fabs(coords[2]));
      
      double hitzdriftlength = t * my_data.tGeometry->get_drift_velocity();
      double hitZ = my_data.tdriftmax * my_data.tGeometry->get_drift_velocity() - hitzdriftlength;
      
      for (auto &iterKey : adcMap)
      {
	if (iterKey.second.first == spechitkey)
	{
	  double adc = iterKey.first;
	  
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
	  
	  break;
	}
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
    
    meanLayer = meanLayer / nHits;
    meanIPhi = meanIPhi / nHits;
    meanIT = meanIT / nHits;
    
    double sigmaLayer = 0.0;
    double sigmaIPhi = 0.0;
    double sigmaIT = 0.0;
    
    double sigmaWeightedLayer = 0.0;
    double sigmaWeightedIPhi = 0.0;
    double sigmaWeightedIT = 0.0;
    
    for (int i = 0; i < (int) clus->getNhits(); i++)
    {
      sigmaLayer += pow(clus->getHitLayer(i) - meanLayer, 2);
      sigmaIPhi += pow(clus->getHitIPhi(i) - meanIPhi, 2);
      sigmaIT += pow(clus->getHitIT(i) - meanIT, 2);
      
      sigmaWeightedLayer += clus->getHitAdc(i) * pow(clus->getHitLayer(i) - (layerSum / adcSum), 2);
      sigmaWeightedIPhi += clus->getHitAdc(i) * pow(clus->getHitIPhi(i) - (iphiSum / adcSum), 2);
      sigmaWeightedIT += clus->getHitAdc(i) * pow(clus->getHitIT(i) - (itSum / adcSum), 2);
    }
    
    clus->setAdc(adcSum);
    clus->setX(clusX);
    clus->setY(clusY);
    clus->setZ(clusZ);
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

    const auto ckey = TrkrDefs::genClusKey(maxKey, laserClusters.size());
    laserClusters.push_back(clus);
    laserClusterKeys.push_back(ckey);
    
  }
  

  void ProcessModuleData(thread_data *my_data)
  {
    
    bgi::rtree<pointKeyLaser, bgi::quadratic<16>> rtree;
    std::multimap<unsigned int, std::pair<std::pair<TrkrDefs::hitkey, TrkrDefs::hitsetkey>, std::array<int, 3>>> adcMap;

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

	if(fabs(it - my_data->peakTimeBin) > 3)
	{
	  continue;
	}

	std::array<int, 3> coords = {(int) layer, iphi, it};

	std::vector<pointKeyLaser> testduplicate;
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
	auto keyCoords = std::make_pair(spechitkey, coords);
	adcMap.insert(std::make_pair(adc, keyCoords));
	rtree.insert(std::make_pair(point(1.0*layer, 1.0*iphi, 1.0*it), spechitkey));

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
      
      auto coords = iterKey->second.second;
      int layer = coords[0];
      int iphi = coords[1];
      int it = coords[2];

      std::vector<pointKeyLaser> clusHits;
      
      rtree.query(bgi::intersects(box(point(layer - my_data->layerMin, iphi - 2, it - 5), point(layer + my_data->layerMax, iphi + 2, it + 5))), std::back_inserter(clusHits));

      calc_cluster_parameter(my_data->cluster_vector, my_data->cluster_key_vector, clusHits, adcMap, *my_data);

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

  if(!m_laserEventInfo->isLaserEvent())
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
	thread_pair.data.cluster_vector = cluster_vector;
	thread_pair.data.cluster_key_vector = cluster_key_vector;
	thread_pair.data.adc_threshold = m_adc_threshold;
	thread_pair.data.peakTimeBin = m_laserEventInfo->getPeakSample(s);
	thread_pair.data.layerMin = (m_lamination ? 2 : 1);
	thread_pair.data.layerMax = (m_lamination ? 2 : 1);
	thread_pair.data.tdriftmax = m_tdriftmax;

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
      
      const auto &data(thread_pair.data);
      
      for(int index = 0; index < (int) data.cluster_vector.size(); ++index)
      {
	auto cluster = data.cluster_vector[index];
	const auto ckey = data.cluster_key_vector[index];
	
	m_clusterlist->addClusterSpecifyKey(ckey, cluster);
      }
    }
  }

  if (Verbosity() > 1)
  {
    std::cout << "LaserClusterizer::process_event " << m_clusterlist->size() << " clusters found" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;

}
