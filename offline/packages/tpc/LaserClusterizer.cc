#include "LaserClusterizer.h"

#include <trackbase/LaserCluster.h>
#include <trackbase/LaserClusterv1.h>
#include <trackbase/LaserClusterContainer.h>
#include <trackbase/LaserClusterContainerv1.h>

#include <trackbase/TpcDefs.h>
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

#include <TFile.h>  

#include <cmath>  // for sqrt, cos, sin
#include <iostream>
#include <map>  // for _Rb_tree_cons...
#include <string>
#include <utility>  // for pair
#include <array>
#include <vector>
#include <limits>

#if !defined(__CINT__) || defined(__CLING__)
//BOOST for combi seeding
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>

#include <boost/geometry/index/rtree.hpp>
#endif



namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;


typedef bg::model::point<float, 3, bg::cs::cartesian> point;
typedef bg::model::box<point> box;
typedef std::pair<TrkrDefs::hitkey, TrkrDefs::hitsetkey> specHitKey;
typedef std::pair<point, specHitKey> pointKeyLaser;


LaserClusterizer::LaserClusterizer(const std::string& name)
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
  auto laserclusters = findNode::getClass<LaserClusterContainerv1>(dstNode, "LASER_CLUSTER");
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
      new PHIODataNode<PHObject>(laserclusters, "LASER_CLUSTER", "PHObject");
    DetNode->addNode(LaserClusterContainerNode);
  }
  

  m_debugFile = new TFile(m_debugFileName.c_str(),"RECREATE");
  m_clusterTree = new TTree("clusterTree","clusterTree");
  m_clusterTree->Branch("event",m_event);
  m_clusterTree->Branch("cluster",&m_currentCluster);
  

  return Fun4AllReturnCodes::EVENT_OK;

}

int LaserClusterizer::process_event(PHCompositeNode *topNode)
{
  ++m_event;

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
  m_clusterlist = findNode::getClass<LaserClusterContainerv1>(topNode, "LASER_CLUSTER");
  if (!m_clusterlist)
    {
      std::cout << PHWHERE << " ERROR: Can't find LASER_CLUSTER." << std::endl;
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


  TrkrHitSetContainer::ConstRange hitsetrange;
  RawHitSetContainer::ConstRange rawhitsetrange;
  //int num_hitsets = 0;

  if(!do_read_raw){
    hitsetrange = m_hits->getHitSets(TrkrDefs::TrkrId::tpcId);
    //num_hitsets = std::distance(hitsetrange.first,hitsetrange.second);
  }else{
    rawhitsetrange = m_rawhits->getHitSets(TrkrDefs::TrkrId::tpcId);
    //num_hitsets = std::distance(rawhitsetrange.first,rawhitsetrange.second);
  }

  m_tdriftmax = AdcClockPeriod * NZBinsSide;  

  bgi::rtree<pointKeyLaser, bgi::quadratic<16> > rtree;
  //std::multimap <unsigned int,std::pair<TrkrDefs::hitkey,TrkrDefs::hitsetkey>> adcMap;
  std::multimap <unsigned int,std::pair<std::pair<TrkrDefs::hitkey,TrkrDefs::hitsetkey>, std::array<float, 3>>> adcMap;
  
  //std::multimap <unsigned int, float *> adcCoords;

  if(!do_read_raw){
    for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
	 hitsetitr != hitsetrange.second;
	 ++hitsetitr)
      {
	//	if(count>0)continue;
	TrkrHitSet *hitset = hitsetitr->second;
	unsigned int layer = TrkrDefs::getLayer(hitsetitr->first);
	int side = TpcDefs::getSide(hitsetitr->first);
	unsigned int sector= TpcDefs::getSectorId(hitsetitr->first);
	PHG4TpcCylinderGeom *layergeom = geom_container->GetLayerCellGeom(layer);
	double r = layergeom->get_radius();

	if(sector != 0 || side != 1) continue;

	//unsigned short phiOffset = (((unsigned short) layergeom->get_phibins())/12) * sector;
	//unsigned short tOffset = 0;

	TrkrDefs::hitsetkey hitsetKey = TpcDefs::genHitSetKey( layer, sector, side );   

	
	TrkrHitSet::ConstRange hitrangei = hitset->getHits();
	
	for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
	     hitr != hitrangei.second;
	     ++hitr)
	  {

	    float_t fadc = (hitr->second->getAdc()) - m_pedestal; // proper int rounding +0.5
	    unsigned short adc = 0;
	    if(fadc>0) adc =  (unsigned short) fadc;
	    if(adc<=0) continue;


	    
	    int iphi = TpcDefs::getPad(hitr->first);
	    int it = TpcDefs::getTBin(hitr->first);
	    double phi = layergeom->get_phi(iphi);
	    double zdriftlength = layergeom->get_zcenter(it) * m_tGeometry->get_drift_velocity();
	    
	    double x = r * cos(phi);
	    double y = r * sin(phi);
	    double z = m_tdriftmax * m_tGeometry->get_drift_velocity() - zdriftlength;
	    if (side == 0) z = -z;
	    

	    std::cout << "hit " << std::distance(hitrangei.first, hitr) << "   sector " << sector << "   side " << side << "   coords=(" << r << "," << phi << "," << z << ")" << std::endl;

	    std::array<float, 3> coords = {(float)x,(float)y,(float)z};
	    
	    std::vector<pointKeyLaser> testduplicate;
	    rtree.query(bgi::intersects(box(point(x-0.001,y-0.001,z-0.001),
					    point(x+0.001,y+0.001,z+0.001))),std::back_inserter(testduplicate));
	    if (!testduplicate.empty()){
	      continue;
	    }

	    TrkrDefs::hitkey hitKey = TpcDefs::genHitKey(iphi, it);
	    
	    auto spechitkey = std::make_pair(hitKey,hitsetKey);
	    auto keyCoords = std::make_pair(spechitkey,coords);
	    //adcMap.insert(std::make_pair(adc,spechitkey));
	    adcMap.insert(std::make_pair(adc,keyCoords));
	    //adcCoords.insert(std::make_pair(adc,coords));
	    rtree.insert(std::make_pair(point(x,y,z), spechitkey));
	    
	  }
      }
  }

  std::cout << "finished looping over hits" << std::endl;
  std::cout << "map size: "  << adcMap.size() << std::endl;
  std::cout << "rtree size: " << rtree.size() << std::endl;

  //done filling rTree
  
  while(adcMap.size() > 0){
    auto iterKey = adcMap.rbegin();
    if(iterKey == adcMap.rend()){
      break;
    }
    
    //    auto iterC = adcCoords.rbegin();
    //    if( iterC == adcCoords.rend()){
    //      break;
    //    }


    //auto spechitkey = iterKey->second;
    auto coords = iterKey->second.second;

    double r = sqrt(coords[0]*coords[0] + coords[1]*coords[1]);
    double phi = atan2(coords[1],coords[0]);
    double z = coords[2];

    std::cout << "current number of clusters: " << m_clusterlist->size() << std::endl;
    std::cout << "max adc=" << iterKey->first << "   coordinates (R,phi,z)=(" << r << "," << phi << "," << z << ")" << std::endl;
    std::cout << "coordinates (X,Y,Z)=(" << coords[0] << "," << coords[1] << "," << coords[2] << ")" << std::endl;

    double half_dr = 1.097047535;
    double half_dphi = 0.015370975;
    double half_dz = 1.0;
    if (r<40.9){
      half_dr = 0.56598627;
      half_dphi = 0.029869262;
    }
    else if(r<58.0){
      half_dr = 1.02069265;
      half_dphi = 0.022094093;
    }

    double ppX = (r + half_dr) * cos(phi + half_dphi);
    double pmX = (r + half_dr) * cos(phi - half_dphi);
    double mmX = (r - half_dr) * cos(phi - half_dphi);
    double mpX = (r - half_dr) * cos(phi + half_dphi);

    double lowX = 1e5;
    if(ppX < lowX) lowX = ppX;
    if(pmX < lowX) lowX = pmX;
    if(mmX < lowX) lowX = mmX;
    if(mpX < lowX) lowX = mpX;

    double highX = -1e5;
    if(ppX > highX) highX = ppX;
    if(pmX > highX) highX = pmX;
    if(mmX > highX) highX = mmX;
    if(mpX > highX) highX = mpX;

    double ppY = (r + half_dr) * sin(phi + half_dphi);
    double pmY = (r + half_dr) * sin(phi - half_dphi);
    double mmY = (r - half_dr) * sin(phi - half_dphi);
    double mpY = (r - half_dr) * sin(phi + half_dphi);

    double lowY = 1e5;
    if(ppY < lowY) lowY = ppY;
    if(pmY < lowY) lowY = pmY;
    if(mmY < lowY) lowY = mmY;
    if(mpY < lowY) lowY = mpY;

    double highY = -1e5;
    if(ppY > highY) highY = ppY;
    if(pmY > highY) highY = pmY;
    if(mmY > highY) highY = mmY;
    if(mpY > highY) highY = mpY;

    double lowZ = z - half_dz;
    double highZ = z + half_dz;

    std::cout << "box low (" << lowX << "," << lowY << "," << lowZ << ")" << std::endl;
    std::cout << "box high (" << highX << "," << highY << "," << highZ << ")" << std::endl;

    bool good_bounds = false;
    if(coords[0]>lowX && coords[0]<highX && coords[1]>lowY && coords[1]<highY && coords[2]>lowZ && coords[2]<highZ) good_bounds = true;

    std::cout << "is max adc point within range? : " << good_bounds << std::endl;

    vector<pointKeyLaser> clusHits;

    if(rtree.empty()) std::cout << "empty tree?" << std::endl;

    rtree.query(bgi::intersects(box(point(lowX,lowY,lowZ),point(highX,highY,highZ))),std::back_inserter(clusHits));
    //rtree.query(bgi::contains(box(point(lowX,lowY,lowZ),point(highX,highY,highZ))),std::back_inserter(clusHits));

    std::cout << "number of clusters in box: " << clusHits.size() << std::endl;

    calc_cluster_parameter(clusHits, adcMap);

    //remove_hits(clusHits, rtree, adcMap, adcCoords);
    remove_hits(clusHits, rtree, adcMap);
    clusHits.clear();

  }

  return Fun4AllReturnCodes::EVENT_OK;


}

int LaserClusterizer::End(PHCompositeNode */*topNode*/)
{
  
  m_debugFile->cd();
  m_clusterTree->Write();

  m_debugFile->Close();

  return Fun4AllReturnCodes::EVENT_OK;
}


//void LaserClusterizer::calc_cluster_parameter(vector<pointKeyLaser> &clusHits, std::multimap<unsigned int,std::pair<TrkrDefs::hitkey,TrkrDefs::hitsetkey>> &adcMap)
void LaserClusterizer::calc_cluster_parameter(vector<pointKeyLaser> &clusHits, std::multimap<unsigned int,std::pair<std::pair<TrkrDefs::hitkey,TrkrDefs::hitsetkey>, std::array<float, 3>>> &adcMap)
{

  double xSum = 0.0;
  double ySum = 0.0;
  double zSum = 0.0;
  
  double adcSum = 0.0;

  double maxAdc = 0.0;
  TrkrDefs::hitsetkey maxKey = 0;

  unsigned int nHits = clusHits.size();

  for(auto iter = clusHits.begin(); iter != clusHits.end(); ++iter){
    float coords[3] = {iter->first.get<0>(), iter->first.get<1>(), iter->first.get<2>()};
    std::pair<TrkrDefs::hitkey,TrkrDefs::hitsetkey> spechitkey = iter->second;

    for(auto iterKey = adcMap.begin(); iterKey != adcMap.end(); ++iterKey){
      if(iterKey->second.first == spechitkey){
	
	double adc = iterKey->first;

	xSum += coords[0] * adc;
	ySum += coords[1] * adc;
	zSum += coords[2] * adc;

	adcSum += adc;

	if(adc > maxAdc){
	  maxAdc = adc;
	  maxKey = spechitkey.second;
	}
	
	break;
      }
    }

  }

  if(nHits == 0) return;

  //std::cout << "making cluster " << m_clusterlist->size() << std::endl;
  //std::cout << "adc=" << adcSum << "   X=" << xSum/adcSum << "   Y=" << ySum/adcSum << "   Z=" << zSum/adcSum << "   nHits=" << nHits << std::endl;

  auto *clus = new LaserClusterv1;
  clus->setAdc(adcSum);
  clus->setX(xSum/adcSum);
  clus->setY(ySum/adcSum);
  clus->setZ(zSum/adcSum);
  clus->setNhits(nHits);

  //std::cout << "made cluster" << std::endl;

  const auto ckey = TrkrDefs::genClusKey( maxKey, m_clusterlist->size());
  
  //std::cout << "made key: " << ckey << std::endl;

  m_clusterlist->addClusterSpecifyKey( ckey, clus);

  m_currentCluster = (LaserClusterv1*)clus->CloneMe();

  //  std::cout << "cloned cluster" << std::endl;

  m_clusterTree->Fill();

  //std::cout << "filled tree" << std::endl;

}

//void LaserClusterizer::remove_hits(std::vector<pointKeyLaser> &clusHits,  bgi::rtree<pointKeyLaser, bgi::quadratic<16> > &rtree, std::multimap <unsigned int, std::pair<TrkrDefs::hitkey,TrkrDefs::hitsetkey>> &adcMap, std::multimap <unsigned int, float*> &adcCoords)
void LaserClusterizer::remove_hits(std::vector<pointKeyLaser> &clusHits,  bgi::rtree<pointKeyLaser, bgi::quadratic<16> > &rtree, std::multimap <unsigned int, std::pair<std::pair<TrkrDefs::hitkey,TrkrDefs::hitsetkey>, std::array<float, 3>>> &adcMap)
{

  for(auto iter = clusHits.begin(); iter != clusHits.end(); ++iter){
    auto spechitkey = iter->second;

    //std::cout << "rtree size before removal: " << rtree.size() << std::endl;

    rtree.remove(*iter);

    //std::cout << "rtree size after removal: " << rtree.size() << std::endl;

    //for(std::multimap <unsigned int, std::pair<TrkrDefs::hitkey,TrkrDefs::hitsetkey>>::iterator iterAdc = adcMap.begin(); iterAdc != adcMap.end();){
    for(auto iterAdc = adcMap.begin(); iterAdc != adcMap.end();){
      if(iterAdc->second.first == spechitkey){
	//std::cout << "current adc " << iterAdc->first << "   hitset=" << iterAdc->second.first.second << "   hitkey=" << iterAdc->second.first.first << "   coords=(" <<iterAdc->second.second[0] << "," << iterAdc->second.second[1] << "," << iterAdc->second.second[2] << ")" << std::endl;
	//std::cout << "tree hitset=" << iter->second.second << "   hitkey=" << iter->second.first << "   coords=(" << iter->first.get<0>() << "," << iter->first.get<1>() << "," << iter->first.get<2>() << ")" << std::endl;
	//std::cout << "adcMap size before removal: " << adcMap.size() << std::endl;
	iterAdc = adcMap.erase(iterAdc);
	//std::cout << "adcMap size after removal: " << adcMap.size() << std::endl;
	break;
      }else{
	++iterAdc;
      }
    }
    
    /*
    for(std::multimap <unsigned int, float*>::iterator iterCoords = adcCoords.begin(); iterCoords != adcCoords.end();){
      double mapX = iterCoords->second[0];
      double mapY = iterCoords->second[1];
      double mapZ = iterCoords->second[2];

      double treeX = iter->first.get<0>();
      double treeY = iter->first.get<1>();
      double treeZ = iter->first.get<2>();
      
      std::cout << "adcCoords index: " << std::distance(adcCoords.begin(),iterCoords) << std::endl;
      std::cout << "map coords: (" << mapX << "," << mapY << "," << mapZ << ")" << std::endl;
      std::cout << "tree coords: (" << treeX << "," << treeY << "," << treeZ << ")" << std::endl;


      if(mapX == treeX && mapY == treeY && mapZ == treeZ){

	std::cout << "adcCoords size before removal: " << adcCoords.size() << std::endl;
	iterCoords = adcCoords.erase(iterCoords);
	std::cout << "adcCoords size after removal: " << adcCoords.size() << std::endl;
	break;
      }else{
	std::cout << "no match" << std::endl;
	++iterCoords;
      }
    }
    */
  }
}
