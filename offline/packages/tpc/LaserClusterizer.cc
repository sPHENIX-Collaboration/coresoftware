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
#include <phool/PHTimer.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <TFile.h>  
#include <TF1.h>

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

  m_itHist_0 = new TH1I("m_itHist_0","side 0;it",360,-0.5,359.5);
  m_itHist_1 = new TH1I("m_itHist_1","side 1;it",360,-0.5,359.5);

  m_clusterTree = new TTree("clusterTree","clusterTree");
  m_clusterTree->Branch("event",&m_event);
  //m_clusterTree->Branch("cluster",&m_currentCluster);
  m_clusterTree->Branch("clusters",&m_eventClusters);
  m_clusterTree->Branch("itHist_0",&m_itHist_0);
  m_clusterTree->Branch("itHist_1",&m_itHist_1);
  m_clusterTree->Branch("nClusters",&m_nClus);
  m_clusterTree->Branch("time_search",&time_search);
  m_clusterTree->Branch("time_clus",&time_clus);
  m_clusterTree->Branch("time_erase",&time_erase);
  m_clusterTree->Branch("time_all",&time_all);
  
  //m_clusterTree->Branch("clusters",&m_clusterlist);

  /*
  m_hitTree = new TTree("hitTree","hitTree");
  m_hitTree->Branch("event",&m_event);
  m_hitTree->Branch("hit",&m_currentHit);
  m_hitTree->Branch("hitH",&m_currentHit_hardware);
  */

  m_tdriftmax = AdcClockPeriod * NZBinsSide;  

  t_all = std::make_unique<PHTimer>("t_all");
  t_all->stop();
  t_search = std::make_unique<PHTimer>("t_search");
  t_search->stop();
  t_clus = std::make_unique<PHTimer>("t_clus");
  t_clus->stop();
  t_erase = std::make_unique<PHTimer>("t_erase");
  t_erase->stop();

  return Fun4AllReturnCodes::EVENT_OK;

}

int LaserClusterizer::process_event(PHCompositeNode *topNode)
{
  ++m_event;

  std::cout << "LaserClusterizer::process_event working on event " << m_event << std::endl;

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


  m_geom_container =
    findNode::getClass<PHG4TpcCylinderGeomContainer>(topNode, "CYLINDERCELLGEOM_SVTX");
  if (!m_geom_container)
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


  bgi::rtree<pointKeyLaser, bgi::quadratic<16> > rtree;
  //std::multimap <unsigned int,std::pair<TrkrDefs::hitkey,TrkrDefs::hitsetkey>> adcMap;
  std::multimap <unsigned int,std::pair<std::pair<TrkrDefs::hitkey,TrkrDefs::hitsetkey>, std::array<int, 3>>> adcMap;
  
  //std::multimap <unsigned int, float *> adcCoords;

  if(!do_read_raw){
    for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
	 hitsetitr != hitsetrange.second;
	 ++hitsetitr)
      {
	//	if(count>0)continue;
	TrkrHitSet *hitset = hitsetitr->second;
	int side = TpcDefs::getSide(hitsetitr->first);

	TrkrHitSet::ConstRange hitrangei = hitset->getHits();
	
	for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
	     hitr != hitrangei.second;
	     ++hitr)
	  {

	    int it = TpcDefs::getTBin(hitr->first);

	    if(side == 0)
	      {
		m_itHist_0->Fill(it);
	      }
	    else
	      {
		m_itHist_1->Fill(it);
	      }

	  }
      }


    double itMeanContent_0 = 0.0;
    double itMeanContent_1 = 0.0;

    for(int i=1; i<=360; i++)
      {
	itMeanContent_0 += m_itHist_0->GetBinContent(i);
	itMeanContent_1 += m_itHist_1->GetBinContent(i);
      }

    itMeanContent_0 = itMeanContent_0 / 360.0;
    itMeanContent_1 = itMeanContent_1 / 360.0;


    m_itHist_0->GetXaxis()->SetRange(150, 360);
    double itMax_0 = m_itHist_0->GetBinCenter(m_itHist_0->GetMaximumBin());
    double itMaxContent_0 = m_itHist_0->GetMaximum();
    m_itHist_0->GetXaxis()->SetRange(0, 0);

    m_itHist_1->GetXaxis()->SetRange(150, 360);
    double itMax_1 = m_itHist_1->GetBinCenter(m_itHist_1->GetMaximumBin());
    double itMaxContent_1 = m_itHist_1->GetMaximum();
    m_itHist_1->GetXaxis()->SetRange(0, 0);

    if(itMaxContent_0 / itMeanContent_0 < 1.5 || itMaxContent_1 / itMeanContent_1 < 1.5)
      {
	return Fun4AllReturnCodes::ABORTEVENT;
      }
    

    for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
	 hitsetitr != hitsetrange.second;
	 ++hitsetitr)
      {
	//	if(count>0)continue;
	TrkrHitSet *hitset = hitsetitr->second;
	unsigned int layer = TrkrDefs::getLayer(hitsetitr->first);
	int side = TpcDefs::getSide(hitsetitr->first);
	unsigned int sector= TpcDefs::getSectorId(hitsetitr->first);
	PHG4TpcCylinderGeom *layergeom = m_geom_container->GetLayerCellGeom(layer);
	double r = layergeom->get_radius();

	//if(sector != 5 || side != 1) continue;

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
	    if(fadc>0)
	      {
		adc =  (unsigned short) fadc;
	      }
	    if(adc<=0){
	      continue;
	    }


	    
	    int iphi = TpcDefs::getPad(hitr->first);
	    int it = TpcDefs::getTBin(hitr->first);

	    if(side == 0 && fabs(it - itMax_0) > 50){
	      continue;
	    }
	    if(side == 1 && fabs(it - itMax_1) > 50){
	      continue;
	    }

	    double phi = layergeom->get_phi(iphi);
	    //iphi -= phiOffset;
	    double zdriftlength = layergeom->get_zcenter(it) * m_tGeometry->get_drift_velocity();

	    float x = r * cos(phi);
	    float y = r * sin(phi);
	    float z = m_tdriftmax * m_tGeometry->get_drift_velocity() - zdriftlength;
	    if (side == 1){
	      z = -z;
	      it = -it;
	    }


	    m_currentHit.clear();
	    m_currentHit.push_back(x);
	    m_currentHit.push_back(y);
	    m_currentHit.push_back(z);
	    m_currentHit.push_back(1.0*adc);

	    m_currentHit_hardware.clear();
	    m_currentHit_hardware.push_back(layer);
	    m_currentHit_hardware.push_back(iphi);
	    m_currentHit_hardware.push_back(it);
	    m_currentHit_hardware.push_back(1.0*adc);

	    //std::cout << "hit " << std::distance(hitrangei.first, hitr) << "   sector " << sector << "   side " << side << "   coords=(" << layer << "," << iphi << "," << it << ")" << std::endl;
	    
	    std::array<int, 3> coords = {(int)layer,iphi,it};
	    
	    std::vector<pointKeyLaser> testduplicate;
	    rtree.query(bgi::intersects(box(point(layer-0.001,iphi-0.001,it-0.001),
					    point(layer+0.001,iphi+0.001,it+0.001))),std::back_inserter(testduplicate));
	    if (!testduplicate.empty()){
	      testduplicate.clear();
	      continue;
	    }

	    //std::cout << "current hit (x,y,z)=(" << m_currentHit[0] << "," << m_currentHit[1] << "," << m_currentHit[2] << ")" << std::endl;
	    //m_hitTree->Fill();

	    TrkrDefs::hitkey hitKey = TpcDefs::genHitKey(iphi, it);
	    
	    auto spechitkey = std::make_pair(hitKey,hitsetKey);
	    auto keyCoords = std::make_pair(spechitkey,coords);
	    //adcMap.insert(std::make_pair(adc,spechitkey));
	    adcMap.insert(std::make_pair(adc,keyCoords));

	    //std::cout << "inserted into adcMap" << std::endl;
	    //adcCoords.insert(std::make_pair(adc,coords));
	    rtree.insert(std::make_pair(point(1.0*layer,1.0*iphi,1.0*it), spechitkey));
	    //std::cout << "inserted into rtree" << std::endl;

	  }
      }
  }


  if(Verbosity() > 1)
    {
      std::cout << "finished looping over hits" << std::endl;
      std::cout << "map size: "  << adcMap.size() << std::endl;
      std::cout << "rtree size: " << rtree.size() << std::endl;
    }

  //done filling rTree

  t_all->restart();

  
  while(adcMap.size() > 0){
    auto iterKey = adcMap.rbegin();
    if(iterKey == adcMap.rend()){
      break;
    }
    
    auto coords = iterKey->second.second;

    int layer = coords[0];
    int iphi = coords[1];
    int it = coords[2];

    int layerMax = layer+1;
    if(layer == 22 || layer == 38 || layer == 54)
      {
	layerMax = layer;
      }
    int layerMin = layer-1;
    if(layer == 7 || layer == 23 || layer == 39)
      {
	layerMin = layer;
      }

    //std::cout << "current number of clusters: " << m_clusterlist->size() << std::endl;
    //std::cout << "max adc=" << iterKey->first << "   coordinates (layer,iphi,it)=(" << layer << "," << iphi << "," << it << ")" << std::endl;


    vector<pointKeyLaser> clusHits;


    //if(rtree.empty()) std::cout << "empty tree?" << std::endl;
    t_search->restart();
    rtree.query(bgi::intersects(box(point(layerMin,iphi-2,it-5),point(layerMax,iphi+2,it+5))),std::back_inserter(clusHits));
    t_search->stop();
    //rtree.query(bgi::contains(box(point(lowX,lowY,lowZ),point(highX,highY,highZ))),std::back_inserter(clusHits));

    //std::cout << "number of hits in box: " << clusHits.size() << std::endl;

    t_clus->restart();
    calc_cluster_parameter(clusHits, adcMap);
    t_clus->stop();

    //remove_hits(clusHits, rtree, adcMap, adcCoords);
    t_erase->restart();
    remove_hits(clusHits, rtree, adcMap);
    t_erase->stop();


    clusHits.clear();

  }

  m_nClus = (int)m_eventClusters.size();
  t_all->stop();

  time_search = t_search->get_accumulated_time() / 1000.;
  time_clus = t_clus->get_accumulated_time() / 1000.;
  time_erase = t_erase->get_accumulated_time() / 1000.;
  time_all = t_all->get_accumulated_time() / 1000.;

  //std::cout << "filling tree" << std::endl;
  m_clusterTree->Fill();
  //std::cout << "tree filled" << std::endl;



  if(Verbosity())
    {
      std::cout << "rtree search time: " << t_search->get_accumulated_time() / 1000. << " sec" << std::endl;
      std::cout << "clustering time: " << t_clus->get_accumulated_time() / 1000. << " sec" << std::endl;
      std::cout << "erasing time: " << t_erase->get_accumulated_time() / 1000. << " sec" << std::endl;
      std::cout << "total time: " << t_all->get_accumulated_time() / 1000. << " sec" << std::endl;
    }

  return Fun4AllReturnCodes::EVENT_OK;


}

int LaserClusterizer::ResetEvent(PHCompositeNode */*topNode*/)
{
  
  m_itHist_0->Reset();
  m_itHist_1->Reset();

  m_eventClusters.clear();
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int LaserClusterizer::End(PHCompositeNode */*topNode*/)
{
  
  m_debugFile->cd();
  //  m_itHist_0->Write();
  //  m_itHist_1->Write();
  m_clusterTree->Write();
  //m_hitTree->Write();

  m_debugFile->Close();

  return Fun4AllReturnCodes::EVENT_OK;
}



//void LaserClusterizer::calc_cluster_parameter(vector<pointKeyLaser> &clusHits, std::multimap<unsigned int,std::pair<TrkrDefs::hitkey,TrkrDefs::hitsetkey>> &adcMap)
void LaserClusterizer::calc_cluster_parameter(vector<pointKeyLaser> &clusHits, std::multimap<unsigned int,std::pair<std::pair<TrkrDefs::hitkey,TrkrDefs::hitsetkey>, std::array<int, 3>>> &adcMap)
{

  double rSum = 0.0;
  double phiSum = 0.0;
  double tSum = 0.0;

  double layerSum = 0.0;
  double iphiSum = 0.0;
  double itSum = 0.0;
  
  double adcSum = 0.0;

  //int maxSide = 0;

  double maxAdc = 0.0;
  TrkrDefs::hitsetkey maxKey = 0;

  unsigned int nHits = clusHits.size();

  auto *clus = new LaserClusterv1;


  for(auto iter = clusHits.begin(); iter != clusHits.end(); ++iter){

    float coords[3] = {iter->first.get<0>(), iter->first.get<1>(), iter->first.get<2>()};
    std::pair<TrkrDefs::hitkey,TrkrDefs::hitsetkey> spechitkey = iter->second;

    //int side = TpcDefs::getSide(spechitkey.second);
    //unsigned int sector= TpcDefs::getSectorId(spechitkey.second);

    PHG4TpcCylinderGeom *layergeom = m_geom_container->GetLayerCellGeom((int)coords[0]);

    //unsigned short phiOffset = (((unsigned short) layergeom->get_phibins())/12) * sector;

    double r = layergeom->get_radius();
    double phi = layergeom->get_phi(coords[1]);
    double t = layergeom->get_zcenter(fabs(coords[2]));

    double hitzdriftlength = t * m_tGeometry->get_drift_velocity();
    double hitZ = m_tdriftmax * m_tGeometry->get_drift_velocity() - hitzdriftlength;

    for(auto iterKey = adcMap.begin(); iterKey != adcMap.end(); ++iterKey){
      if(iterKey->second.first == spechitkey){

	
	double adc = iterKey->first;

	clus->addHit();
	clus->setHitLayer(clus->getNhits() - 1, coords[0]);
	clus->setHitIPhi(clus->getNhits() - 1, coords[1]);
	clus->setHitIT(clus->getNhits() - 1, coords[2]);
	clus->setHitX(clus->getNhits() - 1, r*cos(phi));
	clus->setHitY(clus->getNhits() - 1, r*sin(phi));
	clus->setHitZ(clus->getNhits() - 1, hitZ);
	clus->setHitAdc(clus->getNhits() - 1, (float)adc);


	rSum += r * adc;
	phiSum += phi * adc;
	tSum += t * adc;

	layerSum += coords[0] * adc;
	iphiSum += coords[1] * adc;
	itSum += coords[2] * adc;

	adcSum += adc;

	if(adc > maxAdc){
	  maxAdc = adc;
	  maxKey = spechitkey.second;
	  //maxSide = side;
	}
	
	break;
      }
    }

  }

  if(nHits == 0)
    {
      return;
    }

  // std::cout << "making cluster " << m_clusterlist->size() << std::endl;
  //std::cout << "adc=" << adcSum << "   r=" << rSum/adcSum << "   phi=" << phiSum/adcSum << "   t=" << tSum/adcSum << "   nHits=" << nHits << std::endl;

  double clusR = rSum / adcSum;
  double clusPhi = phiSum / adcSum;
  double clusT = tSum / adcSum;
  double zdriftlength = clusT * m_tGeometry->get_drift_velocity();
  
  double clusX = clusR * cos(clusPhi);
  double clusY = clusR * sin(clusPhi);
  double clusZ = m_tdriftmax * m_tGeometry->get_drift_velocity() - zdriftlength;
  if(itSum<0){
    clusZ = -clusZ;
    for(int i=0; i<(int)clus->getNhits(); i++){
      clus->setHitZ(i, -1*clus->getHitZ(i));
    }
  }


  clus->setAdc(adcSum);
  clus->setX(clusX);
  clus->setY(clusY);
  clus->setZ(clusZ);
  clus->setLayer(layerSum/adcSum);
  clus->setIPhi(iphiSum/adcSum);
  clus->setIT(itSum/adcSum);

  //std::cout << "made cluster" << std::endl;

  const auto ckey = TrkrDefs::genClusKey( maxKey, m_clusterlist->size());
  
  //  std::cout << "made key: " << ckey << std::endl;

  m_clusterlist->addClusterSpecifyKey( ckey, clus);

  m_currentCluster = (LaserClusterv1*)clus->CloneMe();

  //std::cout << "cloned cluster" << std::endl;

  m_eventClusters.push_back((LaserClusterv1*)m_currentCluster->CloneMe());

  //  std::cout << "cloned cluster into m_eventClusters" << std::endl;

  //m_clusterTree->Fill();

  //std::cout << "filled tree" << std::endl;

}

//void LaserClusterizer::remove_hits(std::vector<pointKeyLaser> &clusHits,  bgi::rtree<pointKeyLaser, bgi::quadratic<16> > &rtree, std::multimap <unsigned int, std::pair<TrkrDefs::hitkey,TrkrDefs::hitsetkey>> &adcMap, std::multimap <unsigned int, float*> &adcCoords)
void LaserClusterizer::remove_hits(std::vector<pointKeyLaser> &clusHits,  bgi::rtree<pointKeyLaser, bgi::quadratic<16> > &rtree, std::multimap <unsigned int, std::pair<std::pair<TrkrDefs::hitkey,TrkrDefs::hitsetkey>, std::array<int, 3>>> &adcMap)
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
