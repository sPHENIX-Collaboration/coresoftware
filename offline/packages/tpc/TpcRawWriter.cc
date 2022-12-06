#include "TpcRawWriter.h"

#include <trackbase/TpcDefs.h>
#include <trackbase/MvtxDefs.h>
#include <trackbase/InttDefs.h>
#include <micromegas/MicromegasDefs.h>

#include <trackbase/TrkrClusterContainerv4.h>
#include <trackbase/TrkrClusterv3.h>
#include <trackbase/TrkrClusterv4.h>
#include <trackbase/TrkrClusterHitAssocv3.h>
#include <trackbase/TrkrDefs.h>  // for hitkey, getLayer
#include <trackbase/TrkrHit.h>
#include <trackbase/RawHitv1.h>
#include <trackbase/RawHitTpc.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/RawHitSet.h>
#include <trackbase/RawHitSetv1.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/RawHitSetContainer.h>
#include <trackbase/RawHitSetContainerv1.h>
#include <trackbase/ActsGeometry.h>
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

  

TpcRawWriter::TpcRawWriter(const std::string &name)
  : SubsysReco(name)
{ std::cout << PHWHERE << "Construct TpcRawWriter" << std::endl; }


int TpcRawWriter::InitRun(PHCompositeNode *topNode)
{
  if(topNode)
    std::cout << PHWHERE << "Init TpcRawWriter" << std::endl;
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

  // new containers
  m_rawhits = findNode::getClass<RawHitSetContainerv1>(topNode, "TRKR_RAWHITSET");
  if (!m_rawhits)
  {
    PHNodeIterator dstiter(dstNode);
    auto DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!DetNode)
    {
      DetNode = new PHCompositeNode("TRKR");
      dstNode->addNode(DetNode);
    }

    m_rawhits = new RawHitSetContainerv1;
    auto newNode = new PHIODataNode<PHObject>(m_rawhits, "TRKR_RAWHITSET", "PHObject");
    DetNode->addNode(newNode);
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcRawWriter::process_event(PHCompositeNode *topNode)
{
  //  int print_layer = 18;

  //  if (Verbosity() > 1000)
  if(topNode)
    std::cout << "TpcRawWriter::Process_Event" << std::endl;

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

  // new containers
  m_rawhits = findNode::getClass<RawHitSetContainerv1>(topNode, "TRKR_RAWHITSET");
  if (!m_rawhits)
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
		<< "ActsTrackingGeometry not found on node tree. Exiting"
		<< std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }


  // The hits are stored in hitsets, where each hitset contains all hits in a given TPC readout (layer, sector, side), so clusters are confined to a hitset
  // The TPC clustering is more complicated than for the silicon, because we have to deal with overlapping clusters
  
  // loop over the MVTX HitSet objects
  std::cout << "processing mvtx" << std::endl;
  TrkrHitSetContainer::ConstRange mvtxhitsetrange = m_hits->getHitSets(TrkrDefs::TrkrId::mvtxId);
  //  const int num_hitsets = std::distance(hitsetrange.first,hitsetrange.second);

  int count = 0;
  for (TrkrHitSetContainer::ConstIterator hitsetitr = mvtxhitsetrange.first;
       hitsetitr != mvtxhitsetrange.second;
       ++hitsetitr)
  {
    TrkrHitSet *hitset = hitsetitr->second;

    m_rawhits->findOrAddHitSet(hitsetitr->first);
    RawHitSetv1 *rhitset = dynamic_cast<RawHitSetv1 *>(m_rawhits->findHitSet(hitsetitr->first));
    TrkrHitSet::ConstRange hitrangei = hitset->getHits();

    for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
	 hitr != hitrangei.second;
	 ++hitr){
      unsigned short adc = (unsigned short)(hitr->second->getAdc());
      if(adc>0){
	RawHitv1  *rHit = new RawHitv1;
	unsigned short iphi =  MvtxDefs::getCol(hitr->first);
	unsigned short it = MvtxDefs::getRow(hitr->first);
	rHit->setPhiBin(iphi);
	rHit->setTBin(it);
	rHit->setAdc(adc);
	rhitset->addHit(rHit);
      }
    }

    count++;
  }
  std::cout << "processing intt" << std::endl;
  // loop over the INTT HitSet objects
  TrkrHitSetContainer::ConstRange intt_hitsetrange = m_hits->getHitSets(TrkrDefs::TrkrId::inttId);
  //  const int num_hitsets = std::distance(intt_hitsetrange.first,intt_hitsetrange.second);

  for (TrkrHitSetContainer::ConstIterator hitsetitr = intt_hitsetrange.first;
       hitsetitr != intt_hitsetrange.second;
       ++hitsetitr)
  {
    TrkrHitSet *hitset = hitsetitr->second;

    m_rawhits->findOrAddHitSet(hitsetitr->first);
    RawHitSet *rhitset = m_rawhits->findHitSet(hitsetitr->first);
    TrkrHitSet::ConstRange hitrangei = hitset->getHits();

    for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
	 hitr != hitrangei.second;
	 ++hitr){
      //int sector = InttDefs::getLadderPhiId(hitsetitr->first);
      //int side = InttDefs::getLadderZId(hitsetitr->first);
      //int layer  = TrkrDefs::getLayer(hitsetitr->first);
      unsigned short adc = (unsigned short)(hitr->second->getAdc());
      if(adc>0){
	RawHitv1  *rHit = new RawHitv1;
	unsigned short iphi =  InttDefs::getCol(hitr->first);
	unsigned short it = InttDefs::getRow(hitr->first);
	//std::cout << " layer " << layer << " sector: " << sector << " side " << side << " col: " << iphi << " row " << it << std::endl;
	rHit->setPhiBin(iphi);
	rHit->setTBin(it);
	rHit->setAdc(adc);
	rhitset->addHit(rHit);
      }
    }

    count++;
  }
  std::cout << "processing mm" << std::endl;
  // loop over the micromega HitSet objects
  TrkrHitSetContainer::ConstRange mm_hitsetrange = m_hits->getHitSets(TrkrDefs::TrkrId::micromegasId);
  //  const int num_hitsets = std::distance(mm_hitsetrange.first,mm_hitsetrange.second);

  for (TrkrHitSetContainer::ConstIterator hitsetitr = mm_hitsetrange.first;
       hitsetitr != mm_hitsetrange.second;
       ++hitsetitr)
  {
    TrkrHitSet *hitset = hitsetitr->second;

    m_rawhits->findOrAddHitSet(hitsetitr->first);
    RawHitSet *rhitset = m_rawhits->findHitSet(hitsetitr->first);
    TrkrHitSet::ConstRange hitrangei = hitset->getHits();

    for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
	 hitr != hitrangei.second;
	 ++hitr){

      unsigned short adc = (unsigned short)(hitr->second->getAdc());
      if(adc>0){
	RawHitv1  *rHit = new RawHitv1;
	//	TrkrDefs::hitkey tmp = (hitr->first >> MicromegasDefs::kBitShiftStrip);
	unsigned short iphi = MicromegasDefs::getStrip( hitr->first );
	unsigned short it = 0;
	rHit->setPhiBin(iphi);
	rHit->setTBin(it);
	rHit->setAdc(adc);
	rhitset->addHit(rHit);
      }
    }

    count++;
  }
  std::cout << "processing tpc" << std::endl;
  // loop over the TPC HitSet objects
  TrkrHitSetContainer::ConstRange tpc_hitsetrange = m_hits->getHitSets(TrkrDefs::TrkrId::tpcId);
  //  const int num_hitsets = std::distance(hitsetrange.first,hitsetrange.second);

  int ncheck = 0;
  int allbad = 0;

  for (TrkrHitSetContainer::ConstIterator hitsetitr = tpc_hitsetrange.first;
       hitsetitr != tpc_hitsetrange.second;
       ++hitsetitr)
  //TrkrHitSetContainer::ConstIterator hitsetitr = tpc_hitsetrange.first;
  {
    TrkrHitSet *hitset = hitsetitr->second;
    unsigned int layer = TrkrDefs::getLayer(hitsetitr->first);
    //   if(layer!=7) continue;
    //    int side = TpcDefs::getSide(hitsetitr->first);
    unsigned int sector= TpcDefs::getSectorId(hitsetitr->first);
    PHG4TpcCylinderGeom *layergeom = geom_container->GetLayerCellGeom(layer);
    unsigned short NPhiBins = (unsigned short) layergeom->get_phibins();
    unsigned short NPhiBinsSector = NPhiBins/12;
    unsigned short NZBins = (unsigned short)layergeom->get_zbins();
    //unsigned short NZBinsSide = NZBins/2;
    unsigned short NZBinsSide = NZBins;
    unsigned short NZBinsMin = 0;
    unsigned short PhiOffset = NPhiBinsSector * sector;
    unsigned short ZOffset = NZBinsMin;
    std::vector<int> nhits_thispad;
    std::vector<std::vector<unsigned short>> adcval;
    adcval.resize(NPhiBins);
    for(int nphi= 0; nphi < NPhiBins;nphi++){
      for(int nz= 0; nz < NZBins;nz++){
	adcval[nphi].push_back(0);
	nhits_thispad.push_back(0);
      }
    }
    //std::vector<std::vector<unsigned short>> outadc;
    //outadc.resize(NPhiBins);

    //    RawHitSetContainerv1::Iterator rawhitset_iter = 
    m_rawhits->findOrAddHitSet(hitsetitr->first);
    RawHitSetv1 *rhitset = dynamic_cast<RawHitSetv1 *>(m_rawhits->findHitSet(hitsetitr->first));
    rhitset->setTpcPhiBins(NPhiBins);

    TrkrHitSet::ConstRange hitrangei = hitset->getHits();
    int zbinmax = 498;
    int zbinmin = 0;
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
    for (TrkrHitSet::ConstIterator hitr = hitrangei.first;
	 hitr != hitrangei.second;
	 ++hitr){
      if( TpcDefs::getPad(hitr->first) - PhiOffset < 0 ){
	continue;
      }
      if( TpcDefs::getTBin(hitr->first) - ZOffset < 0 ){
	//std::cout << "WARNING zbin out of range: " << TpcDefs::getTBin(hitr->first) - zoffset  << " | " << zbins <<std::endl;
	continue;
      }
      unsigned short phibin = TpcDefs::getPad(hitr->first) - PhiOffset;
      unsigned short zbin = TpcDefs::getTBin(hitr->first) - ZOffset;
      unsigned short zbinorg = TpcDefs::getTBin(hitr->first);
      if(phibin>=NPhiBinsSector){
	//std::cout << "WARNING phibin out of range: " << phibin << " | " << phibins << std::endl;
	continue;
      }
      if(zbin>=NZBinsSide){
	//std::cout << "WARNING z bin out of range: " << zbin << " | " << zbins << std::endl;
	continue;
      }
      if(zbinorg>zbinmax||zbinorg<zbinmin)
      	continue;
      float_t fadc = (hitr->second->getAdc()) - pedestal; // proper int rounding +0.5
      unsigned short adc = 0;
      if(fadc>0) adc =  (unsigned short) fadc;

      if(adc>0){
	nhits_thispad[phibin]++;
	unsigned short it = TpcDefs::getTBin(hitr->first);
	adcval[phibin][it] = adc;
	//if(layer==20&&sector==4&&phibin==77)
	//std::cout << "input sector: " << sector << " lay: " << layer << " # " << phibin << "|"<< it <<" adc: " << (int)adcval[phibin][it] << " pack: " << std::endl;
	/*RawHitTpc  *rHit = new RawHitTpc;
	unsigned short it = TpcDefs::getTBin(hitr->first);
	rHit->setTBin(it);
	unsigned short iadc = adc;
	rHit->setAdc(iadc);
	rhitset->addTpcHit(phibin,rHit);
	*/
      }
    }
    //    assert(adcval);
    for(int nphi= 0; nphi < NPhiBins;nphi++){
      for(int nz= 0; nz < NZBins;nz++){
	/*if(nphi>=1&&nz<NZBins-1)//Remove isolated single timebin hits
	  if(layer==29&&sector==0&&nphi==0){
	    if(nz>209&&nz<215){
	      std::cout << "single check sector: " << sector << " lay: " << layer << " # " << nphi << "|"<< nz <<" adc: " << (int)adcval[nphi][nz] << " pack: " << std::endl;
	    }
	  }
	*/
	if(adcval[nphi][nz]>0)
	  if((adcval[nphi][nz-1]==0)&&(adcval[nphi][nz+1]==0)){
	    // std::cout << "single hit sector: " << sector << " lay: " << layer << " # " << nphi << "|"<< nz <<" adc: " << (int)adcval[nphi][nz] << " pack: " << std::endl;
	    adcval[nphi][nz] = 0;
	  }
      }
    }
    for(int nphi= 0; nphi < NPhiBins;nphi++){
      int zero_count = 0;
      // int outpos = 0;
      if(nhits_thispad[nphi]==0)continue;
      for(int nz= 0; nz < NZBins;nz++){

	uint8_t thisadc = 0;
	if(nphi>=1&&nz<NZBins-1)//Remove isolated single timebin hits
	  if(adcval[nphi][nz]>0)
	    if((adcval[nphi][nz-1]==0)&&(adcval[nphi][nz+1]==0))
	      adcval[nphi][nz] = 0;

	if(adcval[nphi][nz]>255){//take care of overflows, we want to store 8 bit ADC values
	    thisadc=255;
	    //	    adcval[nphi][nz]=255;
	}
	else{
	  thisadc = adcval[nphi][nz];
	}
	//std::cout << "nz: " << nz << " adc " << (int)thisadc << std::endl;
	if(thisadc>0){//non zero adc, fill adc
	  zero_count = 0;
	  rhitset->m_tpchits[nphi].push_back(thisadc);
	  if(layer==222&&sector==6&&nphi==10)std::cout << "0#"<<nphi<<"nz= " << nz << " filling " << (int)thisadc << " at " << rhitset->m_tpchits[nphi].size() << std::endl;
	}
	else{
	  zero_count++;
	  //	  if(nz==0){
	  //  rhitset->m_tpchits[nphi].push_back(0);
	    //std::cout << " 1filling " << 0 << "|" << rhitset->m_tpchits[nphi].back() << " | " << rhitset->m_tpchits[nphi][outpos++] << std::endl;
	  // if(layer==222&&sector==6&&nphi==10)std::cout << "1#"<<nphi<<"nz= " << nz << " filling " << 0 << " at " << rhitset->m_tpchits[nphi].size() << std::endl;
	  // }else 
	  /*
	  if(adcval[nphi][nz-1]>0){ //trailing edge 0
	    rhitset->m_tpchits[nphi].push_back(0);
	    if(layer==222&&sector==6&&nphi==10)std::cout << "2#"<<nphi<<"nz= " << nz << " filling " << 0 << " at " << rhitset->m_tpchits[nphi].size() << std::endl;
	    //std::cout << " 2filling " << 0 << "|" << rhitset->m_tpchits[nphi].back() << " | " << rhitset->m_tpchits[nphi][outpos++] << std::endl;
	  }
	  */
	}
	if(zero_count>1){
	  if(nz<NZBins-2){
	    if(adcval[nphi][nz+2]>0&&adcval[nphi][nz+1]==0){//fill zero count, end of zero series
	      rhitset->m_tpchits[nphi].push_back(zero_count-1);
	      zero_count = 0;
	      if(layer==222&&sector==6&&nphi==10)
		std::cout << "3#"<<nphi<<"nz= " << nz << " filling " << zero_count-1 << " at " << rhitset->m_tpchits[nphi].size() << " nz+2 " << (nz+2) << " adcval[nz+2]" << (int) adcval[nphi][nz+2] << std::endl;
	    }
	    if(adcval[nphi][nz+1]>0){
	      rhitset->m_tpchits[nphi].push_back(0);//leading edge zero
	      if(layer==222&&sector==6&&nphi==10)std::cout << "4#"<<nphi<<"nz= " << nz << " filling " << 0 << " at " << rhitset->m_tpchits[nphi].size() << std::endl;

	      //std::cout << " 4filling " << 0 << "|" << ((int)rhitset->m_tpchits[nphi].back()) << " | " << ((int)rhitset->m_tpchits[nphi][outpos++]) << std::endl;
	    }
	  }

	}else{
	  if(zero_count==1){
	    rhitset->m_tpchits[nphi].push_back(0);
	    if(layer==222&&sector==6&&nphi==10)std::cout << "5#"<<nphi<<"nz= " << nz << " filling " << 0 << " at " << rhitset->m_tpchits[nphi].size() << std::endl;
	    //std::cout << " 5filling " << 0 << "|" << ((int)rhitset->m_tpchits[nphi].back()) << " | " << ((int)rhitset->m_tpchits[nphi][outpos++]) << std::endl;
	  }
	}
	if(zero_count == 254){
	  zero_count = 0;
	  rhitset->m_tpchits[nphi].push_back(254-2);
	  rhitset->m_tpchits[nphi].push_back(0);
	  if(layer==222&&sector==6&&nphi==10)std::cout << "6#"<<nphi<<"nz= " << nz << " filling " << 254-1 << " at " << rhitset->m_tpchits[nphi].size() << std::endl;
	  //std::cout << "6 filling " << 254-1 << "|" << ((int)rhitset->m_tpchits[nphi].back()) << " | " << ((int)rhitset->m_tpchits[nphi][outpos++]) << std::endl;
	  // rhitset->m_tpchits[nphi].push_back(0);
	  // if(layer==222&&sector==6&&nphi==10)std::cout << "7#"<<nphi<<"nz= " << nz << " filling " << 0 << " at " << rhitset->m_tpchits[nphi].size() << std::endl;
	  //std::cout << " 7filling " << 0 << "|" << ((int)rhitset->m_tpchits[nphi].back()) << " | " << ((int)rhitset->m_tpchits[nphi][outpos++]) << std::endl;
	}
	/*
	  if(layer%15==0){
	  std::cout << "nz #" << nz << " adc " << adcval[nphi][nz] << " n zero " << zero_count << " last (" << rhitset->m_tpchits[nphi].size() << ") " << ((int)rhitset->m_tpchits[nphi].back()) << std::endl;
	if(rhitset->m_tpchits[nphi].size()>=3)
	std::cout << " [ " << rhitset->m_tpchits[nphi].size()-1 << ": "<< ((int)rhitset->m_tpchits[nphi][rhitset->m_tpchits[nphi].size()-1]) << " ] "
		    << " [ " << rhitset->m_tpchits[nphi].size()-2 << ": "<< ((int)rhitset->m_tpchits[nphi][rhitset->m_tpchits[nphi].size()-2]) << " ] "
		    << " [ " << rhitset->m_tpchits[nphi].size()-3<< ": "<< ((int)rhitset->m_tpchits[nphi][rhitset->m_tpchits[nphi].size()-3]) << " ] "
		    << std::endl;
		    }
	*/
      }
      if(zero_count>0&&rhitset->m_tpchits[nphi].back()!=0)rhitset->m_tpchits[nphi].push_back(0); //pad zero at the end
      //      std::cout << "sector: " << sector << " lay: " << layer << " nphi: " << nphi << "|" << NPhiBins<< " nzfilled: " << rhitset->m_tpchits[nphi].size() << std::endl;
    }
    count++;

    //    std::cout << "unpacking..." << std::endl;
    //    unpack and check for correctness
    std::vector<std::vector<uint8_t>> outval;
    outval.resize(NPhiBins);
    for(int nphi= 0; nphi < NPhiBins;nphi++){
      outval[nphi].resize(NZBins,0);
      for(int nz= 0; nz < NZBins;nz++){
	if(outval[nphi][nz]!=0)
	  std::cout << "WARNING!" << std::endl;
      }
    }
    //   now we have a clean output array
    for(int nphi= 0; nphi < NPhiBins;nphi++){
      if(rhitset->m_tpchits[nphi].size()==0) continue;

      int pindex = 0;
      for(unsigned int nzo = 0;nzo<rhitset->m_tpchits[nphi].size();nzo++){
	uint8_t val = rhitset->m_tpchits[nphi][nzo];

	if(val==0)
	  pindex++;
	else{
	  if(nzo==0){
	    outval[nphi][pindex++]=val;
	  }else{
	    if((rhitset->m_tpchits[nphi][nzo-1]==0)&&(rhitset->m_tpchits[nphi][nzo+1]==0))//found zero count
	      pindex+=val;
	    else{
	      outval[nphi][pindex++]=val;
	    }
	  }
	}
      }
    }
    {
      // std::cout << "checking..." << std::endl;
      for(int nphi= 0; nphi <  NPhiBins;nphi++){
	// std::cout << "nphi: " << nphi << "|" << NPhiBins<< std::endl;
	// std::cout << " adcval size "  << adcval[nphi].size() << std::endl;
	// std::cout << " outval size "  << outval[nphi].size() << std::endl;
	int bad = 0;
	for(int nz= 0; nz < NZBins;nz++){
	  ncheck++;
	  if(adcval[nphi][nz]!=outval[nphi][nz]){
	    if(outval[nphi][nz]!=255){
	      bad++;
	      allbad++;
	    }
	  }
	}
	
	if(bad>0){
	  std::cout << "sector: " << sector << " lay: " << layer << " # " << nphi << " bad:" << bad << " packsize: "<<  (int)rhitset->m_tpchits[nphi].size() << std::endl;
	  
	  for(int nz= 0; nz < NZBins;nz++){
	    std::cout << "sector: " << sector << " lay: " << layer << " # " << nphi << "|"<< nz << " bad:" << bad <<" org: " << (int)adcval[nphi][nz] << " pack: "<< (int)outval[nphi][nz] << std::endl;
	  }
	  for(int nz= 0; nz < (int)rhitset->m_tpchits[nphi].size() ;nz++){
	    std::cout << "lay: " << layer << " # " << nphi << "|"<< nz << " bad:" << bad <<" pack: " << (int) rhitset->m_tpchits[nphi][nz] << std::endl;
	  }
	  
	}
	
      }
      
    }
    
    //#ifdef WORK     //    if(layer==7)


      }
  std::cout << " number of hits checked "  << ncheck << " bad: " << allbad << std::endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcRawWriter::End(PHCompositeNode */*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
