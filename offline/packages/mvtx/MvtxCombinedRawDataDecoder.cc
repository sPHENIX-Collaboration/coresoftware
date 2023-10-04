/*!
 * \file MvtxCombinedRawDataDecoder.cc
 * \author Jakub Kvapil <jakub.kvapil@cern.ch>
 */

#include "MvtxCombinedRawDataDecoder.h"

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>

#include <trackbase/MvtxDefs.h>
#include <trackbase/TrkrHitv2.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainerv1.h>

#include <ffarawobjects/MvtxRawHit.h>
#include <ffarawobjects/MvtxRawHitContainer.h>

#include <algorithm>
#include <cassert>

//_________________________________________________________
MvtxCombinedRawDataDecoder::MvtxCombinedRawDataDecoder( const std::string& name ):
  SubsysReco( name )
{}

//_____________________________________________________________________
int MvtxCombinedRawDataDecoder::Init(PHCompositeNode* /*topNode*/ )
{ 
  return Fun4AllReturnCodes::EVENT_OK; 
}

//____________________________________________________________________________..
int MvtxCombinedRawDataDecoder::InitRun(PHCompositeNode *topNode)
{

  // get dst node
  PHNodeIterator iter(topNode);
  auto dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << "MvtxCombinedRawDataDecoder::InitRun - DST Node missing, doing nothing." << std::endl;
    exit(1);
  }

  // create hitset container if needed
  auto hitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!hitsetcontainer)
  {
    // find or create TRKR node
    PHNodeIterator dstiter(dstNode);
    auto trkrnode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
    if (!trkrnode)
    {
      trkrnode = new PHCompositeNode("TRKR");
      dstNode->addNode(trkrnode);
    }

    // create container and add to the tree
    hitsetcontainer = new TrkrHitSetContainerv1;
    auto newNode = new PHIODataNode<PHObject>(hitsetcontainer, "TRKR_HITSET", "PHObject");
    trkrnode->addNode(newNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;

}

//___________________________________________________________________________
int MvtxCombinedRawDataDecoder::process_event(PHCompositeNode *topNode)
{
  // load relevant nodes
  // Get the TrkrHitSetContainer node
  auto trkrhitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  assert(trkrhitsetcontainer);

  MvtxRawHitContainer* mvtx_hit_container = findNode::getClass<MvtxRawHitContainer>(topNode, m_MvtxRawNodeName);
  if (!mvtx_hit_container)
  {
    std::cout << PHWHERE << "::" << __func__ <<  ": Could not get \"" << m_MvtxRawNodeName << "\" from Node Tree" << std::endl;
    std::cout << "Have you built this yet?" << std::endl;
    exit(1);
  }
  
  if (Verbosity() >= VERBOSITY_MORE) mvtx_hit_container->identify();

  int strobe = -1; //Initialise to -1 for debugging
  uint8_t layer = 0; 
  uint8_t stave = 0;
  uint8_t chip = 0;
  uint16_t row = 0;
  uint16_t col = 0;

  for (unsigned int i = 0; i < mvtx_hit_container->get_nhits(); i++)
  {
    MvtxRawHit* mvtx_hit = mvtx_hit_container->get_hit(i);

    strobe = mvtx_hit->get_bco();
    layer = mvtx_hit->get_layer_id();
    stave = mvtx_hit->get_stave_id();
    chip = mvtx_hit->get_chip_id();
    row = mvtx_hit->get_chip_id();
    col = mvtx_hit->get_col();

    if( Verbosity() >= VERBOSITY_A_LOT ){
      std::cout
        << "MVTX raw hit:"
        << " strobe: " << strobe
        << " layer: " << layer
        << " stave: " << stave
        << " chip: " << chip
        << " row: " << row
        << " col: " << col
        << std::endl;
    }
        
    const TrkrDefs::hitsetkey hitsetkey = MvtxDefs::genHitSetKey(layer, stave, chip, strobe);     
    if( !hitsetkey ) continue;

    // get matching hitset
    const auto hitset_it = trkrhitsetcontainer->findOrAddHitSet(hitsetkey);

    // generate hit key
    const TrkrDefs::hitkey hitkey = MvtxDefs::genHitKey(col,row);
    
    // find existing hit, or create
    auto hit = hitset_it->second->getHit(hitkey);
    if( hit ){
      std::cout << PHWHERE << "::" << __func__ << " - duplicated hit, hitsetkey: " << hitsetkey << " hitkey: " << hitkey << std::endl;
      continue;
    }

    // create hit and insert in hitset
    hit = new TrkrHitv2;
    hitset_it->second->addHitSpecificKey(hitkey, hit);
    
    // increment counter
    ++m_hitcounts[hitsetkey];
             
  }

  return Fun4AllReturnCodes::EVENT_OK;

}

//_____________________________________________________________________
int MvtxCombinedRawDataDecoder::End(PHCompositeNode* /*topNode*/ )
{
  if( Verbosity() )
  {
    for( const auto& [hitsetkey, count]:m_hitcounts )
    { std::cout << "MvtxCombinedRawDataDecoder - hitsetkey: " << hitsetkey << ", count: " << count << std::endl; }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
