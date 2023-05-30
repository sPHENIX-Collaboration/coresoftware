/*!
 * \file MvtxRawDataDecoder.cc
 * \author Jakub Kvapil <jakub.kvapil@cern.ch>
 */

#include "MvtxRawDataDecoder.h"

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

#include <algorithm>
#include <cassert>

//_________________________________________________________
MvtxRawDataDecoder::MvtxRawDataDecoder( const std::string& name ):
  SubsysReco( name )
{}

//_____________________________________________________________________
int MvtxRawDataDecoder::Init(PHCompositeNode* /*topNode*/ )
{ 
  return Fun4AllReturnCodes::EVENT_OK; 
}

//____________________________________________________________________________..
int MvtxRawDataDecoder::InitRun(PHCompositeNode *topNode)
{

  // get dst node
  PHNodeIterator iter(topNode);
  auto dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << "MvtxRawDataDecoder::InitRun - DST Node missing, doing nothing." << std::endl;
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
int MvtxRawDataDecoder::process_event(PHCompositeNode *topNode)
{

  // load relevant nodes
  // Get the TrkrHitSetContainer node
  auto trkrhitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  assert(trkrhitsetcontainer);

  // PRDF node
  auto event = findNode::getClass<Event>(topNode, "PRDF");
  assert(event);
  // check event type


//UNCOMMENT THIS ONCE HAVING DECODER!!!!
/*  if(event->getEvtType() >= 8)
  { return Fun4AllReturnCodes::DISCARDEVENT; }

  // get MVTX packet number
  auto packet = event->getPacket(2001);
  if( !packet )
  {
    // no data
    std::cout << "MvtxRawDataDecoder::process_event - event contains no MVTX data" << std::endl;
   // return Fun4AllReturnCodes::EVENT_OK;
  }

*/

  //loop over triggers ADD
  //loop over chips 
  for( int i=0; i<1; i++ )
  {  
    int strobe = 11+i; //get from decoder;
    int layer = 0;//get from decoder;
    int stave = 1;//get from decoder;
    int chip = 2;//get from decoder;

    if( Verbosity() ){
      std::cout
        << "MvtxRawDataDecoder::process_event -"
        << " strobe: " << strobe
        << " layer: " << layer
        << " stave: " << stave
        << " chip: " << chip
        << std::endl;
    }
        
    // map fee and channel to physical hitsetid and physical strip
    // get hitset key matching this fee
    const TrkrDefs::hitsetkey hitsetkey = MvtxDefs::genHitSetKey(layer, stave, chip, strobe);     
    if( !hitsetkey ) continue;

    // get matching hitset
    const auto hitset_it = trkrhitsetcontainer->findOrAddHitSet(hitsetkey);

    for( int ihit=0; ihit<10/*nhit in chip from decoder*/; ++ihit ){
      uint16_t col = 150;//get from decoder;
      uint16_t row = ihit*10;//get from decoder;

      // generate hit key
      const TrkrDefs::hitkey hitkey = MvtxDefs::genHitKey(col,row);
    
      // find existing hit, or create
      auto hit = hitset_it->second->getHit(hitkey);
      if( hit ){
        std::cout << "MvtxRawDataDecoder::process_event - duplicated hit, hitsetkey: " << hitsetkey << " hitkey: " << hitkey << std::endl;
        continue;
      }

      // create hit and insert in hitset
      hit = new TrkrHitv2;
      hitset_it->second->addHitSpecificKey(hitkey, hit);
      
      // increment counter
      ++m_hitcounts[hitsetkey];
             
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;

}

//_____________________________________________________________________
int MvtxRawDataDecoder::End(PHCompositeNode* /*topNode*/ )
{
  if( Verbosity() )
  {
    for( const auto& [hitsetkey, count]:m_hitcounts )
    { std::cout << "MvtxRawDataDecoder - hitsetkey: " << hitsetkey << ", count: " << count << std::endl; }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
