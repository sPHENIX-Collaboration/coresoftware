/*!
 * \file MvtxCombinedRawDataDecoder.cc
 * \author Cameron Dean <cameron.dean@cern.ch@cern.ch>
 * \author Jakub Kvapil <jakub.kvapil@cern.ch>
 */

#include "MvtxCombinedRawDataDecoder.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>

#include <trackbase/TrkrHitv2.h>
#include <trackbase/TrkrHitSet.h>

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
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << "MvtxCombinedRawDataDecoder::InitRun - DST Node missing, doing nothing." << std::endl;
    exit(1);
  }

  // create hitset container if needed
  hit_set_container = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!hit_set_container)
  {
    // find or create TRKR node
    auto trkrNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "TRKR"));
    if (!trkrNode)
    {
      trkrNode = new PHCompositeNode("TRKR");
      dstNode->addNode(trkrNode);
    }

    // create container and add to the tree
    hit_set_container = new TrkrHitSetContainerv1;
    auto newNode = new PHIODataNode<PHObject>(hit_set_container, "TRKR_HITSET", "PHObject");
    trkrNode->addNode(newNode);
  }

  //Check if MVTX event header already exists
  if (m_writeMvtxEventHeader)
  {
    auto mvtxNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "MVTX"));
    if (!mvtxNode)
    {
      mvtxNode = new PHCompositeNode("MVTX");
      dstNode->addNode(mvtxNode);
    }

    mvtx_event_header = findNode::getClass<MvtxEventInfov1>(mvtxNode, "MVTXEVENTHEADER");
    if (!mvtx_event_header)
    {
      mvtx_event_header = new MvtxEventInfov1();
      auto newHeader = new PHIODataNode<PHObject>(mvtx_event_header, "MVTXEVENTHEADER", "PHObject");
      mvtxNode->addNode(newHeader);
    }
  }

  mvtx_raw_event_header = findNode::getClass<MvtxRawEvtHeaderv1>(topNode, m_MvtxRawEvtHeaderNodeName);
  if (!mvtx_raw_event_header)
  {
    std::cout << PHWHERE << "::" << __func__ <<  ": Could not get \"" << m_MvtxRawEvtHeaderNodeName << "\" from Node Tree" << std::endl;
    std::cout << "Have you built this yet?" << std::endl;
    exit(1);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//___________________________________________________________________________
int MvtxCombinedRawDataDecoder::process_event(PHCompositeNode *topNode)
{
  mvtx_raw_event_header = findNode::getClass<MvtxRawEvtHeaderv1>(topNode, m_MvtxRawEvtHeaderNodeName);
  if (Verbosity() >= VERBOSITY_MORE) mvtx_raw_event_header->identify();

  hit_set_container = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");

  mvtx_hit_container = findNode::getClass<MvtxRawHitContainer>(topNode, m_MvtxRawHitNodeName);
  if (!mvtx_hit_container)
  {
    std::cout << PHWHERE << "::" << __func__ <<  ": Could not get \"" << m_MvtxRawHitNodeName << "\" from Node Tree" << std::endl;
    std::cout << "Have you built this yet?" << std::endl;
    exit(1);
  }
  
  if (Verbosity() >= VERBOSITY_MORE) mvtx_hit_container->identify();

  uint64_t strobe = -1; //Initialise to -1 for debugging
  uint8_t layer = 0; 
  uint8_t stave = 0;
  uint8_t chip = 0;
  uint16_t row = 0;
  uint16_t col = 0;
  std::vector<std::pair<uint64_t, uint32_t>> strobe_bc_pairs;

  if (m_writeMvtxEventHeader)
  {
    mvtx_event_header = findNode::getClass<MvtxEventInfov1>(topNode, "MVTXEVENTHEADER");
    assert(mvtx_event_header);
  }

  for (unsigned int i = 0; i < mvtx_hit_container->get_nhits(); i++)
  {
    mvtx_hit = mvtx_hit_container->get_hit(i);

    strobe = mvtx_hit->get_bco();
    layer = mvtx_hit->get_layer_id();
    stave = mvtx_hit->get_stave_id();
    chip = mvtx_hit->get_chip_id();
    row = mvtx_hit->get_row();
    col = mvtx_hit->get_col();

    if( Verbosity() >= VERBOSITY_A_LOT ) mvtx_hit->identify();
        
    const TrkrDefs::hitsetkey hitsetkey = MvtxDefs::genHitSetKey(layer, stave, chip, strobe);     
    if( !hitsetkey ) continue;

    // get matching hitset
    const auto hitset_it = hit_set_container->findOrAddHitSet(hitsetkey);

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
  }

  std::set<uint64_t> l1BCOs = mvtx_raw_event_header->getMvtxLvL1BCO();

  if (m_writeMvtxEventHeader)
  {
    std::set<uint64_t> l1BCOs = mvtx_raw_event_header->getMvtxLvL1BCO();
    for (auto iter = l1BCOs.begin(); iter != l1BCOs.end(); iter++)
    {
      mvtx_event_header->set_strobe_BCO_L1_BCO(strobe, *iter);
    }
    if (Verbosity() >= VERBOSITY_EVEN_MORE) mvtx_event_header->identify();
  } 
  
  return Fun4AllReturnCodes::EVENT_OK;

}

//_____________________________________________________________________
int MvtxCombinedRawDataDecoder::End(PHCompositeNode* /*topNode*/ )
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void MvtxCombinedRawDataDecoder::removeDuplicates(std::vector<std::pair<uint64_t, uint32_t>> &v)
{
  auto end = v.end();
  for (auto it = v.begin(); it != end; ++it)
  {
    end = remove(it + 1, end, *it);
  }
  v.erase(end, v.end());
}
