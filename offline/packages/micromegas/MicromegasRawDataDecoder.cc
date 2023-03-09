/*!
 * \file MicromegasRawDataDecoder.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "MicromegasRawDataDecoder.h"

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>

#include <trackbase/TrkrHitSetContainerv1.h>

#include <cassert>

//_________________________________________________________
MicromegasRawDataDecoder::MicromegasRawDataDecoder( const std::string& name ):
  SubsysReco( name )
{}

//_____________________________________________________________________
int MicromegasRawDataDecoder::Init(PHCompositeNode* /*topNode*/ )
{ 
  // read calibrations
  m_calibration_data.read( m_calibration_filename );
  return Fun4AllReturnCodes::EVENT_OK; 
}

//____________________________________________________________________________..
int MicromegasRawDataDecoder::InitRun(PHCompositeNode *topNode)
{

  // get dst node
  PHNodeIterator iter(topNode);
  auto dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << "MicromegasRawDataDecoder::InitRun - DST Node missing, doing nothing." << std::endl;
    exit(1);
  }

  // create hitset container if needed
  auto hitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!hitsetcontainer)
  {
    std::cout << "MicromegasRawDataDecoder::InitRun - creating TRKR_HITSET." << std::endl;

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
int MicromegasRawDataDecoder::process_event(PHCompositeNode *topNode)
{

  // load relevant nodes
  // Get the TrkrHitSetContainer node
  auto trkrhitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  assert(trkrhitsetcontainer);

  // PRDF node
  auto event = findNode::getClass<Event>(topNode, "PRDF");
  assert( event );

  // check event type
  if(event->getEvtType() >= 8)
  { return Fun4AllReturnCodes::DISCARDEVENT; }

  // get TPOT packet number
  /*
   * for now it is the same packet number as the TPC: 4001.
   * To be fixed at a later stage.
   * check with Martin Purschke
   */
  auto packet = event->getPacket(4001);
  if( !packet )
  {
    // no data
    std::cout << "MicromegasRawDataDecoder::process_event - event contains no TPOT data" << std::endl;
    return Fun4AllReturnCodes::EVENT_OK;
  }

  // get number of datasets (also call waveforms)
  const auto n_waveforms = packet->iValue(0, "NR_WF" );
  if( Verbosity() )
  { std::cout << "MicromegasRawDataDecoder::process_event - n_waveforms: " << n_waveforms << std::endl; }
  for( int i=0; i<n_waveforms; ++i )
  {
    auto channel = packet->iValue( i, "CHANNEL" );
    int fee = packet->iValue(i, "FEE" );
    int samples = packet->iValue( i, "SAMPLES" );
    if( Verbosity() )
    {
      std::cout
        << "MicromegasRawDataDecoder::process_event -"
        << " waveform: " << i
        << " fee: " << fee
        << " channel: " << channel
        << " samples: " << samples
        << std::endl;
    }
    
    // map fee and channel to physical hitsetid and physical strip
    // process signal to get mean ADC
    // save as hit
    
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
