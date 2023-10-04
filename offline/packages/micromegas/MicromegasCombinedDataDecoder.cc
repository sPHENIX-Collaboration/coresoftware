/*!
 * \file MicromegasCombinedDataDecoder.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 */

#include "MicromegasCombinedDataDecoder.h"
#include "MicromegasDefs.h"

#include <ffarawobjects/MicromegasRawHit.h>
#include <ffarawobjects/MicromegasRawHitContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>

#include <trackbase/TrkrHitv2.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainerv1.h>

#include <algorithm>
#include <cassert>
#include <memory>

//_________________________________________________________
MicromegasCombinedDataDecoder::MicromegasCombinedDataDecoder( const std::string& name ):
  SubsysReco( name )
{}

//_____________________________________________________________________
int MicromegasCombinedDataDecoder::Init(PHCompositeNode* /*topNode*/ )
{ 
  // read calibrations
  m_calibration_data.read( m_calibration_filename );
  return Fun4AllReturnCodes::EVENT_OK; 
}

//____________________________________________________________________________..
int MicromegasCombinedDataDecoder::InitRun(PHCompositeNode *topNode)
{

  // get dst node
  PHNodeIterator iter(topNode);
  auto dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << "MicromegasCombinedDataDecoder::InitRun - DST Node missing, doing nothing." << std::endl;
    exit(1);
  }

  // create hitset container if needed
  auto hitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!hitsetcontainer)
  {
    std::cout << "MicromegasCombinedDataDecoder::InitRun - creating TRKR_HITSET." << std::endl;

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
int MicromegasCombinedDataDecoder::process_event(PHCompositeNode *topNode)
{

  // load relevant nodes
  // Get the TrkrHitSetContainer node
  auto trkrhitsetcontainer = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  assert(trkrhitsetcontainer);

  // load raw hits container  
  auto rawhitcontainer = findNode::getClass<MicromegasRawHitContainer>(topNode, m_rawhitnodename);
  assert( rawhitcontainer );

  // loop over raw hits
  if( Verbosity() )
  { std::cout << "MicromegasCombinedDataDecoder::process_event - hits: " << rawhitcontainer->get_nhits() << std::endl; }

  for( unsigned int ihit = 0; ihit < rawhitcontainer->get_nhits(); ++ihit )
  {
    const auto rawhit = rawhitcontainer->get_hit(ihit);
    const auto packet_id = rawhit->get_packetid();

    // make sure packet is valid
    if( std::find( std::begin(MicromegasDefs::m_packet_ids), std::end(MicromegasDefs::m_packet_ids ), packet_id) == std::end(MicromegasDefs::m_packet_ids ) )
    {
      std::cout << "MicromegasCombinedDataDecoder::process_event - invalid packet: " << packet_id << std::endl;
      continue;
    }
      
    const int fee = rawhit->get_fee();
    const auto channel = rawhit->get_channel();
    const int samples = rawhit->get_samples();
          
    // get channel rms and pedestal from calibration data
    const double pedestal = m_calibration_data.get_pedestal( fee, channel );
    const double rms = m_calibration_data.get_rms( fee, channel );
    
    // a rms of zero means the calibration has failed. the data is unusable
    if( rms <= 0 ) continue;
    
    // loop over samples find maximum
    std::vector<int> adc;
    for( int is = std::max( m_sample_min,0 ); is < std::min( m_sample_max,samples ); ++ is )
    { adc.push_back(rawhit->get_adc(is)); }
    
    if( adc.empty() ) continue;
    
    // get max adc value in range
    /* TODO: use more advanced signal processing */
    auto max_adc = *std::max_element( adc.begin(), adc.end() );
            
    // compare to hard min_adc value
    if( max_adc < m_min_adc ) continue;

    // compare to threshold
    if( max_adc < pedestal - m_n_sigma * rms ) continue;
    
    // subtract pedestal
    max_adc -= pedestal;
    
    // map fee and channel to physical hitsetid and physical strip
    // get hitset key matching this fee
    const TrkrDefs::hitsetkey hitsetkey = m_mapping.get_hitsetkey( fee );     
    if( !hitsetkey ) continue;
    
    // get matching physical strip
    int strip = m_mapping.get_physical_strip( fee, channel );
    if( strip < 0 ) continue;
    
    if( Verbosity() )
    {
      const auto bco = rawhit->get_gtm_bco();
      std::cout << "MicromegasCombinedDataDecoder::process_event -"
        << " bco: " << bco
        << " layer: " << int(TrkrDefs::getLayer(hitsetkey)) 
        << " tile: " << int( MicromegasDefs::getTileId( hitsetkey ))
        << " channel: " << channel 
        << " strip: " << strip 
        << " adc: " << max_adc
        << std::endl;
    }
    
    // get matching hitset
    const auto hitset_it = trkrhitsetcontainer->findOrAddHitSet(hitsetkey);
    
    // generate hit key
    const TrkrDefs::hitkey hitkey = MicromegasDefs::genHitKey(strip);
    
    // find existing hit, or create
    auto hit = hitset_it->second->getHit(hitkey);
    if( hit ) 
    {
      // std::cout << "MicromegasCombinedDataDecoder::process_event - duplicated hit, hitsetkey: " << hitsetkey << " strip: " << strip << std::endl;
      continue;
    }
      
    // create hit, assign adc and insert in hitset
    hit = new TrkrHitv2;
    hit->setAdc( max_adc );
    hitset_it->second->addHitSpecificKey(hitkey, hit);
    
    // increment counter
    ++m_hitcounts[hitsetkey];
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int MicromegasCombinedDataDecoder::End(PHCompositeNode* /*topNode*/ )
{
  // if( Verbosity() )
  {
    for( const auto& [hitsetkey, count]:m_hitcounts )
    { std::cout << "MicromegasCombinedDataDecoder::End - hitsetkey: " << hitsetkey << ", count: " << count << std::endl; }
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}
