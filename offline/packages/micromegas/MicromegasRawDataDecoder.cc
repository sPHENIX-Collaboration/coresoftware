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

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>

#include <cassert>

//_________________________________________________________
MicromegasRawDataDecoder::MicromegasRawDataDecoder( const std::string& name ):
  SubsysReco( name )
{}

//_____________________________________________________________________
int MicromegasRawDataDecoder::Init(PHCompositeNode* /*topNode*/ )
{

  // histogram evaluation
  if( m_savehistograms ) create_histograms();
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

  // map fee id to detector index in histogram
  using fee_map_t = std::map<int,int>;
  fee_map_t fee_map = {
    {2, 0},
    {1, 1},
    {3, 2},
    {4, 3},
    {8, 4},
    {9, 5},
    {7, 6},
    {5, 7}
  };

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

    if( m_savehistograms && m_h_fee_id ) m_h_fee_id->Fill(fee);

    // find fee index from map
    const auto iter = fee_map.find( fee );
    if( iter == fee_map.end() )
    {
      std::cout << "MicromegasRawDataDecoder::process_event - unable to find fee " << fee << " in map" << std::endl;
      continue;
    }

    const auto fee_index = iter->second;
    const auto channel_index = fee_index*m_nchannels_fee + channel;

    // loop over samples
    if( m_savehistograms && m_h_adc_channel )
    {
      for( int is = 0; is < samples; ++ is )
      { m_h_adc_channel->Fill( channel_index, packet->iValue(i,is) ); }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
int MicromegasRawDataDecoder::End(PHCompositeNode* /*topNode*/ )
{
  // save histograms
  if( m_savehistograms && m_histogramfile )
  {
    // create mean and rms histograms
    auto profile = m_h_adc_channel->ProfileX("h_adc_channel_profx", 1, -1, "s" );
    auto h_pedestal = new TH1F( "h_pedestal", "pedestal vs channel", m_nchannels_total, 0, m_nchannels_total );
    auto h_rms = new TH1F( "h_rms", "rms vs channel", m_nchannels_total, 0, m_nchannels_total );
    for( int i =0; i<m_nchannels_total; ++i )
    {
      h_pedestal->SetBinContent( i+1, profile->GetBinContent(i+1) );
      h_rms->SetBinContent(i+1, profile->GetBinError(i+1) );
    }

    m_histogramfile->cd();
    m_h_fee_id->Write();
    m_h_adc_channel->Write();
    h_pedestal->Write();
    h_rms->Write();
    m_histogramfile->Close();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//_____________________________________________________________________
void MicromegasRawDataDecoder::create_histograms()
{
  std::cout << "MicromegasRawDataDecoder::create_histograms - writing evaluation histograms to: " << m_histogramfilename << std::endl;
  m_histogramfile.reset( new TFile(m_histogramfilename.c_str(), "RECREATE") );
  m_histogramfile->cd();

  m_h_fee_id = new TH1I( "h_fee_id", "FEE id", 10, 0, 10 );
  m_h_adc_channel = new TH2I( "h_adc_channel", "ADC vs channel", m_nchannels_total, 0, m_nchannels_total, m_max_adc, 0, m_max_adc );
}
