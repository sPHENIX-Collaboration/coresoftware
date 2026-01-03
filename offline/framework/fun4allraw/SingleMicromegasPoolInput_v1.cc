#include "SingleMicromegasPoolInput_v1.h"

#include "Fun4AllStreamingInputManager.h"
#include "InputManagerType.h"

#include <qautils/QAHistManagerDef.h>
#include <qautils/QAUtil.h>

#include <ffarawobjects/MicromegasRawHitContainerv3.h>
#include <ffarawobjects/MicromegasRawHitv3.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/InputFileHandlerReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/getClass.h>
#include <phool/phool.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/Eventiterator.h>
#include <Event/fileEventiterator.h>

#include <TFile.h>
#include <TH1.h>

#include <algorithm>
#include <bitset>

namespace
{

  //! mark invalid ADC values
  static constexpr uint16_t m_adc_invalid = 65000;

  // maximum number of packets
  static constexpr int m_npackets_active = 2;

  // minimum number of requested samples
  static constexpr int m_min_req_samples = 5;

  /* see: https://git.racf.bnl.gov/gitea/Instrumentation/sampa_data/src/branch/fmtv2/README.md */
  enum SampaDataType
  {
    HEARTBEAT_T = 0b000,
    TRUNCATED_DATA_T = 0b001,
    TRUNCATED_TRIG_EARLY_DATA_T = 0b011,
    NORMAL_DATA_T = 0b100,
    LARGE_DATA_T = 0b101,
    TRIG_EARLY_DATA_T = 0b110,
    TRIG_EARLY_LARGE_DATA_T = 0b111,
  };

}  // namespace

//______________________________________________________________
SingleMicromegasPoolInput_v1::SingleMicromegasPoolInput_v1(const std::string& name)
  : SingleStreamingInput(name)
{
  SubsystemEnum(InputManagerType::MICROMEGAS);
  m_rawHitContainerName = "MICROMEGASRAWHIT";
}

//______________________________________________________________
SingleMicromegasPoolInput_v1::~SingleMicromegasPoolInput_v1()
{

  std::cout << "SingleMicromegasPoolInput_v1::~SingleMicromegasPoolInput_v1 - runnumber: " << RunNumber() << std::endl;

  // timer statistics
  m_timer.print_stat();

  for( const auto& [packet,counts]:m_waveform_count_total )
  {
    const auto dropped_bco =  m_waveform_count_dropped_bco[packet];
    const auto dropped_pool =  m_waveform_count_dropped_pool[packet];
    std::cout << "SingleMicromegasPoolInput_v1::~SingleMicromegasPoolInput_v1 -"
      << " packet: " << packet
      << " wf_total: " << counts
      << " wf_dropped_bco: " << dropped_bco
      << " wf_dropped_pool: " << dropped_pool
      << " ratio_bco: " << double(dropped_bco)/counts
      << " ratio_pool: " << double(dropped_pool)/counts
      << std::endl;
  }
  std::cout << std::endl;

  // drop per fee statistics
  for( const auto& [fee,counts]:m_fee_waveform_count_total )
  {
    const auto dropped_bco =  m_fee_waveform_count_dropped_bco[fee];
    std::cout << "SingleMicromegasPoolInput_v1::~SingleMicromegasPoolInput_v1 -"
      << " fee: " << fee
      << " wf_total: " << counts
      << " wf_dropped_bco: " << dropped_bco
      << " ratio_bco: " << double(dropped_bco)/counts
      << std::endl;

  }
  std::cout << std::endl;

  // also printout adjusted multipliers
  for( const auto& [packet_id, bco_matching]:m_bco_matching_information_map )
  {
    const auto old_precision = std::cout.precision();
    std::cout << "SingleMicromegasPoolInput_v1::~SingleMicromegasPoolInput_v1 -"
      << " packet: " << packet_id
      << " multiplier: " <<  std::setprecision(9) << bco_matching.get_adjusted_multiplier() << std::setprecision(old_precision)
      << std::endl;
  }

}

//______________________________________________________________
void SingleMicromegasPoolInput_v1::FillPool(const unsigned int /*nbclks*/)
{
  if (AllDone())  // no more files and all events read
  {
    return;
  }

  while (!GetEventiterator())  // at startup this is a null pointer
  {
    if (OpenNextFile() == InputFileHandlerReturnCodes::FAILURE)
    {
      AllDone(1);
      return;
    }
  }

  while (GetSomeMoreEvents())
  {
    std::unique_ptr<Event> evt(GetEventiterator()->getNextEvent());
    while (!evt)
    {
      fileclose();
      if (OpenNextFile() == InputFileHandlerReturnCodes::FAILURE)
      {
        AllDone(1);
        return;
      }

      // get next event
      evt.reset(GetEventiterator()->getNextEvent());
    }

    if (Verbosity() > 2)
    {
      std::cout << "Fetching next Event" << evt->getEvtSequence() << std::endl;
    }

    if (evt->getEvtType() != DATAEVENT)
    {
      m_NumSpecialEvents++;
      continue;
    }

    RunNumber(evt->getRunNumber());
    if (Verbosity() > 1)
    {
      evt->identify();
    }

    const int EventSequence = evt->getEvtSequence();
    const int npackets = evt->getPacketList(&plist[0], 10);

    if (npackets == 10)
    {
      std::cout << "SingleMicromegasPoolInput_v1::FillPool - too many packets" << std::endl;
      exit(1);
    }

    for (int i = 0; i < npackets; i++)
    {
      // keep pointer to local packet
      std::unique_ptr<Packet> packet(plist[i]);

      // get packet id
      const auto packet_id = packet->getIdentifier();

      if (Verbosity() > 1)
      {
        packet->identify();
      }

      // get relevant bco matching information
      auto& bco_matching_information = m_bco_matching_information_map[packet_id];
      bco_matching_information.set_verbosity(Verbosity());

      // read gtm bco information
      bco_matching_information.save_gtm_bco_information(packet.get());

      // save BCO from tagger internally
      const int n_tagger = packet->lValue(0, "N_TAGGER");
      for (int t = 0; t < n_tagger; ++t)
      {
        const bool is_lvl1 = static_cast<uint8_t>(packet->lValue(t, "IS_LEVEL1_TRIGGER"));
        const bool is_endat = static_cast<uint8_t>(packet->lValue(t, "IS_ENDAT"));
        if (is_lvl1||is_endat)
        {
          const uint64_t gtm_bco = static_cast<uint64_t>(packet->lValue(t, "BCO"));
          m_BeamClockPacket[gtm_bco].insert(packet_id);
          m_BclkStack.insert(gtm_bco);
        }
      }

      // loop over waveforms
      const int nwf = packet->iValue(0, "NR_WF");

      // increment counter and histogram
      m_waveform_count_total[packet_id] += nwf;
      h_waveform_count_total->Fill( std::to_string(packet_id).c_str(), nwf );

      if (Verbosity())
      {
        std::cout
          << "SingleMicromegasPoolInput_v1::FillPool -"
          << " packet: " << packet_id
          << " n_waveform: " << nwf
          << std::endl;
        bco_matching_information.print_gtm_bco_information();
      }

      // find reference from modebits, using BX_COUNTER_SYNC_T
      /*
       * This needs to be done even if the bco matching information is already verified
       * because any BX_COUNTER_SYNC_T event will break past references
       */
      bco_matching_information.find_reference_from_modebits(packet.get());

      // if bco matching information is not verified, try find reference from data
      if (!bco_matching_information.is_verified())
      {
        bco_matching_information.find_reference_from_data(packet.get());
      }

      // if bco matching information is still not verified, drop the packet
      if (!bco_matching_information.is_verified())
      {
        std::cout << "SingleMicromegasPoolInput_v1::FillPool - bco_matching not verified, dropping packet" << std::endl;
        m_waveform_count_dropped_bco[packet_id] += nwf;
        h_waveform_count_dropped_bco->Fill( std::to_string(packet_id).c_str(), nwf );
        continue;
      }

      for (int wf = 0; wf < nwf; wf++)
      {
        // get fee id
        const int fee_id = packet->iValue(wf, "FEE");
        ++m_fee_waveform_count_total[fee_id];

        // get checksum_error and check
        const auto checksum_error = packet->iValue(wf, "CHECKSUMERROR");
        if (checksum_error)
        {
          continue;
        }

        // get fee bco
        const unsigned int fee_bco = packet->iValue(wf, "BCO");

        // find matching gtm bco
        uint64_t gtm_bco = 0;
        const auto result = bco_matching_information.find_gtm_bco(fee_bco);

        if (result)
        {
          // assign gtm bco
          gtm_bco = result.value();
        }
        else
        {
          // increment counter and histogram
          ++m_waveform_count_dropped_bco[packet_id];
          ++m_fee_waveform_count_dropped_bco[fee_id];
          h_waveform_count_dropped_bco->Fill( std::to_string(packet_id).c_str(), 1 );

          // skip the waverform
          continue;
        }

        // get type
        // ignore heartbeat waveforms
        if( packet->iValue(wf, "TYPE" ) == HEARTBEAT_T ) continue;

        // get number of samples and check
        const uint16_t samples = packet->iValue(wf, "SAMPLES");
        if (samples < m_min_req_samples)
        {
          continue;
        }

        // create new hit
        auto newhit = std::make_unique<MicromegasRawHitv3>();
        newhit->set_bco(fee_bco);
        newhit->set_gtm_bco(gtm_bco);

        // packet id, fee id, channel, etc.
        newhit->set_packetid(packet_id);
        newhit->set_fee(fee_id);
        newhit->set_channel(packet->iValue(wf, "CHANNEL"));
        newhit->set_sampaaddress(packet->iValue(wf, "SAMPAADDRESS"));
        newhit->set_sampachannel(packet->iValue(wf, "CHANNEL"));

        // adc values
        for (uint16_t is = 0; is < samples; ++is)
        {
          uint16_t adc = packet->iValue(wf, is);
          if( adc == m_adc_invalid)
          {
            continue;
          }

          uint16_t first = is;
          MicromegasRawHitv3::adc_list_t values;
          for( ;is<samples && (adc = packet->iValue(wf, is)) != m_adc_invalid; ++is )
          { values.push_back(adc); }
          newhit->move_adc_waveform( first, std::move(values));

        }

        m_BeamClockFEE[gtm_bco].insert(fee_id);
        m_FEEBclkMap[fee_id] = gtm_bco;
        if (Verbosity() > 2)
        {
          std::cout << "evtno: " << EventSequence
                    << ", hits: " << wf
                    << ", num waveforms: " << nwf
                    << ", bco: 0x" << std::hex << gtm_bco << std::dec
                    << ", FEE: " << fee_id << std::endl;
        }

        if (StreamingInputManager())
        {
          StreamingInputManager()->AddMicromegasRawHit(gtm_bco, newhit.get());
        }

        m_MicromegasRawHitMap[gtm_bco].push_back(newhit.release());
      }
    }
  }
}

//______________________________________________________________
void SingleMicromegasPoolInput_v1::Print(const std::string& what) const
{
  if (what == "ALL" || what == "FEE")
  {
    for (const auto& bcliter : m_BeamClockFEE)
    {
      std::cout << "Beam clock 0x" << std::hex << bcliter.first << std::dec << std::endl;
      for (auto feeiter : bcliter.second)
      {
        std::cout << "FEM: " << feeiter << std::endl;
      }
    }
  }

  if (what == "ALL" || what == "FEEBCLK")
  {
    for (auto bcliter : m_FEEBclkMap)
    {
      std::cout << "FEE" << bcliter.first << " bclk: 0x"
                << std::hex << bcliter.second << std::dec << std::endl;
    }
  }

  if (what == "ALL" || what == "STORAGE")
  {
    for (const auto& bcliter : m_MicromegasRawHitMap)
    {
      std::cout << "Beam clock 0x" << std::hex << bcliter.first << std::dec << std::endl;
      for (auto feeiter : bcliter.second)
      {
        std::cout
            << "fee: " << feeiter->get_fee()
            << " at " << std::hex << feeiter << std::dec
            << std::endl;
      }
    }
  }

  if (what == "ALL" || what == "STACK")
  {
    for (const auto& iter : m_BclkStack)
    {
      std::cout << "stacked bclk: 0x" << std::hex << iter << std::dec << std::endl;
    }
  }
}

//____________________________________________________________________________
void SingleMicromegasPoolInput_v1::CleanupUsedPackets(const uint64_t bclk, bool dropped)
{

  // delete all raw hits associated to bco smaller than reference, and remove from map
  for(auto iter = m_MicromegasRawHitMap.begin(); iter != m_MicromegasRawHitMap.end() && (iter->first <= bclk); iter = m_MicromegasRawHitMap.erase(iter))
  {
    for (const auto& rawhit : iter->second)
    {
      if( dropped )
      {
        // increment dropped waveform counter and histogram
        ++m_waveform_count_dropped_pool[rawhit->get_packetid()];
        h_waveform_count_dropped_pool->Fill( std::to_string(rawhit->get_packetid()).c_str(), 1 );
      }
      delete rawhit;
    }
  }

  // cleanup bco stacks
  /* it erases all elements for which the bco is no greater than the provided one */
  m_BclkStack.erase(m_BclkStack.begin(), m_BclkStack.upper_bound(bclk));
  m_BeamClockFEE.erase(m_BeamClockFEE.begin(), m_BeamClockFEE.upper_bound(bclk));
  m_BeamClockPacket.erase(m_BeamClockPacket.begin(), m_BeamClockPacket.upper_bound(bclk));

  // cleanup matching information
  for( auto&& bco_matching:m_bco_matching_information_map )
  { bco_matching.second.cleanup(bclk); }

}

//_______________________________________________________
void SingleMicromegasPoolInput_v1::ClearCurrentEvent()
{
  std::cout << "SingleMicromegasPoolInput_v1::ClearCurrentEvent." << std::endl;
  uint64_t currentbclk = *m_BclkStack.begin();
  CleanupUsedPackets(currentbclk);
  return;
}

//_______________________________________________________
bool SingleMicromegasPoolInput_v1::GetSomeMoreEvents()
{
  if (AllDone())
  {
    return false;
  }

  // check minimum pool size
  if (m_MicromegasRawHitMap.size() < m_BcoPoolSize)
  {
    return true;
  }

  // make sure that the latest BCO received by each FEEs is past the current BCO
  std::set<int> toerase;
  uint64_t lowest_bclk = m_MicromegasRawHitMap.begin()->first + m_BcoRange;
  for (auto bcliter : m_FEEBclkMap)
  {
    if (bcliter.second <= lowest_bclk)
    {
      uint64_t highest_bclk = m_MicromegasRawHitMap.rbegin()->first;
      if ((highest_bclk - m_MicromegasRawHitMap.begin()->first) < MaxBclkDiff())
      {
        return true;
      }
      else
      {
        std::cout << PHWHERE << Name() << ": erasing FEE " << bcliter.first
                  << " with stuck bclk: " << std::hex << bcliter.second
                  << " current bco range: 0x" << m_MicromegasRawHitMap.begin()->first
                  << ", to: 0x" << highest_bclk << ", delta: " << std::dec
                  << (highest_bclk - m_MicromegasRawHitMap.begin()->first)
                  << std::dec << std::endl;
        toerase.insert(bcliter.first);
      }
    }
  }
  for(auto& fee: toerase)
  {
    m_FEEBclkMap.erase(fee);
  }

  return false;
}

//_______________________________________________________
void SingleMicromegasPoolInput_v1::CreateDSTNode(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  auto dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }

  PHNodeIterator iterDst(dstNode);
  auto detNode = dynamic_cast<PHCompositeNode*>(iterDst.findFirst("PHCompositeNode", "MICROMEGAS"));
  if (!detNode)
  {
    detNode = new PHCompositeNode("MICROMEGAS");
    dstNode->addNode(detNode);
  }

  auto container = findNode::getClass<MicromegasRawHitContainer>(detNode, m_rawHitContainerName);
  if (!container)
  {
    container = new MicromegasRawHitContainerv3();
    auto newNode = new PHIODataNode<PHObject>(container, m_rawHitContainerName, "PHObject");
    detNode->addNode(newNode);
  }
}

//_______________________________________________________
void SingleMicromegasPoolInput_v1::ConfigureStreamingInputManager()
{
  if (StreamingInputManager())
  {
    StreamingInputManager()->SetMicromegasBcoRange(m_BcoRange);
    StreamingInputManager()->SetMicromegasNegativeBco(m_NegativeBco);
  }
}

//_______________________________________________________
void SingleMicromegasPoolInput_v1::FillBcoQA(uint64_t gtm_bco)
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // set of packets for which the gtm BCO is found at least once
  unsigned int n_waveforms = 0;
  std::set<int> found_packets;
  for (uint64_t gtm_bco_loc = gtm_bco - m_NegativeBco; gtm_bco_loc < gtm_bco + m_BcoRange - m_NegativeBco; ++gtm_bco_loc)
  {
    // packet ids
    const auto packet_iter = m_BeamClockPacket.find(gtm_bco_loc);
    if( packet_iter != m_BeamClockPacket.end() )
    { found_packets.insert( packet_iter->second.begin(), packet_iter->second.end() ); }

    // waveforms
    const auto wf_iter = m_MicromegasRawHitMap.find(gtm_bco_loc);
    if (wf_iter != m_MicromegasRawHitMap.end())
    { n_waveforms += wf_iter->second.size(); }
  }

  // per packet statistics
  h_packet_stat->Fill( "Reference", 1 );

  for( const auto& packet_id:found_packets )
  { h_packet_stat->Fill( std::to_string(packet_id).c_str(), 1 ); }
  h_packet_stat->Fill( "All", found_packets.size()>= m_npackets_active );

  // how many packet_id found for this BCO
  h_packet->Fill(found_packets.size());

  // how many waveforms found for this BCO
  h_waveform->Fill(n_waveforms);
}
//_______________________________________________________
void SingleMicromegasPoolInput_v1::createQAHistos()
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // number of packets found with BCO from felix matching reference BCO
  h_packet = new TH1I( "h_MicromegasBCOQA_npacket_bco", "TPOT Packet Count per GTM BCO; Matching BCO tagger count; GL1 trigger count", 10, 0, 10 );

  // number of waveforms found with BCO from felix matching reference BCO
  h_waveform = new TH1I( "h_MicromegasBCOQA_nwaveform_bco", "TPOT Waveform Count per GTM BCO; Matching Waveform count; GL1 trigger count", 4100, 0, 4100 );

  /*
   * first bin is the number of requested GL1 BCO, for reference
   * next two bins is the number of times the GL1 BCO is found in the taggers list for a given packet_id
   * last bin is the sum
   */
  h_packet_stat = new TH1I( "h_MicromegasBCOQA_packet_stat", "Matching Tagger count per packet; packet id; GL1 trigger count", m_npackets_active+2, 0, m_npackets_active+2 );
  h_packet_stat->GetXaxis()->SetBinLabel(1, "Reference" );
  h_packet_stat->GetXaxis()->SetBinLabel(2, "5001" );
  h_packet_stat->GetXaxis()->SetBinLabel(3, "5002" );
  h_packet_stat->GetXaxis()->SetBinLabel(4, "All" );
  h_packet_stat->GetYaxis()->SetTitle( "trigger count" );

  // total number of waveform per packet
  h_waveform_count_total = new TH1I( "h_MicromegasBCOQA_waveform_count_total", "Total number of waveforms per packet", m_npackets_active, 0, m_npackets_active );

  // number of dropped waveform per packet due to bco mismatch
  h_waveform_count_dropped_bco = new TH1I( "h_MicromegasBCOQA_waveform_count_dropped_bco", "Number of dropped waveforms per packet (bco)", m_npackets_active, 0, m_npackets_active );

  // number of dropped waveform per packet due to fun4all pool mismatch
  h_waveform_count_dropped_pool = new TH1I( "h_MicromegasBCOQA_waveform_count_dropped_pool", "Number of dropped waveforms per packet (pool)", m_npackets_active, 0, m_npackets_active );

  // define axis
  for( const auto& h:std::initializer_list<TH1*>{h_waveform_count_total, h_waveform_count_dropped_bco, h_waveform_count_dropped_pool} )
  {
    h->GetXaxis()->SetBinLabel(1, "5001" );
    h->GetXaxis()->SetBinLabel(2, "5002" );
    h->GetYaxis()->SetTitle( "waveform count" );
    h->GetXaxis()->SetTitle( "packet id" );
  }

  // register all histograms to histogram manager
  for( const auto& h:std::initializer_list<TH1*>{h_packet, h_waveform, h_packet_stat, h_waveform_count_total, h_waveform_count_dropped_bco, h_waveform_count_dropped_pool} )
  {
    h->SetFillStyle(1001);
    h->SetFillColor(kYellow);
    hm->registerHisto(h);
  }
}
