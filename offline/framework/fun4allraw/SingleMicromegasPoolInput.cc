#include "SingleMicromegasPoolInput.h"

#include "Fun4AllStreamingInputManager.h"
#include "InputManagerType.h"

#include <ffarawobjects/MicromegasRawHitContainerv1.h>
#include <ffarawobjects/MicromegasRawHitv1.h>

#include <frog/FROG.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/getClass.h>
#include <phool/phool.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/Eventiterator.h>
#include <Event/fileEventiterator.h>

#include <algorithm>
#include <memory>
#include <set>

namespace
{
  // streamer for lists
  template <class T>
  std::ostream& operator<<(std::ostream& o, const std::list<T>& list)
  {
    if (list.empty())
    {
      o << "{}";
    }
    else
    {
      const bool is_hex = (o.flags()&std::ios_base::hex);
      o << "{ ";
      bool first = true;
      for (const auto& value : list)
      {
        if (!first)
        {
          o << ", ";
        }
        if( is_hex )
	{
	  o << "0x";
	}
        o << value;
        first = false;
      }
      o << " }";
    }
    return o;
  }

  // get the difference of two unsiged numbers
  template<class T>
  inline static constexpr T get_bco_diff( const T& first, const T& second )
  { return first < second ? (second-first):(first-second); }

  // minimum number of requested samples
  static constexpr int m_min_req_samples = 5;

  // define limit for matching two fee_bco
  static constexpr unsigned int max_fee_bco_diff = 10;
  static constexpr unsigned int max_gtm_bco_diff = 100;

  // needed to avoid memory leak. Assumes that we will not be assembling more than 50 events at the same time
  static constexpr unsigned int max_matching_data_size = 50;

}  // namespace

//_________________________________________________________
void SingleMicromegasPoolInput::bco_matching_information_t::truncate( unsigned int maxsize )
{
  while( m_gtm_bco_list.size() > maxsize ) { m_gtm_bco_list.pop_front(); }
  while( m_bco_matching_list.size() > maxsize ) { m_bco_matching_list.pop_front(); }
}

//_________________________________________________________
unsigned int SingleMicromegasPoolInput::bco_matching_information_t::get_predicted_fee_bco( uint64_t gtm_bco ) const
{
  // check proper initialization
  if( !(m_has_gtm_bco_first && m_has_fee_bco_first ) ) { return 0; }

  // this is the clock multiplier from lvl1 to fee clock
  /* todo: should replace with actual rational number for John K. */
  static constexpr double multiplier = 4.2629164;

  // get lvl1 bco difference with proper rollover accounting
  uint64_t gtm_bco_difference = (gtm_bco >= m_gtm_bco_first) ?
    (gtm_bco - m_gtm_bco_first):
    (gtm_bco + (1ULL<<40U) - m_gtm_bco_first);

  // convert to fee bco, and truncate to 20 bits
  uint64_t fee_bco_predicted = m_fee_bco_first + multiplier*(gtm_bco_difference);
  return (unsigned int)(fee_bco_predicted & 0xFFFFFU);
}

//______________________________________________________________
SingleMicromegasPoolInput::SingleMicromegasPoolInput(const std::string& name)
  : SingleStreamingInput(name)
{
  SubsystemEnum(InputManagerType::MICROMEGAS);
  plist = new Packet*[10];
}

//______________________________________________________________
SingleMicromegasPoolInput::~SingleMicromegasPoolInput()
{
  delete[] plist;
}

//______________________________________________________________
void SingleMicromegasPoolInput::FillPool(const unsigned int /*nbclks*/)
{
  if (AllDone())  // no more files and all events read
  {
    return;
  }

  while (!GetEventiterator())  // at startup this is a null pointer
  {
    if (!OpenNextFile())
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
      if (!OpenNextFile())
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

    RunNumber(evt->getRunNumber());
    if (Verbosity() > 1)
    {
      evt->identify();
    }

    if (evt->getEvtType() != DATAEVENT)
    {
      m_NumSpecialEvents++;
      continue;
    }

    const int EventSequence = evt->getEvtSequence();
    const int npackets = evt->getPacketList(plist, 10);

    if (npackets == 10)
    {
      std::cout << "SingleMicromegasPoolInput::FillPool - too many packets" << std::endl;
      exit(1);
    }

    for (int i = 0; i < npackets; i++)
    {
      // keep pointer to local packet
      std::unique_ptr<Packet> packet( plist[i] );

      // get packet id
      const auto packet_id = packet->getIdentifier();

      if (Verbosity() > 1)
      {
        packet->identify();
      }

      // get relevant bco matching information
      auto& bco_matching_information = m_bco_matching_information_map[packet_id];

      // loop over taggers
      const int ntagger = packet->lValue(0, "N_TAGGER");
      for (int t = 0; t < ntagger; t++)
      {
        // only store gtm_bco for level1 type of taggers (not ENDDAT)
        const auto is_lvl1 = static_cast<uint8_t>(packet->lValue(t, "IS_LEVEL1_TRIGGER"));
        if (is_lvl1)
        {
          // get gtm_bco
          const auto gtm_bco = static_cast<uint64_t>(packet->lValue(t, "BCO"));

          // initialize first gtm_bco
          if( !bco_matching_information.m_has_gtm_bco_first )
          {
            bco_matching_information.m_gtm_bco_first = gtm_bco;
            bco_matching_information.m_has_gtm_bco_first = true;

            if( Verbosity() )
            {
              std::cout
                << "SingleMicromegasPoolInput::FillPool -"
                << " packet: " << packet_id
                << std::hex
                << " m_gtm_bco_first: 0x" << gtm_bco
                << std::dec
                << std::endl;
            }
          }

          // store in list of available bco
          bco_matching_information.m_gtm_bco_list.push_back(gtm_bco);
        }
      }

      // loop over waveforms
      const int nwf = packet->iValue(0, "NR_WF");

      if (Verbosity())
      {
        std::cout
          << "SingleMicromegasPoolInput::FillPool -"
          << " packet: " << packet_id
          << " n_gtm_bco: " << bco_matching_information.m_gtm_bco_list.size()
          << " n_waveform: " << nwf
          << std::endl;

        if (!bco_matching_information.m_gtm_bco_list.empty())
        {
          std::cout
            << "SingleMicromegasPoolInput::FillPool -"
            << " packet: " << packet_id
            << " gtm_bco: " << std::hex << bco_matching_information.m_gtm_bco_list << std::dec
            << std::endl;

          // also print predicted fee bco
          std::list<unsigned int> fee_bco_predicted_list;
          std::transform(
            bco_matching_information.m_gtm_bco_list.begin(),
            bco_matching_information.m_gtm_bco_list.end(),
            std::back_inserter(fee_bco_predicted_list),
            [&bco_matching_information](const uint64_t& gtm_bco ){ return bco_matching_information.get_predicted_fee_bco(gtm_bco); } );

          std::cout
            << "SingleMicromegasPoolInput::FillPool -"
            << " packet: " << packet_id
            << " fee_bco_predicted: " << std::hex << fee_bco_predicted_list << std::dec
            << std::endl;
        }
      }

      // keep track of orphans
      using fee_pair_t = std::pair<unsigned int, unsigned int>;
      std::set<fee_pair_t> orphans;

      for (int wf = 0; wf < nwf; wf++)
      {
        // get fee id
        const int fee_id = packet->iValue(wf, "FEE");

        // get checksum_error and check
        const auto checksum_error = packet->iValue(wf, "CHECKSUMERROR");
        if (checksum_error)
        {
          continue;
        }

        // get number of samples and check
        const uint16_t samples = packet->iValue(wf, "SAMPLES");
        if (samples < m_min_req_samples)
        {
          continue;
        }

        // get fee bco
        const unsigned int fee_bco = packet->iValue(wf, "BCO");

        // initialize first fee_bco
        if( !bco_matching_information.m_has_fee_bco_first )
        {
          bco_matching_information.m_fee_bco_first = fee_bco;
          bco_matching_information.m_has_fee_bco_first = true;

          if( Verbosity() )
          {
            std::cout
              << "SingleMicromegasPoolInput::FillPool -"
              << " packet: " << packet_id
              << std::hex
              << " m_fee_bco_first: 0x" << fee_bco
              << std::dec
              << std::endl;
            }
        }

        // find matching gtm bco
        uint64_t gtm_bco = 0;
        const auto bco_matching_iter = std::find_if(
          bco_matching_information.m_bco_matching_list.begin(),
          bco_matching_information.m_bco_matching_list.end(),
          [&fee_bco]( const m_bco_matching_pair_t& pair )
          { return get_bco_diff( pair.first, fee_bco ) < max_fee_bco_diff; } );

        if( bco_matching_iter != bco_matching_information.m_bco_matching_list.end() )
        {
          // assign gtm bco
          gtm_bco = bco_matching_iter->second;
        }
        else
        {

          // find gtm_bco in list that match fee bco
          auto iter = std::find_if(
            bco_matching_information.m_gtm_bco_list.begin(),
            bco_matching_information.m_gtm_bco_list.end(),
            [&fee_bco, &bco_matching_information]( const uint64_t& gtm_bco_local )
            { return get_bco_diff( bco_matching_information.get_predicted_fee_bco(gtm_bco_local), fee_bco ) < max_gtm_bco_diff; } );

          if( iter != bco_matching_information.m_gtm_bco_list.end() )
          {
            gtm_bco = *iter;
            if (Verbosity())
            {
              std::cout << "SingleMicromegasPoolInput::FillPool -"
                << " fee_id: " << fee_id
                << std::hex
                << " fee_bco: 0x" << fee_bco
                << " predicted: 0x" << bco_matching_information.get_predicted_fee_bco(gtm_bco)
                << " gtm_bco: 0x" << gtm_bco
                << std::dec
                << std::endl;
            }

            // fee_bco is new. Assume it corresponds to the first available gtm bco
            // update running fee_bco and gtm_bco pair accordingly
            bco_matching_information.m_bco_matching_list.emplace_back(fee_bco, gtm_bco);

            // remove bco from running list
            bco_matching_information.m_gtm_bco_list.erase(iter);
          }
          else
          {
            if (Verbosity() && orphans.insert(std::make_pair(fee_id, fee_bco)).second)
            {
              std::cout << "SingleMicromegasPoolInput::FillPool -"
                << " fee_id: " << fee_id
                << " fee_bco: 0x" << std::hex << fee_bco << std::dec
                << " gtm_bco: none"
                << std::endl;
            }

            // skip the waverform
            continue;
          }
        }

        // create new hit
        auto newhit = std::make_unique<MicromegasRawHitv1>();
        newhit->set_bco(fee_bco);
        newhit->set_gtm_bco(gtm_bco);

        // packet id, fee id, channel, etc.
        newhit->set_packetid(packet_id);
        newhit->set_fee(fee_id);
        newhit->set_channel(packet->iValue(wf, "CHANNEL"));
        newhit->set_sampaaddress(packet->iValue(wf, "SAMPAADDRESS"));
        newhit->set_sampachannel(packet->iValue(wf, "CHANNEL"));

        // assign samples
        newhit->set_samples(samples);

        // adc values
        for (uint16_t is = 0; is < samples; ++is)
        {
          newhit->set_adc(is, packet->iValue(wf, is));
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
        m_BclkStack.insert(gtm_bco);
      }

      // cleanup
      bco_matching_information.truncate(max_matching_data_size);
    }
  }
}

//______________________________________________________________
void SingleMicromegasPoolInput::Print(const std::string& what) const
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
void SingleMicromegasPoolInput::CleanupUsedPackets(const uint64_t bclk)
{
  std::vector<uint64_t> toclearbclk;
  for (const auto& iter : m_MicromegasRawHitMap)
  {
    if (iter.first <= bclk)
    {
      for (auto pktiter : iter.second)
      {
        delete pktiter;
      }

      toclearbclk.push_back(iter.first);
    }
    else
    {
      break;
    }
  }

  for (auto iter : toclearbclk)
  {
    m_BclkStack.erase(iter);
    m_BeamClockFEE.erase(iter);
    m_MicromegasRawHitMap.erase(iter);
  }
}

//_______________________________________________________
void SingleMicromegasPoolInput::ClearCurrentEvent()
{
  uint64_t currentbclk = *m_BclkStack.begin();
  CleanupUsedPackets(currentbclk);
  return;
}

//_______________________________________________________
bool SingleMicromegasPoolInput::GetSomeMoreEvents()
{
  if (AllDone())
  {
    return false;
  }
  if (m_MicromegasRawHitMap.empty())
  {
    return true;
  }

  uint64_t lowest_bclk = m_MicromegasRawHitMap.begin()->first + m_BcoRange;
  for (auto bcliter : m_FEEBclkMap)
  {
    if (bcliter.second <= lowest_bclk)
    {
      uint64_t highest_bclk = m_MicromegasRawHitMap.rbegin()->first;
      if ((highest_bclk - m_MicromegasRawHitMap.begin()->first) < MaxBclkDiff())
      {
        // std::cout << "FEE " << bcliter.first << " bclk: "
        // 		<< std::hex << bcliter.second << ", req: " << lowest_bclk
        // 		 << " low: 0x" <<  m_MicromegasRawHitMap.begin()->first << ", high: " << highest_bclk << ", delta: " << std::dec << (highest_bclk-m_MicromegasRawHitMap.begin()->first)
        // 		<< std::dec << std::endl;
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
        m_FEEBclkMap.erase(bcliter.first);
      }
    }
  }

  return false;
}

//_______________________________________________________
void SingleMicromegasPoolInput::CreateDSTNode(PHCompositeNode* topNode)
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

  auto container = findNode::getClass<MicromegasRawHitContainer>(detNode, "MICROMEGASRAWHIT");
  if (!container)
  {
    container = new MicromegasRawHitContainerv1();
    auto newNode = new PHIODataNode<PHObject>(container, "MICROMEGASRAWHIT", "PHObject");
    detNode->addNode(newNode);
  }
}

//_______________________________________________________
void SingleMicromegasPoolInput::ConfigureStreamingInputManager()
{
  if (StreamingInputManager())
  {
    StreamingInputManager()->SetMicromegasBcoRange(m_BcoRange);
    StreamingInputManager()->SetMicromegasNegativeBco(m_NegativeBco);
  }
}
