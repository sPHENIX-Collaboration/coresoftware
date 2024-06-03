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
  // minimum number of requested samples
  static constexpr int m_min_req_samples = 5;
}  // namespace

//______________________________________________________________
SingleMicromegasPoolInput::SingleMicromegasPoolInput(const std::string& name)
  : SingleStreamingInput(name)
{
  SubsystemEnum(InputManagerType::MICROMEGAS);
}

//______________________________________________________________
SingleMicromegasPoolInput::~SingleMicromegasPoolInput()
{
  if( Verbosity() )
  {
    std::cout << "SingleMicromegasPoolInput::~SingleMicromegasPoolInput - waveform_count_total: " << m_waveform_count_total << std::endl;
    std::cout << "SingleMicromegasPoolInput::~SingleMicromegasPoolInput - waveform_count_dropped: " << m_waveform_count_dropped << std::endl;
    std::cout << "SingleMicromegasPoolInput::~SingleMicromegasPoolInput - ratio: " << double(m_waveform_count_dropped)/m_waveform_count_total << std::endl;
  }
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
    const int npackets = evt->getPacketList(&plist[0], 10);

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
      bco_matching_information.set_verbosity( Verbosity() );

      // read gtm bco information
      bco_matching_information.save_gtm_bco_information( packet.get() );

      // loop over waveforms
      const int nwf = packet->iValue(0, "NR_WF");
      m_waveform_count_total += nwf;

      if (Verbosity())
      {
        std::cout
          << "SingleMicromegasPoolInput::FillPool -"
          << " packet: " << packet_id
          << " n_waveform: " << nwf
          << std::endl;
      }

      // try find reference
      if( !bco_matching_information.is_verified() )
      { bco_matching_information.find_reference( packet.get() ); }

      // drop packet if not found
      if( !bco_matching_information.is_verified() )
      {
        std::cout << "SingleMicromegasPoolInput::FillPool - bco_matching not verified, dropping packet" << std::endl;
        m_waveform_count_dropped += nwf;
        continue;
      }

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

        // find matching gtm bco
        uint64_t gtm_bco = 0;
        const auto result = bco_matching_information.find_gtm_bco( fee_bco );

        if( result )
        {
          // assign gtm bco
          gtm_bco = result.value();
        }
        else
        {
          // increment count
          ++m_waveform_count_dropped;

          // skip the waverform
          continue;
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
      bco_matching_information.cleanup();
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
