#include "SingleMicromegasPoolInput.h"

#include "Fun4AllStreamingInputManager.h"

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

#include <memory>
#include <set>

namespace
{
  // streamer for lists
  template <class T>
  std::ostream& operator<<(std::ostream& out, const std::list<T>& list)
  {
    if (list.empty())
    {
      out << "{}";
    }
    else
    {
      out << "{ ";
      bool first = true;
      for (const auto& value : list)
      {
        if (!first)
        {
          out << ", ";
        }
        out << value;
        first = false;
      }
      out << " }";
    }
    return out;
  }

  // minimum number of requested samples
  static constexpr int m_min_req_samples = 5;

}  // namespace

//______________________________________________________________
SingleMicromegasPoolInput::SingleMicromegasPoolInput(const std::string& name)
  : SingleStreamingInput(name)
{
  SubsystemEnum(Fun4AllStreamingInputManager::MICROMEGAS);
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
      auto& packet = plist[i];

      // get packet id
      const auto packet_id = packet->getIdentifier();

      if (Verbosity() > 1)
      {
        packet->identify();
      }

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

          // store as available bco for all FEE's bco alignment objects
          for (auto&& bco_alignment : m_bco_alignment_list)
          {
            bco_alignment.gtm_bco_list.push_back(gtm_bco);
          }
        }
      }

      // loop over waveforms
      const int nwf = packet->iValue(0, "NR_WF");

      if (Verbosity())
      {
        std::cout << "SingleMicromegasPoolInput::FillPool -"
                  << " packet: " << packet_id
                  << " n_lvl1_bco: " << m_bco_alignment_list[0].gtm_bco_list.size()
                  << " n_waveform: " << nwf
                  << std::endl;

        std::cout << "SingleMicromegasPoolInput::FillPool -"
                  << " packet: " << packet_id
                  << " bco: " << std::hex << m_bco_alignment_list[0].gtm_bco_list << std::dec
                  << std::endl;
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

        // get matching gtm bco, from bco alignment object
        uint64_t gtm_bco = 0;
        auto&& bco_alignment = m_bco_alignment_list.at(fee_id);
        if (bco_alignment.fee_bco == fee_bco)
        {
          // assign gtm bco
          gtm_bco = bco_alignment.gtm_bco;
        }
        else if (!bco_alignment.gtm_bco_list.empty())
        {
          // get first available gtm_bco from list
          gtm_bco = bco_alignment.gtm_bco_list.front();

          if (Verbosity())
          {
            std::cout << "SingleMicromegasPoolInput::FillPool -"
                      << " fee_id: " << fee_id
                      << " fee_bco: 0x" << std::hex << fee_bco
                      << " gtm_bco: 0x" << gtm_bco
                      << std::dec
                      << std::endl;
          }

          // store as current gtm/fee bco
          bco_alignment.fee_bco = fee_bco;
          bco_alignment.gtm_bco = gtm_bco;

          // remove from list of availables gtm_bco
          bco_alignment.gtm_bco_list.pop_front();
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

      // clear all stored bco list
      for (auto&& bco_alignment : m_bco_alignment_list)
      {
        bco_alignment.gtm_bco_list.clear();
      }

      delete packet;
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
      return true;
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
