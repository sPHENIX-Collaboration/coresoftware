#include "SingleTpcPoolInput.h"

#include "Fun4AllStreamingInputManager.h"
#include "InputManagerType.h"

#include <ffarawobjects/TpcRawHitContainerv2.h>
#include <ffarawobjects/TpcRawHitv2.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/getClass.h>
#include <phool/phool.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/Eventiterator.h>

#include <memory>
#include <set>

SingleTpcPoolInput::SingleTpcPoolInput(const std::string &name)
  : SingleStreamingInput(name)
{
  SubsystemEnum(InputManagerType::TPC);
  m_rawHitContainerName = "TPCRAWHIT";
}

void SingleTpcPoolInput::FillPool(const uint64_t minBCO)
{
  if (AllDone())  // no more files and all events read
  {
    return;
  }
  while (GetEventiterator() == nullptr)  // at startup this is a null pointer
  {
    if (!OpenNextFile())
    {
      AllDone(1);
      return;
    }
  }
  //  std::set<uint64_t> saved_beamclocks;
  while (GetSomeMoreEvents(0))
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
      evt.reset(GetEventiterator()->getNextEvent());
    }
    if (Verbosity() > 2)
    {
      std::cout << "Fetching next Event" << evt->getEvtSequence() << std::endl;
    }
    RunNumber(evt->getRunNumber());
    if (GetVerbosity() > 1)
    {
      evt->identify();
    }
    if (evt->getEvtType() != DATAEVENT)
    {
      m_NumSpecialEvents++;
      if (evt->getEvtType() == ENDRUNEVENT)
      {
        std::cout << "End run flag for " << Name() << " found, remaining TPC data is corrupted" << std::endl;
        AllDone(1);
        return;
      }
      continue;
    }
    int EventSequence = evt->getEvtSequence();
    std::vector<Packet *> pktvec = evt->getPacketVector();
    if (m_skipEarlyEvents)
    {
      for (auto packet : pktvec)
      {
        int numBCOs = packet->lValue(0, "N_TAGGER");
        bool bco_set = false;
        for (int j = 0; j < numBCOs; j++)
        {
          const auto is_lvl1 = static_cast<uint8_t>(packet->lValue(j, "IS_LEVEL1_TRIGGER"));
          const auto is_endat = static_cast<uint8_t>(packet->lValue(j, "IS_ENDAT"));
          if ((is_lvl1 || is_endat) && !bco_set)
          {
            bco_set = true;
            uint64_t bco = packet->lValue(j, "BCO");
            if (bco < minBCO)
            {
              continue;
            }
            m_skipEarlyEvents = false;
          }
        }
      }
    }
    if (m_skipEarlyEvents)
    {
      for (auto packet : pktvec)
      {
        delete packet;
      }
      continue;
    }
    for (auto packet : pktvec)
    {
      // get packet id
      const auto packet_id = packet->getIdentifier();

      if (Verbosity() > 1)
      {
        packet->identify();
      }

      // by default use previous bco clock for gtm bco
      auto &previous_bco = m_packet_bco[packet_id];
      uint64_t gtm_bco = previous_bco;

      uint64_t m_nTaggerInFrame = packet->lValue(0, "N_TAGGER");
      bool skipthis = true;
      uint64_t largest_bco = 0;
      bool bco_set = false;
      for (uint64_t t = 0; t < m_nTaggerInFrame; t++)
      {
        // only store gtm_bco for level1 type of taggers (not ENDDAT)
        const auto is_lvl1 = static_cast<uint8_t>(packet->lValue(t, "IS_LEVEL1_TRIGGER"));

        const auto is_endat = static_cast<uint8_t>(packet->lValue(t, "IS_ENDAT"));
        if ((is_lvl1 || is_endat) && !bco_set)
        {
          bco_set = true;
          gtm_bco = packet->lValue(t, "BCO");
          if (largest_bco < gtm_bco)
          {
            largest_bco = gtm_bco;
          }
          if (gtm_bco < minBCO)
          {
            continue;
          }
          if (Verbosity() > 0)
          {
            std::cout << "bco: 0x" << std::hex << gtm_bco << std::dec << std::endl;
          }
          // store
          skipthis = false;
          previous_bco = gtm_bco;
          if (m_BclkStackPacketMap.find(packet_id) == m_BclkStackPacketMap.end())
          {
            m_BclkStackPacketMap.insert(std::make_pair(packet_id, std::set<uint64_t>()));
          }
          m_BclkStackPacketMap[packet_id].insert(gtm_bco);
        }
      }
      if (skipthis)
      {
        if (Verbosity() > 1)
        {
          std::cout << "Largest bco: 0x" << std::hex << largest_bco << ", minbco 0x"
                    << minBCO << std::dec << ", evtno: " << EventSequence << std::endl;
        }
      }
      else
      {
        int m_nWaveFormInFrame = packet->iValue(0, "NR_WF");
        static int once = 0;
        for (int wf = 0; wf < m_nWaveFormInFrame; wf++)
        {
          if (m_TpcRawHitMap[gtm_bco].size() > 20000)
          {
            if (!once)
            {
              std::cout << "too many hits" << std::endl;
            }
            once++;
            continue;
          }
          else
          {
            if (once)
            {
              std::cout << "many more hits: " << once << std::endl;
            }
            once = 0;
          }
          bool checksumerror = (packet->iValue(wf, "CHECKSUMERROR") > 0);
          if (checksumerror)
          {
            continue;
          }
          bool parityerror = (packet->iValue(wf, "DATAPARITYERROR") > 0);
          auto newhit = std::make_unique<TpcRawHitv2>();
          int FEE = packet->iValue(wf, "FEE");
          newhit->set_bco(packet->iValue(wf, "BCO"));

          // store gtm bco in hit
          newhit->set_gtm_bco(gtm_bco);

          newhit->set_packetid(packet->getIdentifier());
          newhit->set_fee(FEE);
          newhit->set_channel(packet->iValue(wf, "CHANNEL"));
          newhit->set_sampaaddress(packet->iValue(wf, "SAMPAADDRESS"));
          newhit->set_sampachannel(packet->iValue(wf, "SAMPACHANNEL"));
          newhit->set_type(packet->iValue(wf, "TYPE"));
          newhit->set_userword(packet->iValue(wf, "USERWORD"));
          newhit->set_parity(packet->iValue(wf, "DATAPARITY"));
          newhit->set_checksum(packet->iValue(wf, "CHECKSUM"));
          newhit->set_checksumerror(checksumerror);
          newhit->set_parityerror(parityerror);
          //         // checksum and checksum error
          //         newhit->set_checksum( packet->iValue(iwf, "CHECKSUM") );
          //         newhit->set_checksum_error( packet->iValue(iwf, "CHECKSUMERROR") );

          // samples
          // const uint16_t samples = packet->iValue(wf, "SAMPLES");

          // Temp remedy as we set the time window as 425 for now (extended from previous 360
          // due to including of diffused laser flush)
          const uint16_t samples = m_max_tpc_time_samples;

          newhit->set_samples(samples);

          // adc values
          for (uint16_t is = 0; is < samples; ++is)
          {
            uint16_t adval = packet->iValue(wf, is);

            // This is temporary fix for decoder change. Will be changed again for real ZS data decoding.
            // if(adval >= 64000){ newhit->set_samples(is); break;}

            // With this, the hit is unseen from clusterizer
            if (adval <= 64000)
            {
              newhit->set_adc(is, adval);
            }
          }

          m_BeamClockFEE[gtm_bco].insert(FEE);
          m_FEEBclkMap[FEE] = gtm_bco;
          if (Verbosity() > 2)
          {
            std::cout << "evtno: " << EventSequence
                      << ", hits: " << wf
                      << ", num waveforms: " << m_nWaveFormInFrame
                      << ", bco: 0x" << std::hex << gtm_bco << std::dec
                      << ", FEE: " << FEE << std::endl;
          }
          //          packet->convert();
          // if (m_TpcRawHitMap[gtm_bco].size() < 50000)
          // {
          if (StreamingInputManager())
          {
            StreamingInputManager()->AddTpcRawHit(gtm_bco, newhit.get());
          }
          m_TpcRawHitMap[gtm_bco].push_back(newhit.release());
          m_BclkStack.insert(gtm_bco);
          //	}
        }
      }
      delete packet;
    }
  }

  //    Print("HITS");
  //  } while (m_TpcRawHitMap.size() < 10 || CheckPoolDepth(m_TpcRawHitMap.begin()->first));
}

void SingleTpcPoolInput::Print(const std::string &what) const
{
  if (what == "ALL" || what == "FEE")
  {
    for (const auto &bcliter : m_BeamClockFEE)
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
  if (what == "ALL" || what == "HITS")
  {
    const auto bcliter = m_TpcRawHitMap.begin();
    {
      std::cout << Name() << ": Beam clock 0x" << std::hex << bcliter->first
                << std::dec << ", Number of hits: " << bcliter->second.size()
                << std::endl;
    }
  }
  if (what == "ALL" || what == "STORAGE")
  {
    for (const auto &bcliter : m_TpcRawHitMap)
    {
      std::cout << "Beam clock 0x" << std::hex << bcliter.first << std::dec << std::endl;
      for (auto feeiter : bcliter.second)
      {
        std::cout << "fee: " << feeiter->get_fee()
                  << " at " << std::hex << feeiter << std::dec << std::endl;
      }
    }
  }
  if (what == "ALL" || what == "STACK")
  {
    for (auto iter : m_BclkStack)
    {
      std::cout << "stacked bclk: 0x" << std::hex << iter << std::dec << std::endl;
    }
  }
}

void SingleTpcPoolInput::CleanupUsedPackets(const uint64_t bclk)
{
  if (Verbosity() > 2)
  {
    std::cout << "cleaning up bcos < 0x" << std::hex
              << bclk << std::dec << std::endl;
  }
  std::vector<uint64_t> toclearbclk;
  for (const auto &iter : m_TpcRawHitMap)
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
  // for (auto iter :  m_BeamClockFEE)
  // {
  //   iter.second.clear();
  // }
  for (auto iter : toclearbclk)
  {
    m_BclkStack.erase(iter);
    m_BeamClockFEE.erase(iter);
    m_TpcRawHitMap.erase(iter);
  }
}

bool SingleTpcPoolInput::CheckPoolDepth(const uint64_t bclk)
{
  // if (m_FEEBclkMap.size() < 10)
  // {
  //   std::cout << "not all FEEs in map: " << m_FEEBclkMap.size() << std::endl;
  //   return true;
  // }
  for (auto iter : m_FEEBclkMap)
  {
    if (Verbosity() > 2)
    {
      std::cout << "my bclk 0x" << std::hex << iter.second
                << " req: 0x" << bclk << std::dec << std::endl;
    }
    if (iter.second < bclk)
    {
      if (Verbosity() > 1)
      {
        std::cout << "FEE " << iter.first << " beamclock 0x" << std::hex << iter.second
                  << " smaller than req bclk: 0x" << bclk << std::dec << std::endl;
      }
      return false;
    }
  }
  return true;
}

void SingleTpcPoolInput::ClearCurrentEvent()
{
  // called interactively, to get rid of the current event
  uint64_t currentbclk = *m_BclkStack.begin();
  //  std::cout << "clearing bclk 0x" << std::hex << currentbclk << std::dec << std::endl;
  CleanupUsedPackets(currentbclk);
  // m_BclkStack.erase(currentbclk);
  // m_BeamClockFEE.erase(currentbclk);
  return;
}

bool SingleTpcPoolInput::GetSomeMoreEvents(const uint64_t ibclk)
{
  if (AllDone())
  {
    return false;
  }
  if (m_TpcRawHitMap.empty())
  {
    return true;
  }
  uint64_t localbclk = ibclk;
  if (ibclk == 0)
  {
    if (m_TpcRawHitMap.empty())
    {
      return true;
    }
    localbclk = m_TpcRawHitMap.begin()->first;
  }

  std::set<int> toerase;
  for (auto bcliter : m_FEEBclkMap)
  {
    if (bcliter.second <= localbclk)
    {
      uint64_t highest_bclk = m_TpcRawHitMap.rbegin()->first;
      if ((highest_bclk - m_TpcRawHitMap.begin()->first) < MaxBclkDiff())
      {
        // std::cout << "FEE " << bcliter.first << " bclk: "
        // 		<< std::hex << bcliter.second << ", req: " << lowest_bclk
        // 		 << " low: 0x" <<  m_TpcRawHitMap.begin()->first << ", high: " << highest_bclk << ", delta: " << std::dec << (highest_bclk-m_TpcRawHitMap.begin()->first)
        // 		<< std::dec << std::endl;
        return true;
      }
      else
      {
        std::cout << PHWHERE << Name() << ": erasing FEE " << bcliter.first
                  << " with stuck bclk: " << std::hex << bcliter.second
                  << " current bco range: 0x" << m_TpcRawHitMap.begin()->first
                  << ", to: 0x" << highest_bclk << ", delta: " << std::dec
                  << (highest_bclk - m_TpcRawHitMap.begin()->first)
                  << std::dec << std::endl;
        toerase.insert(bcliter.first);
      }
    }
  }
  for (auto iter : toerase)
  {
    m_FEEBclkMap.erase(iter);
  }
  return false;
  // if (CheckPoolDepth(m_TpcRawHitMap.begin()->first))
  // {
  //   if (m_TpcRawHitMap.size() >= 10)
  //   {
  //     return false;
  //   }
  // }
  // return true;
}

void SingleTpcPoolInput::CreateDSTNode(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }
  PHNodeIterator iterDst(dstNode);
  PHCompositeNode *detNode = dynamic_cast<PHCompositeNode *>(iterDst.findFirst("PHCompositeNode", "TPC"));
  if (!detNode)
  {
    detNode = new PHCompositeNode("TPC");
    dstNode->addNode(detNode);
  }
  TpcRawHitContainer *tpchitcont = findNode::getClass<TpcRawHitContainer>(detNode, m_rawHitContainerName);
  if (!tpchitcont)
  {
    tpchitcont = new TpcRawHitContainerv2();
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(tpchitcont, m_rawHitContainerName, "PHObject");
    detNode->addNode(newNode);
  }
}

void SingleTpcPoolInput::ConfigureStreamingInputManager()
{
  if (StreamingInputManager())
  {
    StreamingInputManager()->SetTpcBcoRange(m_BcoRange);
    StreamingInputManager()->SetTpcNegativeBco(m_NegativeBco);
  }
  return;
}
