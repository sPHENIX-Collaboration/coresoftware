#include "SingleTpcTimeFrameInput.h"
#include "TpcTimeFrameBuilder.h"

#include "Fun4AllStreamingInputManager.h"
#include "InputManagerType.h"

#include <ffarawobjects/TpcRawHitContainerv3.h>
#include <ffarawobjects/TpcRawHitv3.h>

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

SingleTpcTimeFrameInput::SingleTpcTimeFrameInput(const std::string &name)
  : SingleStreamingInput(name)
{
  SubsystemEnum(InputManagerType::TPC);
  m_rawHitContainerName = "TPCRAWHIT";
  plist = new Packet *[NTPCPACKETS];
}

SingleTpcTimeFrameInput::~SingleTpcTimeFrameInput()
{
  delete[] plist;

  for (auto iter : m_TpcTimeFrameBuilderMap)
  {
    delete iter.second;
  }
}

void SingleTpcTimeFrameInput::FillPool(const uint64_t targetBCO)
{
  if (Verbosity() > 1)
  {
    std::cout << "SingleTpcTimeFrameInput::FillPool: " << Name()
              << " Entry with targetBCO = 0x"<<std::hex << targetBCO
              <<"("<<std::dec<<targetBCO<<")" << std::endl;
  }

  if (AllDone())  // no more files and all events read
  {
    if (Verbosity() > 1)
    {
      std::cout << "SingleTpcTimeFrameInput::FillPool: " << Name()
                << " AllDone for targetBCO " << targetBCO << std::endl;
    }

    return;
  }
  while (GetEventiterator() == nullptr)  // at startup this is a null pointer
  {
    if (Verbosity() > 1)
    {
      std::cout << "SingleTpcTimeFrameInput::FillPool: " << Name()
                << " GetEventiterator == null for targetBCO " << targetBCO << std::endl;
    }

    if (!OpenNextFile())
    {
      if (Verbosity() > 1)
      {
        std::cout << "SingleTpcTimeFrameInput::FillPool: " << Name()
                  << " OpenNextFile() quit for targetBCO " << targetBCO << std::endl;
      }
      AllDone(1);
      return;
    }
  }
  //  std::set<uint64_t> saved_beamclocks;
  bool require_more_data = true;
  while (require_more_data)
  {
    if (Verbosity() > 3)
    {
      std::cout << "SingleTpcTimeFrameInput::FillPool: " << Name()
                << " require_more_data for targetBCO 0x"
                <<std::hex << targetBCO <<std::dec<< std::endl;
    }

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
      continue;
    }
    // int EventSequence = evt->getEvtSequence();
    int npackets = evt->getPacketList(plist, NTPCPACKETS);

    if (npackets >= NTPCPACKETS)
    {
      std::cout << PHWHERE << " Packets array size " << NTPCPACKETS
                << " too small for " << Name()
                << ", increase NTPCPACKETS and rebuild" << std::endl;
      exit(1);
    }

    for (int i = 0; i < npackets; i++)
    {
      // keep pointer to local packet
      auto &packet = plist[i];
      assert(packet);

      // get packet id
      const auto packet_id = packet->getIdentifier();

      if (m_TpcTimeFrameBuilderMap.find(packet_id) == m_TpcTimeFrameBuilderMap.end())
      {
        if (Verbosity() >= 1)
        {
          std::cout << __PRETTY_FUNCTION__ << ": Creating TpcTimeFrameBuilder for packet id: " << packet_id << std::endl;
        }

        m_TpcTimeFrameBuilderMap[packet_id] = new TpcTimeFrameBuilder(packet_id);
        m_TpcTimeFrameBuilderMap[packet_id]->setVerbosity(Verbosity());
      }

      if (Verbosity() > 1)
      {
        packet->identify();
      }

      assert(m_TpcTimeFrameBuilderMap[packet_id]);
      m_TpcTimeFrameBuilderMap[packet_id]->ProcessPacket(packet);
      require_more_data = require_more_data or m_TpcTimeFrameBuilderMap[packet_id]->isMoreDataRequired(targetBCO);

      delete packet;
      packet = nullptr;
    }  //     for (int i = 0; i < npackets; i++)

    require_more_data = false;
    for (auto & map_builder : m_TpcTimeFrameBuilderMap)
    {
      require_more_data |= map_builder.second->isMoreDataRequired(targetBCO);
    }

    if (not require_more_data)
    {
      for (auto & map_builder : m_TpcTimeFrameBuilderMap)
      {
        auto & timeframe = map_builder.second->getTimeFrame(targetBCO);

        for (auto newhit : timeframe)
        {
          StreamingInputManager()->AddTpcRawHit(targetBCO, newhit); 
        }

      }      
    } 

  }
  //    Print("HITS");
  //  } while (m_TpcRawHitMap.size() < 10 || CheckPoolDepth(m_TpcRawHitMap.begin()->first));
}

void SingleTpcTimeFrameInput::Print(const std::string & /*what*/) const
{
}

void SingleTpcTimeFrameInput::CleanupUsedPackets(const uint64_t bclk)
{
  if (Verbosity() > 2)
  {
    std::cout << "cleaning up bcos < 0x" << std::hex
              << bclk << std::dec << std::endl;
  }

  for (auto &iter : m_TpcTimeFrameBuilderMap)
  {
    iter.second->CleanupUsedPackets(bclk);
  }

}


void SingleTpcTimeFrameInput::ClearCurrentEvent()
{
  std::cout << "SingleTpcTimeFrameInput::ClearCurrentEvent() - deprecated " << std::endl;
  exit(1);

  // // called interactively, to get rid of the current event
  // uint64_t currentbclk = *m_BclkStack.begin();
  // //  std::cout << "clearing bclk 0x" << std::hex << currentbclk << std::dec << std::endl;
  // CleanupUsedPackets(currentbclk);
  // // m_BclkStack.erase(currentbclk);
  // // m_BeamClockFEE.erase(currentbclk);
  return;
}


void SingleTpcTimeFrameInput::CreateDSTNode(PHCompositeNode *topNode)
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
    tpchitcont = new TpcRawHitContainerv3();
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(tpchitcont, m_rawHitContainerName, "PHObject");
    detNode->addNode(newNode);
  }
}

void SingleTpcTimeFrameInput::ConfigureStreamingInputManager()
{
  if (StreamingInputManager())
  {
    StreamingInputManager()->SetTpcBcoRange(m_BcoRange);
    StreamingInputManager()->SetTpcNegativeBco(m_NegativeBco);
  }
  return;
}
