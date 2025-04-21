#include "ClockDiffCheck.h"

#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <ffarawobjects/CaloPacket.h>
#include <ffarawobjects/CaloPacketContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHDataNode.h>
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/getClass.h>

#include <TH1.h>
#include <TSystem.h>

#include <bitset>
#include <iostream>  // for operator<<, endl, basic_ost...
#include <utility>   // for pair
#include <vector>    // for vector

//____________________________________________________________________________..
ClockDiffCheck::ClockDiffCheck(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
int ClockDiffCheck::Init(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ClockDiffCheck::process_event(PHCompositeNode *topNode)
{
  for (auto &iter : m_PacketStuffMap)
  {
    std::get<0>(iter.second) = std::get<1>(iter.second);
    std::get<4>(iter.second) = false;
  }
  OfflinePacket *pkt = findNode::getClass<OfflinePacket>(topNode, "GL1Packet");
  if (pkt)
  {
    FillPacketDiff(pkt);
  }

  std::vector<std::string> nodenames{"CEMCPackets", "HCALPackets", "MBDPackets", "SEPDPackets", "ZDCPackets"};
  for (const auto &iter : nodenames)
  {
    CaloPacketContainer *cemccont = findNode::getClass<CaloPacketContainer>(topNode, iter);
    if (!cemccont)
    {
      //      std::cout << "could not find " << iter << " node" << std::endl;
    }
    else
    {
      FillCaloClockDiff(cemccont);
    }
  }

  std::vector<unsigned int> badPackets;
  uint64_t refdiff = std::numeric_limits<uint64_t>::max();
  auto itergl1 = m_PacketStuffMap.find(14001);
  if (itergl1 != m_PacketStuffMap.end())
  {
    refdiff = std::get<2>(itergl1->second);
  }
  for (auto &iter : m_PacketStuffMap)
  {
    if (!std::get<4>(iter.second))
    {
      std::get<0>(iter.second) = std::numeric_limits<uint64_t>::max();
    }
    else
    {
      if (Verbosity() > 2)
      {
        std::cout << "looking at " << iter.first
                  << ", prev bco: " << std::hex << std::get<0>(iter.second)
                  << ", curr bco: " << std::get<1>(iter.second)
                  << ", clkdiff: " << std::get<2>(iter.second) << std::dec
                  << ", valid: " << std::get<4>(iter.second) << std::endl;
      }
      if (refdiff == std::numeric_limits<uint64_t>::max())
      {
        refdiff = std::get<2>(iter.second);
      }
      else
      {
        if ((refdiff & 0xFFFFFFFFU) != (std::get<2>(iter.second) & 0xFFFFFFFFU))
        {
          badPackets.push_back(iter.first);
          static int nprint = 0;
          if (nprint < 1000 || Verbosity() > 1)
          {
            std::bitset<32> x(refdiff);
            std::bitset<32> y0(std::get<0>(iter.second));
            std::bitset<32> y1(std::get<1>(iter.second));
            std::bitset<32> y2(std::get<2>(iter.second));
            std::cout << "packet " << iter.first << " had different clock diff: 0x" << std::hex
                      << std::get<1>(iter.second) << ", ref diff: 0x" << refdiff << std::dec << std::endl;

            std::cout << "reff: " << x << std::endl;
            std::cout << "this: " << y2 << std::endl;

            std::cout << "prev: " << y0 << std::endl;
            std::cout << "curr: " << y1 << std::endl;
            nprint++;
          }
        }
      }
    }
  }

  for (const auto &nodeiter : nodenames)
  {
    CaloPacketContainer *container = findNode::getClass<CaloPacketContainer>(topNode, nodeiter);
    if (!container)
    {
      continue;
    }
    if (delBadPkts)
    {
      for (unsigned int i = 0; i < container->get_npackets(); i++)
      {
        unsigned int packetID = container->getPacket(i)->getIdentifier();
        for (unsigned int badPacket : badPackets)
        {
          if (badPacket == packetID)
          {
            if (Verbosity() > 1)
            {
              std::cout << "Dropping packet " << container->getPacket(i)->getIdentifier() << " for XMIT clock mismatch" << std::endl;
            }
            container->deletePacket(container->getPacket(i));
            break;
          }
        }
      }
    }

    std::vector<std::vector<int>> EvtCounts;
    std::vector<int> NrAndCount(2);
    NrAndCount[1] = 1;
    int counter = 0;
    int bestEvt = -1;
    int bestEvtCnt = 0;
    unsigned int npacket = container->get_npackets();
    for (unsigned int i = 0; i < npacket; i++)
    {
      CaloPacket *packet = container->getPacket(i);
      if (packet)
      {
        int nrModules = packet->iValue(0, "NRMODULES");
        for (int j = 0; j < nrModules; j++)
        {
          int k;
          for (k = 0; k < counter; k++)
          {
            if (EvtCounts[k][0] == packet->iValue(j, "FEMEVTNR"))
            {
              EvtCounts[k][1]++;
              break;
            }
          }
          if (k >= counter)
          {
            NrAndCount[0] = packet->iValue(j, "FEMEVTNR");
            EvtCounts.push_back(NrAndCount);
            counter++;
          }
        }
      }
    }
    if (counter > 1)
    {
      for (int i = 0; i < counter; i++)
      {
        if (bestEvtCnt < EvtCounts[i][1])
        {
          bestEvtCnt = EvtCounts[i][1];
          bestEvt = EvtCounts[i][0];
        }
      }
      for (unsigned int i = 0; i < npacket; ++i)
      {
        CaloPacket *packet = container->getPacket(i);
        if (packet)
        {
          for (int j = 0; j < packet->iValue(0, "NRMODULES"); j++)
          {
            if (packet->iValue(j, "FEMEVTNR") != bestEvt && bestEvt != -1)
            {
              static int icnt = 0;
              if (icnt < 1000)
              {
                std::cout << "found different FEM clock for packet " << packet->getIdentifier() << std::endl;
                icnt++;
              }
              if (delBadPkts)
              {
                if (Verbosity() > 1)
                {
                  std::cout << "deleting packet " << packet->getIdentifier()
                            << " with fem clock mismatch" << std::endl;
                }
                container->deletePacket(packet);
              }
              break;
            }
          }
        }
      }
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

void ClockDiffCheck::FillCaloClockDiff(CaloPacketContainer *pktcont)
{
  for (unsigned int i = 0; i < pktcont->get_npackets(); i++)
  {
    unsigned int packetid = pktcont->getPacket(i)->getIdentifier();
    if (m_PacketStuffMap.find(packetid) == m_PacketStuffMap.end())
    {
      std::string hname = "clkdiff" + std::to_string(packetid);
      TH1 *h1 = new TH1F(hname.c_str(), hname.c_str(), 100, 0, 99);
      m_PacketStuffMap[packetid] = std::make_tuple(std::numeric_limits<uint64_t>::max(), std::numeric_limits<uint64_t>::max(), std::numeric_limits<uint64_t>::max(), h1, false);
      if (Verbosity() > 3)
      {
        std::cout << "Add tuple for " << packetid << std::endl;
        auto &pktiter = m_PacketStuffMap[packetid];
        std::cout << PHWHERE << "packet init " << packetid << std::hex
                  << ", clk: " << std::get<1>(pktiter)
                  << ", clkdiff: " << std::get<2>(pktiter) << std::dec << ", valid: " << std::get<4>(pktiter)
                  << std::endl;
      }
    }
    else
    {
      auto &pktiter = m_PacketStuffMap[packetid];
      uint64_t clk = pktcont->getPacket(i)->getBCO();
      uint64_t clkdiff = std::numeric_limits<uint64_t>::max();
      std::get<1>(pktiter) = clk;
      // only calculate clk diff and correct clock for rollover if previous clk is set (default is max uint64)
      if (std::get<0>(pktiter) < std::numeric_limits<uint64_t>::max())
      {
        if (clk < std::get<0>(pktiter))
        {
          clk |= 0x100000000U;
        }
        clkdiff = clk - std::get<0>(pktiter);
        clk &= 0xFFFFFFFF;
        std::get<2>(pktiter) = clkdiff;
        std::get<4>(pktiter) = true;
      }
      if (Verbosity() > 2)
      {
	std::cout << "packet " << packetid << ", clk: " << std::hex << clk
		  << ", clk(tup): " << std::get<1>(pktiter) << ", diff: " << clkdiff
		  << ", diff(tup): " << std::get<2>(pktiter) << std::dec << ", valid: " << std::get<4>(pktiter)
		  << std::endl;
      }
    }
  }
}

void ClockDiffCheck::FillPacketDiff(OfflinePacket *pkt)
{
  unsigned int packetid = pkt->getIdentifier();
  uint64_t clk = (pkt->getBCO() & 0xFFFFFFFF);
  if (m_PacketStuffMap.find(packetid) == m_PacketStuffMap.end())
  {
    std::string hname = "clkdiff" + std::to_string(packetid);
    TH1 *h1 = new TH1F(hname.c_str(), hname.c_str(), 100, 0, 99);
    m_PacketStuffMap[packetid] = std::make_tuple(std::numeric_limits<uint64_t>::max(), std::numeric_limits<uint64_t>::max(), std::numeric_limits<uint64_t>::max(), h1, false);
    if (Verbosity() > 3)
    {
      std::cout << "Add tuple for " << packetid << std::endl;
    }
  }
  else
  {
    auto &pktiter = m_PacketStuffMap[packetid];
    uint64_t clkdiff = std::numeric_limits<uint64_t>::max();
    if (std::get<0>(pktiter) < std::numeric_limits<uint64_t>::max())
    {
      if (clk < std::get<0>(pktiter))
      {
        clk |= 0x100000000U;
      }
      clkdiff = clk - std::get<0>(pktiter);
      clk &= 0xFFFFFFFF;
    }
    std::get<1>(pktiter) = clk;
    std::get<2>(pktiter) = clkdiff;
    std::get<4>(pktiter) = true;
    if (Verbosity() > 2)
    {
      std::cout << "packet " << packetid << ", clk: " << std::hex << clk
                << ", clk(tup): " << std::get<1>(pktiter) << ", diff: " << clkdiff
                << ", diff(tup): " << std::get<2>(pktiter) << std::dec << ", valid: " << std::get<4>(pktiter)
                << std::endl;
    }
  }
}
