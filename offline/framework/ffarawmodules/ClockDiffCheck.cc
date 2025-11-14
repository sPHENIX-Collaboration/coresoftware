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
#include <phool/PHPointerListIterator.h>
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

ClockDiffCheck::~ClockDiffCheck()
{
  for (auto &[packetid, tup] : m_PacketStuffMap)
  {
    delete std::get<3>(tup);
    std::get<3>(tup) = nullptr;
  }
}

//____________________________________________________________________________..
int ClockDiffCheck::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << "ClockDiffCheck: " << "Could not find DST Node" << std::endl;
    gSystem->Exit(1);
    exit(1);
  }
  PHNodeIterator iterDst(dstNode);
  PHCompositeNode *pktNode = dynamic_cast<PHCompositeNode *>(iterDst.findFirst("PHCompositeNode", "Packets"));
  PHCompositeNode *pktKeepNode = dynamic_cast<PHCompositeNode *>(iterDst.findFirst("PHCompositeNode", "PacketsKeep"));
  if (pktNode)  // old combined packet containers
  {
    PHNodeIterator iterPkt(pktNode);
    PHPointerListIterator<PHNode> nodeIter(iterPkt.ls());
    PHNode *thisNode;
    while ((thisNode = nodeIter()))
    {
      m_PacketNodeNames.push_back(thisNode->getName());
    }
  }
  if (pktKeepNode)  // old combined packet containers
  {
    PHNodeIterator iterPkt(pktKeepNode);
    PHPointerListIterator<PHNode> nodeIter(iterPkt.ls());
    PHNode *thisNode;
    while ((thisNode = nodeIter()))
    {
      m_PacketNodeNames.push_back(thisNode->getName());
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ClockDiffCheck::process_event(PHCompositeNode *topNode)
{
  count++;
  //  PHNodeIterator
  PHNodeIterator topnodeiter(topNode);
  for (auto &iter : m_PacketStuffMap)
  {
    std::get<0>(iter.second) = std::get<1>(iter.second);
    std::get<4>(iter.second) = false;
  }
  OfflinePacket *pkt = findNode::getClass<OfflinePacket>(topNode, 14001);
  if (pkt)
  {
    FillPacketDiff(pkt);
  }
  pkt = findNode::getClass<OfflinePacket>(topNode, "GL1Packet");
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
      //      std::cout << "ClockDiffCheck: " << "could not find " << iter << " node" << std::endl;
    }
    else
    {
      FillCaloClockDiff(cemccont);
    }
  }

  for (const auto &iter : m_PacketNodeNames)
  {
    if (iter == "14001") { continue;
}
    CaloPacket *calopacket = findNode::getClass<CaloPacket>(topNode, iter);
    if (!calopacket)
    {
      std::cout << "ClockDiffCheck: " << "could not find " << iter << " node" << std::endl;
    }
    else
    {
      if (calopacket->getStatus() != OfflinePacket::PACKET_OK)
      {
        static int npacketprnt = 0;
        if ( npacketprnt < 100)
        {
          std::cout << "packet " << calopacket->getIdentifier() << " status marked as " << calopacket->getStatus() << ". Skipping diff check " << count << std::endl;
          npacketprnt++;
        }
        std::get<1>(m_PacketStuffMap[calopacket->getIdentifier()]) = std::numeric_limits<uint64_t>::max();
        continue;
      }
      FillCaloClockDiffSngl(calopacket);
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
        std::cout << "ClockDiffCheck: " << "looking at " << iter.first
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
            std::cout << "ClockDiffCheck: " << "packet " << iter.first << " had different clock diff: 0x" << std::hex
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
              std::cout << "ClockDiffCheck: " << "Dropping packet " << container->getPacket(i)->getIdentifier() << " for XMIT clock mismatch" << std::endl;
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
                std::cout << "ClockDiffCheck: " << "found different FEM clock for packet " << packet->getIdentifier() << std::endl;
                icnt++;
              }
              if (delBadPkts)
              {
                if (Verbosity() > 1)
                {
                  std::cout << "ClockDiffCheck: " << "deleting packet " << packet->getIdentifier()
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
  for (const auto &iter : m_PacketNodeNames)
  {
    CaloPacket *calopacket = findNode::getClass<CaloPacket>(topNode, iter);
    if (!calopacket)
    {
      continue;
    }
    if (delBadPkts)
    {
      unsigned int packetID = calopacket->getIdentifier();
      for (unsigned int badPacket : badPackets)
      {
        if (badPacket == packetID)
        {
          static int icnt = 0;
          if (icnt < 1000)
          {
            std::cout << "ClockDiffCheck: " << "Dropping packet " << calopacket->getIdentifier() << " for XMIT clock mismatch" << std::endl;
            icnt++;
          }
          calopacket->Reset();
          break;
        }
      }
    }
    bool FemClockOk = CheckFemEventNr(calopacket);
    if (!FemClockOk)
    {
      static int icnt = 0;
      if (icnt < 100)
      {
        std::cout << "FemClockOk failed. Resetting packet " << calopacket->getIdentifier()
          << " with fem event and clock mismatch" << std::endl;
        icnt++;
      }
      calopacket->Reset();
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void ClockDiffCheck::FillCaloClockDiff(CaloPacketContainer *pktcont)
{
  for (unsigned int i = 0; i < pktcont->get_npackets(); i++)
  {
    FillCaloClockDiffSngl(pktcont->getPacket(i));
  }
  return;
}

void ClockDiffCheck::FillCaloClockDiffSngl(CaloPacket *calopkt)
{
  unsigned int packetid = calopkt->getIdentifier();
  if (!m_PacketStuffMap.contains(packetid))
  {
    std::string hname = "clkdiff" + std::to_string(packetid);
    TH1 *h1 = new TH1F(hname.c_str(), hname.c_str(), 100, 0, 99);
    m_PacketStuffMap[packetid] = std::make_tuple(std::numeric_limits<uint64_t>::max(), std::numeric_limits<uint64_t>::max(), std::numeric_limits<uint64_t>::max(), h1, false);
    if (Verbosity() > 3)
    {
      std::cout << "ClockDiffCheck: " << "Add tuple for " << packetid << std::endl;
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
    uint64_t clk = calopkt->getBCO();
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
      std::cout << "ClockDiffCheck: " << "packet " << packetid << ", clk: " << std::hex << clk
                << ", clk(tup): " << std::get<1>(pktiter) << ", diff: " << clkdiff
                << ", diff(tup): " << std::get<2>(pktiter) << std::dec << ", valid: " << std::get<4>(pktiter)
                << std::endl;
    }
  }
}

void ClockDiffCheck::FillPacketDiff(OfflinePacket *pkt)
{
  unsigned int packetid = pkt->getIdentifier();
  uint64_t clk = (pkt->getBCO() & 0xFFFFFFFF);
  if (!m_PacketStuffMap.contains(packetid))
  {
    std::string hname = "clkdiff" + std::to_string(packetid);
    TH1 *h1 = new TH1F(hname.c_str(), hname.c_str(), 100, 0, 99);
    m_PacketStuffMap[packetid] = std::make_tuple(std::numeric_limits<uint64_t>::max(), std::numeric_limits<uint64_t>::max(), std::numeric_limits<uint64_t>::max(), h1, false);
    if (Verbosity() > 3)
    {
      std::cout << "ClockDiffCheck: " << "Add tuple for " << packetid << std::endl;
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
      std::cout << "ClockDiffCheck: " << "packet " << packetid << ", clk: " << std::hex << clk
                << ", clk(tup): " << std::get<1>(pktiter) << ", diff: " << clkdiff
                << ", diff(tup): " << std::get<2>(pktiter) << std::dec << ", valid: " << std::get<4>(pktiter)
                << std::endl;
    }
  }
}

bool ClockDiffCheck::CheckFemEventNr(CaloPacket *calopkt)
{
  int nrModules = calopkt->iValue(0, "NRMODULES");
  std::set<int> EventNoSet;
  for (int j = 0; j < nrModules; j++)
  {
    if (calopkt->getFemStatus(j) == CaloPacket::FEM_OK)
    {
      EventNoSet.insert(calopkt->iValue(j, "FEMEVTNR"));
    }
  }
  if (EventNoSet.size() > 1)
  {
    // at least one packet (6024) has a stuck bit in the fem event nr, check fem clock counter in this case
    // if they are identical FEM is good (not checked if the FEM clock is stuck though)
    std::set<int> FemClockSet;
    for (int j = 0; j < nrModules; j++)
    {
      FemClockSet.insert(calopkt->iValue(j, "FEMCLOCK"));
    }
    if (FemClockSet.size() == 1)
    {
      static int icnt = 0;
      if (icnt < 100)
      {
        icnt++;
        std::cout << "ClockDiffCheck: " << "clk check Packet " << calopkt->getIdentifier() << " has not unique event numbers"
                  << " but FEM Clock counters are identical" << std::endl;
      }
    }
    else  // event nr and fem clock differ
    {
      // now lets find which one is the outlier
      static int icnt = 0;
      if (icnt < 1000)
      {
        icnt++;
        std::cout << "ClockDiffCheck: " << "resetting packet " << calopkt->getIdentifier()
                  << " with fem event and clock mismatch" << std::endl;
        std::map<int, int> EventMap;
        std::map<int, int> ClockMap;
        for (int j = 0; j < nrModules; j++)
        {
          EventMap[calopkt->iValue(j, "FEMEVTNR")]++;
          ClockMap[calopkt->iValue(j, "FEMCLOCK")]++;
        }
        for (const auto iterA : EventMap)
        {
          std::cout << "ClockDiffCheck: " << "Event Nr : " << iterA.first << " shows up " << iterA.second << " times"
                    << std::hex << ", Event Nr 0x" << iterA.first << std::dec << std::endl;
        }
        for (const auto iterA : ClockMap)
        {
          std::cout << "ClockDiffCheck: " << "Clock : 0x" << std::hex << iterA.first << std::dec
                    << " shows up " << iterA.second << " times" << std::endl;
        }
      }
      return false;
    }
  }
  return true;
}
