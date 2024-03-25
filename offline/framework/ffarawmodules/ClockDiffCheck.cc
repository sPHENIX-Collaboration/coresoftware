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
  OfflinePacket *pkt = findNode::getClass<OfflinePacket>(topNode,"GL1Packet");
  if (pkt)
  {
    FillPacketDiff(pkt);
  }

  std::vector<std::string> nodenames {"CEMCPackets", "HCALPackets", "MBDPackets"};
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
  uint64_t refdiff = std::numeric_limits<uint64_t>::max();
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
		  << ", clkdiff: " <<  std::get<2>(iter.second) << std::dec
		  << ", valid: " << std::get<4>(iter.second) << std::endl;
      }
      if (refdiff == std::numeric_limits<uint64_t>::max())
      {
	refdiff = std::get<2>(iter.second);
      }
      else
      {
	if (refdiff != std::get<2>(iter.second))
	{
	  std::bitset<16> x(refdiff);
	  std::bitset<16> y0(std::get<0>(iter.second));
	  std::bitset<16> y1(std::get<1>(iter.second));
	  std::bitset<16> y2(std::get<2>(iter.second));
	  std::cout << "packet " << iter.first << " had different clock diff: 0x" << std::hex
		    <<  std::get<1>(iter.second) << ", ref diff: 0x" << refdiff << std::dec << std::endl;

	  std::cout << "reff: " << x << std::endl;
	  std::cout << "this: " << y2 << std::endl;

	  std::cout << "prev: " << y0 << std::endl;
	  std::cout << "curr: " << y1 << std::endl;
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
    uint64_t clk = pktcont->getPacket(i)->getBCO();
    if (m_PacketStuffMap.find(packetid) == m_PacketStuffMap.end())
    {
      std::string hname = "clkdiff" + std::to_string(packetid);
      TH1 *h1 = new TH1F(hname.c_str(), hname.c_str(),100,0,99);
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
	  clk |= 0x10000;
	}
	clkdiff = clk - std::get<0>(pktiter);
	clk &= 0xFFFF;
      }
      std::get<1>(pktiter) = clk;
      std::get<2>(pktiter) = clkdiff;
      std::get<4>(pktiter) = true;
      if (Verbosity() > 2)
      {
	std::cout << "packet " << packetid << ", clk: " << std::hex << clk
		  << ", tup: " << std::get<1>(pktiter) << ", diff: " << clkdiff
		  << ", tup: " << std::get<2>(pktiter) << std::dec << ", tag: " << std::get<4>(pktiter)
		  << std::endl;
      }
    }
  }
}

void ClockDiffCheck::FillPacketDiff(OfflinePacket *pkt)
{
    unsigned int packetid = pkt->getIdentifier();
    uint64_t clk = (pkt->getBCO()&0xFFFF);
    if (m_PacketStuffMap.find(packetid) == m_PacketStuffMap.end())
    {
      std::string hname = "clkdiff" + std::to_string(packetid);
      TH1 *h1 = new TH1F(hname.c_str(), hname.c_str(),100,0,99);
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
	  clk |= 0x10000;
	}
	clkdiff = clk - std::get<0>(pktiter);
	clk &= 0xFFFF;
      }
      std::get<1>(pktiter) = clk;
      std::get<2>(pktiter) = clkdiff;
      std::get<4>(pktiter) = true;
      if (Verbosity() > 2)
      {
	std::cout << "packet " << packetid << ", clk: " << std::hex << clk
		  << ", tup: " << std::get<1>(pktiter) << ", diff: " << clkdiff
		  << ", tup: " << std::get<2>(pktiter) << std::dec << ", tag: " << std::get<4>(pktiter)
		  << std::endl;
      }
    }
}
