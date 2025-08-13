
#include "CaloPacketSkimmer.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>
#include <ffarawobjects/CaloPacket.h>
#include <ffarawobjects/CaloPacketContainer.h>
#include <variant>

#include <phool/getClass.h>

#include <iostream>  // for operator<<, basic_ostream
#include <map>       // for operator!=, _Rb_tree_iterator
#include <utility>   // for pair

static const std::map<CaloTowerDefs::DetectorSystem, std::string> nodemap{
    {CaloTowerDefs::CEMC, "CEMCPackets"},
    {CaloTowerDefs::HCALIN, "HCALPackets"},
    {CaloTowerDefs::HCALOUT, "HCALPackets"},
    {CaloTowerDefs::ZDC, "ZDCPackets"},
    {CaloTowerDefs::SEPD, "SEPDPackets"},
    {CaloTowerDefs::MBD, "MBDPackets"}};

//____________________________________________________________________________..
CaloPacketSkimmer::CaloPacketSkimmer(const std::string &name)
  : SubsysReco(name)
{
}

int CaloPacketSkimmer::InitRun(PHCompositeNode * /*topNode*/)
{
  // Initialize the detector systems to be used
  if (Verbosity() > 0)
  {
    std::cout << "CaloPacketSkimmer::InitRun - Initializing with detectors: ";
    for (const auto &det : m_CaloDetectors)
    {
      std::cout << getDetectorName(det) << " ";
    }
    std::cout << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloPacketSkimmer::process_event(PHCompositeNode *topNode)
{
  // loop over detectors and return ABORTEVENT if any detector fails
  for (const auto &det : m_CaloDetectors)
  {
    if (Verbosity() > 1)
    {
      std::cout << "CaloPacketSkimmer::process_event - Processing detector: " << getDetectorName(det) << std::endl;
    }
    int ret = processDetector(topNode, det);
    if (ret != Fun4AllReturnCodes::EVENT_OK)
    {
      if (Verbosity() > 1)
      {
        std::cout << "CaloPacketSkimmer: Skipping event due to failure in detector: " << getDetectorName(det) << std::endl;
      }
      return ret;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloPacketSkimmer::processDetector(PHCompositeNode *topNode, CaloTowerDefs::DetectorSystem dettype)
{
  // Process the detector data and filter packets based on the accepted IDs
  auto packetRange = getPacketRange(dettype);
  int startID = packetRange.first;
  int endID = packetRange.second;

  if (Verbosity() > 4)
  {
    std::cout << "Processing " << getDetectorName(dettype) << " from ID " << startID << " to " << endID << std::endl;
  }
  std::variant<CaloPacketContainer *, Event *> event;
  if (m_UseOfflinePacketFlag)
  {
    CaloPacketContainer *cont = findNode::getClass<CaloPacketContainer>(topNode, nodemap.find(dettype)->second);
    if (!cont)
    {
      for (int pid = startID; pid <= endID; pid++)
      {
        if (findNode::getClass<CaloPacket>(topNode, pid))
        {
          m_PacketNodesFlag = true;
          break;
        }
      }
      if (!m_PacketNodesFlag)
      {
        if (Verbosity() > 0)
        {
          std::cout << "CaloPacketSkimmer: unable to find node " << nodemap.find(dettype)->second << "  skiping event" << std::endl;
        }
        return Fun4AllReturnCodes::EVENT_OK;
      }
    }
    event = cont;
  }
  else
  {
    Event *_event = findNode::getClass<Event>(topNode, "PRDF");
    if (_event == nullptr)
    {
      if (Verbosity() > 0)
      {
        std::cout << PHWHERE << " Event node not found, skipping event" << std::endl;
      }
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    if (_event->getEvtType() != DATAEVENT)
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    event = _event;
  }
  // since the function call on Packet and CaloPacket is the same, maybe we can use lambda? ... yes, we can!
  auto process_packet = [&](auto *packet, int pid)
  {
    if (packet)
    {
      int nchannels = packet->iValue(0, "CHANNELS");
      if (nchannels == 0)
      {
        if (Verbosity() > 1)
        {
          std::cout << "CaloPacketSkimmer: No channels in packet " << pid << " for detector " << getDetectorName(dettype) << std::endl;
        }
        return Fun4AllReturnCodes::ABORTEVENT;
      }
    }
    else
    {
      if (Verbosity() > 1)
      {
        std::cout << "CaloPacketSkimmer: Packet " << pid << " not found for detector " << getDetectorName(dettype) << std::endl;
      }
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    return Fun4AllReturnCodes::EVENT_OK;
  };
  // loop over packets
  for (int pid = startID; pid <= endID; pid++)
  {
    if (!m_PacketNodesFlag)
    {
      if (auto *cont = std::get_if<CaloPacketContainer *>(&event))
      {
        if (!*cont)
        {
          std::cout << PHWHERE << " CaloPacketContainer not found: " << nodemap.find(dettype)->second << std::endl;
          return Fun4AllReturnCodes::ABORTEVENT;
        }
        CaloPacket *packet = (*cont)->getPacketbyId(pid);
        if (process_packet(packet, pid) == Fun4AllReturnCodes::ABORTEVENT)
        {
          return Fun4AllReturnCodes::ABORTEVENT;
        }
      }
      else if (auto *_event = std::get_if<Event *>(&event))
      {
        Packet *packet = (*_event)->getPacket(pid);
        if (process_packet(packet, pid) == Fun4AllReturnCodes::ABORTEVENT)
        {
          // I think it is safe to delete a nullptr...yes, tessted
          delete packet;
          return Fun4AllReturnCodes::ABORTEVENT;
        }

        delete packet;
      }
    }
    else
    {
      CaloPacket *calopacket = findNode::getClass<CaloPacket>(topNode, pid);
      if (process_packet(calopacket, pid) == Fun4AllReturnCodes::ABORTEVENT)
      {
        return Fun4AllReturnCodes::ABORTEVENT;
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}
