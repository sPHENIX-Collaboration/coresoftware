#include "CaloPacketSkimmer.h"

#include <qautils/QAHistManagerDef.h>

#include <ffarawobjects/CaloPacket.h>
#include <ffarawobjects/CaloPacketContainer.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>

#include <TH1.h>

#include <cassert>
#include <iostream>  // for operator<<, basic_ostream
#include <map>       // for operator!=, _Rb_tree_iterator
#include <utility>   // for pair
#include <variant>

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

  auto *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  h_aborted_events = new TH1D("h_caloskimmer_aborted_events", "Aborted Events", 1, 0, 1);
  h_aborted_events->SetDirectory(nullptr);
  hm->registerHisto(h_aborted_events);

  h_kept_events = new TH1D("h_caloskimmer_kept_events", "Kept Events", 1, 0, 1);
  h_kept_events->SetDirectory(nullptr);
  hm->registerHisto(h_kept_events);

  h_missing_packets = new TH1D("h_caloskimmer_missing_packets", "Missing Packets", 6000, 6000, 12000);
  h_missing_packets->SetDirectory(nullptr);
  hm->registerHisto(h_missing_packets);

  h_empty_packets = new TH1D("h_caloskimmer_empty_packets", "Empty Packets", 6000, 6000, 12000);
  h_empty_packets->SetDirectory(nullptr);
  hm->registerHisto(h_empty_packets);

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloPacketSkimmer::process_event(PHCompositeNode *topNode)
{
  bool bad_event = false;
  // loop over detectors and return DISCARDEVENT if any detector fails
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
        std::cout << "CaloPacketSkimmer::process_event - Detector " << getDetectorName(det) << " failed" << std::endl;
      }
      bad_event = true;
    }
  }

  h_kept_events->Fill(1);  // Increment histogram for kept events
  if (bad_event)
  {
    h_aborted_events->Fill(1);  // Increment histogram for aborted events
    return Fun4AllReturnCodes::DISCARDEVENT;
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
        h_empty_packets->Fill(pid);
        return Fun4AllReturnCodes::DISCARDEVENT;
      }
    }
    else
    {
      if (Verbosity() > 1)
      {
        std::cout << "CaloPacketSkimmer: Packet " << pid << " not found for detector " << getDetectorName(dettype) << std::endl;
      }
      h_missing_packets->Fill(pid);
      return Fun4AllReturnCodes::DISCARDEVENT;
    }
    return Fun4AllReturnCodes::EVENT_OK;
  };
  // loop over packets
  bool bad_event = false;
  for (int pid = startID; pid <= endID; pid++)
  {
    if (!m_PacketNodesFlag)
    {
      if (auto *cont = std::get_if<CaloPacketContainer *>(&event))
      {
        if (!*cont)
        {
          std::cout << PHWHERE << " CaloPacketContainer not found: " << nodemap.find(dettype)->second << std::endl;
          return Fun4AllReturnCodes::DISCARDEVENT;
        }
        CaloPacket *packet = (*cont)->getPacketbyId(pid);
        if (process_packet(packet, pid) == Fun4AllReturnCodes::DISCARDEVENT)
        {
          bad_event = true;
        }
      }
      else if (auto *_event = std::get_if<Event *>(&event))
      {
        Packet *packet = (*_event)->getPacket(pid);
        if (process_packet(packet, pid) == Fun4AllReturnCodes::DISCARDEVENT)
        {
          // I think it is safe to delete a nullptr...yes, tessted
          delete packet;
          bad_event = true;
        }

        delete packet;
      }
    }
    else
    {
      CaloPacket *calopacket = findNode::getClass<CaloPacket>(topNode, pid);
      if (process_packet(calopacket, pid) == Fun4AllReturnCodes::DISCARDEVENT)
      {
        bad_event = true;
      }
    }
  }

  if (bad_event)
  {
    return Fun4AllReturnCodes::DISCARDEVENT;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloPacketSkimmer::EndRun(const int /*runnumber*/)
{
  std::cout << "CaloPacketSkimmer::EndRun - Analyzed events: " << h_kept_events->GetEntries() << std::endl;
  std::cout << "CaloPacketSkimmer::EndRun - Discarded events: " << h_aborted_events->GetEntries() << std::endl;
  // loop over and print the missing and empty packets
  for (int pid = 6000; pid <= 12000; pid++)
  {
    int bin = h_missing_packets->FindBin(pid);
    if (h_missing_packets->GetBinContent(bin) > 0)
    {
      std::cout << "CaloPacketSkimmer::EndRun - Missing packet: " << pid << " occurences: " << h_missing_packets->GetBinContent(bin) << std::endl;
    }
    if (h_empty_packets->GetBinContent(bin) > 0)
    {
      std::cout << "CaloPacketSkimmer::EndRun - Empty packet: " << pid << " occurences: " << h_empty_packets->GetBinContent(bin) << std::endl;
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
