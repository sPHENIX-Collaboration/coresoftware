#include "CaloFittingQA.h"

// Calo includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <caloreco/CaloTowerDefs.h>

#include <cdbobjects/CDBTTree.h>  // for CDBTTree

#include <ffarawobjects/CaloPacket.h>
#include <ffarawobjects/CaloPacketContainer.h>

#include <ffamodules/CDBInterface.h>

#include <ffaobjects/EventHeader.h>

#include <qautils/QAHistManagerDef.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE


#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>

#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TSystem.h>

#include <boost/format.hpp>

#include <cassert>
#include <cmath>  // for log10, pow, sqrt, abs, M_PI
#include <cstdint>
#include <iostream>  // for operator<<, endl, basic_...
#include <limits>
#include <map>  // for operator!=, _Rb_tree_con...
#include <string>
#include <utility>  // for pair
#include <variant>

static const std::map<CaloTowerDefs::DetectorSystem, std::string> nodemap{
    {CaloTowerDefs::CEMC, "CEMCPackets"},
    {CaloTowerDefs::HCALIN, "HCALPackets"},
    {CaloTowerDefs::HCALOUT, "HCALPackets"},
    {CaloTowerDefs::ZDC, "ZDCPackets"},
    {CaloTowerDefs::SEPD, "SEPDPackets"}};

CaloFittingQA::CaloFittingQA(const std::string& name)
  : SubsysReco(name)
{
}

CaloFittingQA::~CaloFittingQA()
{
  delete cdbttree;
}

int CaloFittingQA::Init(PHCompositeNode* /*unused*/)
{
  if (Verbosity() > 0)
  {
    std::cout << "In CaloFittingQA::Init" << std::endl;
  }

  createHistos();

  if (Verbosity() > 0)
  {
    std::cout << "Leaving CaloFittingQA::Init" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloFittingQA::InitRun(PHCompositeNode* /*unused*/)
{
  if (Verbosity() > 0)
  {
    std::cout << "In CaloFittingQA::InitRun" << std::endl;
  }

  //===============
  // conditions DB flags
  //===============
  m_calibName = "cemc_adcskipmask";
  m_fieldname = "adcskipmask";
  std::string calibdir = CDBInterface::instance()->getUrl(m_calibName);
  if (calibdir.empty())
  {
    std::cout << PHWHERE << "ADC Skip mask not found in CDB, not even in the default... " << std::endl;
    exit(1);
  }
  cdbttree = new CDBTTree(calibdir.c_str());

  if (Verbosity() > 0)
  {
    std::cout << "Leaving CaloFittingQA::InitRun" << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloFittingQA::process_event(PHCompositeNode* topNode)
{
  _eventcounter++;
  //  std::cout << "In process_event" << std::endl;
  process_towers(topNode);

  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloFittingQA::process_towers(PHCompositeNode* topNode)
{
  //---------------------------Event header--------------------------------//
  EventHeader* eventheader =
      findNode::getClass<EventHeader>(topNode, "EventHeader");
  int event_number = 0;
  if (eventheader)
  {
    if (eventheader->isValid())
    {
      event_number = eventheader->get_EvtSequence();
    }
  }
  else
  {
    std::cout << "CaloFittingQA::process_towers()  No event header" << std::endl;
  }

  if (Verbosity() > 0)
  {
    std::cout << _eventcounter << std::endl;
    std::cout << "Event header evtsequence " << event_number << std::endl;
  }
  //  std::cout << "In process_towers" << std::endl;
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  //-------------------------- raw waveforms ------------------------------//
  std::vector<std::vector<float>> cemc_waveforms;
  std::vector<std::vector<float>> ihcal_waveforms;
  std::vector<std::vector<float>> ohcal_waveforms;
  TowerInfoContainer* cemc_sim_waveforms = nullptr;
  TowerInfoContainer* ihcal_sim_waveforms = nullptr;
  TowerInfoContainer* ohcal_sim_waveforms = nullptr;
  if (m_SimFlag) 
  {
    cemc_sim_waveforms = findNode::getClass<TowerInfoContainer>(topNode, "WAVEFORM_CEMC");
    if (!cemc_sim_waveforms) 
    {
      std::cout << PHWHERE << "WAVEFORM_CEMC node missing. Skipping event." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    ihcal_sim_waveforms = findNode::getClass<TowerInfoContainer>(topNode, "WAVEFORM_HCALIN");
    if (!ihcal_sim_waveforms) 
    {
      std::cout << PHWHERE << "WAVEFORM_HCALIN node missing. Skipping event." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    ohcal_sim_waveforms = findNode::getClass<TowerInfoContainer>(topNode, "WAVEFORM_HCALOUT");
    if (!ohcal_sim_waveforms) 
    {
      std::cout << PHWHERE << "WAVEFORM_HCALOUT node missing. Skipping event." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }
  else 
  {
    if(process_data(topNode, CaloTowerDefs::CEMC, cemc_waveforms) == Fun4AllReturnCodes::ABORTEVENT)
    {
      std::cout << PHWHERE << "CEMCPackets node failed unpacking. Skipping event." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    if(process_data(topNode, CaloTowerDefs::HCALIN, ihcal_waveforms) == Fun4AllReturnCodes::ABORTEVENT)
    {
      std::cout << PHWHERE << "HCALPackets node failed unpacking. Skipping event." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    if(process_data(topNode, CaloTowerDefs::HCALOUT, ohcal_waveforms) == Fun4AllReturnCodes::ABORTEVENT)
    {
      std::cout << PHWHERE << "HCALPackets node failed unpacking. Skipping event." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  }

  //-------------------------- raw towers ------------------------------//
  {
    TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_CEMC");
    if (towers)
    {
      int size = towers->size();  // online towers should be the same!
      for (int channel = 0; channel < size; channel++)
      {
        TowerInfo* tower = towers->get_tower_at_channel(channel);
        unsigned int towerkey = towers->encode_key(channel);
        int ieta = towers->getTowerEtaBin(towerkey);
        int iphi = towers->getTowerPhiBin(towerkey);
        float raw_energy = tower->get_energy();
        float zs_energy = -9999;
        if (m_SimFlag)
        {
          zs_energy = cemc_sim_waveforms->get_tower_at_channel(channel)->get_waveform_value(6) - cemc_sim_waveforms->get_tower_at_channel(channel)->get_waveform_value(0);
          h_cemc_etaphi_pedestal->Fill(ieta, iphi, cemc_sim_waveforms->get_tower_at_channel(channel)->get_waveform_value(0));
        }
        else 
        {
          zs_energy = cemc_waveforms.at(channel).at(1) - cemc_waveforms.at(channel).at(0);
          h_cemc_etaphi_pedestal->Fill(ieta, iphi, cemc_waveforms.at(channel).at(0));
        }
        if (Verbosity() > 0) 
        {
          std::cout << "EMCal channel " << channel << " ieta " << ieta << " iphi " << iphi << " template E " << raw_energy << " ZS E " << zs_energy << std::endl;
        }
        if (raw_energy > m_cemc_adc_threshold && raw_energy < m_cemc_high_adc_threshold) 
        {
          h_cemc_etaphi_ZScrosscalib->Fill(ieta, iphi, zs_energy/raw_energy);
        }
      }
    }
  }
  {
    TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_HCALOUT");
    if (towers)
    {
      int size = towers->size();  // online towers should be the same!
      for (int channel = 0; channel < size; channel++)
      {
        TowerInfo* tower = towers->get_tower_at_channel(channel);
        unsigned int towerkey = towers->encode_key(channel);
        int ieta = towers->getTowerEtaBin(towerkey);
        int iphi = towers->getTowerPhiBin(towerkey);
        float raw_energy = tower->get_energy();
        float zs_energy = -9999;
        if (m_SimFlag)
        {
          zs_energy = ohcal_sim_waveforms->get_tower_at_channel(channel)->get_waveform_value(6) - ohcal_sim_waveforms->get_tower_at_channel(channel)->get_waveform_value(0);
          h_ohcal_etaphi_pedestal->Fill(ieta, iphi, ohcal_sim_waveforms->get_tower_at_channel(channel)->get_waveform_value(0));
        }
        else 
        {
          zs_energy = ohcal_waveforms.at(channel).at(1) - ohcal_waveforms.at(channel).at(0);
          h_ohcal_etaphi_pedestal->Fill(ieta, iphi, ohcal_waveforms.at(channel).at(0));
        }
        if (Verbosity() > 0) 
        {
          std::cout << "OHCal channel " << channel << " ieta " << ieta << " iphi " << iphi << " template E " << raw_energy << " ZS E " << zs_energy << std::endl;
        }
        if (raw_energy > m_hcal_adc_threshold && raw_energy < m_hcal_high_adc_threshold)
        {
          h_ohcal_etaphi_ZScrosscalib->Fill(ieta, iphi, zs_energy/raw_energy);
        }
      }
    }
  }
  {
    TowerInfoContainer* towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_HCALIN");
    if (towers)
    {
      int size = towers->size();  // online towers should be the same!
      for (int channel = 0; channel < size; channel++)
      {
        TowerInfo* tower = towers->get_tower_at_channel(channel);
        unsigned int towerkey = towers->encode_key(channel);
        int ieta = towers->getTowerEtaBin(towerkey);
        int iphi = towers->getTowerPhiBin(towerkey);
        float raw_energy = tower->get_energy();
        float zs_energy = -9999;
        if (m_SimFlag)
        {
          zs_energy = ihcal_sim_waveforms->get_tower_at_channel(channel)->get_waveform_value(6) - ihcal_sim_waveforms->get_tower_at_channel(channel)->get_waveform_value(0);
          h_ihcal_etaphi_pedestal->Fill(ieta, iphi, ihcal_sim_waveforms->get_tower_at_channel(channel)->get_waveform_value(0));
        }
        else 
        {
          zs_energy = ihcal_waveforms.at(channel).at(1) - ihcal_waveforms.at(channel).at(0);
          h_ihcal_etaphi_pedestal->Fill(ieta, iphi, ihcal_waveforms.at(channel).at(0));
        }
        if (Verbosity() > 0) 
        {
          std::cout << "IHCal channel " << channel << " ieta " << ieta << " iphi " << iphi << " template E " << raw_energy << " ZS E " << zs_energy << std::endl;
        }
        if (raw_energy > m_hcal_adc_threshold && raw_energy < m_hcal_high_adc_threshold)
        {
          h_ihcal_etaphi_ZScrosscalib->Fill(ieta, iphi, zs_energy/raw_energy);
        }
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloFittingQA::process_data(PHCompositeNode *topNode, CaloTowerDefs::DetectorSystem dettype, std::vector<std::vector<float>> &waveforms)
{

  int packet_low = std::numeric_limits<int>::min();
  int packet_high = std::numeric_limits<int>::min();
  int m_nsamples = 2;
  int m_nchannels = 192;

  if (dettype == CaloTowerDefs::CEMC)
  {
    packet_low = 6001;
    packet_high = 6128;
    m_nchannels = 192;
  }
  else if (dettype == CaloTowerDefs::HCALIN)
  {
    packet_low = 7001;
    packet_high = 7008;
    m_nchannels = 192;
  }
  else if (dettype == CaloTowerDefs::HCALOUT)
  {
    packet_low = 8001;
    packet_high = 8008;
    m_nchannels = 192;
  }
  else if (dettype == CaloTowerDefs::SEPD)
  {
    packet_low = 9001;
    packet_high = 9006;
    m_nchannels = 128;
  }
  else if (dettype == CaloTowerDefs::ZDC)
  {
    packet_low = 12001;
    packet_high = 12001;
    m_nchannels = 128;
  }

  std::variant<CaloPacketContainer*, Event*> event;
  if (m_UseOfflinePacketFlag)
  {
    CaloPacketContainer *hcalcont = findNode::getClass<CaloPacketContainer>(topNode, nodemap.find(dettype)->second);
    if (!hcalcont)
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }
    event = hcalcont;
  }
  else
  {
    Event *_event = findNode::getClass<Event>(topNode, "PRDF");
    if (_event == nullptr)
    {
      std::cout << PHWHERE << " Event not found" << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    if (_event->getEvtType() != DATAEVENT)
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    event = _event;
  }
  //since the function call on Packet and CaloPacket is the same, maybe we can use lambda?
  auto process_packet = [&](auto *packet, int pid)
  {
    if (packet)
    {
      h_packet_events->Fill(pid);
      int nchannels = packet->iValue(0, "CHANNELS");
      unsigned int adc_skip_mask = 0;

      if (dettype == CaloTowerDefs::CEMC)
      {
        adc_skip_mask = cdbttree->GetIntValue(pid, m_fieldname);
      }
      if (dettype == CaloTowerDefs::ZDC)
      {
        nchannels = m_nchannels;
      }
      if (nchannels > m_nchannels)  // packet is corrupted and reports too many channels
      {
         return Fun4AllReturnCodes::ABORTEVENT;
      }
      
      for (int channel = 0; channel < nchannels; channel++)
      {
        if (skipChannel(channel, pid, dettype))
        {
           continue;
        }
        if (dettype == CaloTowerDefs::CEMC)
        {
          if (channel % 64 == 0)
          {
            unsigned int adcboard = (unsigned int) channel / 64;
            if ((adc_skip_mask >> adcboard) & 0x1U)
            {
              for (int iskip = 0; iskip < 64; iskip++)
              {
                std::vector<float> waveform;
                waveform.reserve(m_nsamples);

                for (int samp = 0; samp < m_nsamples; samp++)
                {
                  waveform.push_back(0);
                }
                waveforms.push_back(waveform);
                waveform.clear();
              }
            }
          }
        }

        std::vector<float> waveform;
        waveform.reserve(m_nsamples);
        if (packet->iValue(channel, "SUPPRESSED"))
        {
          waveform.push_back(packet->iValue(channel, "PRE"));
          waveform.push_back(packet->iValue(channel, "POST"));
        }
        else
        {
          waveform.push_back(packet->iValue(0, channel));
          waveform.push_back(packet->iValue(6, channel));
        }
        waveforms.push_back(waveform);
        waveform.clear();
      }

      if (nchannels < m_nchannels && !(dettype == CaloTowerDefs::CEMC && adc_skip_mask < 4))
      {
        for (int channel = 0; channel < m_nchannels - nchannels; channel++)
        {
          if (skipChannel(channel, pid, dettype))
          {
            continue;
          }
          std::vector<float> waveform;
          waveform.reserve(m_nsamples);

          for (int samp = 0; samp < m_nsamples; samp++)
          {
            waveform.push_back(0);
          }
          waveforms.push_back(waveform);
          waveform.clear();
        }
      }
    }
    else  // if the packet is missing treat constitutent channels as zero suppressed
    {
      for (int channel = 0; channel < m_nchannels; channel++)
      {
        if (skipChannel(channel, pid, dettype))
        {
          continue;
        }
        std::vector<float> waveform;
        waveform.reserve(m_nsamples);
        for (int samp = 0; samp < m_nsamples; samp++)
        {
          waveform.push_back(0);
        }
        waveforms.push_back(waveform);
        waveform.clear();
      }
    }
    return Fun4AllReturnCodes::EVENT_OK;
  };


  for (int pid = packet_low; pid <= packet_high; pid++)
  {
    if (auto hcalcont = std::get_if<CaloPacketContainer *>(&event))
    {
      CaloPacket *packet = (*hcalcont)->getPacketbyId(pid);
      if(process_packet(packet, pid) == Fun4AllReturnCodes::ABORTEVENT)
      {
        return Fun4AllReturnCodes::ABORTEVENT;
      }
    }
    else if (auto _event = std::get_if<Event *>(&event))
    {
      Packet *packet = (*_event)->getPacket(pid);
      if(process_packet(packet, pid) == Fun4AllReturnCodes::ABORTEVENT)
      {
        //I think it is safe to delete a nullptr...
        delete packet;
        return Fun4AllReturnCodes::ABORTEVENT;
      }
      else
      {
      delete packet;
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;

}

bool CaloFittingQA::skipChannel(int ich, int pid, CaloTowerDefs::DetectorSystem dettype)
{
  if (dettype == CaloTowerDefs::SEPD)
  {
    int sector = ((ich + 1) / 32);
    int emptych = -999;
    if ((sector == 0) && (pid == 9001))
    {
      emptych = 1;
    }
    else
    {
      emptych = 14 + 32 * sector;
    }
    if (ich == emptych)
    {
      return true;
    }
  }
    
  if (dettype == CaloTowerDefs::ZDC)
  {
     if(((ich > 17) && (ich < 48)) || ((ich > 63) && (ich < 80)) || ((ich > 81) && (ich < 112)))
     {
        return true;
     }
  }
    
  return false;
}

int CaloFittingQA::End(PHCompositeNode* /*unused*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

std::string CaloFittingQA::getHistoPrefix() const { return std::string("h_") + Name() + std::string("_"); }

void CaloFittingQA::createHistos()
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // create and register your histos (all types) here
  h_cemc_etaphi_ZScrosscalib = new TProfile2D(boost::str(boost::format("%scemc_etaphi_ZScrosscalib") % getHistoPrefix()).c_str(), ";eta;phi", 96, 0, 96, 256, 0, 256, -10, 10);
  h_cemc_etaphi_ZScrosscalib->SetDirectory(nullptr);
  hm->registerHisto(h_cemc_etaphi_ZScrosscalib);

  h_ihcal_etaphi_ZScrosscalib = new TProfile2D(boost::str(boost::format("%sihcal_etaphi_ZScrosscalib") % getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, -10, 10);
  h_ihcal_etaphi_ZScrosscalib->SetDirectory(nullptr);
  hm->registerHisto(h_ihcal_etaphi_ZScrosscalib);

  h_ohcal_etaphi_ZScrosscalib = new TProfile2D(boost::str(boost::format("%sohcal_etaphi_ZScrosscalib") % getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, -10, 10);
  h_ohcal_etaphi_ZScrosscalib->SetDirectory(nullptr);
  hm->registerHisto(h_ohcal_etaphi_ZScrosscalib);

  h_cemc_etaphi_pedestal = new TProfile2D(boost::str(boost::format("%scemc_etaphi_pedestal") % getHistoPrefix()).c_str(), ";eta;phi", 96, 0, 96, 256, 0, 256, 0, 16400);
  h_cemc_etaphi_pedestal->SetErrorOption("G");
  h_cemc_etaphi_pedestal->SetDirectory(nullptr);
  hm->registerHisto(h_cemc_etaphi_pedestal);

  h_ihcal_etaphi_pedestal = new TProfile2D(boost::str(boost::format("%sihcal_etaphi_pedestal") % getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, 0, 16400);
  h_ihcal_etaphi_pedestal->SetErrorOption("G");
  h_ihcal_etaphi_pedestal->SetDirectory(nullptr);
  hm->registerHisto(h_ihcal_etaphi_pedestal);

  h_ohcal_etaphi_pedestal = new TProfile2D(boost::str(boost::format("%sohcal_etaphi_pedestal") % getHistoPrefix()).c_str(), ";eta;phi", 24, 0, 24, 64, 0, 64, 0, 16400);
  h_ohcal_etaphi_pedestal->SetErrorOption("G");
  h_ohcal_etaphi_pedestal->SetDirectory(nullptr);
  hm->registerHisto(h_ohcal_etaphi_pedestal);

  h_packet_events = new TH1I(boost::str(boost::format("%spacket_events") % getHistoPrefix()).c_str(), ";packet id", 6010, 6000, 12010);
  h_packet_events->SetDirectory(nullptr);
  hm->registerHisto(h_packet_events);

}
