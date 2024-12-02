#include "CaloTowerBuilder.h"
#include "CaloTowerDefs.h"

#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfoContainerv2.h>
#include <calobase/TowerInfoContainerv3.h>
#include <calobase/TowerInfoContainerv4.h>
#include <calobase/TowerInfoContainerSimv1.h>
#include <calobase/TowerInfoContainerSimv2.h>

#include <ffarawobjects/CaloPacket.h>
#include <ffarawobjects/CaloPacketContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco
#include <phool/recoConsts.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>    // for PHIODataNode
#include <phool/PHNode.h>          // for PHNode
#include <phool/PHNodeIterator.h>  // for PHNodeIterator
#include <phool/PHObject.h>        // for PHObject
#include <phool/getClass.h>

#include <cdbobjects/CDBTTree.h>  // for CDBTTree

#include <ffamodules/CDBInterface.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>

#include <TSystem.h>

#include <climits>
#include <variant>
#include <iostream>  // for operator<<, endl, basic...
#include <memory>    // for allocator_traits<>::val...
#include <vector>    // for vector

static const std::map<CaloTowerDefs::DetectorSystem, std::string> nodemap{
    {CaloTowerDefs::CEMC, "CEMCPackets"},
    {CaloTowerDefs::HCALIN, "HCALPackets"},
    {CaloTowerDefs::HCALOUT, "HCALPackets"},
    {CaloTowerDefs::ZDC, "ZDCPackets"},
    {CaloTowerDefs::SEPD, "SEPDPackets"}};
//____________________________________________________________________________..
CaloTowerBuilder::CaloTowerBuilder(const std::string &name)
  : SubsysReco(name)
{
  WaveformProcessing = new CaloWaveformProcessing();
}

//____________________________________________________________________________..
CaloTowerBuilder::~CaloTowerBuilder()
{
  delete cdbttree;
  delete cdbttree_tbt_zs;
  delete WaveformProcessing;
}

//____________________________________________________________________________..
int CaloTowerBuilder::InitRun(PHCompositeNode *topNode)
{
  WaveformProcessing->set_processing_type(_processingtype);
  WaveformProcessing->set_softwarezerosuppression(m_bdosoftwarezerosuppression, m_nsoftwarezerosuppression);
  if (m_setTimeLim)
  {
     WaveformProcessing->set_timeFitLim(m_timeLim_low,m_timeLim_high);
  }
  if (m_dobitfliprecovery)
  {
    WaveformProcessing->set_bitFlipRecovery(m_dobitfliprecovery);
  }

  if (m_dettype == CaloTowerDefs::CEMC)
  {
    m_detector = "CEMC";
    m_packet_low = 6001;
    m_packet_high = 6128;
    m_nchannels = 192;
    WaveformProcessing->set_template_name("CEMC_TEMPLATE");
    if (_processingtype == CaloWaveformProcessing::NONE)
    {
      WaveformProcessing->set_processing_type(CaloWaveformProcessing::TEMPLATE);
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
  }
  else if (m_dettype == CaloTowerDefs::HCALIN)
  {
    m_packet_low = 7001;
    m_packet_high = 7008;
    m_detector = "HCALIN";
    m_nchannels = 192;
    WaveformProcessing->set_template_name("IHCAL_TEMPLATE");
    if (_processingtype == CaloWaveformProcessing::NONE)
    {
      WaveformProcessing->set_processing_type(CaloWaveformProcessing::TEMPLATE);
    }
  }
  else if (m_dettype == CaloTowerDefs::HCALOUT)
  {
    m_detector = "HCALOUT";
    m_packet_low = 8001;
    m_packet_high = 8008;
    m_nchannels = 192;
    WaveformProcessing->set_template_name("OHCAL_TEMPLATE");
    if (_processingtype == CaloWaveformProcessing::NONE)
    {
      WaveformProcessing->set_processing_type(CaloWaveformProcessing::TEMPLATE);
    }
  }
  else if (m_dettype == CaloTowerDefs::SEPD)
  {
    m_detector = "SEPD";
    m_packet_low = 9001;
    m_packet_high = 9006;
    m_nchannels = 128;
    if (_processingtype == CaloWaveformProcessing::NONE)
    {
      WaveformProcessing->set_processing_type(CaloWaveformProcessing::FAST);  // default the EPD to fast processing
    }
  }
  else if (m_dettype == CaloTowerDefs::ZDC)
  {
    m_detector = "ZDC";
    m_packet_low = 12001;
    m_packet_high = 12001;
    m_nchannels = 128;
    if (_processingtype == CaloWaveformProcessing::NONE)
    {
      WaveformProcessing->set_processing_type(CaloWaveformProcessing::FAST);  // default the ZDC to fast processing
    }
  }
  WaveformProcessing->initialize_processing();
  if(m_dotbtszs)
  {
    cdbttree_tbt_zs = new CDBTTree(m_zsURL.c_str());
  }

  CreateNodeTree(topNode);
  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloTowerBuilder::process_sim()
{
  std::vector<std::vector<float>> waveforms;

  for (int ich = 0; ich < (int) m_CalowaveformContainer->size(); ich++)
  {
    TowerInfo *towerinfo = m_CalowaveformContainer->get_tower_at_channel(ich);
    std::vector<float> waveform;
    waveform.reserve(m_nsamples);
    bool fillwaveform = true;
    //get key
    if(m_dotbtszs)
    {
      unsigned int key = m_CalowaveformContainer->encode_key(ich);
      int zs_threshold = cdbttree_tbt_zs->GetIntValue(key, m_zs_fieldname);
      int pre = towerinfo->get_waveform_value(0);
      //this is always safe since towerinfo v3 has 31 samples
      int post = towerinfo->get_waveform_value(6);
      if((post - pre) <= zs_threshold)
      {
        //zero suppressed
        fillwaveform = false;
        waveform.push_back(pre);
        waveform.push_back(post);
      }
      
    }
    if(fillwaveform)
    {
      for (int samp = 0; samp < m_nsamples; samp++)
      {
        waveform.push_back(towerinfo->get_waveform_value(samp));
      }
    }
    waveforms.push_back(waveform);
    waveform.clear();
  }

  std::vector<std::vector<float>> processed_waveforms = WaveformProcessing->process_waveform(waveforms);
  int n_channels = processed_waveforms.size();
  for (int i = 0; i < n_channels; i++)
  {
    //this is for copying the truth info to the downstream object
    TowerInfo* towerwaveform = m_CalowaveformContainer->get_tower_at_channel(i);
    TowerInfo *towerinfo = m_CaloInfoContainer->get_tower_at_channel(i);
    towerinfo->copy_tower(towerwaveform);
    towerinfo->set_time(processed_waveforms.at(i).at(1));
    towerinfo->set_energy(processed_waveforms.at(i).at(0));
    towerinfo->set_time_float(processed_waveforms.at(i).at(1));
    towerinfo->set_pedestal(processed_waveforms.at(i).at(2));
    towerinfo->set_chi2(processed_waveforms.at(i).at(3));
    bool SZS = isSZS(processed_waveforms.at(i).at(1), processed_waveforms.at(i).at(3));
    if (processed_waveforms.at(i).at(4) == 0) 
    {
      towerinfo->set_isRecovered(false);
    }
    else
    {
      towerinfo->set_isRecovered(true);
    }
    int n_samples = waveforms.at(i).size();
    if (n_samples == m_nzerosuppsamples || SZS)
    {
      towerinfo->set_isZS(true);
    }
    for (int j = 0; j < n_samples; j++)
    {
      towerinfo->set_waveform_value(j, waveforms.at(i).at(j));
    }
  }
  waveforms.clear();

  return Fun4AllReturnCodes::EVENT_OK;
}



int CaloTowerBuilder::process_data(PHCompositeNode *topNode, std::vector<std::vector<float>> &waveforms)
{
  std::variant<CaloPacketContainer*, Event*> event;
  if (m_UseOfflinePacketFlag)
  {
    CaloPacketContainer *hcalcont = findNode::getClass<CaloPacketContainer>(topNode, nodemap.find(m_dettype)->second);
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
      int nchannels = packet->iValue(0, "CHANNELS");
      unsigned int adc_skip_mask = 0;

      if (m_dettype == CaloTowerDefs::CEMC)
      {
        adc_skip_mask = cdbttree->GetIntValue(pid, m_fieldname);
      }
      if (m_dettype == CaloTowerDefs::ZDC)
      {
        nchannels = m_nchannels;
      }
      if (nchannels > m_nchannels)  // packet is corrupted and reports too many channels
      {
         return Fun4AllReturnCodes::ABORTEVENT;
      }
      
      for (int channel = 0; channel < nchannels; channel++)
      {
        if (skipChannel(channel, pid))
        {
           continue;
        }
        if (m_dettype == CaloTowerDefs::CEMC)
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

                for (int samp = 0; samp < m_nzerosuppsamples; samp++)
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
          for (int samp = 0; samp < m_nsamples; samp++)
          {
            waveform.push_back(packet->iValue(samp, channel));
          }
        }
        waveforms.push_back(waveform);
        waveform.clear();
      }

      if (nchannels < m_nchannels && !(m_dettype == CaloTowerDefs::CEMC && adc_skip_mask < 4))
      {
        for (int channel = 0; channel < m_nchannels - nchannels; channel++)
        {
          if (skipChannel(channel, pid))
          {
            continue;
          }
          std::vector<float> waveform;
          waveform.reserve(m_nsamples);

          for (int samp = 0; samp < m_nzerosuppsamples; samp++)
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
        if (skipChannel(channel, pid))
        {
          continue;
        }
        std::vector<float> waveform;
        waveform.reserve(2);
        for (int samp = 0; samp < m_nzerosuppsamples; samp++)
        {
          waveform.push_back(0);
        }
        waveforms.push_back(waveform);
        waveform.clear();
      }
    }
    return Fun4AllReturnCodes::EVENT_OK;
  };


  for (int pid = m_packet_low; pid <= m_packet_high; pid++)
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
  //____________________________________________________________________________..
int CaloTowerBuilder::process_event(PHCompositeNode *topNode)
{
  if (!m_isdata)
  {
    return process_sim();
  }
  std::vector<std::vector<float>> waveforms;
  if(process_data(topNode, waveforms) == Fun4AllReturnCodes::ABORTEVENT)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  // waveform vector is filled here, now fill our output. methods from the base class make sure
  // we only fill what the chosen container version supports
  std::vector<std::vector<float>> processed_waveforms = WaveformProcessing->process_waveform(waveforms);
  int n_channels = processed_waveforms.size();
  for (int i = 0; i < n_channels; i++)
  {
    TowerInfo *towerinfo = m_CaloInfoContainer->get_tower_at_channel(i);
    towerinfo->set_time(processed_waveforms.at(i).at(1));
    towerinfo->set_energy(processed_waveforms.at(i).at(0));
    towerinfo->set_time_float(processed_waveforms.at(i).at(1));
    towerinfo->set_pedestal(processed_waveforms.at(i).at(2));
    towerinfo->set_chi2(processed_waveforms.at(i).at(3));
    bool SZS = isSZS(processed_waveforms.at(i).at(1), processed_waveforms.at(i).at(3));
    if (processed_waveforms.at(i).at(4) == 0) 
    {
      towerinfo->set_isRecovered(false);
    }
    else
    {
      towerinfo->set_isRecovered(true);
    }
    int n_samples = waveforms.at(i).size();
    if (n_samples == m_nzerosuppsamples || SZS)
    {
      if(waveforms.at(i).at(0) == 0)
      {
        towerinfo->set_isNotInstr(true);
      }
      else
      {
        towerinfo->set_isZS(true);
      }
    }
    
    for (int j = 0; j < n_samples; j++)
    {
      towerinfo->set_waveform_value(j, waveforms.at(i).at(j));
    }
  }
  waveforms.clear();

  return Fun4AllReturnCodes::EVENT_OK;
}

bool CaloTowerBuilder::skipChannel(int ich, int pid)
{
  if (m_dettype == CaloTowerDefs::SEPD)
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
    
  if (m_dettype == CaloTowerDefs::ZDC)
  {
     if(((ich > 17) && (ich < 48)) || ((ich > 63) && (ich < 80)) || ((ich > 81) && (ich < 112)))
     {
        return true;
     }
  }
    
  return false;
}

bool CaloTowerBuilder::isSZS(float time, float chi2)
{
  //isfinite
  if (!std::isfinite(time) && !std::isfinite(chi2))
  {
    return true;
  }
  return false;
}

void CaloTowerBuilder::CreateNodeTree(PHCompositeNode *topNode)
{
  PHNodeIterator topNodeItr(topNode);
  // DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(topNodeItr.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << "PHComposite node created: DST" << std::endl;
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }
  if (!m_isdata)
  {
    std::string waveformNodeName = m_inputNodePrefix + m_detector;
    m_CalowaveformContainer = findNode::getClass<TowerInfoContainer>(topNode, waveformNodeName);
    if (!m_CalowaveformContainer)
    {
      std::cout << PHWHERE << "simulation waveform container " << waveformNodeName << " not found" << std::endl;
      gSystem->Exit(1);
      exit(1);
    }
  }

  // towers
  PHNodeIterator nodeItr(dstNode);
  PHCompositeNode *DetNode;
  // enum CaloTowerDefs::DetectorSystem and TowerInfoContainer::DETECTOR are different!!!!
  TowerInfoContainer::DETECTOR DetectorEnum = TowerInfoContainer::DETECTOR::DETECTOR_INVALID;
  std::string DetectorNodeName;
  if (m_dettype == CaloTowerDefs::CEMC)
  {
    DetectorEnum = TowerInfoContainer::DETECTOR::EMCAL;
    DetectorNodeName = "CEMC";
  }
  else if (m_dettype == CaloTowerDefs::SEPD)
  {
    DetectorEnum = TowerInfoContainer::DETECTOR::SEPD;
    DetectorNodeName = "SEPD";
  }
  else if (m_dettype == CaloTowerDefs::ZDC)
  {
    DetectorEnum = TowerInfoContainer::DETECTOR::ZDC;
    DetectorNodeName = "ZDC";
  }
  else if (m_dettype == CaloTowerDefs::HCALIN)
  {
    DetectorEnum = TowerInfoContainer::DETECTOR::HCAL;
    DetectorNodeName = "HCALIN";
  }
  else if (m_dettype == CaloTowerDefs::HCALOUT)
  {
    DetectorEnum = TowerInfoContainer::DETECTOR::HCAL;
    DetectorNodeName = "HCALOUT";
  }
  else
  {
    std::cout << PHWHERE << " Invalid detector type " << m_dettype << std::endl;
    gSystem->Exit(1);
    exit(1);
  }
  DetNode = dynamic_cast<PHCompositeNode *>(nodeItr.findFirst("PHCompositeNode", DetectorNodeName));
  if (!DetNode)
  {
    DetNode = new PHCompositeNode(DetectorNodeName);
    dstNode->addNode(DetNode);
  }
  if (m_buildertype == CaloTowerDefs::kPRDFTowerv1)
  {
    m_CaloInfoContainer = new TowerInfoContainerv1(DetectorEnum);
  }
  else if (m_buildertype == CaloTowerDefs::kPRDFWaveform)
  {
    m_CaloInfoContainer = new TowerInfoContainerv3(DetectorEnum);
  }
  else if (m_buildertype == CaloTowerDefs::kWaveformTowerv2)
  {
    m_CaloInfoContainer = new TowerInfoContainerv2(DetectorEnum);
  }
  else if (m_buildertype == CaloTowerDefs::kPRDFTowerv4)
  {
    m_CaloInfoContainer = new TowerInfoContainerv4(DetectorEnum);
  }
  else if (m_buildertype == CaloTowerDefs::kWaveformTowerSimv1)
  {
    m_CaloInfoContainer = new TowerInfoContainerSimv1(DetectorEnum);
  }
  else
  {
    std::cout << PHWHERE << "invalid builder type " << m_buildertype << std::endl;
    gSystem->Exit(1);
    exit(1);
  }
  TowerNodeName = m_outputNodePrefix + m_detector;
  PHIODataNode<PHObject> *newTowerNode = new PHIODataNode<PHObject>(m_CaloInfoContainer, TowerNodeName, "PHObject");
  DetNode->addNode(newTowerNode);
}
