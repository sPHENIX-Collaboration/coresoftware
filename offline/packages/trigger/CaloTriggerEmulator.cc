#include "CaloTriggerEmulator.h"

#include "LL1Defs.h"
#include "LL1Outv1.h"
#include "TriggerPrimitiveContainerv1.h"
#include "TriggerPrimitivev1.h"

#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv3.h>
#include <calobase/TowerInfoDefs.h>
#include <calobase/TowerInfov3.h>

#include <ffamodules/CDBInterface.h>

#include <ffarawobjects/CaloPacket.h>
#include <ffarawobjects/CaloPacketContainer.h>
#include <ffarawobjects/CaloPacketContainerv1.h>
#include <ffarawobjects/CaloPacketv1.h>

#include <cdbobjects/CDBHistos.h>  // for CDBHistos
#include <cdbobjects/CDBTTree.h>   // for CDBHistos

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>

#include <TFile.h>
#include <TH1.h>
#include <TH1I.h>
#include <TNtuple.h>

#include <bitset>
#include <cassert>
#include <cstdint>
#include <sstream>
#include <string>

// constructor
CaloTriggerEmulator::CaloTriggerEmulator(const std::string &name)
  : SubsysReco(name)
  , m_trigger("NONE")
{
  // initialize all important counters

  m_nevent = 0;
  m_jet_npassed = 0;
  m_pair_npassed = 0;
  m_photon_npassed = 0;
  // is data flag is not used right now
  m_isdata = 1;
  // default nsamples is 16 for mbd, 12 for calos
  m_nsamples = 16;

  // for MBD, this is the peak sample in run-23 data
  m_idx = 12;

  // default values for the lookup tables.
  // TODO: to CDB the LUTs from the database

  // reset variables
  for (unsigned int i = 0; i < 1024; i++)
  {
    m_l1_8x8_table[i] = i;
    if (i >= 0xff)
    {
      m_l1_8x8_table[i] = 0xff;
    }
  }

  for (unsigned int i = 0; i < 1024; i++)
  {
    m_l1_adc_table[i] = (i) &0x3ffU;
  }

  for (unsigned int i = 0; i < 4096; i++)
  {
    m_l1_slewing_table[i] = (i) &0x3ffU;
  }

  for (int i = 0; i < 24576; i++)
  {
    unsigned int key = TowerInfoDefs::encode_emcal(i);
    h_emcal_lut[key] = nullptr;
  }

  for (int i = 0; i < 1536; i++)
  {
    unsigned int key = TowerInfoDefs::encode_hcal(i);
    h_hcalin_lut[key] = nullptr;
  }
  for (int i = 0; i < 1536; i++)
  {
    unsigned int key = TowerInfoDefs::encode_hcal(i);
    h_hcalout_lut[key] = nullptr;
  }

  // Set HCAL LL1 lookup table for the cosmic coincidence trigger.
  if (m_triggerid == TriggerDefs::TriggerId::cosmic_coinTId)
  {
    unsigned int bits1, bits2, sumbits1, sumbits2;
    for (unsigned int i = 0; i < 4096; i++)
    {
      sumbits1 = 0;
      sumbits2 = 0;

      bits1 = (i & 0x3fU);
      bits2 = ((i >> 6U) & 0x3fU);
      for (unsigned int j = 0; j < 3; j++)
      {
        if (((bits1 >> j) & 0x1U) && ((bits2 >> j) & 0x1U))
        {
          sumbits1++;
        }
        if (((bits1 >> (j + 3U)) & 0x1U) && ((bits2 >> (j + 3U)) & 0x1U))
        {
          sumbits2++;
        }
      }

      m_l1_hcal_table[i] = 0;
      if (i == 0)
      {
        continue;
      }
      m_l1_hcal_table[i] |= (sumbits1 ? 0x1U : 0U);
      m_l1_hcal_table[i] |= (sumbits2 ? 0x2U : 0U);
    }
  }

  else if (m_triggerid == TriggerDefs::TriggerId::cosmicTId)
  {
    unsigned int bits1, bits2, sumbits1, sumbits2;
    for (unsigned int i = 0; i < 4096; i++)
    {
      sumbits1 = 0;
      sumbits2 = 0;

      bits1 = (i & 0x3fU);
      bits2 = ((i >> 6U) & 0x3fU);
      for (unsigned int j = 0; j < 6; j++)
      {
        sumbits1 += ((bits1 >> j) & 0x1U);
        sumbits2 += ((bits2 >> j) & 0x1U);
      }

      m_l1_hcal_table[i] = 0;
      if (i == 0)
      {
        continue;
      }
      m_l1_hcal_table[i] |= (sumbits1 ? 0x1U : 0U);
      m_l1_hcal_table[i] |= (sumbits2 ? 0x2U : 0U);
    }
  }

  // point to null for all of these objects to be added to or grabbed from the node tree.
  // this will hold the ll1 info that goes through the emulator
  m_ll1out_pair = nullptr;
  m_ll1out_photon = nullptr;
  m_ll1out_jet = nullptr;

  // waveform containers to be grabbed from node tree.
  // Done int the CaloPacketGetter

  // _waveforms_emcal = nullptr;
  // _waveforms_hcalin = nullptr;
  // _waveforms_hcalout = nullptr;
  // _waveforms_mbd = nullptr;

  // // to hold the primitives constructed from the waveforms.
  // _primitiveprimitives = nullptr;
  // _primitives_emcal = nullptr;
  // _primitives_hcalin = nullptr;
  // _primitives_hcalout = nullptr;
  // _primitives_emcal_ll1 = nullptr;
  // _primitives_hcal_ll1 = nullptr;

  // _primitive = nullptr;

  m_n_primitives = 0;
  m_n_sums = 16;
  m_trig_sample = -1;
  m_trig_sub_delay = 4;

  // define a detector map for detectors included in a trigger
  m_det_map[TriggerDefs::TriggerId::noneTId] = {};
  m_det_map[TriggerDefs::TriggerId::jetTId] = {"EMCAL", "HCALIN", "HCALOUT"};
  m_det_map[TriggerDefs::TriggerId::cosmicTId] = {"HCALIN", "HCALOUT"};
  m_det_map[TriggerDefs::TriggerId::cosmic_coinTId] = {"HCALIN", "HCALOUT"};
  m_det_map[TriggerDefs::TriggerId::pairTId] = {"EMCAL"};
  m_det_map[TriggerDefs::TriggerId::photonTId] = {"EMCAL", "HCALIN", "HCALOUT"};

  // define primitive map as number of primitives in a detector.
  m_prim_map[TriggerDefs::DetectorId::noneDId] = 0;
  m_prim_map[TriggerDefs::DetectorId::emcalDId] = 384;
  m_prim_map[TriggerDefs::DetectorId::hcalinDId] = 24;
  m_prim_map[TriggerDefs::DetectorId::hcaloutDId] = 24;
  m_prim_map[TriggerDefs::DetectorId::hcalDId] = 24;

  m_prim_ll1_map[TriggerDefs::TriggerId::jetTId] = 16;
  m_prim_ll1_map[TriggerDefs::TriggerId::pairTId] = 16;
  m_prim_ll1_map[TriggerDefs::TriggerId::photonTId] = 384;
  // booleans to control the input of detector data
  m_do_emcal = false;
  m_do_hcalin = false;
  m_do_hcalout = false;

  m_masks_channel = {};  //, 70385703};
  m_masks_fiber = {};    //, 70385696};
}

// check whether a channel has been masked
bool CaloTriggerEmulator::CheckChannelMasks(TriggerDefs::TriggerSumKey key)
{
  for (unsigned int &it : m_masks_channel)
  {
    if (key == it)
    {
      return true;
    }
  }
  return false;
}

// cehck whether a fiber has been masked
bool CaloTriggerEmulator::CheckFiberMasks(TriggerDefs::TriggerPrimKey key)
{
  return (std::find(m_masks_fiber.begin(), m_masks_fiber.end(), key) != m_masks_fiber.end());
}

void CaloTriggerEmulator::LoadFiberMasks()
{
  TFile *fin = new TFile(m_optmask_file.c_str(), "r");
  TNtuple *tn_keys = (TNtuple *) fin->Get("tn_optmask");
  float key;
  tn_keys->SetBranchAddress("primkey", &key);
  for (int i = 0; i < tn_keys->GetEntries(); i++)
  {
    tn_keys->GetEntry(i);
    unsigned int primkey = static_cast<unsigned int>((int) key);
    m_masks_fiber.push_back(primkey);
  }

  fin->Close();
  delete fin;
}

// setting the trigger type
void CaloTriggerEmulator::setTriggerType(const std::string &name)
{
  m_trigger = name;
  m_triggerid = TriggerDefs::GetTriggerId(m_trigger);
}

// setting the trigger typ
void CaloTriggerEmulator::setTriggerType(TriggerDefs::TriggerId triggerid)
{
  m_triggerid = triggerid;
}

// make file and histomanager (but nothing goes in the file at the moment)
int CaloTriggerEmulator::Init(PHCompositeNode * /*topNode*/)
{
  return 0;
}

// at the beginning of the run
int CaloTriggerEmulator::InitRun(PHCompositeNode *topNode)
{
  // Get the detectors that are used for a given trigger.

  if (m_triggerid == TriggerDefs::TriggerId::jetTId)
  {
    if (Verbosity() >= 2)
    {
      std::cout << "Using Jet Trigger." << std::endl;
    }
    if (!m_force_emcal)
    {
      m_do_emcal = true;
    }
    if (!m_force_hcalin)
    {
      m_do_hcalin = true;
    }
    if (!m_force_hcalout)
    {
      m_do_hcalout = true;
    }
  }
  else if (m_triggerid == TriggerDefs::TriggerId::photonTId)
  {
    if (Verbosity() >= 2)
    {
      std::cout << "Using Photon Trigger." << std::endl;
    }
    if (!m_force_emcal)
    {
      m_do_emcal = true;
    }
    if (!m_force_hcalin)
    {
      m_do_hcalin = false;
    }
    if (!m_force_hcalout)
    {
      m_do_hcalout = false;
    }

    if (m_do_emcal && (m_do_hcalin || m_do_hcalout))
    {
      std::cout << "Cannot run PHOTON with hcal and emcal" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  }
  else if (m_triggerid == TriggerDefs::TriggerId::cosmicTId)
  {
    if (Verbosity() >= 2)
    {
      std::cout << "Using Cosmic Trigger." << std::endl;
    }
    if (!m_force_emcal)
    {
      m_do_emcal = false;
    }
    if (!m_force_hcalin)
    {
      m_do_hcalin = true;
    }
    if (!m_force_hcalout)
    {
      m_do_hcalout = true;
    }
  }
  else if (m_triggerid == TriggerDefs::TriggerId::cosmic_coinTId)
  {
    if (Verbosity() >= 2)
    {
      std::cout << "Using Cosmic Coincidence Trigger." << std::endl;
    }
    if (!m_force_emcal)
    {
      m_do_emcal = false;
    }
    if (!m_force_hcalin)
    {
      m_do_hcalin = true;
    }
    if (!m_force_hcalout)
    {
      m_do_hcalout = true;
    }
  }
  else if (m_triggerid == TriggerDefs::TriggerId::pairTId)
  {
    if (Verbosity() >= 2)
    {
      std::cout << "Using Pair Trigger." << std::endl;
    }
    if (!m_force_emcal)
    {
      m_do_emcal = true;
    }
    if (!m_force_hcalin)
    {
      m_do_hcalin = false;
    }
    if (!m_force_hcalout)
    {
      m_do_hcalout = false;
    }
  }
  else
  {
    std::cout << __FUNCTION__ << " : No trigger selected using ALL" << std::endl;

    m_trigger = "PHYSICS";
    m_triggerid = TriggerDefs::TriggerId::physicsTId;

    if (!m_force_emcal)
    {
      m_do_emcal = true;
    }
    if (!m_force_hcalin)
    {
      m_do_hcalin = true;
    }
    if (!m_force_hcalout)
    {
      m_do_hcalout = true;
    }
  }

  m_ll1_nodename = "LL1OUT_" + m_trigger;
  m_prim_nodename = "TRIGGERPRIMITIVES_" + m_trigger;

  // Get the calibrations and proroceduce the lookup tables;

  if (Download_Calibrations())
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  CreateNodes(topNode);

  return 0;
}

int CaloTriggerEmulator::Download_Calibrations()
{
  if (m_do_emcal)
  {
    std::string calibName = "cemc_adcskipmask";

    std::string calibdir = CDBInterface::instance()->getUrl(calibName);
    if (calibdir.empty())
    {
      std::cout << PHWHERE << "ADC Skip mask not found in CDB, not even in the default... " << std::endl;
      exit(1);
    }
    cdbttree_adcmask = new CDBTTree(calibdir.c_str());
  }
  if (!m_optmask_file.empty())
  {
    LoadFiberMasks();
  }

  if (m_do_emcal && !m_default_lut_emcal)
  {
    if (!m_emcal_lutname.empty())
    {
      cdbttree_emcal = new CDBHistos(m_emcal_lutname);
    }
    else
    {
      std::string calibdir = CDBInterface::instance()->getUrl("emcal_trigger_lut");
      if (calibdir.empty())
      {
        m_default_lut_emcal = true;
        std::cout << "Could not find and load histograms for EMCAL LUTs! defaulting to the identity table!" << std::endl;
      }
      else
      {
        cdbttree_emcal = new CDBHistos(calibdir);
      }
    }

    if (cdbttree_emcal)
    {
      cdbttree_emcal->LoadCalibrations();

      for (int i = 0; i < 24576; i++)
      {
        std::string histoname = "h_emcal_lut_" + std::to_string(i);
        unsigned int key = TowerInfoDefs::encode_emcal(i);
        h_emcal_lut[key] = (TH1I *) cdbttree_emcal->getHisto(histoname.c_str());
      }
    }
  }
  if (m_do_hcalin && !m_default_lut_hcalin)
  {
    if (!m_hcalin_lutname.empty())
    {
      cdbttree_hcalin = new CDBHistos(m_hcalin_lutname);
    }
    else
    {
      std::string calibdir = CDBInterface::instance()->getUrl("hcalin_trigger_lut");
      if (calibdir.empty())
      {
        m_default_lut_hcalin = true;
        std::cout << "Could not find and load histograms for HCALIN LUTs! defaulting to the identity table!" << std::endl;
      }
      else
      {
        cdbttree_hcalin = new CDBHistos(calibdir);
      }
    }
    if (cdbttree_hcalin)
    {
      cdbttree_hcalin->LoadCalibrations();

      for (int i = 0; i < 1536; i++)
      {
        std::string histoname = "h_hcalin_lut_" + std::to_string(i);
        unsigned int key = TowerInfoDefs::encode_hcal(i);
        h_hcalin_lut[key] = (TH1I *) cdbttree_hcalin->getHisto(histoname.c_str());
      }
    }
  }
  if (m_do_hcalout && !m_default_lut_hcalout)
  {
    if (!m_hcalout_lutname.empty())
    {
      cdbttree_hcalout = new CDBHistos(m_hcalout_lutname);
    }
    else
    {
      std::string calibdir = CDBInterface::instance()->getUrl("hcalout_trigger_lut");
      if (calibdir.empty())
      {
        m_default_lut_hcalout = true;
        std::cout << "Could not find and load histograms for HCALOUT LUTs! defaulting to the identity table!" << std::endl;
      }
      else
      {
        cdbttree_hcalout = new CDBHistos(calibdir);
      }
    }
    if (cdbttree_hcalout)
    {
      cdbttree_hcalout->LoadCalibrations();

      for (int i = 0; i < 1536; i++)
      {
        std::string histoname = "h_hcalout_lut_" + std::to_string(i);
        unsigned int key = TowerInfoDefs::encode_hcal(i);
        h_hcalout_lut[key] = (TH1I *) cdbttree_hcalout->getHisto(histoname.c_str());
      }
    }
  }
  return 0;
}
// process event procedure
int CaloTriggerEmulator::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() >= 1)
  {
    std::cout << __FUNCTION__ << ": event " << m_nevent << std::endl;
  }

  // Get all nodes needed fo
  GetNodes(topNode);

  // process waveforms from the waveform container into primitives
  if (process_waveforms())
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  if (Verbosity())
  {
    std::cout << __FILE__ << __FUNCTION__ << __LINE__ << "::"
              << "done with waveforms" << std::endl;
  }
  // process all the primitives into sums.
  process_primitives();

  // calculate the true LL1 trigger at emcal and hcal.
  if (process_organizer())
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  // calculate the true LL1 trigger algorithm.
  if (process_trigger())
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  m_nevent++;

  if (Verbosity() >= 2)
  {
    identify();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

// RESET event procedure that takes all variables to 0 and clears the primitives.
int CaloTriggerEmulator::ResetEvent(PHCompositeNode * /*topNode*/)
{
  // here, the peak minus pedestal map is cleanly disposed of
  m_peak_sub_ped_emcal.clear();
  m_peak_sub_ped_hcalin.clear();
  m_peak_sub_ped_hcalout.clear();

  return 0;
}
int CaloTriggerEmulator::process_offline()
{
  int sample_start = 1;
  int sample_end = m_nsamples;
  if (m_trig_sample > 0)
  {
    sample_start = m_trig_sample;
    sample_end = m_trig_sample + 1;
  }

  if (m_do_emcal)
  {
    if (Verbosity())
    {
      std::cout << __FILE__ << "::" << __FUNCTION__ << ":: emcal" << std::endl;
    }

    if (!m_emcal_packets)
    {
      return Fun4AllReturnCodes::ABORTRUN;
    }
    unsigned int iwave = 0;

    for (int pid = m_packet_low_emcal; pid <= m_packet_high_emcal; pid++)
    {
      CaloPacket *packet = m_emcal_packets->getPacketbyId(pid);
      if (packet)
      {
        int nchannels = packet->iValue(0, "CHANNELS");
        unsigned int adc_skip_mask = 0;

        adc_skip_mask = cdbttree_adcmask->GetIntValue(pid, m_fieldname);

        for (int channel = 0; channel < nchannels; channel++)
        {
          if (channel % 64 == 0)
          {
            unsigned int adcboard = (unsigned int) channel / 64;
            if ((adc_skip_mask >> adcboard) & 0x1U)
            {
              for (int iskip = 0; iskip < 64; iskip++)
              {
                std::vector<unsigned int> v_peak_sub_ped;
                for (int i = sample_start; i < sample_end; i++)
                {
                  v_peak_sub_ped.push_back(0);
                }
                unsigned int key = TowerInfoDefs::encode_emcal(iwave);
                m_peak_sub_ped_emcal[key] = v_peak_sub_ped;
                iwave++;
              }
            }
          }
          std::vector<unsigned int> v_peak_sub_ped;
          if (packet->iValue(channel, "SUPPRESSED"))
          {
            for (int i = sample_start; i < sample_end; i++)
            {
              v_peak_sub_ped.push_back(0);
            }
          }
          else
          {
            for (int i = sample_start; i < sample_end; i++)
            {
              int16_t maxim = (packet->iValue(i, channel) > packet->iValue(i + 1, channel) ? packet->iValue(i, channel) : packet->iValue(i + 1, channel));
              maxim = (maxim > packet->iValue(i + 2, channel) ? maxim : packet->iValue(i + 2, channel));
              uint16_t sam = 0;
              if (i >= m_trig_sub_delay)
              {
                sam = i - m_trig_sub_delay;
              }
              else
              {
                sam = 0;
              }
              unsigned int sub = 0;
              if (maxim > packet->iValue(sam, channel))
              {
                sub = (((uint16_t) (maxim - packet->iValue(sam, channel))) & 0x3fffU);
              }

              v_peak_sub_ped.push_back(sub);
            }
          }
          unsigned int key = TowerInfoDefs::encode_emcal(iwave);
          m_peak_sub_ped_emcal[key] = v_peak_sub_ped;
          iwave++;
        }
        if (nchannels < 192 && !(adc_skip_mask < 4))
        {
          for (int iskip = 0; iskip < 192 - nchannels; iskip++)
          {
            std::vector<unsigned int> v_peak_sub_ped;
            for (int i = sample_start; i < sample_end; i++)
            {
              v_peak_sub_ped.push_back(0);
            }
            unsigned int key = TowerInfoDefs::encode_emcal(iwave);
            m_peak_sub_ped_emcal[key] = v_peak_sub_ped;
            iwave++;
          }
        }
      }
    }
  }
  if (m_do_hcalout)
  {
    if (Verbosity())
    {
      std::cout << __FILE__ << "::" << __FUNCTION__ << ":: ohcal" << std::endl;
    }

    if (!m_hcal_packets)
    {
      return Fun4AllReturnCodes::ABORTRUN;
    }

    unsigned int iwave = 0;
    for (int pid = m_packet_low_hcalout; pid <= m_packet_high_hcalout; pid++)
    {
      CaloPacket *packet = m_hcal_packets->getPacketbyId(pid);
      if (packet)
      {
        int nchannels = packet->iValue(0, "CHANNELS");

        for (int channel = 0; channel < nchannels; channel++)
        {
          std::vector<unsigned int> v_peak_sub_ped;
          if (packet->iValue(channel, "SUPPRESSED"))
          {
            for (int i = sample_start; i < sample_end; i++)
            {
              v_peak_sub_ped.push_back(0);
            }
          }
          else
          {
            for (int i = sample_start; i < sample_end; i++)
            {
              int16_t maxim = (packet->iValue(i, channel) > packet->iValue(i + 1, channel) ? packet->iValue(i, channel) : packet->iValue(i + 1, channel));
              maxim = (maxim > packet->iValue(i + 2, channel) ? maxim : packet->iValue(i + 2, channel));
              uint16_t sam = 0;
              if (i >= m_trig_sub_delay)
              {
                sam = i - m_trig_sub_delay;
              }
              else
              {
                sam = 0;
              }
              unsigned int sub = 0;
              if (maxim > packet->iValue(sam, channel))
              {
                sub = (((uint16_t) (maxim - packet->iValue(sam, channel))) & 0x3fffU);
              }

              v_peak_sub_ped.push_back(sub);
            }
          }
          unsigned int key = TowerInfoDefs::encode_hcal(iwave);
          m_peak_sub_ped_hcalout[key] = v_peak_sub_ped;
          iwave++;
        }
      }
    }
  }
  if (m_do_hcalin)
  {
    if (Verbosity())
    {
      std::cout << __FILE__ << "::" << __FUNCTION__ << ":: ohcal" << std::endl;
    }
    if (!m_hcal_packets)
    {
      return Fun4AllReturnCodes::ABORTRUN;
    }
    unsigned int iwave = 0;
    for (int pid = m_packet_low_hcalin; pid <= m_packet_high_hcalin; pid++)
    {
      CaloPacket *packet = m_hcal_packets->getPacketbyId(pid);
      if (packet)
      {
        int nchannels = packet->iValue(0, "CHANNELS");

        for (int channel = 0; channel < nchannels; channel++)
        {
          std::vector<unsigned int> v_peak_sub_ped;
          if (packet->iValue(channel, "SUPPRESSED"))
          {
            for (int i = sample_start; i < sample_end; i++)
            {
              v_peak_sub_ped.push_back(0);
            }
          }
          else
          {
            for (int i = sample_start; i < sample_end; i++)
            {
              int16_t maxim = (packet->iValue(i, channel) > packet->iValue(i + 1, channel) ? packet->iValue(i, channel) : packet->iValue(i + 1, channel));
              maxim = (maxim > packet->iValue(i + 2, channel) ? maxim : packet->iValue(i + 2, channel));
              uint16_t sam = 0;
              if (i >= m_trig_sub_delay)
              {
                sam = i - m_trig_sub_delay;
              }
              else
              {
                sam = 0;
              }
              unsigned int sub = 0;
              if (maxim > packet->iValue(sam, channel))
              {
                sub = (((uint16_t) (maxim - packet->iValue(sam, channel))) & 0x3fffU);
              }

              v_peak_sub_ped.push_back(sub);
            }
          }
          unsigned int key = TowerInfoDefs::encode_hcal(iwave);
          m_peak_sub_ped_hcalin[key] = v_peak_sub_ped;
          iwave++;
        }
      }
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
int CaloTriggerEmulator::process_waveforms()
{
  if (!m_isdata)
  {
    return process_sim();
  }

  if (m_useoffline)
  {
    return process_offline();
  }

  if (m_event == nullptr)
  {
    std::cout << PHWHERE << " Event not found" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if (m_event->getEvtType() != DATAEVENT)
  {
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  // Get range of waveforms
  if (Verbosity())
  {
    std::cout << __FILE__ << "::" << __FUNCTION__ << ":: Processing waveforms" << std::endl;
  }

  int sample_start = 1;
  int sample_end = m_nsamples;
  if (m_trig_sample > 0)
  {
    sample_start = m_trig_sample;
    sample_end = m_trig_sample + 1;
  }

  if (m_do_emcal)
  {
    if (Verbosity())
    {
      std::cout << __FILE__ << "::" << __FUNCTION__ << ":: emcal" << std::endl;
    }
    unsigned int iwave = 0;
    for (int pid = m_packet_low_emcal; pid <= m_packet_high_emcal; pid++)
    {
      Packet *packet = m_event->getPacket(pid);
      if (packet)
      {
        int nchannels = packet->iValue(0, "CHANNELS");
        unsigned int adc_skip_mask = 0;

        adc_skip_mask = cdbttree_adcmask->GetIntValue(pid, m_fieldname);

        for (int channel = 0; channel < nchannels; channel++)
        {
          if (channel % 64 == 0)
          {
            unsigned int adcboard = (unsigned int) channel / 64;
            if ((adc_skip_mask >> adcboard) & 0x1U)
            {
              for (int iskip = 0; iskip < 64; iskip++)
              {
                std::vector<unsigned int> v_peak_sub_ped;
                for (int i = sample_start; i < sample_end; i++)
                {
                  v_peak_sub_ped.push_back(0);
                }
                unsigned int key = TowerInfoDefs::encode_emcal(iwave);
                m_peak_sub_ped_emcal[key] = v_peak_sub_ped;
                iwave++;
              }
              continue;
            }
          }
          std::vector<unsigned int> v_peak_sub_ped;
          if (packet->iValue(channel, "SUPPRESSED"))
          {
            for (int i = sample_start; i < sample_end; i++)
            {
              v_peak_sub_ped.push_back(0);
            }
          }
          else
          {
            for (int i = sample_start; i < sample_end; i++)
            {
              int16_t maxim = (packet->iValue(i, channel) > packet->iValue(i + 1, channel) ? packet->iValue(i, channel) : packet->iValue(i + 1, channel));
              maxim = (maxim > packet->iValue(i + 2, channel) ? maxim : packet->iValue(i + 2, channel));
              uint16_t sam = 0;
              if (i >= m_trig_sub_delay)
              {
                sam = i - m_trig_sub_delay;
              }
              else
              {
                sam = 0;
              }
              unsigned int sub = 0;
              if (maxim > packet->iValue(sam, channel))
              {
                sub = (((uint16_t) (maxim - packet->iValue(sam, channel))) & 0x3fffU);
              }

              v_peak_sub_ped.push_back(sub);
            }
          }
          unsigned int key = TowerInfoDefs::encode_emcal(iwave);
          m_peak_sub_ped_emcal[key] = v_peak_sub_ped;
          iwave++;
        }
      }
      delete packet;
    }
  }
  if (m_do_hcalout)
  {
    if (Verbosity())
    {
      std::cout << __FILE__ << "::" << __FUNCTION__ << ":: ohcal" << std::endl;
    }

    unsigned int iwave = 0;
    for (int pid = m_packet_low_hcalout; pid <= m_packet_high_hcalout; pid++)
    {
      Packet *packet = m_event->getPacket(pid);
      if (packet)
      {
        int nchannels = packet->iValue(0, "CHANNELS");

        for (int channel = 0; channel < nchannels; channel++)
        {
          std::vector<unsigned int> v_peak_sub_ped;
          if (packet->iValue(channel, "SUPPRESSED"))
          {
            for (int i = sample_start; i < sample_end; i++)
            {
              v_peak_sub_ped.push_back(0);
            }
          }
          else
          {
            for (int i = sample_start; i < sample_end; i++)
            {
              int16_t maxim = (packet->iValue(i, channel) > packet->iValue(i + 1, channel) ? packet->iValue(i, channel) : packet->iValue(i + 1, channel));
              maxim = (maxim > packet->iValue(i + 2, channel) ? maxim : packet->iValue(i + 2, channel));
              uint16_t sam = 0;
              if (i >= m_trig_sub_delay)
              {
                sam = i - m_trig_sub_delay;
              }
              else
              {
                sam = 0;
              }
              unsigned int sub = 0;
              if (maxim > packet->iValue(sam, channel))
              {
                sub = (((uint16_t) (maxim - packet->iValue(sam, channel))) & 0x3fffU);
              }

              v_peak_sub_ped.push_back(sub);
            }
          }
          unsigned int key = TowerInfoDefs::encode_hcal(iwave);
          m_peak_sub_ped_hcalout[key] = v_peak_sub_ped;
          iwave++;
        }
      }
      delete packet;
    }
  }
  if (m_do_hcalin)
  {
    if (Verbosity())
    {
      std::cout << __FILE__ << "::" << __FUNCTION__ << ":: ohcal" << std::endl;
    }

    unsigned int iwave = 0;
    for (int pid = m_packet_low_hcalin; pid <= m_packet_high_hcalin; pid++)
    {
      Packet *packet = m_event->getPacket(pid);
      if (packet)
      {
        int nchannels = packet->iValue(0, "CHANNELS");

        for (int channel = 0; channel < nchannels; channel++)
        {
          std::vector<unsigned int> v_peak_sub_ped;
          if (packet->iValue(channel, "SUPPRESSED"))
          {
            for (int i = sample_start; i < sample_end; i++)
            {
              v_peak_sub_ped.push_back(0);
            }
          }
          else
          {
            for (int i = sample_start; i < sample_end; i++)
            {
              int16_t maxim = (packet->iValue(i, channel) > packet->iValue(i + 1, channel) ? packet->iValue(i, channel) : packet->iValue(i + 1, channel));
              maxim = (maxim > packet->iValue(i + 2, channel) ? maxim : packet->iValue(i + 2, channel));
              uint16_t sam = 0;
              if (i >= m_trig_sub_delay)
              {
                sam = i - m_trig_sub_delay;
              }
              else
              {
                sam = 0;
              }
              unsigned int sub = 0;
              if (maxim > packet->iValue(sam, channel))
              {
                sub = (((uint16_t) (maxim - packet->iValue(sam, channel))) & 0x3fffU);
              }

              v_peak_sub_ped.push_back(sub);
            }
          }
          unsigned int key = TowerInfoDefs::encode_hcal(iwave);
          m_peak_sub_ped_hcalin[key] = v_peak_sub_ped;
          iwave++;
        }
      }
      delete packet;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int CaloTriggerEmulator::process_sim()
{
  // Get range of waveforms
  if (Verbosity())
  {
    std::cout << __FILE__ << "::" << __FUNCTION__ << ":: Processing waveforms" << std::endl;
  }

  int sample_start = 1;
  int sample_end = m_nsamples;
  if (m_trig_sample > 0)
  {
    sample_start = m_trig_sample;
    sample_end = m_trig_sample + 1;
  }

  if (m_do_emcal)
  {
    if (Verbosity())
    {
      std::cout << __FILE__ << "::" << __FUNCTION__ << ":: emcal" << std::endl;
    }
    if (!m_waveforms_emcal->size())
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }
    // for each waveform, clauclate the peak - pedestal given the sub-delay setting
    for (unsigned int iwave = 0; iwave < (unsigned int) m_waveforms_emcal->size(); iwave++)
    {
      std::vector<unsigned int> v_peak_sub_ped;
      TowerInfo *tower = m_waveforms_emcal->get_tower_at_channel(iwave);
      unsigned int key = TowerInfoDefs::encode_emcal(iwave);
      if (tower->get_isZS())
      {
        for (int i = sample_start; i < sample_end; i++)
        {
          v_peak_sub_ped.push_back(0);
        }
      }
      else
      {
        for (int i = sample_start; i < sample_end; i++)
        {
          int16_t maxim = (tower->get_waveform_value(i) > tower->get_waveform_value(i + 1) ? tower->get_waveform_value(i) : tower->get_waveform_value(i + 1));
          maxim = (maxim > tower->get_waveform_value(i + 2) ? maxim : tower->get_waveform_value(i + 2));
          uint16_t sam = 0;
          if (i >= m_trig_sub_delay)
          {
            sam = i - m_trig_sub_delay;
          }
          else
          {
            sam = 0;
          }
          unsigned int sub = 0;
          if (maxim > tower->get_waveform_value(sam))
          {
            sub = (((uint16_t) (maxim - tower->get_waveform_value(sam))) & 0x3fffU);
          }

          v_peak_sub_ped.push_back(sub);
        }
      }
      // save in global.
      m_peak_sub_ped_emcal[key] = v_peak_sub_ped;
    }
  }
  if (m_do_hcalout)
  {
    if (Verbosity())
    {
      std::cout << __FILE__ << "::" << __FUNCTION__ << ":: ohcal" << std::endl;
    }

    std::vector<int> wave;
    // for each waveform, clauclate the peak - pedestal given the sub-delay setting
    if (!m_waveforms_hcalout->size())
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }

    for (unsigned int iwave = 0; iwave < (unsigned int) m_waveforms_hcalout->size(); iwave++)
    {
      std::vector<unsigned int> v_peak_sub_ped;
      TowerInfo *tower = m_waveforms_hcalout->get_tower_at_channel(iwave);
      unsigned int key = TowerInfoDefs::encode_hcal(iwave);
      if (tower->get_isZS())
      {
        for (int i = sample_start; i < sample_end; i++)
        {
          v_peak_sub_ped.push_back(0);
        }
      }
      else
      {
        for (int i = sample_start; i < sample_end; i++)
        {
          int16_t maxim = (tower->get_waveform_value(i) > tower->get_waveform_value(i + 1) ? tower->get_waveform_value(i) : tower->get_waveform_value(i + 1));
          maxim = (maxim > tower->get_waveform_value(i + 2) ? maxim : tower->get_waveform_value(i + 2));
          uint16_t sam = 0;
          if (i >= m_trig_sub_delay)
          {
            sam = i - m_trig_sub_delay;
          }
          else
          {
            sam = 0;
          }
          unsigned int sub = 0;
          if (maxim > tower->get_waveform_value(sam))
          {
            sub = (((uint16_t) (maxim - tower->get_waveform_value(sam))) & 0x3fffU);
          }

          v_peak_sub_ped.push_back(sub);
        }
      }
      // save in global.
      m_peak_sub_ped_hcalout[key] = v_peak_sub_ped;
    }
  }
  if (m_do_hcalin)
  {
    if (Verbosity())
    {
      std::cout << __FILE__ << "::" << __FUNCTION__ << ":: ihcal" << std::endl;
    }
    if (!m_waveforms_hcalin->size())
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }
    if (Verbosity())
    {
      std::cout << __FILE__ << "::" << __FUNCTION__ << ":: ihcal" << std::endl;
    }

    std::vector<unsigned int> wave;

    // for each waveform, clauclate the peak - pedestal given the sub-delay setting
    for (unsigned int iwave = 0; iwave < (unsigned int) m_waveforms_hcalin->size(); iwave++)
    {
      std::vector<unsigned int> v_peak_sub_ped;
      TowerInfo *tower = m_waveforms_hcalin->get_tower_at_channel(iwave);
      unsigned int key = TowerInfoDefs::encode_hcal(iwave);
      if (tower->get_isZS())
      {
        for (int i = sample_start; i < sample_end; i++)
        {
          v_peak_sub_ped.push_back(0);
        }
      }
      else
      {
        for (int i = sample_start; i < sample_end; i++)
        {
          int16_t maxim = (tower->get_waveform_value(i) > tower->get_waveform_value(i + 1) ? tower->get_waveform_value(i) : tower->get_waveform_value(i + 1));
          maxim = (maxim > tower->get_waveform_value(i + 2) ? maxim : tower->get_waveform_value(i + 2));
          uint16_t sam = 0;
          if (i >= m_trig_sub_delay)
          {
            sam = i - m_trig_sub_delay;
          }
          else
          {
            sam = 0;
          }
          unsigned int sub = 0;
          if (maxim > tower->get_waveform_value(sam))
          {
            sub = (((uint16_t) (maxim - tower->get_waveform_value(sam))) & 0x3fffU);
          }

          v_peak_sub_ped.push_back(sub);
        }
      }
      // save in global.
      m_peak_sub_ped_hcalin[key] = v_peak_sub_ped;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

// procedure to process the peak - pedestal into primitives.
int CaloTriggerEmulator::process_primitives()
{
  int ip;
  int i;
  bool mask;
  int nsample = m_nsamples - 1;
  if (m_trig_sample > 0)
  {
    nsample = 1;
  }

  if (Verbosity())
  {
    std::cout << __FILE__ << "::" << __FUNCTION__ << ":: Processing primitives" << std::endl;
  }

  if (m_do_emcal)
  {
    if (Verbosity())
    {
      std::cout << __FILE__ << "::" << __FUNCTION__ << ":: Processing primitives:: emcal" << std::endl;
    }

    ip = 0;

    // get the number of primitives needed to process
    m_n_primitives = m_prim_map[TriggerDefs::DetectorId::emcalDId];
    for (i = 0; i < m_n_primitives; i++, ip++)
    {
      if (Verbosity())
      {
        std::cout << __FILE__ << "::" << __FUNCTION__ << ":: Processing primitives:: adding " << i << std::endl;
      }
      unsigned int tmp = 0;
      // get the primitive key of what we are making, in order of the packet ID and channel number
      TriggerDefs::TriggerPrimKey primkey = TriggerDefs::getTriggerPrimKey(TriggerDefs::GetTriggerId("NONE"), TriggerDefs::GetDetectorId("EMCAL"), TriggerDefs::GetPrimitiveId("EMCAL"), ip);

      TriggerPrimitive *primitive = m_primitives_emcal->get_primitive_at_key(primkey);
      unsigned int sum = 0;
      // check if masked Fiber;
      mask = CheckFiberMasks(primkey);

      // calculate 16 sums
      for (int isum = 0; isum < m_n_sums; isum++)
      {
        // get sum key
        TriggerDefs::TriggerSumKey sumkey = TriggerDefs::getTriggerSumKey(TriggerDefs::GetTriggerId("NONE"), TriggerDefs::GetDetectorId("EMCAL"), TriggerDefs::GetPrimitiveId("EMCAL"), ip, isum);

        // calculate sums for all samples, hense the vector.
        std::vector<unsigned int> *t_sum = primitive->get_sum_at_key(sumkey);
        t_sum->clear();

        // check to mask channel (if fiber masked, automatically mask the channel)
        bool mask_channel = mask || CheckChannelMasks(sumkey);
        for (int is = 0; is < nsample; is++)
        {
          sum = 0;
          unsigned int temp_sum = 0;

          // if masked, just fill with 0s
          if (!mask_channel)
          {
            for (int j = 0; j < 4; j++)
            {
              // unsigned int iwave = 64*ip + isum*4 + j;
              unsigned int key = TriggerDefs::GetTowerInfoKey(TriggerDefs::GetDetectorId("EMCAL"), ip, isum, j);
              unsigned int lut_input = (m_peak_sub_ped_emcal[key].at(is) >> 4U) & 0x3ffU;
              if (m_default_lut_emcal)
              {
                tmp = (m_l1_adc_table[lut_input] >> 2U);
              }
              else
              {
                unsigned int lut_output = ((unsigned int) h_emcal_lut[key]->GetBinContent(lut_input + 1)) & 0x3ffU;
                tmp = (lut_output >> 2U);
              }
              temp_sum += (tmp & 0xffU);
            }
            sum = ((temp_sum & 0x3ffU) >> 2U) & 0xffU;
            if (Verbosity() >= 10 && sum >= 1)
            {
              std::cout << __FILE__ << "::" << __FUNCTION__ << ":: emcal sum " << sumkey << " = " << sum << std::endl;
            }
          }
          if (Verbosity())
          {
            std::cout << __FILE__ << "::" << __FUNCTION__ << ":: Processing primitives:: adding" << std::endl;
          }

          t_sum->push_back(sum);
        }
      }
    }
  }
  if (m_do_hcalout)
  {
    if (Verbosity())
    {
      std::cout << __FILE__ << "::" << __FUNCTION__ << ":: Processing primitives:: ohcal" << std::endl;
    }

    ip = 0;

    m_n_primitives = m_prim_map[TriggerDefs::DetectorId::hcaloutDId];

    for (i = 0; i < m_n_primitives; i++, ip++)
    {
      TriggerDefs::TriggerPrimKey primkey = TriggerDefs::getTriggerPrimKey(TriggerDefs::GetTriggerId("NONE"), TriggerDefs::GetDetectorId("HCALOUT"), TriggerDefs::GetPrimitiveId("HCALOUT"), ip);
      TriggerPrimitive *primitive = m_primitives_hcalout->get_primitive_at_key(primkey);
      unsigned int sum;
      mask = CheckFiberMasks(primkey);
      for (int isum = 0; isum < m_n_sums; isum++)
      {
        TriggerDefs::TriggerSumKey sumkey = TriggerDefs::getTriggerSumKey(TriggerDefs::GetTriggerId("NONE"), TriggerDefs::GetDetectorId("HCALOUT"), TriggerDefs::GetPrimitiveId("HCALOUT"), ip, isum);
        std::vector<unsigned int> *t_sum = primitive->get_sum_at_key(sumkey);
        mask |= CheckChannelMasks(sumkey);
        for (int is = 0; is < nsample; is++)
        {
          sum = 0;
          unsigned int temp_sum = 0;
          if (!mask)
          {
            for (int j = 0; j < 4; j++)
            {
              unsigned int key = TriggerDefs::GetTowerInfoKey(TriggerDefs::GetDetectorId("HCAL"), ip, isum, j);
              unsigned int lut_input = (m_peak_sub_ped_hcalout[key].at(is) >> 4U) & 0x3ffU;
              unsigned int tmp = 0;
              if (m_default_lut_hcalout)
              {
                tmp = (m_l1_adc_table[lut_input] >> 2U);
              }
              else
              {
                unsigned int lut_output = ((unsigned int) h_hcalout_lut[key]->GetBinContent(lut_input + 1)) & 0x3ffU;
                tmp = (lut_output >> 2U);
              }
              temp_sum += (tmp & 0xffU);
            }
            sum = ((temp_sum & 0x3ffU) >> 2U) & 0xffU;
            if (Verbosity() >= 10 && sum >= 1)
            {
              std::cout << __FILE__ << "::" << __FUNCTION__ << ":: hcalout sum " << sumkey << " = " << sum << std::endl;
            }
          }
          t_sum->push_back(sum);
        }
      }
    }
  }
  if (m_do_hcalin)
  {
    ip = 0;

    if (Verbosity())
    {
      std::cout << __FILE__ << "::" << __FUNCTION__ << ":: Processing primitives:: ihcal" << std::endl;
    }

    m_n_primitives = m_prim_map[TriggerDefs::DetectorId::hcalinDId];

    for (i = 0; i < m_n_primitives; i++, ip++)
    {
      TriggerDefs::TriggerPrimKey primkey = TriggerDefs::getTriggerPrimKey(TriggerDefs::GetTriggerId("NONE"), TriggerDefs::GetDetectorId("HCALIN"), TriggerDefs::GetPrimitiveId("HCALIN"), ip);
      TriggerPrimitive *primitive = m_primitives_hcalin->get_primitive_at_key(primkey);
      unsigned int sum;
      mask = CheckFiberMasks(primkey);
      for (int isum = 0; isum < m_n_sums; isum++)
      {
        TriggerDefs::TriggerSumKey sumkey = TriggerDefs::getTriggerSumKey(TriggerDefs::GetTriggerId("NONE"), TriggerDefs::GetDetectorId("HCALIN"), TriggerDefs::GetPrimitiveId("HCALIN"), ip, isum);
        std::vector<unsigned int> *t_sum = primitive->get_sum_at_key(sumkey);
        mask |= CheckChannelMasks(sumkey);
        for (int is = 0; is < nsample; is++)
        {
          sum = 0;
          unsigned int temp_sum = 0;
          if (!mask)
          {
            for (int j = 0; j < 4; j++)
            {
              unsigned int key = TriggerDefs::GetTowerInfoKey(TriggerDefs::GetDetectorId("HCAL"), ip, isum, j);
              unsigned int lut_input = (m_peak_sub_ped_hcalin[key].at(is) >> 4U) & 0x3ffU;
              unsigned int tmp = 0;
              if (m_default_lut_hcalin)
              {
                tmp = (m_l1_adc_table[lut_input] >> 2U);
              }
              else
              {
                unsigned int lut_output = ((unsigned int) h_hcalin_lut[key]->GetBinContent(lut_input + 1)) & 0x3ffU;
                tmp = (lut_output >> 2U);
              }
              temp_sum += (tmp & 0x3ffU);
            }
            sum = ((temp_sum & 0xfffU) >> 2U) & 0xffU;
            if (Verbosity() >= 10 && sum >= 1)
            {
              std::cout << __FILE__ << "::" << __FUNCTION__ << ":: hcalin sum " << sumkey << " = " << sum << std::endl;
            }
          }
          t_sum->push_back(sum);
        }
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

// Unless this is the MBD or HCAL Cosmics trigger, EMCAL and HCAL will go through here.
// This creates the 8x8 non-overlapping sum and the 4x4 overlapping sum.

int CaloTriggerEmulator::process_organizer()
{
  TriggerPrimitiveContainer::Range range;
  if (Verbosity())
  {
    std::cout << __FILE__ << "::" << __FUNCTION__ << ":: Processing organizer" << std::endl;
  }

  int nsample = m_nsamples - 1;
  // bits are to say whether the trigger has fired. this is what is sent to the GL1
  if (m_trig_sample > 0)
  {
    nsample = 1;
  }
  // 8x8 non-overlapping sums in the EMCAL
  // create the 8x8 non-overlapping sum
  {
    if (Verbosity() >= 2)
    {
      std::cout << __FUNCTION__ << " " << __LINE__ << " processing 8x8 non-overlapping sums" << std::endl;
    }

    m_triggerid = TriggerDefs::TriggerId::jetTId;

    if (!m_primitives_emcal)
    {
      std::cout << "There is no primitive container" << std::endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }

    range = m_primitives_emcal_ll1->getTriggerPrimitives();
    for (TriggerPrimitiveContainerv1::Iter iter = range.first; iter != range.second; ++iter)
    {
      TriggerPrimitivev1::Range sumrange = iter->second->getSums();
      for (TriggerPrimitivev1::Iter siter = sumrange.first; siter != sumrange.second; ++siter)
      {
        for (int is = 0; is < nsample; is++)
        {
          siter->second->push_back(0);
        }
      }
    }

    if (Verbosity())
    {
      std::cout << __FILE__ << "::" << __FUNCTION__ << ":: Processing getting emcal prims" << std::endl;
    }

    // iterate through emcal primitives and organize into the 16 jet primitives each with the 8x8 nonoverlapping sum
    range = m_primitives_emcal->getTriggerPrimitives();

    for (TriggerPrimitiveContainerv1::Iter iter = range.first; iter != range.second; ++iter)
    {
      // get key and see if masked
      TriggerDefs::TriggerPrimKey key = (*iter).first;
      if (CheckFiberMasks(key))
      {
        continue;
      }

      uint16_t sumphi = TriggerDefs::getPrimitivePhiId_from_TriggerPrimKey(key);
      uint16_t sumeta = TriggerDefs::getPrimitiveEtaId_from_TriggerPrimKey(key);

      // based on where the primitive is in the detector, the location of the jet primitive is determined, 0 through 15 in phi.
      uint16_t iprim = sumphi / 2;
      // eta determines the location of the sum within the jet primitive.
      uint16_t isum = (sumeta + (sumphi % 2) * 12);

      TriggerDefs::TriggerPrimKey jet_prim_key = TriggerDefs::getTriggerPrimKey(m_triggerid, TriggerDefs::GetDetectorId("EMCAL"), TriggerDefs::GetPrimitiveId("JET"), iprim);

      TriggerDefs::TriggerPrimKey jet_sum_key = TriggerDefs::getTriggerSumKey(m_triggerid, TriggerDefs::GetDetectorId("EMCAL"), TriggerDefs::GetPrimitiveId("JET"), iprim, isum);

      // add to primitive previously made the sum of the 8x8 non-overlapping sum.
      std::vector<unsigned int> *t_sum = m_primitives_emcal_ll1->get_primitive_at_key(jet_prim_key)->get_sum_at_key(jet_sum_key);

      // get the primitive (16 2x2 sums)
      TriggerPrimitive *primitive = (*iter).second;
      TriggerPrimitivev1::Range sumrange = primitive->getSums();

      // iterate through all 16 sums and add together
      for (TriggerPrimitivev1::Iter iter_sum = sumrange.first; iter_sum != sumrange.second; ++iter_sum)
      {
        if (CheckChannelMasks(iter_sum->first))
        {
          continue;
        }
        int i = 0;
        for (unsigned int &it_s : *(*iter_sum).second)
        {
          t_sum->at(i) += (it_s & 0xffU);
          i++;
        }
      }

      // bit shift by 16 (divide by the 16 towers) to get an 8 bit energy sum.
      for (unsigned int &it_s : *t_sum)
      {
        // unsigned int sumshift = ((it_s >> 0x2U) & 0x3ffU);
        // unsigned int sum_lower = ( sumshift & 0x7fU );
        // unsigned int sum_higher = ( ( sumshift >> 0x7U ) > 0 ? 0x1U : 0x0U );
        // it_s = m_l1_8x8_table[sumshift];
        if (it_s > 0xffU)
        {
          it_s = 0xffU;
        }
      }
    }
  }
  {
    if (Verbosity())
    {
      std::cout << __FILE__ << "::" << __FUNCTION__ << ":: Processing HCAL" << std::endl;
    }

    range = m_primitives_hcal_ll1->getTriggerPrimitives();
    for (TriggerPrimitiveContainerv1::Iter iter = range.first; iter != range.second; ++iter)
    {
      TriggerPrimitivev1::Range sumrange = iter->second->getSums();
      for (TriggerPrimitivev1::Iter siter = sumrange.first; siter != sumrange.second; ++siter)
      {
        for (int is = 0; is < nsample; is++)
        {
          siter->second->push_back(0);
        }
      }
    }

    if (m_primitives_hcalin)
    {
      if (Verbosity())
      {
        std::cout << __FILE__ << "::" << __FUNCTION__ << ":: Processing organizer" << std::endl;
      }

      // iterate through emcal primitives and organize into the 16 jet primitives each with the 8x8 nonoverlapping sum
      range = m_primitives_hcalin->getTriggerPrimitives();

      // for hcalin: there are
      for (TriggerPrimitiveContainerv1::Iter iter = range.first; iter != range.second; ++iter)
      {
        // get key and see if masked
        TriggerDefs::TriggerPrimKey key = (*iter).first;
        if (CheckFiberMasks(key))
        {
          continue;
        }

        // get the primitive (16 2x2 sums)
        TriggerPrimitive *primitive = (*iter).second;
        TriggerPrimitivev1::Range sumrange = primitive->getSums();

        for (TriggerPrimitivev1::Iter iter_sum = sumrange.first; iter_sum != sumrange.second; ++iter_sum)
        {
          TriggerDefs::TriggerSumKey sumkey = (*iter_sum).first;
          uint16_t sumphi = TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(sumkey) * 4 + TriggerDefs::getSumPhiId(sumkey);
          uint16_t sumeta = TriggerDefs::getPrimitiveEtaId_from_TriggerSumKey(sumkey) * 4 + TriggerDefs::getSumEtaId(sumkey);

          int i = 0;
          if (CheckChannelMasks(sumkey))
          {
            continue;
          }
          // based on where the primitive is in the detector, the location of the jet primitive is determined, 0 through 15 in phi.
          uint16_t iprim = sumphi / 2;
          // eta determines the location of the sum within the jet primitive.
          uint16_t isum = (sumeta + (sumphi % 2) * 12);

          TriggerDefs::TriggerPrimKey jet_prim_key = TriggerDefs::getTriggerPrimKey(m_triggerid, TriggerDefs::GetDetectorId("HCAL"), TriggerDefs::GetPrimitiveId("JET"), iprim);

          TriggerDefs::TriggerPrimKey jet_sum_key = TriggerDefs::getTriggerSumKey(m_triggerid, TriggerDefs::GetDetectorId("HCAL"), TriggerDefs::GetPrimitiveId("JET"), iprim, isum);

          // add to primitive previously made the sum of the 8x8 non-overlapping sum.
          std::vector<unsigned int> *t_sum = m_primitives_hcal_ll1->get_primitive_at_key(jet_prim_key)->get_sum_at_key(jet_sum_key);

          for (unsigned int &it_s : *(*iter_sum).second)
          {
            t_sum->at(i) += (it_s & 0xffU);
            i++;
          }
        }
      }
    }

    if (m_primitives_hcalout)
    {
      if (Verbosity())
      {
        std::cout << __FILE__ << "::" << __FUNCTION__ << ":: Processing organizer" << std::endl;
      }

      // iterate through emcal primitives and organize into the 16 jet primitives each with the 8x8 nonoverlapping sum
      range = m_primitives_hcalout->getTriggerPrimitives();

      // for hcalin: there are
      for (TriggerPrimitiveContainerv1::Iter iter = range.first; iter != range.second; ++iter)
      {
        // get key and see if masked
        TriggerDefs::TriggerPrimKey key = (*iter).first;
        if (CheckFiberMasks(key))
        {
          continue;
        }

        // get the primitive (16 2x2 sums)
        TriggerPrimitive *primitive = (*iter).second;
        TriggerPrimitivev1::Range sumrange = primitive->getSums();

        for (TriggerPrimitivev1::Iter iter_sum = sumrange.first; iter_sum != sumrange.second; ++iter_sum)
        {
          TriggerDefs::TriggerSumKey sumkey = (*iter_sum).first;
          uint16_t sumphi = TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(sumkey) * 4 + TriggerDefs::getSumPhiId(sumkey);
          //		      uint16_t sumeta = TriggerDefs::getPrimitiveEtaId_from_TriggerSumKey(sumkey)*4 + TriggerDefs::getSumEtaId(sumkey);
          uint16_t sumeta = TriggerDefs::getPrimitiveEtaId_from_TriggerSumKey(sumkey) * 4 + TriggerDefs::getSumEtaId(sumkey);

          int i = 0;

          if (CheckChannelMasks(sumkey))
          {
            continue;
          }
          // based on where the primitive is in the detector, the location of the jet primitive is determined, 0 through 15 in phi.
          uint16_t iprim = sumphi / 2;
          // eta determines the location of the sum within the jet primitive.
          uint16_t isum = sumeta + (sumphi % 2) * 12;
          TriggerDefs::TriggerPrimKey jet_prim_key = TriggerDefs::getTriggerPrimKey(TriggerDefs::TriggerId::jetTId, TriggerDefs::GetDetectorId("HCAL"), TriggerDefs::GetPrimitiveId("JET"), iprim);

          TriggerDefs::TriggerPrimKey jet_sum_key = TriggerDefs::getTriggerSumKey(TriggerDefs::TriggerId::jetTId, TriggerDefs::GetDetectorId("HCAL"), TriggerDefs::GetPrimitiveId("JET"), iprim, isum);
          if (Verbosity())
          {
            std::cout << __FILE__ << "::" << __FUNCTION__ << "::" << __LINE__ << ":: Processing organizer" << std::endl;
          }

          std::vector<unsigned int> *t_sum = m_primitives_hcal_ll1->get_primitive_at_key(jet_prim_key)->get_sum_at_key(jet_sum_key);
          for (unsigned int &it_s : *(*iter_sum).second)
          {
            t_sum->at(i) += ((it_s) &0xffU);
            i++;
          }
        }
      }
    }

    // iterate through emcal primitives and organize into the 16 jet primitives each with the 8x8 nonoverlapping sum
    if (Verbosity())
    {
      std::cout << __FILE__ << "::" << __FUNCTION__ << "::" << __LINE__ << ":: Processing organizer" << std::endl;
    }
    range = m_primitives_hcal_ll1->getTriggerPrimitives();
    for (TriggerPrimitiveContainerv1::Iter iter = range.first; iter != range.second; ++iter)
    {
      TriggerPrimitive *primitive = (*iter).second;
      TriggerPrimitivev1::Range sumrange = primitive->getSums();

      for (TriggerPrimitivev1::Iter iter_sum = sumrange.first; iter_sum != sumrange.second; ++iter_sum)
      {
        for (unsigned int &it_s : *(*iter_sum).second)
        {
          it_s = (it_s >> 1U) & 0xffU;
        }
      }
    }

    {
      range = m_primitives_jet->getTriggerPrimitives();
      for (TriggerPrimitiveContainerv1::Iter iter = range.first; iter != range.second; ++iter)
      {
        TriggerPrimitivev1::Range sumrange = iter->second->getSums();
        for (TriggerPrimitivev1::Iter siter = sumrange.first; siter != sumrange.second; ++siter)
        {
          for (int is = 0; is < nsample; is++)
          {
            siter->second->push_back(0);
          }
        }
      }
    }
    if (Verbosity())
    {
      std::cout << __FILE__ << "::" << __FUNCTION__ << "::" << __LINE__ << ":: Processing organizer" << std::endl;
    }
    // get jet primitives (after EMCAL and HCAL sum)
    range = m_primitives_jet->getTriggerPrimitives();
    for (TriggerPrimitiveContainerv1::Iter iter = range.first; iter != range.second; ++iter)
    {
      if (Verbosity())
      {
        std::cout << __FILE__ << "::" << __FUNCTION__ << "::" << __LINE__ << ":: Processing organizer" << std::endl;
      }

      TriggerDefs::TriggerPrimKey jet_pkey = (*iter).first;
      TriggerDefs::TriggerPrimKey hcal_pkey = TriggerDefs::getTriggerPrimKey(TriggerDefs::TriggerId::jetTId, TriggerDefs::GetDetectorId("HCAL"), TriggerDefs::GetPrimitiveId("JET"), TriggerDefs::getPrimitiveLocId_from_TriggerPrimKey(jet_pkey));
      TriggerDefs::TriggerPrimKey emcal_pkey = TriggerDefs::getTriggerPrimKey(TriggerDefs::TriggerId::jetTId, TriggerDefs::GetDetectorId("EMCAL"), TriggerDefs::GetPrimitiveId("JET"), TriggerDefs::getPrimitiveLocId_from_TriggerPrimKey(jet_pkey));
      TriggerPrimitive *primitive = (*iter).second;
      TriggerPrimitivev1::Range sumrange = primitive->getSums();
      for (TriggerPrimitivev1::Iter iter_sum = sumrange.first; iter_sum != sumrange.second; ++iter_sum)
      {
        if (Verbosity())
        {
          std::cout << __FILE__ << "::" << __FUNCTION__ << "::" << __LINE__ << ":: Processing organizer" << std::endl;
        }

        TriggerDefs::TriggerSumKey jet_skey = (*iter_sum).first;

        TriggerDefs::TriggerSumKey hcal_skey = TriggerDefs::getTriggerSumKey(TriggerDefs::TriggerId::jetTId, TriggerDefs::GetDetectorId("HCAL"), TriggerDefs::GetPrimitiveId("JET"), TriggerDefs::getPrimitiveLocId_from_TriggerPrimKey(jet_pkey), TriggerDefs::getSumLocId(jet_skey));
        TriggerDefs::TriggerSumKey emcal_skey = TriggerDefs::getTriggerSumKey(TriggerDefs::TriggerId::jetTId, TriggerDefs::GetDetectorId("EMCAL"), TriggerDefs::GetPrimitiveId("JET"), TriggerDefs::getPrimitiveLocId_from_TriggerPrimKey(jet_pkey), TriggerDefs::getSumLocId(jet_skey));

        int i = 0;
        if (Verbosity())
        {
          std::cout << __FILE__ << "::" << __FUNCTION__ << "::" << __LINE__ << ":: Processing organizer" << std::endl;
        }

        for (unsigned int &it_s : *(iter_sum->second))
        {
          unsigned int sum_hcal = m_primitives_hcal_ll1->get_primitive_at_key(hcal_pkey)->get_sum_at_key(hcal_skey)->at(i);
          unsigned int sum_emcal = m_primitives_emcal_ll1->get_primitive_at_key(emcal_pkey)->get_sum_at_key(emcal_skey)->at(i);
          it_s = ((sum_hcal >> 1U) + (sum_emcal >> 1U)) & 0xffU;
          i++;
        }
      }
    }
  }

  return 0;
}

// This is where the LL1 Jet/Pair/Cosmic algorithm is

int CaloTriggerEmulator::process_trigger()
{
  std::vector<unsigned int> bits;

  int nsample = m_nsamples - 1;
  // bits are to say whether the trigger has fired. this is what is sent to the GL1
  if (m_trig_sample > 0)
  {
    nsample = 1;
  }
  bits.reserve(nsample);
  for (int is = 0; is < nsample; is++)
  {
    bits.push_back(0);
  }

  // photon
  // 8x8 non-overlapping sums in the EMCAL
  // create the 8x8 non-overlapping sum

  {
    m_triggerid = TriggerDefs::TriggerId::photonTId;
    std::vector<unsigned int> *trig_bits = m_ll1out_photon->GetTriggerBits();
    if (Verbosity() >= 2)
    {
      std::cout << __FUNCTION__ << " " << __LINE__ << " processing PHOTON trigger , bits before: " << trig_bits->size() << std::endl;
    }

    TriggerPrimitiveContainer::Range range = m_primitives_emcal_ll1->getTriggerPrimitives();

    for (TriggerPrimitiveContainerv1::Iter iter = range.first; iter != range.second; ++iter)
    {
      // get key and see if masked
      TriggerDefs::TriggerPrimKey key = (*iter).first;
      if (CheckFiberMasks(key))
      {
        continue;
      }

      // get the primitive (16 2x2 sums)
      TriggerPrimitive *primitive = (*iter).second;

      TriggerPrimitivev1::Range sumrange = primitive->getSums();

      // iterate through all 24 sums and add together
      for (TriggerPrimitivev1::Iter iter_sum = sumrange.first; iter_sum != sumrange.second; ++iter_sum)
      {
        // check if sum is greater than threshold.
        if (Verbosity() >= 2)
        {
          std::cout << __FUNCTION__ << " " << __LINE__ << " processing PHOTON trigger" << std::endl;
        }

        std::vector<unsigned int> *t_sum = (*iter_sum).second;
        TriggerDefs::TriggerSumKey sumk = (*iter_sum).first;
        for (int is = 0; is < nsample; is++)
        {
          unsigned short bit = getBits(t_sum->at(is), TriggerDefs::TriggerId::photonTId);
          if (bit)
          {
            m_ll1out_photon->addTriggeredSum(sumk, t_sum->at(is));
            m_ll1out_photon->addTriggeredPrimitive(key);
          }
          bits.at(is) |= bit;
        }
      }
    }

    uint16_t pass = 0;
    for (int is = 0; is < nsample; is++)
    {
      pass |= bits.at(is);
      trig_bits->push_back(bits.at(is));
    }

    if (pass)
    {
      m_photon_npassed++;
    }
  }
  {
    if (Verbosity() >= 2)
    {
      std::cout << __FUNCTION__ << " " << __LINE__ << " processing JET trigger" << std::endl;
    }

    // Make the jet primitives
    m_triggerid = TriggerDefs::TriggerId::jetTId;
    std::vector<unsigned int> *trig_bits = m_ll1out_jet->GetTriggerBits();
    std::vector<unsigned int> jet_map[32][9]{};
    for (auto &ie : jet_map)
    {
      for (auto &ip : ie)
      {
        for (int is = 0; is < nsample; is++)
        {
          ip.push_back(0);
        }
      }
    }

    if (!m_primitives_jet)
    {
      std::cout << "There is no primitive container" << std::endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }

    // iterate through emcal primitives and organize into the 16 jet primitives each with the 8x8 nonoverlapping sum
    TriggerPrimitiveContainer::Range range = m_primitives_jet->getTriggerPrimitives();

    for (TriggerPrimitiveContainerv1::Iter iter = range.first; iter != range.second; ++iter)
    {
      // get the primitive (16 2x2 sums)
      TriggerPrimitive *primitive = (*iter).second;
      TriggerPrimitivev1::Range sumrange = primitive->getSums();

      // iterate through all 16 sums and add together
      for (TriggerPrimitivev1::Iter iter_sum = sumrange.first; iter_sum != sumrange.second; ++iter_sum)
      {
        TriggerDefs::TriggerSumKey sumkey = (*iter_sum).first;

        int i = 0;
        int sum_phi = static_cast<int>(TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(sumkey) * 2 + TriggerDefs::getSumPhiId(sumkey));
        int sum_eta = static_cast<int>(TriggerDefs::getSumEtaId(sumkey));
        if (Verbosity() >= 2)
        {
          std::cout << __FUNCTION__ << " " << __LINE__ << " processing JET trigger " << sum_phi << " " << sum_eta << std::endl;
        }

        for (unsigned int &it_s : *(iter_sum->second))
        {
          for (int ijeta = (sum_eta <= 3 ? 0 : sum_eta - 3); ijeta <= (sum_eta > 8 ? 8 : sum_eta); ijeta++)
          {
            for (int ijphi = sum_phi - 3; ijphi <= sum_phi; ijphi++)
            {
              int iphi = (ijphi < 0 ? 32 + ijphi : ijphi);
              jet_map[iphi][ijeta].at(i) += it_s;
            }
          }
          i++;
        }
      }
    }
    if (Verbosity() >= 2)
    {
      std::cout << __FUNCTION__ << " " << __LINE__ << " processing JET trigger" << std::endl;
    }

    int pass = 0;
    for (int ijphi = 0; ijphi < 32; ijphi++)
    {
      for (int ijeta = 0; ijeta < 9; ijeta++)
      {
        if (Verbosity() >= 2)
        {
          std::cout << __FUNCTION__ << " " << __LINE__ << " processing JET trigger " << ijphi << " " << ijeta << std::endl;
        }

        unsigned int sk = ((unsigned int) ijphi & 0xffffU) + (((unsigned int) ijeta & 0xffffU) << 16U);
        std::vector<unsigned int> *sum = m_ll1out_jet->get_word(sk);
        if (Verbosity() >= 2)
        {
          std::cout << __FUNCTION__ << " " << __LINE__ << " processing JET trigger " << ijphi << " " << ijeta << std::endl;
        }

        for (int is = 0; is < nsample; is++)
        {
          if (Verbosity() >= 2)
          {
            std::cout << __FUNCTION__ << " " << __LINE__ << " processing JET trigger " << ijphi << " " << ijeta << std::endl;
          }

          sum->push_back(jet_map[ijphi][ijeta].at(is));
          unsigned short bit = getBits(jet_map[ijphi][ijeta].at(is), TriggerDefs::TriggerId::jetTId);

          if (bit)
          {
            m_ll1out_jet->addTriggeredSum(sk, jet_map[ijphi][ijeta].at(is));
            m_ll1out_jet->addTriggeredPrimitive(sk);
            pass = 1;
          }
          bits.at(is) |= bit;
        }
      }
    }
    if (Verbosity() >= 2)
    {
      std::cout << __FUNCTION__ << " " << __LINE__ << " processing JET trigger" << std::endl;
    }

    for (int is = 0; is < nsample; is++)
    {
      trig_bits->push_back(bits.at(is));
    }

    if (pass)
    {
      m_jet_npassed++;
    }
  }
  // //pair trigger here
  // {
  //   if (Verbosity())
  //     {
  // 	std::cout << __FILE__ << "::" << __FUNCTION__ << "::" << std::dec <<__LINE__ << ":: Processing organizer" << std::endl;
  //     }

  //   m_triggerid = TriggerDefs::TriggerId::pairTId;
  //   if (m_primitives_emcal)
  //     {
  // 	// iterate through emcal primitives and organize into the 16 jet primitives each with the 8x8 nonoverlapping sum
  // 	TriggerPrimitiveContainerv1::Range range = m_primitives_emcal->getTriggerPrimitives();
  // 	std::vector<unsigned int> pair_map[128][47]{};
  // 	for (auto & ie : pair_map)
  // 	  {
  // 	    for (auto & ip : ie)
  // 	      {
  // 		for (int is = 0; is < nsample; is++)
  // 		  {
  // 		    ip.push_back(0);
  // 		  }
  // 	      }
  // 	  }

  // 	for (TriggerPrimitiveContainerv1::Iter iter = range.first; iter != range.second; ++iter)
  // 	  {
  // 	    if (Verbosity())
  // 	      {
  // 		std::cout << __FILE__ << "::" << __FUNCTION__ << "::" << std::dec <<__LINE__ << ":: Processing organizer" << std::endl;
  // 	      }

  // 	    // get key and see if masked
  // 	    TriggerDefs::TriggerPrimKey key = (*iter).first;
  // 	    if (CheckFiberMasks(key))
  // 	      {
  // 		if (Verbosity() >= 2)
  // 		  {
  // 		    std::cout << "masked: " << key << std::endl;
  // 		  }
  // 		continue;
  // 	      }

  // 	    //get eta and phi from emcal
  // 	    uint16_t primphi = TriggerDefs::getPrimitivePhiId_from_TriggerPrimKey(key);
  // 	    uint16_t primeta = TriggerDefs::getPrimitiveEtaId_from_TriggerPrimKey(key);

  // 	    TriggerPrimitivev1::Range sumrange = (iter->second)->getSums();
  // 	    for (TriggerPrimitivev1::Iter iter_sum = sumrange.first; iter_sum != sumrange.second; ++iter_sum)
  // 	      {

  // 		if (CheckChannelMasks(iter_sum->first))
  // 		  {
  // 		    continue;
  // 		  }
  // 		uint16_t sumphi = TriggerDefs::getSumPhiId(iter_sum->first);
  // 		uint16_t sumeta = TriggerDefs::getSumEtaId(iter_sum->first);

  // 		uint16_t globalphi = primphi*4 + sumphi;
  // 		uint16_t globaleta = primeta*4 + sumeta;

  // 		uint16_t lowerphi = 127;
  // 		if (globalphi > 0)
  // 		  {
  // 		    lowerphi = globalphi - 1;
  // 		  }

  // 		int i = 0;
  // 		for (unsigned int &it_s : *(iter_sum->second))
  // 		  {

  // 		    if (globaleta > 0)
  // 		      {
  // 			pair_map[globalphi][globaleta-1].at(i) += it_s;
  // 			pair_map[lowerphi][globaleta-1].at(i) += it_s;
  // 		      }
  // 		    if (globaleta < 47)
  // 		      {
  // 			pair_map[globalphi][globaleta].at(i) += it_s;
  // 			pair_map[lowerphi][globaleta].at(i) += it_s;
  // 		      }
  // 		    i++;
  // 		  }
  // 	      }
  // 	  }

  // 	// no go through and find the peaks.
  // 	//std::vector<unsigned int> pair_map[128][47]{};
  // 	unsigned int maximum = 0;
  // 	int philoc = -1;
  // 	int etaloc = -1;
  // 	for (int ieta = 1; ieta < 47; ieta++)
  // 	  {
  // 	    for (int iphi = 0; iphi < 128; iphi++)
  // 	      {
  // 		if (Verbosity())
  // 		  {
  // 		    std::cout << __FILE__ << "::" << __FUNCTION__ << "::" << std::dec <<__LINE__ << ":: Processing pair "<< iphi << " " << ieta << std::endl;
  // 		  }

  // 		int lowerphi = iphi - 1;
  // 		int higherphi = iphi + 1;
  // 		if (Verbosity())
  // 		  {
  // 		    std::cout << __FILE__ << "::" << __FUNCTION__ << "::" << std::dec <<__LINE__ << ":: Processing pair look"<< lowerphi << " " << higherphi << std::endl;
  // 		  }

  // 		if (higherphi > 127) higherphi = 0;
  // 		if (lowerphi < 0) lowerphi = 127;

  // 		if (Verbosity())
  // 		  {
  // 		    std::cout << __FILE__ << "::" << __FUNCTION__ << "::" << std::dec <<__LINE__ << ":: Processing pair look"<< lowerphi << " " << higherphi << std::endl;
  // 		  }

  // 		if (pair_map[iphi][ieta].at(0) <= maximum) {continue;}

  // 		if (pair_map[iphi][ieta].at(0) < pair_map[iphi][ieta + 1].at(0)) {continue;}
  // 		if (pair_map[iphi][ieta].at(0) < pair_map[iphi + 1][ieta + 1].at(0)) {continue;}
  // 		if (pair_map[iphi][ieta].at(0) < pair_map[iphi - 1][ieta + 1].at(0)) {continue;}
  // 		if (pair_map[iphi][ieta].at(0) < pair_map[iphi + 1][ieta].at(0)) {continue;}
  // 		if (pair_map[iphi][ieta].at(0) <= pair_map[iphi - 1][ieta].at(0)) {continue;}
  // 		if (pair_map[iphi][ieta].at(0) <= pair_map[iphi + 1][ieta - 1].at(0)) {continue;}
  // 		if (pair_map[iphi][ieta].at(0) <= pair_map[iphi][ieta - 1].at(0)) {continue;}
  // 		if (pair_map[iphi][ieta].at(0) <= pair_map[iphi - 1][ieta - 1].at(0)) {continue;}
  // 		maximum = pair_map[iphi][ieta].at(0);
  // 		philoc = iphi;
  // 		etaloc = ieta;
  // 		if (Verbosity())
  // 		  {
  // 		    std::cout << __FILE__ << "::" << __FUNCTION__ << "::" << std::dec <<__LINE__ << ":: Processing pair look"<< maximum << std::endl;
  // 		  }

  // 	      }
  // 	  }

  // 	std::vector<unsigned int> *trig_bits = m_ll1out_pair->GetTriggerBits();
  // 	unsigned int bit = getBits(maximum, TriggerDefs::TriggerId::pairTId);

  // 	if (bit)
  // 	  {
  // 	    unsigned int sumloc = ((unsigned int) philoc) & 0xffffU;
  // 	    sumloc += (((unsigned int) etaloc)  & 0xffffU ) << 16U;
  // 	    unsigned int primloc = philoc/4;
  // 	    m_ll1out_pair->addTriggeredSum(sumloc);
  // 	    m_ll1out_pair->addTriggeredPrimitive(primloc);
  // 	    m_pair_npassed++;
  // 	  }
  // 	trig_bits->push_back(bit);
  //     }

  return Fun4AllReturnCodes::EVENT_OK;
}

void CaloTriggerEmulator::GetNodes(PHCompositeNode *topNode)
{
  if (Verbosity() >= 2)
  {
    std::cout << __FUNCTION__ << std::endl;
  }

  m_event = findNode::getClass<Event>(topNode, "PRDF");

  m_ll1_nodename = "LL1OUT_PAIR";

  m_ll1out_pair = findNode::getClass<LL1Out>(topNode, m_ll1_nodename);

  if (!m_ll1out_pair)
  {
    std::cout << "No LL1Out found... " << std::endl;
    exit(1);
  }

  m_prim_nodename = "TRIGGERPRIMITIVES_PAIR";
  m_primitives_pair = findNode::getClass<TriggerPrimitiveContainer>(topNode, m_prim_nodename);

  if (!m_primitives_pair)
  {
    std::cout << "No TriggerPrimitives found... " << std::endl;
    exit(1);
  }

  m_ll1_nodename = "LL1OUT_JET";
  m_ll1out_jet = findNode::getClass<LL1Out>(topNode, m_ll1_nodename);

  if (!m_ll1out_jet)
  {
    std::cout << "No LL1Out found... " << std::endl;
    exit(1);
  }

  m_prim_nodename = "TRIGGERPRIMITIVES_JET";
  m_primitives_jet = findNode::getClass<TriggerPrimitiveContainer>(topNode, m_prim_nodename);

  if (!m_primitives_jet)
  {
    std::cout << "No TriggerPrimitives found... " << std::endl;
    exit(1);
  }

  m_ll1_nodename = "LL1OUT_PHOTON";
  m_ll1out_photon = findNode::getClass<LL1Out>(topNode, m_ll1_nodename);

  if (!m_ll1out_photon)
  {
    std::cout << "No LL1Out found... " << std::endl;
    exit(1);
  }

  m_prim_nodename = "TRIGGERPRIMITIVES_PHOTON";
  m_primitives_photon = findNode::getClass<TriggerPrimitiveContainer>(topNode, m_prim_nodename);

  if (!m_primitives_photon)
  {
    std::cout << "No TriggerPrimitives found... " << std::endl;
    exit(1);
  }

  bool hcalset = false;
  if (m_do_hcalout)
  {
    if (!m_isdata)
    {
      m_waveforms_hcalout = findNode::getClass<TowerInfoContainer>(topNode, "WAVEFORM_HCALOUT");

      if (!m_waveforms_hcalout)
      {
        std::cout << "No HCALOUT Waveforms found... " << std::endl;
        exit(1);
      }
    }
    m_hcal_packets = findNode::getClass<CaloPacketContainer>(topNode, "HCALPackets");

    if (m_hcal_packets)
    {
      m_useoffline = true;
    }

    m_primitives_hcalout = findNode::getClass<TriggerPrimitiveContainer>(topNode, "TRIGGERPRIMITIVES_HCALOUT");

    if (!m_primitives_hcalout)
    {
      std::cout << "No HCALOUT Primitives found... " << std::endl;
      exit(1);
    }

    m_primitives_hcal_ll1 = findNode::getClass<TriggerPrimitiveContainer>(topNode, "TRIGGERPRIMITIVES_HCAL_LL1");

    if (!m_primitives_hcal_ll1)
    {
      std::cout << "No HCAL Primitives found... " << std::endl;
      exit(1);
    }
    hcalset = true;
  }

  if (m_do_hcalin)
  {
    if (!m_isdata)
    {
      m_waveforms_hcalin = findNode::getClass<TowerInfoContainer>(topNode, "WAVEFORM_HCALIN");

      if (!m_waveforms_hcalin)
      {
        std::cout << "No HCAL Waveforms found... " << std::endl;
        exit(1);
      }
    }

    m_primitives_hcalin = findNode::getClass<TriggerPrimitiveContainer>(topNode, "TRIGGERPRIMITIVES_HCALIN");

    if (!m_primitives_hcalin)
    {
      std::cout << "No HCALIN Primitives found... " << std::endl;
      exit(1);
    }

    if (!hcalset)
    {
      m_hcal_packets = findNode::getClass<CaloPacketContainer>(topNode, "HCALPackets");

      m_primitives_hcal_ll1 = findNode::getClass<TriggerPrimitiveContainer>(topNode, "TRIGGERPRIMITIVES_HCAL_LL1");

      if (!m_primitives_hcal_ll1)
      {
        std::cout << "No HCAL Primitives found... " << std::endl;
        exit(1);
      }
    }
  }
  if (m_do_emcal)
  {
    m_emcal_packets = findNode::getClass<CaloPacketContainer>(topNode, "CEMCPackets");

    if (!m_isdata)
    {
      m_waveforms_emcal = findNode::getClass<TowerInfoContainer>(topNode, "WAVEFORM_CEMC");

      if (!m_waveforms_emcal)
      {
        std::cout << "No EMCAL Waveforms found... " << std::endl;
        exit(1);
      }
    }
    m_primitives_emcal = findNode::getClass<TriggerPrimitiveContainer>(topNode, "TRIGGERPRIMITIVES_EMCAL");

    if (!m_primitives_emcal)
    {
      std::cout << "No EMCAL Primitives found... " << std::endl;
      exit(1);
    }

    m_primitives_emcal_ll1 = findNode::getClass<TriggerPrimitiveContainer>(topNode, "TRIGGERPRIMITIVES_EMCAL_LL1");

    if (!m_primitives_emcal_ll1)
    {
      std::cout << "No EMCAL 8x8 Primitives found... " << std::endl;
      exit(1);
    }

    m_primitives_emcal_2x2_ll1 = findNode::getClass<TriggerPrimitiveContainer>(topNode, "TRIGGERPRIMITIVES_EMCAL_2x2_LL1");

    if (!m_primitives_emcal_2x2_ll1)
    {
      std::cout << "No EMCAL 2x2 Primitives found... " << std::endl;
      exit(1);
    }
  }

  return;
}
void CaloTriggerEmulator::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing doing nothing" << std::endl;
  }

  PHCompositeNode *ll1Node = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "LL1"));
  if (!ll1Node)
  {
    ll1Node = new PHCompositeNode("LL1");
    dstNode->addNode(ll1Node);
  }
  {
    m_ll1_nodename = "LL1OUT_PHOTON";
    LL1Out *ll1out = findNode::getClass<LL1Out>(ll1Node, m_ll1_nodename);
    if (!ll1out)
    {
      ll1out = new LL1Outv1("PHOTON", "NONE");
      PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(ll1out, m_ll1_nodename, "PHObject");
      ll1Node->addNode(LL1OutNode);
    }
    m_prim_nodename = "TRIGGERPRIMITIVES_PHOTON";
    TriggerPrimitiveContainer *ll1out_prim = findNode::getClass<TriggerPrimitiveContainer>(ll1Node, m_prim_nodename);
    if (!ll1out_prim)
    {
      ll1out_prim = new TriggerPrimitiveContainerv1(TriggerDefs::TriggerId::photonTId, TriggerDefs::DetectorId::noneDId);
      PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(ll1out_prim, m_prim_nodename, "PHObject");
      ll1Node->addNode(LL1OutNode);
    }
  }
  {
    m_ll1_nodename = "LL1OUT_JET";
    LL1Out *ll1out = findNode::getClass<LL1Out>(ll1Node, m_ll1_nodename);
    if (!ll1out)
    {
      ll1out = new LL1Outv1("JET", "NONE");
      PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(ll1out, m_ll1_nodename, "PHObject");
      ll1Node->addNode(LL1OutNode);
    }
    m_prim_nodename = "TRIGGERPRIMITIVES_JET";
    TriggerPrimitiveContainer *ll1out_prim = findNode::getClass<TriggerPrimitiveContainer>(ll1Node, m_prim_nodename);
    if (!ll1out_prim)
    {
      ll1out_prim = new TriggerPrimitiveContainerv1(TriggerDefs::TriggerId::jetTId, TriggerDefs::DetectorId::noneDId);
      PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(ll1out_prim, m_prim_nodename, "PHObject");
      ll1Node->addNode(LL1OutNode);
    }
  }
  {
    m_ll1_nodename = "LL1OUT_PAIR";
    LL1Out *ll1out = findNode::getClass<LL1Out>(ll1Node, m_ll1_nodename);
    if (!ll1out)
    {
      ll1out = new LL1Outv1("PAIR", "NONE");
      PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(ll1out, m_ll1_nodename, "PHObject");
      ll1Node->addNode(LL1OutNode);
    }
    m_prim_nodename = "TRIGGERPRIMITIVES_PAIR";
    TriggerPrimitiveContainer *ll1out_prim = findNode::getClass<TriggerPrimitiveContainer>(ll1Node, m_prim_nodename);
    if (!ll1out_prim)
    {
      ll1out_prim = new TriggerPrimitiveContainerv1(TriggerDefs::TriggerId::pairTId, TriggerDefs::DetectorId::noneDId);
      PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(ll1out_prim, m_prim_nodename, "PHObject");
      ll1Node->addNode(LL1OutNode);
    }
  }
  if (m_do_emcal)
  {
    std::string ll1_nodename = "TRIGGERPRIMITIVES_EMCAL";
    TriggerPrimitiveContainer *ll1out_d1 = findNode::getClass<TriggerPrimitiveContainer>(ll1Node, ll1_nodename);
    if (!ll1out_d1)
    {
      ll1out_d1 = new TriggerPrimitiveContainerv1(TriggerDefs::TriggerId::noneTId, TriggerDefs::DetectorId::emcalDId);
      PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(ll1out_d1, ll1_nodename, "PHObject");
      ll1Node->addNode(LL1OutNode);
    }

    ll1_nodename = "TRIGGERPRIMITIVES_EMCAL_LL1";
    ll1out_d1 = findNode::getClass<TriggerPrimitiveContainer>(ll1Node, ll1_nodename);
    if (!ll1out_d1)
    {
      ll1out_d1 = new TriggerPrimitiveContainerv1(TriggerDefs::TriggerId::jetTId, TriggerDefs::DetectorId::emcalDId);
      PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(ll1out_d1, ll1_nodename, "PHObject");
      ll1Node->addNode(LL1OutNode);
    }

    ll1_nodename = "TRIGGERPRIMITIVES_EMCAL_2x2_LL1";
    ll1out_d1 = findNode::getClass<TriggerPrimitiveContainer>(ll1Node, ll1_nodename);
    if (!ll1out_d1)
    {
      ll1out_d1 = new TriggerPrimitiveContainerv1(TriggerDefs::TriggerId::pairTId, TriggerDefs::DetectorId::emcalDId);
      PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(ll1out_d1, ll1_nodename, "PHObject");
      ll1Node->addNode(LL1OutNode);
    }
  }
  if (m_do_hcalout)
  {
    std::string ll1_nodename = "TRIGGERPRIMITIVES_HCALOUT";
    TriggerPrimitiveContainer *ll1out_d = findNode::getClass<TriggerPrimitiveContainer>(ll1Node, ll1_nodename);
    if (!ll1out_d)
    {
      ll1out_d = new TriggerPrimitiveContainerv1(TriggerDefs::TriggerId::noneTId, TriggerDefs::DetectorId::hcaloutDId);
      PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(ll1out_d, ll1_nodename, "PHObject");
      ll1Node->addNode(LL1OutNode);
    }

    ll1_nodename = "TRIGGERPRIMITIVES_HCAL_LL1";
    TriggerPrimitiveContainer *ll1out_d1 = findNode::getClass<TriggerPrimitiveContainer>(ll1Node, ll1_nodename);
    if (!ll1out_d1)
    {
      ll1out_d1 = new TriggerPrimitiveContainerv1(TriggerDefs::TriggerId::jetTId, TriggerDefs::DetectorId::hcalDId);
      PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(ll1out_d1, ll1_nodename, "PHObject");
      ll1Node->addNode(LL1OutNode);
    }
  }
  if (m_do_hcalin)
  {
    std::string ll1_nodename = "TRIGGERPRIMITIVES_HCALIN";
    TriggerPrimitiveContainer *ll1out_d = findNode::getClass<TriggerPrimitiveContainer>(ll1Node, ll1_nodename);
    if (!ll1out_d)
    {
      ll1out_d = new TriggerPrimitiveContainerv1(TriggerDefs::TriggerId::noneTId, TriggerDefs::DetectorId::hcalinDId);
      PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(ll1out_d, ll1_nodename, "PHObject");
      ll1Node->addNode(LL1OutNode);
    }
    ll1_nodename = "TRIGGERPRIMITIVES_HCAL_LL1";
    TriggerPrimitiveContainer *ll1out_d1 = findNode::getClass<TriggerPrimitiveContainer>(ll1Node, ll1_nodename);
    if (!ll1out_d1)
    {
      ll1out_d1 = new TriggerPrimitiveContainerv1(TriggerDefs::TriggerId::jetTId, TriggerDefs::DetectorId::hcalDId);
      PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(ll1out_d1, ll1_nodename, "PHObject");
      ll1Node->addNode(LL1OutNode);
    }
  }
}

int CaloTriggerEmulator::End(PHCompositeNode * /*topNode*/)
{
  delete cdbttree_emcal;
  delete cdbttree_hcalout;
  delete cdbttree_hcalin;

  std::cout << "------------------------" << std::endl;
  std::cout << "Total Jet passed: " << m_jet_npassed << "/" << m_nevent << std::endl;
  std::cout << "Total Photon passed: " << m_photon_npassed << "/" << m_nevent << std::endl;
  std::cout << "Total Pair passed: " << m_pair_npassed << "/" << m_nevent << std::endl;
  std::cout << "------------------------" << std::endl;

  return 0;
}

void CaloTriggerEmulator::identify()
{
  std::cout << " CaloTriggerEmulator: " << m_trigger << std::endl;
  std::cout << " LL1Out: " << std::endl;

  if (m_ll1out_photon)
  {
    m_ll1out_photon->identify();
  }
  if (m_ll1out_jet)
  {
    m_ll1out_jet->identify();
  }
  if (m_ll1out_pair)
  {
    m_ll1out_pair->identify();
  }

  if (m_primitives_photon)
  {
    m_primitives_photon->identify();
  }

  if (m_primitives_jet)
  {
    m_primitives_jet->identify();
  }

  if (m_primitives_pair)
  {
    m_primitives_pair->identify();
  }

  if (m_primitives_emcal_ll1)
  {
    m_primitives_emcal_ll1->identify();
  }

  if (m_primitives_hcal_ll1)
  {
    m_primitives_hcal_ll1->identify();
  }

  if (m_primitives_emcal_2x2_ll1)
  {
    m_primitives_emcal_2x2_ll1->identify();
  }

  if (m_primitives_emcal)
  {
    m_primitives_emcal->identify();
  }

  if (m_primitives_hcalin)
  {
    m_primitives_hcalin->identify();
  }

  if (m_primitives_hcalout)
  {
    m_primitives_hcalout->identify();
  }

  std::cout << " Processed " << m_nevent << std::endl;
}

void CaloTriggerEmulator::useHCALIN(bool use)
{
  m_do_hcalin = use;
  m_force_hcalin = true;
}
void CaloTriggerEmulator::useHCALOUT(bool use)
{
  m_do_hcalout = use;
  m_force_hcalout = true;
}

void CaloTriggerEmulator::useEMCAL(bool use)
{
  m_do_emcal = use;
  m_force_emcal = true;
}

unsigned int CaloTriggerEmulator::getBits(unsigned int sum, TriggerDefs::TriggerId tid)
{
  if (tid == TriggerDefs::TriggerId::jetTId)
  {
    unsigned int bit = 0;
    for (unsigned int i = 0; i < 4; i++)
    {
      bit |= (sum >= m_threshold_jet[i] ? 0x1U << (i) : 0);
    }

    return bit;
  }
  else if (tid == TriggerDefs::TriggerId::pairTId)
  {
    unsigned int bit = 0;
    for (unsigned int i = 0; i < 4; i++)
    {
      bit |= (sum >= m_threshold_pair[i] ? 0x1U << (i) : 0);
    }

    return bit;
  }
  else if (tid == TriggerDefs::TriggerId::photonTId)
  {
    unsigned int bit = 0;
    for (unsigned int i = 0; i < 4; i++)
    {
      bit |= (sum >= m_threshold_photon[i] ? 0x1U << (i) : 0);
    }

    return bit;
  }
  return 0;
}
