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

#include <cdbobjects/CDBHistos.h>  // for CDBHistos

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

#include <TFile.h>
#include <TH1.h>
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
  // initialize all important parameters

  // default values for the lookup tables.
  // TODO: to CDB the LUTs from the database

  for (unsigned int i = 0; i < 1024; i++)
  {
    m_l1_adc_table[i] = (i) &0x3ffU;

    m_l1_adc_table_time[i] = (i) &0x3ffU;
  }
  for (unsigned int i = 0; i < 4096; i++)
  {
    m_l1_slewing_table[i] = (i) &0x1ffU;
  }

  // reset variables
  for (int j = 0; j < 8; j++)
  {
    m_trig_charge[j] = 0;
    for (auto &k : m2_trig_charge)
    {
      k[j] = 0;
    }
  }
  m_trig_nhit = 0;
  for (unsigned int &k : m2_trig_nhit)
  {
    k = 0;
  }
  for (int j = 0; j < 4; j++)
  {
    m_trig_time[j] = 0;
    for (int k = 0; k < 4; k++)
    {
      m2_trig_time[j][k] = 0;
    }
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

  for (int i = 0; i < 2; i++)
  {
    m_out_tsum[i] = 0;
    m_out_nhit[i] = 0;
    m_out_tavg[i] = 0;
    m_out_trem[i] = 0;
  }
  m_out_vtx_sub = 0;
  m_out_vtx_add = 0;

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
  m_ll1out = nullptr;

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
  m_threshold = 1;

  m_nhit1 = 2;
  m_nhit2 = 10;
  m_timediff1 = 10;
  m_timediff2 = 20;
  m_timediff3 = 30;

  // define a detector map for detectors included in a trigger
  m_det_map[TriggerDefs::TriggerId::noneTId] = {};
  m_det_map[TriggerDefs::TriggerId::jetTId] = {"EMCAL", "HCALIN", "HCALOUT"};
  m_det_map[TriggerDefs::TriggerId::mbdTId] = {"MBD"};
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
  m_prim_map[TriggerDefs::DetectorId::mbdDId] = 4;

  m_prim_ll1_map[TriggerDefs::TriggerId::jetTId] = 16;
  m_prim_ll1_map[TriggerDefs::TriggerId::pairTId] = 16;
  m_prim_ll1_map[TriggerDefs::TriggerId::photonTId] = 384;
  // booleans to control the input of detector data
  m_do_emcal = false;
  m_do_hcalin = false;
  m_do_hcalout = false;
  m_do_mbd = false;

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
    m_do_mbd = false;
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
    m_do_mbd = false;

    if (m_do_emcal && (m_do_hcalin || m_do_hcalout))
    {
      std::cout << "Cannot run PHOTON with hcal and emcal" << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  }
  else if (m_triggerid == TriggerDefs::TriggerId::mbdTId)
  {
    if (Verbosity() >= 2)
    {
      std::cout << "Using MBD Trigger." << std::endl;
    }
    if (!m_force_emcal)
    {
      m_do_emcal = false;
    }
    if (!m_force_hcalin)
    {
      m_do_hcalin = false;
    }
    if (!m_force_hcalout)
    {
      m_do_hcalout = false;
    }
    m_do_mbd = true;
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
    m_do_mbd = false;
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

    m_do_mbd = false;
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
    m_do_mbd = false;
  }
  else
  {
    std::cout << __FUNCTION__ << " : No trigger selected " << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_ll1_nodename = "LL1OUT_" + m_trigger;
  m_prim_nodename = "TRIGGERPRIMITIVES_" + m_trigger;

  // Get the calibrations and produce the lookup tables;

  if (Download_Calibrations())
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  CreateNodes(topNode);

  return 0;
}

int CaloTriggerEmulator::Download_Calibrations()
{
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
  m_peak_sub_ped_mbd.clear();

  return 0;
}

int CaloTriggerEmulator::process_waveforms()
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
      unsigned int peak_sub_ped = 0;
      TowerInfo *tower = m_waveforms_emcal->get_tower_at_channel(iwave);
      unsigned int key = TowerInfoDefs::encode_emcal(iwave);
      if (tower->get_nsample() == 2)
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
          int16_t maxim = tower->get_waveform_value(i);
          if (m_use_max)
          {
            int16_t max1 = std::max(tower->get_waveform_value(i), tower->get_waveform_value(i + 1));
            maxim = std::max(max1, tower->get_waveform_value(i + 2));
          }
          int subtraction = maxim - tower->get_waveform_value((i - m_trig_sub_delay > 0 ? i - m_trig_sub_delay : 0));
          // if negative, set to 0
          if (subtraction < 0)
          {
            subtraction = 0;
          }
          peak_sub_ped = (((unsigned int) subtraction) & 0x3fffU);
          if (Verbosity() >= 10 && peak_sub_ped > 16)
          {
            std::cout << __FILE__ << "::" << __FUNCTION__ << ":: emcal peak " << iwave << " = " << peak_sub_ped << std::endl;
          }

          v_peak_sub_ped.push_back(peak_sub_ped);
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
    unsigned int peak_sub_ped;

    std::vector<int> wave;
    // for each waveform, clauclate the peak - pedestal given the sub-delay setting
    if (!m_waveforms_hcalout->size())
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }

    for (unsigned int iwave = 0; iwave < (unsigned int) m_waveforms_hcalout->size(); iwave++)
    {
      std::vector<unsigned int> v_peak_sub_ped;
      peak_sub_ped = 0;
      TowerInfo *tower = m_waveforms_hcalout->get_tower_at_channel(iwave);
      unsigned int key = TowerInfoDefs::encode_hcal(iwave);
      if (tower->get_nsample() == 2)
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
          unsigned int subtraction = 0;
          int16_t maxim = tower->get_waveform_value(i);
          if (m_use_max)
          {
            int16_t max1 = std::max(tower->get_waveform_value(i), tower->get_waveform_value(i + 1));
            maxim = std::max(max1, tower->get_waveform_value(i + 2));
          }
          if (maxim - tower->get_waveform_value((i - m_trig_sub_delay > 0 ? i - m_trig_sub_delay : 0)) > 0)
          {
            subtraction = maxim - tower->get_waveform_value((i - m_trig_sub_delay > 0 ? i - m_trig_sub_delay : 0));
          }
          else
          {
            subtraction = 0;
          }

          peak_sub_ped = (((unsigned int) subtraction) & 0x3fffU);
          if (Verbosity() >= 10 && peak_sub_ped > 16)
          {
            std::cout << __FILE__ << "::" << __FUNCTION__ << ":: hcalout peak " << iwave << " = " << peak_sub_ped << std::endl;
          }

          v_peak_sub_ped.push_back(peak_sub_ped);
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
    unsigned int peak_sub_ped;

    std::vector<unsigned int> wave;

    // for each waveform, clauclate the peak - pedestal given the sub-delay setting
    for (unsigned int iwave = 0; iwave < (unsigned int) m_waveforms_hcalin->size(); iwave++)
    {
      std::vector<unsigned int> v_peak_sub_ped;
      peak_sub_ped = 0;
      TowerInfo *tower = m_waveforms_hcalin->get_tower_at_channel(iwave);
      unsigned int key = TowerInfoDefs::encode_hcal(iwave);

      if (tower->get_nsample() == 2)
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
          int16_t maxim = tower->get_waveform_value(i);
          if (m_use_max)
          {
            int16_t max1 = std::max(tower->get_waveform_value(i), tower->get_waveform_value(i + 1));
            maxim = std::max(max1, tower->get_waveform_value(i + 2));
          }

          int subtraction = maxim - tower->get_waveform_value((i - m_trig_sub_delay > 0 ? i - m_trig_sub_delay : 0));

          if (subtraction < 0)
          {
            subtraction = 0;
          }

          peak_sub_ped = (((unsigned int) subtraction) & 0x3fffU);
          if (Verbosity() >= 10 && peak_sub_ped > 16)
          {
            std::cout << __FILE__ << "::" << __FUNCTION__ << ":: hcalin peak " << iwave << " = " << peak_sub_ped << std::endl;
          }

          v_peak_sub_ped.push_back(peak_sub_ped);
        }
      }
      // save in global.
      m_peak_sub_ped_hcalin[key] = v_peak_sub_ped;
    }
  }

  if (m_do_mbd)
  {
    if (!m_waveforms_mbd->size())
    {
      return Fun4AllReturnCodes::EVENT_OK;
    }

    // for each waveform, clauclate the peak - pedestal given the sub-delay setting
    for (unsigned int iwave = 0; iwave < (unsigned int) m_waveforms_mbd->size(); iwave++)
    {
      std::vector<unsigned int> v_peak_sub_ped;
      unsigned int peak_sub_ped = 0;
      TowerInfo *tower = m_waveforms_mbd->get_tower_at_channel(iwave);
      for (int i = sample_start; i < sample_end; i++)
      {
        int subtraction = tower->get_waveform_value(i) - tower->get_waveform_value((i - i % 6 - 6 + m_trig_sub_delay > 0 ? i - i % 6 - 6 + m_trig_sub_delay : 0));

        if (subtraction < 0)
        {
          subtraction = 0;
        }

        peak_sub_ped = (((unsigned int) subtraction) & 0x3fffU);
        v_peak_sub_ped.push_back(peak_sub_ped);
      }
      // save in global.
      m_peak_sub_ped_mbd[iwave] = v_peak_sub_ped;
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

  if (m_do_mbd)
  {
    // MBD

    ip = 0;

    // get number of primitives
    m_n_primitives = m_prim_map[TriggerDefs::DetectorId::mbdDId];

    for (i = 0; i < m_n_primitives; i++, ip++)
    {
      // make primitive key
      TriggerDefs::TriggerPrimKey primkey = TriggerDefs::getTriggerPrimKey(TriggerDefs::GetTriggerId(m_trigger), TriggerDefs::GetDetectorId("MBD"), TriggerDefs::GetPrimitiveId("MBD"), m_n_primitives - ip);

      // make primitive and check mask;
      TriggerPrimitive *primitive = m_primitives->get_primitive_at_key(primkey);
      mask = CheckFiberMasks(primkey);

      // iterate through samples
      for (int is = 0; is < nsample; is++)
      {
        // reset variables
        for (unsigned int &j : m_trig_charge)
        {
          j = 0;
        }
        m_trig_nhit = 0;
        for (unsigned int &j : m_trig_time)
        {
          j = 0;
        }

        unsigned int tmp, tmp2;
        unsigned int qadd[32];

        // for each section of the board (4 sections of 8 time and 8 charge
        for (int isec = 0; isec < 4; isec++)
        {
          // go through 8 charge channels
          for (int j = 0; j < 8; j++)
          {
            // pass upper 10 bits of charge to get 10 bit LUt outcome
            tmp = m_l1_adc_table[m_peak_sub_ped_mbd[i * 64 + 8 + isec * 16 + j].at(is) >> 4U];

            // put upper 3 bits of the 10 bits into slewing correction later
            qadd[isec * 8 + j] = (tmp & 0x380U) >> 7U;

            // sum up to 11 bits.
            m_trig_charge[isec * 2 + j / 4] += tmp & 0x7ffU;
          }
        }

        // Now the time channels
        for (int isec = 0; isec < 4; isec++)
        {
          // 8 timing channels
          for (int j = 0; j < 8; j++)
          {
            // upper 10 bits go through the LUT
            tmp = m_l1_adc_table[m_peak_sub_ped_mbd[i * 64 + isec * 16 + j].at(is) >> 4U];

            // high bit is the hit bit
            m_trig_nhit += (tmp & 0x200U) >> 9U;

            // get upper 3 bits of charge in the channel, and make it bits 9-11, the time of the chanel is the lower 9 bits from 0-8.
            tmp2 = m_l1_slewing_table[(qadd[isec * 8 + j] << 9U) + (tmp & 0x01ffU)];

            // attribute to the time sum
            m_trig_time[isec] += tmp2;
          }
        }

        // ad in the charge sums

        for (int j = 0; j < 13; j++)
        {
          TriggerDefs::TriggerSumKey sumkey = TriggerDefs::getTriggerSumKey(TriggerDefs::GetTriggerId(m_trigger), TriggerDefs::GetDetectorId("MBD"), TriggerDefs::GetPrimitiveId("MBD"), ip, j);
          if (j < 8)
          {
            primitive->get_sum_at_key(sumkey)->push_back(m_trig_charge[j]);
          }
          else if (j == 8)
          {
            primitive->get_sum_at_key(sumkey)->push_back(m_trig_nhit);
          }
          else
          {
            primitive->get_sum_at_key(sumkey)->push_back(m_trig_time[j - 9]);
          }
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
  std::vector<unsigned int> bits;

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
  if (m_triggerid == TriggerDefs::TriggerId::photonTId || m_triggerid == TriggerDefs::TriggerId::jetTId)
  {
    if (Verbosity() >= 2)
    {
      std::cout << __FUNCTION__ << " " << __LINE__ << " processing 8x8 non-overlapping sums" << std::endl;
    }

    // Make the jet primitives
    if (m_do_emcal)
    {
      if (!m_primitives_emcal)
      {
        std::cout << "There is no primitive container" << std::endl;
        return Fun4AllReturnCodes::EVENT_OK;
      }
      {
        TriggerPrimitiveContainer::Range range = m_primitives_emcal_ll1->getTriggerPrimitives();
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

      // iterate through emcal primitives and organize into the 16 jet primitives each with the 8x8 nonoverlapping sum
      TriggerPrimitiveContainer::Range range = m_primitives_emcal->getTriggerPrimitives();

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
          it_s = ((it_s >> 4U) & 0xffU);
        }
      }
    }
    // Make the jet primitives for hcal
    if (m_do_hcalin || m_do_hcalout)
    {
      // these are the 16 inputs from the HCAL, make them first, with 0s, then add HCALIN and HCALOUT and then divide by 2.
      {
        TriggerPrimitiveContainer::Range range = m_primitives_hcal_ll1->getTriggerPrimitives();
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

      if (m_do_hcalin)
      {
        if (!m_primitives_hcalin)
        {
          std::cout << "There is no primitive container" << std::endl;
          return Fun4AllReturnCodes::EVENT_OK;
        }

        // iterate through emcal primitives and organize into the 16 jet primitives each with the 8x8 nonoverlapping sum
        TriggerPrimitiveContainer::Range range = m_primitives_hcalin->getTriggerPrimitives();

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

      if (m_do_hcalout)
      {
        if (!m_primitives_hcalout)
        {
          std::cout << "There is no primitive container" << std::endl;
          return Fun4AllReturnCodes::EVENT_OK;
        }

        // iterate through emcal primitives and organize into the 16 jet primitives each with the 8x8 nonoverlapping sum
        TriggerPrimitiveContainerv1::Range range = m_primitives_hcalout->getTriggerPrimitives();

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
            TriggerDefs::TriggerPrimKey jet_prim_key = TriggerDefs::getTriggerPrimKey(m_triggerid, TriggerDefs::GetDetectorId("HCAL"), TriggerDefs::GetPrimitiveId("JET"), iprim);

            TriggerDefs::TriggerPrimKey jet_sum_key = TriggerDefs::getTriggerSumKey(m_triggerid, TriggerDefs::GetDetectorId("HCAL"), TriggerDefs::GetPrimitiveId("JET"), iprim, isum);

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

      TriggerPrimitiveContainerv1::Range range = m_primitives_hcal_ll1->getTriggerPrimitives();
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
    }

    if (m_triggerid == TriggerDefs::TriggerId::jetTId)
    {
      {
        TriggerPrimitiveContainer::Range range = m_primitives->getTriggerPrimitives();
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

      TriggerPrimitiveContainerv1::Range range = m_primitives->getTriggerPrimitives();
      for (TriggerPrimitiveContainerv1::Iter iter = range.first; iter != range.second; ++iter)
      {
        TriggerDefs::TriggerPrimKey jet_pkey = (*iter).first;
        TriggerDefs::TriggerPrimKey hcal_pkey = TriggerDefs::getTriggerPrimKey(m_triggerid, TriggerDefs::GetDetectorId("HCAL"), TriggerDefs::GetPrimitiveId("JET"), TriggerDefs::getPrimitiveLocId_from_TriggerPrimKey(jet_pkey));
        TriggerDefs::TriggerPrimKey emcal_pkey = TriggerDefs::getTriggerPrimKey(m_triggerid, TriggerDefs::GetDetectorId("EMCAL"), TriggerDefs::GetPrimitiveId("JET"), TriggerDefs::getPrimitiveLocId_from_TriggerPrimKey(jet_pkey));
        TriggerPrimitive *primitive = (*iter).second;
        TriggerPrimitivev1::Range sumrange = primitive->getSums();
        for (TriggerPrimitivev1::Iter iter_sum = sumrange.first; iter_sum != sumrange.second; ++iter_sum)
        {
          TriggerDefs::TriggerSumKey jet_skey = (*iter_sum).first;
          TriggerDefs::TriggerSumKey hcal_skey = TriggerDefs::getTriggerSumKey(m_triggerid, TriggerDefs::GetDetectorId("HCAL"), TriggerDefs::GetPrimitiveId("JET"), TriggerDefs::getPrimitiveLocId_from_TriggerPrimKey(jet_pkey), TriggerDefs::getSumLocId(jet_skey));
          TriggerDefs::TriggerSumKey emcal_skey = TriggerDefs::getTriggerSumKey(m_triggerid, TriggerDefs::GetDetectorId("EMCAL"), TriggerDefs::GetPrimitiveId("JET"), TriggerDefs::getPrimitiveLocId_from_TriggerPrimKey(jet_pkey), TriggerDefs::getSumLocId(jet_skey));

          int i = 0;
          for (unsigned int &it_s : *(*iter_sum).second)
          {
            unsigned int sum_hcal = m_primitives_hcal_ll1->get_primitive_at_key(hcal_pkey)->get_sum_at_key(hcal_skey)->at(i);
            unsigned int sum_emcal = m_primitives_emcal_ll1->get_primitive_at_key(emcal_pkey)->get_sum_at_key(emcal_skey)->at(i);

            it_s = ((sum_hcal + sum_emcal) >> 1U);
            i++;
          }
        }
      }
    }
  }
  else if (m_triggerid == TriggerDefs::TriggerId::pairTId)
  {
    if (!m_primitives_emcal)
    {
      std::cout << "There is no primitive container" << std::endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }
    // Make the jet primitives
    TriggerPrimitiveContainerv1::Range range;
    TriggerPrimitivev1::Range sumrange;
    // iterate through emcal primitives and organize into the 16 jet primitives each with the 8x8 nonoverlapping sum
    range = m_primitives_emcal->getTriggerPrimitives();

    TriggerDefs::TriggerSumKey temp_sum_key;
    TriggerDefs::TriggerSumKey temp_prim_key;
    for (TriggerPrimitiveContainerv1::Iter iter = range.first; iter != range.second; ++iter)
    {
      // get key and see if masked
      TriggerDefs::TriggerPrimKey key = (*iter).first;
      if (CheckFiberMasks(key))
      {
        if (Verbosity() >= 2)
        {
          std::cout << "masked: " << key << std::endl;
        }
        continue;
      }

      uint16_t primlocid = TriggerDefs::getPrimitiveLocId_from_TriggerPrimKey(key);
      bool prim_right_edge = (primlocid % 12 == 11);
      uint16_t topedge_primlocid = (primlocid / 12 == 31 ? primlocid % 12 : primlocid + 12);

      TriggerDefs::TriggerPrimKey primkey = TriggerDefs::getTriggerPrimKey(m_triggerid, TriggerDefs::GetDetectorId("EMCAL"), TriggerDefs::GetPrimitiveId("PAIR"), primlocid);
      TriggerPrimitive *primitive_photon = m_primitives->get_primitive_at_key(primkey);

      // get the primitive (16 2x2 sums)
      TriggerPrimitive *primitive = (*iter).second;
      // in this primitive we will hold a 4x4 overlapping sum for each 2x2 sum in the original primitive
      unsigned int number_of_sums = primitive->size();
      // iterate through all sums and calculate 4x4 overlapping sum
      for (unsigned int isum = 0; isum < number_of_sums; isum++)
      {
        bool right_edge = (isum % 4 == 3);
        bool top_edge = (isum / 4 == 3);

        if (right_edge && prim_right_edge)
        {
          continue;
        }
        TriggerDefs::TriggerSumKey sumkey = TriggerDefs::getTriggerSumKey(m_triggerid, TriggerDefs::GetDetectorId("EMCAL"), TriggerDefs::GetPrimitiveId("PAIR"), primlocid, isum);
        std::vector<unsigned int> *t_sum = primitive_photon->get_sum_at_key(sumkey);
        ;
        for (int is = 0; is < nsample; is++)
        {
          unsigned int sum = 0;

          temp_sum_key = TriggerDefs::getTriggerSumKey(m_triggerid, TriggerDefs::DetectorId::emcalDId, TriggerDefs::PrimitiveId::calPId, primlocid, isum);
          sum += (primitive->get_sum_at_key(temp_sum_key)->at(is) & 0xffU);
          if (right_edge)
          {
            temp_prim_key = TriggerDefs::getTriggerPrimKey(m_triggerid, TriggerDefs::DetectorId::emcalDId, TriggerDefs::PrimitiveId::calPId, primlocid + 1);
            temp_sum_key = TriggerDefs::getTriggerSumKey(m_triggerid, TriggerDefs::DetectorId::emcalDId, TriggerDefs::PrimitiveId::calPId, primlocid + 1, (isum / 4) * 4);
            sum += (m_primitives_emcal->get_primitive_at_key(temp_prim_key)->get_sum_at_key(temp_sum_key)->at(is) & 0xffU);
          }
          else
          {
            temp_sum_key = TriggerDefs::getTriggerSumKey(m_triggerid, TriggerDefs::DetectorId::emcalDId, TriggerDefs::PrimitiveId::calPId, primlocid, isum + 1);
            sum += (primitive->get_sum_at_key(temp_sum_key)->at(is) & 0xffU);
          }
          if (top_edge)
          {
            temp_prim_key = TriggerDefs::getTriggerPrimKey(m_triggerid, TriggerDefs::DetectorId::emcalDId, TriggerDefs::PrimitiveId::calPId, topedge_primlocid);
            temp_sum_key = TriggerDefs::getTriggerSumKey(m_triggerid, TriggerDefs::DetectorId::emcalDId, TriggerDefs::PrimitiveId::calPId, topedge_primlocid, isum % 4);
            sum += (m_primitives_emcal->get_primitive_at_key(temp_prim_key)->get_sum_at_key(temp_sum_key)->at(is) & 0xffU);
          }
          else
          {
            temp_sum_key = TriggerDefs::getTriggerSumKey(m_triggerid, TriggerDefs::DetectorId::emcalDId, TriggerDefs::PrimitiveId::calPId, primlocid, isum + 4);
            sum += (primitive->get_sum_at_key(temp_sum_key)->at(is) & 0xffU);
          }

          if (top_edge && right_edge)
          {
            temp_prim_key = TriggerDefs::getTriggerPrimKey(m_triggerid, TriggerDefs::DetectorId::emcalDId, TriggerDefs::PrimitiveId::calPId, topedge_primlocid + 1);
            temp_sum_key = TriggerDefs::getTriggerSumKey(m_triggerid, TriggerDefs::DetectorId::emcalDId, TriggerDefs::PrimitiveId::calPId, topedge_primlocid + 1, 0);
            sum += (m_primitives_emcal->get_primitive_at_key(temp_prim_key)->get_sum_at_key(temp_sum_key)->at(is) & 0xffU);
          }
          else if (top_edge)
          {
            temp_prim_key = TriggerDefs::getTriggerPrimKey(m_triggerid, TriggerDefs::DetectorId::emcalDId, TriggerDefs::PrimitiveId::calPId, topedge_primlocid);
            temp_sum_key = TriggerDefs::getTriggerSumKey(m_triggerid, TriggerDefs::DetectorId::emcalDId, TriggerDefs::PrimitiveId::calPId, topedge_primlocid, isum % 4 + 1);
            sum += (m_primitives_emcal->get_primitive_at_key(temp_prim_key)->get_sum_at_key(temp_sum_key)->at(is) & 0xffU);
          }
          else if (right_edge)
          {
            temp_prim_key = TriggerDefs::getTriggerPrimKey(m_triggerid, TriggerDefs::DetectorId::emcalDId, TriggerDefs::PrimitiveId::calPId, primlocid + 1);
            temp_sum_key = TriggerDefs::getTriggerSumKey(m_triggerid, TriggerDefs::DetectorId::emcalDId, TriggerDefs::PrimitiveId::calPId, primlocid + 1, (isum / 4 + 1) * 4);
            sum += (m_primitives_emcal->get_primitive_at_key(temp_prim_key)->get_sum_at_key(temp_sum_key)->at(is) & 0xffU);
          }
          else
          {
            temp_sum_key = TriggerDefs::getTriggerSumKey(m_triggerid, TriggerDefs::DetectorId::emcalDId, TriggerDefs::PrimitiveId::calPId, primlocid, isum + 5);
            sum += (primitive->get_sum_at_key(temp_sum_key)->at(is) & 0xffU);
          }

          sum = (sum >> 2U);
          t_sum->push_back(sum);
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

  std::vector<unsigned int> *trig_bits = m_ll1out->GetTriggerBits();

  // photon
  // 8x8 non-overlapping sums in the EMCAL
  // create the 8x8 non-overlapping sum

  if (m_triggerid == TriggerDefs::TriggerId::photonTId)
  {
    if (Verbosity() >= 2)
    {
      std::cout << __FUNCTION__ << " " << __LINE__ << " processing PHOTON trigger , bits before: " << trig_bits->size() << std::endl;
    }

    // Make the jet primitives
    if (m_do_emcal)
    {
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

          std::vector<unsigned int> *t_sum = (*iter_sum).second;
          TriggerDefs::TriggerSumKey sumk = (*iter_sum).first;
          for (int is = 0; is < nsample; is++)
          {
            unsigned int bit = getBits(t_sum->at(is));
            if (bit)
            {
              m_ll1out->addTriggeredSum(sumk);
              m_ll1out->addTriggeredPrimitive(key);
            }
            bits.at(is) |= bit;
          }
        }
      }
    }
    // Make the jet primitives for hcal
    else
    {
      TriggerPrimitiveContainer::Range range = m_primitives_hcal_ll1->getTriggerPrimitives();

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

          std::vector<unsigned int> *t_sum = (*iter_sum).second;
          TriggerDefs::TriggerSumKey sumk = (*iter_sum).first;
          for (int is = 0; is < nsample; is++)
          {
            unsigned int bit = getBits(t_sum->at(is));
            if (bit)
            {
              m_ll1out->addTriggeredSum(sumk);
              m_ll1out->addTriggeredPrimitive(key);
            }
            bits.at(is) |= bit;
          }
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
      m_npassed++;
    }
  }
  else if (m_triggerid == TriggerDefs::TriggerId::jetTId)
  {
    // Make the jet primitives

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

    if (!m_primitives)
    {
      std::cout << "There is no primitive container" << std::endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }

    // iterate through emcal primitives and organize into the 16 jet primitives each with the 8x8 nonoverlapping sum
    TriggerPrimitiveContainer::Range range = m_primitives->getTriggerPrimitives();

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
        for (unsigned int &it_s : *(*iter_sum).second)
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

    int pass = 0;
    for (int ijphi = 0; ijphi < 32; ijphi++)
    {
      for (int ijeta = 0; ijeta < 9; ijeta++)
      {
        unsigned int sk = ((unsigned int) ijphi & 0xffffU) + (((unsigned int) ijeta & 0xffffU) << 16U);
        std::vector<unsigned int> *sum = m_ll1out->get_word(sk);
        for (int is = 0; is < nsample; is++)
        {
          sum->push_back(jet_map[ijphi][ijeta].at(is));
          unsigned int bit = getBits(jet_map[ijphi][ijeta].at(is));
          if (bit)
          {
            m_ll1out->addTriggeredSum(sk);
            m_ll1out->addTriggeredPrimitive(sk);
            pass = 1;
          }
          bits.at(is) |= bit;
        }
      }
    }

    for (int is = 0; is < nsample; is++)
    {
      trig_bits->push_back(bits.at(is));
    }

    if (pass)
    {
      m_npassed++;
    }
  }
  // pair
  // 4x4 overlapping sum
  else if (m_triggerid == TriggerDefs::TriggerId::pairTId)
  {
  }

  // cosmic trigger (singles)
  else if (m_triggerid == TriggerDefs::TriggerId::cosmicTId)
  {
    if (!m_primitives_hcalout)
    {
      std::cout << "There is no primitive container" << std::endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }

    // iterating through the trigger primitives, and seeing if ANY is above threshold.
    TriggerPrimitiveContainerv1::Range range;
    if (m_do_hcalout)
    {
      range = m_primitives_hcalout->getTriggerPrimitives();

      for (TriggerPrimitiveContainerv1::Iter iter = range.first; iter != range.second; ++iter)
      {
        TriggerDefs::TriggerPrimKey key = (*iter).first;
        if (CheckFiberMasks(key))
        {
          if (Verbosity() >= 2)
          {
            std::cout << "masked: " << key << std::endl;
          }
          continue;
        }

        TriggerPrimitive *primitive = (*iter).second;
        TriggerPrimitivev1::Range sumrange = primitive->getSums();
        for (TriggerPrimitivev1::Iter iter_sum = sumrange.first; iter_sum != sumrange.second; ++iter_sum)
        {
          TriggerDefs::TriggerSumKey sumk = (*iter_sum).first;
          std::vector<unsigned int> *t_sum = (*iter_sum).second;

          if (CheckChannelMasks(sumk))
          {
            continue;
          }

          for (int is = 0; is < nsample; is++)
          {
            unsigned int bit = getBits(t_sum->at(is));
            if (bit)
            {
              m_ll1out->addTriggeredSum(sumk);
              m_ll1out->addTriggeredPrimitive(key);
            }
            bits.at(is) |= bit;
          }
        }
      }
    }

    if (m_do_hcalin)
    {
      range = m_primitives_hcalin->getTriggerPrimitives();
      for (TriggerPrimitiveContainerv1::Iter iter = range.first; iter != range.second; ++iter)
      {
        TriggerDefs::TriggerPrimKey key = (*iter).first;
        if (CheckFiberMasks(key))
        {
          if (Verbosity() >= 2)
          {
            std::cout << "masked: " << key << std::endl;
          }
          continue;
        }
        TriggerPrimitive *primitive = (*iter).second;
        TriggerPrimitivev1::Range sumrange = primitive->getSums();
        for (TriggerPrimitivev1::Iter iter_sum = sumrange.first; iter_sum != sumrange.second; ++iter_sum)
        {
          TriggerDefs::TriggerSumKey sumk = (*iter_sum).first;

          if (CheckChannelMasks(sumk))
          {
            continue;
          }
          std::vector<unsigned int> *t_sum = (*iter_sum).second;
          for (int is = 0; is < nsample; is++)
          {
            unsigned int bit = getBits(t_sum->at(is));

            if (bit)
            {
              m_ll1out->addTriggeredSum(sumk);
              m_ll1out->addTriggeredPrimitive(key);
            }
            bits.at(is) |= bit;
          }
        }
      }
    }
    // check if any sample passes here.

    int pass = 0;
    for (int is = 0; is < nsample; is++)
    {
      std::cout << bits.at(is) << std::endl;
      trig_bits->push_back(bits.at(is));
      if (trig_bits->at(is))
      {
        pass = 1;
      }
    }
    m_npassed += pass;
  }

  // cosmic (coincidence)
  else if (m_triggerid == TriggerDefs::TriggerId::cosmic_coinTId)
  {
    // organize the sums
    unsigned int cosmic_organized_sums[2][12][32];

    if (!m_primitives_hcalout || !m_primitives_hcalin)
    {
      std::cout << "There is no primitive container" << std::endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }
    uint16_t icard, icosmic;
    TriggerPrimitiveContainerv1::Range range;
    for (int isam = 0; isam < nsample; isam++)
    {
      //	  Set everything to 0 to get the sums.
      for (auto &cosmic_organized_sum : cosmic_organized_sums)
      {
        for (auto &iii : cosmic_organized_sum)
        {
          for (unsigned int &iv : iii)
          {
            iv = 0;
          }
        }
      }

      // get all primitives and iterate
      range = m_primitives_hcalout->getTriggerPrimitives();
      for (TriggerPrimitiveContainerv1::Iter iter = range.first; iter != range.second; ++iter)
      {
        // get key
        TriggerDefs::TriggerPrimKey key = (*iter).first;
        // check if primitive is masked
        if (CheckFiberMasks(key))
        {
          if (Verbosity() >= 2)
          {
            std::cout << "masked: " << key << std::endl;
          }
          continue;
        }

        TriggerPrimitive *primitive = (*iter).second;
        // get location in index of phi and eta
        uint16_t primphi = TriggerDefs::getPrimitivePhiId_from_TriggerPrimKey(key);
        uint16_t primeta = TriggerDefs::getPrimitiveEtaId_from_TriggerPrimKey(key);

        // get the card (either 0 or 1);
        icard = (primphi < 4 ? 0 : 1);
        TriggerPrimitivev1::Range sumrange = primitive->getSums();
        // if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<" key: "<<key<<" size: "<<primitive->size()<<std::endl;
        for (TriggerPrimitivev1::Iter iter_sum = sumrange.first; iter_sum != sumrange.second; ++iter_sum)
        {
          // get sum key
          TriggerDefs::TriggerSumKey sumkey = (*iter_sum).first;
          if (CheckChannelMasks(sumkey))
          {
            continue;
          }
          // get the integer index in phi and eta of the usm within the 8x8 area (4x4 sums).
          uint16_t sumphi = TriggerDefs::getSumPhiId(sumkey);
          uint16_t sumeta = TriggerDefs::getSumEtaId(sumkey);

          // get the cosmic area
          icosmic = (primeta * 4 + (3 - sumeta)) / 2;
          int isum = ((primphi % 4) * 8) + sumphi * 2 + (1 - (sumeta % 2));

          cosmic_organized_sums[icard][icosmic][isum] = *((*iter_sum).second->begin() + isam);
        }
      }

      // no the same for inner hcal.
      range = m_primitives_hcalin->getTriggerPrimitives();
      for (TriggerPrimitiveContainerv1::Iter iter = range.first; iter != range.second; ++iter)
      {
        TriggerDefs::TriggerPrimKey key = (*iter).first;
        if (CheckFiberMasks(key))
        {
          if (Verbosity() >= 2)
          {
            std::cout << "masked: " << key << std::endl;
          }
          continue;
        }

        TriggerPrimitive *primitive = (*iter).second;

        uint16_t primphi = TriggerDefs::getPrimitivePhiId_from_TriggerPrimKey(key);
        uint16_t primeta = TriggerDefs::getPrimitiveEtaId_from_TriggerPrimKey(key);
        icard = (primphi < 4 ? 0 : 1);
        TriggerPrimitivev1::Range sumrange = primitive->getSums();
        // if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<" key: "<<key<<" size: "<<primitive->size()<<std::endl;
        for (TriggerPrimitivev1::Iter iter_sum = sumrange.first; iter_sum != sumrange.second; ++iter_sum)
        {
          TriggerDefs::TriggerSumKey sumkey = (*iter_sum).first;

          uint16_t sumphi = TriggerDefs::getSumPhiId(sumkey);
          uint16_t sumeta = TriggerDefs::getSumEtaId(sumkey);

          icosmic = 6 + (primeta * 4 + (3 - sumeta)) / 2;
          int isum = ((primphi % 4) * 8) + sumphi * 2 + (1 - (sumeta % 2));
          if (CheckChannelMasks(sumkey))
          {
            continue;
          }
          cosmic_organized_sums[icard][icosmic][isum] = *((*iter_sum).second->begin() + isam);
        }
      }

      // for the two cards see if there is a coincidence.
      unsigned int hit_cosmic[2] = {0, 0};
      for (unsigned int ica = 0; ica < 2; ica++)
      {
        for (unsigned int ic = 0; ic < 6; ic++)
        {
          for (unsigned int isum = 0; isum < 32; isum++)
          {
            hit_cosmic[ica] |= ((cosmic_organized_sums[ica][ic][isum] > (unsigned int) m_threshold ? 0x1U : 0) << ic);
            hit_cosmic[ica] |= ((cosmic_organized_sums[ica][ic + 6][isum] > (unsigned int) m_threshold ? 0x1U : 0) << (6U + ic));
          }
        }
      }
      if (((m_l1_hcal_table[hit_cosmic[0]] & 0x1U) == 0x1U) && ((m_l1_hcal_table[hit_cosmic[1]] & 0x2U) == 0x2U))
      {
        bits.at(isam) |= 1U;
      }
      if (((m_l1_hcal_table[hit_cosmic[0]] & 0x2U) == 0x2U) && ((m_l1_hcal_table[hit_cosmic[1]] & 0x1U) == 0x1U))
      {
        bits.at(isam) |= 1U;
      }
    }

    int pass = 0;
    for (int is = 0; is < nsample; is++)
    {
      trig_bits->push_back(bits.at(is));
      if (bits.at(is) == 1)
      {
        pass = 1;
      }
    }
    m_npassed += pass;
  }

  // this is the MBD trigger algorithm
  else if (m_triggerid == TriggerDefs::TriggerId::mbdTId)
  {
    {
      std::cout << "There is no primitive container" << std::endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }

    TriggerPrimitiveContainerv1::Range range;
    TriggerPrimitivev1::Range sumrange;
    int ip, isum;

    range = m_primitives->getTriggerPrimitives();

    if (Verbosity() >= 2)
    {
      std::cout << __FUNCTION__ << " " << __LINE__ << " mbd primitives size: " << m_primitives->size() << std::endl;
    }

    std::vector<unsigned int> *word_mbd = nullptr;

    m_word_mbd.clear();
    for (int j = 0; j < 8; j++)
    {
      word_mbd = new std::vector<unsigned int>();
      m_word_mbd.push_back(word_mbd);
    }

    for (int is = 0; is < nsample; is++)
    {
      ip = 0;
      for (TriggerPrimitiveContainerv1::Iter iter = range.first; iter != range.second; ++iter, ip++)
      {
        TriggerPrimitive *primitive = (*iter).second;
        sumrange = primitive->getSums();
        isum = 0;
        for (TriggerPrimitivev1::Iter iter_sum = sumrange.first; iter_sum != sumrange.second; ++iter_sum, isum++)
        {
          if (isum < 8)
          {
            m2_trig_charge[ip][isum] = ((*iter_sum).second)->at(is);
          }
          else if (isum == 8)
          {
            m2_trig_nhit[ip] = ((*iter_sum).second)->at(is);
          }
          else
          {
            m2_trig_time[ip][isum - 9] = ((*iter_sum).second)->at(is);
          }
        }
      }

      if (Verbosity() && is == 11)
      {
        for (int q = 0; q < 8; q++)
        {
          std::cout << "Q" << std::dec << q << ": ";
          for (auto &ipp : m2_trig_charge)
          {
            std::cout << std::hex << ipp[q] << " ";
          }
          std::cout << " " << std::endl;
        }
        std::cout << "NH: ";
        for (unsigned int ipp : m2_trig_nhit)
        {
          std::cout << std::hex << ipp << " ";
        }
        std::cout << " " << std::endl;

        for (int q = 0; q < 4; q++)
        {
          std::cout << "T" << std::dec << q << ": ";
          for (auto &ipp : m2_trig_time)
          {
            std::cout << std::hex << ipp[q] << " ";
          }
          std::cout << " " << std::endl;
        }
      }

      m_out_tsum[0] = 0;
      m_out_tsum[1] = 0;
      m_out_nhit[0] = 0;
      m_out_nhit[1] = 0;
      m_out_tavg[0] = 0;
      m_out_tavg[1] = 0;
      m_out_trem[0] = 0;
      m_out_trem[1] = 0;
      m_out_vtx_sub = 0;
      m_out_vtx_add = 0;

      for (int i = 0; i < 2; i++)
      {
        for (int j = 0; j < 4; j++)
        {
          m_out_tsum[0] += m2_trig_time[i][j];
          m_out_tsum[1] += m2_trig_time[i + 2][j];
        }
        m_out_nhit[0] += m2_trig_nhit[i];
        m_out_nhit[1] += m2_trig_nhit[i + 2];
      }

      if (m_out_nhit[0] != 0)
      {
        m_out_tavg[0] = m_out_tsum[0] / m_out_nhit[0];
        m_out_trem[0] = m_out_tsum[0] % m_out_nhit[0];
      }
      if (m_out_nhit[1] != 0)
      {
        m_out_tavg[1] = m_out_tsum[1] / m_out_nhit[1];
        m_out_trem[1] = m_out_tsum[1] % m_out_nhit[1];
      }

      unsigned int max = m_out_tavg[0];
      unsigned int min = m_out_tavg[1];
      if (min > max)
      {
        max = m_out_tavg[1];
        min = m_out_tavg[0];
      }

      m_out_vtx_sub = (max - min) & 0x1ffU;
      m_out_vtx_add = (m_out_tavg[0] + m_out_tavg[1]) & 0x3ffU;

      m_word_mbd[0]->push_back(m_out_tavg[0]);
      m_word_mbd[1]->push_back(m_out_tavg[1]);
      m_word_mbd[2]->push_back(m_out_nhit[0]);
      m_word_mbd[3]->push_back(m_out_nhit[1]);
      m_word_mbd[4]->push_back(m_out_trem[0]);
      m_word_mbd[5]->push_back(m_out_trem[1]);
      m_word_mbd[6]->push_back(m_out_vtx_sub);
      m_word_mbd[7]->push_back(m_out_vtx_add);

      if (m_out_nhit[0] >= m_nhit1)
      {
        bits.at(is) ^= 1U << 0U;
      }
      if (m_out_nhit[1] >= m_nhit1)
      {
        bits.at(is) ^= 1U << 1U;
      }
      if (m_out_nhit[0] >= m_nhit2)
      {
        bits.at(is) ^= 1U << 2U;
      }
      if (m_out_nhit[1] >= m_nhit2)
      {
        bits.at(is) ^= 1U << 3U;
      }

      if (m_out_nhit[0] >= m_nhit1 && m_out_nhit[1] >= m_nhit1 && m_out_vtx_sub <= m_timediff1)
      {
        bits.at(is) ^= 1U << 4U;
      }
      if (m_out_nhit[0] >= m_nhit1 && m_out_nhit[1] >= m_nhit1 && m_out_vtx_sub <= m_timediff2)
      {
        bits.at(is) ^= 1U << 5U;
      }
      if (m_out_nhit[0] >= m_nhit1 && m_out_nhit[1] >= m_nhit1 && m_out_vtx_sub <= m_timediff3)
      {
        bits.at(is) ^= 1U << 6U;
      }
      if (m_out_nhit[0] >= m_nhit2 && m_out_nhit[1] >= m_nhit2 && m_out_vtx_sub <= m_timediff1)
      {
        bits.at(is) ^= 1U << 7U;
      }
      if (m_out_nhit[0] >= m_nhit2 && m_out_nhit[1] >= m_nhit2 && m_out_vtx_sub <= m_timediff2)
      {
        bits.at(is) ^= 1U << 8U;
      }
      if (m_out_nhit[0] >= m_nhit2 && m_out_nhit[1] >= m_nhit2 && m_out_vtx_sub <= m_timediff3)
      {
        bits.at(is) ^= 1U << 9U;
      }

      if (Verbosity())
      {
        std::cout << "Trigger Word : " << std::bitset<16>(bits.at(is)) << std::dec << std::endl;
      }
    }

    for (int is = 0; is < nsample; is++)
    {
      trig_bits->push_back(bits.at(is));
    }
    if (Verbosity() >= 2)
    {
      std::cout << " " << std::endl;
    }

    for (int iw = 0; iw < 8; iw++)
    {
      std::vector<unsigned int> *sum = m_ll1out->get_word(iw);
      for (int is = 0; is < nsample; is++)
      {
        sum->push_back(m_word_mbd[iw]->at(is));
      }
    }
  }

  else
  {
    std::cout << "Trigger " << m_trigger << " not implemented" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void CaloTriggerEmulator::GetNodes(PHCompositeNode *topNode)
{
  if (Verbosity() >= 2)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  m_ll1out = findNode::getClass<LL1Out>(topNode, m_ll1_nodename);

  if (!m_ll1out)
  {
    std::cout << "No LL1Out found... " << std::endl;
    exit(1);
  }

  m_primitives = findNode::getClass<TriggerPrimitiveContainer>(topNode, m_prim_nodename);

  if (!m_primitives)
  {
    std::cout << "No TriggerPrimitives found... " << std::endl;
    exit(1);
  }

  bool hcalset = false;
  if (m_do_hcalout)
  {
    m_waveforms_hcalout = findNode::getClass<TowerInfoContainer>(topNode, "WAVEFORM_HCALOUT");

    if (!m_waveforms_hcalout)
    {
      std::cout << "No HCALOUT Waveforms found... " << std::endl;
      exit(1);
    }

    m_primitives_hcalout = findNode::getClass<TriggerPrimitiveContainer>(topNode, "TRIGGERPRIMITIVES_HCALOUT");

    if (!m_primitives_hcalout)
    {
      std::cout << "No HCAL Primitives found... " << std::endl;
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
    m_waveforms_hcalin = findNode::getClass<TowerInfoContainer>(topNode, "WAVEFORM_HCALIN");

    if (!m_waveforms_hcalin)
    {
      std::cout << "No HCAL Waveforms found... " << std::endl;
      exit(1);
    }

    m_primitives_hcalin = findNode::getClass<TriggerPrimitiveContainer>(topNode, "TRIGGERPRIMITIVES_HCALIN");

    if (!m_primitives_hcalin)
    {
      std::cout << "No HCAL Primitives found... " << std::endl;
      exit(1);
    }

    if (!hcalset)
    {
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
    m_waveforms_emcal = findNode::getClass<TowerInfoContainer>(topNode, "WAVEFORM_CEMC");

    if (!m_waveforms_emcal)
    {
      std::cout << "No HCAL Waveforms found... " << std::endl;
      exit(1);
    }

    m_primitives_emcal = findNode::getClass<TriggerPrimitiveContainer>(topNode, "TRIGGERPRIMITIVES_EMCAL");

    if (!m_primitives_emcal)
    {
      std::cout << "No HCAL Primitives found... " << std::endl;
      exit(1);
    }

    m_primitives_emcal_ll1 = findNode::getClass<TriggerPrimitiveContainer>(topNode, "TRIGGERPRIMITIVES_EMCAL_LL1");

    if (!m_primitives_emcal_ll1)
    {
      std::cout << "No HCAL Primitives found... " << std::endl;
      exit(1);
    }
  }
  if (m_do_mbd)
  {
    m_waveforms_mbd = findNode::getClass<TowerInfoContainer>(topNode, "WAVEFORM_TOWERS_MBD");

    if (!m_waveforms_mbd)
    {
      std::cout << "No HCAL Waveforms found... " << std::endl;
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

  LL1Out *ll1out = findNode::getClass<LL1Out>(ll1Node, m_ll1_nodename);
  if (!ll1out)
  {
    ll1out = new LL1Outv1(m_trigger, "NONE");
    PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(ll1out, m_ll1_nodename, "PHObject");
    ll1Node->addNode(LL1OutNode);
  }

  TriggerPrimitiveContainer *ll1out_prim = findNode::getClass<TriggerPrimitiveContainer>(ll1Node, m_prim_nodename);
  if (!ll1out_prim)
  {
    ll1out_prim = new TriggerPrimitiveContainerv1(m_triggerid, TriggerDefs::DetectorId::noneDId);
    PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(ll1out_prim, m_prim_nodename, "PHObject");
    ll1Node->addNode(LL1OutNode);
  }

  if (m_do_emcal)
  {
    std::string ll1_nodename = "TRIGGERPRIMITIVES_EMCAL";
    TriggerPrimitiveContainer *ll1out_d = findNode::getClass<TriggerPrimitiveContainer>(ll1Node, ll1_nodename);
    if (!ll1out_d)
    {
      ll1out_d = new TriggerPrimitiveContainerv1(TriggerDefs::TriggerId::noneTId, TriggerDefs::DetectorId::emcalDId);
      PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(ll1out_d, ll1_nodename, "PHObject");
      ll1Node->addNode(LL1OutNode);
    }

    ll1_nodename = "TRIGGERPRIMITIVES_EMCAL_LL1";
    TriggerPrimitiveContainer *ll1out_d1 = findNode::getClass<TriggerPrimitiveContainer>(ll1Node, ll1_nodename);
    if (!ll1out_d1)
    {
      ll1out_d1 = new TriggerPrimitiveContainerv1(m_triggerid, TriggerDefs::DetectorId::emcalDId);
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
      ll1out_d1 = new TriggerPrimitiveContainerv1(m_triggerid, TriggerDefs::DetectorId::hcalDId);
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
      ll1out_d1 = new TriggerPrimitiveContainerv1(m_triggerid, TriggerDefs::DetectorId::hcalDId);
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
  std::cout << "Total passed: " << m_npassed << "/" << m_nevent << std::endl;
  std::cout << "------------------------" << std::endl;

  return 0;
}

void CaloTriggerEmulator::identify()
{
  std::cout << " CaloTriggerEmulator: " << m_trigger << std::endl;
  std::cout << " LL1Out: " << std::endl;

  if (m_ll1out)
  {
    m_ll1out->identify();
  }

  if (m_primitives)
  {
    m_primitives->identify();
  }

  if (m_primitives_emcal_ll1)
  {
    m_primitives_emcal_ll1->identify();
  }

  if (m_primitives_hcal_ll1)
  {
    m_primitives_hcal_ll1->identify();
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

unsigned int CaloTriggerEmulator::getBits(unsigned int sum)
{
  unsigned int bit = 0;
  if (m_single_threshold)
  {
    bit |= (sum >= m_threshold ? 0x1U : 0);
    return bit;
  }
  for (unsigned int i = 0; i < 4; i++)
  {
    bit |= (sum >= m_threshold_calo[i] ? 0x1U << (i) : 0);
  }

  return bit;
}
