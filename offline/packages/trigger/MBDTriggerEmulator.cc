#include "MBDTriggerEmulator.h"

#include "LL1Defs.h"
#include "LL1Outv1.h"
#include "TriggerPrimitiveContainerv1.h"
#include "TriggerPrimitivev1.h"
#include <ffarawobjects/CaloPacketContainer.h>
#include <ffarawobjects/CaloPacket.h>
#include <ffarawobjects/CaloPacketContainerv1.h>
#include <ffarawobjects/CaloPacketv1.h>
#include <ffamodules/CDBInterface.h>
#include <cdbobjects/CDBHistos.h>  // for CDBHistos
#include <Event/Event.h>
#include <Event/EventTypes.h>
#include <Event/packet.h>

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

#include <bitset>
#include <cassert>
#include <cstdint>
#include <sstream>
#include <string>

// constructor
MBDTriggerEmulator::MBDTriggerEmulator(const std::string &name)
  : SubsysReco(name)
{
  // initialize all important counters

  m_nevent = 0;
  m_mbd_npassed = 0;
  m_isdata = false;
  // default nsamples is 16 for mbd, 12 for calos
  m_nsamples = 16;

  m_idx = 8;
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

  for (int i = 0; i < 2; i++)
  {
    m_out_tsum[i] = 0;
    m_out_nhit[i] = 0;
    m_out_tavg[i] = 0;
    m_out_trem[i] = 0;
  }
  m_out_vtx_sub = 0;
  m_out_vtx_add = 0;

  for (uint16_t i = 0; i < 128;i++)
  {
    uint16_t key = ( i & 0x7fffU );
    h_mbd_charge_lut[key] = nullptr;
    key = (0x1U << 0xfU) + ( i & 0x7fffU );
    h_mbd_time_lut[key] = nullptr;
    h_mbd_slewing_lut[key] = nullptr;
  }

  m_ll1out_mbd = nullptr;

  m_n_primitives = 4;
  m_n_sums = 13;
  m_trig_sample = -1;
  m_trig_sub_delay = 4;
  m_threshold = 1;

  m_nhit1 = 1;
  m_nhit2 = 2;
  m_timediff1 = 10;
  m_timediff2 = 20;
  m_timediff3 = 30;

}

// make file and histomanager (but nothing goes in the file at the moment)
int MBDTriggerEmulator::Init(PHCompositeNode * /*topNode*/)
{
  return 0;
}

// at the beginning of the run
int MBDTriggerEmulator::InitRun(PHCompositeNode *topNode)
{
  // Get the detectors that are used for a given trigger.

  m_ll1_nodename = "LL1OUT_MBD";
  m_prim_nodename = "TRIGGERPRIMITIVES_MBD";

  // Get the calibrations and produce the lookup tables;

  if (Download_Calibrations())
  {
    std::cout << "here" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  CreateNodes(topNode);

  return 0;
}

int MBDTriggerEmulator::Download_Calibrations()
{

  if (!m_mbd_charge_lutname.empty())
    {
      cdbttree_mbd_charge = new CDBHistos(m_mbd_charge_lutname);
    }
  else
    {
      std::string calibdir = CDBInterface::instance()->getUrl("mbd_charge_trigger_lut");
      if (calibdir.empty())
	{
	  m_default_lut_mbd = true;
	  std::cout << "Could not find and load histograms for MBD LUTs! defaulting to the identity table!" << std::endl;
	}
      else
	{
	  cdbttree_mbd_charge = new CDBHistos(calibdir);
	}
    }
  if (cdbttree_mbd_charge)
    {
      cdbttree_mbd_charge->LoadCalibrations();
      
      for (uint16_t i = 0; i < 128;i++)
	{
	  std::string histoname = "h_mbd_charge_lut_" + std::to_string(i);
	  uint16_t key = (i & 0x7fffU);
	  h_mbd_charge_lut[key] = (TH1I*) cdbttree_mbd_charge->getHisto(histoname.c_str());
	}
    }
  else
    {
      return Fun4AllReturnCodes::ABORTRUN;
    }  
  // now time channels
  if (!m_mbd_time_lutname.empty())
    {
      cdbttree_mbd_time = new CDBHistos(m_mbd_time_lutname);
    }
  else
    {
      std::string calibdir = CDBInterface::instance()->getUrl("mbd_time_trigger_lut");
      if (calibdir.empty())
	{
	  m_default_lut_mbd = true;
	  std::cout << "Could not find and load histograms for MBD LUTs! defaulting to the identity table!" << std::endl;
	}
      else
	{
	  cdbttree_mbd_time = new CDBHistos(calibdir);
	}
    }
  if (cdbttree_mbd_time)
    {
      cdbttree_mbd_time->LoadCalibrations();
      
      for (uint16_t i = 0; i < 128;i++)
	{
	  std::string histoname = "h_mbd_time_lut_" + std::to_string(i);
	  uint16_t key = (i&0x7fffU) + (0x1U << 0xfU);
	  h_mbd_time_lut[key] = (TH1I*) cdbttree_mbd_time->getHisto(histoname.c_str());
	}
    }
  else
    {
      return Fun4AllReturnCodes::ABORTRUN;
    }  

  if (!m_mbd_slewing_lutname.empty())
    {
      cdbttree_mbd_slewing = new CDBHistos(m_mbd_slewing_lutname);
    }
  else
    {
      std::string calibdir = CDBInterface::instance()->getUrl("mbd_slewing_trigger_lut");
      if (calibdir.empty())
	{
	  m_default_lut_mbd = true;
	  std::cout << "Could not find and load histograms for MBD LUTs! defaulting to the identity table!" << std::endl;
	}
      else
	{
	  cdbttree_mbd_slewing = new CDBHistos(calibdir);
	}
    }
  if (cdbttree_mbd_slewing)
    {
      cdbttree_mbd_slewing->LoadCalibrations();
      
      for (uint16_t i = 0; i < 128;i++)
	{
	  std::string histoname = "h_mbd_slewing_lut_" + std::to_string(i);
	  uint16_t key = (i&0x7fffU) + (0x1U << 0xfU);
	  h_mbd_slewing_lut[key] = (TH1I*) cdbttree_mbd_slewing->getHisto(histoname.c_str());
	}
    }
  else
    {
      return Fun4AllReturnCodes::ABORTRUN;
    }  
  return 0;
}
// process event procedure
int MBDTriggerEmulator::process_event(PHCompositeNode *topNode)
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
int MBDTriggerEmulator::ResetEvent(PHCompositeNode * /*topNode*/)
{
  m_peak_sub_ped_mbd.clear();
  return 0;
}

int MBDTriggerEmulator::process_waveforms()
{
  // Get range of waveforms
  if (Verbosity())
  {
    std::cout << __FILE__ << "::" << __FUNCTION__ << ":: Processing waveforms" << std::endl;
  }

  if (!m_waveforms_mbd)
    {
      return process_raw();
    }
  if (Verbosity())
    {
      std::cout << __LINE__ << std::endl;
    }
  int sample_start = 1;
  int sample_end = m_nsamples;
  if (m_trig_sample > 0)
  {
    sample_start = m_trig_sample;
    sample_end = m_trig_sample + 1;
  }
  if (Verbosity())
    {
      std::cout << __LINE__ << std::endl;
    }

  // for each waveform, clauclate the peak - pedestal given the sub-delay setting
  CaloPacket *dstp[2]{nullptr};

  for (int ipkt = 0; ipkt < 2; ipkt++)
    {
      if (Verbosity())
	{
	  std::cout << __LINE__ << std::endl;
	}

      int pktid = 1001 + ipkt;  // packet id

      dstp[ipkt] = m_waveforms_mbd->getPacketbyId(pktid);

      if (dstp[ipkt])
	{
	  for (int ich = 0; ich < 128; ich++)
	    {

	      //	      int feech = ipkt * 128 + ich;
	      uint16_t isTime = (ich%16 < 8 ? 1 : 0);
	      
	      uint16_t ipmt = ipkt * 64 + (ich/16) * 8 + (ich%8);
	      uint16_t key = ((isTime & 0x1U) << 0xfU) + (ipmt & 0x7fffU);

	      std::vector<unsigned int> v_peak_sub_ped;
	      unsigned int peak_sub_ped = 0;

	      for (int i = sample_start; i < sample_end; i++)
		{
		  int16_t maxim = dstp[ipkt]->iValue(i, ich);// > dstp[ipkt]->iValue(i+1, ich) ? dstp[ipkt]->iValue(i, ich) : dstp[ipkt]->iValue(i+1, ich));
		  //		  maxim = (maxim > dstp[ipkt]->iValue(i+2, ich) ? maxim : dstp[ipkt]->iValue(i+2, ich));

		  if (Verbosity())
		    {
		      std::cout << __LINE__ << std::endl;
		    }

		  int subtraction = maxim - dstp[ipkt]->iValue(((i - m_trig_sub_delay > 0 ? i - m_trig_sub_delay : 0)), ich);

		  if (subtraction < 0)
		    {
		      subtraction = 0;
		    }
		    
		  peak_sub_ped = (((unsigned int) subtraction) & 0x3fffU);
		  v_peak_sub_ped.push_back(peak_sub_ped);
		}
	      // save in global.
	      m_peak_sub_ped_mbd[key] = v_peak_sub_ped;
	    }
	  if (Verbosity())
	    {
	      std::cout << __LINE__ << std::endl;
	    }
	  delete dstp[ipkt];
	}

      if (Verbosity())
	{
	  std::cout << __LINE__ << std::endl;
	}

    }
  if (Verbosity())
    {
      std::cout << __LINE__ << std::endl;
    }

  return Fun4AllReturnCodes::EVENT_OK;
}
int MBDTriggerEmulator::process_raw()
{
  // Get range of waveforms
  if (Verbosity())
  {
    std::cout << __FILE__ << "::" << __FUNCTION__ << ":: Processing waveforms" << std::endl;
  }

  if (m_useoffline)
    {
      return Fun4AllReturnCodes::ABORTEVENT;
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

  int sample_start = 1;
  int sample_end = m_nsamples;
  if (m_trig_sample > 0)
  {
    sample_start = m_trig_sample;
    sample_end = m_trig_sample + 1;
  }

  // for each waveform, clauclate the peak - pedestal given the sub-delay setting
  Packet *dstp[2]{nullptr};

  for (int ipkt = 0; ipkt < 2; ipkt++)
    {
      int pktid = 1001 + ipkt;  // packet id

      dstp[ipkt] = m_event->getPacket(pktid);

      if (dstp[ipkt])
	{
	  for (int ich = 0; ich < 128; ich++)
	    {

	      //	      int feech = ipkt * 128 + ich;
	      uint16_t isTime = (ich%16 < 8 ? 1 : 0);
	      
	      uint16_t ipmt = ipkt * 64 + (ich/16) * 8 + (ich%8);
	      uint16_t key = ((isTime & 0x1U) << 0xfU) + (static_cast<uint16_t>(ipmt) & 0x7fffU);

	      std::vector<unsigned int> v_peak_sub_ped;
	      unsigned int peak_sub_ped = 0;

	      for (int i = sample_start; i < sample_end; i++)
		{
		  int16_t maxim = dstp[ipkt]->iValue(i, ich);// > dstp[ipkt]->iValue(i+1, ich) ? dstp[ipkt]->iValue(i, ich) : dstp[ipkt]->iValue(i+1, ich));
		  //		  maxim = (maxim > dstp[ipkt]->iValue(i+2, ich) ? maxim : dstp[ipkt]->iValue(i+2, ich));


		  int subtraction = maxim - dstp[ipkt]->iValue(((i - m_trig_sub_delay > 0 ? i - m_trig_sub_delay : 0)), ich);

		  if (subtraction < 0)
		    {
		      subtraction = 0;
		    }
		    
		  peak_sub_ped = (((unsigned int) subtraction) & 0x3fffU);
		  v_peak_sub_ped.push_back(peak_sub_ped);
		}
	      // save in global.
	      m_peak_sub_ped_mbd[key] = v_peak_sub_ped;
	    }
	}
      delete dstp[ipkt];
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

// procedure to process the peak - pedestal into primitives.
int MBDTriggerEmulator::process_primitives()
{

  int i = 0;
  int nsample = m_nsamples - 1;
  if (m_trig_sample > 0)
  {
    nsample = 1;
  }

  if (Verbosity())
  {
    std::cout << __FILE__ << "::" << __FUNCTION__ << ":: Processing primitives" << std::endl;
  }
  
  for (i = 0; i < m_n_primitives; i++)
    {
      // make primitive key
      TriggerDefs::TriggerPrimKey primkey = TriggerDefs::getTriggerPrimKey(m_triggerid, TriggerDefs::GetDetectorId("MBD"), TriggerDefs::GetPrimitiveId("MBD"), i);

      // make primitive and check mask;
      TriggerPrimitive *primitive = m_primitives_mbd->get_primitive_at_key(primkey);

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
	// go through 8 charge channels
	for (int j = 0; j < 32; j++)
          {
            // pass upper 10 bits of charge to get 10 bit LUt outcome
            // tmp = m_l1_adc_table[m_peak_sub_ped_mbd[i * 64 + 8 + isec * 16 + j].at(is) >> 4U];
	    uint16_t key = static_cast<uint16_t>(i*32 + j) & 0x7fffU;
	    unsigned int peak_minus_ped = m_peak_sub_ped_mbd[key].at(is) >> 0x4U;
	    tmp = h_mbd_charge_lut[key]->GetBinContent(peak_minus_ped + 1);

            // put upper 3 bits of the 10 bits into slewing correction later
            qadd[j] = (tmp & 0x380U) >> 7U;

            // sum up to 11 bits.
            m_trig_charge[j / 4] += tmp & 0x7fU;
          }
        
	// 8 timing channels
	for (int j = 0; j < 32; j++)
          {

            // pass upper 10 bits of charge to get 10 bit LUt outcome
            // tmp = m_l1_adc_table[m_peak_sub_ped_mbd[i * 64 + 8 + isec * 16 + j].at(is) >> 4U];
	    uint16_t key = (static_cast<uint16_t>(i*32 + j) & 0x7fffU) + (0x1U << 0xfU);
	    unsigned int peak_minus_ped = m_peak_sub_ped_mbd[key].at(is) >> 0x4U;
	    tmp = h_mbd_time_lut[key]->GetBinContent(peak_minus_ped + 1);

            m_trig_nhit += (tmp & 0x200U) >> 9U;
	    unsigned int slewkey = (qadd[j] << 9U) + (tmp & 0x01ffU);
            // get upper 3 bits of charge in the channel, and make it bits 9-11, the time of the chanel is the lower 9 bits from 0-8.
	    tmp2 = h_mbd_slewing_lut[key]->GetBinContent(slewkey + 1);
            // attribute to the time sum
            m_trig_time[j/8] += tmp2;
          }
        // ad in the charge sums

        for (int j = 0; j < 13; j++)
        {
	  TriggerDefs::TriggerSumKey sumkey = TriggerDefs::getTriggerSumKey(m_triggerid, TriggerDefs::GetDetectorId("MBD"), TriggerDefs::GetPrimitiveId("MBD"), i, j);
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
	      primitive->get_sum_at_key(sumkey)->push_back(m_trig_time[j - 9 ]);
	    }

        }
      }
      
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

// Unless this is the MBD or HCAL Cosmics trigger, EMCAL and HCAL will go through here.
// This creates the 8x8 non-overlapping sum and the 4x4 overlapping sum.

// This is where the LL1 Jet/Pair/Cosmic algorithm is

int MBDTriggerEmulator::process_trigger()
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




  std::vector<unsigned int> *trig_bits = m_ll1out_mbd->GetTriggerBits();
  if (!m_primitives_mbd)
    {
      std::cout << "There is no primitive container" << std::endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }

  TriggerPrimitivev1::Range sumrange;
  int ip, isum;

  TriggerPrimitiveContainerv1::Range range = m_primitives_mbd->getTriggerPrimitives();

  if (Verbosity() >= 2)
    {
      std::cout << __FUNCTION__ << " " << __LINE__ << " mbd primitives size: " << m_primitives_mbd->size() << std::endl;
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
  int pass = 0;
  for (int is = 0; is < nsample; is++)
    {
      trig_bits->push_back(bits.at(is));
      if (bits.at(is)) 
	{
	  pass = 1;
	}
    }
  if (pass) 
    {
      m_mbd_npassed++;
    }
  if (Verbosity() >= 2)
    {
      std::cout << " " << std::endl;
    }

  for (int iw = 0; iw < 8; iw++)
    {
      std::vector<unsigned int>* sum = m_ll1out_mbd->get_word(iw);
      for (int is = 0; is < nsample; is++)
	{
	  sum->push_back(m_word_mbd[iw]->at(is));
	}
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

void MBDTriggerEmulator::GetNodes(PHCompositeNode *topNode)
{
  if (Verbosity() >= 2)
  {
    std::cout << __FUNCTION__ << std::endl;
  }
  
m_event = findNode::getClass<Event>(topNode, "PRDF");
  
  m_ll1_nodename="LL1OUT_MBD";
  
  m_ll1out_mbd = findNode::getClass<LL1Out>(topNode, m_ll1_nodename);
  
  if (!m_ll1out_mbd)
    {
      std::cout << "No LL1Out found... " << std::endl;
      exit(1);
    }

  m_prim_nodename = "TRIGGERPRIMITIVES_MBD";
  m_primitives_mbd = findNode::getClass<TriggerPrimitiveContainer>(topNode, m_prim_nodename);

  if (!m_primitives_mbd)
  {
    std::cout << "No TriggerPrimitives found... " << std::endl;
    exit(1);
  }

  m_waveforms_mbd = findNode::getClass<CaloPacketContainer>(topNode, "MBDPackets");
  if (m_waveforms_mbd)
    {
      m_useoffline = true;
    }
  // if (!m_waveforms_mbd)
  //   {
  //     std::cout << "No MBD Waveforms found... " << std::endl;
  //     exit(1);
  //   }

  return;
}
void MBDTriggerEmulator::CreateNodes(PHCompositeNode *topNode)
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

  m_ll1_nodename = "LL1OUT_MBD";
  LL1Out *ll1out = findNode::getClass<LL1Out>(ll1Node, m_ll1_nodename);
  if (!ll1out)
    {
      ll1out = new LL1Outv1("MBD", "NONE");
      PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(ll1out, m_ll1_nodename, "PHObject");
      ll1Node->addNode(LL1OutNode);
    }
  m_prim_nodename = "TRIGGERPRIMITIVES_MBD";
  TriggerPrimitiveContainer *ll1out_prim = findNode::getClass<TriggerPrimitiveContainer>(ll1Node, m_prim_nodename);
  if (!ll1out_prim)
    {
      ll1out_prim = new TriggerPrimitiveContainerv1(TriggerDefs::TriggerId::mbdTId, TriggerDefs::DetectorId::mbdDId);
      PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(ll1out_prim, m_prim_nodename, "PHObject");
      ll1Node->addNode(LL1OutNode);
    }
  return;
}

int MBDTriggerEmulator::End(PHCompositeNode * /*topNode*/)
{
  delete cdbttree_mbd_charge;
  delete cdbttree_mbd_time;
  delete cdbttree_mbd_slewing;

  std::cout << "------------------------" << std::endl;
  std::cout << "Total MBD passed: " << m_mbd_npassed << "/" << m_nevent << std::endl;
  std::cout << "------------------------" << std::endl;

  return 0;
}

void MBDTriggerEmulator::identify()
{
  std::cout << " LL1Out: " << std::endl;
  m_ll1out_mbd->identify();
  std::cout << " TriggerPrimitives: " << std::endl;
  m_primitives_mbd->identify();
  std::cout << " Processed " << m_nevent << std::endl;
}
