#include "CaloTriggerEmulator.h"

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>
#include <calobase/TowerInfoContainerv3.h>
#include <phool/getClass.h>

#include <cdbobjects/CDBTTree.h>  // for CDBTTree
#include <ffamodules/CDBInterface.h>
#include <TFile.h>
#include <cassert>
#include <sstream>
#include <string>
#include <cstdint>
#include <bitset>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include "LL1Defs.h"
#include "LL1Outv1.h"
#include "TriggerPrimitivev1.h"
#include "TriggerPrimitiveContainerv1.h"

using namespace std;


// constructor
CaloTriggerEmulator::CaloTriggerEmulator(const std::string& name, const std::string& filename)
  : SubsysReco(name)
  , outfilename(filename)
  , hm(nullptr)
  , outfile(nullptr)
  , _trigger("NONE")
{
  // initialize all important parameters
 
  _nevent = 0;
  _npassed = 0;
  m_emcal_GeV_per_lut = 0.1;
  m_emcal_lut_offset = 0.0;
  m_emcal_lut_floor = 0.0;
  m_isdata = false;
  // default nsamples is 31;
  m_nsamples = 31;

  // for MBD, this is the peak sample in run-23 data
  _idx = 12;

  // default values for the lookup tables.
  // TODO: to CDB the LUTs from the database

  for (unsigned int i = 0; i < 1024; i++)
    {

      if (i < 1024/4)
	{
	  m_l1_adc_table[i] = (4*i) & 0x3ff;
	}
      else
	{
	  m_l1_adc_table[i] = 0x3ff;
	}
      m_l1_adc_table_time[i] = (i) & 0x3ff;
    }
  for (unsigned int i = 0; i < 4096; i++)
    {
      m_l1_slewing_table[i] = (i) & 0x1ff;
    }

  // reset variables
  for (int j = 0; j < 8; j++)
    {
      m_trig_charge[j] = 0;
      for (int k = 0; k < 4; k++)
	{
	  m2_trig_charge[k][j] = 0;
	}
    }
  m_trig_nhit = 0;
  for (int k = 0; k < 4; k++)
    {
      m2_trig_nhit[k] = 0;
    }
  for (int j = 0; j < 4; j++)
    {
      m_trig_time[j] = 0;
      for (int k = 0; k < 4; k++)
	{
	  m2_trig_time[j][k] = 0;
	}
    }
  for (int i = 0; i < 2; i ++)
    {
      m_out_tsum[i] = 0;
      m_out_nhit[i] = 0;
      m_out_tavg[i] = 0;
      m_out_trem[i] = 0;
    }
  m_out_vtx_sub = 0;
  m_out_vtx_add = 0;

  // Set HCAL LL1 lookup table for the cosmic coincidence trigger.
  if (_triggerid == TriggerDefs::TriggerId::cosmic_coinTId)
    {
      unsigned int bits1, bits2, sumbits1, sumbits2;
      for (int i = 0; i < 4096; i++)
	{
	  sumbits1 = 0;
	  sumbits2 = 0;

	  bits1 = (i & 0x3f);
	  bits2 = ((i>>6) & 0x3f);
	  for (int j = 0; j < 3; j++)
	    {
	      if (((bits1>>j)&0x1) && ((bits2>>j)&0x1)) sumbits1 ++; 
              if (((bits1>>(j+3))&0x1) && ((bits2>>(j+3))&0x1)) sumbits2 ++;
	    }

	  m_l1_hcal_table[i] = 0;
	  if (i == 0) continue;
	  m_l1_hcal_table[i] |= (sumbits1 ? 0x1 : 0);
	  m_l1_hcal_table[i] |= (sumbits2 ? 0x2 : 0);
	}
    }

  else if (_triggerid == TriggerDefs::TriggerId::cosmicTId)
    {
      unsigned int bits1, bits2, sumbits1, sumbits2;
      for (int i = 0; i < 4096; i++)
	{
	  sumbits1 = 0;
	  sumbits2 = 0;

	  bits1 = (i & 0x3f);
	  bits2 = ((i>>6) & 0x3f);
	  for (int j = 0; j < 6; j++)
	    {
	      sumbits1 += ((bits1>>j) & 0x1);
	      sumbits2 += ((bits2>>j) & 0x1);
	    }

	  m_l1_hcal_table[i] = 0;
	  if (i == 0) continue;
	  m_l1_hcal_table[i] |= (sumbits1 ? 0x1 : 0);
	  m_l1_hcal_table[i] |= (sumbits2 ? 0x2 : 0);
	}
    }


  // point to null for all of these objects to be added to or grabbed from the node tree.
  // this will hold the ll1 info that goes through the emulator
  _ll1out = 0;

  // waveform containers to be grabbed from node tree.
  // Done int the CaloPacketGetter

  _waveforms_emcal = 0;
  _waveforms_hcalin = 0;
  _waveforms_hcalout = 0;
  _waveforms_mbd = 0;

  // to hold the primitives constructed from the waveforms.
  _primitives = 0;
  _primitives_emcal = 0;
  _primitives_hcalin = 0;
  _primitives_hcalout = 0;

  _primitive = 0;
  _n_primitives = 0;
  _n_sums = 16;
  _m_trig_sub_delay = 4;
  _m_threshold = 1;

  m_nhit1 = 2;
  m_nhit2 = 10;
  m_timediff1 = 10;
  m_timediff2 = 20;
  m_timediff3 = 30;


  // define a detector map for detectors included in a trigger
  _m_det_map[TriggerDefs::TriggerId::noneTId] = {}; 
  _m_det_map[TriggerDefs::TriggerId::jetTId] = {"EMCAL"};
  _m_det_map[TriggerDefs::TriggerId::mbdTId] = {"MBD"};
  _m_det_map[TriggerDefs::TriggerId::cosmicTId] = {"HCALOUT"};
  _m_det_map[TriggerDefs::TriggerId::cosmic_coinTId] = {"HCALOUT"};
  _m_det_map[TriggerDefs::TriggerId::pairTId] = {"EMCAL"};
  _m_det_map[TriggerDefs::TriggerId::photonTId] = {"EMCAL"};

  // define primitive map as number of primitives in a detector.
  _m_prim_map[TriggerDefs::DetectorId::noneDId] = 0;
  _m_prim_map[TriggerDefs::DetectorId::emcalDId] = 384;
  _m_prim_map[TriggerDefs::DetectorId::hcalinDId] = 24;
  _m_prim_map[TriggerDefs::DetectorId::hcaloutDId] = 24;
  _m_prim_map[TriggerDefs::DetectorId::mbdDId] = 4;

  _m_prim_ll1_map[TriggerDefs::TriggerId::jetTId] = 16;
  _m_prim_ll1_map[TriggerDefs::TriggerId::pairTId] = 16;
  _m_prim_ll1_map[TriggerDefs::TriggerId::photonTId] = 384;
  // booleans to control the input of detector data
  _do_emcal = false;
  _do_hcalin = false;
  _do_hcalout = false;
  _do_mbd = false;

  
  _masks_channel = {};//, 70385703};
  _masks_fiber = {};//, 70385696};
}


// destructr
CaloTriggerEmulator::~CaloTriggerEmulator()
{
  delete hm;
}


// check whether a channel has been masked
bool CaloTriggerEmulator::CheckChannelMasks(TriggerDefs::TriggerSumKey key)
{
  for (auto it = _masks_channel.begin(); it != _masks_channel.end(); ++it)
    {
      if (key == (*it)) return true;
    }
  return false;
}

// cehck whether a fiber has been masked
bool CaloTriggerEmulator::CheckFiberMasks(TriggerDefs::TriggerPrimKey key)
{
  for (auto it = _masks_fiber.begin(); it != _masks_fiber.end(); ++it)
    {
      if (key == (*it)) return true;
    }
  return false;
}


// setting the trigger type
void CaloTriggerEmulator::setTriggerType(const std::string &name)
{
  _trigger = name;
  _triggerid = TriggerDefs::GetTriggerId(_trigger);
  std::cout << "Setting Trigger type: "<<_trigger<< " (" <<_triggerid<<")"<<std::endl;  
}

// make file and histomanager (but nothing goes in the file at the moment)
int CaloTriggerEmulator::Init(PHCompositeNode* /*topNode*/)
{
  hm = new Fun4AllHistoManager(Name());
  outfile = new TFile(outfilename.c_str(), "RECREATE");

  if(Verbosity()>=2) std::cout << __FUNCTION__ << std::endl;

  return 0;
}


// at the beginning of the run
int CaloTriggerEmulator::InitRun(PHCompositeNode* topNode)
{
  if(Verbosity()>=2) std::cout << __FUNCTION__ << std::endl;

  // Get the detectors that are used for a given trigger.

  if (_triggerid == TriggerDefs::TriggerId::jetTId)
    {
      if(Verbosity()>=2) std::cout << "Using Jet Trigger."<<std::endl;
      _do_emcal = true;
      _do_hcalin = false;
      _do_hcalout = false;
      _do_mbd = false;
    }
  else if (_triggerid == TriggerDefs::TriggerId::photonTId)
    {
      if(Verbosity()>=2) std::cout << "Using Photon Trigger."<<std::endl;
      _do_emcal = true;
      _do_hcalin = false;
      _do_hcalout = false;
      _do_mbd = false;
    }
  else if (_triggerid == TriggerDefs::TriggerId::mbdTId)
    {
      if(Verbosity()>=2) std::cout << "Using MBD Trigger."<<std::endl;
      _do_emcal = false;
      _do_hcalin = false;
      _do_hcalout = false;      
      _do_mbd = true;
    }
  else if (_triggerid == TriggerDefs::TriggerId::cosmicTId)
    {
      if(Verbosity()>=2) std::cout << "Using Cosmic Trigger."<<std::endl;
      _do_emcal = false;
      _do_hcalin = false; //true;
      _do_hcalout = true;
      _do_mbd = false;
    }
  else if (_triggerid == TriggerDefs::TriggerId::cosmic_coinTId)
    {
      if(Verbosity()>=2) std::cout << "Using Cosmic Coincidence Trigger."<<std::endl;
      _do_emcal = false;
      _do_hcalin = false;//true;
      _do_hcalout = true;
      _do_mbd = false;
    }
  else if (_triggerid == TriggerDefs::TriggerId::pairTId)
    {
      if(Verbosity()>=2) std::cout << "Using Pair Trigger."<<std::endl;
      _do_emcal = true;
      _do_hcalin = false;
      _do_hcalout = false;
      _do_mbd = false;
    }
  else
    {
      cout << __FUNCTION__ << " : No trigger selected "<<endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

  // Set HCAL LL1 lookup table for the cosmic coincidence trigger.
  if (_triggerid == TriggerDefs::TriggerId::cosmic_coinTId)
    {
      unsigned int bits1, bits2, sumbits1, sumbits2;
      for (int i = 0; i < 4096; i++)
	{
	  sumbits1 = 0;
	  sumbits2 = 0;

	  bits1 = (i & 0x3f);
	  bits2 = ((i>>6) & 0x3f);
	  for (int j = 0; j < 3; j++)
	    {
	      if (((bits1>>j)&0x1) && ((bits2>>j)&0x1)) sumbits1 ++; 
              if (((bits1>>(j+3))&0x1) && ((bits2>>(j+3))&0x1)) sumbits2 ++;
	    }

	  m_l1_hcal_table[i] = 0;
	  if (i == 0) continue;
	  m_l1_hcal_table[i] |= (sumbits1 ? 0x1 : 0);
	  m_l1_hcal_table[i] |= (sumbits2 ? 0x2 : 0);
	}
    }

  // Set the HCAL LL1 lookup table as the singles trigger

  else if (_triggerid == TriggerDefs::TriggerId::cosmicTId)
    {
      unsigned int bits1, bits2, sumbits1, sumbits2;
      for (int i = 0; i < 4096; i++)
	{
	  sumbits1 = 0;
	  sumbits2 = 0;

	  bits1 = (i & 0x3f);
	  bits2 = ((i>>6) & 0x3f);
	  for (int j = 0; j < 6; j++)
	    {
	      sumbits1 += ((bits1>>j) & 0x1);
	      sumbits2 += ((bits2>>j) & 0x1);
	    }

	  m_l1_hcal_table[i] = 0;
	  if (i == 0) continue;
	  m_l1_hcal_table[i] |= (sumbits1 ? 0x1 : 0);
	  m_l1_hcal_table[i] |= (sumbits2 ? 0x2 : 0);
	}
    }

  _ll1_nodename = "LL1OUT_" + _trigger;
  _prim_nodename = "TRIGGERPRIMITIVES_" + _trigger;

  // for each detector in the detector map that is in the trigger.
  for ( auto iter = _m_det_map[_triggerid].begin() ; iter != _m_det_map[_triggerid].end() ; ++iter)
    { 
      // Get the number of primitives for this detector to calculate.
      for (int i = 0; i < _m_prim_map[TriggerDefs::GetDetectorId(*iter)]; i++)
	{

	  // if MBD make the correct number of histograms.
	  if (strcmp((*iter).c_str(), "MBD") == 0)
	    {

	      h_nhit = new TH1D(Form("h_nhit_%d",i),"", 33, -0.5, 32.5);
	      h2_line_up = new TH2D(Form("h_line_up_%d",i),"", m_nsamples, -0.5, static_cast<float>(m_nsamples) - 0.5, 60, -0.5, 59.5);

	      hm->registerHisto(h_nhit);
	      hm->registerHisto(h2_line_up);

	      v_nhit.push_back(h_nhit);
	      v_line_up.push_back(h2_line_up);
	      for (int j = 0; j < 8; j++)
		{

		  h_mbd_charge = new TH1D(Form("h_mbd_charge_%d_%d", j, i), "", 257, 0, 257);
		  v_mbd_charge[Form("h_mbd_charge_%d_%d", j, i)] = h_mbd_charge;
		  hm->registerHisto(h_mbd_charge);
		}	      
	      for (int j = 0; j < 4; j++)
		{
		  h_mbd_time = new TH1D(Form("h_mbd_time_%d_%d", j, i), "", 256, 0, 4096);
		  v_mbd_time[Form("h_mbd_time_%d_%d", j, i)] = h_mbd_time;
		  hm->registerHisto(h_mbd_time);
		}

	      continue;
	    }

	  // make overall histograms
	  peak_primitive = new TH2D(Form("peak_primitive_%s_%d", (*iter).c_str(), i), ";primitive;peak;counts", 16, -0.5, 15.5, m_nsamples - _m_trig_sub_delay, -0.5, m_nsamples - 1 - _m_trig_sub_delay);
	  avg_primitive = new TProfile(Form("avg_primitive_%s_%d",(*iter).c_str(),  i), ";primitive;avg", 16, -0.5, 15.5);
	  primitives = new TH2D(Form("primitives_%s_%d",(*iter).c_str(),  i), ";primitives;", 16, -0.5, 15.5, 64, 0, 256);
	  trigger_fire_map = new TH2D(Form("trigger_fire_map_%s_%d",(*iter).c_str(),  i), ";ch;ch", 4, -0.5, 3.5, 4, -0.5, 3.5);


	  // Make EMCAL and HCAL specific histograms.
	  if (strcmp((*iter).c_str(), "EMCAL") ==0)
	    {
	      v_peak_primitive_emcal.push_back(peak_primitive);
	      v_avg_primitive_emcal.push_back(avg_primitive);
	      v_primitives_emcal.push_back(primitives);
	      v_trigger_fire_map_emcal.push_back(trigger_fire_map);
	    }
	  if (strcmp((*iter).c_str(), "HCALIN") ==0)
	    {
	      v_peak_primitive_hcalin.push_back(peak_primitive);
	      v_avg_primitive_hcalin.push_back(avg_primitive);
	      v_primitives_hcalin.push_back(primitives);
	      v_trigger_fire_map_hcalin.push_back(trigger_fire_map);
	    }
	  if (strcmp((*iter).c_str(), "HCALOUT") ==0)
	    {
	      v_peak_primitive_hcalout.push_back(peak_primitive);
	      v_avg_primitive_hcalout.push_back(avg_primitive);
	      v_primitives_hcalout.push_back(primitives);
	      v_trigger_fire_map_hcalout.push_back(trigger_fire_map);
	    }
	  hm->registerHisto(peak_primitive);
	  hm->registerHisto(avg_primitive);
	  hm->registerHisto(primitives);
	  hm->registerHisto(trigger_fire_map);
	}
    }

  // now jet trigger 
  if (_triggerid == TriggerDefs::TriggerId::jetTId)
    {
      jet_trigger_fire_map = new TH2D(Form("trigger_fire_map_%s_all",_trigger.c_str()), ";ch;ch", 12, -0.5, 11.5, 32, -0.5, 31.5);
      hm->registerHisto(jet_trigger_fire_map);
      for (int i = 0; i < _m_prim_ll1_map[_triggerid]; i++)
	{
	  // make overall histograms
	  peak_primitive = new TH2D(Form("peak_primitive_%s_%d", _trigger.c_str(), i), ";primitive;peak;counts", 24, -0.5, 23.5, m_nsamples - _m_trig_sub_delay, -0.5, m_nsamples - 1 - _m_trig_sub_delay);
	  avg_primitive = new TProfile(Form("avg_primitive_%s_%d",_trigger.c_str(),  i), ";primitive;avg", 24, -0.5, 23.5);
	  primitives = new TH2D(Form("primitives_%s_%d",_trigger.c_str(),  i), ";primitives;", 24, -0.5, 23.5, 64, 0, 256);
	  trigger_fire_map = new TH2D(Form("trigger_fire_map_%s_%d",_trigger.c_str(),  i), ";ch;ch", 12, -0.5, 11.5, 2, -0.5, 1.5);



	  v_peak_primitive_jet.push_back(peak_primitive);
	  v_avg_primitive_jet.push_back(avg_primitive);
	  v_primitives_jet.push_back(primitives);
	  v_trigger_fire_map_jet.push_back(trigger_fire_map);

	  hm->registerHisto(peak_primitive);
	  hm->registerHisto(avg_primitive);
	  hm->registerHisto(primitives);
	  hm->registerHisto(trigger_fire_map);
	}
    }
  if (_triggerid == TriggerDefs::TriggerId::photonTId)
    {
      for (int i = 0; i < _m_prim_ll1_map[_triggerid]; i++)
	{
	  // make overall histograms
	  peak_primitive = new TH2D(Form("peak_primitive_%s_%d", _trigger.c_str(), i), ";primitive;peak;counts", 16, -0.5, 15.5, m_nsamples - _m_trig_sub_delay, -0.5, m_nsamples - 1 - _m_trig_sub_delay);
	  avg_primitive = new TProfile(Form("avg_primitive_%s_%d",_trigger.c_str(),  i), ";primitive;avg", 16, -0.5, 15.5);
	  primitives = new TH2D(Form("primitives_%s_%d",_trigger.c_str(),  i), ";primitives;", 16, -0.5, 15.5, 64, 0, 256);
	  trigger_fire_map = new TH2D(Form("trigger_fire_map_%s_%d",_trigger.c_str(),  i), ";ch;ch", 4, -0.5, 3.5, 4, -0.5, 3.5);

	  v_peak_primitive_photon.push_back(peak_primitive);
	  v_avg_primitive_photon.push_back(avg_primitive);
	  v_primitives_photon.push_back(primitives);
	  v_trigger_fire_map_photon.push_back(trigger_fire_map);

	  hm->registerHisto(peak_primitive);
	  hm->registerHisto(avg_primitive);
	  hm->registerHisto(primitives);
	  hm->registerHisto(trigger_fire_map);
	}
    }


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

  bool m_overrideCalibName = false;
  bool m_overrideFieldName = false;

  if (_do_emcal)
    {
      std::string default_time_independent_calib = "cemc_pi0_twrSlope_v1_default";

      if (!m_overrideCalibName)
	{
	  m_calibName_emcal = "cemc_pi0_twrSlope_v1";
	}
      if (!m_overrideFieldName)
	{
	  m_fieldname_emcal = "Femc_datadriven_qm1_correction";
	}
      std::string calibdir = CDBInterface::instance()->getUrl(m_calibName_emcal);
      if (!calibdir.empty())
	{
	  cdbttree_emcal = new CDBTTree(calibdir);
	}
      else
	{
	  calibdir = CDBInterface::instance()->getUrl(default_time_independent_calib);

	  if (calibdir.empty())
	    {
	      std::cout << "CaloTriggerEmulator::::InitRun No EMCal Calibration NOT even a default" << std::endl;
	      exit(1);
	    }
	  cdbttree_emcal = new CDBTTree(calibdir);
	  std::cout << "CaloTriggerEmulator::::InitRun No specific file for " << m_calibName_emcal << " found, using default calib " << default_time_independent_calib << std::endl;
	}
    }
  else if (_do_hcalin)
    {

      if (!m_overrideCalibName)
	{
	  m_calibName_hcalin = "ihcal_abscalib_cosmic";
	}
      if (!m_overrideFieldName)
	{
	  m_fieldname_hcalin = "ihcal_abscalib_mip";
	}
      std::string calibdir = CDBInterface::instance()->getUrl(m_calibName_hcalin);
      if (!calibdir.empty())
	{
	  cdbttree_hcalin = new CDBTTree(calibdir);
	}
      else
	{
	  std::cout << "CaloTowerCalib::::InitRun No calibration file for domain " << m_calibName_hcalin << " found" << std::endl;
	  exit(1);
	}
    }
  else if (_do_hcalout)
    {

      if (!m_overrideCalibName)
	{
	  m_calibName_hcalout = "ohcal_abscalib_cosmic";
	}
      if (!m_overrideFieldName)
	{
	  m_fieldname_hcalout = "ohcal_abscalib_mip";
	}
      std::string calibdir = CDBInterface::instance()->getUrl(m_calibName_hcalout);
      if (!calibdir.empty())
	{
	  cdbttree_hcalout = new CDBTTree(calibdir);
	}
      else
	{
	  std::cout << "CaloTowerCalib::::InitRun No calibration file for domain " << m_calibName_hcalout << " found" << std::endl;
	  exit(1);
	}
    }

  return 0;
}
// process event procedure
int CaloTriggerEmulator::process_event(PHCompositeNode* topNode)
{

  if(Verbosity()>=1) std::cout << __FUNCTION__ << ": event " <<_nevent<<std::endl;

  // Get all nodes needed fo
  GetNodes(topNode);


  // process waveforms from the waveform container into primitives
  if (process_waveforms())  return Fun4AllReturnCodes::EVENT_OK;
  
  // process all the primitives into sums.
  process_primitives();

  // calculate the true LL1 trigger algorithm.
  if (process_trigger()) return Fun4AllReturnCodes::EVENT_OK;
  
  _nevent++;

  
  if(Verbosity()>=2) identify();
  
  if (Verbosity()>=1) std::cout <<" done .."<<endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

// RESET event procedure that takes all variables to 0 and clears the primitives.
int CaloTriggerEmulator::ResetEvent(PHCompositeNode* /*topNode*/)
{
  if (Verbosity() >=5) std::cout << __FUNCTION__ << __LINE__ <<std::endl;
  // reset the primitives
  if (_primitives_emcal) _primitives_emcal->Reset();
  if (Verbosity() >=5) std::cout << __FUNCTION__ << __LINE__ <<std::endl;
  if ( _primitives_hcalin) _primitives_hcalin->Reset();
  if (Verbosity() >=5) std::cout << __FUNCTION__ << __LINE__ <<std::endl;
  if (_primitives_hcalout) _primitives_hcalout->Reset();
  if (Verbosity() >=5) std::cout << __FUNCTION__ << __LINE__ <<std::endl;
  if (_primitives) _primitives->Reset();
  if (Verbosity() >=5) std::cout << __FUNCTION__ << __LINE__ <<std::endl;
  
  if (Verbosity() >=5) std::cout << __FUNCTION__ << __LINE__ <<std::endl;
  // here, the peak minus pedestal map is cleanly disposed of
  for (auto i = m_peak_sub_ped_emcal.begin(); i  != m_peak_sub_ped_emcal.end(); ++i)
    {
      i->second->clear();
    }
  if (Verbosity() >=5) std::cout << __FUNCTION__ << __LINE__ <<std::endl;
  m_peak_sub_ped_emcal.clear();
  // here, the peak minus pedestal map is cleanly disposed of
  for (auto i = m_peak_sub_ped_hcalin.begin(); i  != m_peak_sub_ped_hcalin.end(); ++i)
    {
      i->second->clear();
    }
  if (Verbosity() >=5) std::cout << __FUNCTION__ << __LINE__ <<std::endl;
  m_peak_sub_ped_hcalin.clear();

  // here, the peak minus pedestal map is cleanly disposed of
  for (auto i = m_peak_sub_ped_hcalout.begin(); i  != m_peak_sub_ped_hcalout.end(); ++i)
    {
      i->second->clear();
    }
  if (Verbosity() >=5) std::cout << __FUNCTION__ << __LINE__ <<std::endl;
  m_peak_sub_ped_hcalout.clear();

  // here, the peak minus pedestal map is cleanly disposed of
  for (auto i = m_peak_sub_ped_mbd.begin(); i  != m_peak_sub_ped_mbd.end(); ++i)
    {
      i->second->clear();
    }
  if (Verbosity() >=5) std::cout << __FUNCTION__ << __LINE__ <<std::endl;
  m_peak_sub_ped_mbd.clear();

  // do the MBD if it is there.


  return 0;
}

int CaloTriggerEmulator::process_waveforms()
{
  // Get range of waveforms
  if (Verbosity()) std::cout << __FILE__<<"::"<<__FUNCTION__<<":: Processing waveforms"<<std::endl;
  if (_do_emcal)
    {  
      
      int peak_sub_ped;
      std::vector<int> *v_peak_sub_ped;
      std::vector<int> wave;
      if (!_waveforms_emcal->size()) return Fun4AllReturnCodes::EVENT_OK;
      // for each waveform, clauclate the peak - pedestal given the sub-delay setting
      for (unsigned int iwave = 0; iwave < (unsigned int)_waveforms_emcal->size(); iwave++)
	{	  
	  
	  v_peak_sub_ped = new std::vector<int>();
	  peak_sub_ped = 0;
	  TowerInfo *tower = _waveforms_emcal->get_tower_at_channel(iwave);
	  unsigned int key = TowerInfoDefs::encode_emcal(iwave);	  
	  for (int i = 0; i < m_nsamples - _m_trig_sub_delay;i++)
	    {
	      peak_sub_ped = tower->get_waveform_value(_m_trig_sub_delay + i) - tower->get_waveform_value(i);
	      // if negative, set to 0
	      if (peak_sub_ped < 0) peak_sub_ped = 0;
	      v_peak_sub_ped->push_back(peak_sub_ped);
	    }
	  // save in global.
	  m_peak_sub_ped_emcal[key] = v_peak_sub_ped;
	}

    }
  if (_do_hcalout)
    {  
      
      if (Verbosity()) std::cout << __FILE__<<"::"<<__FUNCTION__<<":: Processing HCALOUT"<<std::endl;
      int peak_sub_ped;
      std::vector<int> *v_peak_sub_ped;
      std::vector<int> wave;
      // for each waveform, clauclate the peak - pedestal given the sub-delay setting
      if (!_waveforms_hcalout->size()) return Fun4AllReturnCodes::EVENT_OK;

      for (unsigned int iwave = 0; iwave < (unsigned int)_waveforms_hcalout->size(); iwave++)
	{	  

	  v_peak_sub_ped = new std::vector<int>();
	  peak_sub_ped = 0;
	  TowerInfo *tower = _waveforms_hcalout->get_tower_at_channel(iwave);

	  for (int i = 0; i < m_nsamples - _m_trig_sub_delay;i++)
	    {
	      peak_sub_ped = tower->get_waveform_value(_m_trig_sub_delay + i) - tower->get_waveform_value(i);
	      // if negative, set to 0
	      if (peak_sub_ped < 0) peak_sub_ped = 0;
	      v_peak_sub_ped->push_back(peak_sub_ped);
	    }
	  if (Verbosity()) 
	    {
	      std::cout << __FILE__<<"::"<<__FUNCTION__<<":: iwave "<<iwave<<" : ";
	      for (int i = 0; i < m_nsamples - _m_trig_sub_delay;i++)
		{
		  std::cout << v_peak_sub_ped->at(i) << " ";
		}
	      std::cout <<" "<<std::endl;
	    }
	  // save in global.
	  m_peak_sub_ped_hcalout[iwave] = v_peak_sub_ped;
	}
    }
    if (_do_hcalin)
    {  
      if (!_waveforms_hcalout->size()) return Fun4AllReturnCodes::EVENT_OK; 
      int peak_sub_ped;
      std::vector<int> *v_peak_sub_ped;
      std::vector<int> wave;

      // for each waveform, clauclate the peak - pedestal given the sub-delay setting
      for (unsigned int iwave = 0; iwave < (unsigned int)_waveforms_hcalin->size(); iwave++)
	{	  

	  v_peak_sub_ped = new std::vector<int>();
	  peak_sub_ped = 0;
	  TowerInfo *tower = _waveforms_hcalin->get_tower_at_channel(iwave);

	  for (int i = 0; i < m_nsamples - _m_trig_sub_delay;i++)
	    {
	      peak_sub_ped = tower->get_waveform_value(_m_trig_sub_delay + i) - tower->get_waveform_value(i);
	      // if negative, set to 0
	      if (peak_sub_ped < 0) peak_sub_ped = 0;
	      v_peak_sub_ped->push_back(peak_sub_ped);
	    }
	  // save in global.
	  m_peak_sub_ped_hcalin[iwave] = v_peak_sub_ped;
	}
    }
    
  if (_do_mbd)
    { 
      
      
      int peak_sub_ped;
      std::vector<int> *v_peak_sub_ped;
      std::vector<int> wave;
      int ij = 0;

      if (!_waveforms_mbd->size()) return Fun4AllReturnCodes::EVENT_OK;

      // for each waveform, clauclate the peak - pedestal given the sub-delay setting
      for (unsigned int iwave = 0; iwave < (unsigned int)_waveforms_mbd->size(); iwave++)
	{	  
	  v_peak_sub_ped = new std::vector<int>();
	  peak_sub_ped = 0;
	  TowerInfo *tower = _waveforms_mbd->get_tower_at_channel(iwave);

	  for (int i = 0; i < m_nsamples - _m_trig_sub_delay;i++)
	    {
	      peak_sub_ped = tower->get_waveform_value(_m_trig_sub_delay + i) - tower->get_waveform_value(i);
	      // if negative, set to 0
	      if (peak_sub_ped < 0) peak_sub_ped = 0;
	      v_peak_sub_ped->push_back(peak_sub_ped);
	    }
	  // save in global.
	  m_peak_sub_ped_mbd[iwave] = v_peak_sub_ped;
	}
      if(Verbosity()>=1) std::cout << "Processed waves: "<<ij <<std::endl;
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

// procedure to process the peak - pedestal into primitives.  
int CaloTriggerEmulator::process_primitives()
{
  
  unsigned int ip;
  int id_peak;
  unsigned int peak;  
  int i;
  bool mask;
  if (_do_emcal)
    {
      if(Verbosity()>=1) std::cout << __FUNCTION__<<" "<<__LINE__<<"Gathering EMCAL primitives."<<std::endl;
      ip = 0;

      if (Verbosity() > 2) std::cout << __FUNCTION__<<"::"<<__LINE__<<std::endl;
      // get the number of primitives needed to process
      _n_primitives = _m_prim_map[TriggerDefs::DetectorId::emcalDId];
      for (i = 0; i < _n_primitives; i++, ip++)
	{
	  unsigned int tmp;  
	  if (Verbosity() > 2) std::cout << __FUNCTION__<<"::"<<__LINE__<<" prim "<<ip <<std::endl;
	  // get the primitive key of what we are making, in order of the packet ID and channel number
	  TriggerDefs::TriggerPrimKey primkey = TriggerDefs::getTriggerPrimKey(TriggerDefs::GetTriggerId(_trigger), TriggerDefs::GetDetectorId("EMCAL"), TriggerDefs::GetPrimitiveId("EMCAL"), ip);

	  // Make a new primitive, holds all 16 sums
	  _primitive = new TriggerPrimitivev1(primkey);
	 
	  unsigned int sum;
	  // check if masked Fiber;
	  mask = CheckFiberMasks(primkey);

	  // calculate 16 sums
	  for (int isum = 0; isum < _n_sums; isum++)
	    {
	      id_peak = -1;
	      peak = 0;

	      // get sum key
	      TriggerDefs::TriggerSumKey sumkey = TriggerDefs::getTriggerSumKey(TriggerDefs::GetTriggerId(_trigger), TriggerDefs::GetDetectorId("EMCAL"), TriggerDefs::GetPrimitiveId("EMCAL"), ip, isum);

	      // calculate sums for all samples, hense the vector.
	      _sum = new std::vector<unsigned int>();

	      // check to mask channel (if fiber masked, automatically mask the channel)
	      bool mask_channel = mask || CheckChannelMasks(sumkey);
	      for (int is = 0; is < m_nsamples - _m_trig_sub_delay; is++)
		{
		  sum = 0;	
		  
		  // if masked, just fill with 0s
		  if (!mask_channel)
		    {
		      for (int j = 0; j < 4;j++)
			{
			  if (Verbosity() > 2) std::cout << __FUNCTION__<<"::"<<__LINE__<<" tower "<<ip << " "<<isum<<" " << j <<std::endl;
			  unsigned int key = TriggerDefs::GetTowerInfoKey( TriggerDefs::GetDetectorId("EMCAL"), ip, isum, j );
			  if (Verbosity() > 2) std::cout << __FUNCTION__<<"::"<<__LINE__<<" tower "<<key <<std::endl;
			  float calib_value = cdbttree_emcal->GetFloatValue(key, m_fieldname_emcal);
			  float energy_calib = static_cast<float>( m_peak_sub_ped_emcal[key]->at(is) & 0x3fff ) * calib_value;
			  if (energy_calib < m_emcal_lut_floor) tmp = 0;
			  else
			    {
			      tmp = floor((m_emcal_lut_offset + energy_calib)/m_emcal_GeV_per_lut);
			    }
			  if (tmp > 0x3ff) tmp = 0x3ff;
			  else tmp &= 0x3ff;
			  sum += (tmp & 0x3ff);
			}
		      sum = (sum & 0xfff) >> 4;
		      if (peak < sum) 
			{
			  peak = sum;
			  id_peak = is;
			}		      
		    }
		  _sum->push_back(sum);		   
		  
		}

	      if (peak > _m_threshold)
		{
		  v_avg_primitive_emcal.at(i)->Fill(isum, peak);
		  v_primitives_emcal.at(i)->Fill(isum, peak);
			  
		  v_trigger_fire_map_emcal.at(i)->Fill(isum%4, isum/4);
		  v_peak_primitive_emcal.at(i)->Fill(isum, id_peak);
		}

	      // add sum (with sumkey) to the primitive
	      _primitive->add_sum(sumkey, _sum);
	      
	    }
	  // ad primitive to all primitives
	  _primitives_emcal->add_primitive(primkey, _primitive);
	}  

    }
  if (_do_hcalout)
    {
      if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<"Gathering HCALOUT primitives."<<std::endl;
      ip = 0;
      _n_primitives = _m_prim_map[TriggerDefs::DetectorId::hcaloutDId];
      for (i = 0; i < _n_primitives; i++, ip++)
	{
	  unsigned int tmp;  
	  TriggerDefs::TriggerPrimKey primkey = TriggerDefs::getTriggerPrimKey(TriggerDefs::GetTriggerId(_trigger), TriggerDefs::GetDetectorId("HCALOUT"), TriggerDefs::GetPrimitiveId("HCALOUT"), ip);
	  _primitive = new TriggerPrimitivev1(primkey);
	  unsigned int sum;
	  mask = CheckFiberMasks(primkey);
	  for (int isum = 0; isum < _n_sums; isum++)
	    {
	      id_peak = -1;
	      peak = 0;
	      TriggerDefs::TriggerSumKey sumkey = TriggerDefs::getTriggerSumKey(TriggerDefs::GetTriggerId(_trigger), TriggerDefs::GetDetectorId("HCALOUT"), TriggerDefs::GetPrimitiveId("HCALOUT"), ip, isum);
	      _sum = new std::vector<unsigned int>();
	      mask |= CheckChannelMasks(sumkey);
	      for (int is = 0; is < m_nsamples - _m_trig_sub_delay; is++)
		{
		  sum = 0;
		  if (!mask)
		    {
		      for (int j = 0; j < 4;j++)
			{
			  tmp = m_l1_adc_table[m_peak_sub_ped_hcalout[64*ip + isum*4 + j]->at(is) >> 4];
			  sum += (tmp & 0x3ff);
			}
		      sum = (sum & 0x3ff) >> 2;
		      if (peak < sum) 
			{
			  peak = sum;
			  id_peak = is;
			}
		    }
		  _sum->push_back(sum);
		}

	      if (peak > _m_threshold)
		{
		  v_peak_primitive_hcalout.at(i)->Fill(isum, id_peak);
		  v_avg_primitive_hcalout.at(i)->Fill(isum, peak);
		  v_primitives_hcalout.at(i)->Fill(isum, peak);
		  v_trigger_fire_map_hcalout.at(i)->Fill(isum%4, isum/4);
		}

	      _primitive->add_sum(sumkey, _sum);
	      
	    }
	  
	  _primitives_hcalout->add_primitive(primkey, _primitive);
	  if(Verbosity()>=2) std::cout << "Total primitives in hcalout: "<<_primitive->size()<<std::endl;	  
	}  

    }
  if (_do_hcalin)
    {
      if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<"Gathering HCALIN primitives."<<std::endl;
      ip = 0;
      _n_primitives = _m_prim_map[TriggerDefs::DetectorId::hcalinDId];

      for (i = 0; i < _n_primitives; i++, ip++)
	{
	  unsigned int tmp;  
	  TriggerDefs::TriggerPrimKey primkey = TriggerDefs::getTriggerPrimKey(TriggerDefs::GetTriggerId(_trigger), TriggerDefs::GetDetectorId("HCALIN"), TriggerDefs::GetPrimitiveId("HCALIN"), ip);
	  _primitive = new TriggerPrimitivev1(primkey);
	  unsigned int sum;
	  mask = CheckFiberMasks(primkey);
	  for (int isum = 0; isum < _n_sums; isum++)
	    {
	      id_peak = -1;
	      peak = 0;
	      TriggerDefs::TriggerSumKey sumkey = TriggerDefs::getTriggerSumKey(TriggerDefs::GetTriggerId(_trigger), TriggerDefs::GetDetectorId("HCALIN"), TriggerDefs::GetPrimitiveId("HCALIN"), ip, isum);
	      _sum = new std::vector<unsigned int>();
	      mask |= CheckChannelMasks(sumkey);
	      for (int is = 0; is < m_nsamples - _m_trig_sub_delay; is++)
		{
		  sum = 0;
		  if (!mask)
		    {
		      for (int j = 0; j < 4;j++)
			{
			  tmp = m_l1_adc_table[m_peak_sub_ped_hcalin[64*ip + isum*4 + j]->at(is) >> 4];
			  sum += (tmp & 0x3ff);
			}
		      sum = (sum & 0x3ff) >> 2;
		      if (peak < sum) 
			{
			  peak = sum;
			  id_peak = is;
			}
		    
		      if (peak > _m_threshold)
			{
			  v_peak_primitive_hcalin.at(i)->Fill(isum, id_peak);
			  v_avg_primitive_hcalin.at(i)->Fill(isum, peak);
			  v_primitives_hcalin.at(i)->Fill(isum, peak);
			  
			  v_trigger_fire_map_hcalin.at(i)->Fill(isum%4, isum/4);
			}
		    }
		  _sum->push_back(sum);
		}
	      _primitive->add_sum(sumkey, _sum);
	      
	    }
	  if(Verbosity()>=2) std::cout << "Total primitives in hcalin: "<<_primitive->size()<<std::endl;	  
	  _primitives_hcalin->add_primitive(primkey, _primitive);
	  
	}  
    }
  if (_do_mbd)
    {

      // MBD 
      if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<"Gathering MBD primitives."<<std::endl;

      ip = 0;

      // get number of primitives
      _n_primitives = _m_prim_map[TriggerDefs::DetectorId::mbdDId];

      for (i = 0; i < _n_primitives; i++, ip++)
	{

	  // make primitive key
	  TriggerDefs::TriggerPrimKey primkey = TriggerDefs::getTriggerPrimKey(TriggerDefs::GetTriggerId(_trigger), TriggerDefs::GetDetectorId("MBD"), TriggerDefs::GetPrimitiveId("MBD"), _n_primitives - ip);

	  // make primitive and check mask;
	  _primitive = new TriggerPrimitivev1(primkey);
	  mask = CheckFiberMasks(primkey);

	  std::vector<unsigned int> *sum_mbd = nullptr;

	  _sum_mbd.clear();
	  // make 13 sums 
	  // 8 charge, 1 hit, 4 time
	  for (int j = 0 ; j < 13; j++)
	    {
	      sum_mbd = new std::vector<unsigned int>();
	      _sum_mbd.push_back(sum_mbd);
	    }
	      
	  // iterate through samples
	  for (int is = 0; is < m_nsamples - _m_trig_sub_delay; is++) 
	    {

	      // reset variables
	      for (int j = 0; j < 8; j++)
		{
		  m_trig_charge[j] = 0;
		}
	      m_trig_nhit = 0;
	      for (int j = 0; j < 4; j++)
		{
		  m_trig_time[j] = 0;
		}

	      unsigned int tmp, tmp2;  
	      unsigned int qadd[32];

	      // for each section of the board (4 sections of 8 time and 8 charge
	      for (int isec = 0; isec < 4; isec++)
		{
		  // go through 8 charge channels
		  for (int j = 0; j < 8;j++)
		    {

		      // pass upper 10 bits of charge to get 10 bit LUt outcome
		      tmp = m_l1_adc_table[m_peak_sub_ped_mbd[ i*64 + 8 + isec*16 + j ]->at(is) >> 4];
	
		      // put upper 3 bits of the 10 bits into slewing correction later
		      qadd[isec*8+j] = (tmp & 0x380) >> 7;
		      
		      // sum up to 11 bits.
		      m_trig_charge[isec*2 + j/4] += tmp & 0x7ff;
		      
		    }
		}

	      // Now the time channels
	      for (int isec = 0; isec < 4; isec++)
		{
      
		  // 8 timing channels
		  for (int j = 0; j < 8;j++)
		    {
		      // upper 10 bits go through the LUT
		      tmp = m_l1_adc_table[m_peak_sub_ped_mbd[ i*64 + isec*16 + j ]->at(is) >> 4];
		      
		      // high bit is the hit bit
		      m_trig_nhit += (tmp & 0x200) >> 9;
		      
		      // get upper 3 bits of charge in the channel, and make it bits 9-11, the time of the chanel is the lower 9 bits from 0-8. 
		      tmp2 = m_l1_slewing_table[(qadd[isec*8+j] << 9) + (tmp & 0x01ff)];
		      
		      // attribute to the time sum
		      m_trig_time[isec] += tmp2;
		      
		    }
		}

	      // ad in the charge sums
	      for (int j = 0; j < 8; j++)
		{
		  _sum_mbd[j]->push_back(m_trig_charge[j]);
		}
	      // add in the nhits calue
	      _sum_mbd[8]->push_back(m_trig_nhit);

	      // add in the time sums
	      for (int j = 0; j < 4; j++)
		{
		  _sum_mbd[9+j]->push_back(m_trig_time[j]);
		}

	      // if the sample matches the peak given then fill some histograms

	      if (is == _idx)
		{

		  for (int j = 0; j < 8; j++)
		    {
		      v_mbd_charge[Form("h_mbd_charge_%d_%d", j, ip)]->Fill(m_trig_charge[j]);
		    }

		  v_nhit[ip]->Fill(m_trig_nhit);
		  		  
		  for (int j = 0; j < 4; j++)
		    {
		      v_mbd_time[Form("h_mbd_time_%d_%d", j, ip)]->Fill(m_trig_time[j]);
		    }
		}
	    }

	  // add to primitive object
	  for (int j = 0; j < 13; j++)
	    {
	      TriggerDefs::TriggerSumKey sumkey = TriggerDefs::getTriggerSumKey(TriggerDefs::GetTriggerId(_trigger), TriggerDefs::GetDetectorId("MBD"), TriggerDefs::GetPrimitiveId("MBD"), ip, j);
	      _primitive->add_sum(sumkey, _sum_mbd[j]);
	    }
	  _primitives->add_primitive(primkey, _primitive);
	}

      if(Verbosity()>=2) std::cout << "Total primitives in mbd: "<<_primitives->size()<<std::endl; 

    }

  return Fun4AllReturnCodes::EVENT_OK;
}


// This is where the LL1 algorithm is, everything else before was in the ADC

int CaloTriggerEmulator::process_trigger()
{
  std::vector<unsigned int> bits;


  // bits are to say whether the trigger has fired. this is what is sent to the GL1
  
  for (int is = 0; is < m_nsamples - _m_trig_sub_delay; is++)
    {
      bits.push_back(0);
    }

  // jet
  // 8x8 non-overlapping sums in the EMCAL
  // create the 8x8 non-overlapping sum

  if (_triggerid == TriggerDefs::TriggerId::jetTId)
    {

      if(Verbosity()>=2) 
	{
	  std::cout <<__FUNCTION__<<" "<<__LINE__<<" processing JET trigger , bits before: "<< _bits->size();
	}
      if (!_primitives_emcal)
	{
	  std::cout << "There is no primitive container" << std::endl;
	  return Fun4AllReturnCodes::EVENT_OK;
	}

      // Make the jet primitives
      TriggerPrimitiveContainerv1::Range range;      

      int ip;
      _primitives->setTriggerType(_triggerid);
      for (ip = 0; ip < 16; ip++ )
	{
	  TriggerDefs::TriggerPrimKey primkey = TriggerDefs::getTriggerPrimKey(_triggerid, TriggerDefs::GetDetectorId("EMCAL"), TriggerDefs::GetPrimitiveId("JET"), ip);
	  TriggerPrimitive *primitive_jet = new TriggerPrimitivev1(primkey);
	  _primitives->add_primitive(primkey, primitive_jet);
	}

      // iterate through emcal primitives and organize into the 16 jet primitives each with the 8x8 nonoverlapping sum
      range = _primitives_emcal->getTriggerPrimitives();
      if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<" emcal primitives size: "<<_primitives_emcal->size()<<std::endl; 
      
      for (TriggerPrimitiveContainerv1::Iter iter = range.first ; iter != range.second ; ++iter)
	{
	  // get key and see if masked
	  TriggerDefs::TriggerPrimKey key = (*iter).first;
	  if (CheckFiberMasks(key)) {
	    if(Verbosity()>=2) std::cout << "masked: "<<key<<std::endl;
	    continue;
	  }

	  // get the primitive (16 2x2 sums)
	  _primitive  = (*iter).second;
	  TriggerPrimitivev1::Range sumrange = _primitive->getSums(); 
	  if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<" key: "<<key<<" size: "<<_primitive->size()<<std::endl; 

	  // make sum and innitialize to 0.
	  _sum = new std::vector<unsigned int>();

	  for (int isum = 0; isum < m_nsamples - _m_trig_sub_delay; isum++)
	    {
	      _sum->push_back(0);
	    }

	  // iterate through all 16 sums and add together
	  for (TriggerPrimitivev1::Iter iter_sum = sumrange.first; iter_sum != sumrange.second; ++iter_sum)
	    {
	      TriggerDefs::TriggerSumKey sumkey = (*iter_sum).first;
	      int i = 0;
	      if (CheckChannelMasks(sumkey)) continue;
	      if(Verbosity()>=2) cout <<" sum " << sumkey << " size "<<(*iter_sum).second->size()<<std::endl;
	      for (auto it_s = (*iter_sum).second->begin(); it_s != (*iter_sum).second->end(); ++it_s)
		{
		  _sum->at(i) += ((*it_s) & 0xff);
		  i++;
		}
	    }

	  // bit shift by 16 (divide by the 16 towers) to get an 8 bit energy sum.
	  for (auto it_s = _sum->begin(); it_s != _sum->end(); ++it_s)
	    {
	      *it_s = ((*it_s >> 4) & 0xff);
	    }

	  // based on where the primitive is in the detector, the location of the jet primitive is determined, 0 through 15 in phi.
	  uint16_t ijetprim_phi = TriggerDefs::getPrimitivePhiId_from_TriggerPrimKey(key);
	  // eta determines the location of the sum within the jet primitive.
	  uint16_t ijetprim_eta = TriggerDefs::getPrimitiveEtaId_from_TriggerPrimKey(key);
	  
	  TriggerDefs::TriggerPrimKey jet_prim_key = TriggerDefs::getTriggerPrimKey(_triggerid, TriggerDefs::GetDetectorId("EMCAL"), TriggerDefs::GetPrimitiveId("JET"), ijetprim_phi/2);

	  int isum = ijetprim_eta + 12*(ijetprim_phi%2);

	  TriggerDefs::TriggerPrimKey jet_sum_key = TriggerDefs::getTriggerSumKey(_triggerid, TriggerDefs::GetDetectorId("EMCAL"), TriggerDefs::GetPrimitiveId("JET"), ijetprim_phi/2, isum);

	  // add to primitive previously made the sum of the 8x8 non-overlapping sum.
	  _primitives->get_primitive_at_key(jet_prim_key)->add_sum(jet_sum_key, _sum);

	  auto peak = max_element(_sum->begin(), _sum->end());
	  if (*peak > 0)
	    {
	      v_avg_primitive_jet.at(ijetprim_phi/2)->Fill(isum, *peak);
	      v_primitives_jet.at(ijetprim_phi/2)->Fill(isum, *peak);
	      v_trigger_fire_map_jet.at(ijetprim_phi/2)->Fill(isum%12, isum/12);
	      jet_trigger_fire_map->Fill(ijetprim_eta, ijetprim_phi);
	    }
	}

      if (Verbosity() >= 1)
	{
	  std::cout << "TotalTriggerSums: " << _primitives->size() << std::endl;
	}

      if (Verbosity() >= 2)
	{
	  std::cout << " --------- 8x8 Non-overlapping sums over EMCAL --------- " << std::endl;
	  range = _primitives->getTriggerPrimitives();
	  
	  for (TriggerPrimitiveContainerv1::Iter iter = range.first ; iter != range.second ; ++iter)
	    {
	      TriggerDefs::TriggerPrimKey key = (*iter).first;
	      if (CheckFiberMasks(key)) {
		std::cout << "masked: "<<key<<std::endl;
		continue;
	      }

	      // get the primitive (16 2x2 sums)
	      _primitive  = (*iter).second;
	      TriggerPrimitivev1::Range sumrange = _primitive->getSums(); 
	      std::cout << __FUNCTION__<<" "<<__LINE__<<" key: "<<key<<" size: "<<_primitive->size()<<std::endl; 

	      // iterate through all 16 sums and add together
	      for (TriggerPrimitivev1::Iter iter_sum = sumrange.first; iter_sum != sumrange.second; ++iter_sum)
		{
		  TriggerDefs::TriggerSumKey sumkey = (*iter_sum).first;
		  if (CheckChannelMasks(sumkey)) continue;
		  cout <<" sum " << sumkey << ": "<<std::endl;
		  for (auto it_s = (*iter_sum).second->begin(); it_s != (*iter_sum).second->end(); ++it_s)
		    {
		      std::cout << (*it_s) << " ";
		    }
		  std::cout << "" <<std::endl;
		}
	      
	      
	    }
	  
	}
      

    }
  // photon
  // 4x4 overlapping sum
  else if (_triggerid == TriggerDefs::TriggerId::photonTId)
    {

      if(Verbosity()>=2) 
	{
	  std::cout <<__FUNCTION__<<" "<<__LINE__<<" processing PHOTON trigger , bits before: "<< _bits->size();
	}
      if (!_primitives_emcal)
	{
	  std::cout << "There is no primitive container" << std::endl;
	  return Fun4AllReturnCodes::EVENT_OK;
	}

      // Make the jet primitives
      TriggerPrimitiveContainerv1::Range range;      
      TriggerPrimitivev1::Range sumrange;
      _primitives->setTriggerType(_triggerid);

      // iterate through emcal primitives and organize into the 16 jet primitives each with the 8x8 nonoverlapping sum
      range = _primitives_emcal->getTriggerPrimitives();
      if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<" emcal primitives size: "<<_primitives_emcal->size()<<std::endl; 
      TriggerDefs::TriggerSumKey temp_sum_key;            
      TriggerDefs::TriggerSumKey temp_prim_key;            
      for (TriggerPrimitiveContainerv1::Iter iter = range.first ; iter != range.second ; ++iter)
	{
	  // get key and see if masked
	  TriggerDefs::TriggerPrimKey key = (*iter).first;
	  if (CheckFiberMasks(key)) {
	    if(Verbosity()>=2) std::cout << "masked: "<<key<<std::endl;
	    continue;
	  }
	  if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<" int primitive  ";
	  uint16_t primlocid = TriggerDefs::getPrimitiveLocId_from_TriggerPrimKey(key);
	  if(Verbosity()>=2) std::cout << primlocid << std::endl;
	  bool prim_right_edge = ( primlocid%12 == 11 );
	  uint16_t topedge_primlocid = (primlocid/12 == 31 ? primlocid%12 : primlocid + 12);

	  TriggerDefs::TriggerPrimKey primkey = TriggerDefs::getTriggerPrimKey(_triggerid, TriggerDefs::GetDetectorId("EMCAL"), TriggerDefs::GetPrimitiveId("PHOTON"), primlocid);

	  if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<std::endl;
	  TriggerPrimitive *primitive_photon = new TriggerPrimitivev1(primkey);

	  // get the primitive (16 2x2 sums)
	  _primitive  = (*iter).second;
	  if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<std::endl;
	  // in this primitive we will hold a 4x4 overlapping sum for each 2x2 sum in the original primitive
	  unsigned int number_of_sums = _primitive->size();
	  if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<std::endl;
	  // iterate through all sums and calculate 4x4 overlapping sum
	  for (unsigned int isum = 0; isum < number_of_sums; isum++)
	    {
	      if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<std::endl;
	      bool right_edge = (isum%4 == 3);
	      bool top_edge = (isum/4 == 3);

	      if (right_edge && prim_right_edge) continue;
	      if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<std::endl;
	      _sum = new vector<unsigned int>();	
	      if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<std::endl;
	      for (int is = 0; is < m_nsamples - _m_trig_sub_delay; is++)
		{
		  if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<std::endl;
		  unsigned int sum = 0;

		  temp_sum_key = TriggerDefs::getTriggerSumKey(_triggerid, TriggerDefs::DetectorId::emcalDId, TriggerDefs::PrimitiveId::calPId, primlocid, isum);
		  sum += (_primitive->get_sum_at_key(temp_sum_key)->at(is) & 0xff);
		  if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<std::endl;
		  if ( right_edge )
		    {
		      temp_prim_key = TriggerDefs::getTriggerPrimKey(_triggerid, TriggerDefs::DetectorId::emcalDId, TriggerDefs::PrimitiveId::calPId, primlocid + 1);
		      temp_sum_key = TriggerDefs::getTriggerSumKey(_triggerid, TriggerDefs::DetectorId::emcalDId, TriggerDefs::PrimitiveId::calPId, primlocid + 1, (isum/4)*4);
		      sum += (_primitives_emcal->get_primitive_at_key(temp_prim_key)->get_sum_at_key(temp_sum_key)->at(is) & 0xff);
		      if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<std::endl;
		    }
		  else
		    {
		      temp_sum_key = TriggerDefs::getTriggerSumKey(_triggerid, TriggerDefs::DetectorId::emcalDId, TriggerDefs::PrimitiveId::calPId, primlocid, isum + 1);
		      sum += (_primitive->get_sum_at_key(temp_sum_key)->at(is) & 0xff);
		      if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<std::endl;
		    }
		  if ( top_edge )
		    {
		      temp_prim_key = TriggerDefs::getTriggerPrimKey(_triggerid, TriggerDefs::DetectorId::emcalDId, TriggerDefs::PrimitiveId::calPId, topedge_primlocid);
		      temp_sum_key = TriggerDefs::getTriggerSumKey(_triggerid, TriggerDefs::DetectorId::emcalDId, TriggerDefs::PrimitiveId::calPId, topedge_primlocid, isum%4);
		      sum += (_primitives_emcal->get_primitive_at_key(temp_prim_key)->get_sum_at_key(temp_sum_key)->at(is) & 0xff);
		      if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<std::endl;
		    }
		  else
		    {
		      temp_sum_key = TriggerDefs::getTriggerSumKey(_triggerid, TriggerDefs::DetectorId::emcalDId, TriggerDefs::PrimitiveId::calPId, primlocid, isum + 4);
		      sum += (_primitive->get_sum_at_key(temp_sum_key)->at(is) & 0xff);
		      if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<std::endl;
		    }

		  if ( top_edge && right_edge )
		    {
		      temp_prim_key = TriggerDefs::getTriggerPrimKey(_triggerid, TriggerDefs::DetectorId::emcalDId, TriggerDefs::PrimitiveId::calPId, topedge_primlocid + 1);
		      temp_sum_key = TriggerDefs::getTriggerSumKey(_triggerid, TriggerDefs::DetectorId::emcalDId, TriggerDefs::PrimitiveId::calPId, topedge_primlocid + 1, 0);
		      sum += (_primitives_emcal->get_primitive_at_key(temp_prim_key)->get_sum_at_key(temp_sum_key)->at(is) & 0xff);
		      if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<std::endl;
		    }
		  else if ( top_edge )
		    {
		      temp_prim_key = TriggerDefs::getTriggerPrimKey(_triggerid, TriggerDefs::DetectorId::emcalDId, TriggerDefs::PrimitiveId::calPId, topedge_primlocid);
		      temp_sum_key = TriggerDefs::getTriggerSumKey(_triggerid, TriggerDefs::DetectorId::emcalDId, TriggerDefs::PrimitiveId::calPId, topedge_primlocid, isum%4 + 1);
		      sum += (_primitives_emcal->get_primitive_at_key(temp_prim_key)->get_sum_at_key(temp_sum_key)->at(is) & 0xff);
		      if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<std::endl;
		    }
		  else if ( right_edge )
		    {
		      temp_prim_key = TriggerDefs::getTriggerPrimKey(_triggerid, TriggerDefs::DetectorId::emcalDId, TriggerDefs::PrimitiveId::calPId, primlocid + 1);
		      temp_sum_key = TriggerDefs::getTriggerSumKey(_triggerid, TriggerDefs::DetectorId::emcalDId, TriggerDefs::PrimitiveId::calPId, primlocid + 1, (isum/4 + 1)*4);
		      sum += (_primitives_emcal->get_primitive_at_key(temp_prim_key)->get_sum_at_key(temp_sum_key)->at(is) & 0xff);
		      if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<std::endl;
		    }
		  else
		    {
		      temp_sum_key = TriggerDefs::getTriggerSumKey(_triggerid, TriggerDefs::DetectorId::emcalDId, TriggerDefs::PrimitiveId::calPId, primlocid, isum + 5);
		      sum += (_primitive->get_sum_at_key(temp_sum_key)->at(is) & 0xff);
		      if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<std::endl;
		    }

		  sum = (sum >> 2);
		  _sum->push_back(sum);
		}
	      if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<std::endl;
	      primitive_photon->add_sum(TriggerDefs::getTriggerSumKey(_triggerid, TriggerDefs::DetectorId::emcalDId, TriggerDefs::PrimitiveId::photonPId, primlocid, isum), _sum);
	      auto peak = max_element(_sum->begin(), _sum->end());
	      if (*peak > 0)
		{
		  if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<< " " << primlocid << std::endl;
		  v_avg_primitive_photon.at(primlocid)->Fill(isum, *peak);
		  v_primitives_photon.at(primlocid)->Fill(isum, *peak);
		  v_trigger_fire_map_photon.at(primlocid)->Fill(isum%4, isum/4);
		
		}
	      if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<std::endl;
	    }
	  if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<std::endl;
	  _primitives->add_primitive(primkey, primitive_photon);
	}

    }

  // cosmic trigger (singles)
  else  if (_triggerid == TriggerDefs::TriggerId::cosmicTId)
    {

      if(Verbosity()>=2) 
	{
	  std::cout <<__FUNCTION__<<" "<<__LINE__<<" processing COSMIC trigger , bits before: "<< _bits->size();
	}
      if (!_primitives_hcalout)
	{
	  std::cout << "There is no primitive container" << std::endl;
	  return Fun4AllReturnCodes::EVENT_OK;
	}

      // iterating through the trigger primitives, and seeing if ANY is above threshold.
      TriggerPrimitiveContainerv1::Range range;      
      if (_do_hcalout)
	{
	  range = _primitives_hcalout->getTriggerPrimitives();
	  if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<" hcalout primitives size: "<<_primitives_hcalout->size()<<std::endl; 
	  for (TriggerPrimitiveContainerv1::Iter iter = range.first ; iter != range.second ; ++iter)
	    {
	      TriggerDefs::TriggerPrimKey key = (*iter).first;
	      if (CheckFiberMasks(key)) {
		if(Verbosity()>=2) std::cout << "masked: "<<key<<std::endl;
		continue;
	      }

	      _primitive  = (*iter).second;
	      TriggerPrimitivev1::Range sumrange = _primitive->getSums(); 
	      if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<" key: "<<key<<" size: "<<_primitive->size()<<std::endl; 
	      for (TriggerPrimitivev1::Iter iter_sum = sumrange.first; iter_sum != sumrange.second; ++iter_sum)
		{
		  TriggerDefs::TriggerSumKey sumkey = (*iter_sum).first;
		  int i = 0;
		  if (CheckChannelMasks(sumkey)) continue;
		  if(Verbosity()>=2) cout <<" sum " << sumkey << " size "<<(*iter_sum).second->size()<<std::endl;
		  for (auto it_s = (*iter_sum).second->begin(); it_s != (*iter_sum).second->end(); ++it_s)
		    {
		  
		      if ((*it_s) >= _m_threshold) {
			bits.at(i) |= 1;
		      }
		      i++;
		    }
		}
	  
	  
	    }
	}

      if (_do_hcalin)
	{
	  range = _primitives_hcalin->getTriggerPrimitives();
	  if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<" hcalout primitives size: "<<_primitives_hcalin->size()<<std::endl; 
	  for (TriggerPrimitiveContainerv1::Iter iter = range.first ; iter != range.second ; ++iter)
	    {
	      TriggerDefs::TriggerPrimKey key = (*iter).first;
	      if (CheckFiberMasks(key)) {
		if(Verbosity()>=2) std::cout << "masked: "<<key<<std::endl;
		continue;
	      }
	      _primitive = (*iter).second;
	      if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<" key: "<<key<<" size: "<<_primitive->size()<<std::endl; 
	      TriggerPrimitivev1::Range sumrange = _primitive->getSums(); 
	      for (TriggerPrimitivev1::Iter iter_sum = sumrange.first; iter_sum != sumrange.second; ++iter_sum)
		{
		  TriggerDefs::TriggerSumKey sumkey = (*iter_sum).first;
		  int i = 0;
		  if (CheckChannelMasks(sumkey)) {
		    if(Verbosity()>=2) std::cout << "masked: "<<sumkey<<std::endl;
		    continue;
		  }
		  if(Verbosity()>=2) cout <<" sum " << sumkey << " size "<<(*iter_sum).second->size()<<std::endl;
		  for (auto it_s = (*iter_sum).second->begin(); it_s != (*iter_sum).second->end(); ++it_s)
		    {

		      if ((*it_s) >= _m_threshold){
			bits.at(i) |= 1;
		      }
		      i++;
		    }
		}
	  
	  
	    }
	}
      // check if any sample passes here.
      _bits->clear();
      int pass = 0;
      if(Verbosity()>=2) std::cout << "bits after: ";
      for (int is = 0; is < m_nsamples - _m_trig_sub_delay; is++)
	{
	  if(Verbosity()>=2)std::cout <<" "<<bits.at(is);
	  _bits->push_back(bits.at(is));
	  if (bits.at(is) == 1) pass = 1;
	}
      _npassed += pass;
      if(Verbosity()>=2)std::cout <<" "<<std::endl;
    }

  // cosmic (coincidence)
  else if (_triggerid == TriggerDefs::TriggerId::cosmic_coinTId)
    {

      // organize the sums
      unsigned int cosmic_organized_sums[2][12][32];
      if(Verbosity()>=2) 
	{
	  std::cout <<__FUNCTION__<<" "<<__LINE__<<" processing COSMIC COINCIDENCE trigger , bits before: "<< _bits->size();
	}
      if (!_primitives_hcalout || !_primitives_hcalin)
	{
	  std::cout << "There is no primitive container" << std::endl;
	  return Fun4AllReturnCodes::EVENT_OK;
	}
      uint16_t icard, icosmic;
      TriggerPrimitiveContainerv1::Range range;      
      for (int isam = 0; isam < m_nsamples - _m_trig_sub_delay;isam++)
	{

	  //	  Set everything to 0 to get the sums.
	  for (int ii = 0; ii < 2; ii ++)
	    {
	      for (int iii = 0; iii < 12; iii++)
	  	{
	  	  for (int iv = 0; iv < 32; iv ++)
	  	    {
	  	      cosmic_organized_sums[ii][iii][iv] =0;
	  	    }
	  	}
	    }


	  // get all primitives and iterate
	  range = _primitives_hcalout->getTriggerPrimitives();
	  if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<" hcalout primitives size: "<<_primitives_hcalout->size()<<std::endl; 
	  for (TriggerPrimitiveContainerv1::Iter iter = range.first ; iter != range.second ; ++iter)
	    {

	      // get key
	      TriggerDefs::TriggerPrimKey key = (*iter).first;
	      // check if primitive is masked
	      if (CheckFiberMasks(key)) {
		if(Verbosity()>=2) std::cout << "masked: "<<key<<std::endl;
		continue;
	      }
	      
	      
	      _primitive  = (*iter).second;
	      // get location in index of phi and eta
	      uint16_t primphi = TriggerDefs::getPrimitivePhiId_from_TriggerPrimKey(key);
	      uint16_t primeta = TriggerDefs::getPrimitiveEtaId_from_TriggerPrimKey(key);

	      // get the card (either 0 or 1);
	      icard = (primphi < 4 ? 0 : 1);
	      TriggerPrimitivev1::Range sumrange = _primitive->getSums(); 
	      //if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<" key: "<<key<<" size: "<<_primitive->size()<<std::endl; 
	      for (TriggerPrimitivev1::Iter iter_sum = sumrange.first; iter_sum != sumrange.second; ++iter_sum)
		{

		  // get sum key
		  TriggerDefs::TriggerSumKey sumkey = (*iter_sum).first;
		  if (CheckChannelMasks(sumkey)) continue;
		  // get the integer index in phi and eta of the usm within the 8x8 area (4x4 sums).
		  uint16_t sumphi = TriggerDefs::getSumPhiId(sumkey);
		  uint16_t sumeta = TriggerDefs::getSumEtaId(sumkey);

		  // get the cosmic area
		  icosmic = (primeta*4 + (3 - sumeta))/2;
		  int isum = ((primphi%4)*8) + sumphi*2 + (1 - (sumeta%2));

		  // if(Verbosity()>=2) std::cout << "putting sum " << isam << " in " << sumeta <<" - "<<sumphi << " -> " << icard <<" , "<<icosmic<<" , "<<isum<<std::endl; 		  

		  cosmic_organized_sums[icard][icosmic][isum] = *((*iter_sum).second->begin() + isam);		  
		  
		}
	      
	      
	    }

	  // no the same for inner hcal.
	  range = _primitives_hcalin->getTriggerPrimitives();
	  if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<" hcalout primitives size: "<<_primitives_hcalin->size()<<std::endl; 
	  for (TriggerPrimitiveContainerv1::Iter iter = range.first ; iter != range.second ; ++iter)
	    {
	      TriggerDefs::TriggerPrimKey key = (*iter).first;
	      if (CheckFiberMasks(key)) {
		if(Verbosity()>=2) std::cout << "masked: "<<key<<std::endl;
		continue;
	      }
	      
	      _primitive  = (*iter).second;
	      
	      uint16_t primphi = TriggerDefs::getPrimitivePhiId_from_TriggerPrimKey(key);
	      uint16_t primeta = TriggerDefs::getPrimitiveEtaId_from_TriggerPrimKey(key);
	      icard = (primphi < 4 ? 0 : 1);
	      TriggerPrimitivev1::Range sumrange = _primitive->getSums(); 
	      //if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<" key: "<<key<<" size: "<<_primitive->size()<<std::endl; 
	      for (TriggerPrimitivev1::Iter iter_sum = sumrange.first; iter_sum != sumrange.second; ++iter_sum)
		{
		  TriggerDefs::TriggerSumKey sumkey = (*iter_sum).first;

		  uint16_t sumphi = TriggerDefs::getSumPhiId(sumkey);
		  uint16_t sumeta = TriggerDefs::getSumEtaId(sumkey);
		  
		  icosmic = 6 + (primeta*4 + (3 - sumeta))/2;
		  int isum = ((primphi%4)*8) + sumphi*2 + (1 - (sumeta%2));
		  if (CheckChannelMasks(sumkey)) continue;
		  cosmic_organized_sums[icard][icosmic][isum] = *((*iter_sum).second->begin() + isam);		  
		  
		}
	      
	      
	    }


	  // for the two cards see if there is a coincidence.
	  unsigned int hit_cosmic[2] = {0, 0};
	  for (int ica = 0; ica < 2; ica++)
	    {
	      for (int ic = 0; ic < 6; ic++)
		{
		  for (int isum = 0; isum < 32; isum++)
		    {
		      hit_cosmic[ica] |= ((cosmic_organized_sums[ica][ic][isum] > _m_threshold ? 0x1 : 0) << ic);
		      hit_cosmic[ica] |= ((cosmic_organized_sums[ica][ic+6][isum] > _m_threshold ? 0x1 : 0) << (6+ic));
		    }
		}
	    }
	  if (((m_l1_hcal_table[hit_cosmic[0]] & 0x1) == 0x1 ) && ((m_l1_hcal_table[hit_cosmic[1]] & 0x2) == 0x2)) bits.at(isam) |= 1; 
	  if (((m_l1_hcal_table[hit_cosmic[0]] & 0x2) == 0x2 ) && ((m_l1_hcal_table[hit_cosmic[1]] & 0x1) == 0x1)) bits.at(isam) |= 1; 

	  if (Verbosity() >=2 && isam == 6)
	    {
	      for (int ins = 0; ins < 2; ins++)
		{

		  for (int iphii = 0; iphii < 16; iphii++)
		    {
		      std::cout <<" "<<std::endl;
		      for (int ietaa = 0; ietaa < 12; ietaa++)		
			{
			  std:: cout << cosmic_organized_sums[ins][ietaa/2][(iphii%16)*2 + ietaa%2] << " ";
			  if (ietaa % 2) std:: cout << "| ";
			}
		      std::cout << "             ";
		      for (int ietaa = 0; ietaa < 12; ietaa++)		
			{
			  std:: cout << cosmic_organized_sums[ins][6 + ietaa/2][(iphii%16)*2 + ietaa%2] << " ";
			  if (ietaa % 2) std:: cout << "| ";
			}

		    }
		}
	    }

	}

      _bits->clear();
      int pass = 0;
      if(Verbosity()>=2) std::cout << "bits after: ";
      for (int is = 0; is < m_nsamples - _m_trig_sub_delay; is++)
	{
	  if(Verbosity()>=2)std::cout <<" "<<bits.at(is);
	  _bits->push_back(bits.at(is));
	  if (bits.at(is) == 1) pass = 1;
	}
      _npassed += pass;
      if(Verbosity()>=2)std::cout <<" "<<std::endl;
    }

  // this is the MBD trigger algorithm
  else if (_triggerid == TriggerDefs::TriggerId::mbdTId)
    {

      if(Verbosity()>=2) 
	{
	  std::cout <<__FUNCTION__<<" "<<__LINE__<<" processing MBD trigger , bits before: "<< _bits->size();
	}
      if (!_primitives)
	{
	  std::cout << "There is no primitive container" << std::endl;
	  return Fun4AllReturnCodes::EVENT_OK;
	}




      TriggerPrimitiveContainerv1::Range range;      
      TriggerPrimitivev1::Range sumrange;
      int ip, isum;

      range = _primitives->getTriggerPrimitives();
      
      if(Verbosity()>=2) std::cout << __FUNCTION__<<" "<<__LINE__<<" mbd primitives size: "<<_primitives->size()<<std::endl; 

      std::vector<unsigned int> *word_mbd = nullptr;

      _word_mbd.clear();
      for (int j = 0 ; j < 8; j++)
	{
	  word_mbd = new std::vector<unsigned int>();
	  _word_mbd.push_back(word_mbd);
	}


      for (int is = 0; is < m_nsamples - _m_trig_sub_delay; is++)
	{
	  ip = 0;
	  for (TriggerPrimitiveContainerv1::Iter iter = range.first ; iter != range.second ; ++iter, ip++)
	    {
	      _primitive  = (*iter).second;
	      sumrange = _primitive->getSums(); 
	      isum = 0;
	      for (TriggerPrimitivev1::Iter iter_sum = sumrange.first; iter_sum != sumrange.second; ++iter_sum, isum++)
		{
		  if (isum < 8)
		    {
		      m2_trig_charge[ip][isum] = ((*iter_sum).second)->at(is);
		    }
		  else if (isum == 8)  m2_trig_nhit[ip] = ((*iter_sum).second)->at(is);
		  else 
		    {
		      m2_trig_time[ip][isum - 9] = ((*iter_sum).second)->at(is);
		    }
		}
	    }

	  if (Verbosity() && is == 11) {
	    for (int q = 0; q < 8; q++)
	      {

		std::cout <<"Q"<<std::dec<<q<<": "; 
		for (int ipp = 0; ipp < 4; ipp++) std::cout << std::hex<<m2_trig_charge[ipp][q]<<" ";
		std::cout <<" "<<std::endl;
	      }
	    std::cout <<"NH: "; 
	    for (int ipp = 0; ipp < 4; ipp++) std::cout << std::hex<<m2_trig_nhit[ipp]<<" ";
	    std::cout <<" "<<std::endl;

	    for (int q = 0; q < 4; q++)
	      {

		std::cout <<"T"<<std::dec<<q<<": "; 
		for (int ipp = 0; ipp < 4; ipp++) std::cout << std::hex<<m2_trig_time[ipp][q]<<" ";
		std::cout <<" "<<std::endl;
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
	      for (int j = 0 ; j < 4 ; j++)
		{
		  m_out_tsum[0] += m2_trig_time[i][j];
		  m_out_tsum[1] += m2_trig_time[i+2][j];
		}
	      m_out_nhit[0] += m2_trig_nhit[i];
	      m_out_nhit[1] += m2_trig_nhit[i+2];
	    }

	  if (m_out_nhit[0] != 0) 
	    {
	      m_out_tavg[0] = m_out_tsum[0]/m_out_nhit[0];
	      m_out_trem[0] = m_out_tsum[0]%m_out_nhit[0]; 
	    }
	  if (m_out_nhit[1] != 0) 
	    {
	      m_out_tavg[1] = m_out_tsum[1]/m_out_nhit[1];
	      m_out_trem[1] = m_out_tsum[1]%m_out_nhit[1]; 
	    }

	  unsigned int max = m_out_tavg[0];
	  unsigned int min = m_out_tavg[1];
	  if (min > max) 
	    {
	      max = m_out_tavg[1];
	      min = m_out_tavg[0];
	    }
      
	  m_out_vtx_sub = (max - min) & 0x1ff;
	  m_out_vtx_add = (m_out_tavg[0] + m_out_tavg[1]) & 0x3ff;

	  _word_mbd[0]->push_back(m_out_tavg[0]);
	  _word_mbd[1]->push_back(m_out_tavg[1]);
	  _word_mbd[2]->push_back(m_out_nhit[0]);
	  _word_mbd[3]->push_back(m_out_nhit[1]);
	  _word_mbd[4]->push_back(m_out_trem[0]);
	  _word_mbd[5]->push_back(m_out_trem[1]);
	  _word_mbd[6]->push_back(m_out_vtx_sub);
	  _word_mbd[7]->push_back(m_out_vtx_add);

	  if (m_out_nhit[0] >= m_nhit1) bits.at(is) ^= 1 << 0;
	  if (m_out_nhit[1] >= m_nhit1) bits.at(is) ^= 1 << 1;
	  if (m_out_nhit[0] >= m_nhit2) bits.at(is) ^= 1 << 2;
	  if (m_out_nhit[1] >= m_nhit2) bits.at(is) ^= 1 << 3;
      
	  if (m_out_nhit[0] >= m_nhit1 && m_out_nhit[1] >= m_nhit1 && m_out_vtx_sub <= m_timediff1) bits.at(is) ^= 1 << 4;
	  if (m_out_nhit[0] >= m_nhit1 && m_out_nhit[1] >= m_nhit1 && m_out_vtx_sub <= m_timediff2) bits.at(is) ^= 1 << 5;
	  if (m_out_nhit[0] >= m_nhit1 && m_out_nhit[1] >= m_nhit1 && m_out_vtx_sub <= m_timediff3) bits.at(is) ^= 1 << 6;
	  if (m_out_nhit[0] >= m_nhit2 && m_out_nhit[1] >= m_nhit2 && m_out_vtx_sub <= m_timediff1) bits.at(is) ^= 1 << 7;
	  if (m_out_nhit[0] >= m_nhit2 && m_out_nhit[1] >= m_nhit2 && m_out_vtx_sub <= m_timediff2) bits.at(is) ^= 1 << 8;
	  if (m_out_nhit[0] >= m_nhit2 && m_out_nhit[1] >= m_nhit2 && m_out_vtx_sub <= m_timediff3) bits.at(is) ^= 1 << 9;

	  if (Verbosity() )   {
	    std::cout << "Trigger Word : "<<std::bitset<16>(bits.at(is)) << std::dec<<std::endl;  
	  }
	}

      _bits->clear();
      if(Verbosity()>=2) std::cout << "bits after: ";
      for (int is = 0; is < m_nsamples - _m_trig_sub_delay; is++)
	{
	  if(Verbosity()>=2)std::cout <<" "<<bits.at(is);
	  _bits->push_back(bits.at(is));
	}
      if(Verbosity()>=2)std::cout <<" "<<std::endl;
      
      for (int is = 0; is < 8; is++)
	{
	  _ll1out->add_word(is, _word_mbd[is]);
	}
    }

  else 
    {
      std::cout << "Trigger "<<_trigger<< " not implemented"<<std::endl; 
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

void CaloTriggerEmulator::GetNodes(PHCompositeNode* topNode)
{

  if(Verbosity()>=2) std::cout << __FUNCTION__ << std::endl;
  _ll1out = findNode::getClass<LL1Out>(topNode, _ll1_nodename);

  if (!_ll1out) 
    {
      std::cout << "No LL1Out found... " << std::endl;
      exit(1);
    }

  _bits = _ll1out->GetTriggerBits();

  _primitives = findNode::getClass<TriggerPrimitiveContainer>(topNode, _prim_nodename);

  if (!_primitives) 
    {
      std::cout << "No TriggerPrimitives found... " << std::endl;
      exit(1);
    }



  if (_do_hcalout)
    { 
      _waveforms_hcalout = findNode::getClass<TowerInfoContainer>(topNode, "WAVEFORM_HCALOUT");

      if (!_waveforms_hcalout) 
	{
	  std::cout << "No HCALOUT Waveforms found... " << std::endl;
	  exit(1);
	}


      _primitives_hcalout = findNode::getClass<TriggerPrimitiveContainer>(topNode, "TRIGGERPRIMITIVES_HCALOUT");

      if (!_primitives_hcalout) 
	{
	  std::cout << "No HCAL Primitives found... " << std::endl;
	  exit(1);
	}
    }

  if (_do_hcalin)
    { 
      _waveforms_hcalin = findNode::getClass<TowerInfoContainer>(topNode, "WAVEFORM_HCALIN");

      if (!_waveforms_hcalin) 
	{
	  std::cout << "No HCAL Waveforms found... " << std::endl;
	  exit(1);
	}

      _primitives_hcalin = findNode::getClass<TriggerPrimitiveContainer>(topNode, "TRIGGERPRIMITIVES_HCALIN");

      if (!_primitives_hcalin) 
	{
	  std::cout << "No HCAL Primitives found... " << std::endl;
	  exit(1);
	}
    }
  if (_do_emcal)
    { 
      _waveforms_emcal = findNode::getClass<TowerInfoContainer>(topNode, "WAVEFORM_EMCAL");

      if (!_waveforms_emcal) 
	{
	  std::cout << "No HCAL Waveforms found... " << std::endl;
	  exit(1);
	}


      _primitives_emcal = findNode::getClass<TriggerPrimitiveContainer>(topNode, "TRIGGERPRIMITIVES_EMCAL");

      if (!_primitives_emcal) 
	{
	  std::cout << "No HCAL Primitives found... " << std::endl;
	  exit(1);
	}
    }
  if (_do_mbd)
    { 
      _waveforms_mbd = findNode::getClass<TowerInfoContainer>(topNode, "WAVEFORM_TOWERS_MBD");

      if (!_waveforms_mbd) 
	{
	  std::cout << "No HCAL Waveforms found... " << std::endl;
	  exit(1);
	}

    }

  return;

}
void CaloTriggerEmulator::CreateNodes(PHCompositeNode* topNode)
{

  if(Verbosity()>=2) std::cout << __FUNCTION__ << std::endl;

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

  LL1Out *ll1out = findNode::getClass<LL1Out>(ll1Node, _ll1_nodename);
  if (!ll1out)
    {
      ll1out = new LL1Outv1(_trigger, "NONE");
      PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(ll1out, _ll1_nodename, "PHObject");
      ll1Node->addNode(LL1OutNode);
    }
  TriggerPrimitiveContainer *ll1out_prim = findNode::getClass<TriggerPrimitiveContainer>(ll1Node, _prim_nodename);
  if (!ll1out_prim)
    {
      ll1out_prim = new TriggerPrimitiveContainerv1(_trigger);
      PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(ll1out_prim, _prim_nodename, "PHObject");
      ll1Node->addNode(LL1OutNode);
    }
  
  if (_do_emcal)
    {
      std::string ll1_nodename = "TRIGGERPRIMITIVES_EMCAL";
      TriggerPrimitiveContainer *ll1out_d = findNode::getClass<TriggerPrimitiveContainer>(ll1Node, ll1_nodename);
      if (!ll1out_d)
	{
	  ll1out_d = new TriggerPrimitiveContainerv1("EMCAL");
	  PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(ll1out_d, ll1_nodename, "PHObject");
	  ll1Node->addNode(LL1OutNode);
	}

    }
  if (_do_hcalout)
    {
      std::string ll1_nodename = "TRIGGERPRIMITIVES_HCALOUT";
      TriggerPrimitiveContainer *ll1out_d = findNode::getClass<TriggerPrimitiveContainer>(ll1Node, ll1_nodename);
      if (!ll1out_d)
	{
	  ll1out_d = new TriggerPrimitiveContainerv1("HCALOUT");
	  PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(ll1out_d, ll1_nodename, "PHObject");
	  ll1Node->addNode(LL1OutNode);
	}

    }
  if (_do_hcalin)
    {
      std::string ll1_nodename = "TRIGGERPRIMITIVES_HCALIN";
      TriggerPrimitiveContainer *ll1out_d = findNode::getClass<TriggerPrimitiveContainer>(ll1Node, ll1_nodename);
      if (!ll1out_d)
	{
	  ll1out_d = new TriggerPrimitiveContainerv1("HCALIN");
	  PHIODataNode<PHObject> *LL1OutNode = new PHIODataNode<PHObject>(ll1out_d, ll1_nodename, "PHObject");
	  ll1Node->addNode(LL1OutNode);
	}
    }
}

int CaloTriggerEmulator::End(PHCompositeNode* /*topNode*/)
{
  if (Verbosity() >=5) std::cout << __FUNCTION__ << __LINE__ <<std::endl;

  std::cout << "------------------------" <<std::endl;
  std::cout << "Total passed: " <<_npassed<<"/"<<_nevent <<std::endl;
  std::cout << "------------------------" <<std::endl;

  outfile->cd();
  outfile->Write();
  outfile->Close();
  delete outfile;
  hm->dumpHistos(outfilename, "UPDATE");
  if (Verbosity() >=5) std::cout << __FUNCTION__ << __LINE__ <<std::endl;
  return 0;
}

void CaloTriggerEmulator::identify()
{
  std::cout <<  " CaloTriggerEmulator: "<< _trigger << std::endl;
  std::cout <<  " LL1Out: "<< std::endl;
  
  if(_ll1out) _ll1out->identify();

  std::cout << " Processed " << _nevent << std::endl;
}
