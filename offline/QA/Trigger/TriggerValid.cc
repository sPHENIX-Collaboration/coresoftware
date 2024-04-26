#include "TriggerValid.h"

// Trigger includes
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoDefs.h>
#include <calobase/TowerInfoContainer.h>
#include <calotrigger/TriggerDefs.h>
#include <calotrigger/TriggerPrimitiveContainer.h>
#include <calotrigger/TriggerPrimitive.h>
#include <calotrigger/TriggerPrimitiveContainerv1.h>
#include <calotrigger/TriggerPrimitivev1.h>
#include <calotrigger/LL1Out.h>
#include <calotrigger/LL1Outv1.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TNtuple.h>
#include <TProfile2D.h>
#include <TSystem.h>
#include <TTree.h>

#include <cmath>     // for log10, pow, sqrt, abs, M_PI
#include <iostream>  // for operator<<, endl, basic_...
#include <limits>
#include <map>       // for operator!=, _Rb_tree_con...
#include <string>
#include <utility>  // for pair

TriggerValid::TriggerValid(const std::string& name, const std::string& filename)
  : SubsysReco(name)
  , trigger("JET")
  , outfilename(filename)
{
}

TriggerValid::~TriggerValid()
{
  delete hm;
}

int TriggerValid::Init(PHCompositeNode* /*unused*/)
{
  delete hm; // this is a null pointer - make cppcheck happy
  hm = new Fun4AllHistoManager(Name());
  // create and register your histos (all types) here

  if (m_debug)
  {
    std::cout << "In TriggerValid::Init" << std::endl;
  }
  

  outfile = new TFile(outfilename.c_str(), "RECREATE");

  h_emu_emcal_2x2_frequency = new TH2F("h_emu_emcal_2x2_frequency", ";#eta;#phi", 48, 0, 96, 128, 0, 256);
  h_emu_emcal_8x8_frequency = new TH2F("h_emu_emcal_8x8_frequency", ";#eta;#phi", 12, 0, 96, 32, 0, 256);

  h_emu_ihcal_2x2_frequency = new TH2F("h_emu_ihcal_2x2_frequency", ";#eta;#phi", 12, 0, 24, 32, 0, 64);
  h_emu_ohcal_2x2_frequency = new TH2F("h_emu_ohcal_2x2_frequency", ";#eta;#phi", 12, 0, 24, 32, 0, 64);
  h_emu_hcal_2x2_frequency = new TH2F("h_emu_hcal_2x2_frequency", ";#eta;#phi", 12, 0, 24, 32, 0, 64);

  h_emu_jet_frequency = new TH2F("h_emu_jet_frequency", ";#eta;#phi", 9, 0, 9, 32, 0, 32);

  h_emcal_ll1_2x2_frequency = new TH2F("h_emcal_ll1_2x2_frequency", ";#eta;#phi", 48, 0, 96, 128, 0, 256);
  h_emcal_ll1_8x8_frequency = new TH2F("h_emcal_ll1_8x8_frequency", ";#eta;#phi", 12, 0, 96, 32, 0, 256);

  h_hcal_ll1_2x2_frequency = new TH2F("h_hcal_ll1_2x2_frequency", ";#eta;#phi", 12, 0, 24, 64, 0, 64);

  h_jet_ll1_frequency = new TH2F("h_jet_frequency", ";#eta;#phi", 9, 0, 9, 32, 0, 32);

  h_emu_emcal_2x2_avg_out = new TProfile2D("h_emu_emcal_2x2_avg_out", ";#eta;#phi", 48, 0, 96, 128, 0, 256);
  h_emu_emcal_8x8_avg_out = new TProfile2D("h_emu_emcal_8x8_avg_out", ";#eta;#phi", 12, 0, 96, 32, 0, 256);

  h_emu_ihcal_2x2_avg_out = new TProfile2D("h_emu_ihcal_2x2_avg_out", ";#eta;#phi", 12, 0, 24, 32, 0, 64);
  h_emu_ohcal_2x2_avg_out = new TProfile2D("h_emu_ohcal_2x2_avg_out", ";#eta;#phi", 12, 0, 24, 32, 0, 64);
  h_emu_hcal_2x2_avg_out = new TProfile2D("h_emu_hcal_2x2_avg_out", ";#eta;#phi", 12, 0, 24, 64, 0, 64);

  h_emu_jet_avg_out = new TProfile2D("h_emu_jet_avg_out", ";#eta;#phi", 9, 0, 9, 32, 0, 32);
  h_emu_photon_avg_out = new TProfile2D("h_emu_photon_avg_out", ";#eta;#phi", 12, 0, 12, 32, 0, 32);

  h_emcal_ll1_2x2_avg_out = new TProfile2D("h_emcal_ll1_2x2_avg_out", ";#eta;#phi", 48, 0, 96, 128, 0, 256);
  h_emcal_ll1_8x8_avg_out = new TProfile2D("h_emcal_ll1_8x8_avg_out", ";#eta;#phi", 12, 0, 96, 32, 0, 256);

  h_hcal_ll1_2x2_avg_out = new TProfile2D("h_hcal_ll1_2x2_avg_out", ";#eta;#phi", 12, 0, 24, 64, 0, 64);

  h_jet_ll1_avg_out = new TProfile2D("h_jet_ll1_avg_out", ";#eta;#phi", 9, 0, 9, 32, 0, 32);

  h_emcal_2x2_energy_lutsum = new TH2F("h_emcal_2x2_energy_lutsum",";Energy [GeV];LUT output", 100, 0, 10,128, 0, 128);
  h_emcal_8x8_energy_lutsum = new TH2F("h_emcal_8x8_energy_lutsum",";Energy [GeV];LUT output", 100, 0, 10, 64, 0, 64);
  h_hcal_2x2_energy_lutsum = new TH2F("h_hcal_2x2_energy_lutsum",";Energy [GeV];LUT output", 100, 0, 10, 64, 0, 64);
  h_jet_energy_lutsum = new TH2F("h_jet_energy_lutsum",";Energy [GeV];LUT output", 100, 0, 10, 64, 0, 64);

  h_match_emcal = new TProfile2D("h_match_emcal", ";#eta;#phi",48, 0, 48, 128, 0, 128); 
  h_match_emcal_ll1 = new TProfile2D("h_match_emcal_ll1", ";#eta;#phi",12, 0, 12, 32, 0, 32); 
  h_match_hcal_ll1 = new TProfile2D("h_match_hcal_ll1", ";#eta;#phi",12, 0, 12, 32, 0, 32); 
  h_match_jet_ll1 = new TProfile2D("h_match_jet_ll1", ";#eta;#phi",9, 0, 9, 32, 0, 32); 

  if (m_debug)
    {
      std::cout << "Leaving TriggerValid::Init" << std::endl;
    }
  return 0;
}

int TriggerValid::process_event(PHCompositeNode* topNode)
{
  _eventcounter++;

  process_towers(topNode);


  return Fun4AllReturnCodes::EVENT_OK;
}

int TriggerValid::process_towers(PHCompositeNode* topNode)
{
  if (m_debug)
    {
      std::cout << _eventcounter << std::endl;
    }

  
  TowerInfoContainer* towers_emcal = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
  
  TowerInfoContainer* towers_hcalin = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
  
  TowerInfoContainer* towers_hcalout = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");

  LL1Out* ll1out_trigger = findNode::getClass<LL1Out>(topNode, "LL1OUT_JET");

  LL1Out* ll1out_raw_trigger = findNode::getClass<LL1Out>(topNode, "LL1OUT_RAW_JET");

  TriggerPrimitiveContainer* trigger_primitives_raw_emcal = findNode::getClass<TriggerPrimitiveContainer>(topNode, "TRIGGERPRIMITIVES_RAW_EMCAL");

  TriggerPrimitiveContainer* trigger_primitives_raw_emcal_ll1 = findNode::getClass<TriggerPrimitiveContainer>(topNode, "TRIGGERPRIMITIVES_RAW_EMCAL_LL1");

  TriggerPrimitiveContainer* trigger_primitives_raw_trigger = findNode::getClass<TriggerPrimitiveContainer>(topNode, "TRIGGERPRIMITIVES_RAW_JET");

  //  TriggerPrimitiveContainer* trigger_primitives_trigger = findNode::getClass<TriggerPrimitiveContainer>(topNode, "TRIGGERPRIMITIVES_JET");

  TriggerPrimitiveContainer* trigger_primitives_emcal = findNode::getClass<TriggerPrimitiveContainer>(topNode, "TRIGGERPRIMITIVES_EMCAL");

  TriggerPrimitiveContainer* trigger_primitives_emcal_ll1 = findNode::getClass<TriggerPrimitiveContainer>(topNode, "TRIGGERPRIMITIVES_EMCAL_LL1");

  TriggerPrimitiveContainer* trigger_primitives_hcalin = findNode::getClass<TriggerPrimitiveContainer>(topNode, "TRIGGERPRIMITIVES_HCALIN");

  TriggerPrimitiveContainer* trigger_primitives_hcalout = findNode::getClass<TriggerPrimitiveContainer>(topNode, "TRIGGERPRIMITIVES_HCALOUT");

  TriggerPrimitiveContainer* trigger_primitives_hcal_ll1 = findNode::getClass<TriggerPrimitiveContainer>(topNode, "TRIGGERPRIMITIVES_HCAL_LL1");


  std::map<TriggerDefs::TriggerSumKey, unsigned int> v_emcal_ll1_2x2 = {};
  std::map<TriggerDefs::TriggerSumKey, unsigned int> v_emcal_emu_2x2 = {};

  std::map<TriggerDefs::TriggerSumKey, unsigned int> v_emcal_ll1_8x8 = {};
  std::map<TriggerDefs::TriggerSumKey, unsigned int> v_emcal_emu_8x8 = {};

  std::map<TriggerDefs::TriggerSumKey, unsigned int> v_hcal_ll1_2x2 = {};
  std::map<TriggerDefs::TriggerSumKey, unsigned int> v_hcal_emu_2x2 = {};

  std::map<TriggerDefs::TriggerSumKey, unsigned int> v_jet_ll1 = {};
  std::map<TriggerDefs::TriggerSumKey, unsigned int> v_jet_emu = {};

  if (trigger_primitives_raw_emcal)
    {

      TriggerPrimitiveContainerv1::Range range = trigger_primitives_raw_emcal->getTriggerPrimitives();
      for (TriggerPrimitiveContainerv1::Iter iter = range.first ; iter != range.second ; ++iter)
	{
	  TriggerPrimitive* trigger_primitive = (*iter).second;
	  TriggerPrimitivev1::Range srange = trigger_primitive->getSums();
	  for (TriggerPrimitive::Iter siter = srange.first; siter != srange.second; ++siter)
	    {
	      std::vector<unsigned int>* sum = (*siter).second;
	      std::vector<unsigned int>::iterator it = max_element(sum->begin(), sum->end());
	      unsigned int summ = 0;
	      unsigned int sumk = (*siter).first;
		 
	      if (it != sum->end())
		{
		  summ = (*it);
		}
	 
	      uint16_t sum_phi =  TriggerDefs::getSumPhiId(sumk) + 4*TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(sumk);
	      uint16_t sum_eta =  TriggerDefs::getSumEtaId(sumk) + 4*TriggerDefs::getPrimitiveEtaId_from_TriggerSumKey(sumk);
	      
	      if (summ)
		{
		  float i2x2x = h_emcal_ll1_2x2_frequency->GetXaxis()->GetBinCenter(sum_eta+1);
		  float i2x2y = h_emcal_ll1_2x2_frequency->GetYaxis()->GetBinCenter(sum_phi+1);
		  h_emcal_ll1_2x2_frequency->Fill(i2x2x, i2x2y, 1);
		  h_emcal_ll1_2x2_avg_out->Fill(i2x2x, i2x2y, summ);	       
		}
	      v_emcal_ll1_2x2[sumk] = summ; 
	
	    }
	}
    }
   
   
  if (trigger_primitives_raw_emcal_ll1)
    {

      TriggerPrimitiveContainerv1::Range range = trigger_primitives_raw_emcal_ll1->getTriggerPrimitives();
      for (TriggerPrimitiveContainerv1::Iter iter = range.first ; iter != range.second ; ++iter)
	{
	  TriggerPrimitive* trigger_primitive = (*iter).second;
	  TriggerPrimitivev1::Range srange = trigger_primitive->getSums();
	  for (TriggerPrimitive::Iter siter = srange.first; siter != srange.second; ++siter)
	    {
	      std::vector<unsigned int>* sum = (*siter).second;
	      std::vector<unsigned int>::iterator it = max_element(sum->begin(), sum->end());
	      unsigned int summ = 0;
	      unsigned int sumk = (*siter).first;
	 
	      if (it != sum->end())
		{
		  summ = (*it);
		}
	 
	      uint16_t sum_phi =  TriggerDefs::getSumPhiId(sumk) + 2*TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(sumk);
	      uint16_t sum_eta =  TriggerDefs::getSumEtaId(sumk);
	      
	      if (summ)
		{
		  float i2x2x = h_emcal_ll1_8x8_frequency->GetXaxis()->GetBinCenter(sum_eta+1);
		  float i2x2y = h_emcal_ll1_8x8_frequency->GetYaxis()->GetBinCenter(sum_phi+1);
		  h_emcal_ll1_8x8_frequency->Fill(i2x2x, i2x2y, 1);
		  h_emcal_ll1_8x8_avg_out->Fill(i2x2x, i2x2y, summ);	       
		}
	      v_emcal_ll1_8x8[sumk] = summ; 
	    }
	}
    }

  if (trigger_primitives_raw_trigger)
    {

      TriggerPrimitiveContainerv1::Range range = trigger_primitives_raw_trigger->getTriggerPrimitives();
      for (TriggerPrimitiveContainerv1::Iter iter = range.first ; iter != range.second ; ++iter)
	{
	  TriggerPrimitive* trigger_primitive = (*iter).second;
	  TriggerPrimitivev1::Range srange = trigger_primitive->getSums();
	  for (TriggerPrimitive::Iter siter = srange.first; siter != srange.second; ++siter)
	    {
	      std::vector<unsigned int>* sum = (*siter).second;
	      std::vector<unsigned int>::iterator it = max_element(sum->begin(), sum->end());
	      unsigned int summ = 0;
	      unsigned int sumk = (*siter).first;

	 
	 
	      if (it != sum->end())
		{
		  summ = (*it);
		}
	 
	      uint16_t sum_phi =  TriggerDefs::getSumPhiId(sumk) + 4*TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(sumk);
	      uint16_t sum_eta =  TriggerDefs::getSumEtaId(sumk) + 4*TriggerDefs::getPrimitiveEtaId_from_TriggerSumKey(sumk);
	 
	      if (summ)
		{
		  float i2x2x = h_hcal_ll1_2x2_frequency->GetXaxis()->GetBinCenter(sum_eta+1);
		  float i2x2y = h_hcal_ll1_2x2_frequency->GetYaxis()->GetBinCenter(sum_phi+1);
		  h_hcal_ll1_2x2_frequency->Fill(i2x2x, i2x2y, 1);
		  h_hcal_ll1_2x2_avg_out->Fill(i2x2x, i2x2y, summ);	       
		}
	      v_hcal_ll1_2x2[sumk] = summ; 
	    }
	}
    }

  if (trigger_primitives_emcal)
    {

      TriggerPrimitiveContainerv1::Range range = trigger_primitives_emcal->getTriggerPrimitives();
      for (TriggerPrimitiveContainerv1::Iter iter = range.first ; iter != range.second ; ++iter)
	{
	  TriggerPrimitive* trigger_primitive = (*iter).second;
	  TriggerPrimitivev1::Range srange = trigger_primitive->getSums();
	  for (TriggerPrimitive::Iter siter = srange.first; siter != srange.second; ++siter)
	    {
	      std::vector<unsigned int>* sum = (*siter).second;
	      std::vector<unsigned int>::iterator it = max_element(sum->begin(), sum->end());
	      unsigned int summ = 0;
	      unsigned int sumk = (*siter).first;
	 
	      if (it != sum->end())
		{
		  summ = (*it);		  
		}
	 

	      uint16_t sum_phi =  TriggerDefs::getSumPhiId(sumk) + 4*TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(sumk);
	      uint16_t sum_eta =  TriggerDefs::getSumEtaId(sumk) + 4*TriggerDefs::getPrimitiveEtaId_from_TriggerSumKey(sumk);

	 
	      if (summ)
		{
		  float i2x2x = h_emu_emcal_2x2_frequency->GetXaxis()->GetBinCenter(sum_eta+1);
		  float i2x2y = h_emu_emcal_2x2_frequency->GetYaxis()->GetBinCenter(sum_phi+1);
		  h_emu_emcal_2x2_frequency->Fill(i2x2x, i2x2y, 1);
		  h_emu_emcal_2x2_avg_out->Fill(i2x2x, i2x2y, summ);	       
		}
	      v_emcal_emu_2x2[sumk] = summ; 
	    }
	}
    }
   
   
  if (trigger_primitives_emcal_ll1)
    {

      TriggerPrimitiveContainerv1::Range range = trigger_primitives_emcal_ll1->getTriggerPrimitives();
      for (TriggerPrimitiveContainerv1::Iter iter = range.first ; iter != range.second ; ++iter)
	{
	  TriggerPrimitive* trigger_primitive = (*iter).second;
	  TriggerPrimitivev1::Range srange = trigger_primitive->getSums();
	  for (TriggerPrimitive::Iter siter = srange.first; siter != srange.second; ++siter)
	    {
	      std::vector<unsigned int>* sum = (*siter).second;
	      std::vector<unsigned int>::iterator it = max_element(sum->begin(), sum->end());
	      unsigned int summ = 0;
	      unsigned int sumk = (*siter).first;
	 
	      if (it != sum->end())
		{
		  summ = (*it);
		}
	 
	      uint16_t sum_phi =  TriggerDefs::getSumPhiId(sumk) + 2*TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(sumk);
	      uint16_t sum_eta =  TriggerDefs::getSumEtaId(sumk);
	 
	      if (summ)
		{
		  float i2x2x = h_emu_emcal_8x8_frequency->GetXaxis()->GetBinCenter(sum_eta+1);
		  float i2x2y = h_emu_emcal_8x8_frequency->GetYaxis()->GetBinCenter(sum_phi+1);
		  h_emu_emcal_8x8_frequency->Fill(i2x2x, i2x2y, 1);
		  h_emu_emcal_8x8_avg_out->Fill(i2x2x, i2x2y, summ);	       
		}
	      v_emcal_emu_8x8[sumk] = summ; 
	    }
	}
    }

  if (trigger_primitives_hcalin)
    {

      TriggerPrimitiveContainerv1::Range range = trigger_primitives_hcalin->getTriggerPrimitives();
      for (TriggerPrimitiveContainerv1::Iter iter = range.first ; iter != range.second ; ++iter)
	{
	  TriggerPrimitive* trigger_primitive = (*iter).second;
	  TriggerPrimitivev1::Range srange = trigger_primitive->getSums();
	  for (TriggerPrimitive::Iter siter = srange.first; siter != srange.second; ++siter)
	    {
	      std::vector<unsigned int>* sum = (*siter).second;
	      std::vector<unsigned int>::iterator it = max_element(sum->begin(), sum->end());
	      unsigned int summ = 0;
	      unsigned int sumk = (*siter).first;
	 
	 
	      if (it != sum->end())
		{
		  summ = (*it);		  
		}
	 
	      uint16_t sum_phi =  TriggerDefs::getSumPhiId(sumk) + 4*TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(sumk);
	      uint16_t sum_eta =  TriggerDefs::getSumEtaId(sumk) + 4*TriggerDefs::getPrimitiveEtaId_from_TriggerSumKey(sumk);
	 
	      if (summ)
		{
		  float i2x2x = h_emu_ihcal_2x2_frequency->GetXaxis()->GetBinCenter(sum_eta+1);
		  float i2x2y = h_emu_ihcal_2x2_frequency->GetYaxis()->GetBinCenter(sum_phi+1);
		  h_emu_ihcal_2x2_frequency->Fill(i2x2x, i2x2y, 1);
		  h_emu_ihcal_2x2_avg_out->Fill(i2x2x, i2x2y, summ);	       
		}
	    }
	}
    }
  if (trigger_primitives_hcalout)
    {

      TriggerPrimitiveContainerv1::Range range = trigger_primitives_hcalout->getTriggerPrimitives();
      for (TriggerPrimitiveContainerv1::Iter iter = range.first ; iter != range.second ; ++iter)
	{
	  TriggerPrimitive* trigger_primitive = (*iter).second;
	  TriggerPrimitivev1::Range srange = trigger_primitive->getSums();
	  for (TriggerPrimitive::Iter siter = srange.first; siter != srange.second; ++siter)
	    {
	      std::vector<unsigned int>* sum = (*siter).second;
	      std::vector<unsigned int>::iterator it = max_element(sum->begin(), sum->end());
	      unsigned int summ = 0;
	      unsigned int sumk = (*siter).first;
	 
	 
	      if (it != sum->end())
		{
		  summ = (*it);		  
		}
	      uint16_t sum_phi =  TriggerDefs::getSumPhiId(sumk) + 4*TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(sumk);
	      uint16_t sum_eta =  TriggerDefs::getSumEtaId(sumk) + 4*TriggerDefs::getPrimitiveEtaId_from_TriggerSumKey(sumk);
	 
	      if (summ)
		{
		  float i2x2x = h_emu_ohcal_2x2_frequency->GetXaxis()->GetBinCenter(sum_eta+1);
		  float i2x2y = h_emu_ohcal_2x2_frequency->GetYaxis()->GetBinCenter(sum_phi+1);
		  h_emu_ohcal_2x2_frequency->Fill(i2x2x, i2x2y, 1);
		  h_emu_ohcal_2x2_avg_out->Fill(i2x2x, i2x2y, summ);	       
		}
	    }
	}
    }
  if (trigger_primitives_hcal_ll1)
    {

      TriggerPrimitiveContainerv1::Range range = trigger_primitives_hcal_ll1->getTriggerPrimitives();
      for (TriggerPrimitiveContainerv1::Iter iter = range.first ; iter != range.second ; ++iter)
	{
	  TriggerPrimitive* trigger_primitive = (*iter).second;
	  TriggerPrimitivev1::Range srange = trigger_primitive->getSums();
	  for (TriggerPrimitive::Iter siter = srange.first; siter != srange.second; ++siter)
	    {
	      std::vector<unsigned int>* sum = (*siter).second;
	      std::vector<unsigned int>::iterator it = max_element(sum->begin(), sum->end());
	      unsigned int summ = 0;
	      unsigned int sumk = (*siter).first;
	 
	 
	      if (it != sum->end())
		{
		  summ = (*it);
		}

	      uint16_t sum_phi =  TriggerDefs::getSumPhiId(sumk) + 4*TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(sumk);
	      uint16_t sum_eta =  TriggerDefs::getSumEtaId(sumk) + 4*TriggerDefs::getPrimitiveEtaId_from_TriggerSumKey(sumk);
		 
	      if (summ)
		{
		  float i2x2x = h_emu_hcal_2x2_frequency->GetXaxis()->GetBinCenter(sum_eta+1);
		  float i2x2y = h_emu_hcal_2x2_frequency->GetYaxis()->GetBinCenter(sum_phi+1);
		  h_emu_hcal_2x2_frequency->Fill(i2x2x, i2x2y, 1);
		  h_emu_hcal_2x2_avg_out->Fill(i2x2x, i2x2y, summ);	       
		}
	      v_hcal_emu_2x2[sumk] = summ; 
	    }
	}
    }


  // if (trigger_primitives_trigger)
  //   {

  //     TriggerPrimitiveContainerv1::Range range = trigger_primitives_trigger->getTriggerPrimitives();
  //     for (TriggerPrimitiveContainerv1::Iter iter = range.first ; iter != range.second ; ++iter)
  // 	{
  // 	  TriggerPrimitive* trigger_primitive = (*iter).second;
  // 	  TriggerPrimitivev1::Range srange = trigger_primitive->getSums();
  // 	  for (TriggerPrimitive::Iter siter = srange.first; siter != srange.second; ++siter)
  // 	    {
  // 	      std::vector<unsigned int>* sum = (*siter).second;
  // 	      std::vector<unsigned int>::iterator it = max_element(sum->begin(), sum->end());
  // 	      unsigned int summ = 0;
  // 	      unsigned int sumk = 0;
	 
	 
  // 	      if (it != sum->end())
  // 		{
  // 		  summ = (*it);
  // 		  sumk = (*siter).first;
  // 		}
	 
  // 	      uint16_t sum_phi =  TriggerDefs::getSumPhiId(sumk) + 4*TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(sumk);
  // 	      uint16_t sum_eta =  TriggerDefs::getSumEtaId(sumk) + 4*TriggerDefs::getPrimitiveEtaId_from_TriggerSumKey(sumk);
		 
  // 	      if (summ)
  // 		{
  // 		  float i2x2x = h_emu_hcal_2x2_frequency->GetXaxis()->GetBinCenter(sum_eta+1);
  // 		  float i2x2y = h_emu_hcal_2x2_frequency->GetYaxis()->GetBinCenter(sum_phi+1);
  // 		  h_emu_hcal_2x2_frequency->Fill(i2x2x, i2x2y, 1);
  // 		  h_emu_hcal_2x2_avg_out->Fill(i2x2x, i2x2y, summ);	       
  // 		}
  // 	    }
  // 	}
  //   }
   

  if (ll1out_trigger)
    {     
      LL1Outv1::Range range = ll1out_trigger->getTriggerWords();
      for (LL1Outv1::Iter iter = range.first ; iter != range.second ; ++iter)
	{
	  std::vector<unsigned int> *trigger_word = (*iter).second;
       
	  std::vector<unsigned int>::iterator it = max_element(trigger_word->begin(), trigger_word->end());
	  unsigned int summ = 0;
	  unsigned int sumk = 0;
       
	      
	  if (it != trigger_word->end())
	    {
	      summ = (*it);
	      sumk = (*iter).first;
	    }

	  if (summ)
	    {
	      int ijetx = (sumk >> 16U) & 0xffffU;
	      int ijety =(sumk & 0xffffU);

	      h_emu_jet_frequency->Fill(ijetx, ijety, 1);
	      h_emu_jet_avg_out->Fill(ijetx, ijety, summ);
	    }
	  v_jet_emu[sumk] = summ;
	
	}
    }
  if (ll1out_raw_trigger)
    {     
      LL1Outv1::Range range = ll1out_raw_trigger->getTriggerWords();
      for (LL1Outv1::Iter iter = range.first ; iter != range.second ; ++iter)
	{
	  std::vector<unsigned int> *trigger_word = (*iter).second;
       
	  std::vector<unsigned int>::iterator it = max_element(trigger_word->begin(), trigger_word->end());
	  unsigned int summ = 0;
	  unsigned int sumk = 0;
       
	      
	  if (it != trigger_word->end())
	    {
	      summ = (*it);
	      sumk = (*iter).first;
	    }

	  if (summ)
	    {
	      unsigned int ijetx = (sumk >> 16U) & 0xffffU;
	      unsigned int ijety =(sumk & 0xffffU);

	      h_jet_ll1_frequency->Fill(ijetx, ijety, 1);
	      h_jet_ll1_avg_out->Fill(ijetx, ijety, summ);
	      v_jet_ll1[sumk] = summ;
	    }
	
	}
    }
   
  // match emulated with real
  for (auto & it : v_emcal_ll1_2x2)
    {
      unsigned int sumk = it.first;
      uint16_t sum_phi =  TriggerDefs::getSumPhiId(sumk) + 4*TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(sumk);
      uint16_t sum_eta =  TriggerDefs::getSumEtaId(sumk) + 4*TriggerDefs::getPrimitiveEtaId_from_TriggerSumKey(sumk);

      std::map<TriggerDefs::TriggerSumKey, unsigned int>::iterator itt = v_emcal_emu_2x2.find(it.first);
      if (itt != v_emcal_emu_2x2.end())
	{
	  h_match_emcal->Fill(sum_eta, sum_phi, (it.second == (*itt).second ? 1 : 0));
	}
      else
	{
	  std::cout << "emcal 2x2: Trigger Sum " << std::hex<<it.first <<" does not exist -- > " << (*v_emcal_emu_2x2.begin()).first<<std::endl;
	  exit(1);
	}     
    }
  for (auto & it : v_emcal_ll1_8x8)
    {
      unsigned int sumk = it.first;
      uint16_t sum_phi =  TriggerDefs::getSumPhiId(sumk) + 4*TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(sumk);
      uint16_t sum_eta =  TriggerDefs::getSumEtaId(sumk);

      std::map<TriggerDefs::TriggerSumKey, unsigned int>::iterator itt = v_emcal_emu_8x8.find(it.first);
      if (itt != v_emcal_emu_8x8.end())
	{
	  h_match_emcal->Fill(sum_eta, sum_phi, (it.second == (*itt).second ? 1 : 0));
	}
      else
	{
	  std::cout << "emcal 8x8: Trigger Sum " << std::hex<<it.first <<" does not exist --> "<< (*v_emcal_emu_8x8.begin()).first <<std::endl;
	  exit(1);
	} 
    }
  for ( auto & it : v_hcal_ll1_2x2 )
    {
      unsigned int sumk = it.first;
      uint16_t sum_phi =  TriggerDefs::getSumPhiId(sumk) + 4*TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(sumk);
      uint16_t sum_eta =  TriggerDefs::getSumEtaId(sumk);
      std::map<TriggerDefs::TriggerSumKey, unsigned int>::iterator itt = v_hcal_emu_2x2.find(it.first);
      if (itt != v_hcal_emu_2x2.end())
	{
	  h_match_hcal_ll1->Fill(sum_eta, sum_phi, (it.second == (*itt).second ? 1 : 0));
	}
      else
	{
	  std::cout << "hcal 2x2: Trigger Sum " << std::hex<<it.first <<" does not exist -- > "<< (*v_hcal_emu_2x2.begin()).first << " is here" << std::endl;
	  exit(1);
	}
     
    }
  for ( auto & it : v_jet_ll1 )
    {
      uint16_t sum_eta =( it.first >> 16U ) & 0xffffU;
      uint16_t sum_phi = ( it.first ) & 0xffffU;

      std::map<TriggerDefs::TriggerSumKey, unsigned int>::iterator itt = v_jet_emu.find(it.first);
      if (itt != v_jet_emu.end())
	{
	  h_match_jet_ll1->Fill(sum_eta, sum_phi, (it.second == (*itt).second ? 1 : 0));
	}
      else
	{
	  std::cout << "jet: Trigger Sum " << std::hex<<it.first <<" does not exist -->"<< (*v_jet_emu.begin()).first << std::endl;
	  exit(1);
	}
     
    }

  // h_emcal_2x2_energy_lutsum = new TH2F("h_emcal_2x2_energy_lutsum",";Energy [GeV];LUT output", 100, 0, 10, 64, 0, 256);
  // h_emcal_8x8_energy_lutsum = new TH2F("h_emcal_8x8_energy_lutsum",";Energy [GeV];LUT output", 100, 0, 10, 64, 0, 256);
  // h_hcal_2x2_energy_lutsum = new TH2F("h_hcal_2x2_energy_lutsum",";Energy [GeV];LUT output", 100, 0, 10, 64, 0, 256);
  // h_jet_energy_lutsum = new TH2F("h_jet_energy_lutsum",";Energy [GeV];LUT output", 100, 0, 10, 64, 0, 256*16);

  float emcal_energies[12][35]{};
  float hcal_energies[12][35]{};
  if (towers_emcal)
    {
      // go through the emulated 2x2 map for emcal
      for ( auto & it : v_emcal_emu_2x2)
	{
	  unsigned int sumk = it.first;
	  uint16_t sum_phi =  TriggerDefs::getSumPhiId(sumk) + 4*TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(sumk);
	  uint16_t sum_eta =  TriggerDefs::getSumEtaId(sumk) + 4*TriggerDefs::getPrimitiveEtaId_from_TriggerSumKey(sumk);
	  
	  float energy_sum = 0.0;
	  for (int itower = 0; itower < 4; itower++)
	    {
	      int ieta = sum_eta*2 + itower%2;
	      int iphi = sum_phi*2 + itower/2;
	      TowerInfo* tower = towers_emcal->get_tower_at_key(TowerInfoDefs::encode_emcal(ieta, iphi));
	      float offlineenergy = tower->get_energy();
	      energy_sum += offlineenergy;
	    }
	  h_emcal_2x2_energy_lutsum->Fill(energy_sum, it.second);
	}

      // now the 8x8
      for ( auto & it : v_emcal_emu_8x8)
	{
	  unsigned int sumk = it.first;
	  uint16_t sum_phi =  TriggerDefs::getSumPhiId(sumk) + 2*TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(sumk);
	  uint16_t sum_eta =  TriggerDefs::getSumEtaId(sumk);
	  
	  float energy_sum = 0.0;
	  for (int itower = 0; itower < 64; itower++)
	    {
	      int ieta = sum_eta*8 + itower%8;
	      int iphi = sum_phi*8 + itower/8;
	      TowerInfo* tower = towers_emcal->get_tower_at_key(TowerInfoDefs::encode_emcal(ieta, iphi));
	      float offlineenergy = tower->get_energy();
	      energy_sum += offlineenergy;
	    }
	  emcal_energies[sum_eta][sum_phi] = energy_sum;
	  h_emcal_8x8_energy_lutsum->Fill(energy_sum, it.second);
	}

    }

  if (towers_hcalin || towers_hcalout)
    {
      // go through the emulated 2x2 map for emcal
      for ( auto & it : v_hcal_emu_2x2)
	{
	  unsigned int sumk = it.first;
	  uint16_t sum_phi =  TriggerDefs::getSumPhiId(sumk) + 4*TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(sumk);
	  uint16_t sum_eta =  TriggerDefs::getSumEtaId(sumk) + 4*TriggerDefs::getPrimitiveEtaId_from_TriggerSumKey(sumk);
	  
	  float energy_sum = 0.0;
	  for (int itower = 0; itower < 4; itower++)
	    {

	      int ieta = sum_eta*2 + itower%2;
	      int iphi = sum_phi*2 + itower/2;
	      if (towers_hcalin)
		{
		  TowerInfo* tower = towers_hcalin->get_tower_at_key(TowerInfoDefs::encode_hcal(ieta, iphi));
		  float offlineenergy = tower->get_energy();
		  energy_sum += offlineenergy;
		}
	      if (towers_hcalin)
		{
		  TowerInfo* tower = towers_hcalout->get_tower_at_key(TowerInfoDefs::encode_hcal(ieta, iphi));
		  float offlineenergy = tower->get_energy();
		  energy_sum += offlineenergy;
		}

	    }
	  h_hcal_2x2_energy_lutsum->Fill(energy_sum, it.second);
	  hcal_energies[sum_eta][sum_phi] = energy_sum;
	}
    }
  float jet_energies[9][32]{};


  for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 12; j++)
	{
	  emcal_energies[j][i+32] = emcal_energies[j][i];
	  hcal_energies[j][i+32] = hcal_energies[j][i];
	}
    }
  for (int i = 0; i < 9; i++)
    {
      for (int j = 0; j < 32; j++)
	{
	  for (int k = 0; k < 16; k++)
	    {
	      jet_energies[i][j] += emcal_energies[i + k%4][j + k/4];
	      jet_energies[i][j] += hcal_energies[i + k%4][j + k/4];
	    }
	}
    }

  for ( auto & it : v_jet_emu)
    {
      unsigned int sumk = it.first;
      uint16_t sum_phi =  sumk & 0xffffU;
      uint16_t sum_eta =  (sumk >> 16U) & 0xffffU;

      h_jet_energy_lutsum->Fill(jet_energies[sum_eta][sum_phi], it.second);
    }
      
  return Fun4AllReturnCodes::EVENT_OK;
}
   
int TriggerValid::End(PHCompositeNode* /*topNode*/)
{
  outfile->cd();

  outfile->Write();
  outfile->Close();
  delete outfile;
  hm->dumpHistos(outfilename, "UPDATE");
  return 0;
}
