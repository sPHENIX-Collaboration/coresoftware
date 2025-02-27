#include "TriggerValid.h"

// Trigger includes
#include <calobase/TowerInfo.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>

#include <calotrigger/LL1Out.h>
#include <calotrigger/LL1Outv1.h>
#include <calotrigger/TriggerDefs.h>
#include <calotrigger/TriggerPrimitiveContainerv1.h>
#include <calotrigger/TriggerPrimitivev1.h>

#include <qautils/QAHistManagerDef.h>

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

#include <boost/format.hpp>

#include <cmath>     // for log10, pow, sqrt, abs, M_PI
#include <iostream>  // for operator<<, endl, basic_...
#include <limits>
#include <map>  // for operator!=, _Rb_tree_con...
#include <string>
#include <utility>  // for pair

TriggerValid::TriggerValid(const std::string& name)
  : SubsysReco(name)
{
}

int TriggerValid::Init(PHCompositeNode* /*unused*/)
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);
  // create and register your histos (all types) here

  if (m_debug)
  {
    std::cout << "In TriggerValid::Init" << std::endl;
  }
  auto h_gl1_triggers = new TH1D("h_gl1_triggers", ";Trigger Name;#Triggers", 64, -0.5, 63.5);
  hm->registerHisto(h_gl1_triggers);
  for (int i = 0; i < 8; i++)
  {
    auto h_gl1_photon_energy = new TH1D((boost::format("h_gl1_photon_energy_%d") % i).str().c_str(), ";8x8 EMCAL energy [GeV]; counts", 100, 0, 50);
    auto h_gl1_jetpatch_energy = new TH1D((boost::format("h_gl1_jetpatch_energy_%d") % i).str().c_str(), ";.8 x .8 EMCAL+HCAL energy [GeV]", 50, 0, 50);
    hm->registerHisto(h_gl1_photon_energy);
    hm->registerHisto(h_gl1_jetpatch_energy);
  }

  auto h_emu_emcal_2x2_frequency = new TH2F("h_emu_emcal_2x2_frequency", ";#eta;#phi", 48, 0, 96, 128, 0, 256);
  hm->registerHisto(h_emu_emcal_2x2_frequency);
  auto h_emu_emcal_8x8_frequency = new TH2F("h_emu_emcal_8x8_frequency", ";#eta;#phi", 12, 0, 96, 32, 0, 256);
  hm->registerHisto(h_emu_emcal_8x8_frequency);
  auto h_emu_ihcal_2x2_frequency = new TH2F("h_emu_ihcal_2x2_frequency", ";#eta;#phi", 12, 0, 24, 32, 0, 64);
  hm->registerHisto(h_emu_ihcal_2x2_frequency);
  auto h_emu_ohcal_2x2_frequency = new TH2F("h_emu_ohcal_2x2_frequency", ";#eta;#phi", 12, 0, 24, 32, 0, 64);
  hm->registerHisto(h_emu_ohcal_2x2_frequency);
  auto h_emu_hcal_2x2_frequency = new TH2F("h_emu_hcal_2x2_frequency", ";#eta;#phi", 12, 0, 24, 32, 0, 64);
  hm->registerHisto(h_emu_hcal_2x2_frequency);
  for (int i = 0; i < 4; i++)
  {
    auto h_emu_jet_frequency_trig = new TH2F((boost::format("h_emu_jet_frequency_%d") % i).str().c_str(), ";#eta;#phi", 9, 0, 9, 32, 0, 32);
    hm->registerHisto(h_emu_jet_frequency_trig);
    auto h_emu_photon_frequency_trig = new TH2F((boost::format("h_emu_photon_frequency_%d") % i).str().c_str(), ";#eta;#phi", 12, 0, 12, 32, 0, 32);
    hm->registerHisto(h_emu_photon_frequency_trig);
  }

  auto h_emu_emcal_2x2_avg_out = new TProfile2D("h_emu_emcal_2x2_avg_out", ";#eta;#phi", 48, 0, 96, 128, 0, 256);
  hm->registerHisto(h_emu_emcal_2x2_avg_out);
  auto h_emu_emcal_8x8_avg_out = new TProfile2D("h_emu_emcal_8x8_avg_out", ";#eta;#phi", 12, 0, 96, 32, 0, 256);
  hm->registerHisto(h_emu_emcal_8x8_avg_out);
  auto h_emu_ihcal_2x2_avg_out = new TProfile2D("h_emu_ihcal_2x2_avg_out", ";#eta;#phi", 12, 0, 24, 32, 0, 64);
  hm->registerHisto(h_emu_ihcal_2x2_avg_out);
  auto h_emu_ohcal_2x2_avg_out = new TProfile2D("h_emu_ohcal_2x2_avg_out", ";#eta;#phi", 12, 0, 24, 32, 0, 64);
  hm->registerHisto(h_emu_ohcal_2x2_avg_out);
  auto h_emu_hcal_2x2_avg_out = new TProfile2D("h_emu_hcal_2x2_avg_out", ";#eta;#phi", 12, 0, 24, 64, 0, 64);
  hm->registerHisto(h_emu_hcal_2x2_avg_out);
  auto h_emcal_2x2_energy_lutsum = new TH2F("h_emcal_2x2_energy_lutsum", ";LUT output; Energy [GeV]", 256, -0.5, 255.5, 200, 0, 20);
  hm->registerHisto(h_emcal_2x2_energy_lutsum);
  auto h_emcal_8x8_energy_lutsum = new TH2F("h_emcal_8x8_energy_lutsum", ";LUT output; Energy [GeV]", 256, -0.5, 255.5, 200, 0, 20);
  hm->registerHisto(h_emcal_8x8_energy_lutsum);
  auto h_hcal_2x2_energy_lutsum = new TH2F("h_hcal_2x2_energy_lutsum", ";LUT output; Energy [GeV]", 256, -0.5, 255.5, 200, 0, 20);
  hm->registerHisto(h_hcal_2x2_energy_lutsum);
  auto h_hcalin_2x2_energy_lutsum = new TH2F("h_hcalin_2x2_energy_lutsum", ";LUT output; Energy [GeV]", 256, -0.5, 255.5, 200, 0, 20);
  hm->registerHisto(h_hcalin_2x2_energy_lutsum);
  auto h_hcalout_2x2_energy_lutsum = new TH2F("h_hcalout_2x2_energy_lutsum", ";LUT output; Energy [GeV]", 256, -0.5, 255.5, 200, 0, 20);
  hm->registerHisto(h_hcalout_2x2_energy_lutsum);

  auto h_jet_energy_lutsum = new TH2F("h_jet_energy_lutsum", ";LUT output; Energy [GeV]", 4096, -0.5, 4095.5, 200, 0, 20);
  hm->registerHisto(h_jet_energy_lutsum);

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

  LL1Out* ll1out_jet = findNode::getClass<LL1Out>(topNode, "LL1OUT_JET");

  LL1Out* ll1out_photon = findNode::getClass<LL1Out>(topNode, "LL1OUT_PHOTON");

  // y  LL1Out* ll1out_pair = findNode::getClass<LL1Out>(topNode, "LL1OUT_PAIR");

  TriggerPrimitiveContainer* trigger_primitives_emcal = findNode::getClass<TriggerPrimitiveContainer>(topNode, "TRIGGERPRIMITIVES_EMCAL");

  TriggerPrimitiveContainer* trigger_primitives_emcal_ll1 = findNode::getClass<TriggerPrimitiveContainer>(topNode, "TRIGGERPRIMITIVES_EMCAL_LL1");

  TriggerPrimitiveContainer* trigger_primitives_hcalin = findNode::getClass<TriggerPrimitiveContainer>(topNode, "TRIGGERPRIMITIVES_HCALIN");

  TriggerPrimitiveContainer* trigger_primitives_hcalout = findNode::getClass<TriggerPrimitiveContainer>(topNode, "TRIGGERPRIMITIVES_HCALOUT");

  TriggerPrimitiveContainer* trigger_primitives_hcal_ll1 = findNode::getClass<TriggerPrimitiveContainer>(topNode, "TRIGGERPRIMITIVES_HCAL_LL1");

  Gl1Packet* gl1_packet = findNode::getClass<Gl1Packet>(topNode, "GL1Packet");

  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  auto h_emu_emcal_2x2_frequency = dynamic_cast<TH2*>(hm->getHisto("h_emu_emcal_2x2_frequency"));
  auto h_emu_emcal_8x8_frequency = dynamic_cast<TH2*>(hm->getHisto("h_emu_emcal_8x8_frequency"));
  auto h_emu_ihcal_2x2_frequency = dynamic_cast<TH2*>(hm->getHisto("h_emu_ihcal_2x2_frequency"));
  auto h_emu_ohcal_2x2_frequency = dynamic_cast<TH2*>(hm->getHisto("h_emu_ohcal_2x2_frequency"));
  auto h_emu_hcal_2x2_frequency = dynamic_cast<TH2*>(hm->getHisto("h_emu_hcal_2x2_frequency"));

  TH2* h_emu_jet_frequency_trig[4];
  TH2* h_emu_photon_frequency_trig[4];
  for (int i = 0; i < 4; i++)
  {
    h_emu_jet_frequency_trig[i] = dynamic_cast<TH2*>(hm->getHisto((boost::format("h_emu_jet_frequency_%d") % i).str().c_str()));
    h_emu_photon_frequency_trig[i] = dynamic_cast<TH2*>(hm->getHisto((boost::format("h_emu_photon_frequency_%d") % i).str().c_str()));
  }

  auto h_gl1_triggers = dynamic_cast<TH1*>(hm->getHisto("h_gl1_triggers"));
  TH1* h_gl1_jetpatch_energy[8];
  TH1* h_gl1_photon_energy[8];

  for (int i = 0; i < 8; i++)
  {
    h_gl1_jetpatch_energy[i] = dynamic_cast<TH1*>(hm->getHisto((boost::format("h_gl1_jetpatch_energy_%d") % i).str().c_str()));
    h_gl1_photon_energy[i] = dynamic_cast<TH1*>(hm->getHisto((boost::format("h_gl1_photon_energy_%d") % i).str().c_str()));
  }

  auto h_emu_emcal_2x2_avg_out = dynamic_cast<TProfile2D*>(hm->getHisto("h_emu_emcal_2x2_avg_out"));
  auto h_emu_emcal_8x8_avg_out = dynamic_cast<TProfile2D*>(hm->getHisto("h_emu_emcal_8x8_avg_out"));
  auto h_emu_ihcal_2x2_avg_out = dynamic_cast<TProfile2D*>(hm->getHisto("h_emu_ihcal_2x2_avg_out"));
  auto h_emu_ohcal_2x2_avg_out = dynamic_cast<TProfile2D*>(hm->getHisto("h_emu_ohcal_2x2_avg_out"));
  auto h_emu_hcal_2x2_avg_out = dynamic_cast<TProfile2D*>(hm->getHisto("h_emu_hcal_2x2_avg_out"));
  auto h_emcal_2x2_energy_lutsum = dynamic_cast<TH2*>(hm->getHisto("h_emcal_2x2_energy_lutsum"));
  auto h_emcal_8x8_energy_lutsum = dynamic_cast<TH2*>(hm->getHisto("h_emcal_8x8_energy_lutsum"));
  auto h_hcal_2x2_energy_lutsum = dynamic_cast<TH2*>(hm->getHisto("h_hcal_2x2_energy_lutsum"));
  auto h_hcalin_2x2_energy_lutsum = dynamic_cast<TH2*>(hm->getHisto("h_hcalin_2x2_energy_lutsum"));
  auto h_hcalout_2x2_energy_lutsum = dynamic_cast<TH2*>(hm->getHisto("h_hcalout_2x2_energy_lutsum"));
  auto h_jet_energy_lutsum = dynamic_cast<TH2*>(hm->getHisto("h_jet_energy_lutsum"));

  std::map<TriggerDefs::TriggerSumKey, unsigned int> v_emcal_emu_2x2 = {};

  std::map<TriggerDefs::TriggerSumKey, unsigned int> v_emcal_emu_8x8 = {};

  std::map<TriggerDefs::TriggerSumKey, unsigned int> v_hcal_emu_2x2 = {};
  std::map<TriggerDefs::TriggerSumKey, unsigned int> v_hcalin_emu_2x2 = {};
  std::map<TriggerDefs::TriggerSumKey, unsigned int> v_hcalout_emu_2x2 = {};

  std::map<TriggerDefs::TriggerSumKey, unsigned int> v_jet_emu = {};

  std::vector<int> trig_bits{};
  if (gl1_packet)
  {
    ULong64_t gl1_scaledvec = gl1_packet->lValue(0, "ScaledVector");
    for (unsigned int bit = 0; bit < 64; bit++)
    {
      if (((gl1_scaledvec >> bit) & 0x1U) == 0x1U)
      {
        trig_bits.push_back(bit);
        h_gl1_triggers->Fill(bit);
      }
    }
  }

  if (trigger_primitives_emcal)
  {
    TriggerPrimitiveContainerv1::Range range = trigger_primitives_emcal->getTriggerPrimitives();
    for (TriggerPrimitiveContainerv1::Iter iter = range.first; iter != range.second; ++iter)
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

        uint16_t sum_phi = TriggerDefs::getSumPhiId(sumk) + 4 * TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(sumk);
        uint16_t sum_eta = TriggerDefs::getSumEtaId(sumk) + 4 * TriggerDefs::getPrimitiveEtaId_from_TriggerSumKey(sumk);

        if (summ)
        {
          float i2x2x = h_emu_emcal_2x2_frequency->GetXaxis()->GetBinCenter(sum_eta + 1);
          float i2x2y = h_emu_emcal_2x2_frequency->GetYaxis()->GetBinCenter(sum_phi + 1);
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
    for (TriggerPrimitiveContainerv1::Iter iter = range.first; iter != range.second; ++iter)
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

        uint16_t sum_phi = TriggerDefs::getSumPhiId(sumk) + 2 * TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(sumk);
        uint16_t sum_eta = TriggerDefs::getSumEtaId(sumk);

        if (summ)
        {
          float i2x2x = h_emu_emcal_8x8_frequency->GetXaxis()->GetBinCenter(sum_eta + 1);
          float i2x2y = h_emu_emcal_8x8_frequency->GetYaxis()->GetBinCenter(sum_phi + 1);

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
    for (TriggerPrimitiveContainerv1::Iter iter = range.first; iter != range.second; ++iter)
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

        uint16_t sum_phi = TriggerDefs::getSumPhiId(sumk) + 4 * TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(sumk);
        uint16_t sum_eta = TriggerDefs::getSumEtaId(sumk) + 4 * TriggerDefs::getPrimitiveEtaId_from_TriggerSumKey(sumk);

        if (summ)
        {
          float i2x2x = h_emu_ihcal_2x2_frequency->GetXaxis()->GetBinCenter(sum_eta + 1);
          float i2x2y = h_emu_ihcal_2x2_frequency->GetYaxis()->GetBinCenter(sum_phi + 1);
          h_emu_ihcal_2x2_frequency->Fill(i2x2x, i2x2y, 1);
          h_emu_ihcal_2x2_avg_out->Fill(i2x2x, i2x2y, summ);
        }
        v_hcalin_emu_2x2[sumk] = summ;
      }
    }
  }

  if (trigger_primitives_hcalout)
  {
    TriggerPrimitiveContainerv1::Range range = trigger_primitives_hcalout->getTriggerPrimitives();
    for (TriggerPrimitiveContainerv1::Iter iter = range.first; iter != range.second; ++iter)
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
        uint16_t sum_phi = TriggerDefs::getSumPhiId(sumk) + 4 * TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(sumk);
        uint16_t sum_eta = TriggerDefs::getSumEtaId(sumk) + 4 * TriggerDefs::getPrimitiveEtaId_from_TriggerSumKey(sumk);

        if (summ)
        {
          float i2x2x = h_emu_ohcal_2x2_frequency->GetXaxis()->GetBinCenter(sum_eta + 1);
          float i2x2y = h_emu_ohcal_2x2_frequency->GetYaxis()->GetBinCenter(sum_phi + 1);
          h_emu_ohcal_2x2_frequency->Fill(i2x2x, i2x2y, 1);
          h_emu_ohcal_2x2_avg_out->Fill(i2x2x, i2x2y, summ);
        }
        v_hcalout_emu_2x2[sumk] = summ;
      }
    }
  }

  if (trigger_primitives_hcal_ll1)
  {
    TriggerPrimitiveContainerv1::Range range = trigger_primitives_hcal_ll1->getTriggerPrimitives();
    for (TriggerPrimitiveContainerv1::Iter iter = range.first; iter != range.second; ++iter)
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

        uint16_t sum_phi = TriggerDefs::getSumPhiId(sumk) + 2 * TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(sumk);
        uint16_t sum_eta = TriggerDefs::getSumEtaId(sumk);

        if (summ)
        {
          float i2x2x = h_emu_hcal_2x2_frequency->GetXaxis()->GetBinCenter(sum_eta + 1);
          float i2x2y = h_emu_hcal_2x2_frequency->GetYaxis()->GetBinCenter(sum_phi + 1);
          h_emu_hcal_2x2_frequency->Fill(i2x2x, i2x2y, 1);
          h_emu_hcal_2x2_avg_out->Fill(i2x2x, i2x2y, summ);
        }
        v_hcal_emu_2x2[sumk] = summ;
      }
    }
  }

  if (ll1out_photon)
  {
    for (int i = 0; i < 4; i++)
    {
      auto triggered_sums = ll1out_photon->getTriggeredSumKeys(i + 1);
      for (auto& key : triggered_sums)
      {
        unsigned int phi = TriggerDefs::getSumPhiId(key) + TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(key) * 2;
        unsigned int eta = TriggerDefs::getSumEtaId(key);

        h_emu_photon_frequency_trig[i]->Fill(eta, phi);
      }
    }
  }

  if (ll1out_jet)
  {
    for (int i = 0; i < 4; i++)
    {
      auto triggered_sums = ll1out_jet->getTriggeredSumKeys(i + 1);
      for (auto& key : triggered_sums)
      {
        int eta = (key >> 16U) & 0xffffU;
        int phi = (key & 0xffffU);
        h_emu_jet_frequency_trig[i]->Fill(eta, phi);
      }
    }
  }

  // h_emcal_2x2_energy_lutsum = new TH2F("h_emcal_2x2_energy_lutsum",";Energy [GeV];LUT output", 100, 0, 10, 64, 0, 256);
  // h_emcal_8x8_energy_lutsum = new TH2F("h_emcal_8x8_energy_lutsum",";Energy [GeV];LUT output", 100, 0, 10, 64, 0, 256);
  // h_hcal_2x2_energy_lutsum = new TH2F("h_hcal_2x2_energy_lutsum",";Energy [GeV];LUT output", 100, 0, 10, 64, 0, 256);
  // h_jet_energy_lutsum = new TH2F("h_jet_energy_lutsum",";Energy [GeV];LUT output", 100, 0, 10, 64, 0, 256*16);

  float emcal_energies[12][35]{};
  float hcal_energies[12][35]{};
  float max_energy_emcal = 0.0;
  if (towers_emcal)
  {
    // go through the emulated 2x2 map for emcal
    for (auto& it : v_emcal_emu_2x2)
    {
      unsigned int sumk = it.first;
      uint16_t sum_phi = TriggerDefs::getSumPhiId(sumk) + 4 * TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(sumk);
      uint16_t sum_eta = TriggerDefs::getSumEtaId(sumk) + 4 * TriggerDefs::getPrimitiveEtaId_from_TriggerSumKey(sumk);

      float energy_sum = 0.0;
      for (int itower = 0; itower < 4; itower++)
      {
        int ieta = sum_eta * 2 + itower % 2;
        int iphi = sum_phi * 2 + itower / 2;
        TowerInfo* tower = towers_emcal->get_tower_at_key(TowerInfoDefs::encode_emcal(ieta, iphi));
        float offlineenergy = tower->get_energy();
        if (!tower->get_isGood())
        {
          continue;
        }
        energy_sum += offlineenergy;
      }
      h_emcal_2x2_energy_lutsum->Fill(it.second, energy_sum);
    }

    // now the 8x8

    for (auto& it : v_emcal_emu_8x8)
    {
      unsigned int sumk = it.first;
      uint16_t sum_phi = TriggerDefs::getSumPhiId(sumk) + 2 * TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(sumk);
      uint16_t sum_eta = TriggerDefs::getSumEtaId(sumk);

      float energy_sum = 0.0;
      for (int itower = 0; itower < 64; itower++)
      {
        int ieta = sum_eta * 8 + itower % 8;
        int iphi = sum_phi * 8 + itower / 8;
        TowerInfo* tower = towers_emcal->get_tower_at_key(TowerInfoDefs::encode_emcal(ieta, iphi));
        float offlineenergy = tower->get_energy();
        if (!tower->get_isGood())
        {
          continue;
        }
        energy_sum += offlineenergy;
      }
      if (energy_sum > max_energy_emcal)
      {
        max_energy_emcal = energy_sum;
      }
      emcal_energies[sum_eta][sum_phi] = energy_sum;
      h_emcal_8x8_energy_lutsum->Fill(it.second, energy_sum);
    }
  }

  if (towers_hcalin || towers_hcalout)
  {
    // go through the emulated 2x2 map for emcal
    for (auto& it : v_hcal_emu_2x2)
    {
      unsigned int sumk = it.first;
      uint16_t sum_phi = TriggerDefs::getSumPhiId(sumk) + 2 * TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(sumk);
      uint16_t sum_eta = TriggerDefs::getSumEtaId(sumk);

      float energy_sum = 0.0;
      for (int itower = 0; itower < 4; itower++)
      {
        int ieta = sum_eta * 2 + itower % 2;
        int iphi = sum_phi * 2 + itower / 2;
        if (towers_hcalin)
        {
          TowerInfo* tower = towers_hcalin->get_tower_at_key(TowerInfoDefs::encode_hcal(ieta, iphi));
          float offlineenergy = tower->get_energy();

          if (!tower->get_isGood())
          {
            continue;
          }
          energy_sum += offlineenergy;
        }
        if (towers_hcalin)
        {
          TowerInfo* tower = towers_hcalout->get_tower_at_key(TowerInfoDefs::encode_hcal(ieta, iphi));
          float offlineenergy = tower->get_energy();
          if (!tower->get_isGood())
          {
            continue;
          }
          energy_sum += offlineenergy;
        }
      }
      h_hcal_2x2_energy_lutsum->Fill(it.second, energy_sum);
      hcal_energies[sum_eta][sum_phi] = energy_sum;
    }

    // go through the emulated 2x2 map for emcal
    for (auto& it : v_hcalin_emu_2x2)
    {
      unsigned int sumk = it.first;
      uint16_t sum_phi = TriggerDefs::getSumPhiId(sumk) + 2 * TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(sumk);
      uint16_t sum_eta = TriggerDefs::getSumEtaId(sumk);

      float energy_sum = 0.0;
      for (int itower = 0; itower < 4; itower++)
      {
        int ieta = sum_eta * 2 + itower % 2;
        int iphi = sum_phi * 2 + itower / 2;
        if (towers_hcalin)
        {
          TowerInfo* tower = towers_hcalin->get_tower_at_key(TowerInfoDefs::encode_hcal(ieta, iphi));
          float offlineenergy = tower->get_energy();

          if (!tower->get_isGood())
          {
            continue;
          }
          energy_sum += offlineenergy;
        }
      }
      h_hcalin_2x2_energy_lutsum->Fill(it.second, energy_sum);
    }
    // go through the emulated 2x2 map for emcal
    for (auto& it : v_hcalout_emu_2x2)
    {
      unsigned int sumk = it.first;
      uint16_t sum_phi = TriggerDefs::getSumPhiId(sumk) + 2 * TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(sumk);
      uint16_t sum_eta = TriggerDefs::getSumEtaId(sumk);

      float energy_sum = 0.0;
      for (int itower = 0; itower < 4; itower++)
      {
        int ieta = sum_eta * 2 + itower % 2;
        int iphi = sum_phi * 2 + itower / 2;

        if (towers_hcalout)
        {
          TowerInfo* tower = towers_hcalout->get_tower_at_key(TowerInfoDefs::encode_hcal(ieta, iphi));
          float offlineenergy = tower->get_energy();
          if (!tower->get_isGood())
          {
            continue;
          }
          energy_sum += offlineenergy;
        }
      }
      h_hcalout_2x2_energy_lutsum->Fill(it.second, energy_sum);
    }
  }

  float jet_energies[9][32]{};
  float max_energy_jetpatch = 0.0;
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 12; j++)
    {
      emcal_energies[j][i + 32] = emcal_energies[j][i];
      hcal_energies[j][i + 32] = hcal_energies[j][i];
    }
  }
  for (int i = 0; i < 9; i++)
  {
    for (int j = 0; j < 32; j++)
    {
      for (int k = 0; k < 16; k++)
      {
        jet_energies[i][j] += emcal_energies[i + k % 4][j + k / 4];
        jet_energies[i][j] += hcal_energies[i + k % 4][j + k / 4];
      }
    }
  }

  for (auto& it : v_jet_emu)
  {
    unsigned int sumk = it.first;
    uint16_t sum_phi = sumk & 0xffffU;
    uint16_t sum_eta = (sumk >> 16U) & 0xffffU;

    if (max_energy_jetpatch < jet_energies[sum_eta][sum_phi])
    {
      max_energy_jetpatch = jet_energies[sum_eta][sum_phi];
    }
    h_jet_energy_lutsum->Fill(it.second, jet_energies[sum_eta][sum_phi]);
  }

  for (auto& j : trig_bits)
  {
    if (j >= 16 && j < 24)
    {
      h_gl1_jetpatch_energy[j - 16]->Fill(max_energy_jetpatch);
    }
    else if (j >= 24 && j < 32)
    {
      h_gl1_photon_energy[j - 24]->Fill(max_energy_emcal);
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
