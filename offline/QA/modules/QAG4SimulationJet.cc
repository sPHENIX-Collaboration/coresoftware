#include "QAG4SimulationJet.h"

#include "QAHistManagerDef.h"

#include <g4eval/JetEvalStack.h>
#include <g4eval/JetRecoEval.h>
#include <g4eval/JetTruthEval.h>

#include <g4jets/Jet.h>
#include <g4jets/JetMap.h>

#include <g4main/PHG4HitDefs.h>
#include <g4main/PHG4Shower.h>

#include <fun4all/Fun4AllBase.h>
#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>  // for SubsysReco

#include <phool/getClass.h>

#include <TAxis.h>
#include <TH1.h>
#include <TH2.h>
#include <TNamed.h>
#include <TString.h>  // for operator+, TString, Form

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <set>

QAG4SimulationJet::QAG4SimulationJet(const std::string& truth_jet,
                                     enu_flags flags)
  : SubsysReco("QAG4SimulationJet_" + truth_jet)
  , _jetevalstacks()
  , _truth_jet(truth_jet)
  , _reco_jets()
  , _flags(flags)
  , eta_range(-1, 1)
  , _jet_match_dEta(.1)
  , _jet_match_dPhi(.1)
  , _jet_match_dE_Ratio(.5)
{
}

int QAG4SimulationJet::InitRun(PHCompositeNode* topNode)
{
  if (flag(kProcessTruthMatching) || flag(kProcessRecoSpectrum))
  {
    for (std::set<std::string>::const_iterator it_reco_jets = _reco_jets.begin();
         it_reco_jets != _reco_jets.end(); ++it_reco_jets)
    {
      const std::string& reco_jet = *it_reco_jets;

      jetevalstacks_map::iterator it_jetevalstack = _jetevalstacks.find(
          reco_jet);

      if (it_jetevalstack == _jetevalstacks.end())
      {
        _jetevalstacks[reco_jet] = std::shared_ptr<JetEvalStack>(new JetEvalStack(topNode, reco_jet, _truth_jet));
        assert(_jetevalstacks[reco_jet]);
        _jetevalstacks[reco_jet]->set_strict(true);
        _jetevalstacks[reco_jet]->set_verbosity(Verbosity() + 1);
      }
      else
      {
        assert(it_jetevalstack->second);
      }
    }
  }

  if (flag(kProcessTruthSpectrum))
  {
    if (not _jettrutheval)
      _jettrutheval = std::shared_ptr<JetTruthEval>(new JetTruthEval(topNode, _truth_jet));

    assert(_jettrutheval);
    _jettrutheval->set_strict(true);
    _jettrutheval->set_verbosity(Verbosity() + 1);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int QAG4SimulationJet::Init(PHCompositeNode* topNode)
{
  Fun4AllHistoManager* hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  if (flag(kProcessTruthSpectrum))
  {
    if (Verbosity() >= 1)
      std::cout << "QAG4SimulationJet::Init - Process TruthSpectrum " << _truth_jet
                << std::endl;
    Init_Spectrum(topNode, _truth_jet);
  }

  if (flag(kProcessRecoSpectrum))
  {
    for (std::set<std::string>::const_iterator it_reco_jets = _reco_jets.begin();
         it_reco_jets != _reco_jets.end(); ++it_reco_jets)
    {
      const std::string& reco_jet = *it_reco_jets;
      if (Verbosity() >= 1)
        std::cout << "QAG4SimulationJet::Init - Process Reco jet spectrum "
                  << reco_jet << std::endl;
      Init_Spectrum(topNode, reco_jet);
    }
  }

  if (flag(kProcessTruthMatching))
  {
    for (std::set<std::string>::const_iterator it_reco_jets = _reco_jets.begin();
         it_reco_jets != _reco_jets.end(); ++it_reco_jets)
    {
      const std::string& reco_jet = *it_reco_jets;
      if (Verbosity() >= 1)
        std::cout << "QAG4SimulationJet::Init - Process Reco jet spectrum "
                  << reco_jet << std::endl;
      Init_TruthMatching(topNode, reco_jet);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int QAG4SimulationJet::process_event(PHCompositeNode* topNode)
{
  if (Verbosity() > 2)
    std::cout << "QAG4SimulationJet::process_event() entered" << std::endl;

  for (jetevalstacks_map::iterator it_jetevalstack = _jetevalstacks.begin();
       it_jetevalstack != _jetevalstacks.end(); ++it_jetevalstack)
  {
    assert(it_jetevalstack->second);
    it_jetevalstack->second->next_event(topNode);
  }
  if (_jettrutheval)
    _jettrutheval->next_event(topNode);

  if (flag(kProcessTruthSpectrum))
  {
    if (Verbosity() >= 1)
      std::cout << "QAG4SimulationJet::process_event - Process TruthSpectrum "
                << _truth_jet << std::endl;
    process_Spectrum(topNode, _truth_jet, false);
  }

  if (flag(kProcessRecoSpectrum))
  {
    for (std::set<std::string>::const_iterator it_reco_jets = _reco_jets.begin();
         it_reco_jets != _reco_jets.end(); ++it_reco_jets)
    {
      const std::string& reco_jet = *it_reco_jets;
      if (Verbosity() >= 1)
        std::cout
            << "QAG4SimulationJet::process_event - Process Reco jet spectrum "
            << reco_jet << std::endl;
      process_Spectrum(topNode, reco_jet, true);
    }
  }

  if (flag(kProcessTruthMatching))
  {
    for (std::set<std::string>::const_iterator it_reco_jets = _reco_jets.begin();
         it_reco_jets != _reco_jets.end(); ++it_reco_jets)
    {
      const std::string& reco_jet = *it_reco_jets;
      if (Verbosity() >= 1)
        std::cout
            << "QAG4SimulationJet::process_event - Process Reco jet spectrum "
            << reco_jet << std::endl;
      process_TruthMatching(topNode, reco_jet);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//! set eta range
void QAG4SimulationJet::set_eta_range(double low, double high)
{
  if (low > high)
    std::swap(low, high);
  assert(low < high);  // eliminate zero range

  eta_range.first = low;
  eta_range.second = high;
}

//! string description of eta range
//! @return TString as ROOT likes
TString
QAG4SimulationJet::get_eta_range_str(const char* eta_name) const
{
  assert(eta_name);
  return TString(
      Form("%.1f < %s < %.1f", eta_range.first, eta_name, eta_range.second));
}

//! acceptance cut on jet object
bool QAG4SimulationJet::jet_acceptance_cut(const Jet* jet) const
{
  assert(jet);
  bool eta_cut = (jet->get_eta() >= eta_range.first) and (jet->get_eta() <= eta_range.second);
  return eta_cut;
}

std::string
QAG4SimulationJet::get_histo_prefix(const std::string& src_jet_name,
                                    const std::string& reco_jet_name)
{
  std::string histo_prefix = "h_QAG4SimJet_";

  if (src_jet_name.length() > 0)
  {
    histo_prefix += src_jet_name;
    histo_prefix += "_";
  }
  if (reco_jet_name.length() > 0)
  {
    histo_prefix += reco_jet_name;
    histo_prefix += "_";
  }

  return histo_prefix;
}

int QAG4SimulationJet::Init_Spectrum(PHCompositeNode* /*topNode*/,
                                     const std::string& jet_name)
{
  Fun4AllHistoManager* hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // normalization plot with counts
  const int norm_size = 3;

  TH1D* h_norm = new TH1D(
      TString(get_histo_prefix(jet_name)) + "Normalization",
      " Normalization;Item;Count", norm_size, .5, norm_size + .5);
  int i = 1;

  h_norm->GetXaxis()->SetBinLabel(i++, "Event");
  h_norm->GetXaxis()->SetBinLabel(i++, "Inclusive Jets");
  h_norm->GetXaxis()->SetBinLabel(i++, "Leading Jets");
  assert(norm_size >= i - 1);
  h_norm->GetXaxis()->LabelsOption("v");
  hm->registerHisto(h_norm);

  hm->registerHisto(
      new TH1F(
          //
          TString(get_histo_prefix(jet_name)) + "Leading_Et",  //
          TString(jet_name) + " leading jet Et, " + get_eta_range_str() + ";E_{T} (GeV)", 100, 0, 100));
  hm->registerHisto(
      new TH1F(
          //
          TString(get_histo_prefix(jet_name)) + "Leading_eta",  //
          TString(jet_name) + " leading jet #eta, " + get_eta_range_str() + ";#eta", 50, -1, 1));
  hm->registerHisto(
      new TH1F(
          //
          TString(get_histo_prefix(jet_name)) + "Leading_phi",  //
          TString(jet_name) + " leading jet #phi, " + get_eta_range_str() + ";#phi", 50, -M_PI, M_PI));

  TH1F* lcomp = new TH1F(
      //
      TString(get_histo_prefix(jet_name)) + "Leading_CompSize",  //
      TString(jet_name) + " leading jet # of component, " + get_eta_range_str() + ";Number of component;", 100, 1, 3000);
  QAHistManagerDef::useLogBins(lcomp->GetXaxis());
  hm->registerHisto(lcomp);

  hm->registerHisto(
      new TH1F(
          //
          TString(get_histo_prefix(jet_name)) + "Leading_Mass",  //
          TString(jet_name) + " leading jet mass, " + get_eta_range_str() + ";Jet Mass (GeV);", 100, 0, 20));
  hm->registerHisto(
      new TH1F(
          //
          TString(get_histo_prefix(jet_name)) + "Leading_CEMC_Ratio",  //
          TString(jet_name) + " leading jet EMCal ratio, " + get_eta_range_str() + ";Energy ratio CEMC/Total;", 100, 0, 1.01));
  hm->registerHisto(
      new TH1F(
          //
          TString(get_histo_prefix(jet_name)) + "Leading_CEMC_HCalIN_Ratio",  //
          TString(jet_name) + " leading jet EMCal+HCal_{IN} ratio, " + get_eta_range_str() + ";Energy ratio (CEMC + HCALIN)/Total;",
          100, 0, 1.01));

  // reco jet has no definition for leakages, since leakage is not reconstructed as part of jet energy.
  // It is only available for truth jets
  hm->registerHisto(
      new TH1F(
          //
          TString(get_histo_prefix(jet_name)) + "Leading_Leakage_Ratio",  //
          TString(jet_name) + " leading jet leakage ratio, " + get_eta_range_str() + ";Energy ratio, Back leakage/Total;", 100,
          0, 1.01));

  TH1F* h = new TH1F(
      //
      TString(get_histo_prefix(jet_name)) + "Inclusive_E",  //
      TString(jet_name) + " inclusive jet E, " + get_eta_range_str() + ";Total jet energy (GeV)", 100, 1e-3, 100);
  QAHistManagerDef::useLogBins(h->GetXaxis());
  hm->registerHisto(h);
  hm->registerHisto(
      new TH1F(
          //
          TString(get_histo_prefix(jet_name)) + "Inclusive_eta",  //
          TString(jet_name) + " inclusive jet #eta, " + get_eta_range_str() + ";#eta;Jet energy density", 50, -1, 1));
  hm->registerHisto(
      new TH1F(
          //
          TString(get_histo_prefix(jet_name)) + "Inclusive_phi",  //
          TString(jet_name) + " inclusive jet #phi, " + get_eta_range_str() + ";#phi;Jet energy density", 50, -M_PI, M_PI));

  return Fun4AllReturnCodes::EVENT_OK;
}

int QAG4SimulationJet::process_Spectrum(PHCompositeNode* topNode,
                                        const std::string& jet_name, const bool is_reco_jet)
{
  JetMap* jets = findNode::getClass<JetMap>(topNode, jet_name.c_str());
  if (!jets)
  {
    std::cout
        << "QAG4SimulationJet::process_Spectrum - Error can not find DST JetMap node "
        << jet_name << std::endl;
    exit(-1);
  }

  Fun4AllHistoManager* hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  TH1D* h_norm = dynamic_cast<TH1D*>(hm->getHisto(
      get_histo_prefix(jet_name) + "Normalization"));
  assert(h_norm);
  h_norm->Fill("Event", 1);
  h_norm->Fill("Inclusive Jets", jets->size());

  TH1F* ie = dynamic_cast<TH1F*>(hm->getHisto(
      (get_histo_prefix(jet_name)) + "Inclusive_E"  //
      ));
  assert(ie);
  TH1F* ieta = dynamic_cast<TH1F*>(hm->getHisto(
      (get_histo_prefix(jet_name)) + "Inclusive_eta"  //
      ));
  assert(ieta);
  TH1F* iphi = dynamic_cast<TH1F*>(hm->getHisto(
      (get_histo_prefix(jet_name)) + "Inclusive_phi"  //
      ));
  assert(iphi);

  Jet* leading_jet = nullptr;
  double max_et = 0;
  for (JetMap::Iter iter = jets->begin(); iter != jets->end(); ++iter)
  {
    Jet* jet = iter->second;
    assert(jet);

    if (not jet_acceptance_cut(jet))
      continue;

    if (jet->get_et() > max_et)
    {
      leading_jet = jet;
      max_et = jet->get_et();
    }

    ie->Fill(jet->get_e());
    ieta->Fill(jet->get_eta());
    iphi->Fill(jet->get_phi());
  }

  if (leading_jet)
  {
    if (Verbosity())
    {
      std::cout
          << "QAG4SimulationJet::process_Spectrum - processing leading jet with # comp = "
          << leading_jet->size_comp() << std::endl;
      leading_jet->identify();
    }

    h_norm->Fill("Leading Jets", 1);

    TH1F* let = dynamic_cast<TH1F*>(hm->getHisto(
        (get_histo_prefix(jet_name)) + "Leading_Et"  //
        ));
    assert(let);
    TH1F* leta = dynamic_cast<TH1F*>(hm->getHisto(
        (get_histo_prefix(jet_name)) + "Leading_eta"  //
        ));
    assert(leta);
    TH1F* lphi = dynamic_cast<TH1F*>(hm->getHisto(
        (get_histo_prefix(jet_name)) + "Leading_phi"  //
        ));
    assert(lphi);

    TH1F* lcomp = dynamic_cast<TH1F*>(hm->getHisto(
        (get_histo_prefix(jet_name)) + "Leading_CompSize"  //
        ));
    assert(lcomp);
    TH1F* lmass = dynamic_cast<TH1F*>(hm->getHisto(
        (get_histo_prefix(jet_name)) + "Leading_Mass"  //
        ));
    assert(lmass);
    TH1F* lcemcr = dynamic_cast<TH1F*>(hm->getHisto(
        (get_histo_prefix(jet_name)) + "Leading_CEMC_Ratio"  //
        ));
    assert(lcemcr);
    TH1F* lemchcalr = dynamic_cast<TH1F*>(hm->getHisto(
        (get_histo_prefix(jet_name)) + "Leading_CEMC_HCalIN_Ratio"  //
        ));
    assert(lemchcalr);
    TH1F* lleak = dynamic_cast<TH1F*>(hm->getHisto(
        (get_histo_prefix(jet_name)) + "Leading_Leakage_Ratio"  //
        ));
    assert(lleak);

    let->Fill(leading_jet->get_et());
    leta->Fill(leading_jet->get_eta());
    lphi->Fill(leading_jet->get_phi());
    lcomp->Fill(leading_jet->size_comp());
    lmass->Fill(leading_jet->get_mass());

    if (is_reco_jet)
    {  // this is a reco jet
      jetevalstacks_map::iterator it_stack = _jetevalstacks.find(jet_name);
      assert(it_stack != _jetevalstacks.end());
      std::shared_ptr<JetEvalStack> eval_stack = it_stack->second;
      assert(eval_stack);
      JetRecoEval* recoeval = eval_stack->get_reco_eval();
      assert(recoeval);

      if (Verbosity() >= VERBOSITY_A_LOT)
      {
        std::cout << __PRETTY_FUNCTION__ << "Leading Jet " << jet_name << ": ";
        //            leading_jet->identify();
        std::cout << "CEMC_TOWER sum = " << recoeval->get_energy_contribution(leading_jet, Jet::CEMC_TOWER) << std::endl;
        std::cout << "CEMC_CLUSTER sum = " << recoeval->get_energy_contribution(leading_jet, Jet::CEMC_CLUSTER) << std::endl;
        std::cout << "HCALIN_TOWER sum = " << recoeval->get_energy_contribution(leading_jet, Jet::HCALIN_TOWER) << std::endl;
        std::cout << "HCALIN_CLUSTER sum = " << recoeval->get_energy_contribution(leading_jet, Jet::HCALIN_CLUSTER) << std::endl;
        std::cout << "leading_jet->get_e() = " << leading_jet->get_e() << std::endl;
      }

      lcemcr->Fill(                                                         //
          (recoeval->get_energy_contribution(leading_jet, Jet::CEMC_TOWER)  //
           +                                                                //
           recoeval->get_energy_contribution(leading_jet,
                                             Jet::CEMC_CLUSTER)  //
           ) /
          leading_jet->get_e());
      lemchcalr->Fill(                                                      //
          (recoeval->get_energy_contribution(leading_jet, Jet::CEMC_TOWER)  //
           +                                                                //
           recoeval->get_energy_contribution(leading_jet,
                                             Jet::CEMC_CLUSTER)  //
           +                                                     //
           recoeval->get_energy_contribution(leading_jet,
                                             Jet::HCALIN_TOWER)  //
           +                                                     //
           recoeval->get_energy_contribution(leading_jet,
                                             Jet::HCALIN_CLUSTER)  //
           ) /
          leading_jet->get_e());
      // reco jet has no definition for leakages, since leakage is not reconstructed as part of jet energy.
      // It is only available for truth jets
    }
    else
    {  // this is a truth jet
      assert(_jettrutheval);

      double cemc_e = 0;
      double hcalin_e = 0;
      double bh_e = 0;

      std::set<PHG4Shower*> showers = _jettrutheval->all_truth_showers(
          leading_jet);

      for (std::set<PHG4Shower*>::const_iterator it = showers.begin();
           it != showers.end(); ++it)
      {
        if (Verbosity() >= VERBOSITY_A_LOT)
        {
          std::cout << __PRETTY_FUNCTION__ << "Leading Truth Jet shower : ";
          (*it)->identify();
        }

        cemc_e += (*it)->get_edep(
            PHG4HitDefs::get_volume_id("G4HIT_CEMC"));
        cemc_e += (*it)->get_edep(
            PHG4HitDefs::get_volume_id("G4HIT_CEMC_ELECTRONICS"));
        cemc_e += (*it)->get_edep(
            PHG4HitDefs::get_volume_id("G4HIT_ABSORBER_CEMC"));

        hcalin_e += (*it)->get_edep(
            PHG4HitDefs::get_volume_id("G4HIT_HCALIN"));
        hcalin_e += (*it)->get_edep(
            PHG4HitDefs::get_volume_id("G4HIT_ABSORBER_HCALIN"));

        bh_e += (*it)->get_edep(PHG4HitDefs::get_volume_id("G4HIT_BH_1"));
        bh_e += (*it)->get_edep(
            PHG4HitDefs::get_volume_id("G4HIT_BH_FORWARD_PLUS"));
        bh_e += (*it)->get_edep(
            PHG4HitDefs::get_volume_id("G4HIT_BH_FORWARD_NEG"));

        if (Verbosity() >= VERBOSITY_A_LOT)
        {
          //            leading_jet->identify();
          std::cout << "Shower cemc_e sum = "
                    << (*it)->get_edep(PHG4HitDefs::get_volume_id("G4HIT_CEMC"))
                    << " + "
                    << (*it)->get_edep(
                           PHG4HitDefs::get_volume_id("G4HIT_CEMC_ELECTRONICS"))
                    << " + "
                    << (*it)->get_edep(
                           PHG4HitDefs::get_volume_id("G4HIT_ABSORBER_CEMC"))
                    << std::endl;
          std::cout << "Shower hcalin_e sum = "
                    << (*it)->get_edep(
                           PHG4HitDefs::get_volume_id("G4HIT_HCALIN"))
                    << " + "
                    << (*it)->get_edep(
                           PHG4HitDefs::get_volume_id("G4HIT_ABSORBER_HCALIN"))
                    << std::endl;
          std::cout << "Shower bh_e sum = "
                    << (*it)->get_edep(PHG4HitDefs::get_volume_id("G4HIT_BH_1"))
                    << " + "
                    << (*it)->get_edep(
                           PHG4HitDefs::get_volume_id("G4HIT_BH_FORWARD_PLUS"))
                    << " + "
                    << (*it)->get_edep(
                           PHG4HitDefs::get_volume_id("G4HIT_BH_FORWARD_NEG"))
                    << std::endl;
        }
      }

      if (Verbosity() >= VERBOSITY_A_LOT)
      {
        std::cout << "cemc_e sum = " << cemc_e << std::endl;
        std::cout << "hcalin_e sum = " << hcalin_e << std::endl;
        std::cout << "leading_jet->get_e() = " << leading_jet->get_e() << std::endl;
      }

      lcemcr->Fill(cemc_e / leading_jet->get_e());
      lemchcalr->Fill((cemc_e + hcalin_e) / leading_jet->get_e());
      lleak->Fill(bh_e / leading_jet->get_e());
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int QAG4SimulationJet::Init_TruthMatching(PHCompositeNode* /*topNode*/,
                                          const std::string& reco_jet_name)
{
  Fun4AllHistoManager* hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  TH2F* h = new TH2F(
      TString(get_histo_prefix(_truth_jet, reco_jet_name)) + "Matching_Count_Truth_Et",  //
      TString(reco_jet_name) + " Matching Count, " + get_eta_range_str() + ";E_{T, Truth} (GeV)", 20, 0, 100, 3, 0.5, 3.5);
  h->GetYaxis()->SetBinLabel(1, "Total");
  h->GetYaxis()->SetBinLabel(2, "Matched");
  h->GetYaxis()->SetBinLabel(3, "Unique Matched");
  hm->registerHisto(h);

  h = new TH2F(
      //
      TString(get_histo_prefix(_truth_jet, reco_jet_name)) + "Matching_Count_Reco_Et",  //
      TString(reco_jet_name) + " Matching Count, " + get_eta_range_str() + ";E_{T, Reco} (GeV)", 20, 0, 100, 3, 0.5, 3.5);
  h->GetYaxis()->SetBinLabel(1, "Total");
  h->GetYaxis()->SetBinLabel(2, "Matched");
  h->GetYaxis()->SetBinLabel(3, "Unique Matched");
  hm->registerHisto(h);

  h = new TH2F(
      TString(get_histo_prefix(_truth_jet, reco_jet_name)) + "Matching_dEt",  //
      TString(reco_jet_name) + " E_{T} difference, " + get_eta_range_str() + ";E_{T, Truth} (GeV);E_{T, Reco} / E_{T, Truth}", 20, 0, 100, 100,
      0, 2);
  hm->registerHisto(h);

  h = new TH2F(
      TString(get_histo_prefix(_truth_jet, reco_jet_name)) + "Matching_dE",  //
      TString(reco_jet_name) + " Jet Energy Difference, " + get_eta_range_str() + ";E_{Truth} (GeV);E_{Reco} / E_{Truth}", 20, 0, 100, 100, 0, 2);
  hm->registerHisto(h);

  h = new TH2F(
      TString(get_histo_prefix(_truth_jet, reco_jet_name)) + "Matching_dEta",  //
      TString(reco_jet_name) + " #eta difference, " + get_eta_range_str() + ";E_{T, Truth} (GeV);#eta_{Reco} - #eta_{Truth}", 20, 0, 100, 200,
      -.1, .1);
  hm->registerHisto(h);

  h = new TH2F(
      TString(get_histo_prefix(_truth_jet, reco_jet_name)) + "Matching_dPhi",  //
      TString(reco_jet_name) + " #phi difference, " + get_eta_range_str() + ";E_{T, Truth} (GeV);#phi_{Reco} - #phi_{Truth} (rad)", 20, 0, 100,
      200, -.1, .1);
  hm->registerHisto(h);

  return Fun4AllReturnCodes::EVENT_OK;
}

int QAG4SimulationJet::process_TruthMatching(PHCompositeNode* topNode,
                                             const std::string& reco_jet_name)
{
  assert(_jet_match_dPhi > 0);
  assert(_jet_match_dEta > 0);
  assert(_jet_match_dE_Ratio > 0);

  Fun4AllHistoManager* hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  TH2F* Matching_Count_Truth_Et = dynamic_cast<TH2F*>(hm->getHisto(
      (get_histo_prefix(_truth_jet, reco_jet_name)) + "Matching_Count_Truth_Et"  //
      ));
  assert(Matching_Count_Truth_Et);
  TH2F* Matching_Count_Reco_Et = dynamic_cast<TH2F*>(hm->getHisto(
      (get_histo_prefix(_truth_jet, reco_jet_name)) + "Matching_Count_Reco_Et"  //
      ));
  assert(Matching_Count_Reco_Et);
  TH2F* Matching_dEt = dynamic_cast<TH2F*>(hm->getHisto(
      (get_histo_prefix(_truth_jet, reco_jet_name)) + "Matching_dEt"  //
      ));
  assert(Matching_dEt);
  TH2F* Matching_dE = dynamic_cast<TH2F*>(hm->getHisto(
      (get_histo_prefix(_truth_jet, reco_jet_name)) + "Matching_dE"  //
      ));
  assert(Matching_dE);
  TH2F* Matching_dEta = dynamic_cast<TH2F*>(hm->getHisto(
      (get_histo_prefix(_truth_jet, reco_jet_name)) + "Matching_dEta"  //
      ));
  assert(Matching_dEta);
  TH2F* Matching_dPhi = dynamic_cast<TH2F*>(hm->getHisto(
      (get_histo_prefix(_truth_jet, reco_jet_name)) + "Matching_dPhi"  //
      ));
  assert(Matching_dPhi);

  jetevalstacks_map::iterator it_stack = _jetevalstacks.find(reco_jet_name);
  assert(it_stack != _jetevalstacks.end());
  std::shared_ptr<JetEvalStack> eval_stack = it_stack->second;
  assert(eval_stack);
  JetRecoEval* recoeval = eval_stack->get_reco_eval();
  assert(recoeval);

  // iterate over truth jets
  JetMap* truthjets = findNode::getClass<JetMap>(topNode, _truth_jet);
  if (!truthjets)
  {
    std::cout
        << "QAG4SimulationJet::process_TruthMatching - Error can not find DST JetMap node "
        << _truth_jet << std::endl;
    exit(-1);
  }

  // search for leading truth
  Jet* truthjet = nullptr;
  double max_et = 0;
  for (JetMap::Iter iter = truthjets->begin(); iter != truthjets->end(); ++iter)
  {
    Jet* jet = iter->second;
    assert(jet);

    if (not jet_acceptance_cut(jet))
      continue;

    if (jet->get_et() > max_et)
    {
      truthjet = jet;
      max_et = jet->get_et();
    }
  }

  // match leading truth
  if (truthjet)
  {
    if (Verbosity() > 1)
    {
      std::cout << "QAG4SimulationJet::process_TruthMatching - " << _truth_jet
                << " process truth jet ";
      truthjet->identify();
    }

    Matching_Count_Truth_Et->Fill(truthjet->get_et(), "Total", 1);

    {  // inclusive best energy match

      const Jet* recojet = recoeval->best_jet_from(truthjet);
      if (Verbosity() > 1)
      {
        std::cout << "QAG4SimulationJet::process_TruthMatching - " << _truth_jet
                  << " inclusively matched with best reco jet: ";
        recojet->identify();
      }

      if (recojet)
      {
        const double dPhi = recojet->get_phi() - truthjet->get_phi();
        Matching_dPhi->Fill(truthjet->get_et(), dPhi);

        if (fabs(dPhi) < _jet_match_dPhi)
        {
          const double dEta = recojet->get_eta() - truthjet->get_eta();
          Matching_dEta->Fill(truthjet->get_et(), dEta);

          if (fabs(dEta) < _jet_match_dEta)
          {
            const double Et_r = recojet->get_et() / (truthjet->get_et() + 1e-9);
            const double E_r = recojet->get_e() / (truthjet->get_e() + 1e-9);
            Matching_dEt->Fill(truthjet->get_et(), Et_r);
            Matching_dE->Fill(truthjet->get_et(), E_r);

            if (fabs(E_r - 1) < _jet_match_dE_Ratio)
            {
              // matched in eta, phi and energy

              Matching_Count_Truth_Et->Fill(truthjet->get_et(),
                                            "Matched", 1);
            }

          }  //  if (fabs(dEta) < 0.1)

        }  // if (fabs(dPhi) < 0.1)

      }  //       if (recojet)
    }    // inclusive best energy match
    {    // unique match

      const Jet* recojet = recoeval->unique_reco_jet_from_truth(truthjet);
      if (recojet)
      {
        if (Verbosity() > 1)
        {
          std::cout << "QAG4SimulationJet::process_TruthMatching - " << _truth_jet
                    << " uniquely matched with reco jet: ";
          recojet->identify();
        }

        const double dPhi = recojet->get_phi() - truthjet->get_phi();

        if (fabs(dPhi) < _jet_match_dPhi)
        {
          const double dEta = recojet->get_eta() - truthjet->get_eta();

          if (fabs(dEta) < _jet_match_dEta)
          {
            const double E_r = recojet->get_e() / (truthjet->get_e() + 1e-9);
            if (fabs(E_r - 1) < _jet_match_dE_Ratio)
            {
              // matched in eta, phi and energy

              Matching_Count_Truth_Et->Fill(truthjet->get_et(),
                                            "Unique Matched", 1);
            }

          }  //  if (fabs(dEta) < 0.1)

        }  // if (fabs(dPhi) < 0.1)

      }  //       if (recojet)
    }    // unique match

  }  //  if (truthjet)

  // next for reco jets
  JetMap* recojets = findNode::getClass<JetMap>(topNode, reco_jet_name);
  if (!recojets)
  {
    std::cout
        << "QAG4SimulationJet::process_TruthMatching - Error can not find DST JetMap node "
        << reco_jet_name << std::endl;
    exit(-1);
  }

  // search for leading reco jet
  Jet* recojet = nullptr;
  max_et = 0;
  for (JetMap::Iter iter = recojets->begin(); iter != recojets->end(); ++iter)
  {
    Jet* jet = iter->second;
    assert(jet);

    if (not jet_acceptance_cut(jet))
      continue;

    if (jet->get_et() > max_et)
    {
      recojet = jet;
      max_et = jet->get_et();
    }
  }

  // match leading reco jet
  if (recojet)
  {
    if (Verbosity() > 1)
    {
      std::cout << "QAG4SimulationJet::process_TruthMatching - " << reco_jet_name
                << " process reco jet ";
      recojet->identify();
    }

    Matching_Count_Reco_Et->Fill(recojet->get_et(), "Total", 1);

    {  // inclusive best energy match
      Jet* truthjet = recoeval->max_truth_jet_by_energy(recojet);
      if (truthjet)
      {
        const double dPhi = recojet->get_phi() - truthjet->get_phi();
        if (fabs(dPhi) < _jet_match_dPhi)
        {
          const double dEta = recojet->get_eta() - truthjet->get_eta();
          if (fabs(dEta) < _jet_match_dEta)
          {
            const double E_r = recojet->get_e() / (truthjet->get_e() + 1e-9);

            if (fabs(E_r - 1) < _jet_match_dE_Ratio)
            {
              // matched in eta, phi and energy

              Matching_Count_Reco_Et->Fill(recojet->get_et(),
                                           "Unique Matched", 1);
            }

          }  //  if (fabs(dEta) < 0.1)

        }  // if (fabs(dPhi) < 0.1)

      }  //      if (truthjet)
    }    // inclusive best energy match

    {  // unique match
      Jet* truthjet = recoeval->unique_truth_jet_from_reco(recojet);
      if (truthjet)
      {
        const double dPhi = recojet->get_phi() - truthjet->get_phi();
        if (fabs(dPhi) < _jet_match_dPhi)
        {
          const double dEta = recojet->get_eta() - truthjet->get_eta();
          if (fabs(dEta) < _jet_match_dEta)
          {
            const double E_r = recojet->get_e() / (truthjet->get_e() + 1e-9);

            if (fabs(E_r - 1) < _jet_match_dE_Ratio)
            {
              // matched in eta, phi and energy

              Matching_Count_Reco_Et->Fill(recojet->get_et(),
                                           "Matched", 1);
            }

          }  //  if (fabs(dEta) < 0.1)

        }  // if (fabs(dPhi) < 0.1)

      }  //      if (truthjet)
    }    // unique match

  }  // if (recojet)

  return Fun4AllReturnCodes::EVENT_OK;
}
