#include "QAG4SimulationJet.h"

#include <fun4all/SubsysReco.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/PHTFileServer.h>
#include <phool/PHCompositeNode.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>

#include <phool/PHCompositeNode.h>

#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4Particle.h>

#include <g4eval/JetEvalStack.h>

#include <TString.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TVector3.h>
#include <TLorentzVector.h>

#include <exception>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <cassert>
#include <cmath>

#include "QAHistManagerDef.h"

using namespace std;

QAG4SimulationJet::QAG4SimulationJet(const std::string & truth_jet,
    enu_flags flags) :
    SubsysReco("QAG4SimulationJet_" + truth_jet), //
    _jetevalstacks(), //
    _truth_jet(truth_jet), _reco_jets(), _flags(flags), //
    eta_range(-1, 1), //
    _jet_match_dEta(.1), _jet_match_dPhi(.1), _jet_match_dE_Ratio(.5)
{

}

QAG4SimulationJet::~QAG4SimulationJet()
{
}

int
QAG4SimulationJet::InitRun(PHCompositeNode *topNode)
{

  if (flag(kProcessTruthMatching))
    {
      for (set<string>::const_iterator it_reco_jets = _reco_jets.begin();
          it_reco_jets != _reco_jets.end(); ++it_reco_jets)
        {
          const string & reco_jet = *it_reco_jets;

          jetevalstacks_map::iterator it_jetevalstack = _jetevalstacks.find(
              reco_jet);

          if (it_jetevalstack == _jetevalstacks.end())
            {
              _jetevalstacks[reco_jet] = shared_ptr < JetEvalStack
                  > (new JetEvalStack(topNode, reco_jet, _truth_jet));
              assert(_jetevalstacks[reco_jet]);
              _jetevalstacks[reco_jet]->set_strict(true);
              _jetevalstacks[reco_jet]->set_verbosity(verbosity + 1);
            }
          else
            {
              assert(it_jetevalstack->second);
            }
        }
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

int
QAG4SimulationJet::End(PHCompositeNode *topNode)
{

  return Fun4AllReturnCodes::EVENT_OK;
}

int
QAG4SimulationJet::Init(PHCompositeNode *topNode)
{

  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  // normalization plot with counts
  const int norm_size = 2 + _reco_jets.size();

  TH1D * h = new TH1D(TString(get_histo_prefix(_truth_jet)) + "Normalization",
      " Normalization;Z (cm);Count", norm_size, .5, norm_size + .5);
  int i = 1;

  h->GetXaxis()->SetBinLabel(i++, "Event");
  h->GetXaxis()->SetBinLabel(i++, _truth_jet.c_str());

  for (set<string>::const_iterator it_reco_jets = _reco_jets.begin();
      it_reco_jets != _reco_jets.end(); ++it_reco_jets)
    {
      const string & reco_jet = *it_reco_jets;
      h->GetXaxis()->SetBinLabel(i++, reco_jet.c_str());
    }

  assert(norm_size >= i - 1);
  h->GetXaxis()->LabelsOption("v");
  hm->registerHisto(h);

  if (flag(kProcessTruthSpectrum))
    {
      if (verbosity >= 1)
        cout << "QAG4SimulationJet::Init - Process TruthSpectrum " << _truth_jet
            << endl;
      Init_Spectrum(topNode, _truth_jet);
    }

  if (flag(kProcessRecoSpectrum))
    {
      for (set<string>::const_iterator it_reco_jets = _reco_jets.begin();
          it_reco_jets != _reco_jets.end(); ++it_reco_jets)
        {
          const string & reco_jet = *it_reco_jets;
          if (verbosity >= 1)
            cout << "QAG4SimulationJet::Init - Process Reco jet spectrum "
                << reco_jet << endl;
          Init_Spectrum(topNode, reco_jet);
        }
    }

  if (flag(kProcessTruthMatching))
    {
      for (set<string>::const_iterator it_reco_jets = _reco_jets.begin();
          it_reco_jets != _reco_jets.end(); ++it_reco_jets)
        {
          const string & reco_jet = *it_reco_jets;
          if (verbosity >= 1)
            cout << "QAG4SimulationJet::Init - Process Reco jet spectrum "
                << reco_jet << endl;
          Init_TruthMatching(topNode, reco_jet);
        }
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

int
QAG4SimulationJet::process_event(PHCompositeNode *topNode)
{

  if (verbosity > 2)
    cout << "QAG4SimulationJet::process_event() entered" << endl;

  for (jetevalstacks_map::iterator it_jetevalstack = _jetevalstacks.begin();
      it_jetevalstack != _jetevalstacks.end(); ++it_jetevalstack)
    {
      assert(it_jetevalstack->second);
      it_jetevalstack->second->next_event(topNode);
    }

  if (flag(kProcessTruthSpectrum))
    {
      if (verbosity >= 1)
        cout << "QAG4SimulationJet::process_event - Process TruthSpectrum "
            << _truth_jet << endl;
      process_Spectrum(topNode, _truth_jet);
    }

  if (flag(kProcessRecoSpectrum))
    {
      for (set<string>::const_iterator it_reco_jets = _reco_jets.begin();
          it_reco_jets != _reco_jets.end(); ++it_reco_jets)
        {
          const string & reco_jet = *it_reco_jets;
          if (verbosity >= 1)
            cout
                << "QAG4SimulationJet::process_event - Process Reco jet spectrum "
                << reco_jet << endl;
          process_Spectrum(topNode, reco_jet);
        }
    }

  if (flag(kProcessTruthMatching))
    {
      for (set<string>::const_iterator it_reco_jets = _reco_jets.begin();
          it_reco_jets != _reco_jets.end(); ++it_reco_jets)
        {
          const string & reco_jet = *it_reco_jets;
          if (verbosity >= 1)
            cout
                << "QAG4SimulationJet::process_event - Process Reco jet spectrum "
                << reco_jet << endl;
          process_TruthMatching(topNode, reco_jet);
        }
    }

  // at the end, count success events
  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);
  TH1D* h_norm = dynamic_cast<TH1D*>(hm->getHisto(
      get_histo_prefix(_truth_jet) + "Normalization"));
  assert(h_norm);
  h_norm->Fill("Event", 1);

  return Fun4AllReturnCodes::EVENT_OK;
}

//! set eta range
void
QAG4SimulationJet::set_eta_range(double low, double high)
{
  if (low > high)
    swap(low, high);
  assert(low < high); // eliminate zero range

  eta_range.first = low;
  eta_range.second = high;
}

//! string description of eta range
//! @return TString as ROOT likes
TString
QAG4SimulationJet::get_eta_range_str(const char * eta_name) const
{
  assert(eta_name);
  return TString(
      Form("%.1f < %s < %.1f", eta_range.first, eta_name, eta_range.second));
}

//! acceptance cut on jet object
bool
QAG4SimulationJet::jet_acceptance_cut(const Jet * jet) const
{
  assert(jet);
  bool eta_cut = (jet->get_eta() >= eta_range.first)
      and (jet->get_eta() <= eta_range.second);
  return eta_cut;
}

string
QAG4SimulationJet::get_histo_prefix(const std::string & src_jet_name,
    const std::string & reco_jet_name)
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

int
QAG4SimulationJet::Init_Spectrum(PHCompositeNode *topNode,
    const std::string & jet_name)
{

  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);
  hm->registerHisto(
      new TH1F(
          //
          TString(get_histo_prefix(jet_name)) + "Leading_Et", //
          TString(jet_name) + " leading jet Et, " + get_eta_range_str()
              + ";E_{T} (GeV)", 100, 0, 100));
  hm->registerHisto(
      new TH1F(
          //
          TString(get_histo_prefix(jet_name)) + "Leading_eta", //
          TString(jet_name) + " leading jet #eta, " + get_eta_range_str()
              + ";#eta", 50, -1, 1));
  hm->registerHisto(
      new TH1F(
          //
          TString(get_histo_prefix(jet_name)) + "Leading_phi", //
          TString(jet_name) + " leading jet #phi, " + get_eta_range_str()
              + ";#phi", 50, -M_PI, M_PI));

  TH1F * h = new TH1F(
      //
      TString(get_histo_prefix(jet_name)) + "Inclusive_E", //
      TString(jet_name) + " inclusive jet E, " + get_eta_range_str()
          + ";Total jet energy (GeV)", 100, 1e-3, 100);
  QAHistManagerDef::useLogBins(h->GetXaxis());
  hm->registerHisto(h);
  hm->registerHisto(
      new TH1F(
          //
          TString(get_histo_prefix(jet_name)) + "Inclusive_eta", //
          TString(jet_name) + " inclusive jet #eta, " + get_eta_range_str()
              + ";#eta;Jet energy density", 50, -1, 1));
  hm->registerHisto(
      new TH1F(
          //
          TString(get_histo_prefix(jet_name)) + "Inclusive_phi", //
          TString(jet_name) + " inclusive jet #phi, " + get_eta_range_str()
              + ";#phi;Jet energy density", 50, -M_PI, M_PI));

  return Fun4AllReturnCodes::EVENT_OK;
}

int
QAG4SimulationJet::process_Spectrum(PHCompositeNode *topNode,
    const std::string & jet_name)
{
  JetMap* jets = findNode::getClass<JetMap>(topNode, jet_name.c_str());
  if (!jets)
    {
      cout
          << "QAG4SimulationJet::process_Spectrum - Error can not find DST JetMap node "
          << jet_name << endl;
      exit(-1);
    }

  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  TH1D* h_norm = dynamic_cast<TH1D*>(hm->getHisto(
      get_histo_prefix(_truth_jet) + "Normalization"));
  assert(h_norm);
  h_norm->Fill(jet_name.c_str(), jets->size());

  TH1F * ie = dynamic_cast<TH1F*>(hm->getHisto(
      (get_histo_prefix(jet_name)) + "Inclusive_E" //
          ));
  assert(ie);
  TH1F * ieta = dynamic_cast<TH1F*>(hm->getHisto(
      (get_histo_prefix(jet_name)) + "Inclusive_eta" //
          ));
  assert(ieta);
  TH1F * iphi = dynamic_cast<TH1F*>(hm->getHisto(
      (get_histo_prefix(jet_name)) + "Inclusive_phi" //
          ));
  assert(iphi);

  Jet* leading_jet = NULL;
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

      TH1F * let = dynamic_cast<TH1F*>(hm->getHisto(
          (get_histo_prefix(jet_name)) + "Leading_Et" //
              ));
      assert(let);
      TH1F * leta = dynamic_cast<TH1F*>(hm->getHisto(
          (get_histo_prefix(jet_name)) + "Leading_eta" //
              ));
      assert(leta);
      TH1F * lphi = dynamic_cast<TH1F*>(hm->getHisto(
          (get_histo_prefix(jet_name)) + "Leading_phi" //
              ));
      assert(lphi);

      let->Fill(leading_jet->get_et());
      leta->Fill(leading_jet->get_eta());
      lphi->Fill(leading_jet->get_phi());

    }

  return Fun4AllReturnCodes::EVENT_OK;
}

int
QAG4SimulationJet::Init_TruthMatching(PHCompositeNode *topNode,
    const std::string & reco_jet_name)
{
  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  TH2F * h = new TH2F(
      TString(get_histo_prefix(_truth_jet, reco_jet_name))
          + "Matching_Count_Truth_Et", //
      TString(reco_jet_name) + " Matching Count, " + get_eta_range_str()
          + ";E_{T, Truth} (GeV)", 20, 0, 100, 2, 0.5, 2.5);
  h->GetYaxis()->SetBinLabel(1, "Total");
  h->GetYaxis()->SetBinLabel(2, "Matched");
  hm->registerHisto(h);

  h = new TH2F(
      //
      TString(get_histo_prefix(_truth_jet, reco_jet_name))
          + "Matching_Count_Reco_Et", //
      TString(reco_jet_name) + " Matching Count, " + get_eta_range_str()
          + ";E_{T, Reco} (GeV)", 20, 0, 100, 2, 0.5, 2.5);
  h->GetYaxis()->SetBinLabel(1, "Total");
  h->GetYaxis()->SetBinLabel(2, "Matched");
  hm->registerHisto(h);

  h = new TH2F(
      TString(get_histo_prefix(_truth_jet, reco_jet_name)) + "Matching_dEt", //
      TString(reco_jet_name) + " E_{T} difference, " + get_eta_range_str()
          + ";E_{T, Truth} (GeV);E_{T, Reco} / E_{T, Truth}", 20, 0, 100,
      100, 0, 2);
  hm->registerHisto(h);

  h = new TH2F(
      TString(get_histo_prefix(_truth_jet, reco_jet_name)) + "Matching_dE", //
      TString(reco_jet_name) + " Jet Energy Difference, " + get_eta_range_str()
          + ";E_{Truth} (GeV);E_{Reco} / E_{Truth}", 20, 0, 100, 100, 0,
      2);
  hm->registerHisto(h);

  h = new TH2F(
      TString(get_histo_prefix(_truth_jet, reco_jet_name)) + "Matching_dEta", //
      TString(reco_jet_name) + " #eta difference, " + get_eta_range_str()
          + ";E_{T, Truth} (GeV);#eta_{Reco} - #eta_{Truth}", 20, 0, 100,
      200, -.1, .1);
  hm->registerHisto(h);

  h = new TH2F(
      TString(get_histo_prefix(_truth_jet, reco_jet_name)) + "Matching_dPhi", //
      TString(reco_jet_name) + " #phi difference, " + get_eta_range_str()
          + ";E_{T, Truth} (GeV);#phi_{Reco} - #phi_{Truth} (rad)", 20, 0, 100,
      200, -.1, .1);
  hm->registerHisto(h);

  return Fun4AllReturnCodes::EVENT_OK;
}

int
QAG4SimulationJet::process_TruthMatching(PHCompositeNode *topNode,
    const std::string & reco_jet_name)
{
  assert(_jet_match_dPhi > 0);
  assert(_jet_match_dEta > 0);
  assert(_jet_match_dE_Ratio > 0);

  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  TH2F * Matching_Count_Truth_Et = dynamic_cast<TH2F*>(hm->getHisto(
      (get_histo_prefix(_truth_jet, reco_jet_name)) + "Matching_Count_Truth_Et" //
          ));
  assert(Matching_Count_Truth_Et);
  TH2F * Matching_Count_Reco_Et = dynamic_cast<TH2F*>(hm->getHisto(
      (get_histo_prefix(_truth_jet, reco_jet_name)) + "Matching_Count_Reco_Et" //
          ));
  assert(Matching_Count_Reco_Et);
  TH2F * Matching_dEt = dynamic_cast<TH2F*>(hm->getHisto(
      (get_histo_prefix(_truth_jet, reco_jet_name)) + "Matching_dEt" //
          ));
  assert(Matching_dEt);
  TH2F * Matching_dE = dynamic_cast<TH2F*>(hm->getHisto(
      (get_histo_prefix(_truth_jet, reco_jet_name)) + "Matching_dE" //
          ));
  assert(Matching_dE);
  TH2F * Matching_dEta = dynamic_cast<TH2F*>(hm->getHisto(
      (get_histo_prefix(_truth_jet, reco_jet_name)) + "Matching_dEta" //
          ));
  assert(Matching_dEta);
  TH2F * Matching_dPhi = dynamic_cast<TH2F*>(hm->getHisto(
      (get_histo_prefix(_truth_jet, reco_jet_name)) + "Matching_dPhi" //
          ));
  assert(Matching_dPhi);

  jetevalstacks_map::iterator it_stack = _jetevalstacks.find(reco_jet_name);
  assert(it_stack != _jetevalstacks.end());
  shared_ptr<JetEvalStack> eval_stack = it_stack->second;
  assert(eval_stack);
  JetRecoEval* recoeval = eval_stack->get_reco_eval();
  assert(recoeval);

  // iterate over truth jets
  JetMap* truthjets = findNode::getClass<JetMap>(topNode, _truth_jet);
  if (!truthjets)
    {
      cout
          << "QAG4SimulationJet::process_TruthMatching - Error can not find DST JetMap node "
          << _truth_jet << endl;
      exit(-1);
    }
  for (JetMap::Iter iter = truthjets->begin(); iter != truthjets->end(); ++iter)
    {
      Jet* truthjet = iter->second;
      assert(truthjet);

      if (verbosity > 1)
        {
          cout << "QAG4SimulationJet::process_TruthMatching - " << _truth_jet
              << " process truth jet ";
          truthjet->identify();
        }

      if (not jet_acceptance_cut(truthjet))
        continue;

      Matching_Count_Truth_Et->Fill(truthjet->get_et(), "Total", 1);

      const Jet* recojet = recoeval->best_jet_from(truthjet);
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

                  const double Et_r = recojet->get_et()
                      / (truthjet->get_et() + 1e-9);
                  const double E_r = recojet->get_e()
                      / (truthjet->get_e() + 1e-9);
                  Matching_dEt->Fill(truthjet->get_et(), Et_r);
                  Matching_dE->Fill(truthjet->get_et(), E_r);

                  if (fabs(E_r - 1) < _jet_match_dE_Ratio)
                    {
                      // matched in eta, phi and energy

                      Matching_Count_Truth_Et->Fill(truthjet->get_et(),
                          "Matched", 1);
                    }

                } //  if (fabs(dEta) < 0.1)

            } // if (fabs(dPhi) < 0.1)

        } //       if (recojet)

    } //  for (JetMap::Iter iter = truthjets->begin(); iter != truthjets->end(); ++iter)

  // next for reco jets
  JetMap* recojets = findNode::getClass<JetMap>(topNode, reco_jet_name);
  if (!recojets)
    {
      cout
          << "QAG4SimulationJet::process_TruthMatching - Error can not find DST JetMap node "
          << reco_jet_name << endl;
      exit(-1);
    }
  for (JetMap::Iter iter = recojets->begin(); iter != recojets->end(); ++iter)
    {
      Jet* recojet = iter->second;
      assert(recojet);

      if (verbosity > 1)
        {
          cout << "QAG4SimulationJet::process_TruthMatching - " << reco_jet_name
              << " process reco jet ";
          recojet->identify();
        }

      if (not jet_acceptance_cut(recojet))
        continue;

      Matching_Count_Reco_Et->Fill(recojet->get_et(), "Total", 1);

      Jet* truthjet = recoeval->max_truth_jet_by_energy(recojet);
      if (truthjet)
        {

          const double dPhi = recojet->get_phi() - truthjet->get_phi();
          if (fabs(dPhi) < _jet_match_dPhi)
            {

              const double dEta = recojet->get_eta() - truthjet->get_eta();
              if (fabs(dEta) < _jet_match_dEta)
                {

                  const double E_r = recojet->get_e()
                      / (truthjet->get_e() + 1e-9);

                  if (fabs(E_r - 1) < _jet_match_dE_Ratio)
                    {
                      // matched in eta, phi and energy

                      Matching_Count_Reco_Et->Fill(recojet->get_et(), "Matched",
                          1);
                    }

                } //  if (fabs(dEta) < 0.1)

            } // if (fabs(dPhi) < 0.1)

        } //      if (truthjet)

    } //  for (JetMap::Iter iter = recojets->begin(); iter != recojets->end(); ++iter)

  return Fun4AllReturnCodes::EVENT_OK;
}

