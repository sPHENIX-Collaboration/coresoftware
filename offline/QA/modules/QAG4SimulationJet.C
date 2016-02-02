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
    _truth_jet(truth_jet), _reco_jets(), _flags(flags), _ievent(0)
{

}

QAG4SimulationJet::~QAG4SimulationJet()
{
}

int
QAG4SimulationJet::InitRun(PHCompositeNode *topNode)
{
  _ievent = 0;

  if (flag(kProcessComparison))
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
                  > (new JetEvalStack(topNode, _truth_jet, reco_jet));
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

  _ievent = 0;

  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  const int norm_size = 2 + _reco_jets.size();

  TH1D * h = new TH1D(TString(get_histo_prefix(_truth_jet)) + "Normalization",
      " Normalization;Z (cm);Count", norm_size, -.5, norm_size - .5);
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

  if (flag(kProcessComparison))
    {
      for (set<string>::const_iterator it_reco_jets = _reco_jets.begin();
          it_reco_jets != _reco_jets.end(); ++it_reco_jets)
        {
          const string & reco_jet = *it_reco_jets;
          if (verbosity >= 1)
            cout << "QAG4SimulationJet::Init - Process Reco jet spectrum "
                << reco_jet << endl;
          Init_Comparison(topNode, reco_jet);
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

  if (flag(kProcessComparison))
    {
      for (set<string>::const_iterator it_reco_jets = _reco_jets.begin();
          it_reco_jets != _reco_jets.end(); ++it_reco_jets)
        {
          const string & reco_jet = *it_reco_jets;
          if (verbosity >= 1)
            cout
                << "QAG4SimulationJet::process_event - Process Reco jet spectrum "
                << reco_jet << endl;
          process_Comparison(topNode, reco_jet);
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
  hm->registerHisto(new TH1F( //
      TString(get_histo_prefix(jet_name)) + "Leading_Et", //
      TString(jet_name) + " leading jet Et;E_{T} (GeV)", 100, 0, 100));
  hm->registerHisto(new TH1F( //
      TString(get_histo_prefix(jet_name)) + "Leading_eta", //
      TString(jet_name) + " leading jet #eta;#eta", 50, -1, 1));
  hm->registerHisto(new TH1F( //
      TString(get_histo_prefix(jet_name)) + "Leading_phi", //
      TString(jet_name) + " leading jet #phi;#phi", 50, -M_PI, M_PI));

  TH1F * h = new TH1F(
      //
      TString(get_histo_prefix(jet_name)) + "Inclusive_E", //
      TString(jet_name) + " inclusive jet E;Total jet energy (GeV)", 100, 1e-3,
      100);
  QAHistManagerDef::useLogBins(h->GetXaxis());
  hm->registerHisto(h);
  hm->registerHisto(
      new TH1F(
          //
          TString(get_histo_prefix(jet_name)) + "Inclusive_eta", //
          TString(jet_name) + " inclusive jet #eta;#eta;Jet energy density", 50,
          -1, 1));
  hm->registerHisto(
      new TH1F(
          //
          TString(get_histo_prefix(jet_name)) + "Inclusive_phi", //
          TString(jet_name) + " inclusive jet #phi;#phi;Jet energy density", 50,
          -M_PI, M_PI));

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
QAG4SimulationJet::Init_Comparison(PHCompositeNode *topNode,
    const std::string & reco_jet_name)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int
QAG4SimulationJet::process_Comparison(PHCompositeNode *topNode,
    const std::string & reco_jet_name)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

