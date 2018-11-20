// $Id: $

/*!
 * \file JetHepMCLoader.C
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "JetHepMCLoader.h"

#include "JetMapV1.h"
#include "JetV1.h"

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHTypedNodeIterator.h>
#include <phool/getClass.h>

#include <TH1D.h>
#include <TH2F.h>

#include <boost/algorithm/string.hpp>

#include <cassert>
#include <iostream>

using namespace std;

JetHepMCLoader::JetHepMCLoader(const std::string &jetInputCategory)
  : SubsysReco("JetHepMCLoader_" + jetInputCategory)
  , m_jetInputCategory(jetInputCategory)
  , m_saveQAPlots(false)
{
}

JetHepMCLoader::~JetHepMCLoader()
{
}

int JetHepMCLoader::Init(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int JetHepMCLoader::InitRun(PHCompositeNode *topNode)
{
  if (m_saveQAPlots)
  {
    Fun4AllHistoManager *hm = getHistoManager();
    assert(hm);

    const int n_bins = 1 + m_jetSrc.size();

    TH1D *h = new TH1D("hNormalization",  //
                       "Normalization;Items;Summed quantity", n_bins, .5, n_bins + .5);
    int i = 1;
    h->GetXaxis()->SetBinLabel(i++, "Event count");
    for (const hepmc_jet_src &src : m_jetSrc)
    {
      h->GetXaxis()->SetBinLabel(i++, (string("SubEvent ") + src.m_name).c_str());
    }
    h->GetXaxis()->LabelsOption("v");
    hm->registerHisto(h);

    for (const hepmc_jet_src &src : m_jetSrc)
    {
      hm->registerHisto(
          new TH2F((string("hJetEtEta_") + src.m_name).c_str(),  //
                   ("Jet distribution from " + src.m_name + ";Jet #eta;Jet E_{T} [GeV]").c_str(),
                   40, -4., 4.,
                   100, 0, 100));
    }
  }  //   if (m_saveQAPlots)

  return CreateNodes(topNode);
}

int JetHepMCLoader::process_event(PHCompositeNode *topNode)
{
  // For pile-up simulation: define GenEventMap
  PHHepMCGenEventMap *genevtmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");

  if (!genevtmap)
  {
    static bool once = true;

    if (once and Verbosity())
    {
      once = false;

      cout << "HepMCNodeReader::process_event - No PHHepMCGenEventMap node. Do not perform HepMC->Geant4 input" << endl;
    }

    return Fun4AllReturnCodes::DISCARDEVENT;
  }

  if (m_saveQAPlots)
  {
    Fun4AllHistoManager *hm = getHistoManager();
    assert(hm);
    TH1D *h_norm = dynamic_cast<TH1D *>(hm->getHisto("hNormalization"));
    assert(h_norm);
    h_norm->Fill("Event count", 1);
  }  //   if (m_saveQAPlots)

  for (const hepmc_jet_src &src : m_jetSrc)
  {
    JetMap *jets = findNode::getClass<JetMap>(topNode, src.m_name);
    assert(jets);

    jets->set_algo(src.m_algorithmID);
    jets->set_par(src.m_parameter);
    jets->insert_src(Jet::HEPMC_IMPORT);

    PHHepMCGenEvent *genevt =
        genevtmap->get(src.m_embeddingID);

    if (genevt == nullptr) continue;

    HepMC::GenEvent *evt = genevt->getEvent();
    if (!evt)
    {
      cout << PHWHERE << " no evt pointer under HEPMC Node found:";
      genevt->identify();
      return Fun4AllReturnCodes::ABORTEVENT;
    }

    assert(genevt->get_embedding_id() == src.m_embeddingID);

    TH2F *hjet = nullptr;

    if (m_saveQAPlots)
    {
      Fun4AllHistoManager *hm = getHistoManager();
      assert(hm);
      TH1D *h_norm = dynamic_cast<TH1D *>(hm->getHisto("hNormalization"));
      assert(h_norm);
      h_norm->Fill((string("SubEvent ") + src.m_name).c_str(), 1);

      hjet = dynamic_cast<TH2F *>(hm->getHisto(string("hJetEtEta_") + src.m_name));
      assert(hjet);

    }  //   if (m_saveQAPlots)

    const double mom_factor = HepMC::Units::conversion_factor(evt->momentum_unit(), HepMC::Units::GEV);

    for (HepMC::GenEvent::particle_const_iterator p = evt->particles_begin();
         p != evt->particles_end(); ++p)
    {
      HepMC::GenParticle *part = (*p);

      assert(part);
      if (Verbosity() >= VERBOSITY_A_LOT)
        part->print();

      if (part->status() == src.m_tagStatus and part->pdg_id() == src.m_tagPID)
      {
        Jet *jet = new JetV1();

        jet->set_px(part->momentum().px() * mom_factor);
        jet->set_py(part->momentum().py() * mom_factor);
        jet->set_pz(part->momentum().pz() * mom_factor);
        jet->set_e(part->momentum().e() * mom_factor);

        jet->insert_comp(Jet::HEPMC_IMPORT, part->barcode());

        jets->insert(jet);

        if (hjet)
        {
          hjet->Fill(jet->get_eta(), jet->get_et());
        }

      }  //       if (part->status() == src.m_tagStatus and part->pdg_id() == src.m_tagPID)
    }
  }  //  for (const hepmc_jet_src &src : m_jetSrc)

  return Fun4AllReturnCodes::EVENT_OK;
}

int JetHepMCLoader::End(PHCompositeNode *topNode)
{
  if (m_saveQAPlots)
  {
    Fun4AllHistoManager *hm = getHistoManager();
    assert(hm);

    cout << "JetHepMCLoader::End - saving QA histograms to " << Name() + ".root" << endl;
    hm->dumpHistos(Name() + ".root", "RECREATE");
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void JetHepMCLoader::addJet(
    const std::string &name,
    int embeddingID,
    Jet::ALGO algorithm,
    double parameter,
    int tagPID,
    int tagStatus)
{
  string algorithmName = "Undefined_Jet_Algorithm";

  switch (algorithm)
  {
  case Jet::ANTIKT:
    algorithmName = "ANTIKT";
    break;

  case Jet::KT:
    algorithmName = "KT";
    break;

  case Jet::CAMBRIDGE:
    algorithmName = "CAMBRIDGE";
    break;

  default:

    break;
  }

  hepmc_jet_src src{name, embeddingID, algorithmName, algorithm, parameter, tagPID, tagStatus};

  m_jetSrc.push_back(src);
}

int JetHepMCLoader::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    cout << PHWHERE << "DST Node missing, doing nothing." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  for (const hepmc_jet_src &src : m_jetSrc)
  {
    // Create the AntiKt node if required
    PHCompositeNode *AlgoNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", src.m_algorithmName.c_str()));
    if (!AlgoNode)
    {
      AlgoNode = new PHCompositeNode(src.m_algorithmName.c_str());
      dstNode->addNode(AlgoNode);
    }

    // Create the Input node if required
    PHCompositeNode *InputNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", m_jetInputCategory.c_str()));
    if (!InputNode)
    {
      InputNode = new PHCompositeNode(m_jetInputCategory.c_str());
      AlgoNode->addNode(InputNode);
    }

    JetMap *jets = findNode::getClass<JetMap>(topNode, src.m_name);
    if (!jets)
    {
      jets = new JetMapV1();
      PHIODataNode<PHObject> *JetMapNode = new PHIODataNode<PHObject>(jets, src.m_name.c_str(), "PHObject");
      InputNode->addNode(JetMapNode);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

Fun4AllHistoManager *
JetHepMCLoader::getHistoManager()
{
  string histname(Name() + "_HISTOS");

  Fun4AllServer *se = Fun4AllServer::instance();
  Fun4AllHistoManager *hm = se->getHistoManager(histname);

  if (not hm)
  {
    cout
        << "TPCDataStreamEmulator::get_HistoManager - Making Fun4AllHistoManager " << histname
        << endl;
    hm = new Fun4AllHistoManager(histname);
    se->registerHistoManager(hm);
  }

  assert(hm);

  return hm;
}
