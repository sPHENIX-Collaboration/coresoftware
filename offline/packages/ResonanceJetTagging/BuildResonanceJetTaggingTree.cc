#include "BuildResonanceJetTaggingTree.h"

#include "ResonanceJetTagging.h"

#include <phool/phool.h>

/// Jet includes
#include <g4jets/JetMap.h>
#include <g4jets/Jetv1.h>

/// Tracking includes
#include <trackbase_historic/SvtxPHG4ParticleMap_v1.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>

/// HEPMC truth includes
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#pragma GCC diagnostic pop

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

/// Fun4All includes
#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <g4main/PHG4Particle.h>            // for PHG4Particle
#include <g4main/PHG4TruthInfoContainer.h>  // for PHG4TruthInfoContainer
#include <g4main/PHG4VtxPoint.h>            // for PHG4VtxPoint

#include <kfparticle_sphenix/KFParticle_truthAndDetTools.h>

#include <KFParticle.h>
#include <kfparticle_sphenix/KFParticle_Container.h>

/// ROOT includes
#include <TFile.h>
#include <TTree.h>

/// C++ includes
#include <cassert>
#include <sstream>
#include <string>

/**
 * BuildResonanceJetTaggingTree is a class developed to reconstruct jets containing a D-meson
 * The class can be adapted to tag jets using any kind of particle
 * Author: Antonio Silva (antonio.sphenix@gmail.com)
 */

/**
 * Constructor of module
 */
BuildResonanceJetTaggingTree::BuildResonanceJetTaggingTree(const std::string &name, const std::string &filename, const ResonanceJetTagging::TAG tag)
  : SubsysReco(name)
  , m_outfilename(filename)
  , m_tagcontainer_name("")
  , m_jetcontainer_name("")
  , m_truth_jetcontainer_name("")
  , m_dorec(true)
  , m_dotruth(false)
  , m_tag_particle(tag)
  , m_tag_pdg(0)
  , m_outfile(nullptr)
  , m_taggedjettree(nullptr)
{
  /// Initialize variables and trees so we don't accidentally access
  /// memory that was never allocated
  initializeVariables();
  initializeTrees();

  switch (m_tag_particle) {
    case ResonanceJetTagging::TAG::D0:
      m_tag_pdg = 421;
      break;
    case ResonanceJetTagging::TAG::DPLUS:
      m_tag_pdg = 411;
      break;
    case ResonanceJetTagging::TAG::DSTAR:
      m_tag_pdg = 413;
      break;
    case ResonanceJetTagging::TAG::JPSY:
      m_tag_pdg = 433;
      break;
    case ResonanceJetTagging::TAG::K0:
      m_tag_pdg = 311;
      break;
    case ResonanceJetTagging::TAG::GAMMA:
      m_tag_pdg = 22;
      break;
    case ResonanceJetTagging::TAG::ELECTRON:
      m_tag_pdg = 11;
      break;
  }
}

/**
 * Destructor of module
 */
BuildResonanceJetTaggingTree::~BuildResonanceJetTaggingTree()
{
  delete m_taggedjettree;
}

/**
 * Initialize the module and prepare looping over events
 */
int BuildResonanceJetTaggingTree::Init(PHCompositeNode */*topNode*/)
{
  if (Verbosity() > 5)
  {
    std::cout << "Beginning Init in BuildResonanceJetTaggingTree" << std::endl;
  }

  return 0;
}

/**
 * Main workhorse function where each event is looped over and
 * data from each event is collected from the node tree for analysis
 */
int BuildResonanceJetTaggingTree::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 5)
  {
    std::cout << "Beginning process_event in BuildResonanceJetTaggingTree" << std::endl;
  }

  switch (m_tag_particle) {
    case ResonanceJetTagging::TAG::D0:
      return loopD0(topNode);
      break;
    case ResonanceJetTagging::TAG::DPLUS:
      break;
    case ResonanceJetTagging::TAG::DSTAR:
      break;
    case ResonanceJetTagging::TAG::JPSY:
      break;
    case ResonanceJetTagging::TAG::K0:
      break;
    case ResonanceJetTagging::TAG::GAMMA:
      break;
    case ResonanceJetTagging::TAG::ELECTRON:
      break;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

/**
 * End the module and finish any data collection. Clean up any remaining
 * loose ends
 */
int BuildResonanceJetTaggingTree::End(PHCompositeNode */*topNode*/)
{
  if (Verbosity() > 1)
  {
    std::cout << "Ending BuildResonanceJetTaggingTree analysis package" << std::endl;
  }

  /// Change to the outfile
  m_outfile->cd();

  m_taggedjettree->Write();

  /// Write and close the outfile
  m_outfile->Write();
  m_outfile->Close();

  delete m_outfile;

  if (Verbosity() > 1)
  {
    std::cout << "Finished BuildResonanceJetTaggingTree analysis package" << std::endl;
  }

  return 0;
}

int BuildResonanceJetTaggingTree::loopD0(PHCompositeNode *topNode)
{
  KFParticle_Container* kfContainer = nullptr;

  if(m_dorec)
  {
    m_taggedJetMap = getJetMapFromNode(topNode, m_jetcontainer_name);
    if(!m_taggedJetMap) return Fun4AllReturnCodes::ABORTEVENT;

    kfContainer = getKFParticleContainerFromNode(topNode, m_tagcontainer_name);
    if(!kfContainer) return Fun4AllReturnCodes::ABORTEVENT;
  }

  HepMC::GenEvent *hepMCGenEvent = nullptr;

  if(m_dotruth)
  {
    m_truth_taggedJetMap = getJetMapFromNode(topNode, m_truth_jetcontainer_name);
    if(!m_truth_taggedJetMap) return Fun4AllReturnCodes::ABORTEVENT;

    hepMCGenEvent = getGenEventFromNode(topNode, "PHHepMCGenEventMap");
    if(!hepMCGenEvent) return Fun4AllReturnCodes::ABORTEVENT;
  }

  //*********

  Jet *recTagJet = nullptr;
  Jet *genTagJet = nullptr;

  KFParticle *recTag = nullptr;
  HepMC::GenParticle *genTag = nullptr;

  int recDaughtersID[2];

  std::vector<int> recJetIndex;

  if(m_dorec)
  {
    for (JetMapv1::Iter iter = m_taggedJetMap->begin(); iter != m_taggedJetMap->end(); ++iter)
    {
      recTagJet = iter->second;

      if(!recTagJet) continue;

      Jetv1::Iter recTagIter = recTagJet->find(Jet::SRC::VOID);

      if(recTagIter == recTagJet->end_comp()) continue;

      recTag = kfContainer->get(recTagIter->second);

      if(!recTag) continue;
      recDaughtersID[0] = (kfContainer->get(recTagIter->second + 1))->Id();
      recDaughtersID[1] = (kfContainer->get(recTagIter->second + 2))->Id();

      resetTreeVariables();

      m_tagpartpx = recTag->Px();
      m_tagpartpy = recTag->Py();
      m_tagpartpz = recTag->Pz();
      m_tagpartpt = recTag->GetPt();
      m_tagparteta = recTag->GetEta();
      m_tagpartphi = recTag->GetPhi();
      m_tagpartm = recTag->GetMass();

      m_tagjetpx = recTagJet->get_px();
      m_tagjetpy = recTagJet->get_py();
      m_tagjetpz = recTagJet->get_pz();
      m_tagjetpt = recTagJet->get_pt();
      m_tagjeteta = recTagJet->get_eta();
      m_tagjetphi = recTagJet->get_phi();

      genTagJet = nullptr;
      genTag = nullptr;

      if(m_dotruth)
      {
        findMatchedTruthD0(topNode, genTagJet, genTag, recDaughtersID);

        if((!genTagJet) || (!genTag))
        {
          m_taggedjettree->Fill();
        }
        else
        {
          recJetIndex.push_back(genTagJet->get_id());

          m_truth_tagpartpx = genTag->momentum().px();
          m_truth_tagpartpy = genTag->momentum().py();
          m_truth_tagpartpz = genTag->momentum().pz();
          m_truth_tagpartpt = std::sqrt(m_truth_tagpartpx * m_truth_tagpartpx + m_truth_tagpartpy * m_truth_tagpartpy);
          m_truth_tagparteta = atanh(m_truth_tagpartpz / genTag->momentum().e());
          m_truth_tagpartphi = atan(m_truth_tagpartpy / m_truth_tagpartpx);

          m_truth_tagjetpx = genTagJet->get_px();
          m_truth_tagjetpy = genTagJet->get_py();
          m_truth_tagjetpz = genTagJet->get_pz();
          m_truth_tagjetpt = genTagJet->get_pt();
          m_truth_tagjeteta = genTagJet->get_eta();
          m_truth_tagjetphi = genTagJet->get_phi();

          m_taggedjettree->Fill();
        }
      }
      else
      {
        m_taggedjettree->Fill();
      }
    }
  }

  if(m_dotruth)
  {
    resetTreeVariables();

    for (JetMapv1::Iter iter = m_truth_taggedJetMap->begin(); iter != m_truth_taggedJetMap->end(); ++iter)
    {
      genTagJet = iter->second;

      if(!genTagJet) continue;

      //Check if truth was matched to reconstructed
      if(isReconstructed(genTagJet->get_id(), recJetIndex)) continue;

      Jetv1::Iter genTagIter = genTagJet->find(Jet::SRC::VOID);

      genTag = hepMCGenEvent->barcode_to_particle(genTagIter->second);

      m_truth_tagpartpx = genTag->momentum().px();
      m_truth_tagpartpy = genTag->momentum().py();
      m_truth_tagpartpz = genTag->momentum().pz();
      m_truth_tagpartpt = std::sqrt(m_truth_tagpartpx * m_truth_tagpartpx + m_truth_tagpartpy * m_truth_tagpartpy);
      m_truth_tagparteta = atanh(m_truth_tagpartpz / genTag->momentum().e());
      m_truth_tagpartphi = atan(m_truth_tagpartpy / m_truth_tagpartpx);

      m_truth_tagjetpx = genTagJet->get_px();
      m_truth_tagjetpy = genTagJet->get_py();
      m_truth_tagjetpz = genTagJet->get_pz();
      m_truth_tagjetpt = genTagJet->get_pt();
      m_truth_tagjeteta = genTagJet->get_eta();
      m_truth_tagjetphi = genTagJet->get_phi();

      m_taggedjettree->Fill();
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void BuildResonanceJetTaggingTree::findMatchedTruthD0(PHCompositeNode *topNode, Jet *&mcTagJet, HepMC::GenParticle *&mcTag, int decays[])
{
  m_truth_taggedJetMap = getJetMapFromNode(topNode, "D0Jets_Truth_Jet_Container");
  if(!m_truth_taggedJetMap) return;

  HepMC::GenEvent *hepMCGenEvent = getGenEventFromNode(topNode, "PHHepMCGenEventMap");
  if(!hepMCGenEvent) return;

  PHG4Particle *g4particle = nullptr;
  const int nDecays = 2;
  HepMC::GenParticle *mcTags[nDecays];

  PHNodeIterator nodeIter(topNode);
  PHNode *findNode = dynamic_cast<PHNode *>(nodeIter.findFirst("SvtxPHG4ParticleMap"));
  PHG4TruthInfoContainer *truthinfo = nullptr;
  if (findNode)
  {
    findNode = dynamic_cast<PHNode *>(nodeIter.findFirst("G4TruthInfo"));
    if (findNode)
    {
      truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
    }
    else
    {
      std::cout << "KFParticle truth matching: G4TruthInfo does not exist" << std::endl;
      return;
    }
  }

  // Truth map
  SvtxPHG4ParticleMap_v1 *dst_reco_truth_map = findNode::getClass<SvtxPHG4ParticleMap_v1>(topNode, "SvtxPHG4ParticleMap");

  if (!dst_reco_truth_map)
  {
    return;
  }

  for (int idecay = 0; idecay < nDecays; idecay++)
  {
    std::map<float, std::set<int>> truth_set = dst_reco_truth_map->get(decays[idecay]);
    const auto &best_weight = truth_set.rbegin();
    int best_truth_id = *best_weight->second.rbegin();
    g4particle = truthinfo->GetParticle(best_truth_id);

    mcTags[idecay] = getMother(topNode, g4particle);

    if (mcTags[idecay] == nullptr)
    {
      return;
    }
  }
  // check is decays are from the same mother, otherwise it is background
  for (int idecay = 1; idecay < nDecays; idecay++)
  {
    if (mcTags[idecay]->barcode() != mcTags[idecay - 1]->barcode())
    {
      return;
    }
  }

  mcTag = mcTags[0];

  for (JetMapv1::Iter iter = m_truth_taggedJetMap->begin(); iter != m_truth_taggedJetMap->end(); ++iter)
  {
    mcTagJet = iter->second;

    Jetv1::Iter mcTagIter = mcTagJet->find(Jet::SRC::VOID);

    if(int(mcTagIter->second) != mcTag->barcode())
    {
      mcTagJet = nullptr;
      continue;
    }
    else
    {
      break;
    }

  }

  return;
}

HepMC::GenParticle *BuildResonanceJetTaggingTree::getMother(PHCompositeNode *topNode, PHG4Particle *g4daughter)
{
  PHHepMCGenEventMap *hepmceventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");

  /// If the node was not properly put on the tree, return
  if (!hepmceventmap)
  {
    std::cout << PHWHERE
              << "HEPMC event map node is missing, can't collected HEPMC truth particles"
              << std::endl;
    return nullptr;
  }

  PHHepMCGenEvent *hepmcevent = hepmceventmap->get(1);
  if (!hepmcevent)
  {
    return nullptr;
  }

  HepMC::GenEvent *hepMCGenEvent = hepmcevent->getEvent();
  if (!hepMCGenEvent)
  {
    return nullptr;
  }

  HepMC::GenParticle *mcDaughter = nullptr;

  if (g4daughter->get_barcode() > 0)
  {
    mcDaughter = hepMCGenEvent->barcode_to_particle(g4daughter->get_barcode());
  }
  if (!mcDaughter)
  {
    return nullptr;
  }

  HepMC::GenVertex *TagVertex = mcDaughter->production_vertex();
  for (HepMC::GenVertex::particle_iterator it = TagVertex->particles_begin(HepMC::ancestors); it != TagVertex->particles_end(HepMC::ancestors); ++it)
  {
    if (std::abs((*it)->pdg_id()) == m_tag_pdg)
    {
      return (*it);
    }
  }

  return nullptr;
}

bool BuildResonanceJetTaggingTree::isReconstructed(int index, std::vector<int> indexRecVector)
{
  for(auto indexRec : indexRecVector)
  {
    if(index == indexRec) return true;
  }
  return false;
}

void BuildResonanceJetTaggingTree::initializeTrees()
{
  delete m_taggedjettree;
  m_taggedjettree = new TTree("m_taggedjettree", "A tree with Tagged-Jet info");
  m_taggedjettree->Branch("m_tagpartpx", &m_tagpartpx, "m_tagpartpx/D");
  m_taggedjettree->Branch("m_tagpartpy", &m_tagpartpy, "m_tagpartpy/D");
  m_taggedjettree->Branch("m_tagpartpz", &m_tagpartpz, "m_tagpartpz/D");
  m_taggedjettree->Branch("m_tagpartpt", &m_tagpartpt, "m_tagpartpt/D");
  m_taggedjettree->Branch("m_tagparteta", &m_tagparteta, "m_tagparteta/D");
  m_taggedjettree->Branch("m_tagpartphi", &m_tagpartphi, "m_tagpartphi/D");
  m_taggedjettree->Branch("m_tagpartm", &m_tagpartm, "m_tagpartm/D");
  m_taggedjettree->Branch("m_tagjetpx", &m_tagjetpx, "m_tagjetpx/D");
  m_taggedjettree->Branch("m_tagjetpy", &m_tagjetpy, "m_tagjetpy/D");
  m_taggedjettree->Branch("m_tagjetpz", &m_tagjetpz, "m_tagjetpz/D");
  m_taggedjettree->Branch("m_tagjetpt", &m_tagjetpt, "m_tagjetpt/D");
  m_taggedjettree->Branch("m_tagjeteta", &m_tagjeteta, "m_tagjeteta/D");
  m_taggedjettree->Branch("m_tagjetphi", &m_tagjetphi, "m_tagjetphi/D");

  m_taggedjettree->Branch("m_truth_tagpartpx", &m_truth_tagpartpx, "m_truth_tagpartpx/D");
  m_taggedjettree->Branch("m_truth_tagpartpy", &m_truth_tagpartpy, "m_truth_tagpartpy/D");
  m_taggedjettree->Branch("m_truth_tagpartpz", &m_truth_tagpartpz, "m_truth_tagpartpz/D");
  m_taggedjettree->Branch("m_truth_tagpartpt", &m_truth_tagpartpt, "m_truth_tagpartpt/D");
  m_taggedjettree->Branch("m_truth_tagparteta", &m_truth_tagparteta, "m_truth_tagparteta/D");
  m_taggedjettree->Branch("m_truth_tagpartphi", &m_truth_tagpartphi, "m_truth_tagpartphi/D");
  m_taggedjettree->Branch("m_truth_tagjetpx", &m_truth_tagjetpx, "m_truth_tagjetpx/D");
  m_taggedjettree->Branch("m_truth_tagjetpy", &m_truth_tagjetpy, "m_truth_tagjetpy/D");
  m_taggedjettree->Branch("m_truth_tagjetpz", &m_truth_tagjetpz, "m_truth_tagjetpz/D");
  m_taggedjettree->Branch("m_truth_tagjetpt", &m_truth_tagjetpt, "m_truth_tagjetpt/D");
  m_taggedjettree->Branch("m_truth_tagjeteta", &m_truth_tagjeteta, "m_truth_tagjeteta/D");
  m_taggedjettree->Branch("m_truth_tagjetphi", &m_truth_tagjetphi, "m_truth_tagjetphi/D");

}
void BuildResonanceJetTaggingTree::initializeVariables()
{
  delete m_outfile;
  m_outfile = new TFile(m_outfilename.c_str(), "RECREATE");
}

void BuildResonanceJetTaggingTree::resetTreeVariables()
{
  m_tagpartpx = NAN;
  m_tagpartpy = NAN;
  m_tagpartpz = NAN;
  m_tagpartpt = NAN;
  m_tagparteta = NAN;
  m_tagpartphi = NAN;
  m_tagpartm = NAN;
  m_tagjetpx = NAN;
  m_tagjetpy = NAN;
  m_tagjetpz = NAN;
  m_tagjetpt = NAN;
  m_tagjeteta = NAN;
  m_tagjetphi = NAN;
  //Truth info
  m_truth_tagpartpx = NAN;
  m_truth_tagpartpy = NAN;
  m_truth_tagpartpz = NAN;
  m_truth_tagpartpt = NAN;
  m_truth_tagparteta = NAN;
  m_truth_tagpartphi = NAN;
  m_truth_tagjetpx = NAN;
  m_truth_tagjetpy = NAN;
  m_truth_tagjetpz = NAN;
  m_truth_tagjetpt = NAN;
  m_truth_tagjeteta = NAN;
  m_truth_tagjetphi = NAN;
}

JetMapv1* BuildResonanceJetTaggingTree::getJetMapFromNode(PHCompositeNode *topNode, const std::string &name)
{
  JetMapv1 *jetmap = findNode::getClass<JetMapv1>(topNode, name);

  if (!jetmap)
  {
    std::cout << PHWHERE
         << "JetMap node is missing, can't collect jets"
         << std::endl;
    return 0;
  }

  return jetmap;
}
KFParticle_Container* BuildResonanceJetTaggingTree::getKFParticleContainerFromNode(PHCompositeNode *topNode, const std::string &name)
{
  KFParticle_Container *cont = findNode::getClass<KFParticle_Container>(topNode, name);

  if (!cont)
  {
    std::cout << PHWHERE
         << "KFParticle_Container node is missing, can't collect HF particles"
         << std::endl;
    return 0;
  }

  return cont;
}
HepMC::GenEvent* BuildResonanceJetTaggingTree::getGenEventFromNode(PHCompositeNode *topNode, const std::string &name)
{
  PHHepMCGenEventMap *hepmceventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, name);

  /// If the node was not properly put on the tree, return
  if (!hepmceventmap)
  {
    std::cout << PHWHERE
         << "HEPMC event map node is missing, can't collected HEPMC truth particles"
         << std::endl;
    return 0;
  }

  PHHepMCGenEvent *hepmcevent = hepmceventmap->get(1);

  if (!hepmcevent)
  {
    std::cout << PHWHERE
         << "PHHepMCGenEvent node is missing, can't collected HEPMC truth particles"
         << std::endl;
    return 0;
  }

  HepMC::GenEvent *hepMCGenEvent = hepmcevent->getEvent();

  if (!hepmceventmap)
  {
    std::cout << PHWHERE
         << "HepMC::GenEvent node is missing, can't collected HEPMC truth particles"
         << std::endl;
    return 0;
  }

  return hepMCGenEvent;

}
