#include "QAKFParticle.h"

#include <qautils/QAHistManagerDef.h>  // for getHistoManager

#include <kfparticle_sphenix/KFParticle_Container.h>
#include <kfparticle_sphenix/KFParticle_Tools.h>

//#include <g4eval/SvtxClusterEval.h>
//#include <g4eval/SvtxEvalStack.h>  // for SvtxEvalStack

#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>

#include <trackbase_historic/SvtxTrack.h>  // for SvtxTrack
#include <trackbase_historic/SvtxTrackMap.h>

#include <trackbase/TrkrDefs.h>  // for cluskey

#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/getClass.h>

#include <KFParticle.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#pragma GCC diagnostic pop

#include <HepMC/GenParticle.h>
#include <HepMC/IteratorRange.h>

#include <TH1.h>
#include <TH2.h>
#include <TString.h>

#include <CLHEP/Vector/LorentzVector.h>
#include <CLHEP/Vector/ThreeVector.h>

#include <cassert>
#include <cstdlib>  // for abs, NULL
#include <iostream>
#include <map>
#include <string>
#include <utility>  // for pair
#include <vector>

KFParticle_Tools kfpTools;

QAKFParticle::QAKFParticle(const std::string &name, const std::string &mother_name, double min_m, double max_m)
  : SubsysReco(name)
{
  m_min_mass = min_m;
  m_max_mass = max_m;
  m_mother_id = kfpTools.getParticleID(mother_name);
  // m_mother_name = mother_name;
}

int QAKFParticle::InitRun(PHCompositeNode *topNode)
{
  /*
if (!m_svtxEvalStack)
{
  m_svtxEvalStack.reset(new SvtxEvalStack(topNode));
  m_svtxEvalStack->set_strict(false);
  m_svtxEvalStack->set_verbosity(Verbosity());
}
*/
  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
  if (!m_trackMap)
  {
    std::cout << __PRETTY_FUNCTION__ << " Fatal Error : missing " << m_trackMapName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int QAKFParticle::Init(PHCompositeNode * /*topNode*/)
{
  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  TH1 *h(nullptr);

  TH2 *h2(nullptr);

  h = new TH1F(TString(get_histo_prefix()) + "InvMass_KFP",  //
               ";mass [GeV/c^{2}];Entries", 100, m_min_mass, m_max_mass);
  hm->registerHisto(h);

  float eta_min = 0.;
  float eta_max = 1.1;
  float phi_min = -3.14;
  float phi_max = 3.14;
  float pt_min = 0;
  float pt_max = 10;

  h2 = new TH2F(TString(get_histo_prefix()) + "InvMass_KFP_Eta",  //
                ";mass [GeV/c^{2}]; Eta", 100, m_min_mass, m_max_mass, 100, eta_min, eta_max);
  hm->registerHisto(h2);

  h2 = new TH2F(TString(get_histo_prefix()) + "InvMass_KFP_Phi",  //
                ";mass [GeV/c^{2}]; Phi", 100, m_min_mass, m_max_mass, 100, phi_min, phi_max);
  hm->registerHisto(h2);

  h2 = new TH2F(TString(get_histo_prefix()) + "InvMass_KFP_pT",  //
                ";pT [GeV]; mass [GeV/c^{2}] ", 100, pt_min, pt_max, 100, m_min_mass, m_max_mass);
  hm->registerHisto(h2);

  h_mass_KFP = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "InvMass_KFP"));
  assert(h_mass_KFP);

  h_mass_KFP_eta = dynamic_cast<TH2F *>(hm->getHisto(get_histo_prefix() + "InvMass_KFP_Eta"));
  assert(h_mass_KFP_eta);

  h_mass_KFP_phi = dynamic_cast<TH2F *>(hm->getHisto(get_histo_prefix() + "InvMass_KFP_Phi"));
  assert(h_mass_KFP_phi);

  h_mass_KFP_pt = dynamic_cast<TH2F *>(hm->getHisto(get_histo_prefix() + "InvMass_KFP_pT"));
  assert(h_mass_KFP_pt);

  return Fun4AllReturnCodes::EVENT_OK;
}

int QAKFParticle::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 2)
  {
    std::cout << "QAKFParticle::process_event() entered" << std::endl;
  }

  // load relevant nodes from NodeTree
  load_nodes(topNode);

  // histogram manager

  /*
  if (m_svtxEvalStack)
    m_svtxEvalStack->next_event(topNode);

  std::vector<CLHEP::HepLorentzVector> daughters;
  for (auto &[key, track] : *m_trackMap)
  {
    SvtxTrack *thisTrack = getTrack(key, m_trackMap);
    CLHEP::HepLorentzVector *theVector = makeHepLV(topNode, thisTrack->get_id());
    if (theVector) daughters.push_back(*theVector);
  }

  CLHEP::HepLorentzVector mother;
  if (daughters.size() >= 2)
  {
    for (CLHEP::HepLorentzVector daughter : daughters)
    {
      mother += daughter;
    }
  }

  //h_mass->Fill(mother.m());

  daughters.clear();
  */

  // float m_calculated_mother_decaytime = -1;
  // float m_calculated_mother_decaytime_err = -1;
  // const float speedOfLight = 2.99792458e-1;

  // std::map<unsigned int, KFParticle*> Map = m_kfpContainer->returnParticlesByPDGid(m_mother_id);

  /*
  std::vector<int> d_id;
  std::vector<KFParticle> vertex_vec = kfpTools.makeAllPrimaryVertices(topNode, "SvtxVertexMap");

  for (KFParticle_Container::Iter iter = m_kfpContainer->begin(); iter != m_kfpContainer->end(); ++iter)
  {
    if (iter->second->GetPDG() != m_mother_id)
    {
      if (d_id.size() == 0)
      {
        d_id.push_back(abs(iter->second->GetPDG()));
      }
      else
      {
        for (unsigned int j = 0; j < d_id.size(); ++j)
        {
          if (abs(iter->second->GetPDG()) == d_id[j])
          {
            continue;
          }
          else if (j == d_id.size() - 1)
          {
            d_id.push_back(abs(iter->second->GetPDG()));
          }
        }
      }
    }
  }
  */

  for (auto &iter : *m_kfpContainer)
  {
    if (iter.second->GetPDG() == m_mother_id)
    {
      // filling mother histogram information
      h_mass_KFP->Fill(iter.second->GetMass());

      h_mass_KFP_eta->Fill(iter.second->GetMass(), iter.second->GetEta());

      h_mass_KFP_phi->Fill(iter.second->GetMass(), iter.second->GetPhi());

      h_mass_KFP_pt->Fill(iter.second->GetPt(), iter.second->GetMass());
      // h_DecayTime->Fill(part->GetLifeTime());
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

/*
SvtxTrack *QAKFParticle::getTrack(unsigned int track_id, SvtxTrackMap *trackmap)
{
  SvtxTrack *matched_track = NULL;

  for (SvtxTrackMap::Iter iter = trackmap->begin();
       iter != trackmap->end();
       ++iter)
  {
    if (iter->first == track_id) matched_track = iter->second;
  }

  return matched_track;
}
*/

/*
PHG4Particle *QAKFParticle::getTruthTrack(SvtxTrack *thisTrack)
{
  if (!clustereval)
  {
    clustereval = m_svtxEvalStack->get_cluster_eval();
  }

  TrkrDefs::cluskey clusKey = *thisTrack->begin_cluster_keys();
  PHG4Particle *particle = clustereval->max_truth_particle_by_cluster_energy(clusKey);

  return particle;
}
*/

/*
CLHEP::HepLorentzVector *QAKFParticle::makeHepLV(PHCompositeNode *topNode, int track_number)
{
  SvtxTrack *track = getTrack(track_number, m_trackMap);
  PHG4Particle *g4particle = getTruthTrack(track);
  CLHEP::HepLorentzVector *lvParticle = NULL;

  PHHepMCGenEventMap *m_geneventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
  if (!m_geneventmap)
  {
    std::cout << "Missing node PHHepMCGenEventMap" << std::endl;
    std::cout << "You will have no mother information" << std::endl;
    return NULL;
  }

  PHHepMCGenEvent *m_genevt = m_geneventmap->get(1);
  if (!m_genevt)
  {
    std::cout << "Missing node PHHepMCGenEvent" << std::endl;
    std::cout << "You will have no mother information" << std::endl;
    return nullptr;
  }

  HepMC::GenEvent *theEvent = m_genevt->getEvent();

  bool breakOut = false;
  for (HepMC::GenEvent::particle_const_iterator p = theEvent->particles_begin(); p != theEvent->particles_end(); ++p)
  {
    assert((*p));
    if ((*p)->barcode() == g4particle->get_barcode())
    {
      for (HepMC::GenVertex::particle_iterator mother = (*p)->production_vertex()->particles_begin(HepMC::parents);
           mother != (*p)->production_vertex()->particles_end(HepMC::parents); ++mother)
      {
        if (abs((*mother)->pdg_id()) == m_mother_id)
        {
          lvParticle = new CLHEP::HepLorentzVector();
          lvParticle->setVectM(CLHEP::Hep3Vector(track->get_px(), track->get_py(), track->get_pz()), kfpTools.getParticleMass((*p)->pdg_id()));
        }
        else
          continue;
        break;
      }
      breakOut = true;
    }
    if (breakOut) break;
  }

  return lvParticle;
}
*/

int QAKFParticle::load_nodes(PHCompositeNode *topNode)
{
  m_kfpContainer = findNode::getClass<KFParticle_Container>(topNode, "reconstructedParticles_KFParticle_Container");
  if (!m_kfpContainer)
  {
    std::cout << "reconstructedParticles_KFParticle_Container - Fatal Error - "
              << "unable to find DST node "
              << "KFParticle_QA" << std::endl;
    assert(m_kfpContainer);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

std::string QAKFParticle::get_histo_prefix()
{
  return std::string("h_") + Name() + std::string("_") + m_trackMapName + std::string("_");
}