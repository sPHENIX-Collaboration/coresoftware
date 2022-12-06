#include "QAG4SimulationKFParticle.h"

#include "QAHistManagerDef.h"  // for getHistoManager

#include <kfparticle_sphenix/KFParticle_Container.h>
#include <kfparticle_sphenix/KFParticle_Tools.h>

#include <g4eval/SvtxClusterEval.h>
#include <g4eval/SvtxEvalStack.h>  // for SvtxEvalStack

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

QAG4SimulationKFParticle::QAG4SimulationKFParticle(const std::string &name, const std::string &mother_name, double min_m, double max_m)
  : SubsysReco(name)
{
  m_min_mass = min_m;
  m_max_mass = max_m;
  m_mother_id = kfpTools.getParticleID(mother_name);
  m_mother_name = mother_name;
}

int QAG4SimulationKFParticle::InitRun(PHCompositeNode *topNode)
{
  if (!m_svtxEvalStack)
  {
    m_svtxEvalStack.reset(new SvtxEvalStack(topNode));
    m_svtxEvalStack->set_strict(false);
    m_svtxEvalStack->set_verbosity(Verbosity());
  }
  m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
  if (!m_trackMap)
  {
    std::cout << __PRETTY_FUNCTION__ << " Fatal Error : missing " << m_trackMapName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_truthInfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!m_trackMap)
  {
    std::cout << __PRETTY_FUNCTION__ << " Fatal Error : missing G4TruthInfo" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int QAG4SimulationKFParticle::Init(PHCompositeNode * /*topNode*/)
{
  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  TH1 *h(nullptr);

  h = new TH1F(TString(get_histo_prefix()) + "InvMass",  //
               ";mass [GeV/c^{2}];Entries", 100, m_min_mass, m_max_mass);
  hm->registerHisto(h);

  h = new TH1F(TString(get_histo_prefix()) + "InvMass_KFP",  //
               ";mass [GeV/c^{2}];Entries", 100, m_min_mass, m_max_mass);
  hm->registerHisto(h);

  h = new TH1F(TString(get_histo_prefix()) + "DecayTime",  //
               ";Decay Time [ps];Entries", 100, 0, 1.5);
  hm->registerHisto(h);

  h = new TH1F(TString(get_histo_prefix()) + "pT",  //
               ";pT [GeV/c];Entries", 100, 0, 10);
  hm->registerHisto(h);

  h = new TH1F(TString(get_histo_prefix()) + "Chi2_NDF",  //
               ";#chi^{2}/NDF ;Entries", 100, 0, 5);
  hm->registerHisto(h);

  h = new TH1F(TString(get_histo_prefix()) + "Rapidity",  //
               ";y;Entries", 100, -2, 2);
  hm->registerHisto(h);

  h = new TH1F(TString(get_histo_prefix()) + "Mother_DCA_XY",  //
               ";DCA [cm];Entries", 100, -0.05, 0.05);
  hm->registerHisto(h);

  h = new TH1F(TString(get_histo_prefix()) + "Daughter1_pT",  //
               ";pT [GeV/c];Entries", 100, 0, 10);
  hm->registerHisto(h);
  h = new TH1F(TString(get_histo_prefix()) + "Daughter2_pT",  //
               ";pT [GeV/c];Entries", 100, 0, 10);
  hm->registerHisto(h);
  h = new TH1F(TString(get_histo_prefix()) + "Daughter3_pT",  //
               ";pT [GeV/c];Entries", 100, 0, 10);
  hm->registerHisto(h);
  h = new TH1F(TString(get_histo_prefix()) + "Daughter4_pT",  //
               ";pT [GeV/c];Entries", 100, 0, 10);
  hm->registerHisto(h);

  h = new TH1F(TString(get_histo_prefix()) + "Daughter1_DCA_XY_Mother",  //
               ";DCA [cm];Entries", 100, -0.05, 0.05);
  hm->registerHisto(h);
  h = new TH1F(TString(get_histo_prefix()) + "Daughter2_DCA_XY_Mother",  //
               ";DCA [cm];Entries", 100, -0.05, 0.05);
  hm->registerHisto(h);
  h = new TH1F(TString(get_histo_prefix()) + "Daughter3_DCA_XY_Mother",  //
               ";DCA [cm];Entries", 100, -0.05, 0.05);
  hm->registerHisto(h);
  h = new TH1F(TString(get_histo_prefix()) + "Daughter4_DCA_XY_Mother",  //
               ";DCA [cm];Entries", 100, -0.05, 0.05);
  hm->registerHisto(h);

  return Fun4AllReturnCodes::EVENT_OK;
}

int QAG4SimulationKFParticle::process_event(PHCompositeNode *topNode)
{
  if (Verbosity() > 2)
    std::cout << "QAG4SimulationKFParticle::process_event() entered" << std::endl;

  // load relevant nodes from NodeTree
  load_nodes(topNode);

  // histogram manager
  Fun4AllHistoManager *hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  TH1F *h_mass = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "InvMass"));
  assert(h_mass);
  TH1F *h_mass_KFP = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "InvMass_KFP"));
  assert(h_mass_KFP);
  TH1F *h_DecayTime = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "DecayTime"));
  assert(h_DecayTime);
  TH1F *h_pT = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "pT"));
  assert(h_pT);
  TH1F *h_Chi2_NDF = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "Chi2_NDF"));
  assert(h_Chi2_NDF);
  TH1F *h_Rapidity = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "Rapidity"));
  assert(h_Rapidity);
  TH1F *h_Mother_DCA_XY = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "Mother_DCA_XY"));
  assert(h_Mother_DCA_XY);
  TH1F *h_Daughter1_pT = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "Daughter1_pT"));
  assert(h_Daughter1_pT);
  TH1F *h_Daughter1_DCA_XY_Mother = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "Daughter1_DCA_XY_Mother"));
  assert(h_Daughter1_DCA_XY_Mother);
  TH1F *h_Daughter2_pT = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "Daughter2_pT"));
  assert(h_Daughter2_pT);
  TH1F *h_Daughter2_DCA_XY_Mother = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "Daughter2_DCA_XY_Mother"));
  assert(h_Daughter2_DCA_XY_Mother);
  TH1F *h_Daughter3_pT = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "Daughter3_pT"));
  assert(h_Daughter3_pT);
  TH1F *h_Daughter3_DCA_XY_Mother = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "Daughter3_DCA_XY_Mother"));
  assert(h_Daughter3_DCA_XY_Mother);
  TH1F *h_Daughter4_pT = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "Daughter4_pT"));
  assert(h_Daughter4_pT);
  TH1F *h_Daughter4_DCA_XY_Mother = dynamic_cast<TH1F *>(hm->getHisto(get_histo_prefix() + "Daughter4_DCA_XY_Mother"));
  assert(h_Daughter4_DCA_XY_Mother);

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

  h_mass->Fill(mother.m());

  daughters.clear();

  float m_calculated_mother_decaytime = -1;
  float m_calculated_mother_decaytime_err = -1;
  const float speedOfLight = 2.99792458e-1;

  // std::map<unsigned int, KFParticle*> Map = m_kfpContainer->returnParticlesByPDGid(m_mother_id);
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

  for (KFParticle_Container::Iter iter = m_kfpContainer->begin(); iter != m_kfpContainer->end(); ++iter)
  {
    if (iter->second->GetPDG() == m_mother_id)
    {
      // filling mother histogram information
      h_mass_KFP->Fill(iter->second->GetMass());
      // h_DecayTime->Fill(part->GetLifeTime());
      h_pT->Fill(iter->second->GetPt());
      h_Chi2_NDF->Fill(iter->second->Chi2() / iter->second->NDF());
      h_Rapidity->Fill(iter->second->GetRapidity());
      // best PV fit for mother
      int bestCombinationIndex = 0;
      if (vertex_vec.size() > 0)
      {
        for (unsigned int i = 1; i < vertex_vec.size(); ++i)
        {
          if (iter->second->GetDeviationFromVertex(vertex_vec[i]) <
              iter->second->GetDeviationFromVertex(vertex_vec[bestCombinationIndex]))
          {
            bestCombinationIndex = i;
          }
        }
        h_Mother_DCA_XY->Fill(iter->second->GetDistanceFromVertexXY(vertex_vec[bestCombinationIndex]));

        iter->second->SetProductionVertex(vertex_vec[bestCombinationIndex]);
        iter->second->GetLifeTime(m_calculated_mother_decaytime, m_calculated_mother_decaytime_err);
        // part->GetDecayLength(m_calculated_mother_decaylength, m_calculated_mother_decaylength_err);
        m_calculated_mother_decaytime /= speedOfLight;
        // m_calculated_mother_decaytime_err /= speedOfLight;

        h_DecayTime->Fill(m_calculated_mother_decaytime);

        for (unsigned int i = 0; i < d_id.size(); ++i)
        {
          std::map<unsigned int, KFParticle *> D_Map = m_kfpContainer->returnParticlesByPDGid(d_id[i]);
          for (auto &[key, part] : D_Map)
          {
            if (i == 0)
            {
              h_Daughter1_pT->Fill(part->GetPt());
              h_Daughter1_DCA_XY_Mother->Fill(part->GetDistanceFromVertexXY(vertex_vec[bestCombinationIndex]));
            }
            if (i == 1)
            {
              h_Daughter2_pT->Fill(part->GetPt());
              h_Daughter2_DCA_XY_Mother->Fill(part->GetDistanceFromVertexXY(vertex_vec[bestCombinationIndex]));
            }
            if (i == 2)
            {
              h_Daughter3_pT->Fill(part->GetPt());
              h_Daughter3_DCA_XY_Mother->Fill(part->GetDistanceFromVertexXY(vertex_vec[bestCombinationIndex]));
            }
            if (i == 3)
            {
              h_Daughter4_pT->Fill(part->GetPt());
              h_Daughter4_DCA_XY_Mother->Fill(part->GetDistanceFromVertexXY(vertex_vec[bestCombinationIndex]));
            }
          }
        }
      }
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

SvtxTrack *QAG4SimulationKFParticle::getTrack(unsigned int track_id, SvtxTrackMap *trackmap)
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

PHG4Particle *QAG4SimulationKFParticle::getTruthTrack(SvtxTrack *thisTrack)
{
  if (!clustereval)
  {
    clustereval = m_svtxEvalStack->get_cluster_eval();
  }

  TrkrDefs::cluskey clusKey = *thisTrack->begin_cluster_keys();
  PHG4Particle *particle = clustereval->max_truth_particle_by_cluster_energy(clusKey);

  return particle;
}

CLHEP::HepLorentzVector *QAG4SimulationKFParticle::makeHepLV(PHCompositeNode *topNode, int track_number)
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

int QAG4SimulationKFParticle::load_nodes(PHCompositeNode *topNode)
{
  m_truthContainer = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!m_truthContainer)
  {
    std::cout << "QAG4SimulationTracking::load_nodes - Fatal Error - "
              << "unable to find DST node "
              << "G4TruthInfo" << std::endl;
    assert(m_truthContainer);
  }
  m_kfpContainer = findNode::getClass<KFParticle_Container>(topNode, m_mother_name + "_KFParticle_Container");
  if (!m_kfpContainer)
  {
    std::cout << m_mother_name.c_str() << "_KFParticle_Container - Fatal Error - "
              << "unable to find DST node "
              << "G4_QA" << std::endl;
    assert(m_kfpContainer);
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

std::string QAG4SimulationKFParticle::get_histo_prefix()
{
  return std::string("h_") + Name() + std::string("_") + m_trackMapName + std::string("_");
}
