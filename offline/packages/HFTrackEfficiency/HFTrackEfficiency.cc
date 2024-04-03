#include "HFTrackEfficiency.h"

//____________________________________________________________________________..
HFTrackEfficiency::HFTrackEfficiency(const std::string &name)
  : SubsysReco(name)
  , m_triggerOnDecay(false)
  , m_input_track_map_node_name("SvtxTrackMap")
  , m_output_track_map_node_name("HFSelected")
  , m_write_track_map(false)
  , m_outfile_name("outputHFTrackEff.root")
  , m_outfile(nullptr)
  , m_tree(nullptr)
  , m_write_nTuple(false)
  , m_truthRecoMatchPercent(5.)
  , m_nDaughters(2)
{
}

//____________________________________________________________________________..
HFTrackEfficiency::~HFTrackEfficiency() = default;

//____________________________________________________________________________..
int HFTrackEfficiency::Init(PHCompositeNode *topNode)
{
  m_truthRecoMatchPercent /= 100.;

  if (m_write_track_map)
  {
    PHNodeIterator nodeIter(topNode);

    PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(nodeIter.findFirst("PHCompositeNode", "DST"));
    if (!dstNode)
    {
      dstNode = new PHCompositeNode("DST");
      topNode->addNode(dstNode);
      std::cout << "DST node added" << std::endl;
    }

    outputNodeName = m_output_track_map_node_name + "_SvtxTrackMap";
    m_output_trackMap = new SvtxTrackMap_v2();
    PHIODataNode<PHObject> *outputTrackNode = new PHIODataNode<PHObject>(m_output_trackMap, outputNodeName.c_str(), "PHObject");
    dstNode->addNode(outputTrackNode);
    std::cout << outputNodeName << " node added" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int HFTrackEfficiency::process_event(PHCompositeNode *topNode)
{
  PHNodeIterator nodeIter(topNode);
  PHCompositeNode *dstNode = static_cast<PHCompositeNode *>(nodeIter.findFirst("PHCompositeNode", "DST"));
  assert(dstNode);

  std::string df_node_name = m_df_module_name + "_DecayMap";
  m_decayMap = findNode::getClass<DecayFinderContainer_v1>(topNode, df_node_name.c_str());
  if (!m_decayMap)
  {
    std::cout << __FILE__ << ": Missing node " << df_node_name << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_input_trackMap = findNode::getClass<SvtxTrackMap>(topNode, m_input_track_map_node_name.c_str());
  if (!m_input_trackMap)
  {
    std::cout << __FILE__ << ": Missing node " << m_input_track_map_node_name << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_geneventmap = findNode::getClass<PHHepMCGenEventMap>(topNode, "PHHepMCGenEventMap");
  if (!m_geneventmap)
  {
    std::cout << __FILE__ << ": Missing node PHHepMCGenEventMap" << std::endl;
  }

  m_truthInfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!m_truthInfo)
  {
    if (Verbosity() >= VERBOSITY_MORE)
    {
      std::cout << __FILE__ << ": Missing node G4TruthInfo" << std::endl;
    }
  }

  m_dst_truth_reco_map = findNode::getClass<PHG4ParticleSvtxMap_v1>(topNode, "PHG4ParticleSvtxMap");
  if (m_dst_truth_reco_map)
  {
    if (Verbosity() >= VERBOSITY_MORE)
    {
      std::cout << __FILE__ << ": PHG4ParticleSvtxMap found, truth matching will be more accurate" << std::endl;
    }
  }
  else
  {
    if (Verbosity() >= VERBOSITY_MORE)
    {
      std::cout << __FILE__ << ": PHG4ParticleSvtxMap not found, reverting to true matching by momentum relations. Truth matching will be less accurate" << std::endl;
    }
  }

  if (m_decay_descriptor.empty() && !m_decayMap->empty())
  {
    getDecayDescriptor();
    getNDaughters();
  }

  bool reconstructableDecay = false;
  for (auto &iter : *m_decayMap)
  {
    Decay decay = iter.second;
    ++m_counter_allDecays;

    m_all_tracks_reconstructed = findTracks(topNode, decay);
    if (m_all_tracks_reconstructed)
    {
      ++m_counter_acceptedDecays;
      reconstructableDecay = true;
    }

    if (m_write_nTuple)
    {
      if (m_counter_allDecays == 1)
      {
        initializeBranches();
      }
      m_tree->Fill();
      resetBranches();
    }
  }

  if (m_triggerOnDecay && !reconstructableDecay)
  {
    if (Verbosity() >= VERBOSITY_MORE)
    {
      std::cout << "No decays were reconstructable in this event, skipping" << std::endl;
    }
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  else
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }
}

int HFTrackEfficiency::End(PHCompositeNode * /*topNode*/)
{
  PrintEff();

  if (m_write_nTuple && m_counter_allDecays != 0)
  {
    m_outfile->Write();
    m_outfile->Close();
    delete m_outfile;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

bool HFTrackEfficiency::findTracks(PHCompositeNode *topNode, Decay decay)
{
  int trackableParticles[] = {11, 13, 211, 321, 2212};
  bool recoTrackFound = false;
  std::set<SvtxTrack *> selectedTracks;

  CLHEP::HepLorentzVector motherRecoLV;
  CLHEP::HepLorentzVector daughterSumTrueLV;
  CLHEP::HepLorentzVector *motherTrueLV = new CLHEP::HepLorentzVector();
  CLHEP::HepLorentzVector *daughterTrueLV = new CLHEP::HepLorentzVector();

  HepMC::GenEvent *theEvent = nullptr;
  if (m_geneventmap)
  {
    m_genevt = m_geneventmap->get(decay[0].first.first);
    assert(m_genevt);

    theEvent = m_genevt->getEvent();
    HepMC::GenParticle *mother = theEvent->barcode_to_particle(decay[0].first.second);
    assert(mother);
    if (Verbosity() >= VERBOSITY_MORE)
    {
      mother->print();
    }

    m_true_mother_pT = mother->momentum().perp();
    m_true_mother_p = std::sqrt(std::pow(mother->momentum().px(), 2) + std::pow(mother->momentum().py(), 2) + std::pow(mother->momentum().pz(), 2));  // Must have an old HepMC build, no mag function
    m_true_mother_eta = mother->momentum().eta();

    HepMC::GenVertex *thisVtx = mother->production_vertex();
    m_primary_vtx_x = thisVtx->point3d().x();
    m_primary_vtx_y = thisVtx->point3d().y();
    m_primary_vtx_z = thisVtx->point3d().z();
  }

  for (unsigned int i = 1; i < decay.size(); ++i)
  {
    m_dst_track = nullptr;
    int truth_ID = -1;
    if (std::find(std::begin(trackableParticles), std::end(trackableParticles),
                  std::abs(decay[i].second)) != std::end(trackableParticles))
    {
      if (theEvent && decay[i].first.second > -1)
      {
        HepMC::GenParticle *daughterHepMC = theEvent->barcode_to_particle(decay[i].first.second);
        if (Verbosity() >= VERBOSITY_MORE)
        {
          daughterHepMC->print();
        }

        daughterTrueLV->setVectM(CLHEP::Hep3Vector(daughterHepMC->momentum().px(), daughterHepMC->momentum().py(), daughterHepMC->momentum().pz()), getParticleMass(decay[i].second));
        daughterSumTrueLV += *daughterTrueLV;

        m_true_track_PID[i - 1] = daughterHepMC->pdg_id();

        // Now get the decay vertex position
        HepMC::GenVertex *thisVtx = daughterHepMC->production_vertex();
        m_secondary_vtx_x = thisVtx->point3d().x();
        m_secondary_vtx_y = thisVtx->point3d().y();
        m_secondary_vtx_z = thisVtx->point3d().z();

        // We need the G4 ID, not the HepMC ID to use the truth/reco map
        if (m_dst_truth_reco_map)
        {
          PHG4TruthInfoContainer::ConstRange range = m_truthInfo->GetParticleRange();

          for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
          {
            PHG4Particle *daughterG4 = iter->second;

            if (abs(daughterG4->get_px() - daughterTrueLV->x()) <= 5e-3 &&
                abs(daughterG4->get_py() - daughterTrueLV->y()) <= 5e-3 &&
                abs(daughterG4->get_pz() - daughterTrueLV->z()) <= 5e-3 && daughterG4->get_pid() == decay[i].second)
            {
              truth_ID = daughterG4->get_track_id();
              break;
            }
          }
        }
      }
      else
      {
        PHG4TruthInfoContainer::ConstRange range = m_truthInfo->GetParticleRange();

        for (PHG4TruthInfoContainer::ConstIterator iter = range.first; iter != range.second; ++iter)
        {
          PHG4Particle *daughterG4 = iter->second;

          PHG4Particle *motherG4 = nullptr;
          if (daughterG4->get_parent_id() != 0)
          {
            motherG4 = m_truthInfo->GetParticle(daughterG4->get_parent_id());
          }
          else
          {
            continue;
          }

          if (motherG4->get_pid() == decay[0].second && motherG4->get_barcode() == decay[0].first.second && daughterG4->get_pid() == decay[i].second && daughterG4->get_barcode() == decay[i].first.second)
          {
            if (Verbosity() >= VERBOSITY_MORE)
            {
              daughterG4->identify();
            }

            CLHEP::Hep3Vector *mother3Vector = new CLHEP::Hep3Vector(motherG4->get_px(), motherG4->get_py(), motherG4->get_pz());
            motherTrueLV->setVectM((*mother3Vector), getParticleMass(decay[0].second));
            m_true_mother_pT = motherTrueLV->perp();
            m_true_mother_p = mother3Vector->mag();
            m_true_mother_eta = motherTrueLV->pseudoRapidity();

            PHG4VtxPoint *thisVtx = m_truthInfo->GetVtx(motherG4->get_vtx_id());
            m_primary_vtx_x = thisVtx->get_x();
            m_primary_vtx_y = thisVtx->get_y();
            m_primary_vtx_z = thisVtx->get_z();

            daughterTrueLV->setVectM(CLHEP::Hep3Vector(daughterG4->get_px(), daughterG4->get_py(), daughterG4->get_pz()), getParticleMass(decay[i].second));
            daughterSumTrueLV += *daughterTrueLV;

            // Now get the decay vertex position
            thisVtx = m_truthInfo->GetVtx(daughterG4->get_vtx_id());
            m_secondary_vtx_x = thisVtx->get_x();
            m_secondary_vtx_y = thisVtx->get_y();
            m_secondary_vtx_z = thisVtx->get_z();

            m_true_track_PID[i - 1] = daughterG4->get_pid();
            truth_ID = daughterG4->get_track_id();

            delete mother3Vector;
          }
        }
      }

      m_true_track_pT[i - 1] = (float) daughterTrueLV->perp();
      m_true_track_eta[i - 1] = (float) daughterTrueLV->pseudoRapidity();
      m_min_true_track_pT = std::min(m_true_track_pT[i - 1], m_min_true_track_pT);
      m_max_true_track_pT = std::max(m_true_track_pT[i - 1], m_max_true_track_pT);

      if (m_dst_truth_reco_map && truth_ID >= 0)
      {
        std::map<float, std::set<unsigned int>> reco_set = m_dst_truth_reco_map->get(truth_ID);
        if (reco_set.size() == 0)
        {
          continue;
        }
        const auto &best_weight = reco_set.rbegin();
        if (best_weight->second.size() == 0)
        {
          continue;
        }
        unsigned int best_reco_id = *best_weight->second.rbegin();
        m_dst_track = m_input_trackMap->get(best_reco_id);
        if (m_dst_track)
        {
          m_used_truth_reco_map[i - 1] = true;
          recoTrackFound = true;
        }
      }
      else
      {
        for (auto &iter : *m_input_trackMap)
        {
          m_dst_track = iter.second;
          float delta_px = (m_dst_track->get_px() - daughterTrueLV->px()) / daughterTrueLV->px();
          float delta_py = (m_dst_track->get_py() - daughterTrueLV->py()) / daughterTrueLV->py();
          float delta_pz = (m_dst_track->get_pz() - daughterTrueLV->pz()) / daughterTrueLV->pz();

          if (std::abs(delta_px) <= m_truthRecoMatchPercent && std::abs(delta_py) <= m_truthRecoMatchPercent && std::abs(delta_pz) <= m_truthRecoMatchPercent)
          {
            recoTrackFound = true;
            break;
          }
        }
      }

      if (recoTrackFound)
      {
        selectedTracks.insert(m_dst_track);
        if (Verbosity() >= VERBOSITY_MORE)
        {
          m_dst_track->identify();
        }
        m_reco_track_exists[i - 1] = true;
        m_reco_track_pT[i - 1] = m_dst_track->get_pt();
        m_reco_track_eta[i - 1] = m_dst_track->get_eta();
        m_reco_track_chi2nDoF[i - 1] = m_dst_track->get_chisq() / m_dst_track->get_ndf();
        if (m_dst_track->get_silicon_seed())
        {
          m_reco_track_silicon_seeds[i - 1] = static_cast<int>(m_dst_track->get_silicon_seed()->size_cluster_keys());
        }
        else
        {
          m_reco_track_silicon_seeds[i - 1] = 0;
        }
        m_reco_track_tpc_seeds[i - 1] = static_cast<int>(m_dst_track->get_tpc_seed()->size_cluster_keys());
        m_min_reco_track_pT = std::min(m_reco_track_pT[i - 1], m_min_reco_track_pT);
        m_max_reco_track_pT = std::max(m_reco_track_pT[i - 1], m_max_reco_track_pT);

        CLHEP::HepLorentzVector *daughterRecoLV = new CLHEP::HepLorentzVector();
        daughterRecoLV->setVectM(CLHEP::Hep3Vector(m_dst_track->get_px(), m_dst_track->get_py(), m_dst_track->get_pz()), getParticleMass(m_true_track_PID[i - 1]));

        motherRecoLV += *daughterRecoLV;
        delete daughterRecoLV;
      }

      recoTrackFound = false;
    }
  }

  m_true_mother_mass = daughterSumTrueLV.m();
  bool foundDecay = true;
  if (selectedTracks.size() == m_nDaughters)
  {
    m_reco_mother_mass = motherRecoLV.m();
    if (m_write_track_map)
    {
      m_output_trackMap = findNode::getClass<SvtxTrackMap>(topNode, outputNodeName.c_str());
      for (auto &track : selectedTracks)
      {
        m_output_trackMap->insertWithKey(track, track->get_id());
      }
    }
  }
  else
  {
    foundDecay = false;
  }

  selectedTracks.clear();

  delete motherTrueLV;
  delete daughterTrueLV;

  return foundDecay;
}

void HFTrackEfficiency::initializeBranches()
{
  m_outfile = new TFile(m_outfile_name.c_str(), "RECREATE");
  m_tree = new TTree("HFTrackEfficiency", "HFTrackEfficiency");
  m_tree->OptimizeBaskets();
  m_tree->SetAutoSave(-5e6);  // Save the output file every 5MB

  m_tree->Branch("all_tracks_reconstructed", &m_all_tracks_reconstructed, "all_tracks_reconstructed/O");
  m_tree->Branch("true_mother_mass", &m_true_mother_mass, "true_mother_mass/F");
  m_tree->Branch("reco_mother_mass", &m_reco_mother_mass, "reco_mother_mass/F");
  m_tree->Branch("true_mother_pT", &m_true_mother_pT, "true_mother_pT/F");
  m_tree->Branch("true_mother_p", &m_true_mother_p, "true_mother_p/F");
  m_tree->Branch("true_mother_eta", &m_true_mother_eta, "true_mother_eta/F");
  m_tree->Branch("min_true_track_pT", &m_min_true_track_pT, "min_true_track_pT/F");
  m_tree->Branch("min_reco_track_pT", &m_min_reco_track_pT, "min_reco_track_pT/F");
  m_tree->Branch("max_true_track_pT", &m_max_true_track_pT, "max_true_track_pT/F");
  m_tree->Branch("max_reco_track_pT", &m_max_reco_track_pT, "max_reco_track_pT/F");

  for (unsigned int iTrack = 0; iTrack < m_nDaughters; ++iTrack)
  {
    std::string daughter_number = "track_" + std::to_string(iTrack + 1);
    m_tree->Branch("reco_" + TString(daughter_number) + "_exists", &m_reco_track_exists[iTrack], "reco_" + TString(daughter_number) + "_exists/O");
    m_tree->Branch("reco_" + TString(daughter_number) + "_used_truth_reco_map", &m_used_truth_reco_map[iTrack], "reco_" + TString(daughter_number) + "_used_truth_reco_map/O");
    m_tree->Branch("true_" + TString(daughter_number) + "_pT", &m_true_track_pT[iTrack], "true_" + TString(daughter_number) + "_pT/F");
    m_tree->Branch("reco_" + TString(daughter_number) + "_pT", &m_reco_track_pT[iTrack], "reco_" + TString(daughter_number) + "_pT/F");
    m_tree->Branch("true_" + TString(daughter_number) + "_eta", &m_true_track_eta[iTrack], "true_" + TString(daughter_number) + "_eta/F");
    m_tree->Branch("reco_" + TString(daughter_number) + "_eta", &m_reco_track_eta[iTrack], "reco_" + TString(daughter_number) + "_eta/F");
    m_tree->Branch("true_" + TString(daughter_number) + "_PID", &m_true_track_PID[iTrack], "true_" + TString(daughter_number) + "_PID/F");
    m_tree->Branch("reco_" + TString(daughter_number) + "_chi2nDoF", &m_reco_track_chi2nDoF[iTrack], "reco_" + TString(daughter_number) + "_chi2nDoF/F");
    m_tree->Branch("reco_" + TString(daughter_number) + "_silicon_seeds", &m_reco_track_silicon_seeds[iTrack], "reco_" + TString(daughter_number) + "_silicon_seeds/I");
    m_tree->Branch("reco_" + TString(daughter_number) + "_tpc_seeds", &m_reco_track_tpc_seeds[iTrack], "reco_" + TString(daughter_number) + "_tpc_seeds/I");
  }

  m_tree->Branch("true_primary_vertex_x", &m_primary_vtx_x, "true_primary_vertex_x/F");
  m_tree->Branch("true_primary_vertex_y", &m_primary_vtx_y, "true_primary_vertex_y/F");
  m_tree->Branch("true_primary_vertex_z", &m_primary_vtx_z, "true_primary_vertex_z/F");
  m_tree->Branch("true_secondary_vertex_x", &m_secondary_vtx_x, "true_secondary_vertex_x/F");
  m_tree->Branch("true_secondary_vertex_y", &m_secondary_vtx_y, "true_secondary_vertex_y/F");
  m_tree->Branch("true_secondary_vertex_z", &m_secondary_vtx_z, "true_secondary_vertex_z/F");
}

void HFTrackEfficiency::resetBranches()
{
  m_all_tracks_reconstructed = false;
  m_true_mother_mass = std::numeric_limits<float>::quiet_NaN();
  m_reco_mother_mass = std::numeric_limits<float>::quiet_NaN();
  m_true_mother_pT = std::numeric_limits<float>::quiet_NaN();
  m_true_mother_p = std::numeric_limits<float>::quiet_NaN();
  m_true_mother_eta = std::numeric_limits<float>::quiet_NaN();
  m_min_true_track_pT = std::numeric_limits<float>::max();
  m_min_reco_track_pT = std::numeric_limits<float>::max();
  m_max_true_track_pT = -1 * std::numeric_limits<float>::max();
  m_max_reco_track_pT = -1 * std::numeric_limits<float>::max();
  for (unsigned int iTrack = 0; iTrack < m_nDaughters; ++iTrack)
  {
    m_reco_track_exists[iTrack] = false;
    m_used_truth_reco_map[iTrack] = false;
    m_true_track_pT[iTrack] = std::numeric_limits<float>::quiet_NaN();
    m_reco_track_pT[iTrack] = std::numeric_limits<float>::quiet_NaN();
    m_true_track_eta[iTrack] = std::numeric_limits<float>::quiet_NaN();
    m_reco_track_eta[iTrack] = std::numeric_limits<float>::quiet_NaN();
    m_true_track_PID[iTrack] = std::numeric_limits<float>::quiet_NaN();
    m_reco_track_chi2nDoF[iTrack] = std::numeric_limits<float>::quiet_NaN();
    m_reco_track_silicon_seeds[iTrack] = 0;
    m_reco_track_tpc_seeds[iTrack] = 0;
  }

  m_primary_vtx_x = std::numeric_limits<float>::quiet_NaN();
  m_primary_vtx_y = std::numeric_limits<float>::quiet_NaN();
  m_primary_vtx_z = std::numeric_limits<float>::quiet_NaN();
  m_secondary_vtx_x = std::numeric_limits<float>::quiet_NaN();
  m_secondary_vtx_y = std::numeric_limits<float>::quiet_NaN();
  m_secondary_vtx_z = std::numeric_limits<float>::quiet_NaN();
}

void HFTrackEfficiency::getDecayDescriptor()
{
  Decay decay = m_decayMap->begin()->second;
  m_decay_descriptor = getParticleName(decay[0].second);
  m_decay_descriptor += " ->";
  for (unsigned int i = 1; i < decay.size(); ++i)
  {
    m_decay_descriptor += " ";
    m_decay_descriptor += getParticleName(decay[i].second);
  }
}

void HFTrackEfficiency::getNDaughters()
{
  int trackableParticles[] = {11, 13, 211, 321, 2212};
  m_nDaughters = 0;
  Decay decay = m_decayMap->begin()->second;
  for (unsigned int i = 1; i < decay.size(); ++i)
  {
    if (std::find(std::begin(trackableParticles), std::end(trackableParticles),
                  std::abs(decay[i].second)) != std::end(trackableParticles))
    {
      ++m_nDaughters;
    }
  }
}

std::string HFTrackEfficiency::getParticleName(const int PDGID)
{
  return TDatabasePDG::Instance()->GetParticle(PDGID)->GetName();
}

float HFTrackEfficiency::getParticleMass(const int PDGID)
{
  return TDatabasePDG::Instance()->GetParticle(PDGID)->Mass();
}

//____________________________________________________________________________..
void HFTrackEfficiency::PrintEff()
{
  std::cout << "\n--------------- Heavy Flavor Tracking Efficiency ---------------" << std::endl;
  std::cout << "Tracking Efficiency for " << m_decay_descriptor << std::endl;
  std::cout << "Number of decays fully in acceptance: " << m_counter_allDecays << std::endl;
  std::cout << "Number of decays with all tracks reconstructed: " << m_counter_acceptedDecays << std::endl;
  std::cout << "Efficiency: " << std::setw(4) << (float) m_counter_acceptedDecays * 100 / float(m_counter_allDecays) << "%" << std::endl;
  std::cout << "----------------------------------------------------------------\n"
            << std::endl;
}
