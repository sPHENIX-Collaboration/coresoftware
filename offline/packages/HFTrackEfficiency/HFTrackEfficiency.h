// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef HFTRACKEFFICIENCY_H
#define HFTRACKEFFICIENCY_H

#include <decayfinder/DecayFinderContainer_v1.h>  // for DecayFinderContainer_v1
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4VtxPoint.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <trackbase_historic/PHG4ParticleSvtxMap_v1.h>
#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackMap_v2.h>

#include <CLHEP/Vector/LorentzVector.h>
#include <CLHEP/Vector/ThreeVector.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>  // for GenVertex::particle_iterator
#pragma GCC diagnostic pop

#include <HepMC/GenParticle.h>
#include <HepMC/IteratorRange.h>
#include <HepMC/SimpleVector.h>

#include <TBranch.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TTree.h>

#include <string>

class PHCompositeNode;
class PHG4TruthInfoContainer;
class PHG4Particle;
class PHHepMCGenEvent;
class PHHepMCGenEventMap;

namespace CLHEP
{
  class HepLorentzVector;
}

namespace HepMC
{
  class GenParticle;
}

class HFTrackEfficiency : public SubsysReco
{
 public:
  HFTrackEfficiency(const std::string &name = "HFTrackEfficiency");

  ~HFTrackEfficiency() override;

  int Init(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int End(PHCompositeNode *topNode) override;

  void PrintEff();

  void setDFNodeName(const std::string &name) { m_df_module_name = name; }
  void setInputTrackMapName(const std::string &what) { m_input_track_map_node_name = what; }
  void setOutputTrackMapName(const std::string &what) { m_output_track_map_node_name = what; }
  void writeSelectedTrackMap(bool write) { m_write_track_map = write; }
  void writeOutputFile(bool write) { m_write_nTuple = write; }
  void setOutputFileName(const std::string &what) { m_outfile_name = what; }
  void triggerOnDecay(bool trigger) { m_triggerOnDecay = trigger; }
  void setTruthRecoMatchingPercentage(float value) { m_truthRecoMatchPercent = value; };

 private:
  typedef std::vector<std::pair<std::pair<int, int>, int>> Decay;

  bool m_triggerOnDecay;

  PHG4TruthInfoContainer *m_truthInfo = nullptr;
  PHHepMCGenEventMap *m_geneventmap = nullptr;
  PHHepMCGenEvent *m_genevt = nullptr;

  DecayFinderContainer_v1 *m_decayMap = nullptr;
  std::string m_df_module_name;

  SvtxTrackMap *m_input_trackMap = nullptr;
  SvtxTrackMap *m_output_trackMap = nullptr;
  SvtxTrack *m_dst_track = nullptr;
  std::string m_input_track_map_node_name;
  std::string m_output_track_map_node_name;
  std::string outputNodeName;
  bool m_write_track_map;

  std::string m_outfile_name;
  TFile *m_outfile;
  TTree *m_tree;
  bool m_write_nTuple;

  unsigned int m_counter_allDecays = 0;
  unsigned int m_counter_acceptedDecays = 0;
  float m_truthRecoMatchPercent;

  unsigned int m_nDaughters;
  std::string m_decay_descriptor;

  bool findTracks(PHCompositeNode *topNode, Decay decay);
  void initializeBranches();
  void resetBranches();
  void getDecayDescriptor();
  void getNDaughters();
  std::string getParticleName(const int PDGID);
  float getParticleMass(const int PDGID);

  static const int m_maxTracks = 5;
  bool m_all_tracks_reconstructed = false;
  float m_reco_mother_mass = 0.;
  float m_true_mother_pT = 0.;
  float m_true_mother_eta = 0.;
  float m_min_true_track_pT = FLT_MAX;
  float m_min_reco_track_pT = FLT_MAX;
  float m_max_true_track_pT = -1. * FLT_MAX;
  float m_max_reco_track_pT = -1. * FLT_MAX;
  bool m_reco_track_exists[m_maxTracks] = {false};
  float m_true_track_pT[m_maxTracks] = {0.};
  float m_reco_track_pT[m_maxTracks] = {0.};
  float m_true_track_eta[m_maxTracks] = {0.};
  float m_true_track_PID[m_maxTracks] = {0.};
  float m_reco_track_chi2nDoF[m_maxTracks] = {0.};
  int m_reco_track_silicon_seeds[m_maxTracks] = {0};
  int m_reco_track_tpc_seeds[m_maxTracks] = {0};
};

#endif  // HFTRACKEFFICIENCY_H
