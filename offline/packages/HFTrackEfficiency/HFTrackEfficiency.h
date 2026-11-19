// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef HFTRACKEFFICIENCY_H
#define HFTRACKEFFICIENCY_H

#include <fun4all/SubsysReco.h>
#include <g4eval/SvtxEvalStack.h>
#include <g4eval/SvtxTrackEval.h>

#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace CLHEP
{
  class Hep3Vector;
}

class DecayFinderContainerBase;
class PHCompositeNode;
class PHG4Particle;
class PHG4TruthInfoContainer;
class PHG4ParticleSvtxMap;
class PHHepMCGenEvent;
class PHHepMCGenEventMap;
class SvtxTrack;
class SvtxTrackMap;
class TFile;
class TrackSeedContainer;
class TTree;

class HFTrackEfficiency : public SubsysReco
{
 public:
  explicit HFTrackEfficiency(const std::string &name = "HFTrackEfficiency");

  ~HFTrackEfficiency() override = default;

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
  using Decay = std::vector<std::pair<std::pair<int, int>, int>>;

  bool m_triggerOnDecay{false};

  PHG4TruthInfoContainer *m_truthInfo{nullptr};
  PHHepMCGenEventMap *m_geneventmap{nullptr};
  PHHepMCGenEvent *m_genevt{nullptr};

  // Only needed to find the clusters of a track when the track stores seed IDs instead of seed pointers (SvtxTrack_v5)
  TrackSeedContainer *m_siliconSeeds{nullptr};
  TrackSeedContainer *m_tpcSeeds{nullptr};

  SvtxEvalStack *m_svtx_evalstack{nullptr};
  SvtxTrackEval *trackeval{nullptr};
  //PHG4ParticleSvtxMap *m_dst_truth_reco_map{nullptr};

  DecayFinderContainerBase *m_decayMap{nullptr};
  std::string m_df_module_name;

  SvtxTrackMap *m_input_trackMap{nullptr};
  SvtxTrackMap *m_output_trackMap{nullptr};
  SvtxTrack *m_dst_track{nullptr};
  std::string m_input_track_map_node_name;
  std::string m_output_track_map_node_name;
  std::string outputNodeName;
  bool m_write_track_map{false};

  std::string m_outfile_name;
  TFile *m_outfile{nullptr};
  TTree *m_tree{nullptr};
  bool m_write_nTuple{false};

  unsigned int m_counter_allDecays{0};
  unsigned int m_counter_acceptedDecays{0};
  unsigned int m_counter_unmatchedTruth{0};  // Daughters with no Geant4 particle found
  float m_truthRecoMatchPercent;

  unsigned int m_nDaughters;
  std::string m_decay_descriptor;

  bool findTracks(PHCompositeNode *topNode, const Decay &decay);
  PHG4Particle *findGeant4Primary(int embeddingID, int barcode, int pdgID, const CLHEP::Hep3Vector &labMomentum) const;
  bool getGeant4Vertices(const PHG4Particle *motherG4, const PHG4Particle *daughterG4, CLHEP::Hep3Vector &primary, CLHEP::Hep3Vector &production, CLHEP::Hep3Vector &decay) const;
  void countClusters(SvtxTrack *track, PHG4Particle *truthParticle, unsigned int index);
  void initializeBranches();
  void resetBranches();
  void getDecayDescriptor();
  void getNDaughters();
  std::string getParticleName(const int PDGID);
  float getParticleMass(const int PDGID);

  static const int m_maxTracks{5};
  static const int m_nDetectors{4};  // MVTX, INTT, TPC, TPOT
  bool m_all_tracks_reconstructed{false};
  bool m_is_primary{false};
  float m_true_mother_mass{std::numeric_limits<float>::quiet_NaN()};
  float m_reco_mother_mass{std::numeric_limits<float>::quiet_NaN()};
  float m_true_mother_pT{std::numeric_limits<float>::quiet_NaN()};
  float m_reco_mother_pT{std::numeric_limits<float>::quiet_NaN()};
  float m_true_mother_p{std::numeric_limits<float>::quiet_NaN()};
  float m_true_mother_eta{std::numeric_limits<float>::quiet_NaN()};
  float m_true_mother_rapidity{std::numeric_limits<float>::quiet_NaN()};
  float m_true_mother_phi{std::numeric_limits<float>::quiet_NaN()};
  float m_min_true_track_pT{std::numeric_limits<float>::max()};
  float m_min_reco_track_pT{std::numeric_limits<float>::max()};
  float m_max_true_track_pT{std::numeric_limits<float>::lowest()};
  float m_max_reco_track_pT{std::numeric_limits<float>::lowest()};
  bool m_reco_track_exists[m_maxTracks]{false};
  bool m_used_truth_reco_map[m_maxTracks]{false};
  float m_true_track_pT[m_maxTracks]{std::numeric_limits<float>::quiet_NaN()};
  float m_reco_track_pT[m_maxTracks]{std::numeric_limits<float>::quiet_NaN()};
  float m_true_track_eta[m_maxTracks]{std::numeric_limits<float>::quiet_NaN()};
  float m_reco_track_eta[m_maxTracks]{std::numeric_limits<float>::quiet_NaN()};
  float m_true_track_rapidity[m_maxTracks]{std::numeric_limits<float>::quiet_NaN()};
  float m_true_track_phi[m_maxTracks]{std::numeric_limits<float>::quiet_NaN()};
  float m_reco_track_phi[m_maxTracks]{std::numeric_limits<float>::quiet_NaN()};
  float m_true_track_PID[m_maxTracks]{std::numeric_limits<float>::quiet_NaN()};
  float m_reco_track_chi2nDoF[m_maxTracks]{std::numeric_limits<float>::quiet_NaN()};
  // Clusters on the matched reco track (N) and those of them that come from the truth daughter (T), by detector
  // Wrong clusters = N - T, contribution fraction = T / N
  int m_reco_track_nclusters[m_maxTracks][m_nDetectors]{};
  int m_reco_track_ntrue_clusters[m_maxTracks][m_nDetectors]{};
  // True vertices [cm] in the Geant4 lab frame (the frame the tracks are reconstructed in)
  float m_primary_vtx_x{std::numeric_limits<float>::quiet_NaN()};  // Start of the HepMC/Geant4 ancestry of the mother
  float m_primary_vtx_y{std::numeric_limits<float>::quiet_NaN()};
  float m_primary_vtx_z{std::numeric_limits<float>::quiet_NaN()};
  float m_production_vtx_x{std::numeric_limits<float>::quiet_NaN()};  // Where the mother was produced
  float m_production_vtx_y{std::numeric_limits<float>::quiet_NaN()};
  float m_production_vtx_z{std::numeric_limits<float>::quiet_NaN()};
  float m_decay_vtx_x{std::numeric_limits<float>::quiet_NaN()};
  float m_decay_vtx_y{std::numeric_limits<float>::quiet_NaN()};
  float m_decay_vtx_z{std::numeric_limits<float>::quiet_NaN()};
  float m_true_decay_length{std::numeric_limits<float>::quiet_NaN()};     // production -> decay vertex of the mother [cm]
  float m_true_decay_time{std::numeric_limits<float>::quiet_NaN()};       // mother proper time over the same flight [ps]
};

#endif  // HFTRACKEFFICIENCY_H
