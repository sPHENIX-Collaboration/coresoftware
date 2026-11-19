// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef HFTRACKEFFICIENCY_HFTRACKEFFICIENCY_H
#define HFTRACKEFFICIENCY_HFTRACKEFFICIENCY_H

#include <fun4all/SubsysReco.h>

#include <trackbase/TrkrDefs.h>

#include <array>
#include <limits>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace CLHEP
{
  class Hep3Vector;
}

namespace HepMC
{
  class GenParticle;
}

class DecayFinderContainerBase;
class PHCompositeNode;
class PHG4Particle;
class PHG4TruthInfoContainer;
class PHHepMCGenEvent;
class PHHepMCGenEventMap;
class SvtxEvalStack;
class SvtxTrack;
class SvtxTrackEval;
class SvtxTrackMap;
class TFile;
class TrackSeedContainer;
class TTree;

class HFTrackEfficiency : public SubsysReco
{
 public:
  explicit HFTrackEfficiency(const std::string &name = "HFTrackEfficiency");

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
  void setTruthRecoMatchingPercentage(float value) { m_truthRecoMatchPercent = value; }

 private:
  // ((embedding ID, barcode), PDG ID) of the mother followed by its decay products
  using Decay = std::vector<std::pair<std::pair<int, int>, int>>;

  // What is learned about one decay while it is processed, defined in the source file
  struct DecayTruth;
  // What is learned about one daughter of the decay, defined in the source file
  struct DaughterTruth;

  static constexpr unsigned int s_maxTracks{5};
  static constexpr unsigned int s_nDetectors{4};  // MVTX, INTT, TPC, TPOT

  /** @brief Add the output track map to the DST node, only needed if the selected tracks are written. */
  void createOutputTrackMap(PHCompositeNode *topNode);

  /**
   * @brief Find the nodes the event needs.
   * @return Fun4All return code, ABORTRUN if a required node is missing.
   */
  int getNodes(PHCompositeNode *topNode);

  /** @brief Make the evaluation stack on the first event and point it at the new event. */
  void setupEvaluation(PHCompositeNode *topNode);

  /**
   * @brief Work out the decay descriptor and the number of trackable daughters from the first decay found.
   * @return Fun4All return code, ABORTRUN if the decay has more daughters than supported.
   */
  int setupDecay();

  /**
   * @brief Find the tracks of every decay in the event and fill the branches.
   * @return True if at least one decay had all its tracks reconstructed.
   */
  bool processDecays(PHCompositeNode *topNode);

  /**
   * @brief Find the truth and reconstructed information for one decay and fill the branches.
   * @param topNode The top node.
   * @param decay The decay to look at.
   * @return True if every trackable daughter was matched to a different reco track.
   */
  bool findTracks(PHCompositeNode *topNode, const Decay &decay);

  /**
   * @brief Find the truth and reco track of one trackable daughter and fill its branches.
   * @param decay The decay the daughter is part of.
   * @param decayIndex Position of the daughter in the decay.
   * @param trackIndex Position of the daughter amongst the trackable daughters.
   * @param truth Information on the decay, updated with this daughter.
   */
  void processDaughter(const Decay &decay, unsigned int decayIndex, unsigned int trackIndex, DecayTruth &truth);

  /** @brief Fill the mother kinematics and the true vertices from the HepMC record. */
  void fillHepMCMother(const Decay &decay, DecayTruth &truth);

  /** @brief Fill the mother kinematics from HepMC, these are in the generator frame. */
  void fillHepMCMotherKinematics(const HepMC::GenParticle &mother);

  /** @brief Fill the true vertices from HepMC, in the lab frame. */
  void fillHepMCVertices(const HepMC::GenParticle &mother, DecayTruth &truth);

  /** @brief Find the truth information of a daughter in the HepMC record or, if it is not there, in Geant4. */
  void findDaughterTruth(const Decay &decay, unsigned int decayIndex, unsigned int trackIndex, DecayTruth &truth,
                         DaughterTruth &daughter);

  /** @brief Take the truth of a daughter from the HepMC record and find the matching Geant4 particle. */
  void fillHepMCDaughter(const Decay &decay, unsigned int decayIndex, unsigned int trackIndex,
                         const HepMC::GenParticle &daughterHepMC, DecayTruth &truth, DaughterTruth &daughter);

  /**
   * @brief Find how far the Geant4 vertices are shifted from the HepMC ones in the lab frame.
   *
   * HepMCNodeReader can shift and smear the vertex, which cannot be seen from here. The shift is recovered by
   * comparing the production vertex of a HepMC daughter to that of the Geant4 particle made from it.
   */
  void updateVertexOffset(const HepMC::GenParticle &daughterHepMC, const PHG4Particle &daughterG4, DecayTruth &truth);

  /** @brief Find a daughter, and its mother, in the Geant4 truth record. Used when the decay was found there. */
  void findGeant4Daughter(const Decay &decay, unsigned int decayIndex, unsigned int trackIndex, DecayTruth &truth,
                          DaughterTruth &daughter);

  /**
   * @brief Check if a Geant4 particle is the daughter of the decay, and find the mother of the decay if it is.
   * @param decay The decay.
   * @param decayIndex Position of the daughter in the decay.
   * @param candidate The Geant4 particle to check.
   * @return The Geant4 mother of the decay, or nullptr if the candidate is not the daughter.
   */
  PHG4Particle *findGeant4Mother(const Decay &decay, unsigned int decayIndex, const PHG4Particle *candidate) const;

  /** @brief Fill the mother and daughter truth from a Geant4 mother/daughter pair. */
  void fillGeant4Truth(const Decay &decay, PHG4Particle *motherG4, PHG4Particle *daughterG4, unsigned int trackIndex,
                       DecayTruth &truth, DaughterTruth &daughter);

  /**
   * @brief Find the Geant4 primary that HepMCNodeReader made from a HepMC particle.
   * @param embeddingID Embedding ID of the HepMC event the particle is in.
   * @param barcode HepMC barcode of the particle.
   * @param pdgID PDG ID of the particle.
   * @param labMomentum Momentum of the particle in the lab frame [GeV], used if the barcode match fails.
   * @return The Geant4 particle, or nullptr if there is none.
   */
  PHG4Particle *findGeant4Primary(int embeddingID, int barcode, int pdgID, const CLHEP::Hep3Vector &labMomentum) const;

  /**
   * @brief Get the true vertices of a decay from the Geant4 truth record.
   * @param motherG4 The mother, its production vertex is the production vertex of the decay.
   * @param daughterG4 A daughter of the mother, its production vertex is the decay vertex.
   * @param[out] primary Where the parent chain of the mother started [cm].
   * @param[out] production Where the mother was produced [cm].
   * @param[out] decay Where the mother decayed [cm].
   * @return False if any of the vertices is missing.
   */
  bool getGeant4Vertices(const PHG4Particle *motherG4, const PHG4Particle *daughterG4, CLHEP::Hep3Vector &primary,
                         CLHEP::Hep3Vector &production, CLHEP::Hep3Vector &decay) const;

  /** @brief Fill the true kinematics branches of a daughter. */
  void fillTrueDaughterBranches(unsigned int trackIndex, const DaughterTruth &daughter);

  /**
   * @brief Find the reco track of a daughter, from the truth/reco map if there is a Geant4 particle and by
   * momentum otherwise.
   * @return The track, or nullptr if there is none.
   */
  SvtxTrack *findRecoTrack(unsigned int trackIndex, const DaughterTruth &daughter);

  /** @brief Fill the reco branches of a daughter. */
  void fillRecoTrackBranches(SvtxTrack *track, PHG4Particle *truthParticle, unsigned int trackIndex);

  /** @brief Get the keys of the clusters on the seeds of a reco track, also if it only stores the seed IDs. */
  std::set<TrkrDefs::cluskey> collectClusterKeys(SvtxTrack *track) const;

  /**
   * @brief Fill the branches with the number of clusters on a reco track, and how many of those come from a truth
   * particle.
   * @param track The reco track.
   * @param truthParticle The truth daughter the track was matched to, nullptr if there is none.
   * @param trackIndex Index of the daughter.
   */
  void countClusters(SvtxTrack *track, PHG4Particle *truthParticle, unsigned int trackIndex);

  /** @brief Fill the vertex branches, and the decay length and time from the production and decay vertices. */
  void fillVertexBranches(const Decay &decay, DecayTruth &truth);

  /**
   * @brief Fill the reco mother branches and write the selected tracks if every daughter was matched.
   * @return True if every trackable daughter was matched to a different reco track.
   */
  bool storeSelectedTracks(PHCompositeNode *topNode, const DecayTruth &truth);

  void initializeBranches();
  void addMotherBranches();
  void addDaughterBranches(unsigned int trackIndex);
  void addVertexBranches();
  void resetBranches();
  void getDecayDescriptor();
  void getNDaughters();

  bool m_triggerOnDecay{false};

  PHG4TruthInfoContainer *m_truthInfo{nullptr};
  PHHepMCGenEventMap *m_geneventmap{nullptr};
  PHHepMCGenEvent *m_genevt{nullptr};

  // Only needed to find the clusters of a track when the track stores seed IDs instead of seed pointers (SvtxTrack_v5)
  TrackSeedContainer *m_siliconSeeds{nullptr};
  TrackSeedContainer *m_tpcSeeds{nullptr};

  std::unique_ptr<SvtxEvalStack> m_svtx_evalstack;
  SvtxTrackEval *m_trackeval{nullptr};  // Owned by the evaluation stack

  DecayFinderContainerBase *m_decayMap{nullptr};
  std::string m_df_module_name;

  SvtxTrackMap *m_input_trackMap{nullptr};
  SvtxTrackMap *m_output_trackMap{nullptr};
  std::string m_input_track_map_node_name{"SvtxTrackMap"};
  std::string m_output_track_map_node_name{"HFSelected"};
  std::string m_outputNodeName;
  bool m_write_track_map{false};

  std::string m_outfile_name{"outputHFTrackEff.root"};
  std::unique_ptr<TFile> m_outfile;
  TTree *m_tree{nullptr};  // Owned by the output file
  bool m_write_nTuple{false};

  unsigned int m_counter_allDecays{0};
  unsigned int m_counter_acceptedDecays{0};
  unsigned int m_counter_unmatchedTruth{0};  // Daughters with no Geant4 particle found
  float m_truthRecoMatchPercent{5.};

  unsigned int m_nDaughters{2};
  std::string m_decay_descriptor;

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

  std::array<bool, s_maxTracks> m_reco_track_exists{};
  std::array<bool, s_maxTracks> m_used_truth_reco_map{};
  std::array<float, s_maxTracks> m_true_track_pT{};
  std::array<float, s_maxTracks> m_reco_track_pT{};
  std::array<float, s_maxTracks> m_true_track_eta{};
  std::array<float, s_maxTracks> m_reco_track_eta{};
  std::array<float, s_maxTracks> m_true_track_rapidity{};
  std::array<float, s_maxTracks> m_true_track_phi{};
  std::array<float, s_maxTracks> m_reco_track_phi{};
  std::array<float, s_maxTracks> m_true_track_PID{};
  std::array<float, s_maxTracks> m_reco_track_chi2nDoF{};
  // Clusters on the matched reco track (N) and those of them that come from the truth daughter (T), by detector
  // Wrong clusters = N - T, contribution fraction = T / N
  std::array<std::array<int, s_nDetectors>, s_maxTracks> m_reco_track_nclusters{};
  std::array<std::array<int, s_nDetectors>, s_maxTracks> m_reco_track_ntrue_clusters{};

  // True vertices [cm] in the Geant4 lab frame (the frame the tracks are reconstructed in)
  float m_primary_vtx_x{std::numeric_limits<float>::quiet_NaN()};  // Start of the HepMC/Geant4 ancestry of the mother
  float m_primary_vtx_y{std::numeric_limits<float>::quiet_NaN()};
  float m_primary_vtx_z{std::numeric_limits<float>::quiet_NaN()};
  float m_production_vtx_x{std::numeric_limits<float>::quiet_NaN()};  // Where the mother was produced
  float m_production_vtx_y{std::numeric_limits<float>::quiet_NaN()};
  float m_production_vtx_z{std::numeric_limits<float>::quiet_NaN()};
  float m_decay_vtx_x{std::numeric_limits<float>::quiet_NaN()};  // Where the mother decayed
  float m_decay_vtx_y{std::numeric_limits<float>::quiet_NaN()};
  float m_decay_vtx_z{std::numeric_limits<float>::quiet_NaN()};
  float m_true_decay_length{std::numeric_limits<float>::quiet_NaN()};  // production -> decay vertex of the mother [cm]
  float m_true_decay_time{std::numeric_limits<float>::quiet_NaN()};    // mother proper time over the same flight [ps]
};

#endif  // HFTRACKEFFICIENCY_HFTRACKEFFICIENCY_H
