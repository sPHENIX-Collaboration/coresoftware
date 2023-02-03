#ifndef RESONANCEJETTAGGING_H__
#define RESONANCEJETTAGGING_H__

#include <g4jets/Jetv1.h>

#include <fun4all/SubsysReco.h>

#include <fastjet/JetDefinition.hh>

#include <map>      // for map
#include <string>   // for string
#include <utility>  // for pair
#include <vector>

/// Class declarations for use in the analysis module
class PHCompositeNode;
class SvtxTrack;
class PHG4Particlev2;
class PHG4Particle;
class ParticleFlowElement;
class JetMapv1;

namespace CLHEP
{
  class Hep3Vector;
}
namespace HepMC
{
  class GenParticle;
}
namespace fastjet
{
  class PseudoJet;
}

/// Definition of this analysis module class
class ResonanceJetTagging : public SubsysReco
{
 public:
  enum ALGO
  {
    ANTIKT = 0,
    KT = 1,
    CAMBRIDGE = 2
  };

  enum RECOMB
  {
    E_SCHEME = 0,
    PT_SCHEME = 1,
    PT2_SCHEME = 2,
    ET_SCHEME = 3,
    ET2_SCHEME = 4
  };

  enum TAG
  {
    D0 = 0,
    D0TOK3PI = 1,
    DPLUS = 2,
    DSTAR = 3,
    JPSY = 4,
    K0 = 5,
    GAMMA = 6,
    ELECTRON = 7,
    LAMBDAC = 8
  };

  /// Constructor
  ResonanceJetTagging(const std::string &name = "ResonanceJetTagging", const TAG tag = TAG::D0, const std::string &KFparticle_Container_name = "");

  // Destructor
  virtual ~ResonanceJetTagging();

  /// SubsysReco initialize processing method
  int Init(PHCompositeNode *);

  /// SubsysReco event processing method
  int process_event(PHCompositeNode *);

  /// SubsysReco end processing method
  int End(PHCompositeNode *);

  // Particle flow
  void setParticleFlowEtaAcc(double etamin, double etamax)
  {
    m_particleflow_mineta = etamin;
    m_particleflow_maxeta = etamax;
  }
  void setParticleFlowMinEta(double etamin) { m_particleflow_mineta = etamin; }
  void setParticleFlowMaxEta(double etamax) { m_particleflow_maxeta = etamax; }
  double getParticleFlowMinEta() { return m_particleflow_mineta; }
  double getParticleFlowMaxEta() { return m_particleflow_maxeta; }

  // tracks
  void setTrackPtAcc(double ptmin, double ptmax)
  {
    m_track_minpt = ptmin;
    m_track_maxpt = ptmax;
  }
  void setTrackMinPt(double ptmin) { m_track_minpt = ptmin; }
  void setTrackMaxPt(double ptmax) { m_track_maxpt = ptmax; }
  double getTrackMinPt() { return m_track_minpt; }
  double getTrackMaxPt() { return m_track_maxpt; }

  void setTrackEtaAcc(double etamin, double etamax)
  {
    m_track_mineta = etamin;
    m_track_maxeta = etamax;
  }
  void setTrackMinEta(double etamin) { m_track_mineta = etamin; }
  void setTrackMaxEta(double etamax) { m_track_maxeta = etamax; }
  double getTrackMinEta() { return m_track_mineta; }
  double getTrackMaxEta() { return m_track_maxeta; }

  // EMCal clusters
  void setEMCalClusterPtAcc(double ptmin, double ptmax)
  {
    m_EMCal_cluster_minpt = ptmin;
    m_EMCal_cluster_maxpt = ptmax;
  }
  void setEMCalClusterMinPt(double ptmin) { m_EMCal_cluster_minpt = ptmin; }
  void setEMCalClusterMaxPt(double ptmax) { m_EMCal_cluster_maxpt = ptmax; }
  double getEMCalClusterMinPt() { return m_EMCal_cluster_minpt; }
  double getEMCalClusterMaxPt() { return m_EMCal_cluster_maxpt; }

  void setEMCalClusterEtaAcc(double etamin, double etamax)
  {
    m_EMCal_cluster_mineta = etamin;
    m_EMCal_cluster_maxeta = etamax;
  }
  void setEMCalClusterMinEta(double etamin) { m_EMCal_cluster_mineta = etamin; }
  void setEMCalClusterMaxEta(double etamax) { m_EMCal_cluster_maxeta = etamax; }
  double getEMCalClusterMinEta() { return m_EMCal_cluster_mineta; }
  double getEMCalClusterMaxEta() { return m_EMCal_cluster_maxeta; }

  // HCal clusters
  void setHCalClusterPtAcc(double ptmin, double ptmax)
  {
    m_HCal_cluster_minpt = ptmin;
    m_HCal_cluster_maxpt = ptmax;
  }
  void setHCalClusterMinPt(double ptmin) { m_HCal_cluster_minpt = ptmin; }
  void setHCalClusterMaxPt(double ptmax) { m_HCal_cluster_maxpt = ptmax; }
  double getHCalClusterMinPt() { return m_HCal_cluster_minpt; }
  double getHCalClusterMaxPt() { return m_HCal_cluster_maxpt; }

  void setHCalClusterEtaAcc(double etamin, double etamax)
  {
    m_HCal_cluster_mineta = etamin;
    m_HCal_cluster_maxeta = etamax;
  }
  void setHCalClusterMinEta(double etamin) { m_HCal_cluster_mineta = etamin; }
  void setHCalClusterMaxEta(double etamax) { m_HCal_cluster_maxeta = etamax; }
  double getHCalClusterMinEta() { return m_HCal_cluster_mineta; }
  double getHCalClusterMaxEta() { return m_HCal_cluster_maxeta; }

  // Set/Get add Tracks and Clusters
  void setAddParticleFlow(bool b) { m_add_particleflow = b; }
  bool getAddParticleFlow() { return m_add_particleflow; }
  void setAddTracks(bool b) { m_add_tracks = b; }
  bool getAddTracks() { return m_add_tracks; }
  void setAddEMCalClusters(bool b) { m_add_EMCal_clusters = b; }
  bool getAddEMCalClusters() { return m_add_EMCal_clusters; }
  void setAddHCalClusters(bool b) { m_add_HCal_clusters = b; }
  bool getAddHCalClusters() { return m_add_HCal_clusters; }

  // Jet settings
  void setR(double r) { m_jetr = r; }
  double getR(double /*r*/) { return m_jetr; }
  void setJetAlgo(ALGO jetalgo)
  {
    switch (jetalgo)
    {
    case ALGO::ANTIKT:
      m_jetalgo = fastjet::antikt_algorithm;
      break;
    case ALGO::KT:
      m_jetalgo = fastjet::kt_algorithm;
      break;
    case ALGO::CAMBRIDGE:
      m_jetalgo = fastjet::cambridge_algorithm;
      break;
    }
  }
  fastjet::JetAlgorithm getJetAlgo() { return m_jetalgo; }
  void setRecombScheme(RECOMB recomb_scheme)
  {
    switch (recomb_scheme)
    {
    case RECOMB::E_SCHEME:
      m_recomb_scheme = fastjet::E_scheme;
      break;
    case RECOMB::PT_SCHEME:
      m_recomb_scheme = fastjet::pt_scheme;
      break;
    case RECOMB::PT2_SCHEME:
      m_recomb_scheme = fastjet::pt2_scheme;
      break;
    case RECOMB::ET_SCHEME:
      m_recomb_scheme = fastjet::Et_scheme;
      break;
    case RECOMB::ET2_SCHEME:
      m_recomb_scheme = fastjet::Et2_scheme;
      break;
    }
  }
  fastjet::RecombinationScheme getRecombScheme() { return m_recomb_scheme; }
  void setJetParameters(double r, ALGO jetalgo, RECOMB recomb_scheme)
  {
    setR(r);
    setJetAlgo(jetalgo);
    setRecombScheme(recomb_scheme);
  }
  void setJetContainerName(const std::string &n) { m_jetcontainer_name = n; }
  std::string getJetContainerName() { return m_jetcontainer_name; }
  void setDoRecunstructed(bool b) { m_dorec = b; }
  bool getDoRecunstructed() { return m_dorec; }
  void setDoTruth(bool b) { m_dotruth = b; }
  bool getDoTruth() { return m_dotruth; }

 private:
  /// String to contain the jetmap name containing the tagged jets
  std::string m_jetcontainer_name;

  /// Particle Flow selection and acceptance
  double m_particleflow_mineta;
  double m_particleflow_maxeta;

  /// Track selection and acceptance
  double m_track_minpt;
  double m_track_maxpt;
  double m_track_mineta;
  double m_track_maxeta;

  /// EMCal Cluster selection and acceptance
  double m_EMCal_cluster_minpt;
  double m_EMCal_cluster_maxpt;
  double m_EMCal_cluster_mineta;
  double m_EMCal_cluster_maxeta;

  /// HCal Cluster selection and acceptance
  double m_HCal_cluster_minpt;
  double m_HCal_cluster_maxpt;
  double m_HCal_cluster_mineta;
  double m_HCal_cluster_maxeta;

  /// Add Tracks and Clusters
  bool m_add_particleflow;
  bool m_add_tracks;
  bool m_add_EMCal_clusters;
  bool m_add_HCal_clusters;

  // jet settings variables
  double m_jetr;
  fastjet::JetAlgorithm m_jetalgo;
  fastjet::RecombinationScheme m_recomb_scheme;

  JetMapv1 *m_taggedJetMap = nullptr;
  JetMapv1 *m_truth_taggedJetMap = nullptr;

  int m_tag_pdg;
  unsigned int m_jet_id = 0;
  unsigned int m_truth_jet_id = 0;
  bool m_dorec;
  bool m_dotruth;
  int m_nDaughters;
  TAG m_tag_particle;
  std::string m_KFparticle_name;

  /// Methods for grabbing the data
  int tagHFHadronic(PHCompositeNode *topNode);
  void findTaggedJets(PHCompositeNode *topNode, PHG4Particlev2 *Tag, const std::vector<PHG4Particlev2*> &TagDecays);
  void addParticleFlow(PHCompositeNode *topNode, std::vector<fastjet::PseudoJet> &particles, const std::vector<PHG4Particlev2*> &TagDecays, std::map<int, std::pair<Jet::SRC, int>> &fjMap);
  void addTracks(PHCompositeNode *topNode, std::vector<fastjet::PseudoJet> &particles, const std::vector<PHG4Particlev2*> &TagDecays, std::map<int, std::pair<Jet::SRC, int>> &fjMap);
  void addClusters(PHCompositeNode *topNode, std::vector<fastjet::PseudoJet> &particles, std::map<int, std::pair<Jet::SRC, int>> &fjMap);
  void findMCTaggedJets(PHCompositeNode *topNode);

  bool isAcceptableParticleFlow(ParticleFlowElement *pfPart);
  bool isAcceptableTrack(SvtxTrack *track);
  bool isAcceptableEMCalCluster(CLHEP::Hep3Vector &E_vec_cluster);
  bool isAcceptableHCalCluster(CLHEP::Hep3Vector &E_vec_cluster);
  bool isDecay(HepMC::GenParticle *particle, const std::vector<PHG4Particlev2*> &decays);
  bool isDecay(SvtxTrack *track, const std::vector<PHG4Particlev2*> &decays);
  int createJetNode(PHCompositeNode *topNode);
};

#endif
