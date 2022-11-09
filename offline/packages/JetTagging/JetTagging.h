#ifndef JETTAGGING_H__
#define JETTAGGING_H__


#define HomogeneousField


#include <fun4all/SubsysReco.h>

#include <fastjet/ClusterSequence.hh>
#include <fastjet/FunctionOfPseudoJet.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterUtility.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calotrigger/CaloTriggerInfo.h>

#include <HepMC/GenEvent.h>

#include <KFParticle.h>
#include <g4jets/Jetv1.h>
#include <g4jets/JetMapv1.h>

#include <vector>

/// Class declarations for use in the analysis module
class Fun4AllHistoManager;
class PHCompositeNode;
class TFile;
class TTree;
class TH1;
class TH1I;
class PHCompositeNode;
class RawClusterContainer;
class RawCluster;
class SvtxTrackMap;
class GlobalVertex;
class PHHepMCGenEventMap;
class JetRecoEval;
class SvtxTrackEval;
class PHG4TruthInfoContainer;
class PHHepMCGenEvent;
class SvtxTrack;
class PHG4Particle;
class ParticleFlowElement;

/// Definition of this analysis module class
class JetTagging : public SubsysReco
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
    E_SCHEME=0,
    PT_SCHEME=1,
    PT2_SCHEME=2,
    ET_SCHEME=3,
    ET2_SCHEME=4
  };

  /// Constructor
  JetTagging(const std::string &name = "JetTagging",
              const std::string &fname = "JetTagging.root");

  // Destructor
  virtual ~JetTagging();

  /// SubsysReco initialize processing method
  int Init(PHCompositeNode *);

  /// SubsysReco event processing method
  int process_event(PHCompositeNode *);

  /// SubsysReco end processing method
  int End(PHCompositeNode *);

  //Particle flow
  void setParticleFlowEtaAcc(double etamin, double etamax) { m_particleflow_mineta = etamin; m_particleflow_maxeta = etamax; }
  void setParticleFlowMinEta(double etamin) { m_particleflow_mineta = etamin; }
  void setParticleFlowMaxEta(double etamax) { m_particleflow_maxeta = etamax; }
  double getParticleFlowMinEta() { return m_particleflow_mineta; }
  double getParticleFlowMaxEta() { return m_particleflow_maxeta; }

  //tracks
  void setTrackPtAcc(double ptmin, double ptmax) { m_track_minpt = ptmin; m_track_maxpt = ptmax; }
  void setTrackMinPt(double ptmin) { m_track_minpt = ptmin; }
  void setTrackMaxPt(double ptmax) { m_track_maxpt = ptmax; }
  double getTrackMinPt() { return m_track_minpt; }
  double getTrackMaxPt() { return m_track_maxpt; }

  void setTrackEtaAcc(double etamin, double etamax) { m_track_mineta = etamin; m_track_maxeta = etamax; }
  void setTrackMinEta(double etamin) { m_track_mineta = etamin; }
  void setTrackMaxEta(double etamax) { m_track_maxeta = etamax; }
  double getTrackMinEta() { return m_track_mineta; }
  double getTrackMaxEta() { return m_track_maxeta; }

  //EMCal clusters
  void setEMCalClusterPtAcc(double ptmin, double ptmax) { m_EMCal_cluster_minpt = ptmin; m_EMCal_cluster_maxpt = ptmax; }
  void setEMCalClusterMinPt(double ptmin) { m_EMCal_cluster_minpt = ptmin; }
  void setEMCalClusterMaxPt(double ptmax) { m_EMCal_cluster_maxpt = ptmax; }
  double getEMCalClusterMinPt() { return m_EMCal_cluster_minpt; }
  double getEMCalClusterMaxPt() { return m_EMCal_cluster_maxpt; }

  void setEMCalClusterEtaAcc(double etamin, double etamax) { m_EMCal_cluster_mineta = etamin; m_EMCal_cluster_maxeta = etamax; }
  void setEMCalClusterMinEta(double etamin) { m_EMCal_cluster_mineta = etamin; }
  void setEMCalClusterMaxEta(double etamax) { m_EMCal_cluster_maxeta = etamax; }
  double getEMCalClusterMinEta() { return m_EMCal_cluster_mineta; }
  double getEMCalClusterMaxEta() { return m_EMCal_cluster_maxeta; }

  //HCal clusters
  void setHCalClusterPtAcc(double ptmin, double ptmax) { m_HCal_cluster_minpt = ptmin; m_HCal_cluster_maxpt = ptmax; }
  void setHCalClusterMinPt(double ptmin) { m_HCal_cluster_minpt = ptmin; }
  void setHCalClusterMaxPt(double ptmax) { m_HCal_cluster_maxpt = ptmax; }
  double getHCalClusterMinPt() { return m_HCal_cluster_minpt; }
  double getHCalClusterMaxPt() { return m_HCal_cluster_maxpt; }

  void setHCalClusterEtaAcc(double etamin, double etamax) { m_HCal_cluster_mineta = etamin; m_HCal_cluster_maxeta = etamax; }
  void setHCalClusterMinEta(double etamin) { m_HCal_cluster_mineta = etamin; }
  void setHCalClusterMaxEta(double etamax) { m_HCal_cluster_maxeta = etamax; }
  double getHCalClusterMinEta() { return m_HCal_cluster_mineta; }
  double getHCalClusterMaxEta() { return m_HCal_cluster_maxeta; }

  //Set/Get add Tracks and Clusters
  void setAddParticleFlow(bool b) { m_add_particleflow = b; }
  bool getAddParticleFlow() { return m_add_particleflow; }
  void setAddTracks(bool b) { m_add_tracks = b; }
  bool getAddTracks() { return m_add_tracks; }
  void setAddEMCalClusters(bool b) { m_add_EMCal_clusters = b; }
  bool getAddEMCalClusters() { return m_add_EMCal_clusters; }
  void setAddHCalClusters(bool b) { m_add_HCal_clusters = b; }
  bool getAddHCalClusters() { return m_add_HCal_clusters; }

  //Jet settings
  void setR(double r) { m_jetr = r; }
  double getR(double r) { return m_jetr; }
  void setJetAlgo(ALGO jetalgo) {
    switch(jetalgo)
    {
      case ALGO::ANTIKT:
        m_jetalgo = fastjet::antikt_algorithm;
      case ALGO::KT:
        m_jetalgo = fastjet::kt_algorithm;
      case ALGO::CAMBRIDGE:
        m_jetalgo = fastjet::cambridge_algorithm;
    }
  }
  fastjet::JetAlgorithm getJetAlgo() { return m_jetalgo; }
  void setRecombScheme(RECOMB recomb_scheme) {
    switch(recomb_scheme)
    {
      case RECOMB::E_SCHEME:
        m_recomb_scheme = fastjet::E_scheme;
      case RECOMB::PT_SCHEME:
        m_recomb_scheme = fastjet::pt_scheme;
      case RECOMB::PT2_SCHEME:
        m_recomb_scheme = fastjet::pt2_scheme;
      case RECOMB::ET_SCHEME:
        m_recomb_scheme = fastjet::Et_scheme;
      case RECOMB::ET2_SCHEME:
        m_recomb_scheme = fastjet::Et2_scheme;
    }
  }
  fastjet::RecombinationScheme getRecombScheme() { return m_recomb_scheme; }
  void setJetParameters(double r, ALGO jetalgo, RECOMB recomb_scheme)
  {
    setR(r);
    setJetAlgo(jetalgo);
    setRecombScheme(recomb_scheme);
  }

  void setMakeQualityPlots(bool q) { m_qualy_plots = q; }
  bool getMakeQualityPlots() { return m_qualy_plots; }

  void setJetContainerName(std::string n) {m_jetcontainer_name = n;}
  std::string getJetContainerName() {return m_jetcontainer_name;}
  void setSaveDST(bool s) { m_save_dst = s; }
  bool getSaveDST() { return m_save_dst; }
  void setIsMC(bool b) { m_ismc = b; }
  bool getIsMC() { return m_ismc; }
  void setSaveDSTMC(bool s) { m_save_truth_dst = s; }
  bool getSaveDSTMC() { return m_save_truth_dst; }

 private:
  /// String to contain the outfile name containing the trees
  std::string m_outfilename;
  std::string m_jetcontainer_name;

  /// Fun4All Histogram Manager tool
  Fun4AllHistoManager *m_hm;

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

  JetMapv1* m_taggedJetMap;
  JetMapv1* m_truth_taggedJetMap;

  /// TFile to hold the following TTrees and histograms
  TFile *m_outfile;
  TTree *m_taggedjettree;
  TH1 *m_eventcount_h;
  TH1 *m_rec_tagpart_pt;
  TH1 *m_gen_withrec_tagpart_pt; //Distribution of generated particles with a match to a rec particle
  TH1 *m_gen_tagpart_pt;
  TH1 *m_rec_tracks_pt;
  TH1 *m_gen_tracks_pt;
  TH1 *m_rec_emcal_clusters_pt;
  TH1 *m_rec_hcalin_clusters_pt;
  TH1 *m_rec_hcalout_clusters_pt;

  int m_tag_pdg;
  bool m_qualy_plots;
  bool m_save_dst;
  bool m_save_truth_dst;
  unsigned int m_jet_id = 0;
  unsigned int m_truth_jet_id = 0;
  bool m_ismc;

  /// Methods for grabbing the data
  void findTaggedJets(PHCompositeNode *topNode, KFParticle *Tag, KFParticle *TagDecays[], int nDecays);
  void addParticleFlow(PHCompositeNode *topNode, std::vector<fastjet::PseudoJet> &particles, KFParticle *TagDecays[], int nDecays, std::map<int, std::pair<Jet::SRC, int>> &fjMap);
  void addTracks(PHCompositeNode *topNode, std::vector<fastjet::PseudoJet> &particles, KFParticle *TagDecays[], int nDecays, std::map<int, std::pair<Jet::SRC, int>> &fjMap);
  void addClusters(PHCompositeNode *topNode, std::vector<fastjet::PseudoJet> &particles, std::map<int, std::pair<Jet::SRC, int>> &fjMap);
  void getTracks(PHCompositeNode *topNode);
  HepMC::GenParticle* findMCTaggedJets(PHCompositeNode *topNode, KFParticle *decays[], int nDecays);
  HepMC::GenParticle* findMCTag(PHCompositeNode *topNode, KFParticle *decays[], int nDecays, PHG4Particle *mcDaughters[]);
  HepMC::GenParticle* getMother(PHCompositeNode *topNode, PHG4Particle *g4daughter);
  //bool hasMCTagParent(PHG4Particle *g4particle, PHG4TruthInfoContainer *truthinfo, int &parent_id);
  void findNonRecMC(PHCompositeNode *topNode, std::vector<HepMC::GenParticle*> mcTags);
  void doMCLoop(PHCompositeNode *topNode);

  bool isAcceptableParticleFlow(ParticleFlowElement* pfPart);
  bool isAcceptableTrack(SvtxTrack *track);
  bool isAcceptableEMCalCluster(CLHEP::Hep3Vector &E_vec_cluster);
  bool isAcceptableHCalCluster(CLHEP::Hep3Vector &E_vec_cluster);
  bool isDecay(HepMC::GenParticle *particle, PHG4Particle *decays[], int nDecays);
  bool isSameParticle(SvtxTrack *track, KFParticle *particle);
  bool isDecay(SvtxTrack *track, KFParticle *decays[], int nDecays);
  void initializeVariables();
  void initializeTrees();
  int createJetNode(PHCompositeNode* topNode);
  void resetTreeVariables();

  /**
   * Make variables for the relevant trees
   */

  // Tagged-Jet variables
  double m_tagpartpx = -9999.;
  double m_tagpartpy = -9999.;
  double m_tagpartpz = -9999.;
  double m_tagpartpt = -9999.;
  double m_tagparteta = -9999.;
  double m_tagpartphi = -9999.;
  double m_tagpartm = -9999.;
  double m_tagjetpx = -9999.;
  double m_tagjetpy = -9999.;
  double m_tagjetpz = -9999.;
  double m_tagjetpt = -9999.;
  double m_tagjeteta = -9999.;
  double m_tagjetphi = -9999.;
  //Truth info
  double m_truth_tagpartpx = -9999.;
  double m_truth_tagpartpy = -9999.;
  double m_truth_tagpartpz = -9999.;
  double m_truth_tagpartpt = -9999.;
  double m_truth_tagparteta = -9999.;
  double m_truth_tagpartphi = -9999.;
  double m_truth_tagjetpx = -9999.;
  double m_truth_tagjetpy = -9999.;
  double m_truth_tagjetpz = -9999.;
  double m_truth_tagjetpt = -9999.;
  double m_truth_tagjeteta = -9999.;
  double m_truth_tagjetphi = -9999.;
};

#endif
