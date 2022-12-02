#ifndef BUILDRESONANCEJETTAGGINGTREE_H__
#define BUILDRESONANCEJETTAGGINGTREE_H__

#include <fun4all/SubsysReco.h>

#include <g4jets/JetMapv1.h>

#include <HepMC/GenEvent.h>

#include <resonancejettagging/ResonanceJetTagging.h>

#include <vector>

/// Class declarations for use in the analysis module
class PHCompositeNode;
class TFile;
class TTree;
class PHG4Particle;
class KFParticle_Container;

/// Definition of this analysis module class
class BuildResonanceJetTaggingTree : public SubsysReco
{
 public:

  /// Constructor
  BuildResonanceJetTaggingTree(const std::string &name = "BuildResonanceJetTaggingTree", const std::string &fname = "BuildResonanceJetTaggingTree.root", const ResonanceJetTagging::TAG tag = ResonanceJetTagging::TAG::D0);

  // Destructor
  virtual ~BuildResonanceJetTaggingTree();

  /// SubsysReco initialize processing method
  int Init(PHCompositeNode *);

  /// SubsysReco event processing method
  int process_event(PHCompositeNode *);

  /// SubsysReco end processing method
  int End(PHCompositeNode *);
  int loopD0(PHCompositeNode *topNode);
  void findMatchedTruthD0(PHCompositeNode *topNode, Jet *&mcTagJet, HepMC::GenParticle *&mcTag, int decays[]);
  HepMC::GenParticle *getMother(PHCompositeNode *topNode, PHG4Particle *g4daughter);
  bool isReconstructed(int index, std::vector<int> indexRecVector);

  void initializeVariables();
  void initializeTrees();
  void resetTreeVariables();

  void setDoRecunstructed(bool b) { m_dorec = b; }
  bool getDoRecunstructed() { return m_dorec; }
  void setDoTruth(bool b) { m_dotruth = b; }
  bool getDoTruth() { return m_dotruth; }

  void setTagContainerName(const std::string &tagContName) { m_tagcontainer_name = tagContName; }
  std::string getTagContainerName() { return m_tagcontainer_name; }
  void setJetContainerName(const std::string &jetContName) { m_jetcontainer_name = jetContName; }
  std::string getJetContainerName() { return m_jetcontainer_name; }
  void setTruthJetContainerName(const std::string &jetContName) { m_truth_jetcontainer_name = jetContName; }
  std::string getTruthJetContainerName() { return m_truth_jetcontainer_name; }

 private:

  JetMapv1* getJetMapFromNode(PHCompositeNode *topNode, const std::string &name);
  KFParticle_Container* getKFParticleContainerFromNode(PHCompositeNode *topNode, const std::string &name);
  HepMC::GenEvent* getGenEventFromNode(PHCompositeNode *topNode, const std::string &name);
  /// String to contain the outfile name containing the trees
  std::string m_outfilename;
  std::string m_tagcontainer_name;
  std::string m_jetcontainer_name;
  std::string m_truth_jetcontainer_name;
  JetMapv1* m_taggedJetMap;
  JetMapv1* m_truth_taggedJetMap;
  bool m_dorec;
  bool m_dotruth;

  ResonanceJetTagging::TAG m_tag_particle;
  int m_tag_pdg;

  /// TFile to hold the following TTrees and histograms
  TFile *m_outfile = nullptr;
  TTree *m_taggedjettree = nullptr;

  // Tagged-Jet variables
  double m_tagpartpx = NAN;
  double m_tagpartpy = NAN;
  double m_tagpartpz = NAN;
  double m_tagpartpt = NAN;
  double m_tagparteta = NAN;
  double m_tagpartphi = NAN;
  double m_tagpartm = NAN;
  double m_tagjetpx = NAN;
  double m_tagjetpy = NAN;
  double m_tagjetpz = NAN;
  double m_tagjetpt = NAN;
  double m_tagjeteta = NAN;
  double m_tagjetphi = NAN;
  //Truth info
  double m_truth_tagpartpx = NAN;
  double m_truth_tagpartpy = NAN;
  double m_truth_tagpartpz = NAN;
  double m_truth_tagpartpt = NAN;
  double m_truth_tagparteta = NAN;
  double m_truth_tagpartphi = NAN;
  double m_truth_tagjetpx = NAN;
  double m_truth_tagjetpy = NAN;
  double m_truth_tagjetpz = NAN;
  double m_truth_tagjetpt = NAN;
  double m_truth_tagjeteta = NAN;
  double m_truth_tagjetphi = NAN;
};

#endif
