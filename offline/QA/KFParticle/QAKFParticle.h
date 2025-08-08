// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef QAKFPARTICLE_H
#define QAKFPARTICLE_H

#include <fun4all/SubsysReco.h>

#include <g4eval/SvtxEvalStack.h>
#include <calotrigger/TriggerAnalyzer.h>

#include <TH1.h>
#include <TH2.h>
#include <memory>
#include <string>  // for string

class KFParticle_Container;
class PHCompositeNode;
// class PHG4Particle;
// class PHG4TruthInfoContainer;
class SvtxClusterEval;
class SvtxTrackMap;
class SvtxTrack;

/*
namespace CLHEP
{
  class HepLorentzVector;
}
*/

class QAKFParticle : public SubsysReco
{
 public:
  QAKFParticle(const std::string &name, const std::string &mother_name, double min_m, double max_m);

  virtual ~QAKFParticle() = default;

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);

  std::string get_histo_prefix();

  void setTrackMapName(const std::string &name) { m_trackMapName = name; }

  void setKFParticleNodeName(const std::string &name) { m_KFParticleNodeName = name; }

 protected:
  // SvtxClusterEval *clustereval = nullptr;
  int m_mother_id = 0;
  double m_min_mass = 0.;
  double m_max_mass = 10.;
  // std::string m_mother_name;

  TH1F *h_mass_KFP = nullptr;
  TH2F *h_mass_KFP_eta = nullptr;
  TH2F *h_mass_KFP_phi = nullptr;
  TH2F *h_mass_KFP_pt = nullptr;
  TH1F *h_mass_KFP_crossing0 = nullptr;
  TH1F *h_mass_KFP_non_crossing0 = nullptr;
  TH1F *h_mass_KFP_ZDC_Coincidence = nullptr;
  TH1F *h_mass_KFP_MBD_NandS_geq_1_vtx_l_30_cm = nullptr;
  TH1F *h_mass_KFP_Jet_6_GeV_MBD_NandS_geq_1_vtx_l_10_cm = nullptr;  

  TriggerAnalyzer *triggeranalyzer{nullptr};

  int m_ZDC_Coincidence_bit = INT_MAX;
  int m_MBD_NandS_geq_1_vtx_l_30_cm_bit = INT_MAX;
  int m_Jet_6_GeV_MBD_NandS_geq_1_vtx_l_10_cm_bit = INT_MAX; 

 private:
  int load_nodes(PHCompositeNode *);

  void initializeTriggerInfo(PHCompositeNode *);

  // SvtxTrack *getTrack(unsigned int track_id, SvtxTrackMap *trackmap);
  // PHG4Particle *getTruthTrack(SvtxTrack *thisTrack);
  // CLHEP::HepLorentzVector *makeHepLV(PHCompositeNode *topNode, int track_number);

  // PHG4TruthInfoContainer *m_truthContainer = nullptr;

  // std::unique_ptr<SvtxEvalStack> m_svtxEvalStack;

  SvtxTrackMap *m_trackMap = nullptr;
  // PHG4TruthInfoContainer *m_truthInfo = nullptr;
  KFParticle_Container *m_kfpContainer = nullptr;
  std::map<std::string, std::pair<int, float>> particleMasses;
  std::string m_trackMapName = "SvtxTrackMap";
  std::string m_KFParticleNodeName = "reconstructedParticles";

  bool hasTriggerInfo = true;
  static const int nTriggerBits = 64;
  int counter = 0;
};

#endif  // QAKFPARTICLE_H
