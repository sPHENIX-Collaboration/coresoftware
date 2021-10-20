// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef QAG4SIMULATIONKFPARTICLE_H
#define QAG4SIMULATIONKFPARTICLE_H

#include <fun4all/SubsysReco.h>

#include <g4eval/SvtxEvalStack.h>

#include <memory>
#include <string>  // for string

class KFParticle_Container;
class PHCompositeNode;
class PHG4Particle;
class PHG4TruthInfoContainer;
class SvtxClusterEval;
class SvtxTrackMap;
class SvtxTrack;

namespace CLHEP
{
  class HepLorentzVector;
}

class QAG4SimulationKFParticle : public SubsysReco
{
 public:
  QAG4SimulationKFParticle(const std::string &name, const std::string &mother_name, double min_m, double max_m);

  virtual ~QAG4SimulationKFParticle() = default;

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);

  std::string get_histo_prefix();

  void setTrackMapName(const std::string &name) { m_trackMapName = name; }

 protected:
  SvtxClusterEval *clustereval = nullptr;
  int m_mother_id = 0;
  double m_min_mass = 0.;
  double m_max_mass = 10.;
  std::string m_mother_name;

 private:
  int load_nodes(PHCompositeNode *);

  SvtxTrack *getTrack(unsigned int track_id, SvtxTrackMap *trackmap);
  PHG4Particle *getTruthTrack(SvtxTrack *thisTrack);
  CLHEP::HepLorentzVector *makeHepLV(PHCompositeNode *topNode, int track_number);

  PHG4TruthInfoContainer *m_truthContainer = nullptr;

  std::unique_ptr<SvtxEvalStack> m_svtxEvalStack;

  SvtxTrackMap *m_trackMap = nullptr;
  PHG4TruthInfoContainer *m_truthInfo = nullptr;
  KFParticle_Container *m_kfpContainer = nullptr;
  std::map<std::string, std::pair<int, float>> particleMasses;
  std::string m_trackMapName = "SvtxTrackMap";
};

#endif  // QAG4SIMULATIONKFPARTICLE_H
