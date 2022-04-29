#ifndef G4EVAL_PHG4DSTCOMPRESSRECO_H
#define G4EVAL_PHG4DSTCOMPRESSRECO_H

#include <fun4all/SubsysReco.h>

#include <set>
#include <string>

class PHCompositeNode;

class PHG4CellContainer;
class PHG4HitContainer;
class SvtxPHG4ParticleMap;
class PHG4ParticleSvtxMap;
class PHG4TruthInfoContainer;

class RawTowerContainer;

class PHG4DstCompressReco : public SubsysReco
{
 public:
  PHG4DstCompressReco(const std::string &name = "PHG4DstCompressReco");
  ~PHG4DstCompressReco() override {}

  //! run initialization
  int InitRun(PHCompositeNode *topNode) override;

  //! event processing
  int process_event(PHCompositeNode *topNode) override;

  void AddHitContainer(const std::string &name) { _compress_g4hit_names.insert(name); }
  void AddCellContainer(const std::string &name) { _compress_g4cell_names.insert(name); }
  void AddTowerContainer(const std::string &name) { _compress_tower_names.insert(name); }
  //! whether to compress tracker hits based on truth association table produced by PHG4DstCompressReco
  void KeepRecoTrackMatchedParticles(bool b = true) { m_keepRecoTrackMatchedParticles = b; }

 private:
  void SearchG4HitNodes(PHCompositeNode *topNode);

  //! whether to compress tracker hits based on truth association table produced by PHG4DstCompressReco
  //! default to false
  bool m_keepRecoTrackMatchedParticles = false;

  PHG4TruthInfoContainer *_truth_info;
  PHG4ParticleSvtxMap *_truthRecoMap = nullptr;
  SvtxPHG4ParticleMap *_recoTruthMap = nullptr;
  std::set<std::string> _compress_g4hit_names;
  std::set<std::string> _compress_g4cell_names;
  std::set<std::string> _compress_tower_names;

  std::set<PHG4CellContainer *> _g4cells;
  std::set<PHG4HitContainer *> _g4hits;
  std::set<PHG4HitContainer *> _keep_g4hits;

  std::set<RawTowerContainer *> _towers;
};

#endif
