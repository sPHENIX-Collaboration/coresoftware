// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHCOSMICSILICONPROPAGATOR_H
#define PHCOSMICSILICONPROPAGATOR_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class ActsGeometry;
class TrackSeedContainer;
class TrkrClusterContainer;

class PHCosmicSiliconPropagator : public SubsysReco
{
 public:
  PHCosmicSiliconPropagator(const std::string& name = "PHCosmicSiliconPropagator");

  ~PHCosmicSiliconPropagator() override;

  int Init(PHCompositeNode* topNode) override;
  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;
  void set_track_map_name(std::string name) { _track_map_name = name; }
  void set_dca_z_cut(float z) { _dca_z_cut = z; }

 private:
  int createSeedContainer(TrackSeedContainer*& container, const std::string container_name, PHCompositeNode* topNode);

  ActsGeometry* _tgeometry = nullptr;
  TrackSeedContainer* _si_seeds = nullptr;
  TrackSeedContainer* _tpc_seeds = nullptr;
  TrackSeedContainer* _svtx_seeds = nullptr;
  TrkrClusterContainer* _cluster_map = nullptr;

  float _dca_z_cut = 5.;
  std::string _track_map_name = "SvtxTrackSeedContainer";
};

#endif  // PHCOSMICSILICONPROPAGATOR_H
