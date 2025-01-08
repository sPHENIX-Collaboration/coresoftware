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
  void set_track_map_name(const std::string& name) { _track_map_name = name; }
  void set_dca_z_cut(float z) { _dca_z_cut = z; }
  void set_dca_xy_cut(float xy) { _dca_xy_cut = xy; }
  void zero_field() { m_zeroField = true; }
  void resetSvtxSeedContainer() { m_resetContainer = true; }

 private:
  int createSeedContainer(TrackSeedContainer*& container, const std::string& container_name, PHCompositeNode* topNode);

  ActsGeometry* _tgeometry = nullptr;
  TrackSeedContainer* _si_seeds = nullptr;
  TrackSeedContainer* _tpc_seeds = nullptr;
  TrackSeedContainer* _svtx_seeds = nullptr;
  TrkrClusterContainer* _cluster_map = nullptr;

  float m_resetContainer = false;
  float _dca_z_cut = 5.;
  float _dca_xy_cut = 5.;
  bool m_zeroField = false;
  std::string _track_map_name = "SvtxTrackSeedContainer";
};

#endif  // PHCOSMICSILICONPROPAGATOR_H
