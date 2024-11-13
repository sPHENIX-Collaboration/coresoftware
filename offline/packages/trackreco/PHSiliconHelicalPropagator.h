#include <fun4all/SubsysReco.h>
#include <trackbase/ActsGeometry.h>
#include <trackbase/TrackFitUtils.h>
#include <trackbase/TrkrClusterCrossingAssoc.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackSeedContainer.h>

class PHSiliconHelicalPropagator : public SubsysReco
{
 public:
  PHSiliconHelicalPropagator(const std::string& name = "PHSiliconHelicalPropagator");
  ~PHSiliconHelicalPropagator();

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;
  int End(PHCompositeNode* topNode) override;

  void set_track_map_name(const std::string& name) { _track_map_name = name; }
  void zeroField() { m_zeroField = true;  }
  void dca_xy_cut(const float cut) { _dca_cut = cut; }
  void dca_z_cut(const float cut) { _dca_z_cut = cut; }

 private:
  int createSeedContainer(TrackSeedContainer*& container, const std::string& container_name, PHCompositeNode* topNode);

  ActsGeometry* _tgeometry = nullptr;
  TrackSeedContainer* _si_seeds = nullptr;
  TrackSeedContainer* _tpc_seeds = nullptr;
  TrackSeedContainer* _svtx_seeds = nullptr;
  TrkrClusterContainer* _cluster_map = nullptr;
  TrkrClusterCrossingAssoc* _cluster_crossing_map = nullptr;

  float _dca_cut = 1.;
  float _dca_z_cut = 1.;
  bool m_zeroField = false;
  std::string _track_map_name = "SvtxTrackSeedContainer";
};
