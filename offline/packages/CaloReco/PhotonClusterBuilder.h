#ifndef CALORECO_PHOTONCLUSTERBUILDER_H
#define CALORECO_PHOTONCLUSTERBUILDER_H

#include <fun4all/SubsysReco.h>
#include <memory>
#include <string>
#include <vector>

class PHCompositeNode;
class RawClusterContainer;
class RawCluster;
class PhotonClusterv1;
class TowerInfoContainer;
class RawTowerGeomContainer;
class RawTowerGeom;

namespace TMVA
{
  namespace Experimental
  {
    class RBDT;
  }
}  // namespace TMVA

// Simple builder that wraps existing RawClusters above an energy threshold
// into PhotonClusterv1 objects and stores them in a PhotonClusterContainer node.
class PhotonClusterBuilder : public SubsysReco
{
 public:
  explicit PhotonClusterBuilder(const std::string& name = "PhotonClusterBuilder");
  ~PhotonClusterBuilder() override = default;

  int InitRun(PHCompositeNode* topNode) override;
  int process_event(PHCompositeNode* topNode) override;

  void set_input_cluster_node(const std::string& n) { m_input_cluster_node = n; }
  void set_output_photon_node(const std::string& n) { m_output_photon_node = n; }
  void set_ET_threshold(float e) { m_min_cluster_et = e; }
  void set_shower_shape_min_tower_energy(float e) { m_shape_min_tower_E = e; }
  void set_bdt_model_file(const std::string& path) { m_bdt_model_file = path; }
  void set_bdt_feature_list(const std::vector<std::string>& features) { m_bdt_feature_list = features; }
  void set_do_bdt(bool do_bdt) { m_do_bdt = do_bdt; }
  const std::vector<std::string>& get_bdt_feature_list() const { return m_bdt_feature_list; }

 private:
  void CreateNodes(PHCompositeNode* topNode);
  void calculate_shower_shapes(RawCluster* rc, PhotonClusterv1* photon, float eta, float phi);
  void calculate_bdt_score(PhotonClusterv1* photon);
  double getTowerEta(RawTowerGeom* tower_geom, double vx, double vy, double vz);
  std::vector<int> find_closest_hcal_tower(float eta, float phi, RawTowerGeomContainer* geom, TowerInfoContainer* towerContainer, float vertex_z, bool isihcal);
  double deltaR(double eta1, double phi1, double eta2, double phi2);
  bool m_do_bdt{false};

  std::string m_input_cluster_node{"CLUSTERINFO_CEMC"};
  std::string m_output_photon_node{"PHOTONCLUSTER_CEMC"};
  float m_min_cluster_et{5.0f};
  float m_shape_min_tower_E{0.070f};
  std::string m_bdt_model_file{"myBDT_5.root"};
  std::vector<std::string> m_bdt_feature_list;
  float m_vertex{std::numeric_limits<float>::quiet_NaN()};

  RawClusterContainer* m_rawclusters{nullptr};
  RawClusterContainer* m_photon_container{nullptr};
  TowerInfoContainer* m_emc_tower_container{nullptr};
  RawTowerGeomContainer* m_geomEM{nullptr};
  TowerInfoContainer* m_ihcal_tower_container{nullptr};
  RawTowerGeomContainer* m_geomIH{nullptr};
  TowerInfoContainer* m_ohcal_tower_container{nullptr};
  RawTowerGeomContainer* m_geomOH{nullptr};
  std::unique_ptr<TMVA::Experimental::RBDT> m_bdt;
};

#endif  // CALORECO_PHOTONCLUSTERBUILDER_H
