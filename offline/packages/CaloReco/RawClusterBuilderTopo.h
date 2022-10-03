#ifndef CALORECO_RAWCLUSTERTOPO_H
#define CALORECO_RAWCLUSTERTOPO_H

//===========================================================
/// \file RawClusterBuilderTopo.h
/// \brief 3-D topoClustering across calorimeter layers
/// \author Dennis V. Perepelitsa
//===========================================================

#include <fun4all/SubsysReco.h>

#include <map>
#include <string>
#include <utility>  // for pair
#include <vector>

class PHCompositeNode;
class RawClusterContainer;
class RawTowerGeomContainer;

class RawClusterBuilderTopo : public SubsysReco
{
 public:
  explicit RawClusterBuilderTopo(const std::string &name = "RawClusterBuilderTopo");
  ~RawClusterBuilderTopo() override {}

  int InitRun(PHCompositeNode *topNode) override;
  int process_event(PHCompositeNode *topNode) override;
  int End(PHCompositeNode *topNode) override;

  void set_nodename(const std::string &nodename)
  {
    ClusterNodeName = nodename;
  }

  void set_noise(float noise_0 = 0.0025, float noise_1 = 0.006, float noise_2 = 0.03)
  {
    _noise_LAYER[0] = noise_0;
    _noise_LAYER[1] = noise_1;
    _noise_LAYER[2] = noise_2;
  }

  void set_significance(float seed, float grow, float peri)
  {
    _sigma_seed = seed;
    _sigma_grow = grow;
    _sigma_peri = peri;
  }

  void allow_corner_neighbor(bool allow)
  {
    _allow_corner_neighbor = allow;
  }

  void set_enable_HCal(bool enable_HCal)
  {
    _enable_HCal = enable_HCal;
  }

  void set_enable_EMCal(bool enable_EMCal)
  {
    _enable_EMCal = enable_EMCal;
  }

  void set_do_split(bool do_split)
  {
    _do_split = do_split;
  }

  void set_minE_local_max(float minE_0 = 1, float minE_1 = 1, float minE_2 = 1)
  {
    _local_max_minE_LAYER[0] = minE_0;
    _local_max_minE_LAYER[1] = minE_1;
    _local_max_minE_LAYER[2] = minE_2;
  }

  void set_R_shower(float R_shower)
  {
    _R_shower = R_shower;
  }

 private:
  void CreateNodes(PHCompositeNode *topNode);

  std::vector<std::vector<std::vector<float> > > _TOWERMAP_E_LAYER_ETA_PHI;
  std::vector<std::vector<std::vector<int> > > _TOWERMAP_KEY_LAYER_ETA_PHI;
  std::vector<std::vector<std::vector<int> > > _TOWERMAP_STATUS_LAYER_ETA_PHI;

  std::vector<std::vector<float> > _EMTOWERMAP_E_ETA_PHI;
  std::vector<std::vector<int> > _EMTOWERMAP_KEY_ETA_PHI;
  std::vector<std::vector<int> > _EMTOWERMAP_STATUS_ETA_PHI;

  // geometric constants to express IHCal<->EMCal overlap in eta
  static int RawClusterBuilderTopo_constants_EMCal_eta_start_given_IHCal[];

  static int RawClusterBuilderTopo_constants_EMCal_eta_end_given_IHCal[];

  static int RawClusterBuilderTopo_constants_IHCal_eta_given_EMCal[];

  // geometric parameters defined at runtime
  int _EMCAL_NETA;
  int _EMCAL_NPHI;

  int _HCAL_NETA;
  int _HCAL_NPHI;

  // utility functions to express IHCal<->EMCal overlap in phi
  int get_first_matching_EMCal_phi_from_IHCal(int index_hcal_phi)
  {
    return ((68 + 4 * (index_hcal_phi - 32) + _EMCAL_NPHI) % _EMCAL_NPHI);
  }

  int get_matching_HCal_phi_from_EMCal(int index_emcal_phi)
  {
    return ((32 + (index_emcal_phi - 68 + _EMCAL_NPHI) / 4) % _HCAL_NPHI);
  }

  std::vector<int> get_adjacent_towers_by_ID(int ID);

  float calculate_dR(float, float, float, float);

  void export_single_cluster(const std::vector<int> &);

  void export_clusters(const std::vector<int> &, std::map<int, std::pair<int, int> >, unsigned int, std::vector<float>, std::vector<float>, std::vector<float>);

  int get_ID(int ilayer, int ieta, int iphi)
  {
    if (ilayer < 2)
    {
      return ilayer * _HCAL_NETA * _HCAL_NPHI + ieta * _HCAL_NPHI + iphi;
    }
    return _EMCAL_NPHI * _EMCAL_NETA + ieta * _EMCAL_NPHI + iphi;
  }

  int get_ilayer_from_ID(int ID)
  {
    if (ID < _EMCAL_NPHI * _EMCAL_NETA)
    {
      return ((int) (ID / (_HCAL_NETA * _HCAL_NPHI)));
    }
    else
    {
      return 2;
    }
  }

  int get_ieta_from_ID(int ID)
  {
    if (ID < _EMCAL_NPHI * _EMCAL_NETA)
    {
      return ((int) ((ID % (_HCAL_NETA * _HCAL_NPHI)) / (_HCAL_NPHI)));
    }
    else
    {
      return ((int) ((ID - _EMCAL_NPHI * _EMCAL_NETA) / _EMCAL_NPHI));
    }
  }

  int get_iphi_from_ID(int ID)
  {
    if (ID < _EMCAL_NPHI * _EMCAL_NETA)
    {
      return ((int) (ID % _HCAL_NPHI));
    }
    else
    {
      return ((int) ((ID - _EMCAL_NPHI * _EMCAL_NETA) % _EMCAL_NPHI));
    }
  }

  int get_status_from_ID(int ID)
  {
    if (ID < _EMCAL_NPHI * _EMCAL_NETA)
    {
      return _TOWERMAP_STATUS_LAYER_ETA_PHI[get_ilayer_from_ID(ID)][get_ieta_from_ID(ID)][get_iphi_from_ID(ID)];
    }
    return _EMTOWERMAP_STATUS_ETA_PHI[get_ieta_from_ID(ID)][get_iphi_from_ID(ID)];
  }

  float get_E_from_ID(int ID)
  {
    if (ID < _EMCAL_NPHI * _EMCAL_NETA)
    {
      return _TOWERMAP_E_LAYER_ETA_PHI[get_ilayer_from_ID(ID)][get_ieta_from_ID(ID)][get_iphi_from_ID(ID)];
    }
    return _EMTOWERMAP_E_ETA_PHI[get_ieta_from_ID(ID)][get_iphi_from_ID(ID)];
  }

  void set_status_by_ID(int ID, int status)
  {
    if (ID < _EMCAL_NPHI * _EMCAL_NETA)
    {
      _TOWERMAP_STATUS_LAYER_ETA_PHI[get_ilayer_from_ID(ID)][get_ieta_from_ID(ID)][get_iphi_from_ID(ID)] = status;
    }
    else
    {
      _EMTOWERMAP_STATUS_ETA_PHI[get_ieta_from_ID(ID)][get_iphi_from_ID(ID)] = status;
    }
  }

  RawClusterContainer *_clusters = nullptr;

  RawTowerGeomContainer *_geom_containers[3];

  float _noise_LAYER[3];

  float _sigma_seed;
  float _sigma_grow;
  float _sigma_peri;

  bool _allow_corner_neighbor;

  bool _enable_HCal;
  bool _enable_EMCal;

  bool _do_split;
  float _local_max_minE_LAYER[3];
  float _R_shower;

  std::string ClusterNodeName;
};

#endif
