#include "RawClusterBuilderTopo.h"

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawClusterv1.h>
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>  // for abs
#include <exception>
#include <iostream>
#include <list>
#include <memory>  // for allocator_traits<>::valu...
#include <stdexcept>
#include <utility>
#include <vector>

bool sort_by_pair_second(const std::pair<int, float> &a, const std::pair<int, float> &b)
{
  return (a.second > b.second);
}

int RawClusterBuilderTopo::RawClusterBuilderTopo_constants_EMCal_eta_start_given_IHCal[24] =
    {2, 6, 10, 14, 18, 22, 26, 30, 33, 37, 41, 44,
     48, 52, 55, 59, 63, 66, 70, 74, 78, 82, 86, 90};

int RawClusterBuilderTopo::RawClusterBuilderTopo_constants_EMCal_eta_end_given_IHCal[24] =
    {5, 9, 13, 17, 21, 25, 29, 32, 36, 40, 43, 47,
     51, 54, 58, 62, 65, 69, 73, 77, 81, 85, 89, 93};

int RawClusterBuilderTopo::RawClusterBuilderTopo_constants_IHCal_eta_given_EMCal[96] = {
    -1, -1, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2,
    2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5,
    5, 5, 6, 6, 6, 6, 7, 7, 7, 8, 8, 8,
    8, 9, 9, 9, 9, 10, 10, 10, 11, 11, 11, 11,
    12, 12, 12, 12, 13, 13, 13, 14, 14, 14, 14, 15,
    15, 15, 15, 16, 16, 16, 17, 17, 17, 17, 18, 18,
    18, 18, 19, 19, 19, 19, 20, 20, 20, 20, 21, 21,
    21, 21, 22, 22, 22, 22, 23, 23, 23, 23, -1, -1};

float RawClusterBuilderTopo::calculate_dR(float eta1, float eta2, float phi1, float phi2)
{
  float deta = eta1 - eta2;
  float dphi = phi1 - phi2;
  while (dphi > M_PI) dphi -= 2 * M_PI;
  while (dphi < -M_PI) dphi += 2 * M_PI;
  return sqrt(pow(deta, 2) + pow(dphi, 2));
}

std::vector<int> RawClusterBuilderTopo::get_adjacent_towers_by_ID(int ID)
{
  int this_layer = get_ilayer_from_ID(ID);
  int this_eta = get_ieta_from_ID(ID);
  int this_phi = get_iphi_from_ID(ID);

  std::vector<int> adjacent_towers;

  // for both IHCal and OHCal, add adjacent layers in the HCal
  if (this_layer == 0 || this_layer == 1)
  {
    for (int delta_layer = 0; delta_layer <= 1; delta_layer++)
    {
      for (int delta_eta = -1; delta_eta <= 1; delta_eta++)
      {
        for (int delta_phi = -1; delta_phi <= 1; delta_phi++)
        {
          if (delta_layer == 0 && delta_eta == 0 && delta_phi == 0) continue;  // this is the same tower

          int test_eta = this_eta + delta_eta;
          if (test_eta < 0 || test_eta >= _HCAL_NETA)
          {
            continue;
          }  // ignore if at the (eta) edge of calorimeter

          int test_layer = (this_layer + delta_layer) % 2;                  // wrap around in layer
          int test_phi = (this_phi + delta_phi + _HCAL_NPHI) % _HCAL_NPHI;  // wrap around in phi (add 64 to avoid -1)

          // disallow "corner" adjacency (diagonal in eta/phi plane and in different layer) if this option not enabled
          if (!_allow_corner_neighbor && delta_layer == 1 && abs(delta_phi) == 1 && abs(delta_eta) == 1)
          {
            if (Verbosity() > 20) std::cout << "RawClusterBuilderTopo::get_adjacent_towers_by_ID : corner growth not allowed " << std::endl;
            continue;
          }

          // add to list of adjacent towers
          adjacent_towers.push_back(get_ID(test_layer, test_eta, test_phi));
        }
      }
    }
  }

  // for IHCal only, also add 4x4 group of EMCal towers
  if (this_layer == 0 && _enable_EMCal)
  {
    int EMCal_phi_start = get_first_matching_EMCal_phi_from_IHCal(this_phi);
    int EMCal_eta_start = RawClusterBuilderTopo_constants_EMCal_eta_start_given_IHCal[this_eta];
    int EMCal_eta_end = RawClusterBuilderTopo_constants_EMCal_eta_end_given_IHCal[this_eta];

    for (int new_eta = EMCal_eta_start; new_eta <= EMCal_eta_end; new_eta++)
    {
      for (int delta_phi = 0; delta_phi < 4; delta_phi++)
      {
        int new_phi = (EMCal_phi_start + delta_phi + _EMCAL_NPHI) % _EMCAL_NPHI;

        int EMCal_tower = get_ID(2, new_eta, new_phi);
        if (Verbosity() > 20) std::cout << "RawClusterBuilderTopo::get_adjacent_towers_by_ID : HCal tower with eta / phi = " << this_eta << " / " << this_phi << ", adding EMCal tower with eta / phi = " << new_eta << " / " << new_phi << std::endl;
        adjacent_towers.push_back(EMCal_tower);
      }
    }
  }

  // for EMCal, add adjacent EMCal towers and (usually) one IHCal tower
  if (this_layer == 2)
  {
    // EMCal towers first
    for (int delta_eta = -1; delta_eta <= 1; delta_eta++)
    {
      for (int delta_phi = -1; delta_phi <= 1; delta_phi++)
      {
        if (delta_eta == 0 && delta_phi == 0) continue;  // this is the same tower

        int test_eta = this_eta + delta_eta;
        if (test_eta < 0 || test_eta >= _EMCAL_NETA)
        {
          continue;
        }  // ignore if at the (eta) edge of calorimeter

        int test_phi = (this_phi + delta_phi + _EMCAL_NPHI) % _EMCAL_NPHI;  // wrap around in phi (add 256 to avoid -1)

        // add to list of adjacent towers
        adjacent_towers.push_back(get_ID(this_layer, test_eta, test_phi));
      }
    }

    // now add IHCal towers
    if (_enable_HCal)
    {
      int HCal_eta = RawClusterBuilderTopo_constants_IHCal_eta_given_EMCal[this_eta];
      int HCal_phi = get_matching_HCal_phi_from_EMCal(this_phi);

      if (HCal_eta >= 0)
      {
        int IHCal_tower = get_ID(0, HCal_eta, HCal_phi);
        if (Verbosity() > 20) std::cout << "RawClusterBuilderTopo::get_adjacent_towers_by_ID : EMCal tower with eta / phi = " << this_eta << " / " << this_phi << ", adding IHCal tower with eta / phi = " << HCal_eta << " / " << HCal_phi << std::endl;
        adjacent_towers.push_back(IHCal_tower);
      }
      else
      {
        if (Verbosity() > 20) std::cout << "RawClusterBuilderTopo::get_adjacent_towers_by_ID : EMCal tower with eta / phi = " << this_eta << " / " << this_phi << ", does not have matching IHCal due to large eta " << std::endl;
      }
    }
  }

  return adjacent_towers;
}

void RawClusterBuilderTopo::export_single_cluster(const std::vector<int> &original_towers)
{
  if (Verbosity() > 2)
  {
    std::cout << "RawClusterBuilderTopo::export_single_cluster called " << std::endl;
  }

  std::map<int, std::pair<int, int> > tower_ownership;
  for (const int &original_tower : original_towers)
  {
    tower_ownership[original_tower] = std::pair<int, int>(0, -1);  // all towers owned by cluster 0
  }
  export_clusters(original_towers, tower_ownership, 1, std::vector<float>(), std::vector<float>(), std::vector<float>());

  return;
}

void RawClusterBuilderTopo::export_clusters(const std::vector<int> &original_towers, std::map<int, std::pair<int, int> > tower_ownership, unsigned int n_clusters, std::vector<float> pseudocluster_sumE, std::vector<float> pseudocluster_eta, std::vector<float> pseudocluster_phi)
{
  if (n_clusters != 1)  // if we didn't just pass down from export_single_cluster
  {
    if (Verbosity() > 2)
      std::cout << "RawClusterBuilderTopo::export_clusters called on an initial cluster with " << n_clusters << " final clusters " << std::endl;
  }
  // build a RawCluster for output
  std::vector<RawCluster *> clusters;
  std::vector<float> clusters_E;
  std::vector<float> clusters_x;
  std::vector<float> clusters_y;
  std::vector<float> clusters_z;

  for (unsigned int pc = 0; pc < n_clusters; pc++)
  {
    clusters.push_back(new RawClusterv1());
    clusters_E.push_back(0);
    clusters_x.push_back(0);
    clusters_y.push_back(0);
    clusters_z.push_back(0);
  }

  for (int original_tower : original_towers)
  {
    int this_ID = original_tower;
    std::pair<int, int> the_pair = tower_ownership[this_ID];

    if (Verbosity() > 5)
    {
      std::cout << "RawClusterBuilderTopo::export_clusters -> assigning tower " << original_tower << " with ownership ( " << the_pair.first << ", " << the_pair.second << " ) " << std::endl;
    }
    int this_ieta = get_ieta_from_ID(this_ID);
    int this_iphi = get_iphi_from_ID(this_ID);
    int this_layer = get_ilayer_from_ID(this_ID);
    float this_E = get_E_from_ID(this_ID);

    int this_key = 0;
    if (this_layer == 2)
    {
      this_key = _EMTOWERMAP_KEY_ETA_PHI[this_ieta][this_iphi];
    }
    else
    {
      this_key = _TOWERMAP_KEY_LAYER_ETA_PHI[this_layer][this_ieta][this_iphi];
    }

    RawTowerGeom *tower_geom = _geom_containers[this_layer]->get_tower_geometry(this_key);

    if (the_pair.second == -1)
    {
      // assigned only to one cluster, easy
      clusters[the_pair.first]->addTower(this_key, this_E);
      clusters_E[the_pair.first] = clusters_E[the_pair.first] + this_E;
      clusters_x[the_pair.first] = clusters_x[the_pair.first] + this_E * tower_geom->get_center_x();
      clusters_y[the_pair.first] = clusters_y[the_pair.first] + this_E * tower_geom->get_center_y();
      clusters_z[the_pair.first] = clusters_z[the_pair.first] + this_E * tower_geom->get_center_z();

      if (Verbosity() > 5)
      {
        std::cout << " -> tower ID " << this_ID << " fully assigned to pseudocluster " << the_pair.first << std::endl;
      }
    }
    else
    {
      // assigned to two clusters! get energy sharing fraction ...
      float dR1 = calculate_dR(tower_geom->get_eta(), pseudocluster_eta[the_pair.first], tower_geom->get_phi(), pseudocluster_phi[the_pair.first]) / _R_shower;
      float dR2 = calculate_dR(tower_geom->get_eta(), pseudocluster_eta[the_pair.second], tower_geom->get_phi(), pseudocluster_phi[the_pair.second]) / _R_shower;
      float r = std::exp(dR1 - dR2);
      float frac1 = pseudocluster_sumE[the_pair.first] / (pseudocluster_sumE[the_pair.first] + r * pseudocluster_sumE[the_pair.second]);

      if (Verbosity() > 5)
      {
        std::cout << " tower ID " << this_ID << " has dR1 = " << dR1 << " to pseudocluster " << the_pair.first << " , and dR2 = " << dR2 << " to pseudocluster " << the_pair.second << ", so frac1 = " << frac1 << std::endl;
      }
      clusters[the_pair.first]->addTower(this_key, this_E * frac1);
      clusters_E[the_pair.first] = clusters_E[the_pair.first] + this_E * frac1;
      clusters_x[the_pair.first] = clusters_x[the_pair.first] + this_E * tower_geom->get_center_x() * frac1;
      clusters_y[the_pair.first] = clusters_y[the_pair.first] + this_E * tower_geom->get_center_y() * frac1;
      clusters_z[the_pair.first] = clusters_z[the_pair.first] + this_E * tower_geom->get_center_z() * frac1;

      clusters[the_pair.second]->addTower(this_key, this_E * (1 - frac1));
      clusters_E[the_pair.second] = clusters_E[the_pair.second] + this_E * (1 - frac1);
      clusters_x[the_pair.second] = clusters_x[the_pair.second] + this_E * tower_geom->get_center_x() * (1 - frac1);
      clusters_y[the_pair.second] = clusters_y[the_pair.second] + this_E * tower_geom->get_center_y() * (1 - frac1);
      clusters_z[the_pair.second] = clusters_z[the_pair.second] + this_E * tower_geom->get_center_z() * (1 - frac1);
    }
  }

  // iterate through and add to official container

  for (unsigned int cl = 0; cl < n_clusters; cl++)
  {
    clusters[cl]->set_energy(clusters_E[cl]);

    float mean_x = clusters_x[cl] / clusters_E[cl];
    float mean_y = clusters_y[cl] / clusters_E[cl];
    float mean_z = clusters_z[cl] / clusters_E[cl];

    clusters[cl]->set_r(std::sqrt(mean_y * mean_y + mean_x * mean_x));
    clusters[cl]->set_phi(std::atan2(mean_y, mean_x));
    clusters[cl]->set_z(mean_z);

    _clusters->AddCluster(clusters[cl]);

    if (Verbosity() > 1)
    {
      std::cout << "RawClusterBuilderTopo::export_clusters: added cluster with E = " << clusters_E[cl] << ", eta = " << -1 * log(tan(std::atan2(std::sqrt(mean_y * mean_y + mean_x * mean_x), mean_z) / 2.0)) << ", phi = " << std::atan2(mean_y, mean_x) << std::endl;
    }
  }

  return;
}

RawClusterBuilderTopo::RawClusterBuilderTopo(const std::string &name)
  : SubsysReco(name)
{
  // geometry defined at run-time
  _EMCAL_NETA = -1;
  _EMCAL_NPHI = -1;

  _HCAL_NETA = -1;
  _HCAL_NPHI = -1;
  std::fill(std::begin(_geom_containers), std::end(_geom_containers), nullptr);
  _noise_LAYER[0] = 0.0025;
  _noise_LAYER[1] = 0.006;
  _noise_LAYER[2] = 0.03;  // EM

  _sigma_seed = 4.0;
  _sigma_grow = 2.0;
  _sigma_peri = 0.0;

  _allow_corner_neighbor = true;

  _enable_HCal = true;
  _enable_EMCal = true;

  _do_split = true;
  _R_shower = 0.025;

  _local_max_minE_LAYER[0] = 1;
  _local_max_minE_LAYER[1] = 1;
  _local_max_minE_LAYER[2] = 1;

  ClusterNodeName = "TOPOCLUSTER_HCAL";
}

int RawClusterBuilderTopo::InitRun(PHCompositeNode *topNode)
{
  try
  {
    CreateNodes(topNode);
  }
  catch (std::exception &e)
  {
    std::cout << PHWHERE << ": " << e.what() << std::endl;
    throw;
  }

  if (Verbosity() > 0)
  {
    std::cout << "RawClusterBuilderTopo::InitRun: initialized with EMCal enable = " << _enable_EMCal << " and I+OHCal enable = " << _enable_HCal << std::endl;
    std::cout << "RawClusterBuilderTopo::InitRun: initialized with sigma_noise in EMCal / IHCal / OHCal = " << _noise_LAYER[2] << " / " << _noise_LAYER[0] << " / " << _noise_LAYER[1] << std::endl;
    std::cout << "RawClusterBuilderTopo::InitRun: initialized with noise multiples for seeding / growth / perimeter ( S / N / P ) = " << _sigma_seed << " / " << _sigma_grow << " / " << _sigma_peri << std::endl;
    std::cout << "RawClusterBuilderTopo::InitRun: initialized with allow_corner_neighbor = " << _allow_corner_neighbor << " (in HCal)" << std::endl;
    std::cout << "RawClusterBuilderTopo::InitRun: initialized with do_split = " << _do_split << " , R_shower = " << _R_shower << " (angular units) " << std::endl;
    std::cout << "RawClusterBuilderTopo::InitRun: initialized with minE for local max in EMCal / IHCal / OHCal = " << _local_max_minE_LAYER[2] << " / " << _local_max_minE_LAYER[0] << " / " << _local_max_minE_LAYER[1] << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int RawClusterBuilderTopo::process_event(PHCompositeNode *topNode)
{
  RawTowerContainer *towersEM = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_CEMC");
  RawTowerContainer *towersIH = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALIN");
  RawTowerContainer *towersOH = findNode::getClass<RawTowerContainer>(topNode, "TOWER_CALIB_HCALOUT");

  if (!towersEM)
  {
    std::cout << " RawClusterBuilderTopo::process_event : container TOWER_CALIB_CEMC does not exist, aborting " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if (!towersIH)
  {
    std::cout << " RawClusterBuilderTopo::process_event : container TOWER_CALIB_HCALIN does not exist, aborting " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if (!towersOH)
  {
    std::cout << " RawClusterBuilderTopo::process_event : container TOWER_CALIB_HCALOUT does not exist, aborting " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  _geom_containers[0] = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  _geom_containers[1] = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
  _geom_containers[2] = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");

  if (!_geom_containers[0])
  {
    std::cout << " RawClusterBuilderTopo::process_event : container TOWERGEOM_HCALIN does not exist, aborting " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if (!_geom_containers[1])
  {
    std::cout << " RawClusterBuilderTopo::process_event : container TOWERGEOM_HCALOUT does not exist, aborting " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  if (!_geom_containers[2])
  {
    std::cout << " RawClusterBuilderTopo::process_event : container TOWERGEOM_CEMC does not exist, aborting " << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  if (Verbosity() > 10)
  {
    std::cout << "RawClusterBuilderTopo::process_event: " << towersEM->size() << " TOWER_CALIB_CEMC towers" << std::endl;
    std::cout << "RawClusterBuilderTopo::process_event: " << towersIH->size() << " TOWER_CALIB_HCALIN towers" << std::endl;
    std::cout << "RawClusterBuilderTopo::process_event: " << towersOH->size() << " TOWER_CALIB_HCALOUT towers" << std::endl;

    std::cout << "RawClusterBuilderTopo::process_event: pointer to TOWERGEOM_CEMC: " << _geom_containers[2] << std::endl;
    std::cout << "RawClusterBuilderTopo::process_event: pointer to TOWERGEOM_HCALIN: " << _geom_containers[0] << std::endl;
    std::cout << "RawClusterBuilderTopo::process_event: pointer to TOWERGEOM_HCALOUT: " << _geom_containers[1] << std::endl;
  }

  if (_EMCAL_NETA < 0)
  {
    // define geometry only once if it has not been yet
    _EMCAL_NETA = _geom_containers[2]->get_etabins();
    _EMCAL_NPHI = _geom_containers[2]->get_phibins();

    _EMTOWERMAP_STATUS_ETA_PHI.resize(_EMCAL_NETA, std::vector<int>(_EMCAL_NPHI, -2));
    _EMTOWERMAP_KEY_ETA_PHI.resize(_EMCAL_NETA, std::vector<int>(_EMCAL_NPHI, 0));
    _EMTOWERMAP_E_ETA_PHI.resize(_EMCAL_NETA, std::vector<float>(_EMCAL_NPHI, 0));
  }

  if (_HCAL_NETA < 0)
  {
    // define geometry only once if it has not been yet
    _HCAL_NETA = _geom_containers[1]->get_etabins();
    _HCAL_NPHI = _geom_containers[1]->get_phibins();

    _TOWERMAP_STATUS_LAYER_ETA_PHI.resize(2, std::vector<std::vector<int> >(_HCAL_NETA, std::vector<int>(_HCAL_NPHI, -2)));
    _TOWERMAP_KEY_LAYER_ETA_PHI.resize(2, std::vector<std::vector<int> >(_HCAL_NETA, std::vector<int>(_HCAL_NPHI, 0)));
    _TOWERMAP_E_LAYER_ETA_PHI.resize(2, std::vector<std::vector<float> >(_HCAL_NETA, std::vector<float>(_HCAL_NPHI, 0)));
  }

  // reset maps
  // but note -- do not reset keys!
  for (int ieta = 0; ieta < _EMCAL_NETA; ieta++)
  {
    for (int iphi = 0; iphi < _EMCAL_NPHI; iphi++)
    {
      _EMTOWERMAP_STATUS_ETA_PHI[ieta][iphi] = -2;  // set tower does not exist
      _EMTOWERMAP_E_ETA_PHI[ieta][iphi] = 0;        // set zero energy
    }
  }
  for (int ilayer = 0; ilayer < 2; ilayer++)
  {
    for (int ieta = 0; ieta < _HCAL_NETA; ieta++)
    {
      for (int iphi = 0; iphi < _HCAL_NPHI; iphi++)
      {
        _TOWERMAP_STATUS_LAYER_ETA_PHI[ilayer][ieta][iphi] = -2;  // set tower does not exist
        _TOWERMAP_E_LAYER_ETA_PHI[ilayer][ieta][iphi] = 0;        // set zero energy
      }
    }
  }

  // setup
  std::vector<std::pair<int, float> > list_of_seeds;

  // translate towers to our internal representation
  if (_enable_EMCal)
  {
    RawTowerContainer::ConstRange begin_end_EM = towersEM->getTowers();
    for (RawTowerContainer::ConstIterator rtiter = begin_end_EM.first; rtiter != begin_end_EM.second; ++rtiter)
    {
      RawTower *tower = rtiter->second;
      RawTowerGeom *tower_geom = _geom_containers[2]->get_tower_geometry(tower->get_key());

      int ieta = _geom_containers[2]->get_etabin(tower_geom->get_eta());
      int iphi = _geom_containers[2]->get_phibin(tower_geom->get_phi());
      float this_E = tower->get_energy();

      _EMTOWERMAP_STATUS_ETA_PHI[ieta][iphi] = -1;  // change status to unknown
      _EMTOWERMAP_E_ETA_PHI[ieta][iphi] = this_E;
      _EMTOWERMAP_KEY_ETA_PHI[ieta][iphi] = tower->get_key();

      if (this_E > _sigma_seed * _noise_LAYER[2])
      {
        int ID = get_ID(2, ieta, iphi);
        list_of_seeds.emplace_back(ID, this_E);
        if (Verbosity() > 10)
        {
          std::cout << "RawClusterBuilderTopo::process_event: adding EMCal tower at ieta / iphi = " << ieta << " / " << iphi << " with E = " << this_E << std::endl;
          std::cout << " --> ID = " << ID << " , check ilayer / ieta / iphi = " << get_ilayer_from_ID(ID) << " / " << get_ieta_from_ID(ID) << " / " << get_iphi_from_ID(ID) << std::endl;
        };
      }
    }
  }

  // translate towers to our internal representation
  if (_enable_HCal)
  {
    RawTowerContainer::ConstRange begin_end_IH = towersIH->getTowers();
    for (RawTowerContainer::ConstIterator rtiter = begin_end_IH.first; rtiter != begin_end_IH.second; ++rtiter)
    {
      RawTower *tower = rtiter->second;
      RawTowerGeom *tower_geom = _geom_containers[0]->get_tower_geometry(tower->get_key());

      int ieta = _geom_containers[0]->get_etabin(tower_geom->get_eta());
      int iphi = _geom_containers[0]->get_phibin(tower_geom->get_phi());
      float this_E = tower->get_energy();

      _TOWERMAP_STATUS_LAYER_ETA_PHI[0][ieta][iphi] = -1;  // change status to unknown
      _TOWERMAP_E_LAYER_ETA_PHI[0][ieta][iphi] = this_E;
      _TOWERMAP_KEY_LAYER_ETA_PHI[0][ieta][iphi] = tower->get_key();

      if (this_E > _sigma_seed * _noise_LAYER[0])
      {
        int ID = get_ID(0, ieta, iphi);
        list_of_seeds.emplace_back(ID, this_E);
        if (Verbosity() > 10)
        {
          std::cout << "RawClusterBuilderTopo::process_event: adding IHCal tower at ieta / iphi = " << ieta << " / " << iphi << " with E = " << this_E << std::endl;
          std::cout << " --> ID = " << ID << " , check ilayer / ieta / iphi = " << get_ilayer_from_ID(ID) << " / " << get_ieta_from_ID(ID) << " / " << get_iphi_from_ID(ID) << std::endl;
        };
      }
    }

    RawTowerContainer::ConstRange begin_end_OH = towersOH->getTowers();
    for (RawTowerContainer::ConstIterator rtiter = begin_end_OH.first; rtiter != begin_end_OH.second; ++rtiter)
    {
      RawTower *tower = rtiter->second;
      RawTowerGeom *tower_geom = _geom_containers[1]->get_tower_geometry(tower->get_key());

      int ieta = _geom_containers[1]->get_etabin(tower_geom->get_eta());
      int iphi = _geom_containers[1]->get_phibin(tower_geom->get_phi());
      float this_E = tower->get_energy();

      _TOWERMAP_STATUS_LAYER_ETA_PHI[1][ieta][iphi] = -1;  // change status to unknown
      _TOWERMAP_E_LAYER_ETA_PHI[1][ieta][iphi] = this_E;
      _TOWERMAP_KEY_LAYER_ETA_PHI[1][ieta][iphi] = tower->get_key();

      if (this_E > _sigma_seed * _noise_LAYER[1])
      {
        int ID = get_ID(1, ieta, iphi);
        list_of_seeds.emplace_back(ID, this_E);
        if (Verbosity() > 10)
        {
          std::cout << "RawClusterBuilderTopo::process_event: adding OHCal tower at ieta / iphi = " << ieta << " / " << iphi << " with E = " << this_E << std::endl;
          std::cout << " --> ID = " << ID << " , check ilayer / ieta / iphi = " << get_ilayer_from_ID(ID) << " / " << get_ieta_from_ID(ID) << " / " << get_iphi_from_ID(ID) << std::endl;
        };
      }
    }
  }

  if (Verbosity() > 10)
  {
    for (unsigned int n = 0; n < list_of_seeds.size(); n++)
    {
      std::cout << "RawClusterBuilderTopo::process_event: unsorted seed element n = " << n << " , ID / E = " << list_of_seeds.at(n).first << " / " << list_of_seeds.at(n).second << std::endl;
    }
  }

  std::sort(list_of_seeds.begin(), list_of_seeds.end(), sort_by_pair_second);

  if (Verbosity() > 10)
  {
    for (unsigned int n = 0; n < list_of_seeds.size(); n++)
    {
      std::cout << "RawClusterBuilderTopo::process_event: sorted seed element n = " << n << " , ID / E = " << list_of_seeds.at(n).first << " / " << list_of_seeds.at(n).second << std::endl;
    }
  }

  if (Verbosity() > 0)
    std::cout << "RawClusterBuilderTopo::process_event: initialized with " << list_of_seeds.size() << " seeds with E > 4*sigma " << std::endl;

  int cluster_index = 0;  // begin counting clusters

  std::vector<std::vector<int> > all_cluster_towers;  // store final cluster tower lists here

  while (list_of_seeds.size() > 0)
  {
    int seed_ID = list_of_seeds.at(0).first;
    list_of_seeds.erase(list_of_seeds.begin());

    if (Verbosity() > 5)
    {
      std::cout << " RawClusterBuilderTopo::process_event: in seeded loop, current seed has ID = " << seed_ID << " , length of remaining seed vector = " << list_of_seeds.size() << std::endl;
    }

    // if this seed was already claimed by some other seed during its growth, remove it and do nothing
    int seed_status = get_status_from_ID(seed_ID);
    if (seed_status > -1)
    {
      if (Verbosity() > 10)
      {
        std::cout << " --> already owned by cluster # " << seed_status << std::endl;
      }
      continue;  // go onto the next iteration of the loop
    }

    // this seed tower now owned by new cluster
    set_status_by_ID(seed_ID, cluster_index);

    std::vector<int> cluster_tower_ID;
    cluster_tower_ID.push_back(seed_ID);

    std::vector<int> grow_tower_ID;
    grow_tower_ID.push_back(seed_ID);

    // iteratively process growth towers, adding > 2 * sigma neighbors to the list for further checking

    if (Verbosity() > 5)
    {
      std::cout << " RawClusterBuilderTopo::process_event: Entering Growth stage for cluster " << cluster_index << std::endl;
    }

    while (grow_tower_ID.size() > 0)
    {
      int grow_ID = grow_tower_ID.at(0);
      grow_tower_ID.erase(grow_tower_ID.begin());

      if (Verbosity() > 5)
      {
        std::cout << " --> cluster " << cluster_index << ", growth stage, examining neighbors of ID " << grow_ID << ", " << grow_tower_ID.size() << " grow towers left" << std::endl;
      }

      std::vector<int> adjacent_tower_IDs = get_adjacent_towers_by_ID(grow_ID);

      for (int this_adjacent_tower_ID : adjacent_tower_IDs)
      {
        if (Verbosity() > 10)
        {
          std::cout << " --> --> --> checking possible adjacent tower with ID " << this_adjacent_tower_ID << " : ";
        }
        int test_layer = get_ilayer_from_ID(this_adjacent_tower_ID);

        // if tower does not exist, continue
        if (get_status_from_ID(this_adjacent_tower_ID) == -2)
        {
          if (Verbosity() > 10)
          {
            std::cout << "does not exist " << std::endl;
          }
          continue;
        }

        // if tower is owned by THIS cluster already, continue
        if (get_status_from_ID(this_adjacent_tower_ID) == cluster_index)
        {
          if (Verbosity() > 10) std::cout << "already owned by this cluster index " << cluster_index << std::endl;
          continue;
        }

        // if tower has < 2*sigma energy, continue
        if (get_E_from_ID(this_adjacent_tower_ID) < _sigma_grow * _noise_LAYER[test_layer])
        {
          if (Verbosity() > 10) std::cout << "E = " << get_E_from_ID(this_adjacent_tower_ID) << " under 2*sigma threshold " << std::endl;
          continue;
        }

        // if tower is owned by somebody else, continue (although should this really happen?)
        if (get_status_from_ID(this_adjacent_tower_ID) > -1)
        {
          if (Verbosity() > 10) std::cout << "ERROR! in growth stage, encountered >2sigma tower which is already owned?!" << std::endl;
          continue;
        }

        // tower good to be added to cluster and to list of grow towers
        grow_tower_ID.push_back(this_adjacent_tower_ID);
        cluster_tower_ID.push_back(this_adjacent_tower_ID);
        set_status_by_ID(this_adjacent_tower_ID, cluster_index);
        if (Verbosity() > 10) std::cout << "add this tower ( ID " << this_adjacent_tower_ID << " ) to grow list " << std::endl;
      }

      if (Verbosity() > 5) std::cout << " --> after examining neighbors, grow list is now " << grow_tower_ID.size() << ", # of towers in cluster = " << cluster_tower_ID.size() << std::endl;
    }

    // done growing cluster, now add on perimeter towers with E > 0 * sigma
    if (Verbosity() > 5)
    {
      std::cout << " RawClusterBuilderTopo::process_event: Entering Perimeter stage for cluster " << cluster_index << std::endl;
    }
    // we'll be adding on to the cluster list, so get the # of core towers first
    int n_core_towers = cluster_tower_ID.size();

    for (int ic = 0; ic < n_core_towers; ic++)
    {
      int core_ID = cluster_tower_ID.at(ic);

      if (Verbosity() > 5)
      {
        std::cout << " --> cluster " << cluster_index << ", perimeter stage, examining neighbors of ID " << core_ID << ", core cluster # " << ic << " of " << n_core_towers << " total " << std::endl;
      }
      std::vector<int> adjacent_tower_IDs = get_adjacent_towers_by_ID(core_ID);

      for (int this_adjacent_tower_ID : adjacent_tower_IDs)
      {
        if (Verbosity() > 10) std::cout << " --> --> --> checking possible adjacent tower with ID " << this_adjacent_tower_ID << " : ";

        int test_layer = get_ilayer_from_ID(this_adjacent_tower_ID);

        // if tower does not exist, continue
        if (get_status_from_ID(this_adjacent_tower_ID) == -2)
        {
          if (Verbosity() > 10) std::cout << "does not exist " << std::endl;
          continue;
        }

        // if tower is owned by somebody else (including current cluster), continue. ( allowed during perimeter fixing state )
        if (get_status_from_ID(this_adjacent_tower_ID) > -1)
        {
          if (Verbosity() > 10) std::cout << "already owned by other cluster index " << get_status_from_ID(this_adjacent_tower_ID) << std::endl;
          continue;
        }

        // if tower has < 0*sigma energy, continue
        if (get_E_from_ID(this_adjacent_tower_ID) < _sigma_peri * _noise_LAYER[test_layer])
        {
          if (Verbosity() > 10) std::cout << "E = " << get_E_from_ID(this_adjacent_tower_ID) << " under 0*sigma threshold " << std::endl;
          continue;
        }

        // perimeter tower good to be added to cluster
        cluster_tower_ID.push_back(this_adjacent_tower_ID);
        set_status_by_ID(this_adjacent_tower_ID, cluster_index);
        if (Verbosity() > 10) std::cout << "add this tower ( ID " << this_adjacent_tower_ID << " ) to cluster " << std::endl;
      }

      if (Verbosity() > 5) std::cout << " --> after examining perimeter neighbors, # of towers in cluster is now = " << cluster_tower_ID.size() << std::endl;
    }

    // keep track of these
    all_cluster_towers.push_back(cluster_tower_ID);

    // increment cluster index for next one
    cluster_index++;
  }

  if (Verbosity() > 0) std::cout << "RawClusterBuilderTopo::process_event: " << cluster_index << " topo-clusters initially reconstructed, entering splitting step" << std::endl;

  // now entering cluster splitting stage
  int original_cluster_index = cluster_index;  // since it may be updated
  for (int cl = 0; cl < original_cluster_index; cl++)
  {
    std::vector<int> original_towers = all_cluster_towers.at(cl);

    if (!_do_split)
    {
      // don't run splitting, just export entire cluster as it is
      if (Verbosity() > 2) std::cout << "RawClusterBuilderTopo::process_event: splitting step disabled, cluster " << cluster_index << " is final" << std::endl;
      export_single_cluster(original_towers);
      continue;
    }

    std::vector<std::pair<int, float> > local_maxima_ID;

    // iterate through each tower, looking for maxima
    for (int tower_ID : original_towers)
    {
      if (Verbosity() > 10) std::cout << " -> examining tower ID " << tower_ID << " for possible local maximum " << std::endl;

      // check minimum energy
      if (get_E_from_ID(tower_ID) < _local_max_minE_LAYER[get_ilayer_from_ID(tower_ID)])
      {
        if (Verbosity() > 10) std::cout << " -> -> energy E = " << get_E_from_ID(tower_ID) << " < " << _local_max_minE_LAYER[get_ilayer_from_ID(tower_ID)] << " too low" << std::endl;
        continue;
      }

      // examine neighbors
      std::vector<int> adjacent_tower_IDs = get_adjacent_towers_by_ID(tower_ID);
      int neighbors_in_cluster = 0;

      // check for higher neighbox
      bool has_higher_neighbor = false;
      for (int this_adjacent_tower_ID : adjacent_tower_IDs)
      {
        if (get_status_from_ID(this_adjacent_tower_ID) != cl) continue;  // only consider neighbors in cluster, obviously

        neighbors_in_cluster++;

        if (get_E_from_ID(this_adjacent_tower_ID) > get_E_from_ID(tower_ID))
        {
          if (Verbosity() > 10) std::cout << " -> -> has higher-energy neighbor ID / E = " << this_adjacent_tower_ID << " / " << get_E_from_ID(this_adjacent_tower_ID) << std::endl;
          has_higher_neighbor = true;  // at this point we can break -- we won't need to count the number of good neighbors, since we won't even pass the E_neighbor test
          break;
        }
      }

      if (has_higher_neighbor) continue;  // if we broke out, now continue

      // check number of neighbors
      if (neighbors_in_cluster < 4)
      {
        if (Verbosity() > 10) std::cout << " -> -> too few neighbors N = " << neighbors_in_cluster << std::endl;
        continue;
      }

      local_maxima_ID.emplace_back(tower_ID, get_E_from_ID(tower_ID));
    }

    // check for possible EMCal-OHCal seed overlaps

    for (unsigned int n = 0; n < local_maxima_ID.size(); n++)
    {
      // only look at I/OHCal local maxima
      std::pair<int, float> this_LM = local_maxima_ID.at(n);
      if (get_ilayer_from_ID(this_LM.first) == 2) continue;

      float this_phi = _geom_containers[get_ilayer_from_ID(this_LM.first)]->get_phicenter(get_iphi_from_ID(this_LM.first));
      if (this_phi > M_PI) this_phi -= 2 * M_PI;
      float this_eta = _geom_containers[get_ilayer_from_ID(this_LM.first)]->get_etacenter(get_ieta_from_ID(this_LM.first));

      bool has_EM_overlap = false;

      // check all other local maxima for overlaps
      for (unsigned int n2 = 0; n2 < local_maxima_ID.size(); n2++)
      {
        if (n == n2) continue;  // don't check the same one

        // only look at EMCal local mazima
        std::pair<int, float> this_LM2 = local_maxima_ID.at(n2);
        if (get_ilayer_from_ID(this_LM2.first) != 2) continue;

        float this_phi2 = _geom_containers[get_ilayer_from_ID(this_LM2.first)]->get_phicenter(get_iphi_from_ID(this_LM2.first));
        if (this_phi2 > M_PI) this_phi -= 2 * M_PI;
        float this_eta2 = _geom_containers[get_ilayer_from_ID(this_LM2.first)]->get_etacenter(get_ieta_from_ID(this_LM2.first));

        // calculate geometric dR
        float dR = calculate_dR(this_eta, this_eta2, this_phi, this_phi2);

        // check for and report overlaps
        if (dR < 0.15)
        {
          has_EM_overlap = true;
          if (Verbosity() > 2)
          {
            std::cout << "RawClusterBuilderTopo::process_event : removing I/OHal local maximum (ID,E,phi,eta = " << this_LM.first << ", " << this_LM.second << ", " << this_phi << ", " << this_eta << "), ";
            std::cout << "due to EM overlap (ID,E,phi,eta = " << this_LM2.first << ", " << this_LM2.second << ", " << this_phi2 << ", " << this_eta2 << "), dR = " << dR << std::endl;
          }
          break;
        }
      }

      if (has_EM_overlap)
      {
        // remove the I/OHCal local maximum from the list
        local_maxima_ID.erase(local_maxima_ID.begin() + n);
        // make sure to back up one index...
        n = n - 1;
      }  // otherwise, keep this local maximum
    }

    // only now print out full set of local maxima
    if (Verbosity() > 2)
    {
      for (auto this_LM : local_maxima_ID)
      {
        int tower_ID = this_LM.first;
        std::cout << "RawClusterBuilderTopo::process_event in cluster " << cl << ", tower ID " << tower_ID << " is LOCAL MAXIMUM with layer / E = " << get_ilayer_from_ID(tower_ID) << " / " << get_E_from_ID(tower_ID) << ", ";
        float this_phi = _geom_containers[get_ilayer_from_ID(tower_ID)]->get_phicenter(get_iphi_from_ID(tower_ID));
        if (this_phi > M_PI) this_phi -= 2 * M_PI;
        std::cout << " eta / phi = " << _geom_containers[get_ilayer_from_ID(tower_ID)]->get_etacenter(get_ieta_from_ID(tower_ID)) << " / " << this_phi << std::endl;
      }
    }

    // do we have only 1 or 0 local maxima?
    if (local_maxima_ID.size() <= 1)
    {
      if (Verbosity() > 2) std::cout << "RawClusterBuilderTopo::process_event cluster " << cl << " has only " << local_maxima_ID.size() << " local maxima, not splitting " << std::endl;
      export_single_cluster(original_towers);

      continue;
    }

    // engage splitting procedure!

    if (Verbosity() > 2)
    {
      std::cout << "RawClusterBuilderTopo::process_event splitting cluster " << cl << " into " << local_maxima_ID.size() << " according to local maxima!" << std::endl;
    }
    // translate all cluster towers to a map which keeps track of their ownership
    // -1 means unseen
    // -2 means seen and in the seed list now (e.g. don't add it to the seed list again)
    // -3 shared tower, ignore going forward...
    std::map<int, std::pair<int, int> > tower_ownership;
    for (int &original_tower : original_towers)
    {
      tower_ownership[original_tower] = std::pair<int, int>(-1, -1);  // initialize all towers as un-seen
    }
    std::vector<int> seed_list;
    std::vector<int> neighbor_list;
    std::vector<int> shared_list;

    // sort maxima before populating seed list
    std::sort(local_maxima_ID.begin(), local_maxima_ID.end(), sort_by_pair_second);

    // initialize neighbor list
    for (unsigned int s = 0; s < local_maxima_ID.size(); s++)
    {
      tower_ownership[local_maxima_ID.at(s).first] = std::pair<int, int>(s, -1);
      neighbor_list.push_back(local_maxima_ID.at(s).first);
    }

    if (Verbosity() > 100)
    {
      for (int &original_tower : original_towers)
      {
        std::pair<int, int> the_pair = tower_ownership[original_tower];
        std::cout << " Debug Pre-Split: tower_ownership[ " << original_tower << " ] = ( " << the_pair.first << ", " << the_pair.second << " ) ";
        std::cout << " , layer / ieta / iphi = " << get_ilayer_from_ID(original_tower) << " / " << get_ieta_from_ID(original_tower) << " / " << get_iphi_from_ID(original_tower);
        std::cout << std::endl;
      }
    }

    bool first_pass = true;

    do
    {
      if (Verbosity() > 5)
      {
        std::cout << " -> starting split loop with " << seed_list.size() << " seed, " << neighbor_list.size() << " neighbor, and " << shared_list.size() << " shared towers " << std::endl;
      }
      // go through neighbor list, assigning ownership only via the seed list
      std::vector<int> new_ownerships;

      for (unsigned int n = 0; n < neighbor_list.size(); n++)
      {
        int neighbor_ID = neighbor_list.at(n);

        if (Verbosity() > 10)
        {
          std::cout << " -> -> looking at neighbor " << n << " (tower ID " << neighbor_ID << " ) of " << neighbor_list.size() << " total" << std::endl;
        }
        if (first_pass)
        {
          if (Verbosity() > 10)
          {
            std::cout << " -> -> -> special first pass rules, this tower already owned by pseudocluster " << tower_ownership[neighbor_ID].first << std::endl;
          }
          new_ownerships.push_back(tower_ownership[neighbor_ID].first);
        }
        else
        {
          std::map<int, bool> pseudocluster_adjacency;
          for (unsigned int s = 0; s < local_maxima_ID.size(); s++)
          {
            pseudocluster_adjacency[s] = false;
          }
          // look over all towers THIS one is adjacent to, and count up...
          std::vector<int> adjacent_tower_IDs = get_adjacent_towers_by_ID(neighbor_ID);
          for (int this_adjacent_tower_ID : adjacent_tower_IDs)
          {
            if (get_status_from_ID(this_adjacent_tower_ID) != cl) continue;
            if (tower_ownership[this_adjacent_tower_ID].first > -1)
            {
              if (Verbosity() > 20)
              {
                std::cout << " -> -> -> adjacent tower to this one, with ID " << this_adjacent_tower_ID << " , is owned by pseudocluster " << tower_ownership[this_adjacent_tower_ID].first << std::endl;
              }
              pseudocluster_adjacency[tower_ownership[this_adjacent_tower_ID].first] = true;
            }
          }
          int n_pseudocluster_adjacent = 0;
          int last_adjacent_pseudocluster = -1;
          for (unsigned int s = 0; s < local_maxima_ID.size(); s++)
          {
            if (pseudocluster_adjacency[s])
            {
              last_adjacent_pseudocluster = s;
              n_pseudocluster_adjacent++;
              if (Verbosity() > 20)
              {
                std::cout << " -> -> adjacent to pseudocluster " << s << std::endl;
              }
            }
          }

          if (n_pseudocluster_adjacent == 0)
          {
            std::cout << " -> -> ERROR! How can a neighbor tower at this stage be adjacent to no pseudoclusters?? " << std::endl;
            new_ownerships.push_back(9999);
          }
          else if (n_pseudocluster_adjacent == 1)
          {
            if (Verbosity() > 10)
            {
              std::cout << " -> -> neighbor tower " << neighbor_ID << " is ONLY adjacent to one pseudocluster # " << last_adjacent_pseudocluster << std::endl;
            }
            new_ownerships.push_back(last_adjacent_pseudocluster);
          }
          else
          {
            if (Verbosity() > 10)
            {
              std::cout << " -> -> neighbor tower " << neighbor_ID << " is adjacent to " << n_pseudocluster_adjacent << " pseudoclusters, move to shared list " << std::endl;
            }
            new_ownerships.push_back(-3);
          }
        }
      }

      if (Verbosity() > 5)
      {
        std::cout << " -> now updating status of all " << neighbor_list.size() << " original neighbors " << std::endl;
      }
      // transfer neighbor list to seed list or shared list
      for (unsigned int n = 0; n < neighbor_list.size(); n++)
      {
        int neighbor_ID = neighbor_list.at(n);
        if (new_ownerships.at(n) > -1)
        {
          tower_ownership[neighbor_ID] = std::pair<int, int>(new_ownerships.at(n), -1);
          seed_list.push_back(neighbor_ID);
          if (Verbosity() > 20)
          {
            std::cout << " -> -> neighbor ID " << neighbor_ID << " has new status " << new_ownerships.at(n) << std::endl;
          }
        }
        if (new_ownerships.at(n) == -3)
        {
          tower_ownership[neighbor_ID] = std::pair<int, int>(-3, -1);
          shared_list.push_back(neighbor_ID);
          if (Verbosity() > 20)
          {
            std::cout << " -> -> neighbor ID " << neighbor_ID << " has new status " << -3 << std::endl;
          }
        }
      }

      if (Verbosity() > 5)
      {
        std::cout << " producing a new neighbor list ... " << std::endl;
      }
      // populate a new neighbor list from the about-to-be-owned towers before transferring this one
      std::list<int> new_neighbor_list;
      for (unsigned int n = 0; n < neighbor_list.size(); n++)
      {
        int neighbor_ID = neighbor_list.at(n);
        if (new_ownerships.at(n) > -1)
        {
          std::vector<int> adjacent_tower_IDs = get_adjacent_towers_by_ID(neighbor_ID);

          for (int this_adjacent_tower_ID : adjacent_tower_IDs)
          {
            if (get_status_from_ID(this_adjacent_tower_ID) != cl) continue;
            if (tower_ownership[this_adjacent_tower_ID].first == -1)
            {
              new_neighbor_list.push_back(this_adjacent_tower_ID);
              if (Verbosity() > 5)
              {
                std::cout << " -> queueing up to add tower " << this_adjacent_tower_ID << " , neighbor of tower " << neighbor_ID << " to new neighbor list" << std::endl;
              }
            }
          }
        }
      }

      if (Verbosity() > 5)
      {
        std::cout << " new neighbor list has size " << new_neighbor_list.size() << ", but after removing duplicate elements: ";
        new_neighbor_list.sort();
        new_neighbor_list.unique();
        std::cout << new_neighbor_list.size() << std::endl;
      }

      // clear neighbor list!
      neighbor_list.clear();

      // now transfer over new neighbor list
      for (; new_neighbor_list.size() > 0;)
      {
        neighbor_list.push_back(new_neighbor_list.front());
        new_neighbor_list.pop_front();
      }

      first_pass = false;

    } while (neighbor_list.size() > 0);

    if (Verbosity() > 100)
    {
      for (int &original_tower : original_towers)
      {
        std::pair<int, int> the_pair = tower_ownership[original_tower];
        std::cout << " Debug Mid-Split: tower_ownership[ " << original_tower << " ] = ( " << the_pair.first << ", " << the_pair.second << " ) ";
        std::cout << " , layer / ieta / iphi = " << get_ilayer_from_ID(original_tower) << " / " << get_ieta_from_ID(original_tower) << " / " << get_iphi_from_ID(original_tower);
        std::cout << std::endl;
        if (the_pair.first == -1)
        {
          std::vector<int> adjacent_tower_IDs = get_adjacent_towers_by_ID(original_tower);

          for (int this_adjacent_tower_ID : adjacent_tower_IDs)
          {
            if (get_status_from_ID(this_adjacent_tower_ID) != cl) continue;
            std::cout << "    -> adjacent to add tower " << this_adjacent_tower_ID << " , which has status " << tower_ownership[this_adjacent_tower_ID].first << std::endl;
          }
        }
      }
    }

    // calculate pseudocluster energies and positions
    std::vector<float> pseudocluster_sumeta;
    std::vector<float> pseudocluster_sumphi;
    std::vector<float> pseudocluster_sumE;
    std::vector<int> pseudocluster_ntower;
    std::vector<float> pseudocluster_eta;
    std::vector<float> pseudocluster_phi;

    pseudocluster_sumeta.resize(local_maxima_ID.size(), 0);
    pseudocluster_sumphi.resize(local_maxima_ID.size(), 0);
    pseudocluster_sumE.resize(local_maxima_ID.size(), 0);
    pseudocluster_ntower.resize(local_maxima_ID.size(), 0);

    for (int &original_tower : original_towers)
    {
      std::pair<int, int> the_pair = tower_ownership[original_tower];
      if (the_pair.first > -1)
      {
        float this_ID = original_tower;
        pseudocluster_sumE[the_pair.first] += get_E_from_ID(this_ID);
        float this_eta = _geom_containers[get_ilayer_from_ID(this_ID)]->get_etacenter(get_ieta_from_ID(this_ID));
        float this_phi = _geom_containers[get_ilayer_from_ID(this_ID)]->get_phicenter(get_iphi_from_ID(this_ID));
        //float this_phi = ( get_ilayer_from_ID( this_ID ) == 2 ? geomEM->get_phicenter( get_iphi_from_ID( this_ID ) ) : geomOH->get_phicenter( get_iphi_from_ID( this_ID ) ) );
        pseudocluster_sumeta[the_pair.first] += this_eta;
        pseudocluster_sumphi[the_pair.first] += this_phi;
        pseudocluster_ntower[the_pair.first] += 1;
      }
    }

    for (unsigned int pc = 0; pc < local_maxima_ID.size(); pc++)
    {
      pseudocluster_eta.push_back(pseudocluster_sumeta.at(pc) / pseudocluster_ntower.at(pc));
      pseudocluster_phi.push_back(pseudocluster_sumphi.at(pc) / pseudocluster_ntower.at(pc));

      if (Verbosity() > 2)
      {
        std::cout << "RawClusterBuilderTopo::process_event pseudocluster #" << pc << ", E / eta / phi / Ntower = " << pseudocluster_sumE.at(pc) << " / " << pseudocluster_eta.at(pc) << " / " << pseudocluster_phi.at(pc) << " / " << pseudocluster_ntower.at(pc) << std::endl;
      }
    }

    if (Verbosity() > 2)
    {
      std::cout << "RawClusterBuilderTopo::process_event now splitting up shared clusters (including unassigned clusters), initial shared list has size " << shared_list.size() << std::endl;
    }
    // iterate through shared cells, identifying which two they belong to
    while (shared_list.size() > 0)
    {
      // pick the first cell and pop off list
      int shared_ID = shared_list.at(0);
      shared_list.erase(shared_list.begin());

      if (Verbosity() > 5)
      {
        std::cout << " -> looking at shared tower " << shared_ID << ", after this one there are " << shared_list.size() << " shared towers left " << std::endl;
      }
      // look through adjacent pseudoclusters, taking two with highest energies
      std::vector<bool> pseudocluster_adjacency;
      pseudocluster_adjacency.resize(local_maxima_ID.size(), false);

      std::vector<int> adjacent_tower_IDs = get_adjacent_towers_by_ID(shared_ID);

      for (int this_adjacent_tower_ID : adjacent_tower_IDs)
      {
        if (get_status_from_ID(this_adjacent_tower_ID) != cl) continue;
        if (tower_ownership[this_adjacent_tower_ID].first > -1)
        {
          pseudocluster_adjacency[tower_ownership[this_adjacent_tower_ID].first] = true;
        }
        if (tower_ownership[this_adjacent_tower_ID].second > -1)
        {  // can inherit adjacency from shared cluster
          pseudocluster_adjacency[tower_ownership[this_adjacent_tower_ID].second] = true;
        }
        // at the same time, add unowned towers to the list for later examination
        if (tower_ownership[this_adjacent_tower_ID].first == -1)
        {
          shared_list.push_back(this_adjacent_tower_ID);
          tower_ownership[this_adjacent_tower_ID] = std::pair<int, int>(-3, -1);
          if (Verbosity() > 10)
          {
            std::cout << " -> while looking at neighbors, have added un-examined tower " << this_adjacent_tower_ID << " to shared list " << std::endl;
          }
        }
      }

      // now figure out which pseudoclustes this shared tower is adjacent to...
      int highest_pseudocluster_index = -1;
      int second_highest_pseudocluster_index = -1;

      float highest_pseudocluster_E = -1;
      float second_highest_pseudocluster_E = -2;

      for (unsigned int n = 0; n < pseudocluster_adjacency.size(); n++)
      {
        if (!pseudocluster_adjacency[n]) continue;

        if (pseudocluster_sumE[n] > highest_pseudocluster_E)
        {
          second_highest_pseudocluster_E = highest_pseudocluster_E;
          second_highest_pseudocluster_index = highest_pseudocluster_index;

          highest_pseudocluster_E = pseudocluster_sumE[n];
          highest_pseudocluster_index = n;
        }
        else if (pseudocluster_sumE[n] > second_highest_pseudocluster_E)
        {
          second_highest_pseudocluster_E = pseudocluster_sumE[n];
          second_highest_pseudocluster_index = n;
        }
      }

      if (Verbosity() > 5)
      {
        std::cout << " -> highest pseudoclusters its adjacent to are " << highest_pseudocluster_index << " ( E = " << highest_pseudocluster_E << " ) and " << second_highest_pseudocluster_index << " ( E = " << second_highest_pseudocluster_E << " ) " << std::endl;
      }
      // assign these clusters as owners
      tower_ownership[shared_ID] = std::pair<int, int>(highest_pseudocluster_index, second_highest_pseudocluster_index);
    }

    if (Verbosity() > 100)
    {
      for (int &original_tower : original_towers)
      {
        std::pair<int, int> the_pair = tower_ownership[original_tower];
        std::cout << " Debug Post-Split: tower_ownership[ " << original_tower << " ] = ( " << the_pair.first << ", " << the_pair.second << " ) ";
        std::cout << " , layer / ieta / iphi = " << get_ilayer_from_ID(original_tower) << " / " << get_ieta_from_ID(original_tower) << " / " << get_iphi_from_ID(original_tower);
        std::cout << std::endl;
        if (the_pair.first == -1)
        {
          std::vector<int> adjacent_tower_IDs = get_adjacent_towers_by_ID(original_tower);

          for (int this_adjacent_tower_ID : adjacent_tower_IDs)
          {
            if (get_status_from_ID(this_adjacent_tower_ID) != cl) continue;
            std::cout << " -> adjacent to add tower " << this_adjacent_tower_ID << " , which has status " << tower_ownership[this_adjacent_tower_ID].first << std::endl;
          }
        }
      }
    }

    // call helper function
    export_clusters(original_towers, tower_ownership, local_maxima_ID.size(), pseudocluster_sumE, pseudocluster_eta, pseudocluster_phi);
  }

  if (Verbosity() > 1)
  {
    std::cout << "RawClusterBuilderTopo::process_event after splitting (if any) final clusters output to node are: " << std::endl;
    RawClusterContainer::ConstRange begin_end = _clusters->getClusters();
    int ncl = 0;
    for (RawClusterContainer::ConstIterator hiter = begin_end.first; hiter != begin_end.second; ++hiter)
    {
      std::cout << "-> #" << ncl++ << " ";
      hiter->second->identify();
      std::cout << std::endl;
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int RawClusterBuilderTopo::End(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void RawClusterBuilderTopo::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  PHCompositeNode *dstNode = static_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cerr << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    throw std::runtime_error("Failed to find DST node in RawClusterBuilderTopo::CreateNodes");
  }

  PHNodeIterator dstiter(dstNode);
  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TOPOCLUSTER"));
  if (!DetNode)
  {
    DetNode = new PHCompositeNode("TOPOCLUSTER");
    dstNode->addNode(DetNode);
  }

  _clusters = new RawClusterContainer();

  PHIODataNode<PHObject> *clusterNode = new PHIODataNode<PHObject>(_clusters, ClusterNodeName, "PHObject");
  DetNode->addNode(clusterNode);
}
