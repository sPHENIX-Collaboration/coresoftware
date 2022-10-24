#include "RawTowerCombiner.h"

#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerDefs.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeomv1.h>
#include <calobase/RawTowerv1.h>

#include <fun4all/Fun4AllBase.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/SubsysReco.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>

#include <CLHEP/Vector/ThreeVector.h>

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

RawTowerCombiner::RawTowerCombiner(const std::string &name)
  : SubsysReco(name)
  , _tower_node_prefix("SIM")
  , _n_combine_eta(2)
  , _n_combine_phi(2)
  , _towers(nullptr)
  , detector("NONE")
{
}

int RawTowerCombiner::InitRun(PHCompositeNode *topNode)
{
  if (_n_combine_eta == 0)
  {
    std::cout << __PRETTY_FUNCTION__ << " Fatal error _n_combine_eta==0" << std::endl;

    exit(1243);
  }

  if (_n_combine_phi == 0)
  {
    std::cout << __PRETTY_FUNCTION__ << " Fatal error _n_combine_phi==0" << std::endl;

    exit(1243);
  }

  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode",
                                                           "DST"));
  if (!dstNode)
  {
    std::cout << __PRETTY_FUNCTION__ << "DST Node missing, doing nothing."
              << std::endl;
    exit(1);
  }

  try
  {
    CreateNodes(topNode);
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    //exit(1);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int RawTowerCombiner::process_event(PHCompositeNode * /*topNode*/)
{
  assert(_towers);

  const double input_e_sum = _towers->getTotalEdep();
  const double input_n_tower = _towers->size();

  if (Verbosity())
  {
    std::cout << __PRETTY_FUNCTION__ << "Process event entered" << std::endl;
  }

  using tower_id_t = std::pair<int, int>;
  using new_tower_map_t = std::map<tower_id_t, RawTower *>;
  new_tower_map_t new_tower_map;

  RawTowerContainer::ConstRange all_towers = _towers->getTowers();
  for (RawTowerContainer::ConstIterator it = all_towers.first;
       it != all_towers.second; ++it)
  {
    const RawTower *input_tower = it->second;
    assert(input_tower);

    const int intput_eta = input_tower->get_bineta();
    const int intput_phi = input_tower->get_binphi();
    const int output_eta = get_output_bin_eta(intput_eta);
    const int output_phi = get_output_bin_phi(intput_phi);

    RawTower *output_tower = nullptr;
    new_tower_map_t::iterator it_new = new_tower_map.find(
        std::make_pair(output_eta, output_phi));

    if (it_new == new_tower_map.end())
    {
      output_tower = new RawTowerv1(*input_tower);
      assert(output_tower);
      new_tower_map[std::make_pair(output_eta, output_phi)] = output_tower;

      if (Verbosity() >= VERBOSITY_MORE)
      {
        std::cout << __PRETTY_FUNCTION__ << "::" << detector << "::"
                  << " new output tower (prior to tower ID assignment): ";
        output_tower->identify();
      }
    }
    else
    {
      output_tower = it_new->second;

      output_tower->set_energy(
          output_tower->get_energy() + input_tower->get_energy());

      output_tower->set_time(
          (std::abs(output_tower->get_energy()) * output_tower->get_time() + std::abs(input_tower->get_energy()) * input_tower->get_time())  //
          / (std::abs(output_tower->get_energy()) + std::abs(input_tower->get_energy()) + 1e-9)                                              //avoid devide 0
      );

      RawTower::CellConstRange cell_range = input_tower->get_g4cells();

      for (RawTower::CellConstIterator cell_iter = cell_range.first;
           cell_iter != cell_range.second; ++cell_iter)
      {
        output_tower->add_ecell(cell_iter->first, cell_iter->second);
      }

      RawTower::ShowerConstRange shower_range =
          input_tower->get_g4showers();

      for (RawTower::ShowerConstIterator shower_iter = shower_range.first;
           shower_iter != shower_range.second; ++shower_iter)
      {
        output_tower->add_eshower(shower_iter->first,
                                  shower_iter->second);
      }

      if (Verbosity() >= VERBOSITY_MORE)
      {
        std::cout << __PRETTY_FUNCTION__ << "::" << detector << "::"
                  << " merget into output tower (prior to tower ID assignment) : ";
        output_tower->identify();
      }
    }
  }

  // replace content in tower container
  _towers->Reset();

  for (auto & it : new_tower_map)
  {
    const int eta = it.first.first;
    const int phi = it.first.second;
    RawTower *tower = it.second;

    _towers->AddTower(eta, phi, tower);
  }

  if (Verbosity())
  {
    std::cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
              << "input sum energy = " << input_e_sum << " from " << input_n_tower << " towers, merged sum energy = "
              << _towers->getTotalEdep() << " from " << _towers->size() << " towers" << std::endl;
    assert(input_e_sum == _towers->getTotalEdep());
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int RawTowerCombiner::End(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void RawTowerCombiner::CreateNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst(
      "PHCompositeNode", "RUN"));
  if (!runNode)
  {
    std::cerr << __PRETTY_FUNCTION__ << "Run Node missing, doing nothing."
              << std::endl;
    throw std::runtime_error(
        "Failed to find Run node in RawTowerCombiner::CreateNodes");
  }

  const RawTowerDefs::CalorimeterId caloid =
      RawTowerDefs::convert_name_to_caloid(detector);

  const std::string iTowerGeomNodeName = "TOWERGEOM_" + detector;
  RawTowerGeomContainer *towergeom = findNode::getClass<
      RawTowerGeomContainer>(topNode, iTowerGeomNodeName);
  if (!towergeom)
  {
    std::cerr << __PRETTY_FUNCTION__ << " - " << iTowerGeomNodeName
              << " Node missing, doing nothing." << std::endl;
    throw std::runtime_error(
        "Failed to find input tower geometry node in RawTowerCombiner::CreateNodes");
  }
  assert(towergeom);

  const double r = towergeom->get_radius();

  const int new_phibins = ceil(
      double(towergeom->get_phibins()) / double(_n_combine_phi));
  const int new_etabins = ceil(
      double(towergeom->get_etabins()) / double(_n_combine_eta));

  using bound_t = std::pair<double, double>;
  using bound_map_t = std::vector<bound_t>;

  bound_map_t eta_bound_map;
  bound_map_t phi_bound_map;

  for (int ibin = 0; ibin < new_phibins; ibin++)
  {
    const int first_bin = ibin * _n_combine_phi;
    assert(first_bin >= 0 && first_bin < towergeom->get_phibins());

    int last_bin = (ibin + 1) * _n_combine_phi - 1;
    if (last_bin >= towergeom->get_phibins())
      last_bin = towergeom->get_phibins();

    const std::pair<double, double> range1 = towergeom->get_phibounds(
        first_bin);
    const std::pair<double, double> range2 = towergeom->get_phibounds(
        last_bin);

    phi_bound_map.push_back(std::make_pair(range1.first, range2.second));
  }

  for (int ibin = 0; ibin < new_etabins; ibin++)
  {
    const int first_bin = ibin * _n_combine_eta;
    assert(first_bin >= 0 && first_bin < towergeom->get_etabins());

    int last_bin = (ibin + 1) * _n_combine_eta - 1;
    if (last_bin >= towergeom->get_etabins())
      last_bin = towergeom->get_etabins();

    const std::pair<double, double> range1 = towergeom->get_etabounds(
        first_bin);
    const std::pair<double, double> range2 = towergeom->get_etabounds(
        last_bin);

    eta_bound_map.push_back(std::make_pair(range1.first, range2.second));
  }

  // now update the tower geometry object with the new tower structure.
  towergeom->Reset();

  towergeom->set_phibins(new_phibins);
  towergeom->set_etabins(new_etabins);

  for (int ibin = 0; ibin < new_phibins; ibin++)
  {
    towergeom->set_phibounds(ibin, phi_bound_map[ibin]);
  }
  for (int ibin = 0; ibin < new_etabins; ibin++)
  {
    towergeom->set_etabounds(ibin, eta_bound_map[ibin]);
  }

  // setup location of all towers
  for (int iphi = 0; iphi < towergeom->get_phibins(); iphi++)
    for (int ieta = 0; ieta < towergeom->get_etabins(); ieta++)
    {
      RawTowerGeomv1 *tg = new RawTowerGeomv1(
          RawTowerDefs::encode_towerid(caloid, ieta, iphi));

      CLHEP::Hep3Vector tower_pos;
      tower_pos.setRhoPhiEta(r, towergeom->get_phicenter(iphi), towergeom->get_etacenter(ieta));

      tg->set_center_x(tower_pos.x());
      tg->set_center_y(tower_pos.y());
      tg->set_center_z(tower_pos.z());

      towergeom->add_tower_geometry(tg);
    }
  if (Verbosity() >= VERBOSITY_SOME)
  {
    towergeom->identify();
  }

  const std::string input_TowerNodeName = "TOWER_" + _tower_node_prefix + "_" + detector;
  _towers = findNode::getClass<RawTowerContainer>(topNode,
                                                  input_TowerNodeName);
  if (!_towers)
  {
    std::cerr << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
              << " " << input_TowerNodeName << " Node missing, doing bail out!"
              << std::endl;
    throw std::runtime_error(
        "Failed to find " + input_TowerNodeName + " node in RawTowerCombiner::CreateNodes");
  }

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst(
      "PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cerr << __PRETTY_FUNCTION__ << "DST Node missing, doing nothing."
              << std::endl;
    throw std::runtime_error(
        "Failed to find DST node in RawTowerCombiner::CreateNodes");
  }

  return;
}
