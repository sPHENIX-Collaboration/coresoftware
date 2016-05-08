#include "RawTowerCombiner.h"
#include "RawTowerContainer.h"
#include "RawTowerGeomContainer_Cylinderv1.h"
#include "RawTowerGeomv1.h"
#include "RawTowerv1.h"
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellContainer.h>
#include <g4detectors/PHG4CylinderCell.h>
#include <g4detectors/PHG4CylinderCellDefs.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <g4main/PHG4Utils.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/getClass.h>

#include <iostream>
#include <sstream>      // std::stringstream
#include <stdexcept>
#include <map>

using namespace std;

RawTowerCombiner::RawTowerCombiner(const std::string& name) :
    SubsysReco(name), //
    /*std::string*/_tower_node_prefix("SIM"),
    /*std::string*/_output_node_suffix(""), // empty suffix will be automatically assigned later
    /*unsigned int*/_n_combine_eta(2),
    /*unsigned int*/_n_combine_phi(2),
    /*RawTowerContainer**/_intput_towers(NULL),
    /*RawTowerContainer**/_output_towers(NULL),
    /*std::string*/detector("NONE")
{

}

int
RawTowerCombiner::InitRun(PHCompositeNode *topNode)
{
  if (_n_combine_eta == 0)
    {
      cout << __PRETTY_FUNCTION__ << " Fatal error _n_combine_eta==0" << endl;

      exit(1243);
    }

  if (_n_combine_phi == 0)
    {
      cout << __PRETTY_FUNCTION__ << " Fatal error _n_combine_phi==0" << endl;

      exit(1243);
    }

  if (_output_node_suffix.length() == 0)
    {
      // automatic suffix

      stringstream suffix;
      suffix << _n_combine_eta << "x" << _n_combine_phi;
      _output_node_suffix = suffix.str();
    }

  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode",
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
  catch (std::exception& e)
    {
      std::cout << e.what() << std::endl;
      //exit(1);
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

int
RawTowerCombiner::process_event(PHCompositeNode *topNode)
{
  assert(_intput_towers);
  assert(_output_towers);

  if (verbosity)
    {
      std::cout << __PRETTY_FUNCTION__ << "Process event entered" << std::endl;
    }

  RawTowerContainer::ConstRange all_towers = _intput_towers->getTowers();
  for (RawTowerContainer::ConstIterator it = all_towers.first;
      it != all_towers.second; ++it)
    {
      const RawTower *input_tower = it->second;
      assert(input_tower);

      const int intput_eta = input_tower->get_bineta();
      const int intput_phi = input_tower->get_binphi();
      const int output_eta = get_output_bin_eta(intput_eta);
      const int output_phi = get_output_bin_phi(intput_phi);

      RawTower *output_tower = _output_towers->getTower(output_eta, output_phi);

      if (not output_tower)
        {
          output_tower = new RawTowerv1(*input_tower);
          _output_towers->AddTower(intput_eta, intput_phi, output_tower);
          if (verbosity >= VERBOSITY_MORE)
            {
              std::cout << __PRETTY_FUNCTION__ << "::" << detector << "::"
                  << " new output tower " << std::endl;
              _output_towers->identify();
            }
        }
      else
        {
          output_tower->set_energy(
              output_tower->get_energy() + input_tower->get_energy());

          output_tower->set_time(
              (abs(output_tower->get_energy()) * output_tower->get_time()
                  + abs(input_tower->get_energy()) * input_tower->get_time()) //
                  / (abs(output_tower->get_energy())
                      + abs(input_tower->get_energy()) + 1e-9) //avoid devide 0
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

          if (verbosity >= VERBOSITY_MORE)
            {
              std::cout << __PRETTY_FUNCTION__ << "::" << detector << "::"
                  << " merget into output tower " << std::endl;
              _output_towers->identify();
            }
        }

    }

  if (verbosity)
    {
      std::cout << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
          << "input sum energy = " << _intput_towers->getTotalEdep()
          << ", output sum digitalized value = "
          << _output_towers->getTotalEdep() << std::endl;
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

int
RawTowerCombiner::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void
RawTowerCombiner::CreateNodes(PHCompositeNode *topNode)
{
  assert(_output_node_suffix.length() > 0);

  PHNodeIterator iter(topNode);
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst(
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

  const string iTowerGeomNodeName = "TOWERGEOM_" + detector;
  RawTowerGeomContainer *input_towergeom = findNode::getClass<
      RawTowerGeomContainer>(topNode, iTowerGeomNodeName.c_str());
  if (!input_towergeom)
    {
      std::cerr << __PRETTY_FUNCTION__ << " - " << iTowerGeomNodeName
          << " Node missing, doing nothing." << std::endl;
      throw std::runtime_error(
          "Failed to find input tower geometry node in RawTowerCombiner::CreateNodes");
    }

  const string oTowerGeomNodeName = "TOWERGEOM_" + detector
      + _output_node_suffix;
  RawTowerGeomContainer *output_towergeom = findNode::getClass<
      RawTowerGeomContainer>(topNode, oTowerGeomNodeName.c_str());
  if (!output_towergeom)
    {
      output_towergeom = new RawTowerGeomContainer_Cylinderv1(caloid);
      PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(
          output_towergeom, oTowerGeomNodeName.c_str(), "PHObject");
      runNode->addNode(newNode);
    }

  assert(output_towergeom);

  const double r = input_towergeom->get_radius();
  output_towergeom->set_radius(r);
  output_towergeom->set_thickness(input_towergeom->get_thickness());
  output_towergeom->set_phibins(
      ceil(double(input_towergeom->get_phibins()) / double(_n_combine_phi)));
  output_towergeom->set_etabins(
      ceil(double(input_towergeom->get_etabins()) / double(_n_combine_eta)));

  for (int ibin = 0; ibin < output_towergeom->get_phibins(); ibin++)
    {
      const int first_bin = ibin * _n_combine_phi;
      assert(first_bin >= 0 && first_bin < input_towergeom->get_phibins());

      int last_bin = (ibin + 1) * _n_combine_phi - 1;
      if (last_bin >= input_towergeom->get_phibins())
        last_bin = input_towergeom->get_phibins();

      const pair<double, double> range1 = input_towergeom->get_phibounds(
          first_bin);
      const pair<double, double> range2 = input_towergeom->get_phibounds(
          last_bin);

      output_towergeom->set_phibounds(ibin,
          make_pair(range1.first, range2.second));
    }

  for (int ibin = 0; ibin < output_towergeom->get_etabins(); ibin++)
    {
      const int first_bin = ibin * _n_combine_eta;
      assert(first_bin >= 0 && first_bin < input_towergeom->get_etabins());

      int last_bin = (ibin + 1) * _n_combine_eta - 1;
      if (last_bin >= input_towergeom->get_etabins())
        last_bin = input_towergeom->get_etabins();

      const pair<double, double> range1 = input_towergeom->get_etabounds(
          first_bin);
      const pair<double, double> range2 = input_towergeom->get_etabounds(
          last_bin);

      output_towergeom->set_etabounds(ibin,
          make_pair(range1.first, range2.second));

    }

  // setup location of all towers
  for (int iphi = 0; iphi < output_towergeom->get_phibins(); iphi++)
    for (int ieta = 0; ieta < output_towergeom->get_etabins(); ieta++)
      {
        RawTowerGeomv1 * tg = new RawTowerGeomv1(
            RawTowerDefs::encode_towerid(caloid, ieta, iphi));

        tg->set_center_x(r * cos(output_towergeom->get_phicenter(iphi)));
        tg->set_center_y(r * sin(output_towergeom->get_phicenter(iphi)));
        tg->set_center_z(
            r
                / tan(
                    PHG4Utils::get_theta(
                        output_towergeom->get_etacenter(ieta))));
        output_towergeom->add_tower_geometry(tg);
      }
  if (verbosity >= VERBOSITY_SOME)
    {
      output_towergeom->identify();
    }

  const string input_TowerNodeName = "TOWER_" + _tower_node_prefix + "_"
      + detector;
  _intput_towers = findNode::getClass<RawTowerContainer>(topNode,
      input_TowerNodeName.c_str());
  if (!_intput_towers)
    {
      std::cerr << Name() << "::" << detector << "::" << __PRETTY_FUNCTION__
          << " " << input_TowerNodeName << " Node missing, doing bail out!"
          << std::endl;
      throw std::runtime_error(
          "Failed to find " + input_TowerNodeName
              + " node in RawTowerCombiner::CreateNodes");
    }

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst(
      "PHCompositeNode", "DST"));
  if (!dstNode)
    {
      std::cerr << __PRETTY_FUNCTION__ << "DST Node missing, doing nothing."
          << std::endl;
      throw std::runtime_error(
          "Failed to find DST node in RawTowerCombiner::CreateNodes");
    }

  PHNodeIterator dstiter(dstNode);
  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode*>(dstiter.findFirst(
      "PHCompositeNode", detector));
  if (!DetNode)
    {
      DetNode = new PHCompositeNode(detector);
      dstNode->addNode(DetNode);
    }

  // Create the tower nodes on the tree
  _output_towers = new RawTowerContainer(caloid);
  stringstream soutput_towers;
  soutput_towers << "TOWER_";
  if (_tower_node_prefix.length() > 0)
    {
      soutput_towers << _tower_node_prefix << "_";
    }
  soutput_towers << detector << _output_node_suffix;
  PHIODataNode<PHObject> *towerNode = new PHIODataNode<PHObject>(_output_towers,
      soutput_towers.str(), "PHObject");
  DetNode->addNode(towerNode);

  if (verbosity)
    {
      cout << __PRETTY_FUNCTION__ << ": output to " << soutput_towers.str()
          << endl;
    }

  return;
}

