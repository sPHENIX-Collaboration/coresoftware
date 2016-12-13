#include "HcalRawTowerBuilder.h"
#include "RawTowerContainer.h"
#include "RawTowerGeomContainer_Cylinderv1.h"
#include "RawTowerGeomv1.h"
#include "RawTowerv1.h"

#include <g4detectors/PHG4HcalDefs.h>
#include <g4detectors/PHG4Parameters.h>
#include <g4detectors/PHG4ScintillatorSlat.h>
#include <g4detectors/PHG4ScintillatorSlatContainer.h>
#include <g4detectors/PHG4ScintillatorSlatDefs.h>

#include <g4main/PHG4Utils.h>

#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <pdbcalbase/PdbParameterMap.h>
#include <pdbcalbase/PdbParameterMapContainer.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>

#include <TSystem.h>

#include <iostream>
#include <map>
#include <stdexcept>

using namespace std;

HcalRawTowerBuilder::HcalRawTowerBuilder(const std::string& name) :
  SubsysReco(name), 
  PHG4ParameterInterface(name),
  _towers(NULL), 
  rawtowergeom(NULL),
  detector("NONE"), 
  emin(NAN),
  chkenergyconservation(0), 
  _tower_energy_src(enu_tower_energy_src::unknown),
  ncell_to_tower(-1),
  _timer(PHTimeServer::get()->insert_new(name))
{
  InitializeParameters();
}

int
HcalRawTowerBuilder::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode",
      "DST"));
  if (!dstNode)
    {
      std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
      gSystem->Exit(1);
    }
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN" ));
  string paramnodename = "TOWERPARAM_" + detector;


  try
    {
      CreateNodes(topNode);
    }
  catch (std::exception& e)
    {
      std::cout << e.what() << std::endl;
      //exit(1);
    }

  // order first default, 
  // then parameter from g4detector on node tree
   ReadParamsFromNodeTree(topNode);
  // then macro setting
  UpdateParametersWithMacro();
  PHNodeIterator runIter(runNode);
  PHCompositeNode *RunDetNode =  dynamic_cast<PHCompositeNode*>(runIter.findFirst("PHCompositeNode",detector));
  if (! RunDetNode)
    {
      RunDetNode = new PHCompositeNode(detector);
      runNode->addNode(RunDetNode);
    }
  SaveToNodeTree(RunDetNode,paramnodename);
  _tower_energy_src = get_int_param("tower_energy_source");
  emin = get_double_param("emin");
  ncell_to_tower = get_int_param("n_scinti_plates_per_tower");
  if (verbosity >= 1)
    {
      cout << "HcalRawTowerBuilder::InitRun :";
      if (_tower_energy_src == kEnergyDeposition)
	{
        cout << "save Geant4 energy deposition in towers" << endl;
	}
      else if (_tower_energy_src == kLightYield)
	{
        cout << "save light yield in towers" << endl;
	}
      else if (_tower_energy_src == kIonizationEnergy)
	{
        cout << "save ionization energy in towers" << endl;
	}
      else
	{
	  cout << "unknown energy source" << endl;
	}
    }
  TowerGeomNodeName = "TOWERGEOM_" + detector;
  rawtowergeom = findNode::getClass<RawTowerGeomContainer>(topNode,
      TowerGeomNodeName.c_str());
  if (!rawtowergeom)
    {
      rawtowergeom = new RawTowerGeomContainer_Cylinderv1(RawTowerDefs::convert_name_to_caloid(detector));
      PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(rawtowergeom,
          TowerGeomNodeName.c_str(), "PHObject");
      RunDetNode->addNode(newNode);
    }
  rawtowergeom->set_radius(get_double_param(PHG4HcalDefs::innerrad));
  rawtowergeom->set_thickness(get_double_param(PHG4HcalDefs::outerrad)-get_double_param(PHG4HcalDefs::innerrad));
  rawtowergeom->set_phibins(get_int_param(PHG4HcalDefs::n_towers));
  rawtowergeom->set_etabins(get_int_param("etabins"));
  for (int i=0; i<get_int_param(PHG4HcalDefs::n_towers); i++)
    {
      pair<double, double> range = make_pair(i*360./get_int_param(PHG4HcalDefs::n_towers),(i+1)*360./get_int_param(PHG4HcalDefs::n_towers));
      rawtowergeom->set_phibounds(i,range);
    }
  double etalowbound = -1.1;
  for (int i = 0; i < get_int_param("etabins"); i++)
    {
      double etahibound = etalowbound + 2.2 / get_int_param("etabins");
      pair<double, double> range = make_pair(etalowbound, etahibound);
      rawtowergeom->set_etabounds(i, range);
      etalowbound = etahibound;
    }
	 //  if (verbosity > 0)
    {
      rawtowergeom->identify();
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

int
HcalRawTowerBuilder::process_event(PHCompositeNode *topNode)
{
  if (verbosity>3)
    {
      std::cout << PHWHERE << "Process event entered" << std::endl;
    }

  // get cells
  std::string cellnodename = "G4CELL_" + detector;
  PHG4ScintillatorSlatContainer* slats = findNode::getClass<PHG4ScintillatorSlatContainer>(topNode, cellnodename.c_str());
  if (!slats)
    {
      std::cerr << PHWHERE << " " << cellnodename
          << " Node missing, doing nothing." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  // loop over all slats in an event
  PHG4ScintillatorSlatContainer::ConstIterator cell_iter;
  PHG4ScintillatorSlatContainer::ConstRange cell_range = slats->getScintillatorSlats();
  for (cell_iter = cell_range.first; cell_iter != cell_range.second;
      ++cell_iter)
    {
      PHG4ScintillatorSlat *cell = cell_iter->second;

      if (verbosity > 2)
        {
          std::cout << PHWHERE << " print out the cell:" << std::endl;
          cell->identify();
        }
      short twrrow = get_tower_row(cell->get_row());
      // add the energy to the corresponding tower
      // towers are addressed column/row to make the mapping more intuitive
      RawTower *tower = _towers->getTower(cell->get_column(), twrrow);
      if (!tower)
        {
          tower = new RawTowerv1();
          tower->set_energy(0);
          _towers->AddTower(cell->get_column(), twrrow, tower);
        }
      double cell_weight = 0;
      if (_tower_energy_src == kEnergyDeposition)
	{
	  cell_weight = cell->get_edep();
	}
      else if (_tower_energy_src == kLightYield)
	{
	  cell_weight = cell->get_light_yield();
	}
      else if (_tower_energy_src == kIonizationEnergy)
	{
	  cell_weight = cell->get_eion();
	}
      else
	{
	  cout << Name() << ": unknown tower energy source " 
	       << _tower_energy_src << endl;
	  gSystem->Exit(1);
	}

      tower->add_ecell(cell->get_key(), cell_weight);

      tower->set_energy(tower->get_energy() + cell_weight);

    }
  double towerE = 0;
  if (chkenergyconservation)
    {
      double cellE = slats->getTotalEdep();
      towerE = _towers->getTotalEdep();
      if (fabs(cellE - towerE) / cellE > 1e-5)
        {
          cout << "towerE: " << towerE << ", cellE: " << cellE << ", delta: "
              << cellE - towerE << endl;
        }
    }
  if (verbosity)
    {
      towerE = _towers->getTotalEdep();
    }

  _towers->compress(emin);
  if (verbosity)
    {
      cout << "Energy lost by dropping towers with less than " << emin
          << " energy, lost energy: " << towerE - _towers->getTotalEdep()
          << endl;
      _towers->identify();
      RawTowerContainer::ConstRange begin_end = _towers->getTowers();
      RawTowerContainer::ConstIterator iter;
      for (iter = begin_end.first; iter != begin_end.second; ++iter)
        {
          iter->second->identify();
        }
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

void
HcalRawTowerBuilder::CreateNodes(PHCompositeNode *topNode)
{

  PHNodeIterator iter(topNode);
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst(
      "PHCompositeNode", "RUN"));
  if (!runNode)
    {
      std::cerr << PHWHERE << "Run Node missing, doing nothing." << std::endl;
      throw std::runtime_error(
          "Failed to find Run node in HcalRawTowerBuilder::CreateNodes");
    }
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst(
      "PHCompositeNode", "DST"));
  if (!dstNode)
    {
      std::cerr << PHWHERE << "DST Node missing, doing nothing." << std::endl;
      throw std::runtime_error(
          "Failed to find DST node in HcalRawTowerBuilder::CreateNodes");
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
  _towers = new RawTowerContainer(RawTowerDefs::convert_name_to_caloid(detector));
  if (_sim_tower_node_prefix.length() == 0)
    {
      // no prefix, consistent with older convension
      TowerNodeName = "TOWER_" + detector;
    }
  else
    {
      TowerNodeName = "TOWER_" + _sim_tower_node_prefix + "_" + detector;
    }
  PHIODataNode<PHObject> *towerNode = new PHIODataNode<PHObject>(_towers,
      TowerNodeName.c_str(), "PHObject");
  DetNode->addNode(towerNode);
  return;
}

short
HcalRawTowerBuilder::get_tower_row(const short cellrow) const
{
  short twrrow = cellrow/ncell_to_tower;
  return twrrow;
}

void
HcalRawTowerBuilder::SetDefaultParameters()
{
  set_default_int_param(PHG4HcalDefs::scipertwr,5);
  set_default_int_param("tower_energy_source",kLightYield);
  set_default_int_param(PHG4HcalDefs::n_towers,64);
  set_default_int_param("etabins",24);

  set_default_double_param("emin",1.e-6);
  set_default_double_param(PHG4HcalDefs::outerrad,NAN);
  set_default_double_param(PHG4HcalDefs::innerrad,NAN);
}

void
HcalRawTowerBuilder::ReadParamsFromNodeTree(PHCompositeNode *topNode)
{
  PHG4Parameters *pars = new PHG4Parameters("temp");
  // we need the number of scintillator plates per tower
  string geonodename = "G4GEOPARAM_" + detector;
  PdbParameterMapContainer *saveparams = findNode::getClass<PdbParameterMapContainer>(topNode,geonodename);
  if (! saveparams)
    {
      cout << "could not find " << geonodename << endl;
      Fun4AllServer *se = Fun4AllServer::instance();
      se->Print("NODETREE");
      return;
    }
  pars->FillFrom(saveparams,0);
  set_int_param(PHG4HcalDefs::scipertwr,pars->get_int_param(PHG4HcalDefs::scipertwr));
  set_int_param(PHG4HcalDefs::n_towers,pars->get_int_param(PHG4HcalDefs::n_towers));
  set_int_param("etabins",2*pars->get_int_param(PHG4HcalDefs::n_scinti_tiles));
  set_double_param(PHG4HcalDefs::innerrad,pars->get_double_param(PHG4HcalDefs::innerrad));
  set_double_param(PHG4HcalDefs::outerrad,pars->get_double_param(PHG4HcalDefs::outerrad));
  delete pars;
  return;
}
