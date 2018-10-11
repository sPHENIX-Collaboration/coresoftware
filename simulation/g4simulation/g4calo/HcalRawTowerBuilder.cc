#include "HcalRawTowerBuilder.h"
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/RawTowerGeomv1.h>
#include <calobase/RawTowerv1.h>

#include <g4detectors/PHG4HcalDefs.h>
#include <phparameter/PHParameters.h>
#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellContainer.h>
#include <g4detectors/PHG4CellDefs.h>

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
  PHParameterInterface(name),
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
  PHCompositeNode *parNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "PAR" ));
  string geonodename = "TOWERGEO_" + detector;

  PHNodeIterator parIter(parNode);
  PHCompositeNode *ParDetNode =  dynamic_cast<PHCompositeNode*>(parIter.findFirst("PHCompositeNode",detector));
  if (! ParDetNode)
    {
      ParDetNode = new PHCompositeNode(detector);
      parNode->addNode(ParDetNode);
    }
  PutOnParNode(ParDetNode,geonodename);
  _tower_energy_src = get_int_param("tower_energy_source");
  emin = get_double_param("emin");
  ncell_to_tower = get_int_param("n_scinti_plates_per_tower");
  if (Verbosity() >= 1)
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
  double innerrad = get_double_param(PHG4HcalDefs::innerrad);
  double thickness = get_double_param(PHG4HcalDefs::outerrad)-innerrad;
  rawtowergeom->set_radius(innerrad);
  rawtowergeom->set_thickness(thickness);
  rawtowergeom->set_phibins(get_int_param(PHG4HcalDefs::n_towers));
  rawtowergeom->set_etabins(get_int_param("etabins"));
  double geom_ref_radius = innerrad + thickness/2.;
  double phistart = 0;
  for (int i=0; i<get_int_param(PHG4HcalDefs::n_towers); i++)
    {
      double phiend = phistart+2.*M_PI/get_int_param(PHG4HcalDefs::n_towers);
      pair<double, double> range = make_pair(phistart,phiend);
      phistart = phiend;
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
  for (int iphi = 0; iphi < rawtowergeom->get_phibins(); iphi++)
    {
      for (int ieta = 0; ieta < rawtowergeom->get_etabins(); ieta++)
	{
	  RawTowerGeomv1 * tg = new RawTowerGeomv1(RawTowerDefs::encode_towerid(RawTowerDefs::convert_name_to_caloid(detector), ieta, iphi));

	  tg->set_center_x(geom_ref_radius * cos(rawtowergeom->get_phicenter(iphi)));
	  tg->set_center_y(geom_ref_radius * sin(rawtowergeom->get_phicenter(iphi)));
	  tg->set_center_z(geom_ref_radius / tan(PHG4Utils::get_theta(rawtowergeom->get_etacenter(ieta))));
	  rawtowergeom->add_tower_geometry(tg);
	}
    }
  if (Verbosity() > 0)
    {
      rawtowergeom->identify();
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

int
HcalRawTowerBuilder::process_event(PHCompositeNode *topNode)
{
  if (Verbosity()>3)
    {
      std::cout << PHWHERE << "Process event entered" << std::endl;
    }

  // get cells
  std::string cellnodename = "G4CELL_" + detector;
  PHG4CellContainer* slats = findNode::getClass<PHG4CellContainer>(topNode, cellnodename.c_str());
  if (!slats)
    {
      std::cerr << PHWHERE << " " << cellnodename
          << " Node missing, doing nothing." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  // loop over all slats in an event
  PHG4CellContainer::ConstIterator cell_iter;
  PHG4CellContainer::ConstRange cell_range = slats->getCells();
  for (cell_iter = cell_range.first; cell_iter != cell_range.second;
      ++cell_iter)
    {
      PHG4Cell *cell = cell_iter->second;

      if (Verbosity() > 2)
        {
          std::cout << PHWHERE << " print out the cell:" << std::endl;
          cell->identify();
        }
      short twrrow = get_tower_row(PHG4CellDefs::ScintillatorSlatBinning::get_row(cell->get_cellid()));
      // add the energy to the corresponding tower
      // towers are addressed column/row to make the mapping more intuitive
      RawTower *tower = _towers->getTower(PHG4CellDefs::ScintillatorSlatBinning::get_column(cell->get_cellid()), twrrow);
      if (!tower)
        {
          tower = new RawTowerv1();
          tower->set_energy(0);
          _towers->AddTower(PHG4CellDefs::ScintillatorSlatBinning::get_column(cell->get_cellid()), twrrow, tower);
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

      tower->add_ecell(cell->get_cellid(), cell_weight);

      PHG4Cell::ShowerEdepConstRange range = cell->get_g4showers();
      for (PHG4Cell::ShowerEdepConstIterator shower_iter = range.first;
	   shower_iter != range.second;
	   ++shower_iter) 
	{
	  tower->add_eshower(shower_iter->first,shower_iter->second);
	}
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
  if (Verbosity())
    {
      towerE = _towers->getTotalEdep();
    }

  _towers->compress(emin);
  if (Verbosity())
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
  PHParameters *pars = new PHParameters("temp");
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
