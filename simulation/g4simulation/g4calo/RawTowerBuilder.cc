#include "RawTowerBuilder.h"
#include "RawTowerContainer.h"
#include "RawTowerGeomContainer_Cylinderv1.h"
#include "RawTowerGeomv1.h"
#include "RawTowerv1.h"

#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderCellGeom.h>

#include <g4detectors/PHG4CellContainer.h>
#include <g4detectors/PHG4Cell.h>
#include <g4detectors/PHG4CellDefs.h>

#include <g4main/PHG4Utils.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <iostream>
#include <stdexcept>
#include <map>

using namespace std;

RawTowerBuilder::RawTowerBuilder(const std::string& name) :
  SubsysReco(name),
  _towers(nullptr),
  rawtowergeom(nullptr),
  detector("NONE"),
  _cell_binning(PHG4CellDefs::undefined),
  emin(1e-6),
  chkenergyconservation(0),
  _nlayers(-1),
  _nphibins(-1),
  _netabins(-1),
  _etamin(NAN),
  _phimin(NAN),
  _etastep(NAN),
  _phistep(NAN),
  _tower_energy_src(kLightYield),
  _timer(PHTimeServer::get()->insert_new(name))
{}

int
RawTowerBuilder::InitRun(PHCompositeNode *topNode)
{
  std::string geonodename = "CYLINDERCELLGEOM_" + detector;
  PHG4CylinderCellGeomContainer* cellgeos = findNode::getClass <
  PHG4CylinderCellGeomContainer > (topNode, geonodename.c_str());
  if (!cellgeos)
    {
      cout << PHWHERE << " " << geonodename
	   << " Node missing, doing nothing." << std::endl;
      throw std::runtime_error(
			       "Failed to find " + geonodename
			       + " node in RawTowerBuilder::CreateNodes");
    }

  // fill the number of layers in the calorimeter
  _nlayers = cellgeos->get_NLayers();

  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
    {
      std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
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

  if (verbosity >= 1)
    {
      cout << "RawTowerBuilder::InitRun :";
      if (_tower_energy_src == kEnergyDeposition)
	{
	  cout << "save Geant4 energy deposition as the weight of the cells"
	       << endl;
	}
      else if (_tower_energy_src == kLightYield)
	{
	  cout << "save light yield as the weight of the cells" << endl;
	}
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

int
RawTowerBuilder::process_event(PHCompositeNode *topNode)
{
  if (verbosity)
    {
      std::cout << PHWHERE << "Process event entered" << std::endl;
    }

  // get cells
  std::string cellnodename = "G4CELL_" + detector;
  PHG4CellContainer* cells = findNode::getClass<PHG4CellContainer>(topNode, cellnodename.c_str());
  if (!cells)
    {
      std::cerr << PHWHERE << " " << cellnodename
		<< " Node missing, doing nothing." << std::endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  // loop over all cells in an event
  PHG4CellContainer::ConstIterator cell_iter;
  PHG4CellContainer::ConstRange cell_range = cells->getCells();
  for (cell_iter = cell_range.first; cell_iter != cell_range.second; ++cell_iter)
    {
      PHG4Cell *cell = cell_iter->second;

      if (verbosity > 2)
	{
	  std::cout << PHWHERE << " print out the cell:" << std::endl;
	  cell->identify();
	}

      // add the energy to the corresponding tower
      RawTower *tower = nullptr;
      short int firstpar;
      short int secondpar;
      if (cell->has_binning(PHG4CellDefs::sizebinning))
	{
	  firstpar = PHG4CellDefs::SizeBinning::get_zbin(cell->get_cellid());
	  secondpar = PHG4CellDefs::SizeBinning::get_phibin(cell->get_cellid());
	}
      else if (cell->has_binning(PHG4CellDefs::spacalbinning))
	{
	  firstpar = PHG4CellDefs::SpacalBinning::get_etabin(cell->get_cellid());
	  secondpar =  PHG4CellDefs::SpacalBinning::get_phibin(cell->get_cellid());
	}
      else
	{
	  cout << "unknown cell binning, implement 0x" << hex << PHG4CellDefs::get_binning(cell->get_cellid()) << endl;
	  exit(1);
	}
      tower = _towers->getTower(firstpar,secondpar);
      if (!tower)
	{
	  tower = new RawTowerv1();
	  tower->set_energy(0);
	  _towers->AddTower(firstpar,secondpar, tower);
	}
      float cell_weight = 0;
      if (_tower_energy_src == kEnergyDeposition)
	{
	  cell_weight = cell->get_edep();
	}
      else if (_tower_energy_src == kLightYield)
	{
	  cell_weight = cell->get_light_yield();
	}

      tower->add_ecell(cell->get_cellid(), cell_weight);

      PHG4Cell::ShowerEdepConstRange range = cell->get_g4showers();
      for (PHG4Cell::ShowerEdepConstIterator shower_iter = range.first;
	   shower_iter != range.second;
	   ++shower_iter)
	{
	  tower->add_eshower(shower_iter->first, shower_iter->second);
	}

      tower->set_energy(tower->get_energy() + cell_weight);

      if (verbosity > 2)
	{
	  rawtowergeom = findNode::getClass<RawTowerGeomContainer>(topNode, TowerGeomNodeName.c_str());
	  tower->identify();
	}
    }
  double towerE = 0;
  if (chkenergyconservation)
    {
      double cellE = cells->getTotalEdep();
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
RawTowerBuilder::CreateNodes(PHCompositeNode *topNode)
{

  PHNodeIterator iter(topNode);
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
    {
      std::cerr << PHWHERE << "Run Node missing, doing nothing." << std::endl;
      throw std::runtime_error("Failed to find Run node in RawTowerBuilder::CreateNodes");
    }

  PHNodeIterator runIter(runNode);
  PHCompositeNode *RunDetNode =  dynamic_cast<PHCompositeNode*>(runIter.findFirst("PHCompositeNode",detector));
  if (! RunDetNode)
    {
      RunDetNode = new PHCompositeNode(detector);
      runNode->addNode(RunDetNode);
    }

  const RawTowerDefs::CalorimeterId caloid = RawTowerDefs::convert_name_to_caloid(detector);

  // get the cell geometry and build up the tower geometry object
  std::string geonodename = "CYLINDERCELLGEOM_" + detector;
  PHG4CylinderCellGeomContainer* cellgeos = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, geonodename.c_str());
  if (!cellgeos)
    {
      std::cerr << PHWHERE << " " << geonodename
		<< " Node missing, doing nothing." << std::endl;
      throw std::runtime_error(
			       "Failed to find " + geonodename
			       + " node in RawTowerBuilder::CreateNodes");
    }
  TowerGeomNodeName = "TOWERGEOM_" + detector;
  rawtowergeom = findNode::getClass<RawTowerGeomContainer>(topNode, TowerGeomNodeName.c_str());
  if (!rawtowergeom)
    {
      rawtowergeom = new RawTowerGeomContainer_Cylinderv1(caloid);
      PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(rawtowergeom, TowerGeomNodeName.c_str(), "PHObject");
      RunDetNode->addNode(newNode);
    }
  // fill the number of layers in the calorimeter
  _nlayers = cellgeos->get_NLayers();

  // Create the tower nodes on the tree
  PHG4CylinderCellGeomContainer::ConstIterator miter;
  PHG4CylinderCellGeomContainer::ConstRange begin_end =
    cellgeos->get_begin_end();
  int ifirst = 1;
  int first_layer = -1;
  PHG4CylinderCellGeom * first_cellgeo = nullptr;
  double inner_radius = 0;
  double thickness = 0;
  for (miter = begin_end.first; miter != begin_end.second; ++miter)
    {
      PHG4CylinderCellGeom *cellgeo = miter->second;
      first_cellgeo = miter->second;

      if (verbosity)
	{
	  cellgeo->identify();
	}
      thickness += cellgeo->get_thickness();
      if (ifirst)
	{
	  _cell_binning = cellgeo->get_binning();
	  _nphibins = cellgeo->get_phibins();
	  _phimin = cellgeo->get_phimin();
	  _phistep = cellgeo->get_phistep();
	  if (_cell_binning == PHG4CellDefs::etaphibinning
	      || _cell_binning == PHG4CellDefs::etaslatbinning)
	    {
	      _netabins = cellgeo->get_etabins();
	      _etamin = cellgeo->get_etamin();
	      _etastep = cellgeo->get_etastep();
	    }
	  else if (_cell_binning == PHG4CellDefs::sizebinning)
	    {
	      _netabins = cellgeo->get_zbins(); // bin eta in the same number of z bins
	    }
	  else if (_cell_binning == PHG4CellDefs::spacalbinning)
	    {
	      // use eta definiton for each row of towers
	      _netabins = cellgeo->get_etabins();
	    }
	  else
	    {
	      cout << "RawTowerBuilder::CreateNodes::" << Name()
		   << " - Fatal Error - unsupported cell binning method "
		   << _cell_binning << endl;
	    }
	  inner_radius = cellgeo->get_radius();
	  first_layer = cellgeo->get_layer();
	  ifirst = 0;
	}
      else
	{
	  if (_cell_binning != cellgeo->get_binning())
	    {
	      cout << "inconsistent binning method from " << _cell_binning
		   << " layer " << cellgeo->get_layer() << ": "
		   << cellgeo->get_binning() << endl;
	      exit(1);
	    }
	  if (inner_radius > cellgeo->get_radius())
	    {
	      cout << "radius of layer " << cellgeo->get_layer() << " is "
		   << cellgeo->get_radius() << " which smaller than radius "
		   << inner_radius << " of first layer in list " << first_layer
		   << endl;
	    }
	  if (_nphibins != cellgeo->get_phibins())
	    {
	      cout << "mixing number of phibins, fisrt layer: " << _nphibins
		   << " layer " << cellgeo->get_layer() << ": "
		   << cellgeo->get_phibins() << endl;
	      exit(1);
	    }
	  if (_phimin != cellgeo->get_phimin())
	    {
	      cout << "mixing number of phimin, fisrt layer: " << _phimin
		   << " layer " << cellgeo->get_layer() << ": "
		   << cellgeo->get_phimin() << endl;
	      exit(1);
	    }
	  if (_phistep != cellgeo->get_phistep())
	    {
	      cout << "mixing phisteps first layer: " << _phistep << " layer "
		   << cellgeo->get_layer() << ": " << cellgeo->get_phistep()
		   << " diff: " << _phistep - cellgeo->get_phistep() << endl;
	      exit(1);
	    }
	  if (_cell_binning == PHG4CellDefs::etaphibinning
	      || _cell_binning == PHG4CellDefs::etaslatbinning)
	    {
	      if (_netabins != cellgeo->get_etabins())
		{
		  cout << "mixing number of_netabins , fisrt layer: "
		       << _netabins << " layer " << cellgeo->get_layer() << ": "
		       << cellgeo->get_etabins() << endl;
		  exit(1);
		}
	      if (fabs(_etamin - cellgeo->get_etamin()) > 1e-9)
		{
		  cout << "mixing etamin, fisrt layer: " << _etamin << " layer "
		       << cellgeo->get_layer() << ": " << cellgeo->get_etamin()
		       << " diff: " << _etamin - cellgeo->get_etamin() << endl;
		  exit(1);
		}
	      if (fabs(_etastep - cellgeo->get_etastep()) > 1e-9)
		{
		  cout << "mixing eta steps first layer: " << _etastep
		       << " layer " << cellgeo->get_layer() << ": "
		       << cellgeo->get_etastep() << " diff: "
		       << _etastep - cellgeo->get_etastep() << endl;
		  exit(1);
		}
	    }

	  else if (_cell_binning == PHG4CellDefs::sizebinning)
	    {

	      if (_netabins != cellgeo->get_zbins())
		{
		  cout << "mixing number of z bins , fisrt layer: " << _netabins
		       << " layer " << cellgeo->get_layer() << ": "
		       << cellgeo->get_zbins() << endl;
		  exit(1);
		}
	    }
	}
    }
  rawtowergeom->set_radius(inner_radius);
  rawtowergeom->set_thickness(thickness);
  rawtowergeom->set_phibins(_nphibins);
  //  rawtowergeom->set_phistep(_phistep);
  //  rawtowergeom->set_phimin(_phimin);
  rawtowergeom->set_etabins(_netabins);

  if (!first_cellgeo)
    {
      cout << "RawTowerBuilder::CreateNodes - ERROR - can not find first layer of cells "
	   << endl;

      exit(1);
    }

  for (int ibin = 0; ibin < first_cellgeo->get_phibins(); ibin++)
    {

      const pair<double, double> range = first_cellgeo->get_phibounds(ibin);

      rawtowergeom->set_phibounds(ibin, range);
    }

  if (_cell_binning == PHG4CellDefs::etaphibinning
      || _cell_binning == PHG4CellDefs::etaslatbinning
      || _cell_binning == PHG4CellDefs::spacalbinning)
    {
      const double r = inner_radius + thickness / 2.;

      for (int ibin = 0; ibin < first_cellgeo->get_etabins(); ibin++)
	{

	  const pair<double, double> range = first_cellgeo->get_etabounds(ibin);

	  rawtowergeom->set_etabounds(ibin, range);
	}

      // setup location of all towers
      for (int iphi = 0; iphi < rawtowergeom->get_phibins(); iphi++)
	{
	  for (int ieta = 0; ieta < rawtowergeom->get_etabins(); ieta++)
	    {
	      RawTowerGeomv1 * tg = new RawTowerGeomv1(RawTowerDefs::encode_towerid(caloid, ieta, iphi));

	      tg->set_center_x(r * cos(rawtowergeom->get_phicenter(iphi)));
	      tg->set_center_y(r * sin(rawtowergeom->get_phicenter(iphi)));
	      tg->set_center_z(r / tan(PHG4Utils::get_theta(rawtowergeom->get_etacenter(ieta))));
	      rawtowergeom->add_tower_geometry(tg);
	    }
	}
    }
  else if (_cell_binning == PHG4CellDefs::sizebinning)
    {
      const double r = inner_radius + thickness / 2.;

      for (int ibin = 0; ibin < first_cellgeo->get_zbins(); ibin++)
	{

	  const pair<double, double> z_range = first_cellgeo->get_zbounds(ibin);
	  //          const double r = first_cellgeo->get_radius();
	  const double eta1 = -log(tan(atan2(r, z_range.first) / 2.));
	  const double eta2 = -log(tan(atan2(r, z_range.second) / 2.));
	  rawtowergeom->set_etabounds(ibin, make_pair(eta1, eta2));
	}

      // setup location of all towers
      for (int iphi = 0; iphi < rawtowergeom->get_phibins(); iphi++)
	{
	  for (int ibin = 0; ibin < first_cellgeo->get_zbins(); ibin++)
	    {
	      RawTowerGeomv1 * tg = new RawTowerGeomv1(RawTowerDefs::encode_towerid(caloid, ibin, iphi));

	      tg->set_center_x(r * cos(rawtowergeom->get_phicenter(iphi)));
	      tg->set_center_y(r * sin(rawtowergeom->get_phicenter(iphi)));
	      tg->set_center_z(first_cellgeo->get_zcenter(ibin));

	      rawtowergeom->add_tower_geometry(tg);
	    }
	}
    }
  else
    {
      cout << "RawTowerBuilder::CreateNodes - ERROR - unsupported cell geometry "
	   << _cell_binning << endl;
      exit(1);
    }

  if (verbosity >= 1)
    {
      rawtowergeom->identify();
    }

  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
    {
      std::cerr << PHWHERE << "DST Node missing, doing nothing." << std::endl;
      throw std::runtime_error(
			       "Failed to find DST node in RawTowerBuilder::CreateNodes");
    }

  PHNodeIterator dstiter(dstNode);
  PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode*>(dstiter.findFirst("PHCompositeNode", detector));
  if (!DetNode)
    {
      DetNode = new PHCompositeNode(detector);
      dstNode->addNode(DetNode);
    }

  // Create the tower nodes on the tree
  _towers = new RawTowerContainer(caloid);
  if (_sim_tower_node_prefix.length() == 0)
    {
      // no prefix, consistent with older convension
      TowerNodeName = "TOWER_" + detector;
    }
  else
    {
      TowerNodeName = "TOWER_" + _sim_tower_node_prefix + "_" + detector;
    }
  PHIODataNode<PHObject> *towerNode = new PHIODataNode<PHObject>(_towers, TowerNodeName.c_str(), "PHObject");
  DetNode->addNode(towerNode);

  return;
}
