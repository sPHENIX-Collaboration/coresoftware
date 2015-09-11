#include "RawTowerCalibration.h"
#include "RawTowerContainer.h"
#include "RawTowerGeomv2.h"
#include "RawTowerv1.h"
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderCellGeom.h>
#include <g4detectors/PHG4CylinderCellContainer.h>
#include <g4detectors/PHG4CylinderCell.h>
#include <g4detectors/PHG4CylinderCellDefs.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/getClass.h>

#include <iostream>
#include <stdexcept>
#include <map>

using namespace std;

RawTowerCalibration::RawTowerCalibration(const std::string& name):
  SubsysReco(name),
  _towers(NULL),
  rawtowergeom(NULL),
  detector("NONE"),
  _cell_binning(PHG4CylinderCellDefs::undefined),
  emin(1e-6),
  chkenergyconservation(0),
  _nlayers(-1),
  _nphibins(-1),
  _netabins(-1),
  _etamin(NAN),
  _phimin(NAN),
  _etastep(NAN),
  _phistep(NAN),
  _timer( PHTimeServer::get()->insert_new(name) )
{}

int 
RawTowerCalibration::InitRun(PHCompositeNode *topNode)
{
  std::string geonodename = "CYLINDERCELLGEOM_" + detector;
  PHG4CylinderCellGeomContainer* cellgeos = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, geonodename.c_str());
  if (!cellgeos)
    {
      std::cerr << PHWHERE << " " << geonodename << " Node missing, doing nothing." << std::endl;
      throw std::runtime_error("Failed to find " + geonodename + " node in RawTowerCalibration::CreateNodes");
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
  return Fun4AllReturnCodes::EVENT_OK;
}

int 
RawTowerCalibration::process_event(PHCompositeNode *topNode)
{
  if (verbosity)
    {
      std::cout << PHWHERE << "Process event entered" << std::endl;
    }



  return Fun4AllReturnCodes::EVENT_OK;
}

int 
RawTowerCalibration::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void 
RawTowerCalibration::CreateNodes(PHCompositeNode *topNode)
{

  PHNodeIterator iter(topNode);
  PHCompositeNode *runNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
    {
      std::cerr << PHWHERE << "Run Node missing, doing nothing." << std::endl;
      throw std::runtime_error("Failed to find Run node in RawTowerCalibration::CreateNodes");
    }

  // get the cell geometry and build up the tower geometry object
  std::string geonodename = "CYLINDERCELLGEOM_" + detector;
  PHG4CylinderCellGeomContainer* cellgeos = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode, geonodename.c_str());
  if (!cellgeos)
    {
      std::cerr << PHWHERE << " " << geonodename << " Node missing, doing nothing." << std::endl;
      throw std::runtime_error("Failed to find " + geonodename + " node in RawTowerCalibration::CreateNodes");
    }
  TowerGeomNodeName = "TOWERGEOM_" + detector;
  rawtowergeom =  findNode::getClass<RawTowerGeom>(topNode, TowerGeomNodeName.c_str());
  if (! rawtowergeom)
    {

      rawtowergeom = new RawTowerGeomv2();
      PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(rawtowergeom, TowerGeomNodeName.c_str(), "PHObject");
      runNode->addNode(newNode);
    }
  // fill the number of layers in the calorimeter
  _nlayers = cellgeos->get_NLayers();

  // Create the tower nodes on the tree
  std::map<int, PHG4CylinderCellGeom *>::const_iterator miter;
  std::pair <std::map<int, PHG4CylinderCellGeom *>::const_iterator, std::map<int, PHG4CylinderCellGeom *>::const_iterator> begin_end = cellgeos->get_begin_end();
  int ifirst = 1;
  int first_layer = -1;
  PHG4CylinderCellGeom * first_cellgeo = NULL;
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
          if (_cell_binning == PHG4CylinderCellDefs::etaphibinning
              or _cell_binning == PHG4CylinderCellDefs::etaslatbinning)
            {
              _netabins = cellgeo->get_etabins();
              _etamin = cellgeo->get_etamin();
              _etastep = cellgeo->get_etastep();
            }
          else if (_cell_binning == PHG4CylinderCellDefs::sizebinning)
            {
              _netabins = cellgeo->get_zbins();// bin eta in the same number of z bins
            }
          else if (_cell_binning == PHG4CylinderCellDefs::spacalbinning)
            {
              // use eta definiton for each row of towers
              _netabins = cellgeo->get_etabins();
            }
          else
            {
              cout <<"RawTowerCalibration::CreateNodes::"<<Name()
		   <<" - Fatal Error - unsupported cell binning method "<<_cell_binning<<endl;
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
             << " layer " << cellgeo->get_layer() << ": " << cellgeo->get_binning()
             << endl;
              exit(1);
            }
	  if (inner_radius > cellgeo->get_radius())
	    {
	      cout << "radius of layer " << cellgeo->get_layer() 
		   << " is " << cellgeo->get_radius()
		   << " which smaller than radius " << inner_radius
		   << " of first layer in list " << first_layer
		   << endl;
	    }
	  if (_nphibins != cellgeo->get_phibins())
	    {
	      cout << "mixing number of phibins, fisrt layer: " << _nphibins
		   << " layer " << cellgeo->get_layer() << ": " << cellgeo->get_phibins()
		   << endl;
	      exit(1);
	    }
	  if ( _phimin != cellgeo->get_phimin())
	    {
	      cout << "mixing number of phimin, fisrt layer: " << _phimin
		   << " layer " << cellgeo->get_layer() << ": " << cellgeo->get_phimin()
		   << endl;
	      exit(1);
	    }
	  if (_phistep != cellgeo->get_phistep())
	    {
	      cout << "mixing phisteps first layer: " << _phistep
		   << " layer " << cellgeo->get_layer() << ": " << cellgeo->get_phistep()
		   << " diff: " << _phistep - cellgeo->get_phistep()
		   << endl;
	      exit(1);
	    }
          if (_cell_binning == PHG4CylinderCellDefs::etaphibinning
              or _cell_binning == PHG4CylinderCellDefs::etaslatbinning)
            {
              if (_netabins != cellgeo->get_etabins())
                {
                  cout << "mixing number of_netabins , fisrt layer: " << _netabins
                 << " layer " << cellgeo->get_layer() << ": " << cellgeo->get_etabins()
                 << endl;
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

          else if (_cell_binning == PHG4CylinderCellDefs::sizebinning)
            {

              if (_netabins != cellgeo->get_zbins())
                {
                  cout << "mixing number of z bins , fisrt layer: " << _netabins
                 << " layer " << cellgeo->get_layer() << ": " << cellgeo->get_zbins()
                 << endl;
                  exit(1);
                }
            }
	}
    }
  rawtowergeom->set_radius(inner_radius);
  rawtowergeom->set_thickness(thickness);
  rawtowergeom->set_phibins(_nphibins);
  rawtowergeom->set_phistep(_phistep);
  rawtowergeom->set_phimin(_phimin);
  rawtowergeom->set_etabins(_netabins);

  if (first_cellgeo == NULL)
    {
      cout <<"RawTowerCalibration::CreateNodes - ERROR - can not find first layer of cells "<<endl;

      exit(1);
    }

  if (_cell_binning == PHG4CylinderCellDefs::etaphibinning
      or _cell_binning == PHG4CylinderCellDefs::etaslatbinning
      or _cell_binning == PHG4CylinderCellDefs::spacalbinning
      )
    {
//  rawtowergeom->set_etamin(_etamin);
//  rawtowergeom->set_etastep(_etastep);
      for (int ibin = 0; ibin<first_cellgeo->get_etabins(); ibin++)
        {

          const pair<double, double> range = first_cellgeo -> get_etabounds(ibin);

          rawtowergeom->set_etabounds(ibin, range);
        }

    }
  else if (_cell_binning == PHG4CylinderCellDefs::sizebinning)
    {
      for (int ibin = 0; ibin<first_cellgeo->get_zbins(); ibin++)
        {

          const pair<double, double> z_range = first_cellgeo -> get_zbounds(ibin);
          const double r = first_cellgeo->get_radius();
          const double eta1 = -log (tan(atan2(r, z_range.first)/2.) );
          const double eta2 = -log (tan(atan2(r, z_range.second)/2.) );
          rawtowergeom->set_etabounds(ibin, make_pair(eta1, eta2));
        }

    }
  else
    {
      cout <<"RawTowerCalibration::CreateNodes - ERROR - unsupported cell geometry "<<_cell_binning<<endl;
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
      throw std::runtime_error("Failed to find DST node in RawTowerCalibration::CreateNodes");
    }

  // Create the tower nodes on the tree
  _towers = new RawTowerContainer();
  CaliTowerNodeName = "TOWER_SIM_" + detector;
  PHIODataNode<PHObject> *towerNode = new PHIODataNode<PHObject>(_towers, CaliTowerNodeName.c_str(), "PHObject");
  dstNode->addNode(towerNode);

  return;
}




