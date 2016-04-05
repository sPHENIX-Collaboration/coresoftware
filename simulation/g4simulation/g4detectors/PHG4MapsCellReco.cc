#include "PHG4MapsCellReco.h"
#include "PHG4CylinderGeomContainer.h"
#include "PHG4CylinderGeom_MAPS.h"
#include "PHG4CylinderCellGeom_MAPS.h"
#include "PHG4CylinderCellGeomContainer.h"
#include "PHG4CylinderCell_MAPS.h"
#include "PHG4CylinderCellContainer.h"
//#include "PHG4CylinderCellDefs.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4HitContainer.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>


#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>

using namespace std;

PHG4MapsCellReco::PHG4MapsCellReco(const string &name) :
  SubsysReco(name),
  _timer(PHTimeServer::get()->insert_new(name.c_str())),
  chkenergyconservation(0)
{
  memset(nbins, 0, sizeof(nbins));
  Detector(name);
}

int PHG4MapsCellReco::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
    {
      std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
      exit(1);
    }
  hitnodename = "G4HIT_" + detector;
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  if (!g4hit)
    {
      cout << "Could not locate g4 hit node " << hitnodename << endl;
      exit(1);
    }
  cellnodename = "G4CELL_" + detector;
  PHG4CylinderCellContainer *cells = findNode::getClass<PHG4CylinderCellContainer>(topNode , cellnodename);
  if (!cells)
    {
      PHNodeIterator dstiter(dstNode);
      PHCompositeNode *DetNode =
          dynamic_cast<PHCompositeNode*>(dstiter.findFirst("PHCompositeNode",
              detector));
      if (!DetNode)
        {
          DetNode = new PHCompositeNode(detector);
          dstNode->addNode(DetNode);
        }
      cells = new PHG4CylinderCellContainer();
      PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(cells, cellnodename.c_str() , "PHObject");
      DetNode->addNode(newNode);
    
    }

  geonodename = "CYLINDERGEOM_" + detector;
  PHG4CylinderGeomContainer *geo =  findNode::getClass<PHG4CylinderGeomContainer>(topNode , geonodename.c_str());
  if (!geo)
    {
      cout << "Could not locate geometry node " << geonodename << endl;
      exit(1);
    }
  if (verbosity > 0)
    {
      geo->identify();
    }

  return Fun4AllReturnCodes::EVENT_OK;
}


int
PHG4MapsCellReco::process_event(PHCompositeNode *topNode)
{
  _timer.get()->restart();
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  if (!g4hit)
    {
      cout << "Could not locate g4 hit node " << hitnodename << endl;
      exit(1);
    }
  PHG4CylinderCellContainer *cells = findNode::getClass<PHG4CylinderCellContainer>(topNode, cellnodename);
  if (! cells)
    {
      cout << "could not locate cell node " << cellnodename << endl;
      exit(1);
    }

  geonodename = "CYLINDERGEOM_" + detector;
  PHG4CylinderGeomContainer *geo =  findNode::getClass<PHG4CylinderGeomContainer>(topNode , geonodename.c_str());
  if (!geo)
    {
      cout << "Could not locate geometry node " << geonodename << endl;
      exit(1);
    }

  cells->Reset();

  // loop over all of the layers in the hit container
  PHG4HitContainer::LayerIter layer;
  pair<PHG4HitContainer::LayerIter, PHG4HitContainer::LayerIter> layer_begin_end = g4hit->getLayers();
  for (layer = layer_begin_end.first; layer != layer_begin_end.second; ++layer)
    {
      //cout << "---------- MapsReco:  Looping over layers " << endl;

      // loop over the hits in this layer
      PHG4HitContainer::ConstIterator hiter;
      PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits(*layer);

      // we need the geometry object for this layer
      PHG4CylinderGeom_MAPS *layergeom = dynamic_cast<PHG4CylinderGeom_MAPS *> (geo->GetLayerGeom(*layer));
      if(!layergeom)
	exit(1);

      //layergeom->identify();

      for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
	{
	  // get_property_int(const PROPERTY prop_id) const {return INT_MIN;}
	  int stave_number = hiter->second->get_property_int(PHG4Hit::prop_stave_index);
	  int half_stave_number = hiter->second->get_property_int(PHG4Hit::prop_half_stave_index);
	  int module_number = hiter->second->get_property_int(PHG4Hit::prop_module_index);
	  int chip_number = hiter->second->get_property_int(PHG4Hit::prop_chip_index);

	  // As a check, get the positions of the hit strips from the geo object
	  double location[3] = {-1,-1,-1};
	  layergeom->find_sensor_center(stave_number, half_stave_number, module_number, chip_number, location);

	  if(verbosity > 0)
	    {
	      cout << "      Found location from geometry for  " 
		   << " stave number " << stave_number
		   << " half stave number " << half_stave_number 
		   << " module number" << module_number
		   <<  " is:  x = " << location[0]
		   << " y = " << location[1]
		   << " z  = " << location[2]
		   << endl;
	    }

	  // combine ladder index values to get a single key
	  char inkey[1024];
	  sprintf(inkey,"%i-%i_%i_%i",stave_number, half_stave_number, module_number, chip_number);
	  std::string key(inkey);

	  if (celllist.count(key) > 0) {
	    celllist[key]->add_edep(hiter->first, hiter->second->get_edep());
	  } else {
	    celllist[key] = new PHG4CylinderCell_MAPS();
	    celllist[key]->set_layer(*layer);

	    // This encodes the z and phi position of the sensor 
	    celllist[key]->set_stave_index(stave_number);
	    celllist[key]->set_half_stave_index(half_stave_number);
	    celllist[key]->set_module_index(module_number);
	    celllist[key]->set_chip_index(chip_number);
	    
	    celllist[key]->add_edep(hiter->first, hiter->second->get_edep());
	  }
	} // end loop over g4hits

      int numcells = 0;
      for (map<std::string, PHG4CylinderCell_MAPS *>::const_iterator mapiter = celllist.begin();mapiter != celllist.end() ; ++mapiter)
	{	  
	  cells->AddCylinderCell(*layer, mapiter->second);
	  numcells++;

	  if (verbosity > 0)
	    {
	      cout << "Adding cell for sensor_index: " << mapiter->second->get_sensor_index()
		   << ", srip z index: " << mapiter->second->get_binz()
		   << ", strip y index: " << mapiter->second->get_binphi()
		   << ", energy dep: " << mapiter->second->get_edep()
		   << endl;
	    }
	}
      celllist.clear();
      if (verbosity > 0)
	{
	  cout << Name() << ": found " << numcells << " silicon strips with energy deposition" << endl;
	}
    }


if (chkenergyconservation)
  {
    CheckEnergy(topNode);
  }
_timer.get()->stop();
return Fun4AllReturnCodes::EVENT_OK;
}

int
PHG4MapsCellReco::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int
PHG4MapsCellReco::CheckEnergy(PHCompositeNode *topNode)
{
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  PHG4CylinderCellContainer *cells = findNode::getClass<PHG4CylinderCellContainer>(topNode, cellnodename);
  double sum_energy_g4hit = 0.;
  double sum_energy_cells = 0.;
  PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits();
  PHG4HitContainer::ConstIterator hiter;
  for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
    {
      sum_energy_g4hit += hiter->second->get_edep();
    }
  PHG4CylinderCellContainer::ConstRange cell_begin_end = cells->getCylinderCells();
  PHG4CylinderCellContainer::ConstIterator citer;
  for (citer = cell_begin_end.first; citer != cell_begin_end.second; ++citer)
    {
      sum_energy_cells += citer->second->get_edep();

    }
  // the fractional eloss for particles traversing eta bins leads to minute rounding errors
  if (fabs(sum_energy_cells - sum_energy_g4hit) / sum_energy_g4hit > 1e-6)
    {
      cout << "energy mismatch between cells: " << sum_energy_cells
	   << " and hits: " << sum_energy_g4hit
	   << " diff sum(cells) - sum(hits): " << sum_energy_cells - sum_energy_g4hit
	   << endl;
      return -1;
    }
  else
    {
      if (verbosity > 0)
	{
	  cout << Name() << ": total energy for this event: " << sum_energy_g4hit << " GeV" << endl;
	}
    }
  return 0;
}

