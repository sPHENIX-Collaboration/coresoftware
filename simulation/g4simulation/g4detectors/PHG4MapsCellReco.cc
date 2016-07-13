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
  cout << "Creating PHG4MapsCellReco for name = " << name << endl;
  verbosity = 2;
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
      //cout << "---------- PHG4MapsCellReco:  Looping over layers " << endl;

      // loop over the hits in this layer
      PHG4HitContainer::ConstIterator hiter;
      PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits(*layer);

      // we need the geometry object for this layer
      PHG4CylinderGeom_MAPS *layergeom = dynamic_cast<PHG4CylinderGeom_MAPS *> (geo->GetLayerGeom(*layer));
      if(!layergeom)
	exit(1);

      layergeom->identify();

      for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
	{
	  cout << "From PHG4MapsCellReco: Call hit print method: " << endl;
	  hiter->second->print();

	  // get_property_int(const PROPERTY prop_id) const {return INT_MIN;}
	  int stave_number = hiter->second->get_property_int(PHG4Hit::prop_stave_index);
	  int half_stave_number = hiter->second->get_property_int(PHG4Hit::prop_half_stave_index);
	  int module_number = hiter->second->get_property_int(PHG4Hit::prop_module_index);
	  int chip_number = hiter->second->get_property_int(PHG4Hit::prop_chip_index);
	  double xsensor_in = hiter->second->get_property_float(PHG4Hit::prop_local_pos_x_0);
	  double ysensor_in = hiter->second->get_property_float(PHG4Hit::prop_local_pos_y_0);
	  double zsensor_in = hiter->second->get_property_float(PHG4Hit::prop_local_pos_z_0);
	  double xsensor_out = hiter->second->get_property_float(PHG4Hit::prop_local_pos_x_1);
	  double ysensor_out = hiter->second->get_property_float(PHG4Hit::prop_local_pos_y_1);
	  double zsensor_out = hiter->second->get_property_float(PHG4Hit::prop_local_pos_z_1);

	  TVector3 local_in(xsensor_in, ysensor_in, zsensor_in);
	  TVector3 local_out(xsensor_out, ysensor_out, zsensor_out);
      
	  // As a check, get the positions of the hit strips in world coordinates from the geo object 
	  TVector3 location_in = layergeom->get_world_from_local_coords(stave_number, half_stave_number, module_number, chip_number, local_in);
	  TVector3 location_out = layergeom->get_world_from_local_coords(stave_number, half_stave_number, module_number, chip_number, local_out);

	  cout << "      PHG4MapsCellReco:  Found entry location from geometry for  " 
	       << " stave number " << stave_number
	       << " half stave number " << half_stave_number 
	       << " module number" << module_number 
	       << endl
	       <<  " x = " << location_in.X()
	       << " y = " << location_in.Y()
	       << " z  = " << location_in.Z()
	       << " radius " << sqrt( pow(location_in.X(), 2) + pow(location_in.Y(), 2) ) 
	       << " angle " << atan( location_in.Y() / location_in.X() )
	       << endl;
	  cout << "     PHG4MapsCellReco: The entry location from G4 was "
	       << endl
	       << " x = " <<   hiter->second->get_x( 0)
	       << " y " <<  hiter->second->get_y( 0)
	       << " z " << hiter->second->get_z( 0)
	       << " radius " << sqrt( pow(hiter->second->get_x(0), 2) + pow(hiter->second->get_y(0), 2) )
	       << " angle " << atan( hiter->second->get_y(0) / hiter->second->get_x(0) )
	       << endl;
	  cout << " difference in x = " <<   hiter->second->get_x( 0) - location_in.X()  
	       << " difference in y = " <<   hiter->second->get_y( 0) - location_in.Y()  
	       << " difference in z = " <<   hiter->second->get_z( 0) - location_in.Z()  
	       << " difference in radius = " <<  sqrt( pow(hiter->second->get_x(0), 2) + pow(hiter->second->get_y(0), 2) )  - sqrt( pow(location_in.X(),2) + pow(location_in.Y(),2) )    
	       << " in angle = " <<  atan( hiter->second->get_y(0) / hiter->second->get_x(0) )  -  atan( location_in.Y() / location_in.X() )  
	       << endl << endl;

	  cout << "      PHG4MapsCellReco:  Found exit location from geometry for  " 
	       << " stave number " << stave_number
	       << " half stave number " << half_stave_number 
	       << " module number" << module_number 
	       << endl
	       <<  " x = " << location_out.X()
	       << " y = " << location_out.Y()
	       << " z  = " << location_out.Z()
	       << " radius " << sqrt( pow(location_out.X(), 2) + pow(location_out.Y(), 2) ) 
	       << " angle " << atan( location_out.Y() / location_out.X() )
	       << endl;
	  cout << "     PHG4MapsCellReco: The exit location from G4 was "
	       << endl
	       << " x = " <<   hiter->second->get_x( 1)
	       << " y " <<  hiter->second->get_y( 1)
	       << " z " << hiter->second->get_z( 1)
	       << " radius " << sqrt( pow(hiter->second->get_x(1), 2) + pow(hiter->second->get_y(1), 2) )
	       << " angle " << atan( hiter->second->get_y(1) / hiter->second->get_x(1) )
	       << endl;
	  cout << " difference in radius = " <<   sqrt( pow(hiter->second->get_x(1), 2) + pow(hiter->second->get_y(1), 2) )  - sqrt( pow(location_out.X(),2) + pow(location_out.Y(),2) )  
	       << " in angle = " <<  atan( hiter->second->get_y(1) / hiter->second->get_x(1) )  - atan( location_out.Y() / location_out.X() ) 
	       << endl << endl;

	  // Get the pixel number of the input hit
	  int pixel_number = layergeom->get_pixel_from_local_coords(local_in);
	  cout << " CellReco: pixel number = " << pixel_number << endl;
	  
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

	  //if (verbosity > 0)
	  //{
	      cout << "From MapsCellReco: Adding cell for stave: " << mapiter->second->get_stave_index()
		   << " , half stave index: " << mapiter->second->get_half_stave_index()
		   << ", module index: " << mapiter->second->get_module_index()
		   << ", chip index: " << mapiter->second->get_chip_index()
		   << ", energy dep: " << mapiter->second->get_edep()
		   << endl;
	      //}
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

