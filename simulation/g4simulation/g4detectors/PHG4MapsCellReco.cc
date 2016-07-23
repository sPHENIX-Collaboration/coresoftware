#include "PHG4MapsCellReco.h"
#include "PHG4CylinderGeomContainer.h"
#include "PHG4CylinderGeom_MAPS.h"
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

  if(verbosity > 0)  
    cout << "Creating PHG4MapsCellReco for name = " << name << endl;
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

      if(verbosity > 2)
	layergeom->identify();

      for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
	{
	  //cout << "From PHG4MapsCellReco: Call hit print method: " << endl;
	  //hiter->second->print();

	  // get_property_int(const PROPERTY prop_id) const {return INT_MIN;}
	  int stave_number = hiter->second->get_property_int(PHG4Hit::prop_stave_index);
	  int half_stave_number = hiter->second->get_property_int(PHG4Hit::prop_half_stave_index);
	  int module_number = hiter->second->get_property_int(PHG4Hit::prop_module_index);
	  int chip_number = hiter->second->get_property_int(PHG4Hit::prop_chip_index);

	  TVector3 local_in( hiter->second->get_local_x(0),  hiter->second->get_local_y(0),  hiter->second->get_local_z(0) );
 	  TVector3 local_out( hiter->second->get_local_x(1),  hiter->second->get_local_y(1),  hiter->second->get_local_z(1) );

	  /*
	  // The G4 transform to local coordinates for the exit point does not work properly, so the local exit coordinates in the hit object are wrong - now fixed!!
	  // Use the method in the geometry object to convert from world to local coords
	  TVector3 world_out(  hiter->second->get_x( 1),  hiter->second->get_y( 1),  hiter->second->get_z( 1) );
	  TVector3 local_out =  layergeom->get_local_from_world_coords(stave_number, half_stave_number, module_number, chip_number, world_out);
	  */

	  if(verbosity > 4)
	    {
	      cout << endl << "  world entry point position: " << hiter->second->get_x(0) << " " << hiter->second->get_y(0) << " " << hiter->second->get_z(0) << endl;
	      cout << "  world exit point position: " << hiter->second->get_x(1) << " " << hiter->second->get_y(1) << " " << hiter->second->get_z(1) << endl;
	      cout << "  local coords of entry point from G4 " << hiter->second->get_local_x(0)  << " " << hiter->second->get_local_y(0) << " " << hiter->second->get_local_z(0) << endl;      
	      cout << "  local coords of exit point from geom " << local_out.X()  << " " << local_out.Y() << " " << local_out.Z() << endl;      
	      cout << "  local coords of exit point from G4 " << hiter->second->get_local_x(1)  << " " << hiter->second->get_local_y(1) << " " << hiter->second->get_local_z(1) << endl;      
	      TVector3 world_in(  hiter->second->get_x( 0),  hiter->second->get_y( 0),  hiter->second->get_z( 0) );
	      TVector3 local_in_check =  layergeom->get_local_from_world_coords(stave_number, half_stave_number, module_number, chip_number, world_in);
	      cout << "  local coords of entry point from geom (check) " << local_in_check.X()  << " " << local_in_check.Y() << " " << local_in_check.Z() << endl;      	      
	      cout << endl;
	    }

	  if(verbosity > 2)
	    {
	      // As a check, get the positions of the hit strips in world coordinates from the geo object 
	      TVector3 location_in = layergeom->get_world_from_local_coords(stave_number, half_stave_number, module_number, chip_number, local_in);
	      TVector3 location_out = layergeom->get_world_from_local_coords(stave_number, half_stave_number, module_number, chip_number, local_out);
	      
	      cout << endl << "      PHG4MapsCellReco:  Found world entry location from geometry for  " 
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
	      cout << "     PHG4MapsCellReco: The world entry location from G4 was "
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
	      
	      cout << "      PHG4MapsCellReco:  Found world exit location from geometry for  " 
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
	      cout << "     PHG4MapsCellReco: The world exit location from G4 was "
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
	    }

	  /*
	  // Get the pixel number of the entry location
	  int pixel_number_in = layergeom->get_pixel_from_local_coords(pixel_x, pixel_y, local_in);
	  // Get the pixel number of the exit location
	  int pixel_number_out = layergeom->get_pixel_from_local_coords(pixel_x, pixel_y, local_out);

	  cout << "entry pixel number " << pixel_number_in << " exit pixel number " << pixel_number_out << endl;

	  // Are they different?
	  int number_of_pixels = 0;
	  if(pixel_number_out != pixel_number_in)
	    {
	      // There is more than one pixel, have to divide the energy between them?
	      // Get the list of hit pixels

	      // Get the X and Y grid locations of the pixels?
	      int Ngridx_in = layergeom->get_pixel_X_from_pixel_number(pixel_x, pixel_y, pixel_number_in);
	      int Ngridx_out = layergeom->get_pixel_X_from_pixel_number(pixel_x, pixel_y, pixel_number_out);

	      int Ngridy_in = layergeom->get_pixel_Y_from_pixel_number(pixel_x, pixel_y, pixel_number_in);
	      int Ngridy_out = layergeom->get_pixel_Y_from_pixel_number(pixel_x, pixel_y, pixel_number_out);

	      number_of_pixels = abs(Ngridx_in - Ngridx_out) + abs(Ngridy_in - Ngridy_out) + 1;

	      cout << " Ngridx_in " << Ngridx_in
		   << " Ngridy_in " << Ngridy_in
		   << " Ngridx_out " << Ngridx_out
		   << " Ngridy_out " << Ngridy_out
		   << endl; 
	    }
	  else
	    {
	      //There is only one pixel, it gets all of the energy
	      number_of_pixels = 1;

	      // Get the X and Y grid locations of the pixels?
	      int Ngridx_in = layergeom->get_pixel_X_from_pixel_number(pixel_x, pixel_y, pixel_number_in);
	      int Ngridx_out = layergeom->get_pixel_X_from_pixel_number(pixel_x, pixel_y, pixel_number_out);

	      int Ngridy_in = layergeom->get_pixel_Y_from_pixel_number(pixel_x, pixel_y, pixel_number_in);
	      int Ngridy_out = layergeom->get_pixel_Y_from_pixel_number(pixel_x, pixel_y, pixel_number_out);

	      cout << " Ngridx_in " << Ngridx_in
		   << " Ngridy_in " << Ngridy_in
		   << " Ngridx_out " << Ngridx_out
		   << " Ngridy_out " << Ngridy_out
		   << endl; 
	    }
	  cout << "number of pixels  " << number_of_pixels << endl;

	  if(verbosity > 2)
	    {
	      // Testing: get the local coords of the center of the pixel
	      TVector3  pixel_local_coords_in = layergeom->get_local_coords_from_pixel(pixel_x, pixel_y, pixel_number_in);
	      cout << " CellReco: pixel number in = " << pixel_number_in << endl;
	      // check - get pixel number from coords of center of pixel!
	      int pixel_number_check =  layergeom->get_pixel_from_local_coords(pixel_x, pixel_y, local_in);
	      cout << "pixel number in check = " << pixel_number_check << endl;
	      cout << " PHG4MapsCellReco: pixel local coords from pixel number in = " << pixel_local_coords_in.X() << " " << pixel_local_coords_in.Y() << " " << pixel_local_coords_in.Z() << endl;
	      cout << " PHG4MapsCellReco:hit entry point  local coords from G4 = " << local_in.X() << " " << local_in.Y() << " " << local_in.Z() << endl;
	      cout << " PHG4MapsCellReco:hit exit point  local coords from Geom = " << local_out.X() << " " << local_out.Y() << " " << local_out.Z() << endl;
	    }

	  */

	  // This is a temporary make-do to get a single pixel that represents the hit.
	  // Find the mid point between the entry and exit locations in the sensor

	  TVector3 midpoint( (local_in.X() + local_out.X()) / 2.0, (local_in.Y() + local_out.Y()) / 2.0, (local_in.Z() + local_out.Z()) / 2.0 );
	  int pixel_number_mid = layergeom->get_pixel_from_local_coords(midpoint);
	  TVector3  pixel_local_coords_mid = layergeom->get_local_coords_from_pixel(pixel_number_mid);
	  TVector3 pixel_world_coords_mid = layergeom->get_world_from_local_coords(stave_number, half_stave_number, module_number, chip_number, pixel_local_coords_mid);

	  if(verbosity > 4)
	    cout << " pixel number " << pixel_number_mid 
		 << " pixel local coords " << pixel_local_coords_mid.X() << " " << pixel_local_coords_mid.Y() << " " << pixel_local_coords_mid.Z() 
		 << " pixel world coords " << pixel_world_coords_mid.X() << " " << pixel_world_coords_mid.Y() << " " << pixel_world_coords_mid.Z() 
		 << endl << endl;
	  
	  // combine ladder index values to get a single key
	  char inkey[1024];
	  sprintf(inkey,"%i-%i_%i_%i",stave_number, half_stave_number, module_number, chip_number);
	  std::string key(inkey);

	  if (celllist.count(key) > 0) {
	    celllist[key]->add_edep(hiter->first, hiter->second->get_edep());
	  } else {
	    celllist[key] = new PHG4CylinderCell_MAPS();
	    celllist[key]->set_layer(*layer);

	    // This encodes the z and phi position of the sensor, and pixel number within the sensor
	    celllist[key]->set_stave_index(stave_number);
	    celllist[key]->set_half_stave_index(half_stave_number);
	    celllist[key]->set_module_index(module_number);
	    celllist[key]->set_chip_index(chip_number);
	    celllist[key]->set_pixel_index(pixel_number_mid);
	    
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
	      cout << "From MapsCellReco: Adding cell for stave: " << mapiter->second->get_stave_index()
		   << " , half stave index: " << mapiter->second->get_half_stave_index()
		   << ", module index: " << mapiter->second->get_module_index()
		   << ", chip index: " << mapiter->second->get_chip_index()
		   << ", pixel index: " << mapiter->second->get_pixel_index()
		   << ", energy dep: " << mapiter->second->get_edep()
		   << endl;
	    }
	}
      celllist.clear();
      if (verbosity > 0)
	{
	  cout << Name() << ": found " << numcells << " silicon pixels with energy deposition" << endl;
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

