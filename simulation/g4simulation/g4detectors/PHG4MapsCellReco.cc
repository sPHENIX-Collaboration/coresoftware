#include "PHG4MapsCellReco.h"
#include "PHG4CylinderGeomContainer.h"
#include "PHG4CylinderGeom_MAPS.h"
#include "PHG4CylinderCell_MAPS.h"
#include "PHG4CylinderCellContainer.h"

#include <g4main/PHG4Hit.h>
#include <g4main/PHG4Hitv1.h>
#include <g4main/PHG4HitContainer.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>

#include <TSystem.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>

using namespace std;

PHG4MapsCellReco::PHG4MapsCellReco(const string &name) :
  SubsysReco(name),
  detector(name),
  _timer(PHTimeServer::get()->insert_new(name)),
  chkenergyconservation(0)
{
  memset(nbins, 0, sizeof(nbins));

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
	  if(verbosity >4)
	    hiter->second->print();

	  // get_property_int(const PROPERTY prop_id) const {return INT_MIN;}
	  int stave_number = hiter->second->get_property_int(PHG4Hit::prop_stave_index);
	  int half_stave_number = hiter->second->get_property_int(PHG4Hit::prop_half_stave_index);
	  int module_number = hiter->second->get_property_int(PHG4Hit::prop_module_index);
	  int chip_number = hiter->second->get_property_int(PHG4Hit::prop_chip_index);

	  TVector3 local_in( hiter->second->get_local_x(0),  hiter->second->get_local_y(0),  hiter->second->get_local_z(0) );
 	  TVector3 local_out( hiter->second->get_local_x(1),  hiter->second->get_local_y(1),  hiter->second->get_local_z(1) );

	  if(verbosity > 4)
	    {
	      cout << endl << "  world entry point position: " << hiter->second->get_x(0) << " " << hiter->second->get_y(0) << " " << hiter->second->get_z(0) << endl;
	      cout << "  world exit point position: " << hiter->second->get_x(1) << " " << hiter->second->get_y(1) << " " << hiter->second->get_z(1) << endl;
	      cout << "  local coords of entry point from G4 " << hiter->second->get_local_x(0)  << " " << hiter->second->get_local_y(0) << " " << hiter->second->get_local_z(0) << endl;   
	      TVector3 world_in(  hiter->second->get_x( 0),  hiter->second->get_y( 0),  hiter->second->get_z( 0) );
	      TVector3 local_in_check =  layergeom->get_local_from_world_coords(stave_number, half_stave_number, module_number, chip_number, world_in);
	      cout << "  local coords of entry point from geom (a check) " << local_in_check.X()  << " " << local_in_check.Y() << " " << local_in_check.Z() << endl;      	         
	      cout << "  local coords of exit point from G4 " << hiter->second->get_local_x(1)  << " " << hiter->second->get_local_y(1) << " " << hiter->second->get_local_z(1) << endl;      
	      cout << "  local coords of exit point from geom (a check) " << local_out.X()  << " " << local_out.Y() << " " << local_out.Z() << endl;      
	      cout << endl;
	    }

	  if(verbosity > 2)
	    {
	      // As a check, get the positions of the hit pixels in world coordinates from the geo object 
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


	  // Get the pixel number of the entry location
	  int pixel_number_in = layergeom->get_pixel_from_local_coords(local_in);
	  // Get the pixel number of the exit location
	  int pixel_number_out = layergeom->get_pixel_from_local_coords(local_out);

	  if(verbosity > 0)
	    cout << "entry pixel number " << pixel_number_in << " exit pixel number " << pixel_number_out << endl;

	  vector<int> vpixel;
	  vector<int> vxbin;
	  vector<int> vzbin;
	  vector<double> vlen;
	  double trklen = 0.0;

	  // Are they different?
	  bool test_one_pixel = false;  // normally false!

	  if(pixel_number_out != pixel_number_in)
	    {
	      if(test_one_pixel)
		{
		  // For testing, assign the hit to 1 pixel, the one encompassing the center of the track line through the sensor
		  // All of the energy will assigned to this pixel

		  TVector3 pathvec = local_in - local_out;
		  trklen = sqrt( pow( local_in.X() - local_out.X(), 2 ) + pow( local_in.Z() - local_out.Z(), 2) );     // only the length in x and z plane

		  // Find the mid point between the entry and exit locations in the sensor
		  TVector3 midpoint( (local_in.X() + local_out.X()) / 2.0, (local_in.Y() + local_out.Y()) / 2.0, (local_in.Z() + local_out.Z()) / 2.0 );
		  int pixel_number_mid = layergeom->get_pixel_from_local_coords(midpoint);
		  // get the phi and Z index for this pixel, it is needed later by the clustering
		  int xbin = layergeom->get_pixel_X_from_pixel_number(pixel_number_mid);
		  int zbin = layergeom->get_pixel_Z_from_pixel_number(pixel_number_mid);
		  
		  vpixel.push_back(pixel_number_mid);
		  vxbin.push_back(xbin);
		  vzbin.push_back(zbin);
		  vlen.push_back(trklen);

		  if(verbosity > 0)
		    cout << " Test one pixel: pixel number mid " << pixel_number_mid 
			 << " pixel_number_in " << pixel_number_in
			 << " pixel_number_out " << pixel_number_out
			 << " xbin " << xbin << " zbin " << zbin << " trklen " << trklen 
			 << " edep " <<  hiter->second->get_edep()
			 << endl << endl; 
		}
	      else
		{
		  // This is the correct way
		  // There is more than one pixel, have to divide the energy between them
		  
		  // Get the X and Y bin locations of the pixels
		  int xbin_in = layergeom->get_pixel_X_from_pixel_number(pixel_number_in);
		  int zbin_in = layergeom->get_pixel_Z_from_pixel_number(pixel_number_in);
		  int xbin_out = layergeom->get_pixel_X_from_pixel_number(pixel_number_out);
		  int zbin_out = layergeom->get_pixel_Z_from_pixel_number(pixel_number_out);
		  
		  // make sure that the entry and exit bins are in increasing order, needed by line_and_rectangle_intersect
		  if(xbin_in > xbin_out)
		    {
		      int tmp = xbin_out;
		      xbin_out = xbin_in;
		      xbin_in = tmp;
		    }
		  if(zbin_in > zbin_out)
		    {
		      int tmp = zbin_out;
		      zbin_out = zbin_in;
		      zbin_in = tmp;
		    }
		  
		  // get the line connecting the entry and exit points
		  double ax = local_in.X();
		  double az = local_in.Z();
		  double bx = local_out.X();
		  double bz = local_out.Z();
		  
		  TVector3 pathvec = local_in - local_out;
		  trklen = sqrt( pow( local_in.X() - local_out.X(), 2 ) + pow( local_in.Z() - local_out.Z(), 2) );     // only the length in x and z plane
		  
		  for(int xbin=xbin_in; xbin<=xbin_out;xbin++)
		    for(int zbin=zbin_in;zbin<=zbin_out;zbin++)
		      {
			// get the pixel center
			int pixnum = layergeom->get_pixel_number_from_xbin_zbin( xbin, zbin );
			TVector3 tmp = layergeom->get_local_coords_from_pixel( pixnum );
			// note that we need cx < dx and cz < dz
			double cx = tmp.X() - layergeom->get_pixel_x() / 2.0;
			double dx = tmp.X() + layergeom->get_pixel_x() / 2.0;
			double cz = tmp.Z() - layergeom->get_pixel_z() / 2.0;
			double dz = tmp.Z() + layergeom->get_pixel_z() / 2.0;
			double rr = 0.0;

			bool yesno = line_and_rectangle_intersect(ax, az, bx, bz, cx, cz, dx, dz, &rr);
			
			if (yesno)
			  {
			    if (verbosity > 0) cout << "CELL FIRED: "  << " xbin " << xbin << " zbin  " << zbin << " rr " << rr  << endl;
			    if(verbosity>2)
			      {
				cout << "##### line: ax " << ax << " az " << az << " bx " << bx << " bz " << bz << endl;
				cout << "####### cell:  cx " << cx << " cz " << cz << " dx " << dx << " dz " << dz << endl;
			      }

			    vpixel.push_back(pixnum);
			    vxbin.push_back(xbin);
			    vzbin.push_back(zbin);
			    vlen.push_back(rr);
			  }
		      }
		}
	    }
	  else
	    {
	      //There is only one pixel, it gets all of the energy
	      int xbin = layergeom->get_pixel_X_from_pixel_number(pixel_number_in);
	      int zbin = layergeom->get_pixel_Z_from_pixel_number(pixel_number_in);
	 
	      vpixel.push_back(pixel_number_in);
	      vxbin.push_back(xbin);
	      vzbin.push_back(zbin);
	      vlen.push_back(trklen);
	    }

	  // loop over all fired cells for this hit and add them to the celllist
	  for (unsigned int i1 = 0; i1 < vpixel.size(); i1++)   // loop over all fired cells
	    {
	      int pixel_number = vpixel[i1];

	      // combine ladder index and pixel values to get a single unique key for this pixel
	      // layers:     0 - 2
	      // Stave index:   0 - 47   = 6 bits
	      // Half stave index:  0 - 1 2 bits
	      // Module index: 0 - 6 3 bits
	      // Chip index:  0 - 13
	      // Pixel index:   0 - 1.14E+06  // yes, that is 1.14 million
	      // check validity (if values are within their assigned number of bits)
	      unsigned long long tmp = pixel_number;
	      unsigned long long inkey = tmp << 32;
	      static unsigned int stave_number_bits = 0x8;
	      static unsigned int stave_number_max = pow(2,stave_number_bits);
	      static unsigned int half_stave_number_bits = 0x2;
	      static unsigned int half_stave_number_max = pow(2,half_stave_number_bits);
	      static unsigned int module_number_bits = 0x2;
	      static unsigned int module_number_max = pow(2,module_number_bits);
	      static unsigned int chip_number_bits = 0x4;
	      static unsigned int chip_number_max = pow(2,chip_number_bits);

	      if (static_cast<unsigned int> (stave_number) > stave_number_max)
		{
		  cout << "stave number " << stave_number << " exceeds valid value " << stave_number_max << endl;
		  gSystem->Exit(1);
		  exit(1); // make coverity happy which does not know about gSystem->Exit()
		}
	      if (static_cast<unsigned int> (half_stave_number) > half_stave_number_max)
		{
		  cout << "half stave number " << half_stave_number << " exceeds valid value " << half_stave_number_max << endl;
		  gSystem->Exit(1);
		  exit(1); // make coverity happy which does not know about gSystem->Exit()
		}
	      if (static_cast<unsigned int> (module_number) > module_number_max)
		{
		  cout << "module_number " << module_number << " exceeds valid value " << module_number_max << endl;
		  gSystem->Exit(1);
		  exit(1); // make coverity happy which does not know about gSystem->Exit()
		}
	      if (static_cast<unsigned int> (chip_number) > chip_number_max)
		{
		  cout << "chip_number " << chip_number << " exceeds valid value " << chip_number_max << endl;
		  gSystem->Exit(1);
		  exit(1); // make coverity happy which does not know about gSystem->Exit()
		}
	      inkey += stave_number;
	      inkey += (half_stave_number << stave_number_bits);
	      inkey += (module_number << (stave_number_bits+half_stave_number_bits));
	      inkey += (chip_number << (stave_number_bits+half_stave_number_bits+module_number_bits));
	      PHG4CylinderCell *cell = nullptr;
	      map<unsigned long long, PHG4CylinderCell*>::iterator it;
	      it = celllist.find(inkey);
	      if (it != celllist.end())
		{
		  cell = it->second;
		}
	      else
		{
		  cell = new PHG4CylinderCell_MAPS();
		  celllist[inkey] = cell;
		  cell->set_layer(*layer);
		  cell->set_stave_index(stave_number);
		  cell->set_half_stave_index(half_stave_number);
		  cell->set_module_index(module_number);
		  cell->set_chip_index(chip_number);
		  cell->set_pixel_index(pixel_number);
		  cell->set_phibin(vxbin[i1]);
		  cell->set_zbin(vzbin[i1]);
		}
	      double edep;
	      if(trklen > 0.0)
		{
		  edep = hiter->second->get_edep() * vlen[i1] / trklen;
		}
	      else
		{
		  edep = hiter->second->get_edep();
		}
              cell->add_edep(hiter->first, edep);

		  
	      if(verbosity > 1)
		{
		  cout << " looping over fired cells: cell " << i1 << " inkey 0x" << hex << inkey << dec 
		       << " cell length " << vlen[i1] << " trklen " << trklen
		       << " cell edep " << edep << " total edep " << hiter->second->get_edep() << endl;
	
		}
	    }
	} // end loop over g4hits
      
      int numcells = 0;
      for (map<unsigned long long, PHG4CylinderCell *>::const_iterator mapiter = celllist.begin();mapiter != celllist.end() ; ++mapiter)
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

bool PHG4MapsCellReco::lines_intersect(
				       double ax,
				       double ay,
				       double bx,
				       double by,
				       double cx,
				       double cy,
				       double dx,
				       double dy,
				       double* rx, // intersection point (output)
				       double* ry
				       )
{
  
  // Find if a line segment limited by points A and B
  // intersects line segment limited by points C and D.
  // First check if an infinite line defined by A and B intersects
  // segment (C,D). If h is from 0 to 1 line and line segment intersect
  // Then check in intersection point is between C and D
  
  double ex = bx - ax; // E=B-A
  double ey = by - ay;
  double fx = dx - cx; // F=D-C
  double fy = dy - cy;
  double px = -ey;     // P
  double py = ex;

  double bottom = fx * px + fy * py; // F*P
  double gx = ax - cx; // A-C
  double gy = ay - cy;
  double top = gx * px + gy * py; // G*P

  double h = 99999.;
  if (bottom != 0.)
    {
      h = top / bottom;
    }

//intersection point R = C + F*h
  if (h > 0. && h < 1.)
    {
      *rx = cx + fx * h;
      *ry = cy + fy * h;
      //cout << "      line/segment intersection coordinates: " << *rx << " " << *ry << endl;
      if ((*rx > ax && *rx > bx) || (*rx < ax && *rx < bx) || (*ry < ay && *ry < by) || (*ry > ay && *ry > by))
        {
          //cout << "       NO segment/segment intersection!" << endl;
          return false;
        }
      else
        {
          //cout << "       segment/segment intersection!" << endl;
          return true;
        }
    }

  return false;
}

bool  PHG4MapsCellReco::line_and_rectangle_intersect(
						     double ax,
						     double ay,
						     double bx,
						     double by,
						     double cx,
						     double cy,
						     double dx,
						     double dy,
						     double* rr // length of the line segment inside the rectangle (output)
						     )
{
  
  // find if a line isegment limited by points (A,B)
  // intersects with a rectangle defined by two
  // corner points (C,D) two other points are E and F
  //   E--------D
  //   |        |
  //   |        |
  //   C--------F
  
  if (cx > dx || cy > dy)
    {
      cerr << "ERROR: Bad rectangle definition!" << endl;
      return false;
    }

  double ex = cx;
  double ey = dy;
  double fx = dx;
  double fy = cy;
  double rx = 99999.;
  double ry = 99999.;

  vector<double> vx;
  vector<double> vy;

  bool i1 = lines_intersect(ax, ay, bx, by, cx, cy, fx, fy, &rx, &ry);
  if (i1)
    {
      vx.push_back(rx);
      vy.push_back(ry);
    }
  bool i2 = lines_intersect(ax, ay, bx, by, fx, fy, dx, dy, &rx, &ry);
  if (i2)
    {
      vx.push_back(rx);
      vy.push_back(ry);
    }
  bool i3 = lines_intersect(ax, ay, bx, by, ex, ey, dx, dy, &rx, &ry);
  if (i3)
    {
      vx.push_back(rx);
      vy.push_back(ry);
    }
  bool i4 = lines_intersect(ax, ay, bx, by, cx, cy, ex, ey, &rx, &ry);
  if (i4)
    {
      vx.push_back(rx);
      vy.push_back(ry);
    }

//cout << "Rectangle intersections: " << i1 << " " << i2 << " " << i3 << " " << i4 << endl;
//cout << "Number of intersections = " << vx.size() << endl;

  *rr = 0.;
  if (vx.size() == 2)
    {
      *rr = sqrt( (vx[0] - vx[1]) * (vx[0] - vx[1]) + (vy[0] - vy[1]) * (vy[0] - vy[1]) );
//  cout << "Length of intersection = " << *rr << endl;
    }
  if (vx.size() == 1)
    {
      // find which point (A or B) is within the rectangle
      if (ax > cx && ay > cy && ax < dx && ay < dy)   // point A is inside the rectangle
        {
          //cout << "Point A is inside the rectangle." << endl;
          *rr = sqrt((vx[0] - ax) * (vx[0] - ax) + (vy[0] - ay) * (vy[0] - ay));
        }
      if (bx > cx && by > cy && bx < dx && by < dy)   // point B is inside the rectangle
        {
          //cout << "Point B is inside the rectangle." << endl;
          *rr = sqrt((vx[0] - bx) * (vx[0] - bx) + (vy[0] - by) * (vy[0] - by));
        }
    }

  if (i1 || i2 || i3 || i4)
    {
      return true;
    }
  return false;
}
