#include "PHG4MapsCellReco.h"
#include "PHG4CylinderGeomContainer.h"
#include "PHG4CylinderGeom_MAPS.h"
#include "PHG4CylinderCell_MAPS.h"
#include "PHG4CylinderCellContainer.h"

#include "PHG4Cellv1.h"
#include "PHG4CellContainer.h"
#include "PHG4CellDefs.h"

#include <phparameter/PHParametersContainer.h>
#include <phparameter/PHParameterContainerInterface.h>


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
  PHParameterContainerInterface(name),
  detector(name),
  _timer(PHTimeServer::get()->insert_new(name)),
  chkenergyconservation(0)
{
  SetDefaultParameters();   // sets default timing window
  memset(nbins, 0, sizeof(nbins));

  if(Verbosity() > 0)  
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
  PHG4CellContainer *cells = findNode::getClass<PHG4CellContainer>(topNode , cellnodename);
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
      cells = new PHG4CellContainer();
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
  if (Verbosity() > 0)
    {
      geo->identify();
    }
  
  return Fun4AllReturnCodes::EVENT_OK;
}


int
PHG4MapsCellReco::process_event(PHCompositeNode *topNode)
{
  //cout << PHWHERE << "Entering process_event for PHG4MapsCellreco" << endl;
  
  _timer.get()->restart();
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  if (!g4hit)
    {
      cout << "Could not locate g4 hit node " << hitnodename << endl;
      exit(1);
    }
  PHG4CellContainer *cells = findNode::getClass<PHG4CellContainer>(topNode, cellnodename);
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

      if(Verbosity() > 2)
	layergeom->identify();

      // Get some layer parameters for later use
      double xpixw = layergeom->get_pixel_x();
      double xpixw_half = xpixw/2.0;
      double zpixw = layergeom->get_pixel_z();
      double zpixw_half = zpixw/2.0;
      int maxNX = layergeom->get_NX();
      int maxNZ = layergeom->get_NZ();

      for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
	{
	  //cout << "From PHG4MapsCellReco: Call hit print method: " << endl;
	  if(Verbosity() >4)
	    hiter->second->print();

	  // checking ADC timing integration window cut
	  if(*layer > 2)
	    {
	      cout  << PHWHERE << "Maps layers only go up to three! Quit." << endl;
	      exit(1);
	    }
	  if(Verbosity() > 1)
	    cout << " layer " << *layer << " t0 " << hiter->second->get_t(0) << " t1 " << hiter->second->get_t(1)
		 << " tmin " <<  tmin_max[*layer].first << " tmax " <<  tmin_max[*layer].second
		 << endl;
	  if (hiter->second->get_t(0) > tmin_max[*layer].second) continue;
	  if (hiter->second->get_t(1) < tmin_max[*layer].first) continue;

	  // get_property_int(const PROPERTY prop_id) const {return INT_MIN;}
	  int stave_number = hiter->second->get_property_int(PHG4Hit::prop_stave_index);
	  int half_stave_number = hiter->second->get_property_int(PHG4Hit::prop_half_stave_index);
	  int module_number = hiter->second->get_property_int(PHG4Hit::prop_module_index);
	  int chip_number = hiter->second->get_property_int(PHG4Hit::prop_chip_index);

	  TVector3 local_in( hiter->second->get_local_x(0),  hiter->second->get_local_y(0),  hiter->second->get_local_z(0) );
 	  TVector3 local_out( hiter->second->get_local_x(1),  hiter->second->get_local_y(1),  hiter->second->get_local_z(1) );
	  TVector3 midpoint( (local_in.X() + local_out.X()) / 2.0, (local_in.Y() + local_out.Y()) / 2.0, (local_in.Z() + local_out.Z()) / 2.0 );

	  if(Verbosity() > 4)
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

	  if(Verbosity() > 2)
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

	  if(pixel_number_in < 0 || pixel_number_out < 0)
	    {
	      cout << "Oops!  got negative pixel number in layer " << layergeom->get_layer()
		   << " pixel_number_in " << pixel_number_in
		   << " pixel_number_out " << pixel_number_out
		   << " local_in = " << local_in.X() << " " << local_in.Y() << " " << local_in.Z()
		   << " local_out = " << local_out.X() << " " << local_out.Y() << " " << local_out.Z()
		   << endl;

	    }

	  if(Verbosity() > 0)
	    cout << "entry pixel number " << pixel_number_in << " exit pixel number " << pixel_number_out << endl;

	  vector<int> vpixel;
	  vector<int> vxbin;
	  vector<int> vzbin;
	  vector<double> vlen;
	  vector< pair <double, double> > venergy;
	  //double trklen = 0.0;

	  //===================================================
	  // OK, now we have found which sensor the hit is in, extracted the hit
	  // position in local sensor coordinates,  and found the pixel numbers of the 
	  // entry point and exit point

	  //====================================================
	  // Beginning of charge sharing implementation
	  //    Find tracklet line inside sensor
	  //    Divide tracklet line into n segments (vary n until answer stabilizes) 
	  //    Find centroid of each segment
	  //    Diffuse charge at each centroid
	  //    Apportion charge between neighboring pixels
	  //    Add the pixel energy contributions from different track segments together
	  //====================================================

	  TVector3 pathvec = local_in - local_out;

	  // See figure 7.3 of the thesis by  Lucasz Maczewski (arXiv:10053.3710) for diffusion simulations in a MAPS epitaxial layer
	  // The diffusion widths below were inspired by those plots, corresponding to where the probability drops off to 1/3 of the peak value
	  // However note that we make the simplifying assumption that the probability distribution is flat within this diffusion width,
	  // while in the simulation it is not
	  //double diffusion_width_max = 35.0e-04;   // maximum diffusion radius 35 microns, in cm
	  //double diffusion_width_min = 12.0e-04;   // minimum diffusion radius 12 microns, in cm
	  double diffusion_width_max = 25.0e-04;   // maximum diffusion radius 35 microns, in cm
	  double diffusion_width_min = 8.0e-04;   // minimum diffusion radius 12 microns, in cm

	  double ydrift_max = pathvec.Y();
	  int nsegments = 4;

	  // we want to make a list of all pixels possibly affected by this hit
	  // we take the entry and exit locations in local coordinates, and build
	  // a rectangular array of pixels that encompasses both, with "nadd" pixels added all around

	  int xbin_in = layergeom->get_pixel_X_from_pixel_number(pixel_number_in);
	  int zbin_in = layergeom->get_pixel_Z_from_pixel_number(pixel_number_in);
	  int xbin_out = layergeom->get_pixel_X_from_pixel_number(pixel_number_out);
	  int zbin_out = layergeom->get_pixel_Z_from_pixel_number(pixel_number_out);
	  
	  int xbin_max, xbin_min;
	  int nadd = 2;
	  if(xbin_in > xbin_out)
	    {
	      xbin_max = xbin_in + nadd;
	      xbin_min = xbin_out - nadd;
	    }
	  else
	    {
	      xbin_max = xbin_out + nadd;
	      xbin_min = xbin_in - nadd;
	    }

	  int zbin_max, zbin_min;
	  if(zbin_in > zbin_out)
	    {
	      zbin_max = zbin_in + nadd;
	      zbin_min = zbin_out -nadd;
	    }
	  else
	    {
	      zbin_max = zbin_out + nadd;
	      zbin_min = zbin_in -nadd;
	    }


	  // need to check that values of xbin and zbin are within the valid range
	  if(xbin_min < 0) xbin_min = 0;
	  if(zbin_min < 0) zbin_min = 0;
	  if(xbin_max > maxNX) xbin_max = maxNX;
	  if(zbin_max > maxNZ) xbin_max = maxNZ;

	  if(Verbosity() > 1)
	    {
	      cout << " xbin_in " << xbin_in << " xbin_out " << xbin_out << " xbin_min " << xbin_min << " xbin_max " << xbin_max << endl;
	      cout << " zbin_in " << zbin_in << " zbin_out " << zbin_out << " zbin_min " << zbin_min << " zbin_max " << zbin_max << endl;
	    }

	  // skip this hit if it involves an unreasonable  number of pixels
	  // this skips it if either the xbin or ybin range traversed is greater than 8 (for 8 adding two pixels at each end makes the range 12) 
	  if(xbin_max - xbin_min > 12 || zbin_max - zbin_min > 12)
	    continue;

	  // this hit is skipped earlier if this dimensioning would be exceeded
	  double pixenergy[12][12] = {}; // init to 0
	  double pixeion[12][12] = {}; // init to 0

	  // Loop over track segments and diffuse charge at each segment location, collect energy in pixels
	  for(int i=0;i<nsegments;i++)
	    {
	      // Find the tracklet segment location
	      // If there are n segments of equal length, we want 2*n intervals
	      // The 1st segment is centered at interval 1, the 2nd at interval 3, the nth at interval 2n -1
	      double interval = 2 * (double ) i  + 1;
	      double frac = interval / (double) (2 * nsegments);
	      TVector3 segvec(pathvec.X() * frac, pathvec.Y() * frac, pathvec.Z() * frac);
	      segvec = segvec + local_out;
	      
	      //  Find the distance to the back of the sensor from the segment location
	      // That projection changes only the value of y
	      double ydrift = segvec.Y()  - local_out.Y();
	      
	      // Caculate the charge diffusion over this drift distance
	      // increases from diffusion width_min to diffusion_width_max 
	      double ydiffusion_radius = diffusion_width_min + (ydrift / ydrift_max) * (diffusion_width_max - diffusion_width_min);  
	      
	      if(Verbosity() > 5)
		cout << " segment " << i 
		     << " interval " << interval
		     << " frac " << frac
		     << " local_in.X " << local_in.X()
		     << " local_in.Z " << local_in.Z()
		     << " local_in.Y " << local_in.Y()
		     << " pathvec.X " << pathvec.X()
		     << " pathvec.Z " << pathvec.Z()
		     << " pathvec.Y " << pathvec.Y()
		     << " segvec.X " << segvec.X()
		     << " segvec.Z " << segvec.Z()
		     << " segvec.Y " << segvec.Y()
		     << " ydrift " << ydrift
		     << " ydrift_max " << ydrift_max
		     << " ydiffusion_radius " << ydiffusion_radius
		     << endl;
	      
	      // Now find the area of overlap of the diffusion circle with each pixel and apportion the energy
	      for(int ix = xbin_min;ix<=xbin_max;ix++)
		{
		  for(int iz = zbin_min;iz<=zbin_max;iz++)
		    {
		      // Find the pixel corners for this pixel number
		      int pixnum = layergeom->get_pixel_number_from_xbin_zbin( ix, iz );		  

		      if(pixnum < 0)
			{
			  cout << " pixnum < 0 , pixnum = " << pixnum << endl;
			  cout << " ix " << ix << " iz " << iz << endl;
			  cout << " xbin_min " << xbin_min << " zbin_min " << zbin_min
			       << " xbin_max " << xbin_max << " zbin_max " << zbin_max 
			       << endl;
			  cout << " maxNX " << maxNX << " maxNZ " << maxNZ
			       << endl;
			}

		      TVector3 tmp = layergeom->get_local_coords_from_pixel( pixnum );
		      // note that (x1,z1) is the top left corner, (x2,z2) is the bottom right corner of the pixel - circle_rectangle_intersection expects this ordering
		      double x1 = tmp.X() - xpixw_half;
		      double z1 = tmp.Z() + zpixw_half;
		      double x2 = tmp.X() + xpixw_half;
		      double z2 = tmp.Z() - zpixw_half;
		      
		      // here segvec.X and segvec.Z are the center of the circle, and diffusion_radius is the circle radius
		      // circle_rectangle_intersection returns the overlap area of the circle and the pixel. It is very fast if there is no overlap.
		      double pixarea_frac = circle_rectangle_intersection(x1, z1, x2, z2, segvec.X(), segvec.Z(), ydiffusion_radius) / (M_PI * pow(ydiffusion_radius,2) );
		      // assume that the energy is deposited uniformly along the tracklet length, so that this segment gets the fraction 1/nsegments of the energy
		      pixenergy[ix-xbin_min][iz-zbin_min] += pixarea_frac * hiter->second->get_edep() / (float) nsegments;
		      if (hiter->second->has_property(PHG4Hit::prop_eion))
			{
                          pixeion[ix-xbin_min][iz-zbin_min] += pixarea_frac * hiter->second->get_eion() / (float) nsegments;
			}
		      if(Verbosity() > 5)
			{
			  cout << "    pixnum " << pixnum << " xbin " << ix << " zbin " << iz
			       << " pixel_area fraction of circle " << pixarea_frac << " accumulated pixel energy " << pixenergy[ix-xbin_min][iz-zbin_min]
			       << endl;
			}
		    }
		}
	    }  // end loop over segments

	  // now we have the energy deposited in each pixel, summed over all tracklet segments. We make a vector of all pixels with non-zero energy deposited
	  for(int ix=xbin_min;ix<=xbin_max;ix++)
	    {
	      for(int iz=zbin_min;iz<=zbin_max;iz++)
		{
		  if( pixenergy[ix-xbin_min][iz-zbin_min] > 0.0 )
		    {	      
		      int pixnum = layergeom->get_pixel_number_from_xbin_zbin( ix, iz ); 
		      vpixel.push_back(pixnum);
		      vxbin.push_back(ix);
		      vzbin.push_back(iz);
		      pair <double,double> tmppair = make_pair(pixenergy[ix-xbin_min][iz-zbin_min],pixeion[ix-xbin_min][iz-zbin_min]);
		      venergy.push_back(tmppair);  	  
		      if(Verbosity() > 1)
			cout << " Added pixel number " << pixnum << " xbin " << ix << " zbin " << iz << " to vectors with energy " << pixenergy[ix-xbin_min][iz-zbin_min] << endl;
		    }		    	
		}
	    }
	  
	  //===================================
	  // End of charge sharing implementation
	  //===================================


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
	      PHG4Cell *cell = nullptr;
	      map<unsigned long long, PHG4Cell*>::iterator it;
	      it = celllist.find(inkey);
	      if (it != celllist.end())
		{
		  cell = it->second;
		}
	      else
		{

		  unsigned int index = celllist.size();
		  index++;
		  PHG4CellDefs::keytype key = PHG4CellDefs::MapsBinning::genkey(*layer,index);
		  cell = new PHG4Cellv1(key);
		  celllist[inkey] = cell;
		  cell->set_stave_index(stave_number);
		  cell->set_half_stave_index(half_stave_number);
		  cell->set_module_index(module_number);
		  cell->set_chip_index(chip_number);
		  cell->set_pixel_index(pixel_number);
		  cell->set_phibin(vxbin[i1]);
		  cell->set_zbin(vzbin[i1]);
		}
	      cell->add_edep(hiter->first, venergy[i1].first);
	      if (venergy[i1].second > 0)
		{
                  cell->add_eion(venergy[i1].second);
		}
	      cell->add_edep( venergy[i1].first);
	      
	      if(Verbosity() > 1)
		{
		  cout << " looping over fired cells: cell " << i1 << " inkey 0x" << hex << inkey << dec 
		       << " cell energy " << venergy[i1].first
		       << " cell edep " << cell->get_edep() << " total edep " << hiter->second->get_edep() << endl;
		}
	      
	    }
	} // end loop over g4hits
      
      int numcells = 0;
      for (map<unsigned long long, PHG4Cell *>::const_iterator mapiter = celllist.begin();mapiter != celllist.end() ; ++mapiter)
	{	  
	  cells->AddCell(mapiter->second);
	  numcells++;
	  
	  if (Verbosity() > 0)
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
      if (Verbosity() > 0)
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
  PHG4CellContainer *cells = findNode::getClass<PHG4CellContainer>(topNode, cellnodename);
  double sum_energy_g4hit = 0.;
  double sum_energy_cells = 0.;
  PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits();
  PHG4HitContainer::ConstIterator hiter;
  for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
    {
      sum_energy_g4hit += hiter->second->get_edep();
    }
  PHG4CellContainer::ConstRange cell_begin_end = cells->getCells();
  PHG4CellContainer::ConstIterator citer;
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
      if (Verbosity() > 0)
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

double  PHG4MapsCellReco::circle_rectangle_intersection( double x1, double y1,  double x2,  double y2,  double mx,  double my,  double r )
{
  // Find the area of overlap of a circle and rectangle 
  // Calls sA, which uses an analytic formula to determine the integral of the circle between limits set by the corners of the rectangle

  // move the rectangle to the frame where the circle is at (0,0)
  x1 -= mx; 
  x2 -= mx; 
  y1 -= my; 
  y2 -= my;

  if(Verbosity() > 7)
    {
      cout << " mx " << mx << " my " << my << " r " << r << " x1 " << x1 << " x2 " << x2 << " y1 " << y1 << " y2 " << y2 << endl;
      cout << " sA21 " << sA(r,x2,y1)
	   << " sA11 " << sA(r,x1,y1)
	   << " sA22 " << sA(r,x2,y2)
	   << " sA12 " << sA(r,x1,y2)
	   << endl;
    }

  return sA(r, x2, y1) - sA(r, x1, y1) - sA(r, x2, y2) + sA(r, x1, y2);
  
}

double  PHG4MapsCellReco::sA(double r, double x, double y) 
{
  // Uses analytic formula for the integral of a circle between limits set by the corner of a rectangle
  // It is called repeatedly to find the overlap area between the circle and rectangle
  // I found this code implementing the integral on a web forum called "ars technica",
  // https://arstechnica.com/civis/viewtopic.php?t=306492
  // posted by "memp"

  double a;

  if (x < 0) 
    {
      return -sA(r, -x, y);
    }
  
	if (y < 0) 
	  {
	    return -sA(r, x, -y);
	  }

	if (x > r) 
	  {
	    x = r;
	  }
	
	if (y > r) 
	  {
	    y = r;
	  }
	
	if (x*x + y*y > r*r) 
	  {
	    a = r*r*asin(x/r) + x*sqrt(r*r-x*x)
	      + r*r*asin(y/r) + y*sqrt(r*r-y*y)
	      - r*r*M_PI_2;
	    
	    a *= 0.5;
	  } 
	else 
	  {
	    a = x*y;
	  }
	
	return a;
}

void
PHG4MapsCellReco::set_timing_window(const int detid, const double tmin, const double tmax)
{
  if (Verbosity())
    cout << "PHG4MapsCellReco: Setting MAPS timing window parameters from macro for detid = " << detid << " to tmin = " << tmin << " tmax = " << tmax << endl;
  tmin_max.insert(std::make_pair(detid, std::make_pair(tmin, tmax)));

  return;
}

void
PHG4MapsCellReco::SetDefaultParameters()
{
  cout << "PHG4MapsCellReco: Setting MAPS timing window defaults to tmin = -2000 and  tmax = 2000 " << endl;
  for(int ilayer = 0;ilayer<3;ilayer++)
    tmin_max.insert(std::make_pair(ilayer, std::make_pair(-2000, 2000)));

  return;
}
