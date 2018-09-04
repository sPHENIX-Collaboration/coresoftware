#include "PHG4SiliconTrackerCellReco.h"
#include "PHG4CellContainer.h"
#include "PHG4Cellv1.h"
#include "PHG4CylinderCellGeom.h"
#include "PHG4CylinderCellGeomContainer.h"
#include "PHG4CylinderGeomContainer.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hitv1.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>

#include <PHG4CylinderGeom_Siladders.h>
#include <boost/format.hpp>
#include <TF1.h>
#include <TVector3.h>

#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace std;

PHG4SiliconTrackerCellReco::PHG4SiliconTrackerCellReco(const std::string &name)
  : SubsysReco(name)
  , _timer(PHTimeServer::get()->insert_new(name.c_str()))
  , chkenergyconservation(0)
  , tmin_default(-20.0)  // FVTX NIM paper Fig 32, collision has a timing spread around the triggered event. Accepting negative time too.
  ,  // ns
  tmax_default(80.0) // FVTX NIM paper Fig 32
  ,  // ns
  tmin_max()
{
  memset(nbins, 0, sizeof(nbins));
  Detector(name);

  hitnodename = "G4HIT_" + detector;
  cellnodename = "G4CELL_" + detector;
  geonodename = "CYLINDERGEOM_" + detector;
}

int PHG4SiliconTrackerCellReco::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
    exit(1);
  }

  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  if (!g4hit)
  {
    std::cout << "Could not locate g4 hit node " << hitnodename << std::endl;
    exit(1);
  }

  PHG4CellContainer *cells = findNode::getClass<PHG4CellContainer>(topNode, cellnodename);
  if (!cells)
  {
    PHNodeIterator dstiter(dstNode);

    PHCompositeNode *DetNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", detector));

    if (!DetNode)
    {
      DetNode = new PHCompositeNode(detector);
      dstNode->addNode(DetNode);
    }

    cells = new PHG4CellContainer();
    PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(cells, cellnodename.c_str(), "PHObject");
    DetNode->addNode(newNode);
  }

  PHG4CylinderGeomContainer *geo = findNode::getClass<PHG4CylinderGeomContainer>(topNode, geonodename.c_str());
  if (!geo)
  {
    std::cout << "Could not locate geometry node " << geonodename << std::endl;
    exit(1);
  }

  if (verbosity > 0)
    geo->identify();

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4SiliconTrackerCellReco::process_event(PHCompositeNode *topNode)
{
  _timer.get()->restart();

  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  if (!g4hit)
  {
    std::cout << "Could not locate g4 hit node " << hitnodename << std::endl;
    exit(1);
  }

  PHG4CellContainer *cells = findNode::getClass<PHG4CellContainer>(topNode, cellnodename);
  if (!cells)
  {
    std::cout << "could not locate cell node " << cellnodename << std::endl;
    exit(1);
  }
  cells->Reset();

  PHG4CylinderGeomContainer *geo = findNode::getClass<PHG4CylinderGeomContainer>(topNode, geonodename.c_str());
  if (!geo)
  {
    std::cout << "Could not locate geometry node " << geonodename << std::endl;
    exit(1);
  }

  // loop over all of the layers in the hit container
  // we need the geometry object for this layer
  if(verbosity > 2) cout << " PHG4SiliconTrackerCellReco: Loop over hits" << endl;
  PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits();
  for (PHG4HitContainer::ConstIterator hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
  {
    const int sphxlayer = hiter->second->get_detid();
    PHG4CylinderGeom_Siladders *layergeom = (PHG4CylinderGeom_Siladders*) geo->GetLayerGeom(sphxlayer);

    // checking ADC timing integration window cut
    // uses default values for now
    // these should depend on layer radius
    if (hiter->second->get_t(0) > tmax_default)
      continue;
    if (hiter->second->get_t(1) < tmin_default)
      continue;

    // I made this (small) diffusion up for now, we will get actual values for the INTT later
    double diffusion_width = 5.0e-04;   // diffusion radius 5 microns, in cm

    const int ladder_z_index = hiter->second->get_ladder_z_index();
    const int ladder_phi_index = hiter->second->get_ladder_phi_index();

     // What we have is a hit in the sensor. We have not yet assigned the strip(s) that were hit, we do that here
    //========================================================================

    // Get the entry point of the hit in sensor local coordinates
    TVector3 local_in( hiter->second->get_local_x(0),  hiter->second->get_local_y(0),  hiter->second->get_local_z(0) );
    TVector3 local_out( hiter->second->get_local_x(1),  hiter->second->get_local_y(1),  hiter->second->get_local_z(1) );
    TVector3 pathvec = local_in - local_out;

    int strip_y_index_in, strip_z_index_in, strip_y_index_out, strip_z_index_out;    
    layergeom->find_strip_index_values(ladder_z_index, local_in.y(), local_in.z(), strip_y_index_in, strip_z_index_in);
    layergeom->find_strip_index_values(ladder_z_index, local_out.y(), local_out.z(), strip_y_index_out, strip_z_index_out);

    if(verbosity > 5)
      {    
	// check to see if we get back the positions from these strip index values
	double check_location[3] = {-1,-1,-1};
	layergeom->find_strip_center_localcoords(ladder_z_index, strip_y_index_in, strip_z_index_in, check_location);
	cout << " G4 entry location = " << local_in.X() << "  " << local_in.Y() << "  " << local_in.Z() << endl; 
	cout << " Check entry location = " << check_location[0] << "  " << check_location[1] << "  " << check_location[2] << endl; 
	layergeom->find_strip_center_localcoords(ladder_z_index, strip_y_index_out, strip_z_index_out, check_location);
	cout << " G4 exit location = " << local_out.X() << "  " << local_out.Y() << "  " << local_out.Z() << endl; 
	cout << " Check exit location = " << check_location[0] << "  " << check_location[1] << "  " << check_location[2] << endl; 
      }

    // Now we find how many strips were crossed by this track, and divide the energy between them
    int  minstrip_z = strip_z_index_in;
    int maxstrip_z = strip_z_index_out;
    if(minstrip_z > maxstrip_z) swap(minstrip_z, maxstrip_z);

    int minstrip_y = strip_y_index_in;
    int maxstrip_y = strip_y_index_out;
    if(minstrip_y > maxstrip_y) swap(minstrip_y, maxstrip_y);

    // Use an algorithm similar to the one for the MVTX pixels, since it facilitates adding charge diffusion
    // for now we assume small charge diffusion
    
    vector<int> vybin;
    vector<int> vzbin;
    //vector<double> vlen;
    vector< pair <double, double> > venergy;
    
    //====================================================
    // Beginning of charge sharing implementation
    //    Find tracklet line inside sensor
    //    Divide tracklet line into n segments (vary n until answer stabilizes) 
    //    Find centroid of each segment
    //    Diffuse charge at each centroid
    //    Apportion charge between neighboring pixels
    //    Add the pixel energy contributions from different track segments together
    //====================================================
    
    // skip this hit if it involves an unreasonable  number of pixels
    // this skips it if either the xbin or ybin range traversed is greater than 8 (for 8 adding two pixels at each end makes the range 12) 
    if(maxstrip_y - minstrip_y > 12 || maxstrip_z - minstrip_z > 12)
      continue;
    
    // this hit is skipped above if this dimensioning would be exceeded
    double stripenergy[12][12] = {}; // init to 0
    double stripeion[12][12] = {}; // init to 0
    
    int nsegments = 10;
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
	
	// Caculate the charge diffusion over this drift distance
	// increases from diffusion width_min to diffusion_width_max 
	double diffusion_radius = diffusion_width;   
	
	if(verbosity > 5)
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
	       << " diffusion_radius " << diffusion_radius
	       << endl;
	
	// Now find the area of overlap of the diffusion circle with each pixel and apportion the energy
	for(int iz=minstrip_z; iz <= maxstrip_z; iz++)
	  {
	    for(int iy = minstrip_y; iy <= maxstrip_y; iy++)
	      {
		// Find the pixel corners for this pixel number
		double location[3] = {-1,-1,-1};
		layergeom->find_strip_center_localcoords(ladder_z_index, iy, iz, location);
		// note that (y1,z1) is the top left corner, (y2,z2) is the bottom right corner of the pixel - circle_rectangle_intersection expects this ordering		
		double y1 = location[1] - layergeom->get_strip_y_spacing() / 2.0;
		double y2 = location[1] + layergeom->get_strip_y_spacing() / 2.0;
		double z1 = location[2] + layergeom->get_strip_z_spacing() / 2.0;
		double z2 = location[2] - layergeom->get_strip_z_spacing() / 2.0;
		
		// here segvec.Y and segvec.Z are the center of the circle, and diffusion_radius is the circle radius
		// circle_rectangle_intersection returns the overlap area of the circle and the pixel. It is very fast if there is no overlap.
		double striparea_frac = circle_rectangle_intersection(y1, z1, y2, z2, segvec.Y(), segvec.Z(), diffusion_radius) / (M_PI * pow(diffusion_radius,2) );
		// assume that the energy is deposited uniformly along the tracklet length, so that this segment gets the fraction 1/nsegments of the energy
		stripenergy[iy-minstrip_y][iz-minstrip_z] += striparea_frac * hiter->second->get_edep() / (float) nsegments;
		if (hiter->second->has_property(PHG4Hit::prop_eion))
		  {
		    stripeion[iy-minstrip_y][iz-minstrip_z] += striparea_frac * hiter->second->get_eion() / (float) nsegments;
		  }
		if(verbosity > 5)
		  {
		    cout << "    strip y index " << iy <<  " strip z index  " << iz
			 << " strip area fraction of circle " << striparea_frac << " accumulated pixel energy " << stripenergy[iy-minstrip_y][iz-minstrip_z]
			 << endl;
		  }
	      }
	  }
      }  // end loop over segments
    
    // now we have the energy deposited in each pixel, summed over all tracklet segments. We make a vector of all pixels with non-zero energy deposited
    for(int iz=minstrip_z; iz <= maxstrip_z; iz++)
      {
	for(int iy = minstrip_y; iy <= maxstrip_y; iy++)
	  {
	    if( stripenergy[iy-minstrip_y][iz-minstrip_z] > 0.0 )
	      {	      
		vybin.push_back(iy);
		vzbin.push_back(iz);
		pair <double,double> tmppair = make_pair(stripenergy[iy-minstrip_y][iz-minstrip_z],stripeion[iy-minstrip_y][iz-minstrip_z]);
		venergy.push_back(tmppair);  	  
		if(verbosity > 1)
		  cout << " Added ybin " << iy << " zbin " << iz << " to vectors with energy " << stripenergy[iy-minstrip_y][iz-minstrip_z] << endl;
	      }		    	
	  }
      }
    
    //===================================
    // End of charge sharing implementation
    //===================================
    
    // Add the strips fired by this hit to the cell list
    //===============================
    
    for (unsigned int i1 = 0; i1 < vybin.size(); i1++)   // loop over all fired cells
      {
	// this string must be unique - it needs the layer too, or in high multiplicity events it will add g4 hits in different layers with the same key together
	std::string key = boost::str(boost::format("%d-%d-%d-%d-%d") % sphxlayer % ladder_z_index % ladder_phi_index % vzbin[i1] % vybin[i1]).c_str();
	PHG4Cell *cell = nullptr;
	map<string, PHG4Cell *>::iterator it;
	
	it = celllist.find(key);
	// If there is an existing cell to add this hit to, find it    
	if (it != celllist.end())
	  {
	    cell = it->second;
	    if(verbosity > 2)  
	      cout << " found existing cell with key " << key << endl;
	  }
	
	// There is not an existing cell to add this hit to, start a new cell    
	if(!cell)
	  {
	    if(verbosity > 2) cout << " did not find existing cell with key " << key << " start a new one" << endl;
	    unsigned int index = celllist.size();
	    index++;
	    PHG4CellDefs::keytype cellkey = PHG4CellDefs::MapsBinning::genkey(sphxlayer, index);
	    cell = new PHG4Cellv1(cellkey);
	    celllist[key] = cell;
	    // This encodes the z and phi position of the sensor
	    //          celllist[key]->set_sensor_index(boost::str(boost::format("%d_%d") %ladder_z_index %ladder_phi_index).c_str());
	    
	    cell->set_ladder_z_index(ladder_z_index);
	    cell->set_ladder_phi_index(ladder_phi_index);
	    
	    // The z and phi position of the hit strip within the sensor
	    cell->set_zbin(vzbin[i1]);
	    cell->set_phibin(vybin[i1]);
	  }
	
	// One way or another we have a cell pointer - add this hit to the cell
	cell->add_edep(venergy[i1].first);
	cell->add_edep(hiter->first, venergy[i1].first);  // adds the g4hit association to the cell 
	cell->add_eion(venergy[i1].second);
      }

  }  // end loop over g4hits

  
  int numcells = 0;
  for (std::map<std::string, PHG4Cell *>::const_iterator mapiter = celllist.begin(); mapiter != celllist.end(); ++mapiter)
  {
    cells->AddCell(mapiter->second);
    numcells++;

    if (verbosity > 0)
      {
	std::cout << "Adding cell for "
		  << " layer " << mapiter->second->get_layer()
		  << " ladder z index: " << mapiter->second->get_ladder_z_index()
		  << ", ladder phi index: " << mapiter->second->get_ladder_phi_index()
		  << ", srip z: " << mapiter->second->get_zbin()
		  << ", strip y: " << mapiter->second->get_phibin()
		  << ", energy dep: " << mapiter->second->get_edep()
		  << std::endl;
      }
  }
  celllist.clear();
  
  if (verbosity > 0)
    std::cout << Name() << ": found " << numcells << " silicon strips with energy deposition" << std::endl;

  if (chkenergyconservation)
  {
    CheckEnergy(topNode);
  }
  _timer.get()->stop();
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4SiliconTrackerCellReco::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4SiliconTrackerCellReco::CheckEnergy(PHCompositeNode *topNode)
{
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  PHG4CellContainer *cells = findNode::getClass<PHG4CellContainer>(topNode, cellnodename);
  double sum_energy_g4hit = 0.;
  double sum_energy_cells = 0.;

  PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits();
  PHG4HitContainer::ConstIterator hiter;
  for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
    sum_energy_g4hit += hiter->second->get_edep();

  PHG4CellContainer::ConstRange cell_begin_end = cells->getCells();
  PHG4CellContainer::ConstIterator citer;
  for (citer = cell_begin_end.first; citer != cell_begin_end.second; ++citer)
    sum_energy_cells += citer->second->get_edep();

  // the fractional eloss for particles traversing eta bins leads to minute rounding errors
  if (fabs(sum_energy_cells - sum_energy_g4hit) / sum_energy_g4hit > 1e-6)
  {
    std::cout << "energy mismatch between cells: " << sum_energy_cells
              << " and hits: " << sum_energy_g4hit
              << " diff sum(cells) - sum(hits): " << sum_energy_cells - sum_energy_g4hit
              << std::endl;

    return -1;
  }
  else
  {
    if (verbosity > 0)
      std::cout << Name() << ": total energy for this event: " << sum_energy_g4hit << " GeV" << std::endl;
  }
  return 0;
}

double  PHG4SiliconTrackerCellReco::circle_rectangle_intersection( double x1, double y1,  double x2,  double y2,  double mx,  double my,  double r )
{
  // Find the area of overlap of a circle and rectangle 
  // Calls sA, which uses an analytic formula to determine the integral of the circle between limits set by the corners of the rectangle

  // move the rectangle to the frame where the circle is at (0,0)
  x1 -= mx; 
  x2 -= mx; 
  y1 -= my; 
  y2 -= my;

  if(verbosity > 7)
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

double  PHG4SiliconTrackerCellReco::sA(double r, double x, double y) 
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
