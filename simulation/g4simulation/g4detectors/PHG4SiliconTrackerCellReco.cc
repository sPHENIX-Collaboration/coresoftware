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
    const int sphxlayer = hiter->second->get_layer();
    PHG4CylinderGeom_Siladders *layergeom = (PHG4CylinderGeom_Siladders*) geo->GetLayerGeom(sphxlayer);

    // checking ADC timing integration window cut
    // uses default values for now
    // these should depend on layer radius
    if (hiter->second->get_t(0) > tmax_default)
      continue;
    if (hiter->second->get_t(1) < tmin_default)
      continue;

    const int ladder_z_index = hiter->second->get_ladder_z_index();
    const int ladder_phi_index = hiter->second->get_ladder_phi_index();

     // What we have here is a hit in the sensor. We have not yet assigned the strip(s) that were hit, we do that here.
    //===========================================================================

    // Get the entry point of the hit in sensor local coordinates
    TVector3 local_in( hiter->second->get_local_x(0),  hiter->second->get_local_y(0),  hiter->second->get_local_z(0) );
    TVector3 local_out( hiter->second->get_local_x(1),  hiter->second->get_local_y(1),  hiter->second->get_local_z(1) );
    //TVector3 midpoint( (local_in.X() + local_out.X()) / 2.0, (local_in.Y() + local_out.Y()) / 2.0, (local_in.Z() + local_out.Z()) / 2.0 );
    TVector3 pathvec = local_in - local_out;
    double path_length = pathvec.Mag();

    int strip_y_index_in, strip_z_index_in, strip_y_index_out, strip_z_index_out;    
    layergeom->find_strip_index_values(ladder_z_index, local_in.y(), local_in.z(), strip_y_index_in, strip_z_index_in);
    layergeom->find_strip_index_values(ladder_z_index, local_out.y(), local_out.z(), strip_y_index_out, strip_z_index_out);
    
    // Now we find how many strips were crossed by this track, and divide the energy between them

    int  minstrip_z = strip_z_index_in;
    int maxstrip_z = strip_z_index_out;
    if(minstrip_z > maxstrip_z) swap(minstrip_z, maxstrip_z);

    int minstrip_y = strip_y_index_in;
    int maxstrip_y = strip_y_index_out;
    if(minstrip_y > maxstrip_y) swap(minstrip_y, maxstrip_y);

    cout << " minstrip_y " << minstrip_y << " maxstrip_y " << maxstrip_y << " minstrip_z " << minstrip_z << " maxstrip_z " << maxstrip_z << endl;

    // define the (assumed) straight line of the path
    double dzdy = ( local_out.Z() - local_in.Z() ) / ( local_out.Y() - local_in.Y() );
    double z0 =  local_in.Z()  - dzdy * local_in.Y();
    TF1* fpath1 = new TF1("fpath1", "[0] + [1]*x");
    fpath1->SetParameter(0, z0);
    fpath1->SetParameter(1, dzdy);

    double dydz = ( local_out.Y() - local_in.Y() ) / ( local_out.Z() - local_in.Z() );
    double y0 =  local_in.Y()  - dydz * local_in.Z();
    TF1* fpath2 = new TF1("fpath2", "[0] + [1]*x");
    fpath2->SetParameter(0, y0);
    fpath2->SetParameter(1, dydz);

    // loop over the strips touched by this track
    for(int iz=minstrip_z; iz <= maxstrip_z; iz++)
      {
	for(int iy = minstrip_y; iy <= maxstrip_y; iy++)
	  {
	    double point[2][2] = {-10000, -10000, -10000, -10000};

	    double location[3] = {-1,-1,-1};
	    layergeom->find_strip_center_localcoords(ladder_z_index, iy, iz, location);
	    // define the rectangular box corresponding to this strip
	    double strip_ymin = location[1] - layergeom->get_strip_y_spacing() / 2.0;
	    double strip_ymax = location[1] + layergeom->get_strip_y_spacing() / 2.0;
	    double strip_zmin = location[2] - layergeom->get_strip_z_spacing() / 2.0;
	    double strip_zmax = location[2] + layergeom->get_strip_z_spacing() / 2.0;

	    // does the line pass through this box?
	    int got_it = 0;
	    double z1 = fpath1->Eval(strip_ymin);
	    if( z1 > strip_zmin && z1 < strip_zmax) 
	      {
		// crosses the boundary at ymin
		point[got_it][0] = strip_ymin;
		point[got_it][1] = z1; 
		got_it++;
	      }
	      double z2 = fpath1->Eval(strip_ymax);
	      if (z2 > strip_zmin && z2 < strip_zmax)
	      {
		// crosses the boundary at ymax
		point[got_it][0] = strip_ymax;
		point[got_it][1] = z1; 
		got_it++;
	      }

	    double y1 = fpath2->Eval(strip_zmin);
	    if( (y1 > strip_ymin && y1 < strip_ymax))
	      {
		// crosses the boundary at zmin
		point[got_it][0] = y1;
		point[got_it][1] = strip_zmin; 
		got_it++;
	      }
	    double y2 = fpath2->Eval(strip_zmax);
	    if( (y2 > strip_ymin && y2 < strip_ymax) )
	      {
		// crosses the boundary at zmax
		point[got_it][0] = y2;
		point[got_it][1] = strip_zmax; 
		got_it++;
	      }

	    // the path does not enter the volume
	    if(got_it == 0)
	      {
		cout << " track did not intersect z or y boundary of strip with iy " << iy << " iz " << iz << "  ---  got_it = " << got_it << endl; 
		continue;
	      }
	      if(got_it != 2)
		{
		  cout << "Oops! That's not right: got_it = " << got_it << endl;
		  continue;
		}

	      // now we can get the energy deposit in this strip
	      double plength = sqrt( pow( (point[0][0] - point[1][0]), 2) + pow( (point[0][1] - point[1][1]), 2)  );
	      double efrac = plength/path_length;
	      cout << " Make cell for layer " << sphxlayer << " ladder_z_index " << ladder_z_index << " ladder_phi_index " << ladder_phi_index << endl;
	      cout << "         strip_y_index " << iy << " strip_z_index " << iz << endl; 
	      cout << "         path_length " << path_length << " plength " << plength << " efrac " << efrac << " strip energy " << efrac * hiter->second->get_edep()
		   << endl;

	    // Add this strip to the cell list
	    //====================

	    // this string must be unique - it needs the layer too, or in high multiplicity events it will add g4 hits in different layers with the same key together
	    std::string key = boost::str(boost::format("%d-%d-%d-%d-%d") % sphxlayer % ladder_z_index % ladder_phi_index % iz % iy).c_str();
	    PHG4Cell *cell = nullptr;
	    map<string, PHG4Cell *>::iterator it;
	    
	    it = celllist.find(key);
	    // If there is an existing cell to add this hit to, find it    
	    if (it != celllist.end())
	      {
		cell = it->second;
		if(verbosity > 2)  cout << " found existing cell with key " << key << endl;
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
		cell->set_zbin(iz);
		cell->set_phibin(iy);
	      }
	    
	    // One way or another we have a cell pointer - add this hit to the cell
	    cell->add_edep(hiter->first, efrac * hiter->second->get_edep());
	    cell->add_edep(efrac * hiter->second->get_edep());
	  }
      }
    
    


    //=============================================================================
    /*    
    // As a check, get the positions of the hit strips from the geo object
    double location[3] = {-1, -1, -1};
    layergeom->find_strip_center(ladder_z_index, ladder_phi_index, strip_z_index, strip_y_index, location);

    if(verbosity > 2) 
    {
	cout << endl << "  g4 hit:  layer " <<  hiter->second->get_layer() << " edep " <<  hiter->second->get_edep() << endl;
	cout << "   Hit entry point x,y,z = " << hiter->second->get_x(0) << "  " << hiter->second->get_y(0) << "  " << hiter->second->get_z(0) << endl;
	cout << "   Hit exit point x,y,z = " << hiter->second->get_x(1) << "  " << hiter->second->get_y(1) << "  " << hiter->second->get_z(1) << endl;
	cout << "  ladder z index " <<  hiter->second->get_ladder_z_index() << " ladder phi index " <<  hiter->second->get_ladder_phi_index() 
	     << " strip z index " <<  hiter->second->get_strip_z_index() << " strip y index " <<   hiter->second->get_strip_y_index() << endl;
	cout << "   strip x,y,z from geometry object = " << location[0] << "  " << location[1] << "  " << location[2] << endl;
	cout << endl;
    }

    // this string must be unique - it needs the layer too, or in high multiplicity events it will add g4 hits in different layers with the same key together
    std::string key = boost::str(boost::format("%d-%d-%d-%d-%d") % sphxlayer % ladder_z_index % ladder_phi_index % strip_z_index % strip_y_index).c_str();
    PHG4Cell *cell = nullptr;
    map<string, PHG4Cell *>::iterator it;

    it = celllist.find(key);
    // If there is an existing cell to add this hit to, find it    
    if (it != celllist.end())
      {
	cell = it->second;
	if(verbosity > 2)  cout << " found existing cell with key " << key << endl;
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
	cell->set_zbin(strip_z_index);
	cell->set_phibin(strip_y_index);
      }

    // One way or another we have a cell pointer - add this hit to the cell
    cell->add_edep(hiter->first, hiter->second->get_edep());
    cell->add_edep(hiter->second->get_edep());
    */

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
                << ", srip z index: " << mapiter->second->get_zbin()
                << ", strip y index: " << mapiter->second->get_phibin()
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
