#include "PHG4SiliconTrackerCellReco.h"
#include "PHG4CylinderGeomContainer.h"
#include "PHG4CylinderGeomv4.h"
#include "PHG4CylinderCellGeomContainer.h"
#include "PHG4CylinderCellGeom.h"
#include "PHG4CylinderCellv2.h"
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
#include <fun4all/getClass.h>


#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>

using namespace std;

PHG4SiliconTrackerCellReco::PHG4SiliconTrackerCellReco(const string &name) :
  SubsysReco(name),
  _timer(PHTimeServer::get()->insert_new(name.c_str())),
  chkenergyconservation(0)
{
  memset(nbins, 0, sizeof(nbins));
  Detector(name);
}

int PHG4SiliconTrackerCellReco::InitRun(PHCompositeNode *topNode)
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

  /*
  seggeonodename = "CYLINDERCELLGEOM_" + detector;
  PHG4CylinderCellGeomContainer *seggeo = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode , seggeonodename.c_str());
  if (!seggeo)
    {
      seggeo = new PHG4CylinderCellGeomContainer();
      PHCompositeNode *runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN" ));
      PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(seggeo, seggeonodename.c_str() , "PHObject");
      runNode->addNode(newNode);
    }

  map<int, PHG4CylinderGeom *>::const_iterator miter;
  pair <map<int, PHG4CylinderGeom *>::const_iterator, map<int, PHG4CylinderGeom *>::const_iterator> begin_end = geo->get_begin_end();
  map<int, std::pair <double, int> >::iterator sizeiter;
  for (miter = begin_end.first; miter != begin_end.second; ++miter)
    {
      PHG4CylinderGeom *layergeom = miter->second;
      int layer = layergeom->get_layer();
      int numslats = layergeom->get_nscint();
      // create geo object and fill with variables common to all binning methods
      PHG4CylinderCellGeom *layerseggeo = new PHG4CylinderCellGeom();
      layerseggeo->set_layer(layergeom->get_layer());
      layerseggeo->set_radius(layergeom->get_radius());
      layerseggeo->set_thickness(layergeom->get_thickness());
      int nslatbins = numslats / nslatscombined;
	  // check if the number of slats is a multiple of our number of slats read out
	  // if not print out a warning and drop the remaining slats 
	  // from reconstruction since
	  // we don't want weird differently sized towers 
          if (numslats % nslatscombined)
            {
	      //              nslatbins++;
              cout << Name() << ": total number of slats " << numslats 
                   << " not multiple of readout combined slats "
                   <<  nslatscombined << " dropping the last " 
                   << numslats % nslatscombined
                   << " slats from reconstruction" << endl;
	      cout << "This will affect clusters which span the rollover phi range" << endl;
            }
	  if (verbosity > 1)
	    {
	      layergeom->identify();
	    }
          layerseggeo->set_binning(phg4cylindercelldefs::etaslatbinning);
          layerseggeo->set_etabins(22);
          layerseggeo->set_etamin(-1.1);
          layerseggeo->set_etastep((1.1+1.1)/22.);
          layerseggeo->set_phimin(layergeom->get_phi_slat_zero());
	  double deltaphi = 2.*M_PI/layergeom->get_nscint();
	  layerseggeo->set_phistep(deltaphi*nslatscombined);
	  layerseggeo->set_phibins(nslatbins);
          etastep[layer] = (1.1+1.1)/22.;
      // add geo object filled by different binning methods
      seggeo->AddLayerCellGeom(layerseggeo);
      if (verbosity > 1)
	{
	  layerseggeo->identify();
	}
    }
  */

  return Fun4AllReturnCodes::EVENT_OK;
}


int
PHG4SiliconTrackerCellReco::process_event(PHCompositeNode *topNode)
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
      //cout << "---------- SiliconTrackerReco:  Looping over layers " << endl;

      // loop over the hits in this layer
      PHG4HitContainer::ConstIterator hiter;
      PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits(*layer);

      // we need the geometry object for this layer
      PHG4CylinderGeom *layergeom = geo->GetLayerGeom(*layer);
      //layergeom->identify();

      for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
	{
	  int ladder_z_index = hiter->second->get_ladder_z_index();
	  int ladder_phi_index = hiter->second->get_ladder_phi_index();
	  int strip_z_index = hiter->second->get_strip_z_index();
	  int strip_y_index = hiter->second->get_strip_y_index();

	  // As a check, get the positions of the hit strips from the geo object
	  double location[3] = {-1,-1,-1};
	  layergeom->find_strip_center(ladder_z_index, ladder_phi_index, strip_z_index, strip_y_index, location);

	  // cout << "      Found location from geometry for ladder_z_index " << ladder_z_index 
	  //      << " ladder_phi_index " << ladder_phi_index 
	  //      << " strip_z_index " << strip_z_index 
	  //      << " strip_y_index " << strip_y_index
	  //      <<  " is:  x = " << location[0]
	  //      << " y = " << location[1]
	  //      << " z  = " << location[2]
	  //      << endl;

	  // combine ladder index values to get a single integer
	  //     A ladder consists of ladder segments butted end to end in z
	  //     Each ladder segment contains one sensor, each sensor has multiple columns of tracking strips
	  // ladder_z index is the z bin for the ladder segment containing the sensor with this hit strip
	  // ladder_phi index is the phi bin for the ladder segment containing the sensor with this hit strip
	  // strip_z_index is the strip column inside the sensor
	  // strip_y_index is the number of the strip i the column
	  char inkey[1024];
	  sprintf(inkey,"%i-%i_%i_%i",ladder_z_index, ladder_phi_index, strip_z_index, strip_y_index);
	  std::string key(inkey);

	  if (celllist.count(key) > 0) {
	    celllist[key]->add_edep(hiter->first, hiter->second->get_edep());
	  } else {
	    celllist[key] = new PHG4CylinderCellv2();
	    celllist[key]->set_layer(*layer);

	    // This encodes the z and phi position of the sensor 
	    ostringstream name;
	    name << ladder_z_index << "_" << ladder_phi_index;
	    std::string sensor_index = name.str();
	    celllist[key]->set_sensor_index(sensor_index);

	    celllist[key]->set_ladder_z_index(ladder_z_index);
	    celllist[key]->set_ladder_phi_index(ladder_phi_index);
	    
	    // The z and phi position of the hit strip within the sensor
	    celllist[key]->set_zbin(strip_z_index);
	    celllist[key]->set_phibin(strip_y_index);	  

	    celllist[key]->add_edep(hiter->first, hiter->second->get_edep());
	  }
	} // end loop over g4hits

      int numcells = 0;
      for (map<std::string, PHG4CylinderCell *>::const_iterator mapiter = celllist.begin();mapiter != celllist.end() ; ++mapiter)
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
PHG4SiliconTrackerCellReco::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int
PHG4SiliconTrackerCellReco::CheckEnergy(PHCompositeNode *topNode)
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

