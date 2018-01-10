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

  // ADF: removed this because it hard-codes the layer numbers!
  // binning.insert(std::make_pair(2, 0));
  // binning.insert(std::make_pair(3, 1));
  // binning.insert(std::make_pair(4, 2));
  // binning.insert(std::make_pair(5, 3));
  // for (std::map<int,int>::iterator iter = binning.begin(); iter != binning.end(); ++iter)
  //   {
  //     // if the user doesn't set an integration window, set the default
  //     tmin_max.insert(std::make_pair(/*layer*/iter->first, std::make_pair(tmin_default, tmax_default)));
  //   }

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
    PHG4CylinderGeom *layergeom = geo->GetLayerGeom(sphxlayer);

    // checking ADC timing integration window cut
    // uses default values for now
    // these should depend on layer radius
    if (hiter->second->get_t(0) > tmax_default)
      continue;
    if (hiter->second->get_t(1) < tmin_default)
      continue;

    const int ladder_z_index = hiter->second->get_ladder_z_index();
    const int ladder_phi_index = hiter->second->get_ladder_phi_index();
    const int strip_z_index = hiter->second->get_strip_z_index();
    const int strip_y_index = hiter->second->get_strip_y_index();

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
