#include "PHG4SvtxAddConnectedCells.h"
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>
#include <g4detectors/PHG4CylinderCellContainer.h>
#include <g4detectors/PHG4CylinderCellGeomContainer.h>
#include <g4detectors/PHG4CylinderCell.h>
#include <g4detectors/PHG4CylinderCellv1.h>
#include <g4detectors/PHG4CylinderCellGeom.h>

#include <boost/tuple/tuple.hpp>
#include <boost/format.hpp>

#define BOOST_NO_HASH // Our version of boost.graph is incompatible with GCC-4.3 w/o this flag
#include <boost/bind.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
using namespace boost;

#include <iostream>
#include <stdexcept>
#include <cmath>

using namespace std;

static const float twopi = 2.0*M_PI;

bool PHG4SvtxAddConnectedCells::lessthan(const PHG4CylinderCell* lhs, 
				   const PHG4CylinderCell* rhs) {

  if( lhs->get_binphi() < rhs->get_binphi() ) return true;
  else if( lhs->get_binphi() == rhs->get_binphi() ){
    if( lhs->get_binz() < rhs->get_binz() ) return true;
  }

  return false;
}

PHG4SvtxAddConnectedCells::PHG4SvtxAddConnectedCells(const string &name) :
  SubsysReco(name),
  connected_phi_offset(128), // offset in phi bins for connected cells
  _cells(NULL),
  _geom_container(NULL),
  _timer(PHTimeServer::get()->insert_new(name)) {

  memset(ncells_connected,0,sizeof(ncells_connected));

  // Default values. Override using the setters if desired
  ncells_connected[4] = 2;
  ncells_connected[5] = 2;
  ncells_connected[6] = 5;

  if (verbosity > 0) {
    cout << "PHG4SvtxAddConnectedCells: Default settings are: " << endl;
    cout << "phi bins offset for connected cells " << connected_phi_offset << endl;
    for(int i=0; i<7;i++)
      cout << " number of connected cells for layer " << i << " is " << ncells_connected[i] << endl;
  }
}

int PHG4SvtxAddConnectedCells::InitRun(PHCompositeNode* topNode) {

  // get the SVX geometry object
  _geom_container = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode,"CYLINDERCELLGEOM_SVTX");
  if (!_geom_container)
    {
      cout << PHWHERE
	   << " ERROR: Can't find PHG4CylinderCellGeomContainerNode."
	   << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

  return Fun4AllReturnCodes::EVENT_OK;
}

int PHG4SvtxAddConnectedCells::process_event(PHCompositeNode *topNode) {

  _timer.get()->restart();

  //----------
  // Get Nodes
  //----------

  _cells = findNode::getClass<PHG4CylinderCellContainer>(topNode,"G4CELL_SVTX");
  if (!_cells) {
    cout << PHWHERE << " ERROR: Can't find G4CELL_SVTX." << endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  // Since the G4 layers don't necessarily correspond to the silicon
  // layers, and don't necessarily start from zero (argh), we create
  // our own layers numbers that are consecutive starting from zero.
  //
  PHG4CylinderCellGeomContainer::ConstRange layerrange = _geom_container->get_begin_end();
  for(PHG4CylinderCellGeomContainer::ConstIterator layeriter = layerrange.first;
      layeriter != layerrange.second;
      ++layeriter) {

    int layer = layeriter->second->get_layer();

    // Different layers have different numbers of cells connected together
    // nconnected is the number of additional cells to add
    int nconnected = ncells_connected[layer];

    vector<PHG4CylinderCell*> cell_list;
    PHG4CylinderCellContainer::ConstRange cellrange = _cells->getCylinderCells(layer);
    for(PHG4CylinderCellContainer::ConstIterator celliter = cellrange.first;
        celliter != cellrange.second;
        ++celliter) {

      PHG4CylinderCell* cell = celliter->second;
      cell_list.push_back(cell);
    }

    if (cell_list.size() == 0) continue; // if no cells, go to the next layer

    sort(cell_list.begin(), cell_list.end(), PHG4SvtxAddConnectedCells::lessthan);

    // Get the zbin and phibin for each hit cell, and add connected cells

    for(unsigned int i=0; i<cell_list.size(); i++) 
      {
	int ibinz = cell_list[i]->get_binz();
	int ibinphi = cell_list[i]->get_binphi();

	PHG4HitDefs::keytype g4hitid = ~0;
	int foundit = 0;
	PHG4CylinderCell::EdepConstRange maprange = cell_list[i]->get_g4hits();
	for (PHG4CylinderCell::EdepConstIterator hititer = maprange.first;
	     hititer != maprange.second;
	     hititer++) {
	  g4hitid = hititer->first;
	  foundit = 1;
	}

	// Should be impossible
	if(!foundit)
	  {
	    cout << "Oops! G4 hitid not found" << endl;
	    continue;
	  }

	if(verbosity > 0)	
	  cout << "   Found cell: layer " << cell_list[i]->get_layer() 
	       << " cell binz " << cell_list[i]->get_binz() 
	       << " binphi " << cell_list[i]->get_binphi()
	       << " cellid " << cell_list[i]->get_cell_id()
	       << " edep " << cell_list[i]->get_edep()
	       << " g4hitid " << g4hitid 
	       << endl;
	
	for(int ic=0; ic<nconnected;ic++)
	  {
	    int newphibin = ibinphi - (ic+1) * connected_phi_offset;
	    if(newphibin < 0)
	      newphibin = ibinphi + (ic+1) * connected_phi_offset;
	    
	    if (verbosity > 1)
	      {
		cout << "Adding connected cell in bin phi: " << newphibin
		     << ", z bin: " << ibinz
		     << ", energy dep: " << cell_list[i]->get_edep()
		     << endl;
	      }
	    
	    // Note - this does not check that it is not already a hit cell! 
	    PHG4CylinderCell *newcell = new PHG4CylinderCellv1();
	    newcell->set_layer(layer);
	    newcell->set_phibin(newphibin);
	    newcell->set_zbin(ibinz);
	    
	    newcell->add_edep(g4hitid, cell_list[i]->get_edep());
	    
	    _cells->AddCylinderCell(layer, newcell);

	  } // end loop over new cells to be added

      } // end loop over cell list

    if(verbosity > 1)
      {
	cout << "New list of cells for layer " << layer << endl;
	
	PHG4CylinderCellContainer::ConstRange cellrange_now = _cells->getCylinderCells(layer);
	for(PHG4CylinderCellContainer::ConstIterator celliter = cellrange_now.first;
	    celliter != cellrange_now.second;
	    ++celliter) {
	  
	  PHG4CylinderCell* cell = celliter->second;
	  
	  // Get the G4 hitid
	  PHG4HitDefs::keytype g4hitid = ~0;
	  PHG4CylinderCell::EdepConstRange maprange = cell->get_g4hits();
	  for (PHG4CylinderCell::EdepConstIterator hititer = maprange.first;
	       hititer != maprange.second;
	       hititer++) {
	    g4hitid = hititer->first;
	  }
	  
	  cout << " layer " << cell->get_layer()
	       << " binz " << cell->get_binz()
	       << " binphi " << cell->get_binphi()
	       << " edep " << cell->get_edep()
	       << " cell_id " << cell->get_cell_id()
	       << " g4hitid " << g4hitid
	       << endl;
	} 
      }  // end of verbose output

  }  // end of loop over layers

  _timer.get()->stop();
  return Fun4AllReturnCodes::EVENT_OK;
}
