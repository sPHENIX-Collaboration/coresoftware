#include "PHG4ForwardCalCellReco.h"
#include "PHG4CylinderCellv3.h"
#include "PHG4CylinderCellContainer.h"
#include "PHG4CylinderCellDefs.h"

#include <g4main/PHG4Hit.h>
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
#include <limits>       // std::numeric_limits

using namespace std;

PHG4ForwardCalCellReco::PHG4ForwardCalCellReco(const string &name) :
  SubsysReco(name),
  _timer(PHTimeServer::get()->insert_new(name.c_str())),
  chkenergyconservation(0),
  tmin_default(0.0),  // ns
  tmax_default(60.0), // ns
  tmin_max()
{
  memset(nbins, 0, sizeof(nbins));  
}

int PHG4ForwardCalCellReco::InitRun(PHCompositeNode *topNode)
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
  
  return Fun4AllReturnCodes::EVENT_OK;
}


int
PHG4ForwardCalCellReco::process_event(PHCompositeNode *topNode)
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

  PHG4HitContainer::LayerIter layer;
  pair<PHG4HitContainer::LayerIter, PHG4HitContainer::LayerIter> layer_begin_end = g4hit->getLayers();
  for (layer = layer_begin_end.first; layer != layer_begin_end.second; ++layer)
    {
      PHG4HitContainer::ConstIterator hiter;
      PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits(*layer);
      for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
	{
	  // checking ADC timing integration window cut
	  if (hiter->second->get_t(0)>tmax_default) continue;
	  if (hiter->second->get_t(1)<tmin_default) continue;

          // only hits that deposited energy (or geantinos)
          if (hiter->second->get_edep()<=0)
            continue;

	  unsigned int key = (hiter->second->get_index_j()<<16) + hiter->second->get_index_k();
	  if (celllist.find(key) == celllist.end())
	    {
	      celllist[key] = new PHG4CylinderCellv3();
	      celllist[key]->set_layer(*layer);
	      celllist[key]->set_j_index(hiter->second->get_index_j());
	      celllist[key]->set_k_index(hiter->second->get_index_k());
	      celllist[key]->set_l_index(hiter->second->get_index_l());
	    }

	  celllist[key]->add_edep(hiter->first, hiter->second->get_edep(), hiter->second->get_light_yield());
	  celllist[key]->add_shower_edep(hiter->second->get_shower_id(), hiter->second->get_edep());

	}

      int numcells = 0;
      for (map<unsigned int, PHG4CylinderCell *>::const_iterator mapiter = celllist.begin();mapiter != celllist.end() ; ++mapiter)
	{
	  cells->AddCylinderCellSpecifyKey(mapiter->first, mapiter->second);
	  numcells++;
	}
      celllist.clear();
      if (Verbosity() > 0)
	{
	  cout << Name() << ": found " << numcells << " eta/slat cells with energy deposition" << endl;
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
PHG4ForwardCalCellReco::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}


int
PHG4ForwardCalCellReco::CheckEnergy(PHCompositeNode *topNode)
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
      if (Verbosity() > 0)
	{
	  cout << Name() << ": total energy for this event: " << sum_energy_g4hit << " GeV" << endl;
	}
    }
  return 0;
}


