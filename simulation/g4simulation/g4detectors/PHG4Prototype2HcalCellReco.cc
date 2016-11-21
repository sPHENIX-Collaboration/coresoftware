#include "PHG4Prototype2HcalCellReco.h"
#include "PHG4ScintillatorSlatv1.h"
#include "PHG4ScintillatorSlatContainer.h"

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
// for hcal dimension
#define ROWDIM 320
#define COLUMNDIM 22
static PHG4ScintillatorSlat *slatarray[ROWDIM][COLUMNDIM];

PHG4Prototype2HcalCellReco::PHG4Prototype2HcalCellReco(const string &name) :
  SubsysReco(name),
  _timer(PHTimeServer::get()->insert_new(name.c_str())),
  chkenergyconservation(0),
  tmin(0.0),  // ns
  tmax(60.0) // ns
{
  memset(slatarray, 0, sizeof(slatarray));
}

int PHG4Prototype2HcalCellReco::InitRun(PHCompositeNode *topNode)
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
  PHG4ScintillatorSlatContainer *slats = findNode::getClass<PHG4ScintillatorSlatContainer>(topNode , cellnodename);
  if (!slats)
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
      slats = new PHG4ScintillatorSlatContainer();
      PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(slats, cellnodename.c_str() , "PHObject");
      DetNode->addNode(newNode);
    }

  return Fun4AllReturnCodes::EVENT_OK;
}


int
PHG4Prototype2HcalCellReco::process_event(PHCompositeNode *topNode)
{
  _timer.get()->restart();
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  if (!g4hit)
    {
      cout << "Could not locate g4 hit node " << hitnodename << endl;
      exit(1);
    }
  PHG4ScintillatorSlatContainer *slats = findNode::getClass<PHG4ScintillatorSlatContainer>(topNode, cellnodename);
  if (! slats)
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
	  if (hiter->second->get_t(0)>tmax) continue;
	  if (hiter->second->get_t(1)<tmin) continue;
	  short icolumn = hiter->second->get_scint_id();
	  short irow = hiter->second->get_row();
	  if ( irow >= ROWDIM || irow < 0)
	    {
	      cout << "row " << irow
		   << " exceed array size: " << ROWDIM
		   << " adjust ROWDIM and recompile" << endl;
	      exit(1);
	    }

	  if (icolumn >= COLUMNDIM || icolumn < 0)
	    {
	      cout << "column: " << icolumn
		   << " exceed array size: " << COLUMNDIM
		   << " adjust COLUMNDIM and recompile" << endl;
	      exit(1);
	    }


	  if (!slatarray[irow][icolumn])
	    {
	      slatarray[irow][icolumn] = new PHG4ScintillatorSlatv1();
	    }
	  slatarray[irow][icolumn]->add_edep(hiter->second->get_edep(),
					     hiter->second->get_eion(),
					     hiter->second->get_light_yield());
	  slatarray[irow][icolumn]->add_hit_key(hiter->first);
	  // cout << "row: " << hiter->second->get_row() 
	  // 	   << ", column: " << hiter->second->get_scint_id() << endl;
	  // checking ADC timing integration window cut
	} // end loop over g4hits
      int nslathits = 0;
      for (int irow = 0; irow<ROWDIM; irow++)
	{
	  for (int icolumn = 0; icolumn<COLUMNDIM; icolumn++)
	    {
	      if (slatarray[irow][icolumn])
		{
		  PHG4ScintillatorSlatDefs::keytype key = PHG4ScintillatorSlatDefs::genkey(irow,icolumn);
		  slats->AddScintillatorSlat(key,slatarray[irow][icolumn]);
		  slatarray[irow][icolumn] = NULL;
		  nslathits++;
		}
	    }
	}
      if (verbosity > 0)
	{
	  cout << Name() << ": found " << nslathits << " slats with energy deposition" << endl;
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
PHG4Prototype2HcalCellReco::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}


int
PHG4Prototype2HcalCellReco::CheckEnergy(PHCompositeNode *topNode)
{
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  PHG4ScintillatorSlatContainer *slats = findNode::getClass<PHG4ScintillatorSlatContainer>(topNode, cellnodename);
  double sum_energy_g4hit = 0.;
  double sum_energy_cells = 0.;
  PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits();
  PHG4HitContainer::ConstIterator hiter;
  for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
    {
      sum_energy_g4hit += hiter->second->get_edep();
    }
  PHG4ScintillatorSlatContainer::ConstRange cell_begin_end = slats->getScintillatorSlats();
  PHG4ScintillatorSlatContainer::ConstIterator citer;
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
