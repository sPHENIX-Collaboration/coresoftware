#include "PHG4HcalCellReco.h"
#include "PHG4ScintillatorSlatv1.h"
#include "PHG4ScintillatorSlatContainer.h"
#include "PHG4Parameters.h"

#include <g4main/PHG4Hit.h>
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
#include <limits>       // std::numeric_limits

using namespace std;
// for hcal dimension
#define ROWDIM 320
#define COLUMNDIM 24

static PHG4ScintillatorSlat *slatarray[ROWDIM][COLUMNDIM];

PHG4HcalCellReco::PHG4HcalCellReco(const string &name) :
  SubsysReco(name),
  PHG4ParameterInterface(name),
  _timer(PHTimeServer::get()->insert_new(name.c_str())),
  chkenergyconservation(0),
  tmin(NAN),  // ns
  tmax(NAN) // ns
{
  memset(slatarray, 0, sizeof(slatarray));
  InitializeParameters();
}

int
PHG4HcalCellReco::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
    {
      std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
      exit(1);
    }
  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN" ));
  PHCompositeNode *parNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "PAR" ));

  string paramnodename = "G4CELLPARAM_" + detector;
  string geonodename = "G4CELLGEO_" + detector;
  hitnodename = "G4HIT_" + detector;
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  if (!g4hit)
    {
      cout << Name() << " Could not locate G4HIT node " << hitnodename << endl;
      gSystem->Exit(1);
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
  UpdateParametersWithMacro();
  // save this to the run wise tree to store on DST
  PHNodeIterator runIter(runNode);
  PHCompositeNode *RunDetNode =  dynamic_cast<PHCompositeNode*>(runIter.findFirst("PHCompositeNode",detector));
  if (! RunDetNode)
    {
      RunDetNode = new PHCompositeNode(detector);
      runNode->addNode(RunDetNode);
    }
  SaveToNodeTree(RunDetNode,paramnodename);
  // save this to the parNode for use
  PHNodeIterator parIter(parNode);
  PHCompositeNode *ParDetNode =  dynamic_cast<PHCompositeNode*>(parIter.findFirst("PHCompositeNode",detector));
  if (! ParDetNode)
    {
      ParDetNode = new PHCompositeNode(detector);
      parNode->addNode(ParDetNode);
    }
  PutOnParNode(ParDetNode,geonodename);
  tmin = get_double_param("tmin");
  tmax = get_double_param("tmax");
  return Fun4AllReturnCodes::EVENT_OK;
}


int
PHG4HcalCellReco::process_event(PHCompositeNode *topNode)
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

  PHG4HitContainer::ConstIterator hiter;
  PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits();
  for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
    {
      if (hiter->second->get_t(0) > tmax) continue;
      if (hiter->second->get_t(1) < tmin) continue;
      short icolumn = hiter->second->get_scint_id();
      int introw = (hiter->second->get_hit_id() >> PHG4HitDefs::hit_idbits);
      if ( introw >= ROWDIM || introw < 0)
	{
	  cout << "row " << introw
	       << " exceed array size: " << ROWDIM
	       << " adjust ROWDIM and recompile" << endl;
	  exit(1);
	}
      // after checking for size of introw so we do not run into
      // overflow issues, put this into the short we want later
      short irow = introw;
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
      // cout << "row: " << irow
      //  	   << ", column: " << hiter->second->get_scint_id() << endl;
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

  if (chkenergyconservation)
    {
      CheckEnergy(topNode);
    }
  _timer.get()->stop();
  return Fun4AllReturnCodes::EVENT_OK;
}

int
PHG4HcalCellReco::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}


int
PHG4HcalCellReco::CheckEnergy(PHCompositeNode *topNode)
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
      cout << "hint: timing cuts tmin/tmax will do this to you" << endl;
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

void
PHG4HcalCellReco::SetDefaultParameters()
{
  set_default_double_param("tmax",60.0);
  set_default_double_param("tmin",0.0);
  return;
}

void
PHG4HcalCellReco::set_timing_window(const double tmi, const double tma)
{
  set_double_param("tmin",tmi);
  set_double_param("tmax",tma);
}
