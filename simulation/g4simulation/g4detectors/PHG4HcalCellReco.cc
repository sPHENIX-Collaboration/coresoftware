#include "PHG4HcalCellReco.h"
#include "PHG4CylinderGeomContainer.h"
#include "PHG4CylinderGeom.h"
#include "PHG4CylinderCellGeomContainer.h"
#include "PHG4CylinderCellGeom.h"
#include "PHG4CylinderCellv1.h"
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

PHG4HcalCellReco::PHG4HcalCellReco(const string &name) :
  SubsysReco(name),
  _timer(PHTimeServer::get()->insert_new(name.c_str())),
  nslatscombined(1),
  netabins(24), // that is our default
  chkenergyconservation(0),
  timing_min(0.0),
  timing_max(numeric_limits<double>::max())
{
  memset(nbins, 0, sizeof(nbins));  
}

int PHG4HcalCellReco::InitRun(PHCompositeNode *topNode)
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
          layerseggeo->set_binning(PHG4CylinderCellDefs::etaslatbinning);
          layerseggeo->set_etabins(netabins);
          layerseggeo->set_etamin(-1.1);
          layerseggeo->set_etastep((1.1+1.1)/netabins);
          layerseggeo->set_phimin(layergeom->get_phi_slat_zero());
	  double deltaphi = 2.*M_PI/layergeom->get_nscint();
	  layerseggeo->set_phistep(deltaphi*nslatscombined);
	  layerseggeo->set_phibins(nslatbins);
          etastep[layer] = (1.1+1.1)/netabins;
      // add geo object filled by different binning methods
      seggeo->AddLayerCellGeom(layerseggeo);
      if (verbosity > 1)
	{
	  layerseggeo->identify();
	}
    }
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
  PHG4CylinderCellContainer *cells = findNode::getClass<PHG4CylinderCellContainer>(topNode, cellnodename);
  if (! cells)
    {
      cout << "could not locate cell node " << cellnodename << endl;
      exit(1);
    }

  PHG4CylinderCellGeomContainer *seggeo = findNode::getClass<PHG4CylinderCellGeomContainer>(topNode , seggeonodename.c_str());
  if (! seggeo)
    {
      cout << "could not locate geo node " << seggeonodename << endl;
      exit(1);
    }

  PHG4HitContainer::LayerIter layer;
  pair<PHG4HitContainer::LayerIter, PHG4HitContainer::LayerIter> layer_begin_end = g4hit->getLayers();
  for (layer = layer_begin_end.first; layer != layer_begin_end.second; ++layer)
    {
      PHG4HitContainer::ConstIterator hiter;
      PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits(*layer);
      PHG4CylinderCellGeom *geo = seggeo->GetLayerCellGeom(*layer);
      int nslatbins = geo->get_phibins();
      for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
	{
	  // checking ADC timing integration window cut
	  if (hiter->second->get_t(0)>timing_max) continue;
	  if (hiter->second->get_t(1)<timing_min) continue;

	  int slatno = hiter->second->get_scint_id();
	  int slatbin;
	  slatbin = hiter->second->get_layer() / nslatscombined;
	  if (slatbin < 0 || slatbin > nslatbins)
	    {
	      if (slatbin + 1 > nslatbins)
		{
		  if (verbosity > 0)
		    {
		      cout << "dealing with non fitting slat binning, this one will be dropped" << endl;
		    }
		  continue;
		}
	      cout << "slatbin out of range: " << slatbin
		   << ", slatbins 0 - " << nslatbins
		   << ", slat no: " << hiter->second->get_scint_id()
		   << endl;
	    }

	  unsigned int key = (slatbin<<16) + slatno;
	  if (celllist.find(key) == celllist.end())
	    {
	      celllist[key] = new PHG4CylinderCellv1();
	      celllist[key]->set_layer(*layer);
	      celllist[key]->set_phibin(slatbin);
	      celllist[key]->set_etabin(slatno);
	    }

    celllist[key]->add_edep(hiter->first, hiter->second->get_edep(),hiter->second->get_light_yield());
    celllist[key]->add_shower_edep(hiter->second->get_shower_id(), hiter->second->get_edep());
	} // end loop over g4hits
      int numcells = 0;
      for (map<unsigned int, PHG4CylinderCell *>::const_iterator mapiter = celllist.begin();mapiter != celllist.end() ; ++mapiter)
	{
	  cells->AddCylinderCell(*layer, mapiter->second);
	  numcells++;
	  if (verbosity > 1)
	    {
	      cout << "Adding cell in bin slat: " << (mapiter->first >> 16)
		   << ", eta: " << (mapiter->first & 0xFFFF)
		   << ", energy dep: " << mapiter->second->get_edep()
		   << endl;
	    }
	}
      celllist.clear();
      if (verbosity > 0)
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
PHG4HcalCellReco::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void
PHG4HcalCellReco::etasize_nslat(const int i, const double deltaeta, const int nslat)
{
  if (nslat >= 1)
    {
      nslatscombined = nslat;
      set_size(i, deltaeta, nslat, PHG4CylinderCellDefs::etaslatbinning);
    }
  else
    {
      cout << PHWHERE << ": bad number of slats to combine: " << nslat << endl;
      exit(1);
    }
  return;
}

void
PHG4HcalCellReco::set_size(const int i, const double sizeA, const int sizeB, const int what)
{
  if (binning.find(i) != binning.end())
    {
      cout << "size for layer " << i << " already set" << endl;
      return;
    }
  binning[i] = what;
  cell_size[i] = std::make_pair(sizeA, sizeB);
  return;
}

int
PHG4HcalCellReco::CheckEnergy(PHCompositeNode *topNode)
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


double
PHG4HcalCellReco::get_phi_slat_zero_low(const double radius, const double thickness, const double tiltangle)
{
  // A/sin(alpha) = C/sin(gamma)
  // beta = 90-gamma

  double sinalpha = ((radius+thickness/2.)/radius) * sin(tiltangle);
  double beta =  asin(sinalpha) - tiltangle;
  cout << "beta: " << beta * 180./M_PI << endl;
  return beta;
}

double
PHG4HcalCellReco::get_phi_slat_zero_up(const double radius, const double thickness, const double tiltangle)
{
  // A/sin(alpha) = C/sin(gamma)
  // beta = 90-gamma
  double a = radius+thickness;
  double c = radius+thickness/2.;
  double alpha = M_PI - tiltangle;
  double singamma = c/a * sin(alpha);
  double beta =  M_PI - asin(singamma) - alpha;
  cout << "beta: " << beta * 180./M_PI << endl;
  return beta;
}
