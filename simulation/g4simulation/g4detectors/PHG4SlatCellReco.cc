#include "PHG4SlatCellReco.h"
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


#include <boost/tuple/tuple.hpp>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <limits>       // std::numeric_limits

using namespace std;

#define ARRAYDIM1 1000
#define ARRAYDIM2 1000
static PHG4CylinderCell *cellptarray[ARRAYDIM1][ARRAYDIM2];

PHG4SlatCellReco::PHG4SlatCellReco(const string &name) :
  SubsysReco(name),
  _timer(PHTimeServer::get()->insert_new(name.c_str())),
  nslatscombined(1),
  chkenergyconservation(0),
  tmin_default(0.0),  // ns
  tmax_default(60.0), // ns
  tmin_max()
{
  memset(nbins, 0, sizeof(nbins));
  memset(cellptarray, 0, sizeof(cellptarray));
}

int PHG4SlatCellReco::InitRun(PHCompositeNode *topNode)
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
      //     double circumference = layergeom->get_radius() * 2 * M_PI;
      //     double length_in_z = layergeom->get_zmax() - layergeom->get_zmin();
      sizeiter = cell_size.find(layer);
      if (sizeiter == cell_size.end())
        {
          cout << Name() << ": no cell sizes for layer " << layer << endl;
          exit(1);
        }
      // create geo object and fill with variables common to all binning methods
      PHG4CylinderCellGeom *layerseggeo = new PHG4CylinderCellGeom();
      layerseggeo->set_layer(layergeom->get_layer());
      layerseggeo->set_radius(layergeom->get_radius());
      layerseggeo->set_thickness(layergeom->get_thickness());
      if (binning[layer] == PHG4CylinderCellDefs::etaslatbinning)
        {
          // calculate eta at radius+ thickness (outer radius)
          // length via eta coverage is calculated using the outer radius
          double etamin = get_eta(layergeom->get_radius() + layergeom->get_thickness(), layergeom->get_zmin());
          double etamax = get_eta(layergeom->get_radius() + layergeom->get_thickness(), layergeom->get_zmax());
          zmin_max[layer] = make_pair(etamin, etamax);
          double etastepsize = (sizeiter->second).first;
          double d_etabins;
          double fract = modf((etamax - etamin) / etastepsize, &d_etabins);
          if (fract != 0)
            {
              d_etabins++;
            }
          etastepsize = (etamax - etamin) / d_etabins;
          (sizeiter->second).first = etastepsize;
          int etabins = d_etabins;
          double etahi = etamin + etastepsize;
          for (int i = 0; i < etabins; i++)
            {
              if (etahi > (etamax + 1.e-6)) // etahi is a tiny bit larger due to numerical uncertainties
                {
                  cout << "etahi: " << etahi << ", etamax: " << etamax << endl;
                }
              etahi +=  etastepsize;
            }
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
          n_phi_z_bins[layer] = make_pair(nslatbins, etabins);
	  if (verbosity > 1)
	    {
	      layergeom->identify();
	    }
          layerseggeo->set_binning(PHG4CylinderCellDefs::etaslatbinning);
          layerseggeo->set_etabins(etabins);
          layerseggeo->set_etamin(etamin);
          layerseggeo->set_etastep(etastepsize);
          layerseggeo->set_phimin(layergeom->get_phi_slat_zero());
	  double deltaphi = 2.*M_PI/layergeom->get_nscint();
	  layerseggeo->set_phistep(deltaphi*nslatscombined);
	  layerseggeo->set_phibins(nslatbins);
          etastep[layer] = etastepsize;
        }
      // add geo object filled by different binning methods
      seggeo->AddLayerCellGeom(layerseggeo);
      if (verbosity > 1)
	{
	  layerseggeo->identify();
	}
    }

  for (std::map<int,int>::iterator iter = binning.begin(); 
       iter != binning.end(); ++iter) {
    int layer = iter->first;
    // if the user doesn't set an integration window, set the default
    tmin_max.insert(std::make_pair(layer,std::make_pair(tmin_default,tmax_default)));    
  }
  
  return Fun4AllReturnCodes::EVENT_OK;
}


int
PHG4SlatCellReco::process_event(PHCompositeNode *topNode)
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
      int nslatbins = n_phi_z_bins[*layer].first;
      int nzbins = n_phi_z_bins[*layer].second;
      if ( nslatbins >= ARRAYDIM1)
        {
          cout << "required bins in numslats " << nslatbins
               << " exceed array size: " << ARRAYDIM1
               << " adjust ARRAYDIM1 and recompile" << endl;
          exit(1);
        }

      if (nzbins >= ARRAYDIM2)
        {
          cout << "required bins in z: " << nzbins
               << " exceed array size: " << ARRAYDIM2
               << " adjust ARRAYDIM2 and recompile" << endl;
          exit(1);
        }


      // ------- eta/phi binning ------------------------------------------------------------------------
      if (binning[*layer] == PHG4CylinderCellDefs::etaslatbinning)
        {
          for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
            {
              // checking ADC timing integration window cut
              if (hiter->second->get_t(0)>tmin_max[*layer].second) continue;
	      if (hiter->second->get_t(1)<tmin_max[*layer].first) continue;

              double etaphi[2];
              int slatbin;
              double etabin[2];
              for (int i = 0; i < 2; i++)
                {
                  etaphi[i] = get_eta(sqrt(hiter->second->get_x(i) * hiter->second->get_x(i) + hiter->second->get_y(i) * hiter->second->get_y(i)), hiter->second->get_z(i));
                  etabin[i] = geo->get_etabin( etaphi[i]);
                }
              slatbin = hiter->second->get_scint_id() / nslatscombined;
              if (etabin[0] < 0 || etabin[0] >= nzbins   || etabin[1] < 0 || etabin[1] >= nzbins)
                {
		  if (verbosity > 1)
		    {
		      cout << "eta out of range: etabin[0]: " << etabin[0]
			   << ", etabin[1]: " << etabin[1]
			   << ", nzbins: 0 - " << nzbins
			   << ", etaphi[0]: " << etaphi[0]
			   << ", etaphi[1]: " << etaphi[1]
			   << endl;
		    }
		  continue;
		}

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
	      if (etabin[0] < 0)

                {
                  if (verbosity > 0)
                    {
                      hiter->second->identify();
                    }
                  continue;
                }

              int intetabin = etabin[0];
	      int intetabinout = etabin[1];
	      if (intetabin != intetabinout)
		{
		  vector< boost::tuple<int, int, double> > bins_fraction;
		  double etaboundary = NAN;
		  double boundlow = etaphi[0];
		  double boundhi = etaphi[1];
		  int binstart = intetabin;
		  int binend = intetabinout;
		  if (intetabin > intetabinout)
		    {
		      binstart = intetabinout;
		      binend = intetabin;
		      boundlow = etaphi[1];
		      boundhi = etaphi[0];
		    }
		  for (int i = binstart; i <= binend; i++)
		    {
		      etaboundary = geo->get_etamin() +  geo->get_etastep() * (i + 1);
		      double fraction = fabs((boundlow - etaboundary) / (etaphi[0] - etaphi[1]));
		      if (etaboundary > boundhi)
			{
			  fraction = fabs((boundlow - boundhi) / (etaphi[0] - etaphi[1]));
			}
		      bins_fraction.push_back(boost::make_tuple(slatbin, i, fraction));
		      boundlow = etaboundary;
		    }
		  while (bins_fraction.size() > 0)
		    {
		      slatbin = bins_fraction.back().get<0>();
		      intetabin = bins_fraction.back().get<1>();
		      if (!cellptarray[slatbin][intetabin])
			{
			  cellptarray[slatbin][intetabin] = new PHG4CylinderCellv1();
			  cellptarray[slatbin][intetabin]->set_layer(*layer);
			  cellptarray[slatbin][intetabin]->set_phibin(slatbin);
			  cellptarray[slatbin][intetabin]->set_etabin(intetabin);
			}
		      cellptarray[slatbin][intetabin]->add_edep(hiter->first, hiter->second->get_edep()*bins_fraction.back().get<2>(), hiter->second->get_light_yield()*bins_fraction.back().get<2>());
		      bins_fraction.pop_back();
		    }
		}
	      else // all in same cell this is easy (keep from previous code)
		{
		  if (!cellptarray[slatbin][intetabin])
		    {
		      cellptarray[slatbin][intetabin] = new PHG4CylinderCellv1();
		      cellptarray[slatbin][intetabin]->set_layer(*layer);
		      cellptarray[slatbin][intetabin]->set_phibin(slatbin);
		      cellptarray[slatbin][intetabin]->set_etabin(intetabin);
		    }
      cellptarray[slatbin][intetabin]->add_edep(hiter->first, hiter->second->get_edep(), hiter->second->get_light_yield());
      cellptarray[slatbin][intetabin]->add_shower_edep(hiter->second->get_shower_id(), hiter->second->get_edep());
		}
	    } // end loop over g4hits
          int numcells = 0;
          for (int islat = 0; islat < nslatbins; islat++)
            {
              for (int iz = 0; iz < nzbins; iz++)
                {
                  if (cellptarray[islat][iz])
                    {
                      cells->AddCylinderCell(*layer, cellptarray[islat][iz]);
                      numcells++;
                      if (verbosity > 1)
                        {
                          cout << "Adding cell in bin slat: " << islat
                               << ", eta: " << iz
                               << ", energy dep: " << cellptarray[islat][iz]->get_edep()
                               << endl;
                        }
                      cellptarray[islat][iz] = 0;
                    }

                }
            }

          if (verbosity > 0)
            {
              cout << Name() << ": found " << numcells << " eta/slat cells with energy deposition" << endl;
            }
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
PHG4SlatCellReco::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void
PHG4SlatCellReco::cellsize(const int i, const double sr, const double sz)
{
  set_size(i, sr, sz, PHG4CylinderCellDefs::sizebinning);
}

void
PHG4SlatCellReco::etaphisize(const int i, const double deltaeta, const double deltaphi)
{
  set_size(i, deltaeta, deltaphi, PHG4CylinderCellDefs::etaphibinning);
  return;
}

void
PHG4SlatCellReco::etasize_nslat(const int i, const double deltaeta, const int nslat)
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
PHG4SlatCellReco::set_size(const int i, const double sizeA, const int sizeB, const int what)
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

pair<double, double>
PHG4SlatCellReco::get_etaphi(const double x, const double y, const double z)
{
  double eta;
  double phi;
  double radius;
  double theta;
  radius = sqrt(x * x + y * y);
  phi = atan2(y, x);
  theta = atan2(radius, z);
  eta = -log(tan(theta / 2.));
  return make_pair(eta, phi);
}

double
PHG4SlatCellReco::get_eta(const double radius, const double z)
{
  double eta;
  double theta;
  theta = atan2(radius, fabs(z));
  eta = -log(tan(theta / 2.));
  if (z < 0)
    {
      eta = -eta;
    }
  return eta;
}

int
PHG4SlatCellReco::CheckEnergy(PHCompositeNode *topNode)
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
PHG4SlatCellReco::get_phi_slat_zero_low(const double radius, const double thickness, const double tiltangle)
{
  // A/sin(alpha) = C/sin(gamma)
  // beta = 90-gamma

  double sinalpha = ((radius+thickness/2.)/radius) * sin(tiltangle);
  double beta =  asin(sinalpha) - tiltangle;
  cout << "beta: " << beta * 180./M_PI << endl;
  return beta;
}

double
PHG4SlatCellReco::get_phi_slat_zero_up(const double radius, const double thickness, const double tiltangle)
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
