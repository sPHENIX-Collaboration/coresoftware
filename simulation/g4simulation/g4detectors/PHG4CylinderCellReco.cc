#include "PHG4CylinderCellReco.h"
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

#include <TROOT.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <limits>       // std::numeric_limits


using namespace std;

PHG4CylinderCellReco::PHG4CylinderCellReco(const string &name) :
  SubsysReco(name),
  _timer(PHTimeServer::get()->insert_new("PHG4CylinderCellReco")),
  chkenergyconservation(0),
  tmin_default(-0.0),  // ns
  tmax_default(100.0), // ns
  tmin_max()
{
  memset(nbins, 0, sizeof(nbins));
}

int PHG4CylinderCellReco::InitRun(PHCompositeNode *topNode)
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
  cellnodename = "G4CELL_" + outdetector;
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
  seggeonodename = "CYLINDERCELLGEOM_" + outdetector;
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
  map<int, std::pair <double, double> >::iterator sizeiter;
  for (miter = begin_end.first; miter != begin_end.second; ++miter)
    {
      PHG4CylinderGeom *layergeom = miter->second;
      int layer = layergeom->get_layer();
      double circumference = layergeom->get_radius() * 2 * M_PI;
      double length_in_z = layergeom->get_zmax() - layergeom->get_zmin();
      sizeiter = cell_size.find(layer);
      if (sizeiter == cell_size.end())
        {
          cout << "no cell sizes for layer " << layer << endl;
          exit(1);
        }
      // create geo object and fill with variables common to all binning methods
      PHG4CylinderCellGeom *layerseggeo = new PHG4CylinderCellGeom();
      layerseggeo->set_layer(layergeom->get_layer());
      layerseggeo->set_radius(layergeom->get_radius());
      layerseggeo->set_thickness(layergeom->get_thickness());
      if (binning[layer] == PHG4CylinderCellDefs::etaphibinning)
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
          double etalow = etamin;
          double etahi = etalow + etastepsize;
          for (int i = 0; i < etabins; i++)
            {
              if (etahi > (etamax + 1.e-6)) // etahi is a tiny bit larger due to numerical uncertainties
                {
                  cout << "etahi: " << etahi << ", etamax: " << etamax << endl;
                }
              etahi +=  etastepsize;
            }

          double phimin = -M_PI;
          double phimax = M_PI;
          double phistepsize = (sizeiter->second).second;
          double d_phibins;
          fract = modf((phimax - phimin) / phistepsize, &d_phibins);
          if (fract != 0)
            {
              d_phibins++;
            }
          phistepsize = (phimax - phimin) / d_phibins;
          (sizeiter->second).second = phistepsize;
          int phibins = d_phibins;
          double philow = phimin;
          double phihi = philow + phistepsize;
          for (int i = 0; i < phibins; i++)
            {
              if (phihi > phimax)
                {
                  cout << "phihi: " << phihi << ", phimax: " << phimax << endl;
                }
              phihi +=  phistepsize;
            }
          pair<int, int> phi_z_bin = make_pair(phibins, etabins);
          n_phi_z_bins[layer] = phi_z_bin;
          layerseggeo->set_binning(PHG4CylinderCellDefs::etaphibinning);
          layerseggeo->set_etabins(etabins);
          layerseggeo->set_etamin(etamin);
          layerseggeo->set_etastep(etastepsize);
          layerseggeo->set_phimin(phimin);
          layerseggeo->set_phibins(phibins);
          layerseggeo->set_phistep(phistepsize);
          phistep[layer] = phistepsize;
          etastep[layer] = etastepsize;
        }
      else if (binning[layer] == PHG4CylinderCellDefs::sizebinning)
        {
          zmin_max[layer] = make_pair(layergeom->get_zmin(), layergeom->get_zmax());
          double size_z = (sizeiter->second).second;
          double size_r = (sizeiter->second).first;
          double bins_r;
          // unlikely but if the circumference is a multiple of the cell size
          // use result of division, if not - add 1 bin which makes the
          // cells a tiny bit smaller but makes them fit
          double fract = modf(circumference / size_r, &bins_r);
          if (fract != 0)
            {
              bins_r++;
            }
          nbins[0] = bins_r;
          size_r = circumference / bins_r;
          (sizeiter->second).first = size_r;
          double phistepsize = 2 * M_PI / bins_r;
          double phimin = -M_PI;
          double phimax = phimin + phistepsize;
          phistep[layer] = phistepsize;
          for (int i = 0 ; i < nbins[0]; i++)
            {
              if (phimax > (M_PI + 1e-9))
                {
                  cout << "phimax: " << phimax << ", M_PI: " << M_PI
		       << "phimax-M_PI: " << phimax-M_PI << endl;
                }
              phimax += phistepsize;
            }
          // unlikely but if the length is a multiple of the cell size
          // use result of division, if not - add 1 bin which makes the
          // cells a tiny bit smaller but makes them fit
          fract = modf(length_in_z / size_z, &bins_r);
          if (fract != 0)
            {
              bins_r++;
            }
          nbins[1] = bins_r;
          pair<int, int> phi_z_bin = make_pair(nbins[0], nbins[1]);
          n_phi_z_bins[layer] = phi_z_bin;
          // update our map with the new sizes
          size_z = length_in_z / bins_r;
          (sizeiter->second).second = size_z;
          double zlow = layergeom->get_zmin();
          double zhigh = zlow + size_z;;
          for (int i = 0 ; i < nbins[1]; i++)
            {
              if (zhigh > (layergeom->get_zmax()+1e-9))
                {
                  cout << "zhigh: " << zhigh << ", zmax " 
                       << layergeom->get_zmax()
		       << ", zhigh-zmax: " <<  zhigh-layergeom->get_zmax()
                       << endl;
                }
              zhigh += size_z;
            }
          layerseggeo->set_binning(PHG4CylinderCellDefs::sizebinning);
          layerseggeo->set_zbins(nbins[1]);
          layerseggeo->set_zmin(layergeom->get_zmin());
          layerseggeo->set_zstep(size_z);
          layerseggeo->set_phibins(nbins[0]);
          layerseggeo->set_phistep(phistepsize);
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
  
  // print out settings
  if (verbosity > 0) {
    cout << "===================== PHG4CylinderCellReco::InitRun() =====================" << endl;
    cout << " " << outdetector << " Segmentation Description: " << endl;
    for (std::map<int,int>::iterator iter = binning.begin(); 
	 iter != binning.end(); ++iter) {
      int layer = iter->first;

      if (binning[layer] == PHG4CylinderCellDefs::etaphibinning) {
	// phi & eta bin is usually used to make projective towers
	// so just print the first layer
	cout << " Layer #" << binning.begin()->first << "-" << binning.rbegin()->first << endl;
	cout << "   Nbins (phi,eta): (" << n_phi_z_bins[layer].first << ", " << n_phi_z_bins[layer].second << ")" << endl;
	cout << "   Cell Size (phi,eta): (" << cell_size[layer].first << " rad, " << cell_size[layer].second << " units)" << endl;
	break;
      } else if (binning[layer] == PHG4CylinderCellDefs::sizebinning) {	
	cout << " Layer #" << layer << endl;
	cout << "   Nbins (phi,z): (" << n_phi_z_bins[layer].first << ", " << n_phi_z_bins[layer].second << ")" << endl;
	cout << "   Cell Size (phi,z): (" << cell_size[layer].first << " cm, " << cell_size[layer].second << " cm)" << endl;
      }
    }
    cout << "===========================================================================" << endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}


int
PHG4CylinderCellReco::process_event(PHCompositeNode *topNode)
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

  map<int, std::pair <double, double> >::iterator sizeiter;
  PHG4HitContainer::LayerIter layer;
  pair<PHG4HitContainer::LayerIter, PHG4HitContainer::LayerIter> layer_begin_end = g4hit->getLayers();
  //   cout << "number of layers: " << g4hit->num_layers() << endl;
  //   cout << "number of hits: " << g4hit->size() << endl;
  //   for (layer = layer_begin_end.first; layer != layer_begin_end.second; layer++)
  //     {
  //       cout << "layer number: " << *layer << endl;
  //     }
  for (layer = layer_begin_end.first; layer != layer_begin_end.second; layer++)
    {
      PHG4HitContainer::ConstIterator hiter;
      PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits(*layer);
      PHG4CylinderCellGeom *geo = seggeo->GetLayerCellGeom(*layer);
      int nphibins = n_phi_z_bins[*layer].first;
      int nzbins = n_phi_z_bins[*layer].second;

      // ------- eta/phi binning ------------------------------------------------------------------------
      if (binning[*layer] == PHG4CylinderCellDefs::etaphibinning)
        {
          for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; hiter++)
            {
              // checking ADC timing integration window cut
              if (hiter->second->get_t(0)>tmin_max[*layer].second) continue;
	      if (hiter->second->get_t(1)<tmin_max[*layer].first) continue;
	      
              pair<double, double> etaphi[2];
              double phibin[2];
              double etabin[2];
              for (int i = 0; i < 2; i++)
                {
                  etaphi[i] = get_etaphi(hiter->second->get_x(i), hiter->second->get_y(i), hiter->second->get_z(i));
                  etabin[i] = geo->get_etabin( etaphi[i].first );
                  phibin[i] = geo->get_phibin( etaphi[i].second );
                }
              // check bin range
              if (phibin[0] < 0 || phibin[0] >= nphibins || phibin[1] < 0 || phibin[1] >= nphibins)
                {
                  continue;
                }
              if (etabin[0] < 0 || etabin[0] >= nzbins   || etabin[1] < 0 || etabin[1] >= nzbins)
                {
                  continue;
                }

              if (etabin[0] < 0)
                {
                  if (verbosity > 0)
                    {
                      hiter->second->identify();
                    }
                  continue;
                }

              int intphibin = phibin[0];
              int intetabin = etabin[0];
              int intphibinout = phibin[1];
              int intetabinout = etabin[1];

              // Determine all fired cells

              double ax = (etaphi[0]).second; // phi
              double ay = (etaphi[0]).first;  // eta
              double bx = (etaphi[1]).second;
              double by = (etaphi[1]).first;
              if (intphibin > intphibinout)
                {
                  int tmp = intphibin;
                  intphibin = intphibinout;
                  intphibinout = tmp;
                }
              if (intetabin > intetabinout)
                {
                  int tmp = intetabin;
                  intetabin = intetabinout;
                  intetabinout = tmp;
                }

              double trklen = sqrt((ax - bx) * (ax - bx) + (ay - by) * (ay - by));
	      // if entry and exit hit are the same (seems to happen rarely), trklen = 0
	      // which leads to a 0/0 and an NaN in edep later on
	      // this code does for particles in the same cell a trklen/trklen (vdedx[ii]/trklen)
	      // so setting this to any non zero number will do just fine
	      // I just pick -1 here to flag those strange hits in case I want t oanalyze them
	      // later on
	      if (trklen == 0)
		{
		  trklen = -1.;
		}
              vector<int> vphi;
              vector<int> veta;
              vector<double> vdedx;

              if (intphibin == intphibinout && intetabin == intetabinout)   // single cell fired
                {
                  if (verbosity > 0) cout << "SINGLE CELL FIRED: " << intphibin << " " << intetabin << endl;
                  vphi.push_back(intphibin);
                  veta.push_back(intetabin);
                  vdedx.push_back(trklen);
                }
              else
                {
                  for (int ibp = intphibin; ibp <= intphibinout; ibp++)
                    {
                      for (int ibz = intetabin; ibz <= intetabinout; ibz++)
                        {
                          double cx = geo->get_phicenter(ibp) - geo->get_phistep() / 2.;
                          double dx = geo->get_phicenter(ibp) + geo->get_phistep() / 2.;
                          double cy = geo->get_etacenter(ibz) - geo->get_etastep() / 2.;
                          double dy = geo->get_etacenter(ibz) + geo->get_etastep() / 2.;
                          double rr = 0.;
                          //cout << "##### line: " << ax << " " << ay << " " << bx << " " << by << endl;
                          //cout << "####### cell: " << cx << " " << cy << " " << dx << " " << dy << endl;
                          bool yesno = line_and_rectangle_intersect(ax, ay, bx, by, cx, cy, dx, dy, &rr);
                          if (yesno)
                            {
                              if (verbosity > 0) cout << "CELL FIRED: " << ibp << " " << ibz << " " << rr << endl;
                              vphi.push_back(ibp);
                              veta.push_back(ibz);
                              vdedx.push_back(rr);
                            }
                        }
                    }
                }
              if (verbosity > 0) cout << "NUMBER OF FIRED CELLS = " << vphi.size() << endl;

              double tmpsum = 0.;
              for (unsigned int ii = 0; ii < vphi.size(); ii++)
                {
                  tmpsum += vdedx[ii];
                  vdedx[ii] = vdedx[ii] / trklen;
                  if (verbosity > 0) cout << "  CELL " << ii << "  dE/dX = " <<  vdedx[ii] << endl;
                }
              if (verbosity > 0) cout << "    TOTAL TRACK LENGTH = " << tmpsum << " " << trklen << endl;

              for (unsigned int i1 = 0; i1 < vphi.size(); i1++)   // loop over all fired cells
                {

                  int iphibin = vphi[i1];
                  int ietabin = veta[i1];

		  // This is the key for cellptmap
		  // It is constructed using the phi and z (or eta) bin index values
		  // It will be unique for a given phi and z (or eta) bin combination
		  char inkey[1024];
		  sprintf(inkey,"%i-%i",iphibin,ietabin);
		  std::string key(inkey);

		  if(verbosity > 1)
		    cout << " iphibin " << iphibin << " ietabin " << ietabin << " key " << key << endl;

		  if(cellptmap.count(key) > 0)
		    {
		      if(verbosity > 1)
			cout << "  add energy to existing cell " << endl;

		      cellptmap.find(key)->second->add_edep(hiter->first, hiter->second->get_edep()*vdedx[i1], hiter->second->get_light_yield()*vdedx[i1]);
		      cellptmap.find(key)->second->add_shower_edep(hiter->second->get_shower_id(), hiter->second->get_edep()*vdedx[i1]);
		    }
		  else
		    {
		      if(verbosity > 1)
			cout << "    did not find a previous entry for key = " << key << " add a new one" << endl;

		      cellptmap[key] = new PHG4CylinderCellv1();
		      it = cellptmap.find(key);
                      it->second->set_layer(*layer);
                      it->second->set_phibin(iphibin);
                      it->second->set_etabin(ietabin);		      
		      it->second->add_edep(hiter->first, hiter->second->get_edep()*vdedx[i1], hiter->second->get_light_yield()*vdedx[i1]);
		      it->second->add_shower_edep(hiter->second->get_shower_id(), hiter->second->get_edep()*vdedx[i1]);
		    }

		  // just a sanity check - we don't want to mess up by having Nan's or Infs in our energy deposition
		  if (! isfinite(hiter->second->get_edep()*vdedx[i1]))
		    {
		      cout << PHWHERE << " invalid energy dep " << hiter->second->get_edep()
			   << " or path length: " << vdedx[i1] << endl;
		    }
                }
              vphi.clear();
              veta.clear();

            } // end loop over g4hits

          int numcells = 0;

	  for(it = cellptmap.begin(); it != cellptmap.end(); ++it)
	    {
	      cells->AddCylinderCell(*layer, it->second);
	      numcells++;
	      if (verbosity > 1)
		{
		  cout << "Adding cell in bin phi: " << it->second->get_binphi()
		       << " phi: " << geo->get_phicenter(it->second->get_binphi()) * 180./M_PI
		       << ", z bin: " << it->second->get_bineta()
		       << ", z: " <<  geo->get_etacenter(it->second->get_bineta())
		       << ", energy dep: " << it->second->get_edep()
		       << endl;
		}
	    }

          if (verbosity > 0)
            {
              cout << Name() << ": found " << numcells << " eta/phi cells with energy deposition" << endl;
            }
        }



      else // ------ size binning ---------------------------------------------------------------
        {
          sizeiter = cell_size.find(*layer);
          if (sizeiter == cell_size.end())
            {
              cout << "logical screwup!!! no sizes for layer " << *layer << endl;
              exit(1);
            }
          double zstepsize = (sizeiter->second).second;
          double phistepsize = phistep[*layer];

          for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; hiter++)
            {
              // checking ADC timing integration window cut
              if (hiter->second->get_t(0)>tmin_max[*layer].second) continue;
	      if (hiter->second->get_t(1)<tmin_max[*layer].first) continue;

              double xinout[2];
              double yinout[2];
              double px[2];
              double py[2];
              double phi[2];
              double z[2];
              double phibin[2];
              double zbin[2];
              if (verbosity > 0) cout << "--------- new hit in layer # " << *layer << endl;

              for (int i = 0; i < 2; i++)
                {
                  xinout[i] = hiter->second->get_x(i);
                  yinout[i] = hiter->second->get_y(i);
                  px[i] = hiter->second->get_px(i);
                  py[i] = hiter->second->get_py(i);
                  phi[i] = atan2(hiter->second->get_y(i), hiter->second->get_x(i));
                  z[i] =  hiter->second->get_z(i);
                  phibin[i] = geo->get_phibin( phi[i] );
                  zbin[i] = geo->get_zbin( hiter->second->get_z(i) );

                  if (verbosity > 0) cout << " " << i << "  phibin: " << phibin[i] << ", phi: " << phi[i] << ", stepsize: " << phistepsize << endl;
                  if (verbosity > 0) cout << " " << i << "  zbin: " << zbin[i] << ", z = " << hiter->second->get_z(i) << ", stepsize: " << zstepsize << " offset: " <<  zmin_max[*layer].first << endl;
                }
              // check bin range
              if (phibin[0] < 0 || phibin[0] >= nphibins || phibin[1] < 0 || phibin[1] >= nphibins)
                {
                  continue;
                }
              if (zbin[0] < 0 || zbin[0] >= nzbins   || zbin[1] < 0 || zbin[1] >= nzbins)
                {
                  continue;
                }


              if (zbin[0] < 0)
                {
                  hiter->second->identify();
                  continue;
                }

              int intphibin = phibin[0];
              int intzbin = zbin[0];
              int intphibinout = phibin[1];
              int intzbinout = zbin[1];

              if (verbosity > 0)
                {
                  cout << "    phi bin range: " << intphibin << " to " << intphibinout << " phi: " << phi[0] << " to " << phi[1] << endl;
                  cout << "    Z bin range: " << intzbin << " to " << intzbinout << " Z: " << z[0] << " to " << z[1] << endl;
                  cout << "    phi difference: " << (phi[1] - phi[0])*1000. << " milliradians." << endl;
                  cout << "    phi difference: " << 2.5*(phi[1] - phi[0])*10000. << " microns." << endl;
                  cout << "    path length = " << sqrt((xinout[1] - xinout[0])*(xinout[1] - xinout[0]) + (yinout[1] - yinout[0])*(yinout[1] - yinout[0])) << endl;
                  cout << "       px = " << px[0] << " " << px[1] << endl;
                  cout << "       py = " << py[0] << " " << py[1] << endl;
                  cout << "       x = " << xinout[0] << " " << xinout[1] << endl;
                  cout << "       y = " << yinout[0] << " " << yinout[1] << endl;
                }

              // Determine all fired cells

              double ax = phi[0];
              double ay = z[0];
              double bx = phi[1];
              double by = z[1];
              if (intphibin > intphibinout)
                {
                  int tmp = intphibin;
                  intphibin = intphibinout;
                  intphibinout = tmp;
                }
              if (intzbin > intzbinout)
                {
                  int tmp = intzbin;
                  intzbin = intzbinout;
                  intzbinout = tmp;
                }

              double trklen = sqrt((ax - bx) * (ax - bx) + (ay - by) * (ay - by));
              // if entry and exit hit are the same (seems to happen rarely), trklen = 0
              // which leads to a 0/0 and an NaN in edep later on
              // this code does for particles in the same cell a trklen/trklen (vdedx[ii]/trklen)
              // so setting this to any non zero number will do just fine
              // I just pick -1 here to flag those strange hits in case I want t oanalyze them
              // later on
              if (trklen == 0)
          {
            trklen = -1.;
          }
              vector<int> vphi;
              vector<int> vz;
              vector<double> vdedx;

              if (intphibin == intphibinout && intzbin == intzbinout)   // single cell fired
                {
                  if (verbosity > 0) cout << "SINGLE CELL FIRED: " << intphibin << " " << intzbin << endl;
                  vphi.push_back(intphibin);
                  vz.push_back(intzbin);
                  vdedx.push_back(trklen);
                }
              else
                {
                  for (int ibp = intphibin; ibp <= intphibinout; ibp++)
                    {
                      for (int ibz = intzbin; ibz <= intzbinout; ibz++)
                        {
                          double cx = geo->get_phicenter(ibp) - geo->get_phistep() / 2.;
                          double dx = geo->get_phicenter(ibp) + geo->get_phistep() / 2.;
                          double cy = geo->get_zcenter(ibz) - geo->get_zstep() / 2.;
                          double dy = geo->get_zcenter(ibz) + geo->get_zstep() / 2.;
                          double rr = 0.;
                          //cout << "##### line: " << ax << " " << ay << " " << bx << " " << by << endl;
                          //cout << "####### cell: " << cx << " " << cy << " " << dx << " " << dy << endl;
                          bool yesno = line_and_rectangle_intersect(ax, ay, bx, by, cx, cy, dx, dy, &rr);
                          if (yesno)
                            {
                              if (verbosity > 0) cout << "CELL FIRED: " << ibp << " " << ibz << " " << rr << endl;
                              vphi.push_back(ibp);
                              vz.push_back(ibz);
                              vdedx.push_back(rr);
                            }
                        }
                    }
                }
              if (verbosity > 0) cout << "NUMBER OF FIRED CELLS = " << vz.size() << endl;

              double tmpsum = 0.;
              for (unsigned int ii = 0; ii < vz.size(); ii++)
                {
                  tmpsum += vdedx[ii];
                  vdedx[ii] = vdedx[ii] / trklen;
                  if (verbosity > 0) cout << "  CELL " << ii << "  dE/dX = " <<  vdedx[ii] << endl;
                }
              if (verbosity > 0) cout << "    TOTAL TRACK LENGTH = " << tmpsum << " " << trklen << endl;

              for (unsigned int i1 = 0; i1 < vphi.size(); i1++)   // loop over all fired cells
                {
                  int iphibin = vphi[i1];
                  int izbin = vz[i1];

		  char inkey[1024];
		  sprintf(inkey,"%i-%i",iphibin,izbin);
		  std::string key(inkey);

		  if(verbosity > 1)
		    cout << " iphibin " << iphibin << " izbin " << izbin << " key " << key << endl;
	
		  // check to see if there is already an entry for this cell
		  if(cellptmap.count(key) > 0)
		    {
		      if(verbosity > 1)
			cout << "  add energy to existing cell for key = " << cellptmap.find(key)->first << endl;

		      cellptmap.find(key)->second->add_edep(hiter->first, hiter->second->get_edep()*vdedx[i1], hiter->second->get_light_yield()*vdedx[i1]);
		      cellptmap.find(key)->second->add_shower_edep(hiter->second->get_shower_id(), hiter->second->get_edep()*vdedx[i1]);

		      if(verbosity > 1 && std::isnan(hiter->second->get_light_yield()*vdedx[i1]))
            {

              cout << "    NAN lighy yield with vdedx[i1] = "<<vdedx[i1]
              <<" and hiter->second->get_light_yield() = "<<hiter->second->get_light_yield() << endl;

            }
		    }
		  else
		    {
		      if(verbosity > 1)
			cout << "    did not find a previous entry for key = " << key << " create a new one" << endl;

		      cellptmap[key] = new PHG4CylinderCellv1();
		      it = cellptmap.find(key);
		      it->second->set_layer(*layer);
                      it->second->set_phibin(iphibin);
                      it->second->set_zbin(izbin);
		      it->second->add_edep(hiter->first, hiter->second->get_edep()*vdedx[i1], hiter->second->get_light_yield()*vdedx[i1]);
		      it->second->add_shower_edep(hiter->second->get_shower_id(), hiter->second->get_edep()*vdedx[i1]);

		      if(verbosity > 1 && std::isnan(hiter->second->get_light_yield()*vdedx[i1]))
            {

              cout << "    NAN lighy yield with vdedx[i1] = "<<vdedx[i1]
              <<" and hiter->second->get_light_yield() = "<<hiter->second->get_light_yield() << endl;

            }

		    }
		}
              vphi.clear();
              vz.clear();
	      
            } // end loop over hits

          int numcells = 0;

	  for(it = cellptmap.begin(); it != cellptmap.end(); ++it)
	    {
	      cells->AddCylinderCell(*layer, it->second);
	      numcells++;
	      if (verbosity > 1)
		{
		  cout << "Adding cell for key " << it->first << " in bin phi: " << it->second->get_binphi()
		       << " phi: " << geo->get_phicenter(it->second->get_binphi()) * 180./M_PI
		       << ", z bin: " << it->second->get_binz()
		       << ", z: " <<  geo->get_zcenter(it->second->get_binz())
		       << ", energy dep: " << it->second->get_edep()
		       << endl;
		}
	    }
	  
          if (verbosity > 0)
            {
              cout << "found " << numcells << " z/phi cells with energy deposition" << endl;
            }
        }

      //==========================================================
      // now reset the cell map before moving on to the next layer
      if(verbosity > 1)
	cout << "cellptmap for layer " << *layer << " has final length " << cellptmap.size();
      while(cellptmap.begin() != cellptmap.end())
	{
	  // Assumes that mmmory is freed by the cylinder cell container when it is destroyed
	  cellptmap.erase(cellptmap.begin());
	} 
      if(verbosity > 1)
	cout << " reset it to " << cellptmap.size() << endl;
    }
  if (chkenergyconservation)
    {
      CheckEnergy(topNode);
    }
  _timer.get()->stop();
  return Fun4AllReturnCodes::EVENT_OK;
}

int
PHG4CylinderCellReco::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void
PHG4CylinderCellReco::cellsize(const int i, const double sr, const double sz)
{
  set_size(i, sr, sz, PHG4CylinderCellDefs::sizebinning);
}

void
PHG4CylinderCellReco::etaphisize(const int i, const double deltaeta, const double deltaphi)
{
  set_size(i, deltaeta, deltaphi, PHG4CylinderCellDefs::etaphibinning);
  return;
}

void
PHG4CylinderCellReco::set_size(const int i, const double sizeA, const double sizeB, const int what)
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
PHG4CylinderCellReco::get_etaphi(const double x, const double y, const double z)
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
PHG4CylinderCellReco::get_eta(const double radius, const double z)
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

//---------------------------------------------------------------

bool PHG4CylinderCellReco::lines_intersect(
  double ax,
  double ay,
  double bx,
  double by,
  double cx,
  double cy,
  double dx,
  double dy,
  double* rx, // intersection point (output)
  double* ry
)
{

// Find if a line segment limited by points A and B
// intersects line segment limited by points C and D.
// First check if an infinite line defined by A and B intersects
// segment (C,D). If h is from 0 to 1 line and line segment intersect
// Then check in intersection point is between C and D

  double ex = bx - ax; // E=B-A
  double ey = by - ay;
  double fx = dx - cx; // F=D-C
  double fy = dy - cy;
  double px = -ey;     // P
  double py = ex;

  double bottom = fx * px + fy * py; // F*P
  double gx = ax - cx; // A-C
  double gy = ay - cy;
  double top = gx * px + gy * py; // G*P

  double h = 99999.;
  if (bottom != 0.)
    {
      h = top / bottom;
    }

//intersection point R = C + F*h
  if (h > 0. && h < 1.)
    {
      *rx = cx + fx * h;
      *ry = cy + fy * h;
      //cout << "      line/segment intersection coordinates: " << *rx << " " << *ry << endl;
      if ((*rx > ax && *rx > bx) || (*rx < ax && *rx < bx) || (*ry < ay && *ry < by) || (*ry > ay && *ry > by))
        {
          //cout << "       NO segment/segment intersection!" << endl;
          return false;
        }
      else
        {
          //cout << "       segment/segment intersection!" << endl;
          return true;
        }
    }

  return false;
}

//---------------------------------------------------------------

bool  PHG4CylinderCellReco::line_and_rectangle_intersect(
  double ax,
  double ay,
  double bx,
  double by,
  double cx,
  double cy,
  double dx,
  double dy,
  double* rr // length of the line segment inside the rectangle (output)
)
{

// find if a line isegment limited by points (A,B)
// intersects with a rectangle defined by two
// corner points (C,D) two other points are E and F
//   E--------D
//   |        |
//   |        |
//   C--------F

  if (cx > dx || cy > dy)
    {
      cerr << "ERROR: Bad rectangle definition!" << endl;
      return false;
    }

  double ex = cx;
  double ey = dy;
  double fx = dx;
  double fy = cy;
  double rx = 99999.;
  double ry = 99999.;

  vector<double> vx;
  vector<double> vy;

  bool i1 = lines_intersect(ax, ay, bx, by, cx, cy, fx, fy, &rx, &ry);
  if (i1)
    {
      vx.push_back(rx);
      vy.push_back(ry);
    }
  bool i2 = lines_intersect(ax, ay, bx, by, fx, fy, dx, dy, &rx, &ry);
  if (i2)
    {
      vx.push_back(rx);
      vy.push_back(ry);
    }
  bool i3 = lines_intersect(ax, ay, bx, by, ex, ey, dx, dy, &rx, &ry);
  if (i3)
    {
      vx.push_back(rx);
      vy.push_back(ry);
    }
  bool i4 = lines_intersect(ax, ay, bx, by, cx, cy, ex, ey, &rx, &ry);
  if (i4)
    {
      vx.push_back(rx);
      vy.push_back(ry);
    }

//cout << "Rectangle intersections: " << i1 << " " << i2 << " " << i3 << " " << i4 << endl;
//cout << "Number of intersections = " << vx.size() << endl;

  *rr = 0.;
  if (vx.size() == 2)
    {
      *rr = sqrt( (vx[0] - vx[1]) * (vx[0] - vx[1]) + (vy[0] - vy[1]) * (vy[0] - vy[1]) );
//  cout << "Length of intersection = " << *rr << endl;
    }
  if (vx.size() == 1)
    {
      // find which point (A or B) is within the rectangle
      if (ax > cx && ay > cy && ax < dx && ay < dy)   // point A is inside the rectangle
        {
          //cout << "Point A is inside the rectangle." << endl;
          *rr = sqrt((vx[0] - ax) * (vx[0] - ax) + (vy[0] - ay) * (vy[0] - ay));
        }
      if (bx > cx && by > cy && bx < dx && by < dy)   // point B is inside the rectangle
        {
          //cout << "Point B is inside the rectangle." << endl;
          *rr = sqrt((vx[0] - bx) * (vx[0] - bx) + (vy[0] - by) * (vy[0] - by));
        }
    }

  if (i1 || i2 || i3 || i4)
    {
      return true;
    }
  return false;
}


int
PHG4CylinderCellReco::CheckEnergy(PHCompositeNode *topNode)
{
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode, hitnodename.c_str());
  PHG4CylinderCellContainer *cells = findNode::getClass<PHG4CylinderCellContainer>(topNode, cellnodename);
  double sum_energy_g4hit = 0.;
  double sum_energy_cells = 0.;
  PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits();
  PHG4HitContainer::ConstIterator hiter;
  for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; hiter++)
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
  if (fabs(sum_energy_cells - sum_energy_g4hit)/sum_energy_g4hit > 1e-6) 
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
      cout << Name() << ":total energy for this event: " << sum_energy_g4hit << " GeV" << endl;
	}
    }
  return 0;
}

void
PHG4CylinderCellReco::Detector(const std::string &d)
{
  detector = d;
  // only set the outdetector nodename if it wasn't set already
  // in case the order in the macro is outdetector(); detector();
  if (outdetector.size() == 0)
    {
      OutputDetector(d);
    }
  return;
}
