#include "PHG4FullProjSpacalCellReco.h"

#include "PHG4CylinderGeomContainer.h"
#include "PHG4CylinderGeom_Spacalv3.h"
#include "PHG4CylinderCellGeomContainer.h"
#include "PHG4CylinderCellGeom.h"
#include "PHG4CylinderCellGeom_Spacalv1.h"
#include "PHG4CylinderCellContainer.h"
#include "PHG4CylinderCellDefs.h"

#include <g4main/PHG4Hit.h>
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
#include <cassert>
#include <boost/foreach.hpp>

using namespace std;

PHG4FullProjSpacalCellReco::PHG4FullProjSpacalCellReco(const string &name) :
    SubsysReco(name), _timer(PHTimeServer::get()->insert_new(name.c_str())), chkenergyconservation(
        0)
{
}

int
PHG4FullProjSpacalCellReco::InitRun(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);

  // Looking for the DST node
  PHCompositeNode *dstNode;
  dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode",
      "DST"));
  if (!dstNode)
    {
      std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
      exit(1);
    }
  hitnodename = "G4HIT_" + detector;
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode,
      hitnodename.c_str());
  if (!g4hit)
    {
      cout << "Could not locate g4 hit node " << hitnodename << endl;
      exit(1);
    }
  cellnodename = "G4CELL_" + detector;
  PHG4CylinderCellContainer *cells = findNode::getClass<
      PHG4CylinderCellContainer>(topNode, cellnodename);
  if (!cells)
    {
      cells = new PHG4CylinderCellContainer();
      PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(cells,
          cellnodename.c_str(), "PHObject");
      dstNode->addNode(newNode);
    }

  geonodename = "CYLINDERGEOM_" + detector;
  PHG4CylinderGeomContainer *geo =
      findNode::getClass<PHG4CylinderGeomContainer>(topNode,
          geonodename.c_str());
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
  PHG4CylinderCellGeomContainer *seggeo = findNode::getClass<
      PHG4CylinderCellGeomContainer>(topNode, seggeonodename.c_str());
  if (!seggeo)
    {
      seggeo = new PHG4CylinderCellGeomContainer();
      PHCompositeNode *runNode = dynamic_cast<PHCompositeNode*>(iter.findFirst(
          "PHCompositeNode", "RUN"));
      PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(seggeo,
          seggeonodename.c_str(), "PHObject");
      runNode->addNode(newNode);
    }

  map<int, PHG4CylinderGeom *>::const_iterator miter;
  pair<map<int, PHG4CylinderGeom *>::const_iterator,
      map<int, PHG4CylinderGeom *>::const_iterator> begin_end =
      geo->get_begin_end();
  map<int, std::pair<double, int> >::iterator sizeiter;
  for (miter = begin_end.first; miter != begin_end.second; ++miter)
    {
      const PHG4CylinderGeom *layergeom_raw = miter->second;
      assert(layergeom_raw);

      // a special implimentation of PHG4CylinderGeom is required here.
      const PHG4CylinderGeom_Spacalv3 *layergeom =
          dynamic_cast<const PHG4CylinderGeom_Spacalv3 *>(layergeom_raw);

      if (!layergeom)
        {
          cout << "PHG4FullProjSpacalCellReco::InitRun - Fatal Error -"
              << " require to work with a version of SPACAL geometry of PHG4CylinderGeom_Spacalv3 or higher. "
              << "However the incoming geometry has version "
              << layergeom_raw->ClassName() << endl;
          exit(1);
        }

      int layer = layergeom->get_layer();

      // create geo object and fill with variables common to all binning methods
      PHG4CylinderCellGeom_Spacalv1 *layerseggeo =
          new PHG4CylinderCellGeom_Spacalv1();
      layerseggeo->set_layer(layergeom->get_layer());
      layerseggeo->set_radius(layergeom->get_radius());
      layerseggeo->set_thickness(layergeom->get_thickness());

      if (verbosity > 1)
        {
          layergeom->identify();
        }

      layerseggeo->set_binning(PHG4CylinderCellDefs::spacalbinning);

      // construct a map to convert tower_ID into the older eta bins.

      const PHG4CylinderGeom_Spacalv3::tower_map_t & tower_map =
          layergeom->get_sector_tower_map();
      const PHG4CylinderGeom_Spacalv3::sector_map_t & sector_map =
          layergeom-> get_sector_map();
      const int nphibin = layergeom->get_azimuthal_n_sec() * layergeom->get_max_phi_bin_in_sec();
      const double deltaphi = 2. * M_PI / nphibin;

      map<double, int> map_z_tower_z_ID;
      double phi_min = NAN;

      BOOST_FOREACH(const PHG4CylinderGeom_Spacalv3::tower_map_t::value_type& tower_pair, tower_map)
        {

        const int & tower_ID = tower_pair.first;
        const PHG4CylinderGeom_Spacalv3::geom_tower & tower = tower_pair.second;

        // inspect index in sector 0
        std::pair<int,int> tower_z_phi_ID = layergeom->get_tower_z_phi_ID(tower_ID, 0);

        const int & tower_ID_z = tower_z_phi_ID.first;
        const int & tower_ID_phi = tower_z_phi_ID.second;

        if (tower_ID_phi == 0)
          {
            //assign phi min according phi bin 0
            const double phi_min = atan2( tower.centralY, tower.centralX ) + sector_map[0];
          }

        if (tower_ID_phi == layergeom->get_max_phi_bin_in_sec()/2)
          {
            //assign eta min according phi bin 0
            map_z_tower_z_ID[tower.centralZ] = tower_ID_z;

            //TODO: assign eta ranges.

          }
          // ...
        }//       BOOST_FOREACH(const PHG4CylinderGeom_Spacalv3::tower_map_t::value_type& tower_pair, tower_map)


      assert(not isnan(phi_min));

      PHG4CylinderCellGeom_Spacalv1::tower_z_ID_eta_bin_map_t tower_z_ID_eta_bin_map;
      int eta_bin=0;
      BOOST_FOREACH( map<double, int>::value_type& z_tower_z_ID, map_z_tower_z_ID)
      {
        tower_z_ID_eta_bin_map[z_tower_z_ID.second] = eta_bin;
        eta_bin++;
      }
      layerseggeo->set_tower_z_ID_eta_bin_map(tower_z_ID_eta_bin_map);

      if (verbosity > 1)
        {
          cout << "PHG4FullProjSpacalCellReco::InitRun - Tower mapping:"<<endl;

          BOOST_FOREACH(const PHG4CylinderCellGeom_Spacalv1::tower_z_ID_eta_bin_map_t::value_type&  tower_z_ID_eta_bin,
              layerseggeo->get_tower_z_ID_eta_bin_map())
          {
            cout <<"\t"<<"Tower Z ID["<<tower_z_ID_eta_bin.first<<"] \t-> Eta Bin "<<tower_z_ID_eta_bin.second<<endl;
          }
        }

      layerseggeo->set_etabins(22);
      layerseggeo->set_etamin(NAN);
      layerseggeo->set_etastep(NAN);

      layerseggeo->set_phimin(phi_min);
      layerseggeo->set_phistep(deltaphi);
      layerseggeo->set_phibins(nphibin);

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
PHG4FullProjSpacalCellReco::process_event(PHCompositeNode *topNode)
{
  _timer.get()->restart();
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode,
      hitnodename.c_str());
  if (!g4hit)
    {
      cout << "Could not locate g4 hit node " << hitnodename << endl;
      exit(1);
    }
  PHG4CylinderCellContainer *cells = findNode::getClass<
      PHG4CylinderCellContainer>(topNode, cellnodename);
  if (!cells)
    {
      cout << "could not locate cell node " << cellnodename << endl;
      exit(1);
    }

  PHG4CylinderCellGeomContainer *seggeo = findNode::getClass<
      PHG4CylinderCellGeomContainer>(topNode, seggeonodename.c_str());
  if (!seggeo)
    {
      cout << "could not locate geo node " << seggeonodename << endl;
      exit(1);
    }

  PHG4HitContainer::LayerIter layer;
  pair<PHG4HitContainer::LayerIter, PHG4HitContainer::LayerIter> layer_begin_end =
      g4hit->getLayers();
  for (layer = layer_begin_end.first; layer != layer_begin_end.second; ++layer)
    {
      PHG4HitContainer::ConstIterator hiter;
      PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits(*layer);
      PHG4CylinderCellGeom *geo = seggeo->GetLayerCellGeom(*layer);
      int nslatbins = geo->get_phibins();
      for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
        {
          int slatno = hiter->second->get_scint_id();
          int slatbin;
          slatbin = hiter->second->get_layer() / nslatscombined;
          if (slatbin < 0 || slatbin > nslatbins)
            {
              if (slatbin + 1 > nslatbins)
                {
                  if (verbosity > 0)
                    {
                      cout
                          << "dealing with non fitting slat binning, this one will be dropped"
                          << endl;
                    }
                  continue;
                }
              cout << "slatbin out of range: " << slatbin << ", slatbins 0 - "
                  << nslatbins << ", slat no: " << hiter->second->get_scint_id()
                  << endl;
            }

          unsigned int key = (slatbin << 16) + slatno;
          if (celllist.find(key) == celllist.end())
            {
              celllist[key] = new PHG4CylinderCellv1();
              celllist[key]->set_layer(*layer);
              celllist[key]->set_phibin(slatbin);
              celllist[key]->set_etabin(slatno);
            }

          celllist[key]->add_edep(hiter->first, hiter->second->get_edep());
        } // end loop over g4hits
      int numcells = 0;
      for (map<unsigned int, PHG4CylinderCell *>::const_iterator mapiter =
          celllist.begin(); mapiter != celllist.end(); ++mapiter)
        {
          cells->AddCylinderCell(*layer, mapiter->second);
          numcells++;
          if (verbosity > 1)
            {
              cout << "Adding cell in bin slat: " << (mapiter->first >> 16)
                  << ", eta: " << (mapiter->first & 0xFFFF) << ", energy dep: "
                  << mapiter->second->get_edep() << endl;
            }
        }
      celllist.clear();
      if (verbosity > 0)
        {
          cout << Name() << ": found " << numcells
              << " eta/slat cells with energy deposition" << endl;
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
PHG4FullProjSpacalCellReco::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

void
PHG4FullProjSpacalCellReco::etasize_nslat(const int i, const double deltaeta,
    const int nslat)
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
PHG4FullProjSpacalCellReco::set_size(const int i, const double sizeA,
    const int sizeB, const int what)
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
PHG4FullProjSpacalCellReco::CheckEnergy(PHCompositeNode *topNode)
{
  PHG4HitContainer *g4hit = findNode::getClass<PHG4HitContainer>(topNode,
      hitnodename.c_str());
  PHG4CylinderCellContainer *cells = findNode::getClass<
      PHG4CylinderCellContainer>(topNode, cellnodename);
  double sum_energy_g4hit = 0.;
  double sum_energy_cells = 0.;
  PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits();
  PHG4HitContainer::ConstIterator hiter;
  for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
    {
      sum_energy_g4hit += hiter->second->get_edep();
    }
  PHG4CylinderCellContainer::ConstRange cell_begin_end =
      cells->getCylinderCells();
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
          << " diff sum(cells) - sum(hits): "
          << sum_energy_cells - sum_energy_g4hit << endl;
      return -1;
    }
  else
    {
      if (verbosity > 0)
        {
          cout << Name() << ": total energy for this event: "
              << sum_energy_g4hit << " GeV" << endl;
        }
    }
  return 0;
}

double
PHG4FullProjSpacalCellReco::get_phi_slat_zero_low(const double radius,
    const double thickness, const double tiltangle)
{
  // A/sin(alpha) = C/sin(gamma)
  // beta = 90-gamma

  double sinalpha = ((radius + thickness / 2.) / radius) * sin(tiltangle);
  double beta = asin(sinalpha) - tiltangle;
  cout << "beta: " << beta * 180. / M_PI << endl;
  return beta;
}

double
PHG4FullProjSpacalCellReco::get_phi_slat_zero_up(const double radius,
    const double thickness, const double tiltangle)
{
  // A/sin(alpha) = C/sin(gamma)
  // beta = 90-gamma
  double a = radius + thickness;
  double c = radius + thickness / 2.;
  double alpha = M_PI - tiltangle;
  double singamma = c / a * sin(alpha);
  double beta = M_PI - asin(singamma) - alpha;
  cout << "beta: " << beta * 180. / M_PI << endl;
  return beta;
}
