#include "PHG4FullProjSpacalCellReco.h"

#include "PHG4CylinderGeomContainer.h"
#include "PHG4CylinderGeom_Spacalv3.h"
#include "PHG4CylinderCell_Spacalv1.h"
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
#include <phool/getClass.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <cassert>
#include <boost/foreach.hpp>
#include <exception>
#include <limits>       // std::numeric_limits
using namespace std;

PHG4FullProjSpacalCellReco::PHG4FullProjSpacalCellReco(const string &name) :
    SubsysReco(name), _timer(PHTimeServer::get()->insert_new(name.c_str())), chkenergyconservation(
        0), timing_window_size(numeric_limits<double>::max())
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
      cout
          << "PHG4FullProjSpacalCellReco::InitRun - Could not locate g4 hit node "
          << hitnodename << endl;
      topNode->print();
      exit(1);
    }
  cellnodename = "G4CELL_" + detector;
  PHG4CylinderCellContainer *cells = findNode::getClass<
      PHG4CylinderCellContainer>(topNode, cellnodename);
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
      PHIODataNode<PHObject> *newNode = new PHIODataNode<PHObject>(cells,
          cellnodename.c_str(), "PHObject");
      DetNode->addNode(newNode);
    }

  geonodename = "CYLINDERGEOM_" + detector;
  PHG4CylinderGeomContainer *geo =
      findNode::getClass<PHG4CylinderGeomContainer>(topNode,
          geonodename.c_str());
  if (!geo)
    {
      cout
          << "PHG4FullProjSpacalCellReco::InitRun - Could not locate geometry node "
          << geonodename << endl;
      topNode->print();
      exit(1);
    }
  if (verbosity > 0)
    {
      cout << "PHG4FullProjSpacalCellReco::InitRun - incoming geometry:"
          << endl;
      geo->identify();
      assert(geo->get_NLayers()>0);
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
      if (verbosity > 1)
        {
          layergeom->identify();
        }

      layergeom-> subtower_consistency_check();

//      int layer = layergeom->get_layer();

// create geo object and fill with variables common to all binning methods
      PHG4CylinderCellGeom_Spacalv1 *layerseggeo =
          new PHG4CylinderCellGeom_Spacalv1();
      layerseggeo->set_layer(layergeom->get_layer());
      layerseggeo->set_radius(layergeom->get_radius());
      layerseggeo->set_thickness(layergeom->get_thickness());
      layerseggeo->set_binning(PHG4CylinderCellDefs::spacalbinning);

      // construct a map to convert tower_ID into the older eta bins.

      const PHG4CylinderGeom_Spacalv3::tower_map_t & tower_map =
          layergeom->get_sector_tower_map();
      const PHG4CylinderGeom_Spacalv3::sector_map_t & sector_map =
          layergeom->get_sector_map();
      const int nphibin = layergeom->get_azimuthal_n_sec() // sector
      * layergeom->get_max_phi_bin_in_sec() // blocks per sector
          * layergeom->get_n_subtower_phi(); // subtower per block
      const double deltaphi = 2. * M_PI / nphibin;

      typedef map<double, int> map_z_tower_z_ID_t;
      map_z_tower_z_ID_t map_z_tower_z_ID;
      double phi_min = NAN;

      BOOST_FOREACH(const PHG4CylinderGeom_Spacalv3::tower_map_t::value_type& tower_pair, tower_map)
        {

          const int & tower_ID = tower_pair.first;
          const PHG4CylinderGeom_Spacalv3::geom_tower & tower =
              tower_pair.second;

          // inspect index in sector 0
          std::pair<int, int> tower_z_phi_ID = layergeom->get_tower_z_phi_ID(
              tower_ID, 0);

          const int & tower_ID_z = tower_z_phi_ID.first;
          const int & tower_ID_phi = tower_z_phi_ID.second;

          if (tower_ID_phi == 0)
            {
              //assign phi min according phi bin 0
              phi_min = atan2(tower.centralY,
                  tower.centralX
                      + 0.25
                          * (tower.pDx1 + tower.pDx2 + tower.pDx3 + tower.pDx4))
                  + sector_map.begin()->second;
            }

          if (tower_ID_phi == layergeom->get_max_phi_bin_in_sec() / 2)
            {
              //assign eta min according phi bin 0
              map_z_tower_z_ID[tower.centralZ] = tower_ID_z;
            }
          // ...
        } //       BOOST_FOREACH(const PHG4CylinderGeom_Spacalv3::tower_map_t::value_type& tower_pair, tower_map)

      assert(! std::isnan(phi_min));
      layerseggeo->set_phimin(phi_min);
      layerseggeo->set_phistep(deltaphi);
      layerseggeo->set_phibins(nphibin);

      PHG4CylinderCellGeom_Spacalv1::tower_z_ID_eta_bin_map_t tower_z_ID_eta_bin_map;
      int eta_bin = 0;
      BOOST_FOREACH( map_z_tower_z_ID_t::value_type& z_tower_z_ID, map_z_tower_z_ID)
        {
          tower_z_ID_eta_bin_map[z_tower_z_ID.second] = eta_bin;
          eta_bin++;
        }
      layerseggeo->set_tower_z_ID_eta_bin_map(tower_z_ID_eta_bin_map);
      layerseggeo->set_etabins(eta_bin * layergeom->get_n_subtower_eta());
      layerseggeo->set_etamin(NAN);
      layerseggeo->set_etastep(NAN);

      //build eta bin maps
      BOOST_FOREACH(const PHG4CylinderGeom_Spacalv3::tower_map_t::value_type& tower_pair, tower_map)
        {

          const int & tower_ID = tower_pair.first;
          const PHG4CylinderGeom_Spacalv3::geom_tower & tower =
              tower_pair.second;

          // inspect index in sector 0
          std::pair<int, int> tower_z_phi_ID = layergeom->get_tower_z_phi_ID(
              tower_ID, 0);
          const int & tower_ID_z = tower_z_phi_ID.first;
          const int & tower_ID_phi = tower_z_phi_ID.second;
          const int & etabin = tower_z_ID_eta_bin_map[tower_ID_z];

          if (tower_ID_phi == layergeom->get_max_phi_bin_in_sec() / 2)
            {
              // half z-range
              const double dz = fabs(
                  0.5 * (tower.pDy1 + tower.pDy2) * sin(tower.pRotationAngleX));

              const double eta_central = -log(
                  tan(0.5 * atan2(tower.centralY, tower.centralZ)));
              // half eta-range
              const double deta =
                  fabs(
                      eta_central
                          - (-log(
                              tan(
                                  0.5
                                      * atan2(tower.centralY,
                                          tower.centralZ + dz)))));

              for (int sub_tower_ID_y = 0; sub_tower_ID_y < tower.NSubtowerY;
                  ++sub_tower_ID_y)
                {
                  assert(tower.NSubtowerY <=layergeom->get_n_subtower_eta());
                  // do not overlap to the next bin.
                  const int sub_tower_etabin = etabin
                      * layergeom->get_n_subtower_eta() + sub_tower_ID_y;

                  layerseggeo->set_etabounds(sub_tower_etabin,
                      make_pair<double, double>(
                          eta_central - deta
                              + sub_tower_ID_y * 2 * deta / tower.NSubtowerY,
                          eta_central - deta
                              + (sub_tower_ID_y + 1) * 2 * deta
                                  / tower.NSubtowerY));
                  layerseggeo->set_zbounds(sub_tower_etabin,
                      make_pair<double, double>(
                          tower.centralZ - dz
                              + sub_tower_ID_y * 2 * dz / tower.NSubtowerY,
                          tower.centralZ - dz
                              + (sub_tower_ID_y + 1) * 2 * dz
                                  / tower.NSubtowerY));
                }

            }
          // ...
        } //       BOOST_FOREACH(const PHG4CylinderGeom_Spacalv3::tower_map_t::value_type& tower_pair, tower_map)

      // add geo object filled by different binning methods
      seggeo->AddLayerCellGeom(layerseggeo);
      if (verbosity >= VERBOSITY_SOME)
        {
          cout << "PHG4FullProjSpacalCellReco::InitRun::" << Name()
              << " - Done layer" << (layergeom->get_layer())
              << ". Print out the cell geometry:" << endl;
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
      cout
          << "PHG4FullProjSpacalCellReco::process_event - Fatal Error - Could not locate g4 hit node "
          << hitnodename << endl;
      exit(1);
    }
  PHG4CylinderCellContainer *cells = findNode::getClass<
      PHG4CylinderCellContainer>(topNode, cellnodename);
  if (!cells)
    {
      cout
          << "PHG4FullProjSpacalCellReco::process_event - Fatal Error - could not locate cell node "
          << cellnodename << endl;
      exit(1);
    }

  PHG4CylinderCellGeomContainer *seggeo = findNode::getClass<
      PHG4CylinderCellGeomContainer>(topNode, seggeonodename.c_str());
  if (!seggeo)
    {
      cout
          << "PHG4FullProjSpacalCellReco::process_event - Fatal Error - could not locate cell geometry node "
          << seggeonodename << endl;
      exit(1);
    }

  PHG4CylinderGeomContainer *layergeo = findNode::getClass<
      PHG4CylinderGeomContainer>(topNode, geonodename.c_str());
  if (!layergeo)
    {
      cout
          << "PHG4FullProjSpacalCellReco::process_event - Fatal Error - Could not locate sim geometry node "
          << geonodename << endl;
      exit(1);
    }

  PHG4HitContainer::LayerIter layer;
  pair<PHG4HitContainer::LayerIter, PHG4HitContainer::LayerIter> layer_begin_end =
      g4hit->getLayers();
  for (layer = layer_begin_end.first; layer != layer_begin_end.second; ++layer)
    {
      // layer loop
      PHG4HitContainer::ConstIterator hiter;
      PHG4HitContainer::ConstRange hit_begin_end = g4hit->getHits(*layer);

      PHG4CylinderCellGeom *geo_raw = seggeo->GetLayerCellGeom(*layer);
      PHG4CylinderCellGeom_Spacalv1 * geo =
          dynamic_cast<PHG4CylinderCellGeom_Spacalv1 *>(geo_raw);
      assert(geo);

      const PHG4CylinderGeom *layergeom_raw = layergeo->GetLayerGeom(*layer);
      assert(layergeom_raw);
      // a special implimentation of PHG4CylinderGeom is required here.
      const PHG4CylinderGeom_Spacalv3 *layergeom =
          dynamic_cast<const PHG4CylinderGeom_Spacalv3 *>(layergeom_raw);
      assert(layergeom);

      for (hiter = hit_begin_end.first; hiter != hit_begin_end.second; ++hiter)
        {
          // checking ADC timing integration window cut
          if (hiter->second->get_t(0) > timing_window_size)
            continue;

          // hit loop
          int scint_id = hiter->second->get_scint_id();

          unsigned int key = static_cast<unsigned int>(scint_id);
          if (celllist.find(key) == celllist.end())
            {
              // decode scint_id
              PHG4CylinderGeom_Spacalv3::scint_id_coder decoder(scint_id);

              // convert to z_ID, phi_ID
              std::pair<int, int> tower_z_phi_ID =
                  layergeom->get_tower_z_phi_ID(decoder.tower_ID,
                      decoder.sector_ID);
              const int & tower_ID_z = tower_z_phi_ID.first;
              const int & tower_ID_phi = tower_z_phi_ID.second;

              PHG4CylinderGeom_Spacalv3::tower_map_t::const_iterator it_tower =
                  layergeom->get_sector_tower_map().find(decoder.tower_ID);
              assert(it_tower != layergeom->get_sector_tower_map().end());

              // convert tower_ID_z to to eta bin number
              int etabin = -1;
              try
                {
                  etabin = geo->get_etabin_block(tower_ID_z); // block eta bin
                }
              catch (exception & e)
                {
                  cout << "Print cell geometry:" << endl;
                  geo->identify();
                  cout << "Print scint_id_coder:" << endl;
                  decoder.identify();
                  cout << "Print the hit:" << endl;
                  hiter->second->print();
                  cout << "PHG4FullProjSpacalCellReco::process_event::"
                      << Name() << " - Fatal Error - " << e.what() << endl;
                  exit(1);
                }

              const int sub_tower_ID_x = it_tower->second.get_sub_tower_ID_x(
                  decoder.fiber_ID);
              const int sub_tower_ID_y = it_tower->second.get_sub_tower_ID_y(
                  decoder.fiber_ID);

              celllist[key] = new PHG4CylinderCell_Spacalv1();
              celllist[key]->set_layer(*layer);
              celllist[key]->set_phibin(
                  tower_ID_phi * layergeom->get_n_subtower_phi()
                      + sub_tower_ID_x);
              celllist[key]->set_etabin(
                  etabin * layergeom->get_n_subtower_eta() + sub_tower_ID_y);
              celllist[key]->set_fiber_ID(decoder.fiber_ID);
            }

          celllist[key]->add_edep(hiter->first, hiter->second->get_edep(),
              hiter->second->get_light_yield());
          celllist[key]->add_shower_edep(hiter->second->get_shower_id(),
              hiter->second->get_edep());

        } // end loop over g4hits
      int numcells = 0;
      for (map<unsigned int, PHG4CylinderCell *>::const_iterator mapiter =
          celllist.begin(); mapiter != celllist.end(); ++mapiter)
        {
          cells->AddCylinderCell(*layer, mapiter->second);
          numcells++;
          if (verbosity > 1)
            {
              cout << "PHG4FullProjSpacalCellReco::process_event::" << Name()
                  << " - " << "Adding cell in bin eta "
                  << (mapiter->second->get_bineta()) << " phi "
                  << (mapiter->second->get_binphi()) << " fiber "
                  << (mapiter->second->get_fiber_ID()) << ", energy dep: "
                  << mapiter->second->get_edep() << ", light yield: "
                  << mapiter->second->get_light_yield() << endl;
            }
        }
      celllist.clear();
      if (verbosity > 0)
        {
          cout << "PHG4FullProjSpacalCellReco::process_event::" << Name()
              << " - " << " found " << numcells
              << " fibers with energy deposition" << endl;
        }
    }

  if (chkenergyconservation or verbosity > 4)
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
      cout
          << "PHG4FullProjSpacalCellReco::CheckEnergy - energy mismatch between cells: "
          << sum_energy_cells << " and hits: " << sum_energy_g4hit
          << " diff sum(cells) - sum(hits): "
          << sum_energy_cells - sum_energy_g4hit << endl;
      return -1;
    }
  else
    {
      if (verbosity > 0)
        {
          cout << "PHG4FullProjSpacalCellReco::CheckEnergy::" << Name()
              << " - total energy for this event: " << sum_energy_g4hit
              << " GeV. Passed CheckEnergy" << endl;
        }
    }
  return 0;
}

